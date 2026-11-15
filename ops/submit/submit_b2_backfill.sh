#!/usr/bin/env bash
# submit_b2_backfill.sh — create + submit B2 refold leaves for EXISTING midrun_v2
# designs that lack a refold_b2 leaf. Mirrors submit_esm_backfill.sh.
# Default partition gpu.bulk (PUBLIC Turing pool, preemptible) + gpumem:8g so the
# 11G 2080 Ti qualify (B2 needs ~2.7G). Self-compresses (worker packs designs/).
# Idempotent: skips uuids that already have a refold_b2 leaf.
#
#   submit_b2_backfill.sh                                          # dry-run, all stacks
#   submit_b2_backfill.sh --stack halludesign__ligandmpnn --go     # SUBMIT HD leaves
#   PART=gpu.4h GRES=gpumem:24g submit_b2_backfill.sh --go --limit 1
set -u
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
WORKER=$RB/code/scripts/stages/refold_b2.sh
RNASEQ=ACUGCCGUAUAGGCAG
PART="${PART:-gpu.bulk}"
GRES="${GRES:-gpumem:8g}"
WALL="${WALL:-24:00:00}"
GO=0; LIMIT=0; ONLY_UUID=""
STACKS="rfd3__ligandmpnn boltzdesign__ligandmpnn boltzgen__boltzif boltzgen__ligandmpnn halludesign__ligandmpnn"
while [ $# -gt 0 ]; do case "$1" in
  --go) GO=1;; --limit) LIMIT=$2; shift;; --uuid) ONLY_UUID=$2; shift;;
  --stack) STACKS=$2; shift;; --part) PART=$2; shift;;
  *) echo "unknown arg: $1" >&2; exit 2;; esac; shift; done
mk_sbatch(){ # $1=RDIR $2=SRC $3=UUID
  cat <<SB
#!/bin/bash
#SBATCH --job-name=refold_b2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --partition=$PART
#SBATCH --time=$WALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=$GRES
$([ -n "${NODELIST:-}" ] && echo "#SBATCH --nodelist=$NODELIST")
#SBATCH --output=$1/slurm.log
#SBATCH --error=$1/slurm.err
set -euo pipefail
export B2_SIF=/cluster/work/wicky/containers/boltz2.sif
export BOLTZ_CACHE=/cluster/work/wicky/params/.boltz
export B2_DIFFUSION_SAMPLES=1
export B2_OUTPUT_FORMAT=mmcif
export B2_WRITE_FULL_PAE=true
export RNA_MSA_PATH=""
export RNA_TEMPLATE_PATH=""
export RNA_TEMPLATE_CHAIN="B"
bash $WORKER $1 $2 $RNASEQ "$3"
SB
}
n=0; sub=0
for S in $STACKS; do
  M=$CAMP/$S/${LABEL:-midrun_v2}; [ -d "$M" ] || continue
  for slice in "$M"/designs/*.tar.zst; do
    [ -e "$slice" ] || continue
    u=$(basename "$slice" .tar.zst)
    [ -n "$ONLY_UUID" ] && [ "$u" != "$ONLY_UUID" ] && continue
    RDIR=$M/refold_b2/msaoff_$u
    [ -e "$RDIR" ] && continue
    n=$((n+1))
    if [ "$GO" = 1 ]; then
      mkdir -p "$RDIR"
      mk_sbatch "$RDIR" "$M" "$u" > "$RDIR/run.sbatch"
      jid=$(sbatch --parsable "$RDIR/run.sbatch" 2>/dev/null)
      echo "$S $u -> ${jid:-FAILED}"; sub=$((sub+1))
      [ "$LIMIT" -gt 0 ] && [ "$sub" -ge "$LIMIT" ] && break 2
    else
      echo "would submit: $S/$u"
      [ "$LIMIT" -gt 0 ] && [ "$n" -ge "$LIMIT" ] && break 2
    fi
  done
done
echo "---"
echo "$([ "$GO" = 1 ] && echo SUBMITTED || echo would-submit): $([ "$GO" = 1 ] && echo $sub || echo $n) B2 refold leaves   (partition=$PART gres=$GRES wall=$WALL)"
