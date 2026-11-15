#!/usr/bin/env bash
# submit_af3_backfill.sh — create + submit AF3 refold leaves for EXISTING midrun_v2
# designs that lack a refold_af3 leaf. Mirrors submit_esm_backfill.sh.
# Default partition gpuhe.bulk (wicky Ampere pool, preemptible) + gpumem:24g + the
# Ampere nodelist (AF3 needs sm_80; stage auto-picks xla on Turing but we keep it
# on Ampere here per ops policy). Self-compresses. Idempotent.
#
#   submit_af3_backfill.sh                                          # dry-run, all stacks
#   submit_af3_backfill.sh --stack halludesign__ligandmpnn --go     # SUBMIT HD leaves
#   PART=gpu.4h GRES=gpumem:24g NODELIST="" submit_af3_backfill.sh --go --limit 1   # Turing xla trial
set -u
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
WORKER=$RB/code/scripts/stages/refold_af3.sh
RNASEQ=ACUGCCGUAUAGGCAG
PART="${PART:-gpuhe.bulk}"
GRES="${GRES:-gpumem:24g}"
WALL="${WALL:-24:00:00}"
NODELIST="${NODELIST-eu-g4-[001-031],eu-g6-[001-080]}"
GO=0; LIMIT=0; ONLY_UUID=""
STACKS="rfd3__ligandmpnn boltzdesign__ligandmpnn boltzgen__boltzif boltzgen__ligandmpnn halludesign__ligandmpnn"
while [ $# -gt 0 ]; do case "$1" in
  --go) GO=1;; --limit) LIMIT=$2; shift;; --uuid) ONLY_UUID=$2; shift;;
  --stack) STACKS=$2; shift;; --part) PART=$2; shift;;
  *) echo "unknown arg: $1" >&2; exit 2;; esac; shift; done
mk_sbatch(){ # $1=RDIR $2=SRC $3=UUID
  local NL=""
  [ -n "$NODELIST" ] && NL="#SBATCH --nodelist=$NODELIST"
  cat <<SB
#!/bin/bash
#SBATCH --job-name=refold_af3
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --partition=$PART
#SBATCH --time=$WALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=$GRES
$NL
#SBATCH --output=$1/slurm.log
#SBATCH --error=$1/slurm.err
set -euo pipefail
export AF3_SIF=/cluster/work/wicky/containers/alphafold3.sif
export AF3_PARAMS=/cluster/work/wicky/params/af3_params
export MODEL_SEEDS_JSON='[1]'
export AF3_NUM_DIFFUSION_SAMPLES=1
export RNA_MSA_PATH=""
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
    RDIR=$M/refold_af3/msaoff_$u
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
echo "$([ "$GO" = 1 ] && echo SUBMITTED || echo would-submit): $([ "$GO" = 1 ] && echo $sub || echo $n) AF3 refold leaves   (partition=$PART gres=$GRES wall=$WALL nodelist=${NODELIST:-none})"
