#!/usr/bin/env bash
# submit_esm_backfill.sh — create + submit ESM refold leaves for EXISTING midrun_v2
# designs, mirroring the per-uuid refold_b2/msaoff_<uuid> leaves. Self-compresses
# (refold_esm.sh packs designs/ -> designs.tar.zst). Idempotent: skips uuids that
# already have a refold_esm leaf.
#
#   submit_esm_backfill.sh                      # dry-run, all stacks (counts only)
#   submit_esm_backfill.sh --stack boltzgen__ligandmpnn --uuid 0000df81 --go   # one leaf
#   submit_esm_backfill.sh --go                 # SUBMIT all missing leaves
#   submit_esm_backfill.sh --go --limit 50      # submit at most 50
#   PART=gpuhe.bulk submit_esm_backfill.sh --go # choose partition (default gpuhe.24h)
set -u
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
WORKER=$RB/code/scripts/stages/refold_esm.sh
RNASEQ=ACUGCCGUAUAGGCAG
PART="${PART:-gpuhe.24h}"
NODELIST="eu-g4-[001-031],eu-g6-[001-080]"   # ampere_24g (sm_80+); ESM Triton kernel needs Ampere/Ada
GO=0; LIMIT=0; ONLY_UUID=""
STACKS="rfd3__ligandmpnn boltzdesign__ligandmpnn boltzgen__boltzif boltzgen__ligandmpnn halludesign__ligandmpnn"
while [ $# -gt 0 ]; do case "$1" in
  --go) GO=1;; --limit) LIMIT=$2; shift;; --uuid) ONLY_UUID=$2; shift;;
  --stack) STACKS=$2; shift;; --part) PART=$2; shift;;
  *) echo "unknown arg: $1" >&2; exit 2;; esac; shift; done

mk_sbatch(){ # $1=RDIR $2=SRC $3=UUID
  cat <<SB
#!/bin/bash
#SBATCH --job-name=refold_esm
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --partition=$PART
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:24g
#SBATCH --nodelist=$NODELIST
#SBATCH --output=$1/slurm.log
#SBATCH --error=$1/slurm.err
set -euo pipefail
export ESM_SIF=/cluster/work/wicky/containers/esmc.sif
export ESM_HF_CACHE=$RB/params/esmfold2_hf_cache/huggingface
export ESM_NUM_LOOPS=20
export ESM_NUM_SAMPLING_STEPS=100
export ESM_NUM_DIFFUSION_SAMPLES=1
export ESM_SEED=0
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
    RDIR=$M/refold_esm/msaoff_$u
    [ -e "$RDIR" ] && continue          # already created -> skip (idempotent)
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
echo "$([ "$GO" = 1 ] && echo SUBMITTED || echo would-submit): $([ "$GO" = 1 ] && echo $sub || echo $n) ESM refold leaves   (partition=$PART)"
