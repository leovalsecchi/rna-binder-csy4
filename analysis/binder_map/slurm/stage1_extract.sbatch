#!/bin/bash
# Stage 1: extract backbones + fingerprints. One array task per shard.
#   sbatch --job-name=bm_s1_esm --array=0-63%16 slurm/stage1_extract.sbatch esm
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

LAYER="${1:?usage: stage1_extract.sbatch <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

cd "$FM" || exit 1
echo "=== stage1 layer=$LAYER shard=$SLURM_ARRAY_TASK_ID on $(hostname) $(date) ==="

"$PY" scripts/stage1_extract.py \
  --layer "$LAYER" \
  --shard "$SLURM_ARRAY_TASK_ID" \
  --n-shards "${N_SHARDS:-64}" \
  --build "$FM/${BUILD:-build}" \
  --results "$FM/${RESULTS:-results}"
rc=$?
echo "=== exit $rc $(date) ==="
exit $rc
