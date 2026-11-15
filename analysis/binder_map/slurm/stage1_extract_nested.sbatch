#!/bin/bash
# Stage 1 (nested): one array task per outer midrun archive.
#   sbatch --array=0-34%35 slurm/stage1_extract_nested.sbatch
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

cd "$FM" || exit 1
echo "=== stage1-nested shard=$SLURM_ARRAY_TASK_ID on $(hostname) $(date) ==="
command -v zstd >/dev/null || { echo "FATAL: zstd not on PATH"; exit 1; }

"$PY" scripts/stage1_extract_nested.py \
  --shard "$SLURM_ARRAY_TASK_ID" \
  --n-shards "${N_SHARDS:-35}" \
  --build "$FM/${BUILD:-build_shared}" \
  --results "$FM/${RESULTS:-results_shared}"
rc=$?
echo "=== exit $rc $(date) ==="
exit $rc
