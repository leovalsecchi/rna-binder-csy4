#!/bin/bash
# Parameterised single-node stage runner.
#   sbatch --job-name=bm_s4_esm slurm/stage_single.sbatch <script> [args...]
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%j.log
set -uo pipefail            # NOT -e: report failures rather than vanishing

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

mkdir -p "$FM/logs"
echo "=== $SLURM_JOB_NAME on $(hostname) at $(date) ==="
echo "cmd: $PY $*"

cd "$FM" || exit 1
"$PY" "$@"
rc=$?
echo "=== exit $rc at $(date) ==="
exit $rc
