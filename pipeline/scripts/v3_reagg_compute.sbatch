#!/bin/bash
# v3_reagg_compute.sbatch — stage 2: array task computes whole-set structural
# metrics on ONE chunk of collected_clean (compute_metrics --input/--output).
#SBATCH --job-name=v3_reagg_compute
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/midrun_v3/_metrics_chunks/logs/%A_%a.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"                                # CWD=$RB so apptainer binds code/pylibs (pydssp) + pyrosetta
export PATH="$RB/code/bin:$PATH"
CK="$RB/analyses/midrun_v3/_metrics_chunks"
IDX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_${IDX}.csv"
OUT="$CK/chunk_${IDX}_metrics.csv"
[ -s "$IN" ] || { echo "no chunk $IN"; exit 0; }
echo "=== chunk $IDX: $(( $(wc -l < "$IN") - 1 )) designs on $(hostname) $(date) ==="
python3 "$RB/code/scripts/compute_metrics.py" \
    --input "$IN" --output "$OUT" \
    --jobs "${SLURM_CPUS_PER_TASK:-8}" --progress
echo "=== chunk $IDX done -> $OUT ($(date)) ==="
