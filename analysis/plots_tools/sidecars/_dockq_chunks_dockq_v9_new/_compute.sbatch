#!/bin/bash
#SBATCH --job-name=dqp_compute
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new/logs/%A_%a.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder; export PATH="/cluster/work/wicky/lvalsecchi/rna_binder/code/bin:$PATH"
IDX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
IN="analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new/chunk_${IDX}.csv"; OUT="analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new/chunk_${IDX}_dockq.csv"
[ -s "$IN" ] || { echo "empty $IN"; exit 0; }
[ -f "$OUT" ] && { echo "$OUT exists, skip (resume)"; exit 0; }
echo "chunk $IDX: $(( $(wc -l < "$IN") - 1 )) designs on $(hostname)"
python3 /cluster/work/wicky/lvalsecchi/rna_binder/analyses/plots_tools/dockq_refolds.py --input "$IN" --out "$OUT" --jobs "${SLURM_CPUS_PER_TASK:-8}"
