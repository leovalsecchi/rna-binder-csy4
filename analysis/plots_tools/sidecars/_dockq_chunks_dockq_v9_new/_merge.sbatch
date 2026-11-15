#!/bin/bash
#SBATCH --job-name=dqp_merge
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=analyses/plots_tools/sidecars/_dqp_merge_%j.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder
FIRST=$(ls analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new/chunk_*_dockq.csv 2>/dev/null | head -1)
[ -n "$FIRST" ] || { echo "FATAL: no chunk dockq csvs"; exit 1; }
head -1 "$FIRST" > "analyses/plots_tools/sidecars/dockq_v9_new.csv"
n=0; for f in analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new/chunk_*_dockq.csv; do tail -n +2 "$f" >> "analyses/plots_tools/sidecars/dockq_v9_new.csv"; n=$((n+1)); done
echo "merged $n chunks -> $(( $(wc -l < "analyses/plots_tools/sidecars/dockq_v9_new.csv") - 1 )) designs -> analyses/plots_tools/sidecars/dockq_v9_new.csv"
rm -rf "analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new"
echo "=== dqp merge DONE $(date) ==="
