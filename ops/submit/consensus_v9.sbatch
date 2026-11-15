#!/bin/bash
# consensus_v9.sbatch — 3-backend worst-case consensus for midrun_v9, mirroring the v2v7
# run (job 8939470: -c 8 --mem-per-cpu=16G). Euler REJECTS --mem; use --mem-per-cpu.
#SBATCH --job-name=cons_v9
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_cons_v9_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
python3 code/scripts/consensus_integrate.py \
  --master   analyses/midrun_v9/collected_clean.csv \
  --chunkdir analyses/cross_campaign/fullmetrics_chunks_v9 \
  --out      analyses/midrun_v9/collected_clean_fullmetrics.csv
echo "=== DONE $(date) ==="
ls -la analyses/midrun_v9/
