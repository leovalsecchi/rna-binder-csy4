#!/bin/bash
# dockq_rna_v9.sbatch — recombine the RNA-aligned LRMSD into DockQ for midrun_v9,
# then swap it in as the default via the same safety-ordered driver used for v2v7/v8.
#SBATCH --job-name=dockqrna_v9
#SBATCH --partition=normal.4h
#SBATCH --time=03:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_dockqrna_v9_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
CC=analyses/cross_campaign
M=analyses/midrun_v9/collected_clean_fullmetrics.csv

echo "=== (1) recombine $(date) ==="
python3 code/scripts/dockq_rna_recombine.py \
  --sidecar-dir $CC/lrmsd_rna_v9 \
  --master      $M \
  --out         $CC/dockq_rna_v9.csv \
  --label       v9

echo "=== (2) apply (archives original, verifies row count) $(date) ==="
bash ops/submit/fold_dockq_rna.sh "$RB/$M" "$RB/$CC/dockq_rna_v9.csv" \
     "$RB/$CC/archive_pre_dockq_rna_v9"

echo "=== DONE $(date) ==="
echo "cols now: $(head -1 $M | tr "," "\n" | wc -l)  (v8 master = 144)"
