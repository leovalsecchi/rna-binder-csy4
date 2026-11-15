#!/bin/bash
# gapfill_merge.sbatch — replace 4 (campaign,stack) slices with fresh clean sets ->
# combined_clean_v2v7_filled.csv. Directives mirror merge_v2v7.sbatch (known-good).
#   sbatch gapfill_merge.sbatch
#SBATCH --job-name=gfmerge
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_gapfill_merge.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
G=analyses/gapfill
python3 code/scripts/gapfill_merge.py \
  analyses/cross_campaign/combined_clean_v2v7.csv \
  analyses/cross_campaign/combined_clean_v2v7_filled.csv \
  v2:$G/midrun_v2__rfd3__ligandmpnn/collected_clean.csv \
  v3:$G/midrun_v3__halludesign__ligandmpnn/collected_clean.csv \
  v4:$G/midrun_v4__halludesign__ligandmpnn/collected_clean.csv \
  v5:$G/midrun_v5__halludesign__ligandmpnn/collected_clean.csv
echo "=== gapfill merge DONE $(date) ==="
