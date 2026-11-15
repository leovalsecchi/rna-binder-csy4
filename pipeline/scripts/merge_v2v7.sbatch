#!/bin/bash
# merge_v2v7.sbatch — merge v6 + v7 clean into the v2-v5 master -> combined_clean_v2v7.csv.
# Directives copied from concat_clean_version.sbatch / v3_reaggregate.sbatch (known-good):
# system python3 + code/bin on PATH, no module load.
#   sbatch merge_v2v7.sbatch
#SBATCH --job-name=mergev2v7
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_merge_v2v7.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
python3 code/scripts/merge_combined_v2v7.py \
  analyses/cross_campaign/combined_clean.csv \
  analyses/cross_campaign/combined_clean_v2v7.csv \
  v6:analyses/midrun_v6/collected_clean.csv \
  v7:analyses/midrun_v7/collected_clean.csv
echo "=== merge DONE $(date) ==="
