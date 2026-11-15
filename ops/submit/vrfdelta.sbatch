#!/bin/bash
#SBATCH --job-name=vrfdelta
#SBATCH --partition=normal.24h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfdelta_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
python3 verify_merge_delta.py \
  --pre  combined_clean_v2v9_fullmetrics_bf_final.csv \
  --post csy4_master_foot.csv
