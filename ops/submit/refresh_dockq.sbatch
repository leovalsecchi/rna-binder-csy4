#!/bin/bash
#SBATCH --job-name=rfdockq
#SBATCH --partition=normal.4h
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/rfdockq_%j.log
set -uo pipefail
CC=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
cd "$CC"
python3 refresh_dockq_store.py --master csy4_master_foot.csv \
  --store filter_store_v2v9 --gate 0.7 ${RF_ARGS:-}
