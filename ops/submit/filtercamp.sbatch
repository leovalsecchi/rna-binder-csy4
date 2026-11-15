#!/bin/bash
#SBATCH --job-name=filtcamp
#SBATCH --partition=normal.4h
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/filtcamp_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
python3 filter_campaign.py \
  --master csy4_master_foot.csv \
  --store  filter_store_v2v9 \
  --out    csy4_filtered.csv \
  --t-dockq 0.7 --t-rmsd 3.0 --t-iptm 0.6 --t-repeat 6 --t-foot 5
