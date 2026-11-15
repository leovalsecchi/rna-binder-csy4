#!/bin/bash
#SBATCH --job-name=vrfdockq2
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfdockq2_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
python3 analyses/cross_campaign/verify_dockq2.py \
  analyses/cross_campaign/csy4_master_foot.csv
