#!/bin/bash
#SBATCH --job-name=vrfdockq
#SBATCH --partition=normal.4h
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfdockq_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
python3 verify_dockq.py csy4_master_foot.csv
