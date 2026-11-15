#!/bin/bash
#SBATCH --job-name=appendv8top
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/appendv8top_%j.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder
python3 ops/submit/append_v8hdtopup.py
