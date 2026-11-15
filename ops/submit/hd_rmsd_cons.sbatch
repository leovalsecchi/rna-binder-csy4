#!/bin/bash
#SBATCH --job-name=hd_rmsdcons
#SBATCH --partition=normal.4h
#SBATCH --time=03:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_hd_rmsdcons_%j.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder
python3 scratch_read/hd_rmsd_consistency.py
