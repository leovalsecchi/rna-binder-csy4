#!/usr/bin/env bash
#SBATCH --job-name=rfd3_calib_finalize
#SBATCH --partition=normal.4h
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaign_kit/plans/_calib/finalize_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/campaign_kit
module load stack/2024-06 python/3.12.8 2>/dev/null || true
# reads the calibration job's sacct VRAM/runtime, derives, writes compat/euler.compat.yaml
python3 plans/_calib/finalize.py "${CALIB_JID:-7310716}" gpu
