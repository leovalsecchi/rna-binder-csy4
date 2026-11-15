#!/bin/bash
#SBATCH --job-name=vrftruth
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrftruth_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"                     # archive refs in the master are relative to RB, not to cross_campaign
S=$(mktemp -d /scratch/vrf.XXXXXX); export TMPDIR="$S"; export ARCH_READ_CACHE_MAX=2
python3 analyses/cross_campaign/verify_dockq_truth.py \
  --master analyses/cross_campaign/csy4_master_foot.csv \
  --backend ${VT_BE:-af3} --n-per-stack ${VT_N:-6}
rm -rf "$S"
