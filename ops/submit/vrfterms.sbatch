#!/bin/bash
#SBATCH --job-name=vrfterms
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfterms_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
S=$(mktemp -d /scratch/vrf.XXXXXX); export TMPDIR="$S"; export ARCH_READ_CACHE_MAX=2
echo "################ rows the merge REWROTE (fnat_factor present) ################"
python3 analyses/cross_campaign/verify_dockq_terms.py \
  --master analyses/cross_campaign/csy4_master_foot.csv \
  --backend ${VT_BE:-af3} --n-per-stack ${VT_N:-3} --scan 999999999 --require-factor
echo
echo "################ rows the merge LEFT ALONE (any) ################"
python3 analyses/cross_campaign/verify_dockq_terms.py \
  --master analyses/cross_campaign/csy4_master_foot.csv \
  --backend ${VT_BE:-af3} --n-per-stack ${VT_N:-3} --scan 999999999
rm -rf "$S"
