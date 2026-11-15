#!/bin/bash
#SBATCH --job-name=sasadiv
#SBATCH --partition=normal.24h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/sasadiv_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
S=$(mktemp -d /scratch/sd.XXXXXX); export TMPDIR="$S"; export ARCH_READ_CACHE_MAX=2
python3 analyses/cross_campaign/sasa_diverge_probe.py \
  --master analyses/cross_campaign/csy4_master_foot.csv \
  --sidecar-dir analyses/cross_campaign/sasa_sidecars_af3 \
  --backend af3 --per-stack 4
rm -rf "$S"
