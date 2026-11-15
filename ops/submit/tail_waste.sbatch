#!/bin/bash
#SBATCH --job-name=tail_waste
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --tmp=200000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_tail_waste_%j.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder
export ARCH_READ_CACHE_MAX=3
apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" \
  --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=3 \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 scratch_read/tail_waste.py
