#!/usr/bin/env bash
#SBATCH --job-name=calib_rfd3_gpupr
#SBATCH --partition=gpupr.4h
#SBATCH --qos=es_wicky/gpuhe/4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:80g
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/campaign_kit/_calib/calib_rfd3_gpupr_%j.log

set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true
export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID}}"; mkdir -p "$TMPDIR"
export CK_UUID="calib" CK_SEED="1" CK_INDEX="0"
export CK_TARGET_PDB="" CK_RNA=""
export CK_OUT="$TMPDIR/calib"; mkdir -p "$CK_OUT"
echo "[calibration] model=rfd3 pool=gpupr units=2 card=80G"

# run 2 unit(s) of the model — the probe reads this job's peak VRAM + runtime from sacct
apptainer run --nv --bind /cluster/work/wicky:/cluster/work/wicky  /cluster/work/wicky/lvalsecchi/rna_binder/containers/rfdiffusion3.sif python run_inference.py --diffusion_batch_size 1 --seed "$CK_SEED" --out "$CK_OUT"

echo "[calibration] done; produced $(ls -A "$CK_OUT" 2>/dev/null | wc -l) output item(s)"
