#!/usr/bin/env bash
#SBATCH --job-name=calib_rfd3_gpu
#SBATCH --partition=gpu.bulk
#SBATCH --qos=es_wicky/gpu/bulk
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:8g
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/campaign_kit/_calib/calib_rfd3_gpu_%j.log

set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true
export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID}}"; mkdir -p "$TMPDIR"
export CK_UUID="calib" CK_SEED="1" CK_INDEX="0"
export CK_TARGET_PDB="" CK_RNA=""
export CK_OUT="$TMPDIR/calib"; mkdir -p "$CK_OUT"
echo "[calibration] model=rfd3 pool=gpu units=1 card=24G"

# run 1 unit(s) of the model — the probe reads this job's peak VRAM + runtime from sacct
source /cluster/work/wicky/lvalsecchi/rfdiffusion3/activate/foundry.activate.sh && cd /cluster/work/wicky/lvalsecchi/rfdiffusion3/examples/rna && frun rfd3 design out_dir="$CK_OUT" inputs=csy4_crrna.json n_batches=${CK_N_BATCHES:-1} diffusion_batch_size=${CK_DBS:-1} prevalidate_inputs=true skip_existing=false

echo "[calibration] done; produced $(ls -A "$CK_OUT" 2>/dev/null | wc -l) output item(s)"
