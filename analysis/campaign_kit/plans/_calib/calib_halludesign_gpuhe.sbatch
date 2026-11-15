#!/usr/bin/env bash
#SBATCH --job-name=calib_halludesign_gpuhe
#SBATCH --partition=gpuhe.4h
#SBATCH --qos=es_wicky/gpuhe/4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:24g
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/campaign_kit/_calib/calib_halludesign_gpuhe_%j.log

set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true
export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID}}"; mkdir -p "$TMPDIR"
export CK_UUID="calib" CK_SEED="1" CK_INDEX="0"
export CK_TARGET_PDB="" CK_RNA=""
export CK_OUT="$TMPDIR/calib"; mkdir -p "$CK_OUT"
echo "[calibration] model=halludesign pool=gpuhe units=2 card=24G"

# run 2 unit(s) of the model — the probe reads this job's peak VRAM + runtime from sacct
cd /cluster/work/wicky/lvalsecchi/halludesign_upgrade && apptainer run --nv --bind /cluster/work/wicky:/cluster/work/wicky --env PYTHONPATH="/cluster/work/wicky/lvalsecchi/rna_binder/vendor/protenix_pkg" --env CUDA_HOME="/usr/local/cuda-real" /cluster/work/wicky/lvalsecchi/rna_binder/containers/af3_protenix_esm.sif python ./HalluDesign_run.py --HalluDesign_model cross_model --random_init_chain_spec A:80 --num_designs 1 --num_seqs 2 --num_recycles 3 --ref_time_steps 60 --mpnn ligand_mpnn --fix_chain_index B --seed "$CK_SEED" --out "$CK_OUT"

echo "[calibration] done; produced $(ls -A "$CK_OUT" 2>/dev/null | wc -l) output item(s)"
