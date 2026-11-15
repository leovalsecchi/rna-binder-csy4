#!/bin/bash
#SBATCH --job-name=rf_esm_scratch_test
#SBATCH --partition=gpu.4h
#SBATCH --qos=es_wicky/gpu/4
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=1
#SBATCH --gres=gpu:1
#SBATCH --nodelist=eu-g3-[071-077,083-084]
#SBATCH --tmp=30G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_esm_scratchtest/msaoff_04fcaa8c/slurm.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_esm_scratchtest/msaoff_04fcaa8c/slurm.err
set -euo pipefail
export ESM_SIF=/cluster/work/wicky/containers/esmc.sif
export ESM_HF_CACHE=/cluster/work/wicky/lvalsecchi/rna_binder/params/esmfold2_hf_cache/huggingface
export ESM_NUM_LOOPS=20
export ESM_NUM_SAMPLING_STEPS=100
export ESM_NUM_DIFFUSION_SAMPLES=1
export ESM_SEED=0
bash /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/refold_esm_scratch.sh "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_esm_scratchtest/msaoff_04fcaa8c" "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4" ACUGCCGUAUAGGCAG "04fcaa8c"
