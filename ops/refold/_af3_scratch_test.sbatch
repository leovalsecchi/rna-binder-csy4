#!/bin/bash
#SBATCH --job-name=rf_af3_scratch_test
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
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_af3_scratchtest/msaoff_04fcaa8c/slurm.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_af3_scratchtest/msaoff_04fcaa8c/slurm.err
set -euo pipefail
export AF3_SIF=/cluster/work/wicky/containers/alphafold3.sif
export AF3_PARAMS=/cluster/work/wicky/params/af3_params
export MODEL_SEEDS_JSON="[1]"
export AF3_NUM_DIFFUSION_SAMPLES=1
export RNA_MSA_PATH=""
bash /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/refold_af3_scratch.sh "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_af3_scratchtest/msaoff_04fcaa8c" "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4" ACUGCCGUAUAGGCAG "04fcaa8c"
