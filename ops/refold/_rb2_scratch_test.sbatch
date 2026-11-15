#!/bin/bash
#SBATCH --job-name=rb2_scratch_test
#SBATCH --partition=gpu.120h
#SBATCH --qos=es_wicky/gpu/120
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=1
#SBATCH --gres=gpu:1
#SBATCH --nodelist=eu-g3-[071-077,083-084]
#SBATCH --exclude=eu-g6-064
#SBATCH --tmp=30G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_b2_scratchtest/msaoff_04fcaa8c/slurm.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_b2_scratchtest/msaoff_04fcaa8c/slurm.err
set -euo pipefail
export B2_SIF=/cluster/work/wicky/containers/boltz2.sif
export BOLTZ_CACHE=/cluster/work/wicky/params/.boltz
export B2_DIFFUSION_SAMPLES=1
export B2_OUTPUT_FORMAT=mmcif
export B2_WRITE_FULL_PAE=true
export RNA_MSA_PATH=""
export RNA_TEMPLATE_PATH=""
export RNA_TEMPLATE_CHAIN="B"
echo "TMPDIR=$TMPDIR ; node=$(hostname)"
bash /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/refold_b2_scratch.sh "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4/refold_b2_scratchtest/msaoff_04fcaa8c" "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v4" ACUGCCGUAUAGGCAG "04fcaa8c"
