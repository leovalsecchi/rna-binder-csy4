#!/bin/bash
#SBATCH --job-name=refold_af3
#SBATCH --array=0-439%50
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --signal=B:USR1@600
#SBATCH --partition=gpuhe.bulk
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:24g
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/midrun_v3/refold_af3/_array/%A_%a.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/midrun_v3/refold_af3/_array/%A_%a.err
set -uo pipefail

# --- array task -> uuid mapping (line N+1 of the list = task N) ---
CAMPAIGN=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/midrun_v3
LIST=/cluster/work/wicky/lvalsecchi/rna_binder/.cache/af3_array_uuids.txt
RNA=ACUGCCGUAUAGGCAG
UUID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")
[ -z "${UUID:-}" ] && { echo "no uuid for array task $SLURM_ARRAY_TASK_ID — exiting clean"; exit 0; }
RUNDIR=$CAMPAIGN/refold_af3/msaoff_$UUID
echo "=== array task $SLURM_ARRAY_TASK_ID -> uuid $UUID ($(date)) on $(hostname) ==="

# --- identical env + invocation to the per-uuid run.sbatch (verbatim) ---
export AF3_SIF=/cluster/work/wicky/containers/alphafold3.sif
export AF3_PARAMS=/cluster/work/wicky/params/af3_params
export MODEL_SEEDS_JSON='[1]'
export AF3_NUM_DIFFUSION_SAMPLES=1
export RNA_MSA_PATH=""
bash /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/refold_af3.sh "$RUNDIR" "$CAMPAIGN" "$RNA" "$UUID"
