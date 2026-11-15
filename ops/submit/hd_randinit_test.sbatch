#!/bin/bash
#SBATCH --job-name=hd_randinit_test
#SBATCH --partition=gpuhe.4h
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:20g
#SBATCH --nodelist=eu-g4-[001-031],eu-g6-[001-080]
#SBATCH --time=04:00:00
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/randinit_test/slurm.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/randinit_test/slurm.err
set -euo pipefail
echo "=== Job started: $(date) ==="
echo "=== Host: $(hostname) ==="
nvidia-smi || true
export CONTAINER=/cluster/work/wicky/lvalsecchi/rna_binder/containers/af3_protenix_esm.sif
OUTDIR=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn/randinit_test/designs
rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"

apptainer run --nv \
  --home /cluster/work/wicky/lvalsecchi/rna_binder/containers/runtime_home \
  --bind /cluster/work/wicky:/cluster/work/wicky \
  --bind /cluster/software:/cluster/software:ro \
  --bind /cluster/software/manual/cuda/13.0.2/x86_64:/usr/local/cuda-real:ro \
  --env PYTHONPATH=/cluster/work/wicky/lvalsecchi/rna_binder/containers/protenix_pkg \
  --env XDG_CACHE_HOME=/cluster/work/wicky/lvalsecchi/rna_binder/containers/runtime_cache \
  --env XDG_CONFIG_HOME=/cluster/work/wicky/lvalsecchi/rna_binder/containers/runtime_config \
  --env MPLCONFIGDIR=/cluster/work/wicky/lvalsecchi/rna_binder/containers/runtime_cache/matplotlib \
  --env WANDB_MODE=offline \
  --env CUTLASS_PATH=/cluster/work/wicky/lvalsecchi/rna_binder/containers/cutlass \
  --env CUDA_HOME=/usr/local/cuda-real \
  --env DS_SKIP_CUDA_CHECK=1 \
  --env TMPDIR=/cluster/work/wicky/lvalsecchi/rna_binder/containers/runtime_home/tmp \
  --env LD_LIBRARY_PATH=/usr/local/cuda-real/lib64:/.singularity.d/libs \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  $CONTAINER \
  bash -c '
    export PATH=/usr/local/cuda-real/bin:$PATH
    cd /cluster/work/wicky/lvalsecchi/halludesign && \
    python ./HalluDesign_run.py \
      --HalluDesign_model cross_model \
      --template_path /cluster/work/wicky/lvalsecchi/rna_binder/targets/csy4_crrna_free/template_protenix.json \
      --extra_json_path /cluster/work/wicky/lvalsecchi/rna_binder/targets/csy4_crrna_free/template_af3.json \
      --random_init_chain_spec "A:80" \
      --num_designs 4 \
      --mpnn ligand_mpnn \
      --mpnn_num_packs 1 \
      --fix_chain_index B \
      --rna ACUGCCGUAUAGGCAG \
      --num_seqs 10 \
      --num_recycles 2 \
      --design_epoch_begin 0 \
      --ref_time_steps 60 \
      --output_dir '"$OUTDIR"'
  '
echo "=== done: $(date) ==="
