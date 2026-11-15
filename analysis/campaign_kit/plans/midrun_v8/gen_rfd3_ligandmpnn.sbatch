#!/usr/bin/env bash
#SBATCH --job-name=midrun_v8_gen_rfd3_ligandmpnn
#SBATCH --partition=gpu.4h
#SBATCH --qos=es_wicky/gpu/4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=1
#SBATCH --gres=gpumem:8g
#SBATCH --array=0-151%200
#SBATCH --exclude=eu-g3-071
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/campaign_kit/midrun_v8/midrun_v8_gen_rfd3_ligandmpnn_%A_%a.log

set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true

# node-local scratch — extraction/outputs NEVER touch /cluster (inode-flood guardrail)
export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID}}"; mkdir -p "$TMPDIR"

# --- per-unit identity + inputs ---
export CK_INDEX="$SLURM_ARRAY_TASK_ID"
export CK_UUID="$(printf 'midrun_v8-gen_rfd3_ligandmpnn-%06d' "$CK_INDEX")"
export CK_SEED="$(( 80001 + CK_INDEX ))"
export CK_TARGET_PDB="/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/csy4_rna_only.pdb"
export CK_RNA="ACUGCCGUAUAGGCAG"
export CK_OUT="$TMPDIR/$CK_UUID"; mkdir -p "$CK_OUT"
PUBLISH="/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/gen_rfd3_ligandmpnn/midrun_v8"; mkdir -p "$PUBLISH"
[ -f "$PUBLISH/$CK_UUID.tar.zst" ] && { echo "already done: $CK_UUID"; exit 0; }   # idempotent/resumable

# --- run ONE unit (model run contract; must write into $CK_OUT) ---
apptainer run --nv --bind /cluster/work/wicky:/cluster/work/wicky  /cluster/work/wicky/lvalsecchi/rna_binder/containers/rfdiffusion3.sif python run_inference.py --diffusion_batch_size 8 --seed "$CK_SEED" --out "$CK_OUT"

# --- pack + verify + publish (generic): node-local -> Lustre as ONE verified file ---
if [ -n "$(ls -A "$CK_OUT" 2>/dev/null)" ]; then
  ( cd "$CK_OUT" && tar -cf - . | zstd -q -19 -T0 -o "$TMPDIR/$CK_UUID.tar.zst" )
  zstd -tq "$TMPDIR/$CK_UUID.tar.zst"
  mv -f "$TMPDIR/$CK_UUID.tar.zst" "$PUBLISH/$CK_UUID.tar.zst"
  echo "published $PUBLISH/$CK_UUID.tar.zst"
else
  echo "ERROR: $CK_OUT is empty — run_contract produced no output" >&2; exit 1
fi
