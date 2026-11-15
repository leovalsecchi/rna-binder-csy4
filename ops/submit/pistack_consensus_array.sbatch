#!/bin/bash
#SBATCH --job-name=pistack_be
#SBATCH --array=0-19
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/pistack_be_%A_%a.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder

# ---- HARD GUARD: extraction scratch MUST be node-local, never /cluster ------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/pistack.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ] || case "$(readlink -f "$SCRATCH")" in /cluster/*) true;; *) false;; esac; then
  echo "ABORT: no node-local scratch (SCRATCH='$SCRATCH')"; exit 3
fi
export TMPDIR="$SCRATCH"
CH=$(printf 'analyses/cross_campaign/pistack_chunks/chunk_%02d.csv' "$SLURM_ARRAY_TASK_ID")
echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID -> $CH  TMPDIR=$TMPDIR  host=$(hostname)"

SIF=/cluster/work/wicky/containers/boltzgen.sif
run() {  # $1=path-col  $2=col-suffix ; each a SEPARATE process so scratch frees between backends
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
    "$SIF" env PYTHONPATH=code/src:code/scripts python3 code/scripts/pi_stack_sidecar.py \
    --csv "$CH" --path-col "$1" --col-suffix "$2" --jobs 8 --progress
}
run b2_cif_path  _b2
run esm_cif_path _esm
echo "[$(date '+%F %T')] DONE $CH"
