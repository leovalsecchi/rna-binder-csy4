#!/bin/bash
#SBATCH --job-name=pistack100k
#SBATCH --array=0-19
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/pistack_%A_%a.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder

# ---- HARD GUARD: extraction scratch MUST be node-local, never /cluster ------
# Pick a node-local scratch; refuse to run if it would land on shared Lustre
# (that is what flooded inodes on 2026-07-10). TMPDIR on Euler compute nodes is
# per-job node-local (/scratch/...) and auto-purged at job end.
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in
  /cluster/*|"") SCRATCH="" ;;                     # unusable -> try /scratch
esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/pistack.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ] || case "$(readlink -f "$SCRATCH")" in /cluster/*) true;; *) false;; esac; then
  echo "ABORT: no node-local scratch (SCRATCH='$SCRATCH' resolves under /cluster or unwritable)"; exit 3
fi
export TMPDIR="$SCRATCH"
CH=$(printf 'analyses/cross_campaign/pistack_chunks/chunk_%02d.csv' "$SLURM_ARRAY_TASK_ID")
echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID -> $CH  TMPDIR=$TMPDIR  host=$(hostname)"

# Pre-flight: confirm archive_read will extract OFF /cluster, else abort before any work.
apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env TMPDIR="$TMPDIR" \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 -c \
  'import os,sys; sys.path.insert(0,"code/scripts"); from archive_read import _node_local_base; b=_node_local_base() or ""; import os.path as p; ok=b and not p.realpath(b).startswith("/cluster"); print("scratch base:",b,"node-local:",bool(ok)); sys.exit(0 if ok else 4)'

# archive_read extracts to node-local $TMPDIR (our 2026-07-10 fix): no Lustre inodes.
apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
  --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/pi_stack_sidecar.py --csv "$CH" --jobs 8 --progress
echo "[$(date '+%F %T')] DONE $CH"
