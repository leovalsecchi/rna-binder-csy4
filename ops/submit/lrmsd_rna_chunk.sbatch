#!/bin/bash
# lrmsd_rna_chunk.sbatch — RNA-aligned LRMSD (target held fixed, binder measured) vs the
# designed backbone, for PACK chunks of designs per array task.
#
# Reuses the shared fullmetrics chunk inputs (design_id + 3 refold path cols +
# designed_cif_path). Archives are decompressed to NODE-LOCAL scratch by archive_read, so
# nothing lands on Lustre and the project inode quota is untouched.
#
#   sbatch --array=0-<N>%<throttle> \
#     --export=ALL,LR_CHUNKDIR=...,LR_OUTDIR=...,LR_PACK=20 lrmsd_rna_chunk.sbatch
#SBATCH --job-name=lrmsdrna
#SBATCH --partition=normal.4h
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/lrmsdrna_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF=/cluster/work/wicky/containers/boltzgen.sif
J="${SLURM_CPUS_PER_TASK:-8}"
: "${LR_CHUNKDIR:?set LR_CHUNKDIR}"; : "${LR_OUTDIR:?set LR_OUTDIR}"

# ---- HARD GUARD: extraction MUST be node-local, never /cluster Lustre ----------------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/lr.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=5      # design archive + 3 refold archives live per work item

PACK="${LR_PACK:-1}"
T="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$LR_OUTDIR"
for k in $(seq 0 $((PACK-1))); do
  IDX=$(printf "%05d" $((T*PACK+k)))
  CIN="$LR_CHUNKDIR/chunk_${IDX}.csv"
  OUT="$LR_OUTDIR/lrmsd_${IDX}.csv"
  [ -s "$CIN" ] || continue
  if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then
    echo "[$(date '+%F %T')] $OUT complete — skip"; continue
  fi
  echo "[$(date '+%F %T')] chunk=$CIN host=$(hostname) TMPDIR=$TMPDIR jobs=$J"
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
    --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
    --env ARCH_READ_CACHE_MAX=5 --env PYTHONWARNINGS=ignore \
    "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
    code/scripts/lrmsd_rna_aligned.py --csv "$CIN" --output "$OUT" --jobs "$J" \
    || echo "WARN rc=$? $CIN"
done
echo "[$(date '+%F %T')] DONE task $T (pack=$PACK)"
