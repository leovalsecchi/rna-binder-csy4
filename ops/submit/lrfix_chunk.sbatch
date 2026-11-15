#!/bin/bash
# lrfix_chunk.sbatch — recompute lrmsd_rna_* / lrmsd_prot_* for ALL THREE backends with the
# register-corrected _match, over chunks cut from the lrfix copy of the v2v9 master.
#
#   sbatch --array=0-<N>%<throttle> \
#     --export=ALL,LR_CHUNKDIR=...,LR_OUTDIR=...,LR_PACK=5 lrfix_chunk.sbatch
#SBATCH --job-name=lrfix
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/lrfix_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF=/cluster/work/wicky/containers/boltzgen.sif
J="${SLURM_CPUS_PER_TASK:-8}"
: "${LR_CHUNKDIR:?set LR_CHUNKDIR}"; : "${LR_OUTDIR:?set LR_OUTDIR}"
SCRIPT="${LR_SCRIPT:-code/scripts/lrmsd_rna_aligned.py}"

# ---- HARD GUARD: extraction MUST be node-local, never /cluster Lustre ----------------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/lrfix.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=5

PACK="${LR_PACK:-5}"
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
  echo "[$(date '+%F %T')] chunk=$CIN host=$(hostname) jobs=$J"
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
    --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
    --env ARCH_READ_CACHE_MAX=5 --env PYTHONWARNINGS=ignore \
    "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
    "$SCRIPT" --csv "$CIN" --output "$OUT" --jobs "$J" \
    || echo "WARN rc=$? $CIN"
done
echo "[$(date '+%F %T')] DONE task $T (pack=$PACK)"
