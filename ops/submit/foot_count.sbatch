#!/bin/bash
# foot_chunk.sbatch — hairpin-foot clash (A-form + transplant) for ONE chunk, af3 + esm.
# Both backends share one node-local archive-extract cache. Continuous values only.
#   sbatch --array=0-2137%500 foot_chunk.sbatch
#SBATCH --job-name=footcount
#SBATCH --partition=normal.4h
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/footclash_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
FC=$RB/analyses/cross_campaign/foot_clash
cd "$RB"

# ---- HARD GUARD: archive extraction MUST be node-local, never /cluster Lustre ----
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/fc.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=3
export FOOT_EXT_MODEL=$FC/extended_csy4_model_0.cif

IDX=$(printf "%05d" "${SLURM_ARRAY_TASK_ID:-0}")
CIN="$FC/chunks/chunk_${IDX}.csv"
[ -s "$CIN" ] || { echo "missing $CIN"; exit 0; }
SIF=/cluster/work/wicky/containers/boltzgen.sif

for BE in af3:af3_msaoff_cif_path esm:esm_cif_path; do
  TAG="${BE%%:*}"; COL="${BE##*:}"
  OUT="$FC/out2/foot_${TAG}_${IDX}.csv"
  [ -s "$OUT" ] && { echo "skip existing $OUT"; continue; }   # idempotent on requeue
  singularity exec -B /cluster "$SIF" \
    python "$FC/foot_clash_sidecar2.py" "$CIN" "$COL" "$OUT.tmp" || { echo "FAIL $TAG $IDX"; continue; }
  mv -f "$OUT.tmp" "$OUT"
done
echo "done idx=$IDX on $(hostname)"
