#!/bin/bash
# pairrmsd_chunk.sbatch — pairwise refolder-agreement RMSD (af3/b2/esm, 3 pairs x
# {prot, rna, pose}) for ONE chunk of designs. Reuses the SHARED per-chunk inputs
# already written for the full-metrics run (chunk_<NN>.csv carries design_id + all
# three refold path columns), so no re-chunking of the 6 GB master is needed.
#
# Driven by env (set at submit): PR_CHUNKDIR (input chunks), PR_OUTDIR (sidecars).
#   sbatch --array=0-<N>%<throttle> \
#     --export=ALL,PR_CHUNKDIR=...,PR_OUTDIR=... pairrmsd_chunk.sbatch
#SBATCH --job-name=pairrmsd
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/pairrmsd_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF_BG=/cluster/work/wicky/containers/boltzgen.sif
J="${SLURM_CPUS_PER_TASK:-8}"
: "${PR_CHUNKDIR:?set PR_CHUNKDIR}"; : "${PR_OUTDIR:?set PR_OUTDIR}"

# ---- HARD GUARD: archive extraction MUST be node-local, never /cluster Lustre --------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/pr.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ] || { case "$(readlink -f "$SCRATCH")" in /cluster/*) true;; *) false;; esac; }; then
  echo "ABORT: no node-local scratch (SCRATCH='$SCRATCH')"; exit 3
fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=4      # 4 x ~100 MB x 8 workers; work item = one archive TRIPLE

PACK="${PR_PACK:-1}"                 # chunks handled per array task (scheduler overhead amortiser)
T="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$PR_OUTDIR"
for k in $(seq 0 $((PACK-1))); do
  IDX=$(printf "%05d" $((T*PACK+k)))
  CIN="$PR_CHUNKDIR/chunk_${IDX}.csv"
  OUT="$PR_OUTDIR/pairrmsd_${IDX}.csv"
  [ -s "$CIN" ] || continue
  # idempotent on requeue: a complete sidecar (same row count as its input) is kept
  if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then
    echo "[$(date +%F %T)] $OUT complete - skip"; continue
  fi
  echo "[$(date +%F %T)] chunk=$CIN host=$(hostname) TMPDIR=$TMPDIR jobs=$J"
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
    --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
    --env ARCH_READ_CACHE_MAX=4 --env PYTHONWARNINGS=ignore \
    "$SIF_BG" env PYTHONPATH=code/src:code/scripts python3 \
    code/scripts/refold_pair_rmsd.py --csv "$CIN" --output "$OUT" --jobs "$J" || echo "WARN rc=$? $CIN"
done
echo "[$(date +%F %T)] DONE task $T (pack=$PACK)"
