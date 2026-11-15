#!/bin/bash
# sasa_chunk.sbatch — RNA-side dSASA moiety/groove decomposition (dsasa_rna_*, frac_dsasa_*)
# for one backend, over the lrfix chunks (which already carry all three refold path columns).
#
# WORKLIST MODE (preferred — Euler backfill caps starts at 50 jobs/user/partition per cycle,
# so scale by task SIZE, not task COUNT; a worklist also stops tasks burning a start on skips):
#   sbatch --array=0-<NTASK-1> \
#     --export=ALL,SA_CHUNKDIR=...,SA_OUTDIR=...,SA_PATHCOL=af3_msaoff_cif_path,\
# SA_WORKLIST=...,SA_NTASK=<NTASK> sasa_chunk.sbatch
#
# LEGACY PACK MODE (positional, processes SA_PACK consecutive indices per task):
#   sbatch --array=0-<N>%<throttle> --export=ALL,...,SA_PACK=5 sasa_chunk.sbatch
#SBATCH --job-name=sasa
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M
#SBATCH --tmp=100000
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/sasa_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF=/cluster/work/wicky/containers/boltzgen.sif
J="${SLURM_CPUS_PER_TASK:-8}"
: "${SA_CHUNKDIR:?}"; : "${SA_OUTDIR:?}"; : "${SA_PATHCOL:?}"

SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/sasa.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=2

T="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$SA_OUTDIR"

# Build this task's index list: strided slice of the worklist, or a positional pack.
if [ -n "${SA_WORKLIST:-}" ]; then
  NT="${SA_NTASK:?SA_NTASK required with SA_WORKLIST}"
  IDXS=$(awk -v t="$T" -v n="$NT" 'NR % n == t' "$SA_WORKLIST")
else
  PACK="${SA_PACK:-5}"
  IDXS=$(for k in $(seq 0 $((PACK-1))); do printf "%05d\n" $((T*PACK+k)); done)
fi
NIDX=$(printf "%s\n" "$IDXS" | grep -c . || true)
echo "[$(date '+%F %T')] task $T: $NIDX chunks, col=$SA_PATHCOL host=$(hostname) jobs=$J"

for IDX in $IDXS; do
  CIN="$SA_CHUNKDIR/chunk_${IDX}.csv"
  OUT="$SA_OUTDIR/sasa_${IDX}.csv"
  [ -s "$CIN" ] || continue
  if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then
    echo "[$(date '+%F %T')] $OUT complete — skip"; continue
  fi
  echo "[$(date '+%F %T')] chunk=$CIN"
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
    --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
    --env ARCH_READ_CACHE_MAX=2 --env PYTHONWARNINGS=ignore \
    "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
    code/scripts/contacts_sidecar.py --csv "$CIN" --output "$OUT" \
    --path-col "$SA_PATHCOL" --jobs "$J" --sasa \
    || echo "WARN rc=$? $CIN"
done
echo "[$(date '+%F %T')] DONE task $T ($NIDX chunks)"
