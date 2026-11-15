#!/bin/bash
# fnat_chunk.sbatch — per-design fnat denominator correction factor (nat_full/nat_trim).
# Reads ONE structure per design (the design itself); no refolds, no trimmed files written.
#
# WORKLIST MODE (preferred — see sasa_chunk.sbatch header for why):
#   sbatch --array=0-<NTASK-1> \
#     --export=ALL,FF_CHUNKDIR=...,FF_OUTDIR=...,FF_WORKLIST=...,FF_NTASK=<NTASK> fnat_chunk.sbatch
#SBATCH --job-name=fnatf
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M
#SBATCH --tmp=100000
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/fnatf_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
J="${SLURM_CPUS_PER_TASK:-8}"
: "${FF_CHUNKDIR:?}"; : "${FF_OUTDIR:?}"
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/ff.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"; export ARCH_READ_CACHE_MAX=2
T="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$FF_OUTDIR"

if [ -n "${FF_WORKLIST:-}" ]; then
  NT="${FF_NTASK:?FF_NTASK required with FF_WORKLIST}"
  IDXS=$(awk -v t="$T" -v n="$NT" 'NR % n == t' "$FF_WORKLIST")
else
  PACK="${FF_PACK:-10}"
  IDXS=$(for k in $(seq 0 $((PACK-1))); do printf "%05d\n" $((T*PACK+k)); done)
fi
NIDX=$(printf "%s\n" "$IDXS" | grep -c . || true)
echo "[$(date '+%F %T')] task $T: $NIDX chunks, host=$(hostname) jobs=$J"

for IDX in $IDXS; do
  CIN="$FF_CHUNKDIR/chunk_${IDX}.csv"
  OUT="$FF_OUTDIR/fnatf_${IDX}.csv"
  [ -s "$CIN" ] || continue
  if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then
    echo "[$(date '+%F %T')] $OUT complete — skip"; continue
  fi
  echo "[$(date '+%F %T')] chunk=$CIN"
  PYTHONPATH=code/src:code/scripts python3 code/scripts/fnat_factor.py \
    --csv "$CIN" --output "$OUT" --jobs "$J" || echo "WARN rc=$? $CIN"
done
echo "[$(date '+%F %T')] DONE task $T ($NIDX chunks)"
