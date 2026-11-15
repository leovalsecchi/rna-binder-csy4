#!/bin/bash
# v9_sidecar_chunk.sbatch — run a per-backend sidecar script (contacts/dsasa, pi-stack, ...)
# over PACK chunks of the shared fullmetrics chunk inputs, for af3 AND esm.
#
# One runner instead of a bespoke array per metric: the sidecars all take
# --csv/--output/--path-col and all live in boltzgen.sif. b2 is deliberately NOT run — the
# filter tool must never gain a Boltz-2 column.
#
#   sbatch --array=0-<N>%<throttle> --export=ALL,\
#     SD_CHUNKDIR=...,SD_OUTDIR=...,SD_SCRIPT=code/scripts/contacts_sidecar.py,\
#     SD_PREFIX=dsasa,SD_PACK=40,SD_EXTRA=--sasa v9_sidecar_chunk.sbatch
#SBATCH --job-name=v9side
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/v9side_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF="${SD_SIF:-/cluster/work/wicky/containers/boltzgen.sif}"
J="${SLURM_CPUS_PER_TASK:-8}"
: "${SD_CHUNKDIR:?set SD_CHUNKDIR}"; : "${SD_OUTDIR:?set SD_OUTDIR}"
: "${SD_SCRIPT:?set SD_SCRIPT}";     : "${SD_PREFIX:?set SD_PREFIX}"
EXTRA="${SD_EXTRA:-}"
PACK="${SD_PACK:-40}"

# ---- HARD GUARD: archive extraction MUST be node-local, never /cluster Lustre --------
# Extraction onto Lustre is what flooded the project inode quota on 2026-07-10.
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/sd.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=4

T="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$SD_OUTDIR"
for k in $(seq 0 $((PACK-1))); do
  IDX=$(printf "%05d" $((T*PACK+k)))
  CIN="$SD_CHUNKDIR/chunk_${IDX}.csv"
  [ -s "$CIN" ] || continue
  for BE in af3:af3_msaoff_cif_path esm:esm_cif_path; do
    TAG="${BE%%:*}"; COL="${BE##*:}"
    OUT="$SD_OUTDIR/${SD_PREFIX}_${TAG}_${IDX}.csv"
    # idempotent on requeue: a complete sidecar (row count matches its input) is kept
    if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then continue; fi
    echo "[$(date '+%F %T')] $TAG chunk=$IDX host=$(hostname)"
    # --col-suffix only exists on pi_stack_sidecar.py, NOT on contacts_sidecar.py —
    # passing it unconditionally would error out the latter. SD_SUFFIX=1 opts in.
    SFX=""; [ "${SD_SUFFIX:-0}" = "1" ] && SFX="--col-suffix _$TAG"
    apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
        --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" \
        --env ARCH_READ_CACHE_MAX=4 --env PYTHONWARNINGS=ignore \
        "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
        "$SD_SCRIPT" --csv "$CIN" --output "$OUT.tmp" --path-col "$COL" \
        $SFX --jobs "$J" $EXTRA \
        && mv -f "$OUT.tmp" "$OUT" || echo "WARN rc=$? $TAG $IDX"
  done
done
echo "[$(date '+%F %T')] DONE task $T (pack=$PACK)"
