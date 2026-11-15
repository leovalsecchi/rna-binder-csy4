#!/bin/bash
# lrmsd_convention.sbatch — both LRMSD directions (protein-aligned = DockQ's, vs
# RNA-aligned = the binder-campaign convention) for the seedtest structures.
# One array task = one (backend, chunk) unit from vsdesign_manifest.tsv.
#   sbatch --array=0-28%29 lrmsd_convention.sbatch
#SBATCH --job-name=lrmsdconv
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=800M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/lrmsdconv_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
D=$RB/seedtest_2026-07-30
cd "$RB"
SIF=/cluster/work/wicky/containers/boltzgen.sif
J="${SLURM_CPUS_PER_TASK:-16}"
OUTDIR=$D/lrmsd_convention
mkdir -p "$OUTDIR"

SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/lc.XXXXXX);; esac
fi
[ -n "$SCRATCH" ] && [ -w "$SCRATCH" ] || { echo "ABORT: no node-local scratch"; exit 3; }
export TMPDIR="$SCRATCH"

T="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$(sed -n "$((T+1))p" "$D/vsdesign_manifest.tsv")
[ -n "$LINE" ] || { echo "no manifest line $T"; exit 0; }
BE=$(echo "$LINE" | cut -f1)
CHUNK=$(echo "$LINE" | cut -f2)
TAG="${BE}_${CHUNK%.tar.zst}"
OUT="$OUTDIR/${TAG}.csv"
[ -s "$OUT" ] && { echo "[$(date '+%F %T')] $OUT exists — skip"; exit 0; }

WORK="$TMPDIR/$TAG"
mkdir -p "$WORK"
echo "[$(date '+%F %T')] $TAG host=$(hostname) jobs=$J"
tar --zstd -xf "$D/outputs/$BE/$CHUNK" -C "$WORK" || { echo "extract failed"; exit 1; }

apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
  --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" --env PYTHONWARNINGS=ignore \
  "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/lrmsd_convention.py --dir "$WORK/designs" \
  --refs $RB/analyses/cross_campaign/designed_cifs_index.csv \
  --backend "$BE" --output "$OUT" --jobs "$J"
rc=$?
rm -rf "$WORK"
echo "[$(date '+%F %T')] DONE rc=$rc $OUT"
exit $rc
