#!/bin/bash
# hbbb_array.sbatch — base-vs-backbone split of the interface H-bonds, all three refolders.
#
# Fragmented deliberately: one array task = one (chunk, backend) pair, ~2000 structures, so
# tasks are ~2 h and backfill readily rather than sitting as a few long jobs. PyRosetta is
# the cost here (~4 s/structure across 8 workers), not IO.
#
# Verified before launch: hb_total reproduces the stored n_h_bonds_interface_af3 EXACTLY on
# 199/199 probe designs, so the split refines the production metric rather than replacing it.
#
#   sbatch --array=0-149%80 hbbb_array.sbatch
#SBATCH --job-name=hbbb
#SBATCH --partition=normal.4h
#SBATCH --time=03:50:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hbbb_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF=/cluster/work/wicky/containers/pyrosetta.sif
MAN=analyses/cross_campaign/hbbb_manifest.tsv
OUT=analyses/cross_campaign/hbbb_out
mkdir -p "$OUT"

SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/hbbb.XXXXXX);; esac
fi
[ -n "$SCRATCH" ] && [ -w "$SCRATCH" ] || { echo "ABORT: no node-local scratch"; exit 3; }
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=3

T="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$(sed -n "$((T+1))p" "$MAN")
[ -n "$LINE" ] || { echo "no manifest line $T"; exit 0; }
CH=$(echo "$LINE" | cut -f1)
BE=$(echo "$LINE" | cut -f2)
COL=$(echo "$LINE" | cut -f3)
O="$OUT/$(basename "${CH%.csv}")_${BE}.csv"
[ -s "$O" ] && { echo "$O exists — skip"; exit 0; }

echo "[$(date '+%F %T')] task $T  chunk=$CH  backend=$BE  host=$(hostname)  TMPDIR=$TMPDIR"
apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" \
  --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=3 --env PYTHONWARNINGS=ignore \
  "$SIF" env PYTHONPATH=code/src:code/scripts python3 code/scripts/hbond_base_backbone.py \
  --csv "$CH" --output "$O" --path-col "$COL" --jobs 8
echo "[$(date '+%F %T')] DONE rc=$? $O"
