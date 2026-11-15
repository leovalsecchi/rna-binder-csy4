#!/bin/bash
#SBATCH --job-name=inodepar
#SBATCH --partition=normal.4h
#SBATCH --time=03:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/inodepar_%A_%a.log
# One array task per subtree — the census is filesystem-walk bound, so splitting the
# tree across nodes is the only thing that speeds it up.
set -uo pipefail
B=/cluster/work/wicky/lvalsecchi
C=$B/rna_binder/campaigns/csy4_crrna_free
OUT=$B/_inode_census; mkdir -p "$OUT"
TARGETS=(
  "$C/rfd3__ligandmpnn"
  "$C/boltzgen__ligandmpnn"
  "$C/boltzgen__boltzif"
  "$C/boltzdesign__ligandmpnn"
  "$C/halludesign__ligandmpnn"
  "$B/rna_binder/analyses"
  "$B/rna_binder/scratch_folds"
  "$B/rna_binder/seedtest_2026-07-30"
)
T=${SLURM_ARRAY_TASK_ID:-0}
D="${TARGETS[$T]:-}"
[ -n "$D" ] && [ -d "$D" ] || { echo "no target $T"; exit 0; }
N=$(basename "$D")
echo "[$(date '+%F %T')] task $T -> $D"
TOT=$(find "$D" 2>/dev/null | wc -l)
echo "TOTAL_INODES $N $TOT"
du --inodes -d 3 "$D" 2>/dev/null | sort -rn | head -40 > "$OUT/${T}_${N}.txt"
echo "--- top of $N ---"; head -18 "$OUT/${T}_${N}.txt"
echo "[$(date '+%F %T')] DONE $N total=$TOT"
