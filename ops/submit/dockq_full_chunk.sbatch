#!/bin/bash
# dockq_full_chunk.sbatch — DockQ with ALL FOUR components per backend, one array task per
# chunk. Reuses the chunk inputs already written by dockq_parallel.sh (merge_chunks.sbatch
# deliberately did not delete them).
#
# Exists because dockq_refolds.py keeps only dockq_self and discards dockq_LRMSD /
# dockq_iRMSD / dockq_fnat — and dockq_LRMSD is the filter tool's sc_RMSD axis, while
# dockq_rna_recombine.py needs all three to rescore into the RNA-aligned frame.
#
#   sbatch --array=0-127%40 dockq_full_chunk.sbatch
#SBATCH --job-name=dqfull
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/dqfull_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CK=$RB/analyses/plots_tools/sidecars/_dockq_chunks_dockq_v9_new
cd "$RB"; export PATH="$RB/code/bin:$PATH"

IDX=$(printf "%03d" "${SLURM_ARRAY_TASK_ID:-0}")
IN="$CK/chunk_${IDX}.csv"
OUT="$CK/chunk_${IDX}_dqfull.csv"
[ -s "$IN" ] || { echo "empty $IN"; exit 0; }

# idempotent on requeue: keep a complete output (row count matches its input)
if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$IN")" ]; then
  echo "$OUT complete - skip"; exit 0
fi

echo "chunk $IDX: $(( $(wc -l < "$IN") - 1 )) designs on $(hostname) $(date)"
python3 code/scripts/dockq_refolds_full.py --input "$IN" --out "$OUT.tmp" \
        --jobs "${SLURM_CPUS_PER_TASK:-8}" \
  && mv -f "$OUT.tmp" "$OUT" \
  || { echo "FAIL chunk $IDX rc=$?"; exit 1; }
echo "done $IDX $(date)"
