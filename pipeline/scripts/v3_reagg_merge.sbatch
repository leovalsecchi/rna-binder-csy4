#!/bin/bash
# v3_reagg_merge.sbatch — stage 3: concat the per-chunk metrics into the full
# collected_clean_metrics.csv, then sample10k metrics + cross-campaign refresh.
#SBATCH --job-name=v3_reagg_merge
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v3_reagg_merge_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
S=code/scripts
O=analyses/midrun_v3
CK="$O/_metrics_chunks"
J="${SLURM_CPUS_PER_TASK:-16}"

echo "=== (1) merge chunk metrics -> collected_clean_metrics.csv $(date) ==="
FULL="$O/collected_clean_metrics.csv"
FIRST=$(ls "$CK"/chunk_*_metrics.csv 2>/dev/null | head -1)
[ -n "$FIRST" ] || { echo "FATAL: no chunk metrics found"; exit 1; }
head -1 "$FIRST" > "$FULL"
n=0
for f in "$CK"/chunk_*_metrics.csv; do
    tail -n +2 "$f" >> "$FULL"; n=$((n+1))
done
echo "  merged $n chunks -> $(( $(wc -l < "$FULL") - 1 )) rows -> $FULL"

echo "=== (2) sample10k metrics (for plots) $(date) ==="
python3 "$S/compute_metrics.py" --dir "$O" --jobs "$J" --progress

echo "=== (3) refresh cross-campaign pool (v3 + v2) $(date) ==="
python3 "$S/combine_and_rank.py" \
    --v3 "$O" --v2 analyses/cross_v2_bg_bd \
    --out analyses/cross_campaign --top 10 --per-stack 10
python3 "$S/plot_cross.py" \
    analyses/cross_campaign/combined_sample10k_metrics.csv \
    analyses/cross_campaign/plots

echo "=== (4) cleanup chunks $(date) ==="
rm -rf "$CK"
echo "=== DONE $(date) — full CSV: $FULL + refreshed cross_campaign ==="
