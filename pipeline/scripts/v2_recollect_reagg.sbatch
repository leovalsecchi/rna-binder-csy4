#!/bin/bash
# v2_recollect_reagg.sbatch — re-run the v2 aggregation with the fixed aggregate.py
# (seq_key now includes gen+invfold, so BG+LiMPNN keeps its OWN sequence). Steps:
# re-collect(3 v2 stacks) -> clean -> sample10k metrics, then dispatch the PARALLEL
# full metrics + combine (full_metrics_parallel.sh with COMBINE_AFTER). This corrects
# the false cysteine + false global-dedup that were collapsing BG+LiMPNN into BG+BIF.
#SBATCH --job-name=v2_recollect_reagg
#SBATCH --partition=normal.24h
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v2_recollect_reagg_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"; export PATH="$RB/code/bin:$PATH"
S=code/scripts; C=campaigns/csy4_crrna_free; O=analyses/cross_v2_bg_bd
J="${SLURM_CPUS_PER_TASK:-32}"

echo "=== (1) re-collect v2 with FIXED aggregate.py $(date) ==="
python3 "$S/collect.py" \
    "$C/boltzgen__boltzif/midrun_v2" \
    "$C/boltzgen__ligandmpnn/midrun_v2" \
    "$C/boltzdesign__ligandmpnn/midrun_v2" \
    -o "$O/collected_raw.csv" --jobs "$J" --progress

echo "=== (2) clean + sample (BG+LiMPNN should stop collapsing into BG+BIF) $(date) ==="
python3 "$S/clean.py" --input "$O/collected_raw.csv" --outdir "$O"

echo "=== (3) sample10k metrics (for plots) $(date) ==="
python3 "$S/compute_metrics.py" --dir "$O" --jobs "$J" --progress

echo "=== (4) dispatch parallel FULL metrics + combine $(date) ==="
export COMBINE_AFTER="--v3 analyses/midrun_v3 --v2 $O --out analyses/cross_campaign"
sbatch --export=ALL "$S/full_metrics_parallel.sh" "$O" 32
echo "=== recollect+reagg prep DONE $(date); full-metrics array + combine dispatched ==="
