#!/bin/bash
# v3_reaggregate.sbatch — re-run the midrun_v3 analysis (rfd3 + halludesign) after
# the RFD3 v3 ESM gap-fill lands, so analyses/midrun_v3 shows ESM closed, then
# refresh the cross-campaign pool that consumes v3. Submitted with
# --dependency=afterany on the ESM gap-fill arrays.
#SBATCH --job-name=v3_reaggregate
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v3_reaggregate_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
S=code/scripts
C=campaigns/csy4_crrna_free
J="${SLURM_CPUS_PER_TASK:-16}"

echo "=== (1) re-collect midrun_v3 (rfd3 + halludesign; picks up filled ESM/AF3) $(date) ==="
python3 "$S/collect.py" \
    "$C/rfd3__ligandmpnn/midrun_v3" \
    "$C/halludesign__ligandmpnn/midrun_v3" \
    -o analyses/midrun_v3/collected_raw.csv --jobs "$J" --progress

echo "=== (2) clean + sample $(date) ==="
python3 "$S/clean.py" --input analyses/midrun_v3/collected_raw.csv --outdir analyses/midrun_v3

echo "=== (3) sample10k structural metrics $(date) ==="
python3 "$S/compute_metrics.py" --dir analyses/midrun_v3 --jobs "$J" --progress

echo "=== (4) refresh cross-campaign pool (v3 + v2) $(date) ==="
python3 "$S/combine_and_rank.py" \
    --v3 analyses/midrun_v3 --v2 analyses/cross_v2_bg_bd \
    --out analyses/cross_campaign --top 10 --per-stack 10
python3 "$S/plot_cross.py" \
    analyses/cross_campaign/combined_sample10k_metrics.csv \
    analyses/cross_campaign/plots

echo "=== DONE $(date) — refreshed analyses/midrun_v3/ + analyses/cross_campaign/ ==="
