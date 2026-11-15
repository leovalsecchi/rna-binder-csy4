#!/bin/bash
# v2_finalize.sbatch — after the v2 AF3/ESM gap refolds drain, re-collect the v2
# stacks (picks up the new _gap leaves), then rebuild the cross-campaign pool +
# top-10 + plots. Submitted with --dependency=afterany on the refold arrays.
#SBATCH --job-name=v2_finalize_crosscampaign
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v2_finalize_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"          # bundled zstd for packed runs
S=code/scripts
C=campaigns/csy4_crrna_free
J="${SLURM_CPUS_PER_TASK:-16}"

echo "=== (1) re-collect v2 (picks up new AF3/ESM _gap leaves) $(date) ==="
python3 "$S/collect.py" \
    "$C/boltzgen__boltzif/midrun_v2" \
    "$C/boltzgen__ligandmpnn/midrun_v2" \
    "$C/boltzdesign__ligandmpnn/midrun_v2" \
    -o analyses/cross_v2_bg_bd/collected_raw.csv --jobs "$J" --progress

echo "=== (2) clean + sample $(date) ==="
python3 "$S/clean.py" --input analyses/cross_v2_bg_bd/collected_raw.csv --outdir analyses/cross_v2_bg_bd

echo "=== (3) sample10k structural metrics $(date) ==="
python3 "$S/compute_metrics.py" --dir analyses/cross_v2_bg_bd --jobs "$J" --progress

echo "=== (4) combine v3+v2, re-rank consensus, re-export top-N $(date) ==="
python3 "$S/combine_and_rank.py" \
    --v3 analyses/midrun_v3 --v2 analyses/cross_v2_bg_bd \
    --out analyses/cross_campaign --top 10 --per-stack 10

echo "=== (5) re-plot $(date) ==="
python3 "$S/plot_cross.py" \
    analyses/cross_campaign/combined_sample10k_metrics.csv \
    analyses/cross_campaign/plots

echo "=== DONE $(date) — refreshed analyses/cross_campaign/ ==="
