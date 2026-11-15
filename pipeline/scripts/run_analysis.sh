#!/bin/bash
# run_analysis.sh — drive the collect -> clean -> compute_metrics pipeline into one
# tight analyses/<name>/ dir. Run on a compute node or `sbatch` it (collect reads a
# CIF per design for the sequence, so the full campaign is I/O-heavy).
#
#   sbatch code/scripts/run_analysis.sh <name> <run_dir> [run_dir ...]
#   bash   code/scripts/run_analysis.sh <name> <run_dir> [run_dir ...]
#
# Produces analyses/<name>/:
#   collected_raw.csv  collected_clean.csv  sample10k.csv  sample10k_metrics.csv  README.md
# Then for the whole set:  python3 code/scripts/compute_metrics.py --dir analyses/<name> --full
#SBATCH --job-name=rb_analysis
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_rb_analysis_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"          # bundled zstd for packed (.tar.zst) runs

NAME="${1:?usage: run_analysis.sh <name> <run_dir> [run_dir ...]}"; shift
[ "$#" -ge 1 ] || { echo "need at least one run dir" >&2; exit 1; }
OUT="analyses/$NAME"; S=code/scripts; JOBS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$OUT"

echo "=== (1) collect $(date) ==="
python3 "$S/collect.py" "$@" -o "$OUT/collected_raw.csv" --jobs "$JOBS" --progress

echo "=== (2) clean + sample $(date) ==="
python3 "$S/clean.py" --input "$OUT/collected_raw.csv" --outdir "$OUT"

echo "=== (3) compute_metrics on sample10k $(date) ==="
python3 "$S/compute_metrics.py" --dir "$OUT" --jobs "$JOBS" --progress

cat > "$OUT/README.md" <<EOF
# Analysis: $NAME

Generated $(date) by code/scripts/run_analysis.sh

## Inputs (campaign run dirs)
$(for d in "$@"; do echo "- $d"; done)

## Files
- \`collected_raw.csv\`   — one row per design; AF3-msaoff / B2 / ESM iPTM+PTM, refold CIF paths
- \`collected_clean.csv\` — cys-free, all-three-refolds, sequence-deduped (keep best AF3 iPTM)
- \`sample10k.csv\`        — random 600/stack (seed 0) of the clean set, for fast iteration
- \`sample10k_metrics.csv\`— sample10k + structural / RNA-contact / H-bond / tip-weighted-contact metrics

## Whole-set metrics (slower)
\`\`\`
python3 code/scripts/compute_metrics.py --dir $OUT --full --jobs $JOBS --progress
\`\`\`
EOF
echo "=== DONE $(date) -> $OUT ==="
ls -la "$OUT"
