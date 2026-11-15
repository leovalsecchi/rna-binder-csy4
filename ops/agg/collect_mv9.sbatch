#!/bin/bash
# collect_mv9.sbatch — midrun_v9 collect + clean (steps 1+2 of run_analysis.sh).
# Step 3 (compute_metrics --dir) is DELIBERATELY SKIPPED: it looks for sample3k.csv
# while clean.py writes sample10k.csv, which is what made the v8 job (8959232) exit
# FAILED after producing perfectly good raw+clean CSVs.
# Snapshot semantics: clean.py requires all 3 refolds, so uuids whose af3/esm are
# still in flight simply drop out; re-run this later to pick them up.
#SBATCH --job-name=collect_mv9
#SBATCH --partition=normal.24h
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G   # was 2G: the 08-05 run peaked at 63.99 GiB of 64 (99.98%), and v9 has grown ~900 HD uuids since
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_collect_mv9_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"          # bundled zstd for packed (.tar.zst) runs

OUT=analyses/midrun_v9
C=campaigns/csy4_crrna_free
mkdir -p "$OUT"

RUNS=(
  "$C/rfd3__ligandmpnn/midrun_v9"
  "$C/boltzgen__ligandmpnn/midrun_v9"
  "$C/boltzgen__boltzif/midrun_v9"
  "$C/boltzdesign__ligandmpnn/midrun_v9"
  "$C/halludesign__ligandmpnn/midrun_v9"
)
for d in "${RUNS[@]}"; do [ -d "$d" ] || { echo "MISSING run dir: $d" >&2; exit 1; }; done

echo "=== (1) collect $(date) ==="
python3 code/scripts/collect.py "${RUNS[@]}" \
    -o "$OUT/collected_raw.csv" --jobs "${SLURM_CPUS_PER_TASK:-16}" --progress

echo "=== (2) clean + sample $(date) ==="
python3 code/scripts/clean.py --input "$OUT/collected_raw.csv" --outdir "$OUT"

echo "=== DONE $(date) ==="
ls -la "$OUT"
