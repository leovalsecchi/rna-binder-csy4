#!/bin/bash
# collect_mv8.sbatch — midrun_v8 RE-collect + clean (steps 1+2 of run_analysis.sh).
# Cloned verbatim from ops/agg/collect_mv9.sbatch (known-good, 3 COMPLETED runs).
#
# WHY: analyses/midrun_v8/collected_raw.csv was built 2026-07-29 15:19, but v8 HD's
# af3/esm refold leaves kept landing until 2026-08-06 17:08. clean.py requires all 3
# refolds, so the 07-29 snapshot dropped every HD design whose af3/esm arrived later:
# HD raw 71,400 -> clean 12,150 (17.0%) vs 95.6-100% for every other campaign.
# Spot-check confirms the refolds are COMPLETE on disk (gen 150 / af3 300 / b2 150 /
# esm 150 per uuid), so this is a stale-snapshot problem, NOT missing compute.
# v9 hit the identical issue and was fixed by exactly this re-run ("RE-COLLECT 08-07").
#
# NON-DESTRUCTIVE: writes to analyses/midrun_v8_recollect/, leaving the existing
# analyses/midrun_v8/ (the current master's source) untouched for comparison.
#
# Step 3 (compute_metrics --dir) stays DELIBERATELY SKIPPED: it looks for sample3k.csv
# while clean.py writes sample10k.csv — that is what made the earlier v8 job (8959232)
# exit FAILED despite producing perfectly good raw+clean CSVs.
#SBATCH --job-name=collect_mv8
#SBATCH --partition=normal.24h
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G   # v9 OOMed at 2G (peaked 63.99/64 GiB); v8 has MORE slices (3546 vs 2603)
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_collect_mv8_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"          # bundled zstd for packed (.tar.zst) runs

OUT=analyses/midrun_v8_recollect
C=campaigns/csy4_crrna_free
mkdir -p "$OUT"

RUNS=(
  "$C/rfd3__ligandmpnn/midrun_v8"
  "$C/boltzgen__ligandmpnn/midrun_v8"
  "$C/boltzgen__boltzif/midrun_v8"
  "$C/boltzdesign__ligandmpnn/midrun_v8"
  "$C/halludesign__ligandmpnn/midrun_v8"
)
for d in "${RUNS[@]}"; do [ -d "$d" ] || { echo "MISSING run dir: $d" >&2; exit 1; }; done

echo "=== (1) collect $(date) ==="
python3 code/scripts/collect.py "${RUNS[@]}" \
    -o "$OUT/collected_raw.csv" --jobs "${SLURM_CPUS_PER_TASK:-16}" --progress

echo "=== (2) clean + sample $(date) ==="
python3 code/scripts/clean.py --input "$OUT/collected_raw.csv" --outdir "$OUT"

# Report-only. Every command is guarded: under `set -e` a bare [ -f ] or a grep -c
# that finds 0 returns 1 and would exit the job FAILED *after* the CSVs were written
# correctly -- precisely the failure mode this script's step-3 comment warns about.
set +e
echo "=== (3) HD before/after $(date) ==="
for f in analyses/midrun_v8/collected_raw.csv "$OUT/collected_raw.csv" \
         analyses/midrun_v8/collected_clean.csv "$OUT/collected_clean.csv"; do
  if [ -f "$f" ]; then
    n=$(grep -c csy4free_HD_LiMPNN "$f")
    printf "  %-56s HD=%s\n" "$f" "${n:-0}"
  else
    printf "  %-56s (absent)\n" "$f"
  fi
done
set -e

echo "=== DONE $(date) ==="
ls -la "$OUT"
