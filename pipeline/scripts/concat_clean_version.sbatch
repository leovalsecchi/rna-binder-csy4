#!/bin/bash
# concat_clean_version.sbatch — concat the per-stack collected_raw_<stack>.csv into one
# collected_raw.csv, then clean.py (cys-free + all-3-refolds + dedup ACROSS stacks).
# Run afterok on the 5 collect_stack jobs.
#   sbatch concat_clean_version.sbatch <midrun_vN>
#SBATCH --job-name=ccclean
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_cc_%x_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
S=code/scripts
V="${1:?usage: concat_clean_version.sbatch <midrun_vN>}"
OUT="analyses/$V"
PARTS=$(ls "$OUT"/collected_raw_*.csv 2>/dev/null || true)
[ -n "$PARTS" ] || { echo "FATAL: no collected_raw_*.csv in $OUT"; exit 1; }
echo "=== concat parts $(date) ==="; echo "$PARTS" | sed 's/^/  /'
# verify identical headers across parts (schema safety)
FIRST=$(echo "$PARTS" | head -1)
H=$(head -1 "$FIRST")
for f in $PARTS; do
  [ "$(head -1 "$f")" = "$H" ] || { echo "FATAL: header mismatch in $f"; exit 1; }
done
head -1 "$FIRST" > "$OUT/collected_raw.csv"
for f in $PARTS; do tail -n +2 "$f" >> "$OUT/collected_raw.csv"; done
echo "concat: $(( $(wc -l < "$OUT/collected_raw.csv") - 1 )) rows"
python3 "$S/clean.py" --input "$OUT/collected_raw.csv" --outdir "$OUT"
echo "=== clean DONE: $(( $(wc -l < "$OUT/collected_clean.csv") - 1 )) rows $(date) ==="
