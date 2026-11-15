#!/bin/bash
# build_filter_store.sbatch — one streaming pass over the 9.7 GB master CSV to build the
# columnar store the binder-filter tool loads. Single-threaded and IO-bound (~40 min for
# 5.34M rows). Runs as a batch job because a long single-thread pass over Lustre gets killed
# on a login node.
#SBATCH --job-name=filter_store
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/filter_store_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"

SRC="$RB/analyses/cross_campaign/combined_clean_v2v8_fullmetrics.csv"
OUT="$RB/analyses/cross_campaign/filter_store_v2v8"

echo "[$(date '+%F %T')] host=$(hostname)"
echo "  src=$SRC"
echo "  out=$OUT"
rm -rf "$OUT.tmp"
PYTHONPATH=code/src:code/scripts python3 code/scripts/build_filter_store.py \
    --csv "$SRC" --out "$OUT.tmp" --n-hint 6000000
rc=$?
if [ $rc -ne 0 ]; then echo "BUILD FAILED rc=$rc"; exit $rc; fi

rm -rf "$OUT" && mv "$OUT.tmp" "$OUT"
echo "[$(date '+%F %T')] DONE -> $OUT"
du -sh "$OUT"

# tarball for transfer to the Mac
tar -C "$(dirname "$OUT")" -cf "$OUT.tar" "$(basename "$OUT")"
ls -la "$OUT.tar"
