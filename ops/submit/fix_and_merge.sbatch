#!/bin/bash
# fix_and_merge.sbatch — recompute the ipTM consensus on BOTH masters, then concatenate.
#
# Order matters: the aggregates are fixed per-table FIRST (each table's own legacy values are
# preserved in place), then the merge is a pure concat onto the v2-v7 schema. Doing it the
# other way round would bake the campaign-dependent old rule into the merged table.
#
# Every swap is: write .fixed -> compare line counts -> move original to the archive -> swap.
#SBATCH --job-name=fixmerge
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/fixmerge_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CC=$RB/analyses/cross_campaign
SIF=/cluster/work/wicky/containers/boltzgen.sif
run(){ apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp "$SIF" python3 "$@"; }

A=$CC/combined_clean_v2v7_fullmetrics.csv
B=$RB/analyses/midrun_v8/collected_clean_fullmetrics.csv
M=$CC/combined_clean_v2v8_fullmetrics.csv
ARCH=$CC/archive_pre_iptm_fix_20260805
mkdir -p "$ARCH"

echo "=== 1. fix ipTM aggregates on v2-v7 ==="
run "$RB/code/scripts/fix_iptm_master.py" --table "$A" --out "$A.fixed" || exit 1
echo "=== 2. fix ipTM aggregates on v8 ==="
run "$RB/code/scripts/fix_iptm_master.py" --table "$B" --out "$B.fixed" || exit 1

for f in "$A" "$B"; do
  o=$(wc -l < "$f"); n=$(wc -l < "$f.fixed")
  if [ "$o" -ne "$n" ]; then echo "ABORT: row mismatch $f ($o vs $n)"; exit 4; fi
  mv "$f" "$ARCH/$(basename "$f")"
  mv "$f.fixed" "$f"
  echo "swapped $(basename "$f") ($n lines incl header); original -> $ARCH"
done

echo "=== 3. merge onto the v2-v7 schema ==="
run "$RB/code/scripts/merge_masters.py" --a "$A" --b "$B" --out "$M" || exit 1

echo "=== verify ==="
ra=$(wc -l < "$A"); rb=$(wc -l < "$B"); rm_=$(wc -l < "$M")
echo "rows: v2v7 $((ra-1)) + v8 $((rb-1)) = $((ra+rb-2)) ; merged $((rm_-1))"
if [ $((ra+rb-2)) -ne $((rm_-1)) ]; then echo "ABORT: merged row count wrong"; exit 5; fi
echo "cols: v2v7 $(head -1 "$A" | tr ',' '\n' | wc -l) ; merged $(head -1 "$M" | tr ',' '\n' | wc -l)"
ls -la "$M"
echo DONE
