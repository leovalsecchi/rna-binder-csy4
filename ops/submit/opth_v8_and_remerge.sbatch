#!/bin/bash
# opth_v8_and_remerge.sbatch — give v8 the optimised-H columns, then rebuild the merged master.
#
# NOTHING IS RECOMPUTED. The optH sidecars (opth_full_{af3,b2,esm}.csv) already hold all
# 5,343,524 designs — the manifest builder unioned v2-v7 with v8 — but opth_integrate.py was
# only ever pointed at the v2-v7 master, which is why v8's 19 optH columns were blank.
# This just runs the same integration against the v8 master and redoes the concat.
#SBATCH --job-name=opthv8
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/opthv8_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CC=$RB/analyses/cross_campaign
SIF=/cluster/work/wicky/containers/boltzgen.sif
run(){ apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp "$SIF" python3 "$@"; }

A=$CC/combined_clean_v2v7_fullmetrics.csv
B=$RB/analyses/midrun_v8/collected_clean_fullmetrics.csv
M=$CC/combined_clean_v2v8_fullmetrics.csv

echo "=== v8 before: $(head -1 "$B" | tr ',' '\n' | wc -l) cols ==="
echo "=== 1. integrate optH into the v8 master ==="
run "$RB/code/scripts/opth_integrate.py" "$B" || exit 1
echo "=== v8 after:  $(head -1 "$B" | tr ',' '\n' | wc -l) cols ==="

echo "=== 2. re-merge onto the v2-v7 schema ==="
run "$RB/code/scripts/merge_masters.py" --a "$A" --b "$B" --out "$M.new" || exit 1

ra=$(wc -l < "$A"); rb=$(wc -l < "$B"); rn=$(wc -l < "$M.new")
echo "rows: v2v7 $((ra-1)) + v8 $((rb-1)) = $((ra+rb-2)) ; new merged $((rn-1))"
if [ $((ra+rb-2)) -ne $((rn-1)) ]; then echo "ABORT: merged row count wrong"; rm -f "$M.new"; exit 5; fi
mv "$M.new" "$M"
echo "cols: v2v7 $(head -1 "$A" | tr ',' '\n' | wc -l) ; merged $(head -1 "$M" | tr ',' '\n' | wc -l)"
ls -la "$M"
echo DONE
