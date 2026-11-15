#!/bin/bash
# merge_final.sbatch — fold the dSASA-moiety (af3/b2/esm) and fnat-correction sidecars into the
# master, then rewrite the DockQ columns the fnat fix touches. Single streaming pass.
#
#   sbatch --dependency=afterany:<sasa/fnatf array ids> merge_final.sbatch
#
# If a sidecar set is still incomplete when this fires, merge_final.py refuses to write and this
# script re-queues ITSELF behind whatever sasa/fnatf arrays are still in the queue (capped by
# MG_GEN) — so it lands on its own once the last chunk is in, without clobbering anything early.
#SBATCH --job-name=mergefinal
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/mergefinal_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CC="$RB/analyses/cross_campaign"
SELF="$RB/ops/submit/merge_final.sbatch"
cd "$CC"
GEN="${MG_GEN:-0}"
echo "[$(date '+%F %T')] merge_final gen=$GEN host=$(hostname)"

python3 merge_final.py \
  --master combined_clean_v2v9_fullmetrics_bf_final.csv \
  --out    csy4_master_complete.csv \
  --chunk-dir lrfix_chunks \
  --sasa-dir 'sasa_sidecars_{be}' \
  --fnat-dir fnatf_sidecars
RC=$?
if [ $RC -eq 0 ]; then
  echo "[$(date '+%F %T')] MERGE COMPLETE"
  exit 0
fi

echo "[$(date '+%F %T')] merge did not run (rc=$RC)"
if [ "$GEN" -ge 20 ]; then echo "gen cap reached — not re-queuing"; exit $RC; fi
DEPS=$(squeue -u "$USER" -r -h -o '%F %j' | awk '$2=="sasa"||$2=="fnatf"{print $1}' | sort -u | paste -sd: -)
if [ -n "$DEPS" ]; then
  echo "re-queuing behind $DEPS (gen $((GEN+1)))"
  sbatch --dependency=afterany:"$DEPS" --export=ALL,MG_GEN=$((GEN+1)) "$SELF"
else
  echo "no sasa/fnatf jobs left in the queue — sidecars are genuinely incomplete, resubmit them."
fi
exit $RC
