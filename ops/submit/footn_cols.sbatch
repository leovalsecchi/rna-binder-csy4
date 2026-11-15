#!/bin/bash
# Add footn_* columns to the shared filter store. Chained afterok to the merge job so it can
# only run on a COMPLETE merge — a partial CSV would leave designs NaN, which the filter treats
# as fail-closed, silently shrinking every result set.
#SBATCH --job-name=footncols
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000M
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/footncols_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
FC=$RB/analyses/cross_campaign/foot_clash
python3 "$RB/code/scripts/build_footn_cols.py" \
  "$RB/analyses/cross_campaign/filter_store_v2v9" \
  "$FC/foot_counts_merged.csv"
echo "exit=$?"
