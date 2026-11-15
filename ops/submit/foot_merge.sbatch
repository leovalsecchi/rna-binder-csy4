#!/bin/bash
# Merge the foot-count sidecar output into one CSV. Run as a JOB, not a login-node process:
# each `ssh euler` lands on a DIFFERENT login node, so a nohup'd process there is invisible
# (and unkillable/uncheckable) from the next session. Shared storage + sbatch removes that.
#SBATCH --job-name=footmerge
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/footmerge_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
FC=$RB/analyses/cross_campaign/foot_clash
cd "$RB"
python3 "$FC/merge_foot_counts.py" "$FC/foot_counts_merged.csv"
echo "exit=$?"
ls -la "$FC/foot_counts_merged.csv"
