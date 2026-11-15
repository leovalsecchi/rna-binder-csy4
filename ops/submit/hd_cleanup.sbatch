#!/bin/bash
#SBATCH --job-name=hd_cleanup
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%A_%a.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
# Runs hd_cleanup.sh for one label per array task. Login nodes kill the big
# Lustre walks (v7/v8/v9 died silently there), so the walk belongs on a
# compute node. ARGS is passed through: empty = dry run, "--apply" = real.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
LABELS=(midrun_v2 midrun_v3 midrun_v4 midrun_v5 midrun_v6 midrun_v7 midrun_v8 midrun_v9)
L=${LABELS[$SLURM_ARRAY_TASK_ID]}
echo "task $SLURM_ARRAY_TASK_ID -> $L   args='${ARGS:-}'   $(date)"
bash "$RB/ops/submit/hd_cleanup.sh" ${ARGS:-} "$L"
echo "done $(date)"
