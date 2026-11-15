#!/bin/bash
#SBATCH --job-name=hd_shard
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%A_%a.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
# Sharded re-run of hd_cleanup for the two pack-heavy labels. One label per uuid dir
# takes ~2 min to tar+zstd (~500 small files each, Lustre metadata bound), so v8 (446
# dirs) and v9 (1236) are hours of serial work. Splitting the candidate list N ways
# turns that into one wall-clock shard.
#
# Safe to overlap with an already-running pass: pack_dir takes a non-blocking flock per
# uuid and no-ops on any dir that already has designs.tar.zst, so a shard that reaches a
# dir another task finished (or is holding) just moves on.
#
#   tasks  0-11  -> midrun_v8, shard i of 12
#   tasks 12-35  -> midrun_v9, shard i-12 of 24
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
T=$SLURM_ARRAY_TASK_ID
if [ "$T" -lt 12 ]; then L=midrun_v8; I=$T; N=12; else L=midrun_v9; I=$((T - 12)); N=24; fi
echo "task $T -> $L shard $I/$N   $(date)"
bash "$RB/ops/submit/hd_cleanup.sh" --apply --shard "$I/$N" "$L"
echo "done $(date)"
