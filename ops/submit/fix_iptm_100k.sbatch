#!/bin/bash
#SBATCH --job-name=fixiptm
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/fixiptm_%j.log
RB=/cluster/work/wicky/lvalsecchi/rna_binder
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif python3 $RB/code/scripts/fix_iptm_consensus_100k.py \
  --table $RB/analyses/cross_campaign/app_metrics_100k.csv \
  --backup-dir $RB/analyses/cross_campaign/archive_pre_iptm_fix_20260805
