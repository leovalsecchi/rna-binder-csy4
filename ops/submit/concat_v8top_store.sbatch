#!/bin/bash
#SBATCH --job-name=concatstore
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/concatstore_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; CC=$RB/analyses/cross_campaign
cd $RB
python3 code/scripts/concat_stores.py \
  --a $CC/filter_store_v2v9 \
  --b $CC/filter_store_v8hdtopup \
  --out $CC/filter_store_v2v9_hdtop --force
