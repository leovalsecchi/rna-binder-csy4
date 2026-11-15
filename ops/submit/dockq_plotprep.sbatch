#!/bin/bash
#SBATCH --job-name=dockqprep
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/dockqprep_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
SIF=/cluster/work/wicky/containers/boltzgen.sif
P=$RB/analyses/cross_campaign/plotdata
cd "$RB"
run() {
  apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
    "$SIF" python3 code/scripts/dockq_plot_prep.py "$@"
}
run --master $RB/analyses/cross_campaign/combined_clean_v2v7_fullmetrics.csv \
    --campaign v2v7 --out $P/v2v7 --ids $P/seed2k_ids.txt
run --master $RB/analyses/midrun_v8/collected_clean_fullmetrics.csv \
    --campaign v8 --out $P/v8 --ids $P/seed2k_ids.txt
echo DONE
