#!/bin/bash
#SBATCH --job-name=sse_bb
#SBATCH --array=0-28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:25:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/sse_bb/logs/%A_%a.log
SH=$(printf "shard_%02d" $SLURM_ARRAY_TASK_ID)
cd /cluster/work/wicky/lvalsecchi
./binder_map_venv/bin/python sse_worker.py $SH
