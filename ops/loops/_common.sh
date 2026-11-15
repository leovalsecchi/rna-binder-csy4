#!/usr/bin/env bash
# _common.sh — shared paths + submit helper for the generalized loop runner.
# Sourced by loop.slurm (the runner) and loopctl (the control CLI).
RB=/cluster/work/wicky/lvalsecchi/rna_binder
LOOPS="$RB/ops/loops"
STATE="$LOOPS/state"
CONFD="$LOOPS/loops.d"
LOGDIR="$RB/ops/logs"
mkdir -p "$STATE"
ts(){ date '+%F %T'; }

# submit_loop <conf> : (re)submit loop.slurm for this conf. Every SBATCH option
# is derived from the conf, so the self-chain and `loopctl start` submit identically.
submit_loop(){
  local conf="$1"
  ( set -a; source "$conf"; set +a
    local log="$LOGDIR/_${NAME}_loop.log"
    local opts=(-J "${NAME}_loop" -o "$log" --open-mode=append \
                --partition=normal.4h --time=04:00:00 --mem-per-cpu=2G)
    [ "${CPUS:-1}" != "1" ] && opts+=(-c "$CPUS")
    [ "${REQUEUE:-1}" = "1" ] && opts+=(--requeue)
    sbatch "${opts[@]}" "$LOOPS/loop.slurm" "$conf"
  )
}
