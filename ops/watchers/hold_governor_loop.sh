#!/bin/bash
#SBATCH --job-name=hold_governor_loop
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/watchers/_hold_governor_loop.log
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Self-chaining babysitter for hold_governor.py. Runs a pass every PERIOD
# seconds for ~3h, then resubmits itself (survives the 4h walltime). Stop the
# whole chain by:  touch /cluster/work/wicky/lvalsecchi/rna_binder/_hold_governor.STOP
set -u
RB=/cluster/work/wicky/lvalsecchi/rna_binder/ops/watchers
TARGET=${TARGET:-200}
PERIOD=${PERIOD:-90}
cd "$RB" || exit 1

if [ -f "$RB/_hold_governor.STOP" ]; then
  echo "[$(date +%T)] STOP sentinel present -> not starting"; exit 0
fi

END=$(( $(date +%s) + 3*3600 ))   # work ~3h, then chain
while [ "$(date +%s)" -lt "$END" ]; do
  if [ -f "$RB/_hold_governor.STOP" ]; then
    echo "[$(date +%T)] STOP sentinel -> exiting (no re-chain)"; exit 0
  fi
  python3 "$RB/hold_governor.py" --target "$TARGET" --apply
  sleep "$PERIOD"
done

# self-chain (only reached if no STOP)
sbatch "$RB/hold_governor_loop.sh"
