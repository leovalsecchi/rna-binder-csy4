#!/usr/bin/env bash
# midrun_v6 GEN — Turing stacks (RFD3/BoltzGen/BD1). Generates per-uuid run.sbatch via
# launch.py --num-runs N (creates N nested run dirs), then CANCELS the broken native
# --ntasks=4 array launch.py submits (we re-bucket with --cpus-per-task afterward).
# HD is generated separately (needs unique seed per uuid).
#
# v6 change vs v5: route the stacks that FIT 4h to the high-cap gpu.4h tier
# (es_wicky/gpu/4 = 96 GPU / 1770G, ~2x the 24h cap). Measured mv5 gen runtimes:
#   RFD3 24m, BoltzGen 2.7h max -> gpu.4h (04:00:00, comfortable margin).
#   BD1 6.4h median / 11h max -> STAYS gpu.24h + gpu.120h (cannot fit 4h).
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v6
CLOG="$RB/ops/submit/mv6_gen_cancelled.txt"; : > "$CLOG"

gen(){  # cfg N part nodelist gres wall
  local cfg="$1" N="$2" part="$3" nodelist="$4" gres="$5" wall="$6"
  echo ">>> $cfg  N=$N  part=$part  nodes=$nodelist  gres=$gres  wall=$wall"
  local out aj
  out=$(python3 scripts/launch.py "$cfg" --label "$LABEL" --num-runs "$N" \
    --override slurm.gpus=1 --override slurm.gres="$gres" \
    --override slurm.partition="$part" --override slurm.time="$wall" \
    --override slurm.nodelist="$nodelist" 2>&1)
  aj=$(printf '%s\n' "$out" | grep -oE "jobid=[0-9]+" | grep -oE "[0-9]+" | tail -1)
  if [ -n "$aj" ]; then scancel "$aj" && echo "   cancelled broken array $aj"; echo "$aj" >> "$CLOG"
  else echo "   !! no array jobid parsed — output tail:"; printf '%s\n' "$out" | tail -3; fi
}

# RFD3 (141): gpu.4h, all 18 public 24G nodes (~24min jobs, huge 4h cap)
gen configs/rfd3__multif/default.yaml        141 gpu.4h   "eu-g3-[062-077,083-084]" gpu:1 04:00:00
# BoltzGen (152): gpu.4h, all 18 nodes (~2.3h median / 2.7h max -> fits 4h)
gen configs/boltzgen__multif/default.yaml    152 gpu.4h   "eu-g3-[062-077,083-084]" gpu:1 04:00:00
# BD1 (150): CANNOT fit 4h (6.4h median) -> gpu.24h x75 + gpu.120h x75 (as v5)
gen configs/boltzdesign__multif/default.yaml  75 gpu.24h  "eu-g3-[062-070]"         gpu:1 23:30:00
gen configs/boltzdesign__multif/default.yaml  75 gpu.120h "eu-g3-[071-077,083-084]" gpu:1 48:00:00
echo "=== Turing gen done. Cancelled arrays: $(cat "$CLOG" | tr '\n' ' ') ==="
