#!/usr/bin/env bash
# midrun_v5 GEN — Turing stacks (RFD3/BoltzGen/BD1). Generates per-uuid run.sbatch via
# launch.py --num-runs N (creates N nested run dirs), then CANCELS the broken native
# --ntasks=4 array launch.py submits (we re-bucket with --cpus-per-task afterward).
# HD is generated separately (needs unique seed per uuid).
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v5
CLOG="$RB/ops/submit/mv5_gen_cancelled.txt"; : > "$CLOG"

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

# RFD3 (141): public gpu.24h, all 18 public 24G nodes (short ~50min jobs)
gen configs/rfd3__multif/default.yaml        141 gpu.24h  "eu-g3-[062-077,083-084]" gpu:1 23:30:00
# BoltzGen (152): gpu.24h x76 + gpu.120h x76
gen configs/boltzgen__multif/default.yaml     76 gpu.24h  "eu-g3-[062-070]"         gpu:1 23:30:00
gen configs/boltzgen__multif/default.yaml     76 gpu.120h "eu-g3-[071-077,083-084]" gpu:1 48:00:00
# BD1 (150): gpu.24h x75 + gpu.120h x75
gen configs/boltzdesign__multif/default.yaml  75 gpu.24h  "eu-g3-[062-070]"         gpu:1 23:30:00
gen configs/boltzdesign__multif/default.yaml  75 gpu.120h "eu-g3-[071-077,083-084]" gpu:1 48:00:00
echo "=== Turing gen done. Cancelled arrays: $(cat "$CLOG" | tr '\n' ' ') ==="
