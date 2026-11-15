#!/usr/bin/env bash
# midrun_v6 GEN — HalluDesign (de-novo). UNCHANGED from v5 per user: gpuhe.24h + gpupr.24h,
# 23:30 walltime (HD is 3.6h median / 8.8h max -> does NOT fit 4h). MUST loop --num-runs 1
# with a UNIQUE seed per uuid (--num-runs N reuses one seed -> identical designs). Each call
# creates a flat run dir + submits a single job which we CANCEL (we only want run.sbatch;
# the campaign runs via the re-bucketed --cpus-per-task arrays). Ampere/Ada + A100 only.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v6
CFG=configs/halludesign/random_init.yaml
SEED=70001
MAN="$RB/ops/submit/mv6_hd_seed_manifest.txt"; : > "$MAN"

hd_lane(){  # N part nodelist gres
  local N="$1" part="$2" nodelist="$3" gres="$4"
  echo ">>> HD lane: N=$N part=$part nodes=$nodelist gres=$gres seed_start=$SEED"
  local k out uuid jid
  for ((k=0;k<N;k++)); do
    out=$(python3 scripts/launch.py "$CFG" --label "$LABEL" --num-runs 1 \
      --override slurm.gpus=1 --override slurm.gres="$gres" \
      --override slurm.partition="$part" --override slurm.time=23:30:00 \
      --override slurm.nodelist="$nodelist" \
      --override generator_params.params.seed="$SEED" 2>&1)
    uuid=$(printf '%s\n' "$out" | grep -oE "uuid8=[0-9a-f]+|UUID8: [0-9a-f]+" | grep -oE "[0-9a-f]{8}" | head -1)
    jid=$(printf '%s\n'  "$out" | grep -oE "jobid=[0-9]+" | grep -oE "[0-9]+" | head -1)
    [ -n "$jid" ] && scancel "$jid" 2>/dev/null
    printf "%s\t%s\t%s\t%s\n" "$SEED" "${uuid:-NA}" "${jid:-NA}" "$part" >> "$MAN"
    SEED=$((SEED+1))
  done
  echo "   lane done (seed now $SEED)"
}

hd_lane 400 gpuhe.24h "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g
hd_lane 200 gpupr.24h "eu-a65-[02-07],eu-ts-[01-02]"    gpumem:40g
echo "=== HD gen done: 600 run dirs, seeds 70001-$((SEED-1)), single jobs cancelled ==="
