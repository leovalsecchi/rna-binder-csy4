#!/usr/bin/env bash
# midrun_v8 GEN — HalluDesign (de-novo). 2400 uuids, seeds 90001-92400.
# Same recipe as v7 (submit_mv7_gen_hd.sh): unique seed per uuid via --num-runs 1 loop;
# each launch.py creates a flat run dir + submits a single job which we CANCEL (we only
# want run.sbatch; the campaign runs via the re-bucketed arrays from bucket_midrun_v8.py).
# run.sbatch carries --resume (bulk-safe under preemption).
#
# Routing (user 2026-07-19): 80/20 gpuhe/gpupr, spread across 24h/120h/bulk in each pool.
# Tier is baked into each run.sbatch's #SBATCH partition/time/nodelist/gres, which the
# bucketer reads to form one array per tier. gpuhe.* partitions ALSO contain Turing eu-g3
# (sm_75) which HD CANNOT use -> pin eu-g4/eu-g6 (Ampere/Ada) only, exactly as v7.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v9
CFG=configs/halludesign/random_init.yaml
SEED=100001
MAN="$RB/ops/submit/mv9_hd_seed_manifest.txt"; : > "$MAN"

hd_lane(){  # N part nodelist gres wall
  local N="$1" part="$2" nodelist="$3" gres="$4" wall="$5"
  echo ">>> HD lane: N=$N part=$part nodes=$nodelist gres=$gres wall=$wall seed_start=$SEED"
  local k out uuid jid
  for ((k=0;k<N;k++)); do
    out=$(python3 scripts/launch.py "$CFG" --label "$LABEL" --num-runs 1 \
      --override slurm.gpus=1 --override slurm.gres="$gres" \
      --override slurm.partition="$part" --override slurm.time="$wall" \
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

# gpuhe 80% = 1920 (Ampere/Ada eu-g4/eu-g6, gpumem:24g)
hd_lane 216 gpuhe.24h  "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 23:30:00
hd_lane 216 gpuhe.120h "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 47:30:00
hd_lane 1728 gpuhe.bulk "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 23:30:00
# gpupr 20% = 480 (A100 eu-a65/eu-ts, gpumem:40g; gpupr.120h EXCLUDES eu-a65-04)
hd_lane 24 gpupr.24h  "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00
hd_lane 24 gpupr.120h "eu-a65-[02-03,05-07],eu-ts-[01-02]" gpumem:40g 47:30:00
hd_lane 192 gpupr.bulk "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00
echo "=== HD gen done: 2400 run dirs, seeds 90001-$((SEED-1)) ==="
