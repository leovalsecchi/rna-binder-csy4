#!/usr/bin/env bash
# midrun_v9 HD TOPUP — 1000 more uuids = ~150,000 designs (measured 150 designs/uuid,
# median==mean over 12 sampled v9 packs). Seeds 102401-103400: v9 HD used 100001-102400,
# so this is collision-free on the SAME label.
#
# Isolated topup per the established recipe: SEPARATE manifest (never `: >` the original)
# + an ALLOWLIST bucketer, so the 2400 existing v9 uuids are never resubmitted.
#
# Routing copied verbatim from submit_mv9_gen_hd.sh, scaled 2400 -> 1000 at the same
# proportions: 90/10 gpuhe/gpupr, and 80/10/10 bulk/24h/120h inside each pool.
# gpuhe.* also contains Turing eu-g3 (sm_75) which HD CANNOT use -> pin eu-g4/eu-g6 only.
# gpupr needs gpumem:40g (NOT 24g) - it is a discrete node-match key, not a VRAM floor.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v9
CFG=configs/halludesign/random_init.yaml
SEED=102401
MAN="$RB/ops/submit/mv9_hd_topup_seed_manifest.txt"; : > "$MAN"

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

# gpuhe 90% = 900 (Ampere/Ada eu-g4/eu-g6, gpumem:24g)
hd_lane  90 gpuhe.24h  "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 23:30:00
hd_lane  90 gpuhe.120h "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 47:30:00
hd_lane 720 gpuhe.bulk "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 23:30:00
# gpupr 10% = 100 (A100 eu-a65/eu-ts, gpumem:40g; gpupr.120h EXCLUDES eu-a65-04)
hd_lane  10 gpupr.24h  "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00
hd_lane  10 gpupr.120h "eu-a65-[02-03,05-07],eu-ts-[01-02]" gpumem:40g 47:30:00
hd_lane  80 gpupr.bulk "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00
echo "=== HD topup gen done: 1000 run dirs, seeds 102401-$((SEED-1)) ==="
echo "=== manifest: $MAN ==="
