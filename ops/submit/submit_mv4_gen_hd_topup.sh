#!/usr/bin/env bash
# midrun_v4 HD GEN TOPUP — 15 more de-novo uuids (seeds 50601-50615) to lift HD from
# 88,180 to ~90k, matching the other stacks. EXACT same proven recipe as
# submit_mv4_gen_hd.sh: loop --num-runs 1 with a UNIQUE seed per uuid, cancel the single
# job (keep only run.sbatch), Ampere/Ada + A100 only, 23:30 walltime. Writes a SEPARATE
# manifest so the original 600-uuid manifest is never touched.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v4
CFG=configs/halludesign/random_init.yaml
SEED=50601
MAN="$RB/ops/submit/mv4_hd_topup_seed_manifest.txt"; : > "$MAN"

hd_lane(){  # N part nodelist gres
  local N="$1" part="$2" nodelist="$3" gres="$4"
  echo ">>> HD topup lane: N=$N part=$part nodes=$nodelist gres=$gres seed_start=$SEED"
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

hd_lane 10 gpuhe.24h "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g
hd_lane 5  gpupr.24h "eu-a65-[02-07],eu-ts-[01-02]"    gpumem:40g
echo "=== HD topup gen done: 15 run dirs, seeds 50601-$((SEED-1)), single jobs cancelled ==="
