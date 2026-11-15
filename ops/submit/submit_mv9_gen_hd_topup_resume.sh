#!/usr/bin/env bash
#SBATCH --job-name=hd_topup_gen
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/mv9_hd_topup_gen_resume.log
#
# RESUME of submit_mv9_gen_hd_topup.sh — the original died on the LOGIN NODE after 645/1000
# uuids (manifest froze at 12:10, seed 103045). Login nodes reap long-running loops, so this
# half runs as a BATCH JOB instead.
#
# Two deliberate differences from the original:
#   1. NO `: > "$MAN"` — this APPENDS to the existing manifest. Truncating would orphan the
#      645 run dirs already created (the allowlist bucketer reads this file).
#   2. Only the lanes that did not finish: gpuhe.bulk 255 of 720, then the three gpupr lanes.
# Seeds resume at 103046 and end at 103400, continuing the 102401- block with no reuse.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB/code"
LABEL=midrun_v9
CFG=configs/halludesign/random_init.yaml
SEED=103046
MAN="$RB/ops/submit/mv9_hd_topup_seed_manifest.txt"      # APPEND — never truncate

before=$(wc -l < "$MAN")
echo "resuming with $before uuids already in the manifest; first new seed $SEED"

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

hd_lane 255 gpuhe.bulk "eu-g4-[001-031],eu-g6-[001-080]" gpumem:24g 23:30:00
hd_lane  10 gpupr.24h  "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00
hd_lane  10 gpupr.120h "eu-a65-[02-03,05-07],eu-ts-[01-02]" gpumem:40g 47:30:00
hd_lane  80 gpupr.bulk "eu-a65-[02-07],eu-ts-[01-02]"       gpumem:40g 23:30:00

after=$(wc -l < "$MAN")
echo "=== resume done: manifest $before -> $after (last seed $((SEED-1))) ==="
