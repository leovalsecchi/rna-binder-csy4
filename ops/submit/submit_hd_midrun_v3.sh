#!/usr/bin/env bash
# HD de-novo midrun_v3 — 4090 pool (eu-g6) across 4h/24h/120h + a100 (gpupr) precautioned 24h/120h.
# Per-uuid unique seed (de-novo determinism), --num-runs 1 loop, NO auto-refold.
# Logs a manifest (seed/uuid/jobid) for the later gap-sweep.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB/code"
MANIFEST="$RB/hd_midrun_v3_manifest.tsv"
[ -f "$MANIFEST" ] || printf "ts\tpool\tpartition\tseed\tuuid8\tjobid\n" > "$MANIFEST"
SEED=10001

sub () {  # nodelist gres partition wall count pool
  local nodelist="$1" gres="$2" part="$3" wall="$4" n="$5" pool="$6" k out jobid uuid
  for ((k=0; k<n; k++)); do
    out=$(python3 scripts/launch.py configs/halludesign/random_init.yaml \
      --label midrun_v3 --num-runs 1 \
      --override slurm.gpus=1 \
      --override slurm.gres="$gres" \
      --override slurm.partition="$part" \
      --override slurm.time="$wall" \
      --override slurm.nodelist="$nodelist" \
      --override generator_params.params.seed="$SEED" 2>&1)
    jobid=$(printf '%s' "$out" | grep -oE "jobid=[0-9]+" | head -1 | cut -d= -f2)
    uuid=$(printf '%s'  "$out" | grep -oE "uuid8=[0-9a-f]+" | head -1 | cut -d= -f2)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$(date +%FT%T)" "$pool" "$part" "$SEED" "${uuid:-NA}" "${jobid:-FAIL}" >> "$MANIFEST"
    [ -z "$jobid" ] && echo "[WARN seed=$SEED $pool $part] submit FAILED: $(printf '%s' "$out" | tail -1)"
    SEED=$((SEED+1))
  done
  echo ">>> done $pool $part : $n submitted (seed now $SEED)"
}

# ---- 4090 pool (eu-g6) : fast (2h17) + biggest pool ; fits 4h ----
sub "eu-g6-[001-080]" gpumem:24g gpuhe.4h   04:00:00  300 4090
sub "eu-g6-[001-080]" gpumem:24g gpuhe.24h  23:30:00  400 4090
sub "eu-g6-[001-080]" gpumem:24g gpuhe.120h 119:00:00 200 4090

# ---- a100 pool (gpupr) : PRECAUTIONS = 24h/120h only (3h49 too close to 4h cap),
#      pinned nodelist, gpumem:40g forces a100, modest count (slow + reservation-gated) ----
sub "eu-a65-[02-07],eu-ts-[01-02]" gpumem:40g gpupr.24h  23:30:00 12 a100
sub "eu-a65-[02-07],eu-ts-[01-02]" gpumem:40g gpupr.120h 119:00:00 12 a100

echo "=== DONE: 4090=900 (300/400/200), a100=24 (12/12); seeds 10001-$((SEED-1)). Manifest: $MANIFEST ==="
