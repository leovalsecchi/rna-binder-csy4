#!/usr/bin/env bash
# HD de-novo midrun_v3 — 3090 pool (eu-g4, EXCL flaky eu-g4-012) on 24h + 120h queues.
# Completes the 1200-uuid campaign: seeds 10925-11200 (276 uuids).
# 3090 is ~2-3x slower than 4090 (de-novo ~5-7h) -> generous walltimes (no timeout).
# eu-g4 has idle nodes so jobs start immediately regardless of walltime.
# Config already at mem_per_cpu 7G (28G). NO auto-refold. Appends shared manifest.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB/code"
MANIFEST="$RB/hd_midrun_v3_manifest.tsv"
[ -f "$MANIFEST" ] || printf "ts\tpool\tpartition\tseed\tuuid8\tjobid\n" > "$MANIFEST"
SEED=10925

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

# 3090 pool (eu-g4, exclude flaky 012). Generous walltime: 3090 is slow + eu-g4 has idle nodes.
sub "eu-g4-[001-011,013-031]" gpumem:24g gpuhe.24h  16:00:00   180 3090
sub "eu-g4-[001-011,013-031]" gpumem:24g gpuhe.120h 1-00:30:00 96  3090

echo "=== DONE: 3090 = 276 (24h=180 / 120h=96); seeds 10925-$((SEED-1)). Manifest: $MANIFEST ==="
