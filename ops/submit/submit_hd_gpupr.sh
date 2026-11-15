#!/usr/bin/env bash
# HD de-novo midrun_v3 — gpupr (a100) extra capacity (user: gpupr cleared, ~50 free a100 GPUs).
# Separate QOS (es_wicky/gpupr/*) -> net capacity beyond the gpuhe caps. Seeds 11201-11260.
# a100 slow (~3h49) -> generous walltimes. Config at 28G. NO auto-refold. Appends manifest.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB/code"
MANIFEST="$RB/hd_midrun_v3_manifest.tsv"
SEED=11201

sub () {  # nodelist gres partition wall count pool
  local nodelist="$1" gres="$2" part="$3" wall="$4" n="$5" pool="$6" k out jobid uuid
  for ((k=0; k<n; k++)); do
    out=$(python3 scripts/launch.py configs/halludesign/random_init.yaml \
      --label midrun_v3 --num-runs 1 \
      --override slurm.gpus=1 --override slurm.gres="$gres" \
      --override slurm.partition="$part" --override slurm.time="$wall" \
      --override slurm.nodelist="$nodelist" \
      --override generator_params.params.seed="$SEED" 2>&1)
    jobid=$(printf '%s' "$out" | grep -oE "jobid=[0-9]+" | head -1 | cut -d= -f2)
    uuid=$(printf '%s'  "$out" | grep -oE "uuid8=[0-9a-f]+" | head -1 | cut -d= -f2)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$(date +%FT%T)" "$pool" "$part" "$SEED" "${uuid:-NA}" "${jobid:-FAIL}" >> "$MANIFEST"
    [ -z "$jobid" ] && echo "[WARN seed=$SEED $pool $part] FAILED: $(printf '%s' "$out" | tail -1)"
    SEED=$((SEED+1))
  done
  echo ">>> done $pool $part : $n (seed now $SEED)"
}

sub "eu-a65-[02-07],eu-ts-[01-02]" gpumem:40g gpupr.24h  08:00:00   40 a100
sub "eu-a65-[02-07],eu-ts-[01-02]" gpumem:40g gpupr.120h 1-00:30:00 20 a100

echo "=== DONE gpupr: 60 (24h=40/120h=20); seeds 11201-$((SEED-1)). ==="
