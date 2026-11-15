#!/bin/bash
# apo_submit_v9.sh — submit the v9 apo arrays across both GPU pools.
#
# Chunk assignment is deliberate: chunks 0-322 are the STACK-BALANCED phase (200 of
# every stack per chunk). Those go to the TIMED tiers, sized to exactly each tier's
# memory-cap concurrency (mem cap / 48 G per task) so they start immediately and are
# never preempted. The surplus goes to the bulk tiers, starting with the remaining
# balanced chunks. Every chunk is internally balanced, so any completion order still
# yields equal per-stack coverage — the ordering just maximises how much lands early.
#
# Directives copied verbatim from the v8 apo jobs (9453681/9454196/9454204/9454206/
# 9454210) and the running gpuhe apo array: cpu=8, mem-per-cpu=6G (48 G/task), gpus=1.
#   gpu.*   pool -> --nodelist=eu-g3-[062-077,083-084], NO gpumem
#   gpuhe.* pool -> --gres=gpumem:24g, NO nodelist  (11 G 2080 Tis live only in gpu.*)
# Walltime ROUTES to tier and overrides --partition, so: 4h=04:00:00, 24h=06:00:00,
# 120h=1-06:00:00, bulk=06:00:00 + a scontrol move while still PENDING (a plain
# --partition=*.bulk submit silently remaps to *.24h).
# NO %-throttles: a pending array is one scheduler record, so a throttle buys zero
# queue-window protection and only caps throughput below the QOS entitlement.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
S=apo_fold_array_v9.sh
PIN="eu-g3-[062-077,083-084]"
IDS=apo_chunks_v9/.apo_v9_jobids
: > "$IDS"

sub_timed () {   # name range partition qos time extra...
  local name=$1 range=$2 part=$3 qos=$4 tl=$5; shift 5
  local j
  j=$(sbatch --parsable --array="$range" --partition="$part" --qos="$qos" \
        --time="$tl" "$@" "$S")
  echo "$name  $j  array=$range  $part/$qos  $tl"
  echo "$name $j $range $part" >> "$IDS"
}

sub_bulk () {    # name range partition qos extra...
  local name=$1 range=$2 part=$3 qos=$4; shift 4
  local j
  j=$(sbatch --parsable --hold --array="$range" --time=06:00:00 "$@" "$S")
  scontrol update jobid="$j" partition="$part" qos="$qos"
  scontrol release "$j"
  echo "$name  $j  array=$range  -> $part/$qos (moved while pending)"
  echo "$name $j $range $part" >> "$IDS"
}

echo "=== timed tiers: stack-balanced chunks 0-203 ==="
sub_timed gpu.4h     0-71    gpu.4h    es_wicky/gpu/4     04:00:00   --nodelist="$PIN"
sub_timed gpuhe.4h   72-143  gpuhe.4h  es_wicky/gpuhe/4   04:00:00   --gres=gpumem:24g
sub_timed gpu.24h    144-162 gpu.24h   es_wicky/gpu/24    06:00:00   --nodelist="$PIN"
sub_timed gpuhe.24h  163-181 gpuhe.24h es_wicky/gpuhe/24  06:00:00   --gres=gpumem:24g
sub_timed gpu.120h   182-192 gpu.120h  es_wicky/gpu/120   1-06:00:00 --nodelist="$PIN"
sub_timed gpuhe.120h 193-203 gpuhe.120h es_wicky/gpuhe/120 1-06:00:00 --gres=gpumem:24g

echo "=== bulk tiers: chunks 204-1192 (204-322 still balanced) ==="
sub_bulk  gpu.bulk   204-698  gpu.bulk   es_wicky/gpu/bulk   --nodelist="$PIN"
sub_bulk  gpuhe.bulk 699-1192 gpuhe.bulk es_wicky/gpuhe/bulk --gres=gpumem:24g

echo
echo "job ids -> $IDS"
cat "$IDS"
