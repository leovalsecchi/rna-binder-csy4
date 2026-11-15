#!/usr/bin/env bash
# Distribute pending refold jobs across walltime tiers, respecting GPU fit + walltime.
# B2 -> gpu.* (2080 Ti);  ESM/AF3 -> gpuhe.* (Ampere 24G);  AF3 never in 4h (too slow).
set -uo pipefail
U=lvalsecchi

move_chunks() {
  # args: jobname  then repeating: partition count timelimit(- = keep)
  local name=$1; shift
  mapfile -t ids < <(squeue -u "$U" -h -t PENDING -n "$name" -o "%i" 2>/dev/null)
  local total=${#ids[@]} i=0
  echo "### $name: $total pending"
  while [ "$#" -ge 3 ]; do
    local part=$1 cnt=$2 tl=$3; shift 3
    [ "$cnt" = "rest" ] && cnt=$((total - i))
    local end=$((i + cnt)); [ "$end" -gt "$total" ] && end=$total
    local chunk="" k=$i
    while [ "$k" -lt "$end" ]; do chunk="$chunk,${ids[$k]}"; k=$((k+1)); done
    chunk=${chunk#,}
    local moved=$((end - i))
    if [ -n "$chunk" ]; then
      if [ "$tl" = "-" ]; then
        scontrol update jobid="$chunk" partition="$part" 2>/dev/null \
          && echo "  -> $part : $moved" || echo "  -> $part : FAILED ($moved)"
      else
        scontrol update jobid="$chunk" partition="$part" TimeLimit="$tl" 2>/dev/null \
          && echo "  -> $part (t=$tl) : $moved" || echo "  -> $part (t=$tl) : FAILED ($moved)"
      fi
    fi
    i=$end
  done
}

move_chunks refold_b2  gpu.4h   100 04:00:00   gpu.24h   70  -   gpu.bulk   rest -
move_chunks refold_esm gpuhe.4h 100 04:00:00   gpuhe.24h 100 -   gpuhe.bulk rest -
move_chunks refold_af3 gpuhe.24h 80 -          gpuhe.120h 100 -  gpuhe.bulk rest -

echo
echo "=== RESULT: my refold jobs by partition ==="
squeue -u "$U" -h -n refold_b2,refold_esm,refold_af3 -o "%P %T" 2>/dev/null | sort | uniq -c | sort -k2
