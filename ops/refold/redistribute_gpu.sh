#!/usr/bin/env bash
# Fill the idle Turing-24G gpu.* pool (quadro/titan) with AF3/ESM pulled from the
# overloaded gpuhe.*. AF3 -> Turing-24G via xla; ESM -> via PyTorch fallback. Both need
# 24G so we pin the 24G Turing nodes (NOT the 11G 2080 Ti, which stay for B2).
set -uo pipefail
U=lvalsecchi
TUR24='eu-g3-[062-077,083-084],eu-lo-g3-[049-061]'   # quadro+titan 24G

# move_to_gpu NAME DST_PART COUNT TIMELIMIT(- keep) SRC_PARTS...
move_to_gpu() {
  local name=$1 dst=$2 cnt=$3 tl=$4; shift 4
  local srcs; srcs=$(printf '%s,' "$@"); srcs=${srcs%,}
  local m=0
  for j in $(squeue -u "$U" -h -t PENDING -n "$name" -p "$srcs" -o "%i" 2>/dev/null | head -"$cnt"); do
    if [ "$tl" = "-" ]; then
      scontrol update jobid="$j" partition="$dst" nodelist="$TUR24" 2>/dev/null && m=$((m+1))
    else
      scontrol update jobid="$j" partition="$dst" nodelist="$TUR24" TimeLimit="$tl" 2>/dev/null && m=$((m+1))
    fi
  done
  echo "  $name -> $dst (Turing-24G): $m"
}

echo "=== ESM onto gpu.* Turing-24G ==="
move_to_gpu refold_esm gpu.4h   80 04:00:00 gpuhe.24h gpuhe.bulk
move_to_gpu refold_esm gpu.bulk 40 -        gpuhe.bulk gpuhe.24h
echo "=== AF3 onto gpu.* Turing-24G (no 4h) ==="
move_to_gpu refold_af3 gpu.24h  60 -        gpuhe.24h gpuhe.bulk
move_to_gpu refold_af3 gpu.120h 60 -        gpuhe.bulk gpuhe.24h

echo
echo "=== gpu.* Turing-24G reason check (fit OK?) ==="
squeue -u "$U" -h -p gpu.4h,gpu.24h,gpu.120h,gpu.bulk -o "%R" 2>/dev/null | grep -ivE "^\(None\)|Resources|Priority" | sort | uniq -c | head
echo "=== ALL refolds by partition ==="
squeue -u "$U" -h -n refold_b2,refold_esm,refold_af3 -o "%P %T" 2>/dev/null | sort | uniq -c | sort -k2
echo "=== by TIER ==="
squeue -u "$U" -h -n refold_b2,refold_esm,refold_af3 -o "%P" 2>/dev/null | sed -E 's/.*\.//' | sort | uniq -c
