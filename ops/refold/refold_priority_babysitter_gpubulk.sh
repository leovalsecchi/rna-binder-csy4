#!/bin/bash
# gpu.bulk progress-priority for B2 refolds on the PUBLIC gpu.bulk tier.
# SEPARATE from refold_priority_babysitter.sh (which manages the wicky gpuhe.bulk
# tier) because the two pools have DIFFERENT fairshare bases -> each needs its own
# NICE_MAX probe. Only B2 runs on gpu.bulk (AF3/ESM need Ampere, gpu.* is Turing).
# Nice low-progress jobs DOWN so near-done leaves jump ahead.
# Usage: refold_priority_babysitter_gpubulk.sh [--go]   (dry-run default)
set -uo pipefail
NICE_MAX="${NICE_MAX:-20}"  # [nice-cap] pinned small (was base*0.9~6183 -> ranked us below other users)
GO=0; [ "${1:-}" = "--go" ] && GO=1
[ "$GO" = 1 ] && echo "[applying nice]" || echo "[dry-run]"
# my pending B2 jobs on the PUBLIC gpu.bulk partition (exact match, not gpuhe.bulk)
JOBS=$(scontrol show job -o 2>/dev/null | grep 'UserId=lvalsecchi' | grep -E 'JobName=refold_b2' | grep 'JobState=PENDING' | grep -E 'Partition=gpu\.bulk')
# auto-probe the gpu.bulk base priority unless NICE_MAX is pinned via env
probe=""; base=""
if [ -z "${NICE_MAX:-}" ]; then
  probe=$(printf '%s\n' "$JOBS" | grep -oE 'JobId=[0-9]+' | head -1 | cut -d= -f2)
  if [ -n "$probe" ]; then
    [ "$GO" = 1 ] && scontrol update jobid="$probe" nice=0 >/dev/null 2>&1
    base=$(scontrol show job -o "$probe" 2>/dev/null | grep -oE 'Priority=[0-9]+' | head -1 | cut -d= -f2)
  fi
  if [ -n "$base" ] && [ "$base" -gt 100 ] 2>/dev/null; then
    NICE_MAX=$(( base * 9 / 10 ))
  else
    NICE_MAX=1800   # fallback: probe unavailable/clamped
  fi
fi
echo "[gpu.bulk base-probe job=${probe:-none} base=${base:-?} -> NICE_MAX=$NICE_MAX]"
cnt=0
while read -r line; do
  [ -z "$line" ] && continue
  jid=$(echo "$line" | grep -oE 'JobId=[0-9]+' | head -1 | cut -d= -f2)
  leaf=$(echo "$line" | grep -oE 'StdOut=\S+' | head -1 | sed 's#StdOut=##; s#/[^/]*$##')
  { [ -z "$jid" ] || [ -z "$leaf" ]; } && continue
  exp=$(grep -oE '[0-9]+' "$leaf/expected.json" 2>/dev/null | head -1); exp=${exp:-600}
  done=$(find "$leaf/outputs" -name '*.cif' 2>/dev/null | wc -l)
  if [ "$done" -eq 0 ] && [ -f "$leaf/outputs.tar.zst.list" ]; then done=$(grep -c '\.cif$' "$leaf/outputs.tar.zst.list" 2>/dev/null); fi
  if [ "$done" -eq 0 ]; then done=$(find "$leaf/designs" -name '*.cif' 2>/dev/null | wc -l); fi
  if [ "$done" -eq 0 ] && [ -f "$leaf/designs.tar.zst.list" ]; then done=$(grep -c '\.cif$' "$leaf/designs.tar.zst.list" 2>/dev/null); fi
  [ "$done" -gt "$exp" ] && done=$exp
  nice=$(( (exp-done)*NICE_MAX/exp ))
  [ "$GO" = 1 ] && scontrol update jobid="$jid" nice="$nice" >/dev/null 2>&1
  cnt=$((cnt+1))
done < <(printf '%s\n' "$JOBS")
echo "niced $cnt pending gpu.bulk B2 refolds by progress (NICE_MAX=$NICE_MAX)"
