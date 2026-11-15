#!/bin/bash
# Progress-priority for BULK refolds (af3+b2): nice low-progress DOWN so near-done jump ahead
# (Euler: scontrol top disabled, users can't raise priority). Scale-safe + reliable:
# ONE `scontrol show job -o` (comma/space lists are "Invalid"; %o is null), filter my pending
# bulk refolds, map leaf via StdOut, read progress from outputs/ or the packed outputs.tar.zst.list.
# Usage: refold_priority_babysitter.sh [--go]   (dry-run default)
#
# NICE_MAX MUST stay below the bulk base priority, else SLURM clamps every job's
# (base - nice) to the floor of 1, they all tie, and the scheduler runs them oldest-first
# (ignoring progress). The base is fairshare-driven (~1981 on 2026-06-16) and drifts over time,
# so each run we AUTO-PROBE it: set one pending bulk job to nice=0, read its Priority (= base),
# and use 90% of that as NICE_MAX. Pin NICE_MAX via env to skip the probe.
set -uo pipefail
NICE_MAX="${NICE_MAX:-20}"  # [nice-cap] pinned small (was base*0.9~6183 -> ranked us below other users)
GO=0; [ "${1:-}" = "--go" ] && GO=1
[ "$GO" = 1 ] && echo "[applying nice]" || echo "[dry-run]"

# list my pending bulk refolds once (reused for probe + main loop)
JOBS=$(scontrol show job -o 2>/dev/null | grep 'UserId=lvalsecchi' | grep -E 'JobName=refold_(af3|b2|esm)' | grep 'JobState=PENDING' | grep -iE 'Partition=gpuhe.bulk')

# auto-probe the bulk base priority unless NICE_MAX is pinned via env
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
    NICE_MAX=1800   # fallback: probe unavailable/clamped (dry-run, no pending jobs, fresh fairshare)
  fi
fi
echo "[base-probe job=${probe:-none} base=${base:-?} -> NICE_MAX=$NICE_MAX]"

cnt=0
while read -r line; do
  [ -z "$line" ] && continue
  jid=$(echo "$line" | grep -oE 'JobId=[0-9]+' | head -1 | cut -d= -f2)
  leaf=$(echo "$line" | grep -oE 'StdOut=\S+' | head -1 | sed 's#StdOut=##; s#/[^/]*$##')
  { [ -z "$jid" ] || [ -z "$leaf" ]; } && continue
  exp=$(grep -oE '[0-9]+' "$leaf/expected.json" 2>/dev/null | head -1); exp=${exp:-600}
  done=$(find "$leaf/outputs" -name '*.cif' 2>/dev/null | wc -l)
  if [ "$done" -eq 0 ] && [ -f "$leaf/outputs.tar.zst.list" ]; then done=$(grep -c '\.cif$' "$leaf/outputs.tar.zst.list" 2>/dev/null); fi
  # ESMFold2 writes canonical CIFs straight into designs/ (no outputs/ stage).
  if [ "$done" -eq 0 ]; then done=$(find "$leaf/designs" -name '*.cif' 2>/dev/null | wc -l); fi
  if [ "$done" -eq 0 ] && [ -f "$leaf/designs.tar.zst.list" ]; then done=$(grep -c '\.cif$' "$leaf/designs.tar.zst.list" 2>/dev/null); fi
  [ "$done" -gt "$exp" ] && done=$exp
  nice=$(( (exp-done)*NICE_MAX/exp ))
  [ "$GO" = 1 ] && scontrol update jobid="$jid" nice="$nice" >/dev/null 2>&1
  cnt=$((cnt+1))
done < <(printf '%s\n' "$JOBS")
echo "niced $cnt pending bulk refolds by progress (NICE_MAX=$NICE_MAX)"
