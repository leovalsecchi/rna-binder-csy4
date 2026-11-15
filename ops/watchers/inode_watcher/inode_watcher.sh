#!/bin/bash
# inode_watcher.sh — ONE tick of the project-inode brake.
#
# Reads the wicky PROJECT inode usage (lfs quota -p 510551) and, with hysteresis,
# HOLDS all schedulable PENDING jobs when inodes get high, RELEASES them once
# inodes fall back. Running jobs are NEVER touched (they pack down -> reclaim inodes).
#
#   DRY-RUN (default): prints what it would do, changes nothing.
#   --apply          : actually scontrol hold / release.
#
# Safety:
#   * only holds jobs that are PENDING *and schedulable* (reason !~ /Held/), so it
#     never re-touches a job you (or an admin) already held manually.
#   * on RELEASE it frees BOTH the job-IDs it recorded in _held_jobids AND any
#     orphaned JobHeldUser managed jobs (a prior run may have lost its record),
#     so inode-holds self-heal once inodes recover. Infra loops (EXCLUDE_RE) are
#     never released.
#   * never holds the infra loops/babysitters (EXCLUDE_RE).
#   * fail-safe: if the quota read fails, it does NOTHING.
set -uo pipefail

IW="${IW:-/cluster/work/wicky/lvalsecchi/rna_binder/ops/watchers/inode_watcher}"
PROJ="${INODE_PROJ:-510551}"
FS="${INODE_FS:-/cluster/work/wicky}"
HIGH="${INODE_HIGH:-7400000}"     # brake just under the 7.5M soft quota
LOW="${INODE_LOW:-7200000}"       # release once we fall back below this
UNAME="${INODE_USER:-lvalsecchi}"

HELD_FILE="$IW/_held_jobids"
STATE_FILE="$IW/_state"           # "held" | "normal"
ENABLED_FILE="$IW/enabled"        # "1" = active, anything else = paused
LOG="$IW/watcher.log"

EXCLUDE_RE='gpu_watcher_loop|inode_watcher_loop|gpubulk_prio_loop|refold_prio_loop|mv4_refold_loop|bif_fill_refold_loop|hd_refold_bulk_loop'

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

ts(){ date '+%F %T'; }
log(){ echo "$(ts) $*" | tee -a "$LOG"; }

mkdir -p "$IW"

# --- enabled toggle ---
en=$(cat "$ENABLED_FILE" 2>/dev/null || echo 1)
if [ "$en" != "1" ]; then log "disabled (enabled=$en) -> no action"; exit 0; fi

# --- read current project inode count (robust to lfs line-wrapping) ---
# numeric tokens after the header, in order: kbytes bquota blimit files iquota ilimit
# (bgrace/igrace are '-' and dropped; this survives the path-on-its-own-line wrap)
nums=$(lfs quota -p "$PROJ" "$FS" 2>/dev/null | tail -n +3 | tr -s ' ' '\n' | grep -E '^[0-9]+$')
files=$(echo "$nums" | sed -n '4p')
iquota=$(echo "$nums" | sed -n '5p')
ilimit=$(echo "$nums" | sed -n '6p')
if ! [[ "$files" =~ ^[0-9]+$ ]]; then
  log "WARN could not read inode count (got '$files') -> fail-safe, no action"
  exit 0
fi

state=$(cat "$STATE_FILE" 2>/dev/null || echo normal)
mode=$([ "$APPLY" = 1 ] && echo APPLY || echo DRY-RUN)

# list of schedulable pending non-infra job ids (id name reason)
schedulable_ids(){
  squeue -u "$UNAME" -h -t PD -o "%i %j %r" \
    | awk '$3 !~ /Held/' \
    | grep -vE "$EXCLUDE_RE" \
    | awk '{print $1}'
}

if [ "$files" -ge "$HIGH" ]; then
  # ---- HOLD phase (covers initial trip AND newcomers while already held) ----
  ids=$(schedulable_ids)
  n=$(echo -n "$ids" | grep -c . || true)
  if [ "$n" -eq 0 ]; then
    log "[$mode] files=$files >= HIGH=$HIGH ; nothing schedulable to hold (state=$state)"
  else
    log "[$mode] files=$files >= HIGH=$HIGH -> HOLD $n schedulable pending job(s)"
    if [ "$APPLY" = 1 ]; then
      scontrol hold $(echo "$ids" | tr '\n' ' ') 2>>"$LOG" || true
      { cat "$HELD_FILE" 2>/dev/null; echo "$ids"; } | grep -E '^[0-9]' | sort -un > "$HELD_FILE.tmp" && mv "$HELD_FILE.tmp" "$HELD_FILE"
      echo held > "$STATE_FILE"
    fi
  fi

elif [ "$files" -le "$LOW" ]; then
  # ---- RELEASE phase (self-healing) ----
  # Free the IDs we recorded PLUS any orphaned JobHeldUser managed jobs: a
  # prior run may have held jobs whose _held_jobids record was later reset,
  # leaving them stuck (held forever) even after inodes recovered. Infra loops
  # (EXCLUDE_RE) are never touched; releasing a non-held id is a harmless no-op.
  recorded=$(cat "$HELD_FILE" 2>/dev/null | grep -E '^[0-9]' || true)
  orphaned=$(squeue -u "$UNAME" -h -t PD -o "%i %j %r" | awk '$3 ~ /Held/' | grep -vE "$EXCLUDE_RE" | awk '{print $1}')
  rel=$(printf '%s\n%s\n' "$recorded" "$orphaned" | grep -E '^[0-9]' | sort -u)
  n=$(echo -n "$rel" | grep -c . || true)
  if [ "$n" -gt 0 ]; then
    log "[$mode] files=$files <= LOW=$LOW -> RELEASE $n job(s) (recorded+orphaned)"
    if [ "$APPLY" = 1 ]; then
      echo "$rel" | xargs -r -n 200 scontrol release 2>>"$LOG" || true
      : > "$HELD_FILE"
      echo normal > "$STATE_FILE"
    fi
  else
    [ "$APPLY" = 1 ] && echo normal > "$STATE_FILE"
    log "[$mode] files=$files <= LOW=$LOW -> nothing held to release (state=$state)"
  fi

else
  log "[$mode] files=$files (HIGH=$HIGH LOW=$LOW) state=$state -> steady, no change"
fi
