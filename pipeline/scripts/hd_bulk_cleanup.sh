#!/bin/bash
# hd_bulk_cleanup.sh — salvage-then-purge for orphaned HalluDesign run dirs.
#   (v7-ready + ARRAY-SAFE rewrite, 2026-07-13)
#
# An "orphan" = a uuid run dir with NO designs.tar.zst (never packed) AND NO
# active job writing it. Caused mainly by bulk preemption / CUDA device-asserts
# that kill HD after it has written some (or all) of its ligand_mpnn designs but
# before/while the pack step runs. The bulk (--requeue) strategy produces these
# continuously, which is why this is meant to run on a loop (see hd_cleanup_loop.slurm).
#
# Each orphan is classified by how many *ligand_mpnn_*.pdb it holds — the exact
# file postprocess packs:
#   * >=1 design  -> SALVAGE: run postprocess_halludesign.sh, which canonicalizes
#                    pdb->cif, packs designs.tar.zst, and removes the designs/
#                    scratch. Recovers real designs AND reclaims inodes.
#   * 0 designs   -> PURGE: rm -rf designs/ (pure prediction scratch). Inodes only.
# The lightweight shell (run.sbatch / config / logs) is always kept, so any uuid
# stays trivially resubmittable.
#
# SAFETY (why this is race-free even under the array model):
#   The OLD version protected live uuids by grepping squeue for job names
#   "hd_<uuid8>". Bucketed v7 runs under the ARRAY name (e.g. "mv7_hd"), so that
#   grep matches nothing -> the old guard was blind to array jobs. This version
#   instead uses a FRESHNESS guard: a live HD run appends designs/processing_results.csv
#   EVERY cycle and the array wrapper writes slurm.{log,err} live, so any run dir
#   touched within $HD_STALE_MIN minutes is treated as active and never touched.
#   A preempted/dead job stops writing instantly, so its dir goes stale and becomes
#   reapable after the window. (The legacy hd_<uuid8> squeue check is also kept as a
#   belt-and-suspenders for the single-job model.) HD's own run.sbatch rm -rf's
#   OUTDIR at start, so even a resubmit after salvage is safe.
#
# Config (env overrides):
#   HD_LABEL      campaign label to scan          (default: midrun_v7)
#   HD_STALE_MIN  minutes a dir must be idle       (default: 30) before it is an orphan
#
# Usage:
#   hd_bulk_cleanup.sh            # DRY-RUN (default): classify + report, no writes
#   hd_bulk_cleanup.sh --apply    # salvage packable orphans, purge empty ones
set -uo pipefail

HD=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn
PP=/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/postprocess_halludesign.sh
LOG=/cluster/work/wicky/lvalsecchi/rna_binder/.cache/hd_bulk_cleanup.log
LABEL="${HD_LABEL:-midrun_v7}"
STALE_MIN="${HD_STALE_MIN:-30}"
APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1
ts() { date +"%Y-%m-%d %H:%M:%S"; }
mkdir -p "$(dirname "$LOG")"

now=$(date +%s)
stale_s=$((STALE_MIN * 60))

# legacy protected set: uuids with a running/pending SINGLE-job HD gen ("hd_<uuid8>")
ACT=$(mktemp)
squeue -u "$USER" --noheader -o "%j" 2>/dev/null \
  | grep -oE "hd_[0-9a-f]{8}" | sed "s/hd_//" | sort -u > "$ACT" || true

# newest mtime (epoch s) among a run dir's live-writer signals; 0 if none exist.
# Keyed on files a LIVE HD run touches continuously: per-cycle processing_results.csv,
# the array-wrapper slurm.{log,err}, and any freshly written design file.
freshest() {
  local d="$1" newest=0 m
  for probe in "$d/designs/processing_results.csv" "$d/slurm.err" "$d/slurm.log"; do
    [ -f "$probe" ] || continue
    m=$(stat -c %Y "$probe" 2>/dev/null) || continue
    [ -n "$m" ] && [ "$m" -gt "$newest" ] && newest=$m
  done
  # newest actual file inside designs/ (catches a run writing pdbs but not yet the csv)
  m=$(find "$d/designs" -type f -printf '%T@\n' 2>/dev/null | sort -rn | head -1 | cut -d. -f1)
  [ -n "$m" ] && [ "$m" -gt "$newest" ] && newest=$m
  echo "$newest"
}

done=0; active=0; salv=0; salv_designs=0; purged=0; purged_ino=0; salv_fail=0
for d in "$HD"/${LABEL}_*/; do
  [ -d "$d" ] || continue
  u=$(basename "$d"); u=${u#${LABEL}_}
  [ -f "$d/designs.tar.zst" ] && { done=$((done+1)); continue; }
  grep -qx "$u" "$ACT" && { active=$((active+1)); continue; }
  [ -d "$d/designs" ] || continue
  # FRESHNESS guard: recently-written => a job is writing it now => never touch.
  fresh=$(freshest "$d")
  if [ "$fresh" -gt 0 ] && [ $((now - fresh)) -lt "$stale_s" ]; then
    active=$((active+1)); continue
  fi
  # ORPHAN (stale, unpacked, jobless)
  n=$(find "$d/designs" -name "*ligand_mpnn_*.pdb" -type f 2>/dev/null | wc -l)
  if [ "$n" -ge 1 ]; then
    if [ "$APPLY" = "1" ]; then
      if timeout 900 bash "$PP" "$d" "" >/dev/null 2>&1 && [ -f "$d/designs.tar.zst" ]; then
        got=$(wc -l < "$d/designs.tar.zst.list" 2>/dev/null || echo "$n")
        salv=$((salv+1)); salv_designs=$((salv_designs+got))
        echo "$(ts) SALVAGED $u ($got designs)" >> "$LOG"
      else
        salv_fail=$((salv_fail+1))
        echo "$(ts) SALVAGE_FAILED $u (mpnn=$n, left intact)" >> "$LOG"
      fi
    else
      salv=$((salv+1)); salv_designs=$((salv_designs+n))
    fi
  else
    ino=$(find "$d/designs" 2>/dev/null | wc -l)
    purged_ino=$((purged_ino+ino))
    if [ "$APPLY" = "1" ]; then
      rm -rf "$d/designs" && purged=$((purged+1))
      echo "$(ts) PURGED $u ($ino inodes, 0 designs)" >> "$LOG"
    else
      purged=$((purged+1))
    fi
  fi
done
rm -f "$ACT"

mode=$([ "$APPLY" = "1" ] && echo APPLY || echo DRY-RUN)
echo "[$mode label=$LABEL stale=${STALE_MIN}m] done=$done active=$active(protected) | SALVAGE: dirs=$salv designs=$salv_designs failed=$salv_fail | PURGE: dirs=$purged inodes=$purged_ino"
[ "$APPLY" = "1" ] && echo "$(ts) SUMMARY $mode label=$LABEL salvaged=$salv designs=$salv_designs purged=$purged purge_inodes=$purged_ino salv_fail=$salv_fail" >> "$LOG"
exit 0
