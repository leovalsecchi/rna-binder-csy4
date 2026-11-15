#!/bin/bash
# hd_combine.sh — consolidate FINISHED HalluDesign per-uuid leaves into the shared
# <run>/designs/<uuid8>.tar.zst pool (flat .cif — the format the refolders expect),
# so hd_refold_bulk.py can pick them up. Two cases:
#
#   CLEAN  (per-uuid <run>_<uuid>/designs.tar.zst exists; recycle_* already
#           canonicalized+purged by the gen finalize): unpack -> mv the canonical
#           .cif to the shared pool -> pack_slice (flat per-uuid slice). Leaves the
#           per-uuid designs.tar.zst in place as a backup.
#   LOOSE  (timed-out: <run>_<uuid>/designs/recycle_* present, never finalized):
#           run postprocess_halludesign.sh RUN_DIR SHARED  — flattens the recycle
#           trees, canonicalizes pdb->cif, pack_slices into the shared pool, then
#           removes the per-uuid staging dir. This RECOVERS the partial (completed)
#           recycles of a timed-out job — HD has no resume, but the finished recycles
#           are valid designs and this packs exactly those.
#
# Skips uuids still generating (live hd_<uuid8> job) and uuids already in the shared
# pool. Idempotent. DRY-RUN unless --go.
#   Usage: hd_combine.sh [RUN] [--go] [--limit N]   (RUN default midrun_v3)
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
HDBASE=$RB/campaigns/csy4_crrna_free/halludesign__ligandmpnn
PP=$RB/code/scripts/stages/postprocess_halludesign.sh

RUN=midrun_v3; GO=0; LIMIT=0; FULL_MIN=150   # timed-out leaf must have >=FULL_MIN designs (all 15 recycles)
while [ $# -gt 0 ]; do
  case "$1" in
    --go) GO=1;;
    --limit) LIMIT="${2:-0}"; shift;;
    midrun_v*) RUN="$1";;
    *) echo "unknown arg: $1" >&2;;
  esac
  shift
done

SHARED="$HDBASE/$RUN/designs"
# shellcheck disable=SC1090
source "$RB/code/scripts/stages/archive_pipeline.sh" 2>/dev/null || true
declare -F pack_slice >/dev/null || { echo "FATAL: pack_slice not available" >&2; exit 1; }

squeue -u lvalsecchi -h -o '%j' 2>/dev/null | grep -oE 'hd_[0-9a-f]{8}' | sed 's/hd_//' | sort -u > /tmp/hd_active.txt
mkdir -p "$SHARED"

clean=0; loose_full=0; partial=0; empty=0; already=0; ok_clean=0; ok_loose=0
for d in "$HDBASE/${RUN}_"*/; do
  [ -d "$d" ] || continue
  u=$(basename "$d" | sed "s/^${RUN}_//")
  [[ "$u" =~ ^[0-9a-f]{8}$ ]] || continue
  grep -qx "$u" /tmp/hd_active.txt && continue              # still generating
  [ -f "$SHARED/$u.tar.zst" ] && { already=$((already+1)); continue; }

  if [ -f "$d/designs.tar.zst" ]; then
    clean=$((clean+1))
    [ "$GO" = 1 ] || continue
    sc=$(mktemp -d)
    if zstd -dq -c "$d/designs.tar.zst" | ( cd "$sc" && tar -x ); then
      mv "$sc"/designs/*_"${u}"_* "$SHARED"/ 2>/dev/null || true
      if pack_slice "$SHARED" "$u" >/dev/null 2>&1 && [ -f "$SHARED/$u.tar.zst" ]; then
        ok_clean=$((ok_clean+1))
      fi
    fi
    rm -rf "$sc"
  elif ls -d "$d"/designs/recycle_* >/dev/null 2>&1; then
    # timed-out leaf: recover ONLY if FULL (all 15 recycles -> >=150 designs).
    # Partial leaves only hold early/un-converged cycles (plddt ~60) and are
    # missing their best (recycle ~14-15) designs, so we skip them.
    dc=$(find "$d/designs" -name '*ligand_mpnn_*.pdb' -type f 2>/dev/null | wc -l)
    if [ "$dc" -lt "$FULL_MIN" ]; then partial=$((partial+1)); continue; fi
    loose_full=$((loose_full+1))
    [ "$GO" = 1 ] || continue
    bash "$PP" "$d" "$SHARED" >/dev/null 2>&1 || true
    [ -f "$SHARED/$u.tar.zst" ] && ok_loose=$((ok_loose+1))
  else
    empty=$((empty+1))
  fi
  [ "$LIMIT" -gt 0 ] && [ $((ok_clean+ok_loose)) -ge "$LIMIT" ] && break
done
echo "[hd_combine $RUN GO=$GO] candidates: clean=$clean loose_FULL=$loose_full loose_partial(skip)=$partial empty=$empty already=$already | combined: clean=$ok_clean loose_full=$ok_loose"
