#!/bin/bash
# hd_purge_dryrun.sh — DRY RUN. Classify HalluDesign flat per-uuid staging dirs and report
# how many inodes each class would reclaim. Deletes NOTHING.
#
# HD writes one flat dir per uuid instead of pooling like the other stacks. When the
# canonicalize->pack tail completes, the dir collapses to ~9 inodes (designs.tar.zst +
# .list + config/logs). When it does NOT, a loose designs/ of ~3,200 files is left behind.
# Those are the inode hog.
#
# The loose dirs are NOT uniformly safe to delete: the designs are only redundant if that
# uuid also landed in the SHARED pool <label>/designs/<uuid>.tar.zst. So each dir is put in
# one of:
#   PACKED        has designs.tar.zst, no loose designs/      -> nothing to do
#   LOOSE_POOLED  loose designs/, uuid IS in the shared pool  -> loose copy is redundant, safe to drop
#   LOOSE_ORPHAN  loose designs/, uuid NOT in the pool        -> MUST pack, never delete
#   OTHER         anything else                               -> inspect by hand
#
# Usage: hd_purge_dryrun.sh [label ...]      (default: midrun_v6 midrun_v7 midrun_v8)
set -uo pipefail
HD=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn
LABELS=("$@"); [ ${#LABELS[@]} -eq 0 ] && LABELS=(midrun_v6 midrun_v7 midrun_v8)
OUT=/cluster/work/wicky/lvalsecchi/_inode_census
mkdir -p "$OUT"

for L in "${LABELS[@]}"; do
  POOL="$HD/$L/designs"
  echo "=================================================================="
  echo "$L   shared pool: $(ls "$POOL"/*.tar.zst 2>/dev/null | wc -l) slices"
  # one filesystem walk for all per-uuid dirs of this label
  du --inodes -d 1 "$HD"/${L}_* 2>/dev/null | awk '$2 ~ /'"${L}"'_[0-9a-f]+$/ {print $1, $2}' \
    > "$OUT/${L}_counts.txt"
  n=$(wc -l < "$OUT/${L}_counts.txt")
  echo "  per-uuid dirs: $n"

  : > "$OUT/${L}_classified.tsv"
  while read -r cnt d; do
    u=$(basename "$d" | sed "s/^${L}_//")
    has_tar=0; [ -s "$d/designs.tar.zst" ] && has_tar=1
    has_loose=0; [ -d "$d/designs" ] && has_loose=1
    in_pool=0; [ -s "$POOL/${u}.tar.zst" ] && in_pool=1
    if   [ $has_loose -eq 0 ] && [ $has_tar -eq 1 ]; then cls=PACKED
    elif [ $has_loose -eq 1 ] && [ $in_pool -eq 1 ]; then cls=LOOSE_POOLED
    elif [ $has_loose -eq 1 ] && [ $in_pool -eq 0 ]; then cls=LOOSE_ORPHAN
    else cls=OTHER; fi
    printf '%s\t%s\t%s\t%s\n' "$cls" "$cnt" "$u" "$d" >> "$OUT/${L}_classified.tsv"
  done < "$OUT/${L}_counts.txt"

  awk -F'\t' '{c[$1]++; s[$1]+=$2}
    END{printf "  %-13s %7s %12s\n","CLASS","DIRS","INODES";
        for(k in c) printf "  %-13s %7d %12d\n", k, c[k], s[k]}' \
    "$OUT/${L}_classified.tsv" | sort -k3 -rn
  # reclaimable = loose inodes above the ~9 a packed dir keeps, POOLED only
  awk -F'\t' '$1=="LOOSE_POOLED"{n++; s+=$2-9}
    END{printf "  >> SAFE TO RECLAIM (LOOSE_POOLED, keeping metadata): %d dirs, %d inodes\n", n+0, s+0}' \
    "$OUT/${L}_classified.tsv"
  awk -F'\t' '$1=="LOOSE_ORPHAN"{n++; s+=$2}
    END{printf "  >> MUST PACK FIRST (LOOSE_ORPHAN, do NOT delete):    %d dirs, %d inodes\n", n+0, s+0}' \
    "$OUT/${L}_classified.tsv"
done
echo "=================================================================="
echo "DRY RUN — nothing deleted. Per-dir classification in $OUT/<label>_classified.tsv"
