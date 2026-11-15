#!/bin/bash
# hd_cleanup.sh — reclaim inodes from HalluDesign per-uuid staging dirs.
#
# HD writes one flat dir per uuid, <label>_<uuid8>/, whose designs/ subtree is the RAW
# working tree: designs/recycle_N/{*.pdb, *protenix_eval/, *_mpnn_eval/}. Thousands of
# inodes per uuid. The canonicalize->pack tail is supposed to copy the finished designs
# into the SHARED pool <label>/designs/<uuid8>.tar.zst and delete the working tree; when
# the job times out mid-run that never happens and the tree is left behind.
#
# Two disjoint actions, decided per uuid dir:
#
#   PURGE   the uuid IS in the shared pool (slice >=200 B and a non-empty .list)
#           -> the working tree is redundant intermediates. Delete designs/*/ ,
#              KEEP designs/processing_results.csv (the per-recycle metrics table)
#              and everything at the dir root (config, run_meta, slurm logs).
#
#   PACK    the uuid is NOT in the pool -> the working tree is the ONLY copy of those
#           designs. tar+zstd designs/ -> designs.tar.zst (+ .list sidecar), verified
#           before the originals are removed. Lossless; unpack with untarz.
#
# Guards: skips any dir touched in the last $FRESH_HOURS hours (a running HD job owns
# it), refuses paths outside campaigns/, and never deletes a pool slice.
#
# Usage:  hd_cleanup.sh [--apply] [--only-pack|--only-purge] [--fresh-hours N] LABEL...
#         default is a DRY RUN that only prints what it would do.
set -uo pipefail

HD=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn
OUT=/cluster/work/wicky/lvalsecchi/_inode_census
APPLY=0; DO_PACK=1; DO_PURGE=1; FRESH_HOURS=6; ZL=${ZL:-3}
SHARD=0; NSHARD=1   # --shard I/N : handle only the dirs whose position % N == I
LABELS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1 ;;
    --only-pack) DO_PURGE=0 ;;
    --only-purge) DO_PACK=0 ;;
    --fresh-hours) FRESH_HOURS="$2"; shift ;;
    --shard) SHARD="${2%%/*}"; NSHARD="${2##*/}"; shift ;;
    *) LABELS+=("$1") ;;
  esac; shift
done
[ ${#LABELS[@]} -eq 0 ] && { echo "usage: $0 [--apply] LABEL..." >&2; exit 2; }
mkdir -p "$OUT"
TAG=$(date +%Y%m%d_%H%M%S)

pooled_ok() {  # <label> <uuid8>
  local s="$HD/$1/designs/$2.tar.zst"
  [ -s "$s" ] || return 1
  [ "$(stat -c %s "$s")" -ge 200 ] || return 1
  [ -s "$s.list" ] || return 1
  return 0
}

fresh() {      # <dir> : true if anything at depth<=2 was modified recently
  [ -n "$(find "$1" -maxdepth 2 -newermt "-${FRESH_HOURS} hours" -print -quit 2>/dev/null)" ]
}

pack_dir() {   # <uuidroot> : designs/ -> designs.tar.zst, verified, then rm -rf designs/
  local root="$1" d="$1/designs"
  local arch="$root/designs.tar.zst" tmp="$root/designs.tar.zst.tmp"
  local list="$root/designs.tar.zst.list" lock="$root/designs.lock"
  (
    flock -n -x 200 || { echo "PACK-SKIP  locked        $root"; exit 3; }
    [ -d "$d" ] || exit 0
    [ -f "$arch" ] && { echo "PACK-SKIP  already packed $root"; exit 0; }
    ( cd "$root" && find designs ! -type d | sort ) > "$list.tmp" || { rm -f "$list.tmp"; exit 1; }
    if [ ! -s "$list.tmp" ]; then   # designs/ holds only empty dirs — nothing to archive
      rm -f "$list.tmp"; rm -rf "$d"; echo "EMPTIED    (no files)     $root"; exit 0
    fi
    rm -f "$tmp"
    ( cd "$root" && tar -c --files-from "$list.tmp" | zstd -q -${ZL} -T2 -o "$tmp" ) \
      || { echo "PACK-FAIL  tar/zstd      $root" >&2; rm -f "$tmp" "$list.tmp"; exit 1; }
    zstd -tq "$tmp" || { echo "PACK-FAIL  verify        $root" >&2; rm -f "$tmp" "$list.tmp"; exit 1; }
    # member count must match the listing, else something changed under us
    local want have
    want=$(wc -l < "$list.tmp")
    have=$(zstd -dc "$tmp" | tar -t | grep -vc '/$')
    [ "$want" -eq "$have" ] || { echo "PACK-FAIL  count $want!=$have $root" >&2; rm -f "$tmp" "$list.tmp"; exit 1; }
    mv -f "$tmp" "$arch" && mv -f "$list.tmp" "$list" || exit 1
    rm -rf "$d"
    echo "PACKED     $want files       $root"
  ) 200>"$lock"
}

for L in "${LABELS[@]}"; do
  MAN="$OUT/${L}_cleanup_${TAG}_s${SHARD}of${NSHARD}.tsv"; : > "$MAN"
  np=0; nu=0; ip=0; iu=0; nskip=0; nfresh=0; idx=0
  for root in "$HD/${L}"_*; do
    [ -d "$root" ] || continue
    case "$root" in *"/campaigns/"*) ;; *) echo "refusing non-campaign $root" >&2; continue ;; esac
    u=$(basename "$root"); u=${u#"${L}_"}
    [[ "$u" =~ ^[0-9a-f]{6,}$ ]] || continue
    [ -d "$root/designs" ] || { nskip=$((nskip+1)); continue; }
    idx=$((idx + 1))                                   # position among candidates
    [ $(( (idx - 1) % NSHARD )) -eq "$SHARD" ] || continue
    if fresh "$root"; then nfresh=$((nfresh+1)); continue; fi
    # the inode count is a second full walk of the tree — only worth it for the dry
    # run's manifest; when applying, pack_dir's own listing reports the real number.
    if [ $APPLY -eq 1 ]; then n=0; else n=$(find "$root/designs" | wc -l); fi
    if pooled_ok "$L" "$u"; then
      [ $DO_PURGE -eq 1 ] || continue
      nu=$((nu+1)); iu=$((iu+n))
      printf 'PURGE\t%s\t%s\t%s\n' "$n" "$u" "$root" >> "$MAN"
      if [ $APPLY -eq 1 ]; then
        find "$root/designs" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} + \
          && echo "PURGED     $n inodes      $root"
      fi
    else
      [ $DO_PACK -eq 1 ] || continue
      np=$((np+1)); ip=$((ip+n))
      printf 'PACK\t%s\t%s\t%s\n' "$n" "$u" "$root" >> "$MAN"
      [ $APPLY -eq 1 ] && pack_dir "$root"
    fi
  done
  printf '%s  PURGE %d dirs / %d inodes   PACK %d dirs / %d inodes   (skipped %d no-designs, %d fresh<%dh)  manifest %s\n' \
    "$L" "$nu" "$iu" "$np" "$ip" "$nskip" "$nfresh" "$FRESH_HOURS" "$MAN"
done
[ $APPLY -eq 1 ] || echo "DRY RUN — nothing changed. Re-run with --apply."
