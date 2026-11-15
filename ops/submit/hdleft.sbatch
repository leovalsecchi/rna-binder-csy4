#!/bin/bash
#SBATCH --job-name=hdleft
#SBATCH --partition=normal.24h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdleft_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
HD=$RB/campaigns/csy4_crrna_free/halludesign__ligandmpnn
CC=$RB/analyses/cross_campaign

echo "=== HD packed slices per LABEL (excluding per-uuid working dirs) ==="
printf "%-44s %8s %14s\n" label slices bytes
tot=0
for d in "$HD"/*/; do
  l=$(basename "$d")
  # per-uuid working dirs are <label>_<8hex>; keep only real pool labels
  case "$l" in *_[0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]) continue;; esac
  [ -d "$d/designs" ] || continue
  n=$(find "$d/designs" -maxdepth 1 -name '*.tar.zst' 2>/dev/null | wc -l)
  [ "$n" -eq 0 ] && continue
  b=$(du -sb "$d/designs" 2>/dev/null | cut -f1)
  printf "%-44s %8s %14s\n" "$l" "$n" "$b"
  tot=$((tot+n))
done
echo "TOTAL HD slices across all labels: $tot"

echo
echo "=== HD uuids in the MASTER ==="
grep -o 'csy4free_HD_LiMPNN_len80_[0-9a-f]\{8\}' "$CC/csy4_master_foot.csv" \
  | sed 's/.*_//' | sort -u > /tmp/hd_master_uuids.txt
echo "distinct HD uuids in master: $(wc -l < /tmp/hd_master_uuids.txt)"

echo
echo "=== HD uuids PACKED on disk, per label, and how many are MISSING from master ==="
printf "%-44s %8s %8s %8s\n" label packed in_master MISSING
for d in "$HD"/*/; do
  l=$(basename "$d")
  case "$l" in *_[0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]) continue;; esac
  [ -d "$d/designs" ] || continue
  find "$d/designs" -maxdepth 1 -name '*.tar.zst' -printf '%f\n' 2>/dev/null \
    | sed 's/\.tar\.zst$//' | grep -oE '^[0-9a-f]{8}' | sort -u > /tmp/hd_lbl.txt || true
  p=$(wc -l < /tmp/hd_lbl.txt)
  [ "$p" -eq 0 ] && continue
  i=$(comm -12 /tmp/hd_lbl.txt /tmp/hd_master_uuids.txt | wc -l)
  m=$((p-i))
  printf "%-44s %8s %8s %8s\n" "$l" "$p" "$i" "$m"
done
