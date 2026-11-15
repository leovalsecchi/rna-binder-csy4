#!/bin/bash
# Post-hoc AF3 canonicalize for leaves whose outputs.tar.zst is populated but
# designs.tar.zst is the empty placeholder. Unpacks outputs to node-local /tmp,
# flattens to canonical designs/ (Stage-3 of refold_af3.sh), packs designs.tar.zst
# into the leaf (inode-neutral: replaces the empty placeholder). Optionally drops
# outputs.tar.zst (arg2=1) to reclaim inodes. No loose files ever touch shared FS.
# 2026-06-19 FIX: iterate ALL inner dirs per design (was: first only via -print -quit),
# so designs with a timestamped AND a non-timestamped inner dir keep every sample.
set -euo pipefail
LEAF="${1%/}"; DROP_OUTPUTS="${2:-0}"
source /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/archive_pipeline.sh
export PATH="/cluster/work/wicky/lvalsecchi/rna_binder/code/bin:$PATH"
[ -f "$LEAF/outputs.tar.zst" ] || { echo "SKIP(no-outputs): $LEAF"; exit 0; }
[ "$(stat -c %s "$LEAF/outputs.tar.zst")" -gt 1024 ] || { echo "SKIP(empty-outputs): $LEAF"; exit 0; }
SCRATCH=$(mktemp -d /tmp/canon_af3.XXXXXX); trap "rm -rf \"$SCRATCH\"" EXIT
unpack_unit "$LEAF/outputs.tar.zst" "$SCRATCH" >/dev/null
OUTPUTS_DIR="$SCRATCH/outputs"; CANONICAL_DIR="$SCRATCH/designs"; mkdir -p "$CANONICAL_DIR"; cd "$CANONICAL_DIR"
n=0; shopt -s nullglob
for design_dir in "$OUTPUTS_DIR"/*; do
  [ -d "$design_dir" ] || continue
  design_name=$(basename "$design_dir")
  for inner_dir in "$design_dir"/*/; do                 # ALL inner dirs (timestamped + non-timestamped)
    inner_dir="${inner_dir%/}"
    [ -d "$inner_dir" ] || continue
    for sample_dir in "$inner_dir"/seed-*_sample-*; do
      [ -d "$sample_dir" ] || continue
      sub=$(basename "$sample_dir")
      seed=$(echo "$sub" | sed -n "s/^seed-\([0-9]*\)_sample-.*/\1/p")
      sample=$(echo "$sub" | sed -n "s/^seed-[0-9]*_sample-\([0-9]*\)/\1/p")
      [ -f "$sample_dir/model.cif" ] || continue
      cp "$sample_dir/model.cif" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}.cif"
      [ -f "$sample_dir/summary_confidences.json" ] && cp "$sample_dir/summary_confidences.json" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}_summary_confidences.json"
      n=$((n+1))
    done
    inner_name=$(basename "$inner_dir")
    [ -f "$inner_dir/${inner_name}_model.cif" ] && cp "$inner_dir/${inner_name}_model.cif" "$CANONICAL_DIR/${design_name}_AF3_top.cif"
  done
done
shopt -u nullglob
[ "$n" -gt 0 ] || { echo "WARN(0-samples): $LEAF"; exit 0; }
rm -f "$LEAF/designs.tar.zst" "$LEAF/designs.tar.zst.list"
pack_unit "$CANONICAL_DIR" >/dev/null
mv -f "$SCRATCH/designs.tar.zst" "$LEAF/designs.tar.zst"
mv -f "$SCRATCH/designs.tar.zst.list" "$LEAF/designs.tar.zst.list"
echo "OK($n samples): $LEAF -> designs.tar.zst ($(wc -l < "$LEAF/designs.tar.zst.list") members)"
if [ "$DROP_OUTPUTS" = "1" ]; then rm -f "$LEAF/outputs.tar.zst" "$LEAF/outputs.tar.zst.list" "$LEAF/outputs.lock"; echo "  dropped outputs.tar.zst"; fi
