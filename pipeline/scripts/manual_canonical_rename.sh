#!/bin/bash
# manual_canonical_rename.sh — extracted Stage 3 of refold_af3.sh
# Usage: ./manual_canonical_rename.sh REFOLD_DIR

set -euo pipefail

REFOLD_DIR="${1:?Usage: $0 REFOLD_DIR}"
OUTPUTS_DIR="$REFOLD_DIR/outputs"
CANONICAL_DIR="$REFOLD_DIR/designs"

mkdir -p "$CANONICAL_DIR"

echo "=== Manual canonical rename ==="
echo " Outputs dir:    $OUTPUTS_DIR"
echo " Canonical dir:  $CANONICAL_DIR"
echo ""

cd "$CANONICAL_DIR"
n_canonical=0

shopt -s nullglob
for design_dir in "$OUTPUTS_DIR"/*; do
    [ -d "$design_dir" ] || continue
    design_name=$(basename "$design_dir")
    inner_dir=$(find "$design_dir" -mindepth 1 -maxdepth 1 -type d | head -1)
    if [ -z "$inner_dir" ] || [ ! -d "$inner_dir" ]; then
        echo "  WARN: no inner dir under $design_dir, skipping"
        continue
    fi
    for sample_dir in "$inner_dir"/seed-*_sample-*; do
        [ -d "$sample_dir" ] || continue
        sub=$(basename "$sample_dir")
        seed=$(echo "$sub" | sed -n 's/^seed-\([0-9]*\)_sample-.*/\1/p')
        sample=$(echo "$sub" | sed -n 's/^seed-[0-9]*_sample-\([0-9]*\)/\1/p')
        src_cif="$sample_dir/model.cif"
        [ -f "$src_cif" ] || continue
        new_name="${design_name}_AF3_seed${seed}_sample${sample}.cif"
        if [ -f "$CANONICAL_DIR/$new_name" ]; then
            continue
        fi
        cp "$src_cif" "$CANONICAL_DIR/$new_name"
        for sidecar in summary_confidences.json confidences.json; do
            src_side="$sample_dir/$sidecar"
            if [ -f "$src_side" ]; then
                cp "$src_side" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}_${sidecar}"
            fi
        done
        n_canonical=$((n_canonical + 1))
    done
done

echo ""
echo "================================================================"
echo " Done. New canonical .cif files written: $n_canonical"
echo " Canonical dir: $CANONICAL_DIR"
echo "================================================================"
