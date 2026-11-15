#!/bin/bash
# postprocess_boltzdesign.sh — flatten BD1's nested results_final/ into canonical names.
#
# Usage:
#   bash postprocess_boltzdesign.sh <RUN_DIR>
#
# BD1 writes designs to:
#   outputs/rna_<target>_<suffix>/results_final/boltz_results_<...>_itr<N>_length<L>/
#       predictions/<...>_itr<N>_length<L>/<...>_itr<N>_length<L>_model_0.cif
#
# Canonical target:
#   designs/<target>_BD1_LiMPNN_len<L>_<uuid8>_design<idx>.cif
#
# Where idx is the sequential index across sorted itr numbers.
# Idempotent: re-running won't duplicate (uses cp -n logic).

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 RUN_DIR" >&2
    exit 2
fi

RUN_DIR="$1"
DESIGNS_DIR="$RUN_DIR/designs"
OUTPUTS_DIR="$RUN_DIR/outputs"

if [ ! -d "$OUTPUTS_DIR" ]; then
    echo "ERROR: $OUTPUTS_DIR does not exist (BD1 hasn't produced outputs yet)" >&2
    exit 1
fi

META="$RUN_DIR/run_meta.json"
CONFIG="$RUN_DIR/resolved_config.yaml"

if [ ! -f "$META" ] || [ ! -f "$CONFIG" ]; then
    echo "ERROR: missing run_meta.json or resolved_config.yaml in $RUN_DIR" >&2
    exit 1
fi

mkdir -p "$DESIGNS_DIR"

HEADERS=$(python3 << PYEOF
import json, re, sys
import yaml
with open("$META") as f:
    meta = json.load(f)
with open("$CONFIG") as f:
    cfg = yaml.safe_load(f)
run_basename = "$(basename "$RUN_DIR")"
uuid_match = re.search(r"_([0-9a-f]{8})\$", run_basename)
uuid8 = uuid_match.group(1) if uuid_match else meta.get("uuid8", meta.get("uuid", "unknown"))
target_full = cfg.get("target_id", "csy4_crrna_free")
target_short = "csy4free" if target_full == "csy4_crrna_free" else target_full.replace("_", "")
length = 80
gp = cfg.get("generator_params", {})
for k in ("length_min", "length", "binder_length"):
    if k in gp:
        try:
            length = int(gp[k]); break
        except Exception:
            pass
print(f"{target_short}|{uuid8}|{length}")
PYEOF
)

TARGET_SHORT=$(echo "$HEADERS" | cut -d'|' -f1)
UUID8=$(echo "$HEADERS" | cut -d'|' -f2)
LENGTH=$(echo "$HEADERS" | cut -d'|' -f3)
PREFIX="${TARGET_SHORT}_BD1_LiMPNN_len${LENGTH}_${UUID8}"

echo "================================================================"
echo " BD1 postprocess"
echo " Run dir:     $RUN_DIR"
echo " Prefix:      $PREFIX"
echo "================================================================"

RESULTS_FINAL=""
for d in "$OUTPUTS_DIR"/*/results_final; do
    if [ -d "$d" ]; then
        RESULTS_FINAL="$d"; break
    fi
done

if [ -z "$RESULTS_FINAL" ]; then
    echo "ERROR: no results_final/ dir found under $OUTPUTS_DIR" >&2
    exit 1
fi
echo "Results root: $RESULTS_FINAL"

mapfile -t ITR_DIRS < <(
    ls -d "$RESULTS_FINAL"/boltz_results_*_itr*_length* 2>/dev/null \
    | awk '{ match($0, /_itr([0-9]+)_/, m); if (m[1] != "") print m[1] "\t" $0; }' \
    | sort -n \
    | cut -f2
)

if [ ${#ITR_DIRS[@]} -eq 0 ]; then
    echo "WARNING: no boltz_results_*_itr*_length* dirs in $RESULTS_FINAL" >&2
    exit 0
fi
echo "Found ${#ITR_DIRS[@]} iteration dirs"

design_idx=0
n_copied=0
for itr_dir in "${ITR_DIRS[@]}"; do
    cif=$(find "$itr_dir" -name '*_model_0.cif' -type f | head -1)
    if [ -z "$cif" ]; then
        echo "  SKIP (no model_0.cif): $(basename "$itr_dir")" >&2
        design_idx=$((design_idx + 1)); continue
    fi
    dst="$DESIGNS_DIR/${PREFIX}_design${design_idx}.cif"
    if [ ! -f "$dst" ]; then
        cp "$cif" "$dst"
        n_copied=$((n_copied + 1))
    fi
    conf=$(find "$itr_dir" -name 'confidence_*_model_0.json' -type f | head -1)
    if [ -n "$conf" ]; then
        conf_dst="$DESIGNS_DIR/${PREFIX}_design${design_idx}_confidence.json"
        [ ! -f "$conf_dst" ] && cp "$conf" "$conf_dst"
    fi
    design_idx=$((design_idx + 1))
done

echo ""
echo "Canonical-named files added: $n_copied"
echo "Total canonical .cif in $DESIGNS_DIR:"
ls "$DESIGNS_DIR"/${PREFIX}_design*.cif 2>/dev/null | wc -l
echo "================================================================"
