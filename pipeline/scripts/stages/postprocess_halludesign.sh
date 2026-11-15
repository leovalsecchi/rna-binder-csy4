#!/bin/bash
# postprocess_halludesign.sh — flatten HD's nested recycle_N/.../ligand_mpnn_S_K.pdb
# files into canonical names at top-level designs/.
#
# Usage:
#   bash postprocess_halludesign.sh <RUN_DIR>
#
# Where RUN_DIR is the HD run dir (contains designs/ subdir with recycle_N/ inside).
#
# Source pattern:
#   designs/recycle_<N>/<scaffold>_recycle_<N>_mpnn_eval/ligand_mpnn/<scaffold>_recycle_<N>ligand_mpnn_<S>_<K>.pdb
#
# Canonical target:
#   designs/<target>_HD_LiMPNN_len<L>_<uuid8>_recycle<N>_sample<S>_pack<K>.pdb
#
# Where:
#   <target>, <uuid8>, <L> come from run_meta.json + resolved_config.yaml
#   N = cycle (1-indexed by HD)
#   S = sequence index (0-indexed: 0..num_seqs-1)
#   K = sidechain pack index (1-indexed: 1..number_of_packs_per_design, default 1)
#
# This script is idempotent: re-running won't duplicate (uses cp -n).

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 RUN_DIR" >&2
    exit 2
fi

RUN_DIR="$1"
SHARED_DESIGNS_DIR="${2:-}"   # --num-runs>1: combine target (empty = single-run)
DESIGNS_DIR="$RUN_DIR/designs"

if [ ! -d "$DESIGNS_DIR" ]; then
    echo "ERROR: $DESIGNS_DIR does not exist" >&2
    exit 1
fi

META="$RUN_DIR/run_meta.json"
CONFIG="$RUN_DIR/resolved_config.yaml"

if [ ! -f "$META" ] || [ ! -f "$CONFIG" ]; then
    echo "ERROR: missing run_meta.json or resolved_config.yaml in $RUN_DIR" >&2
    exit 1
fi

# Extract uuid8, target, length from provenance.
# Use python for robust parsing.
HEADERS=$(python3 << PYEOF
import json, re, sys
import yaml

with open("$META") as f:
    meta = json.load(f)
with open("$CONFIG") as f:
    cfg = yaml.safe_load(f)

# uuid8: try meta["uuid"] or meta["uuid8"] or derive from run dir name
run_basename = "$(basename "$RUN_DIR")"
uuid_match = re.search(r"_([0-9a-f]{8})\$", run_basename)
uuid8 = uuid_match.group(1) if uuid_match else meta.get("uuid8", meta.get("uuid", "unknown"))

target_full = cfg.get("target_id", "csy4_crrna_free")
# Map to short form: csy4_crrna_free → csy4free
target_short = "csy4free" if target_full == "csy4_crrna_free" else target_full.replace("_", "")

# Length: try several places. Default 80 if not found.
length = 80
gp = cfg.get("generator_params", {})
# Look in any nested dict for "length", "len", "binder_length"
def find_length(d, depth=0):
    if depth > 4:
        return None
    if isinstance(d, dict):
        for k, v in d.items():
            if k in ("length", "binder_length", "len"):
                try:
                    return int(v)
                except:
                    pass
            r = find_length(v, depth + 1)
            if r is not None:
                return r
    elif isinstance(d, list):
        for v in d:
            r = find_length(v, depth + 1)
            if r is not None:
                return r
    return None
maybe = find_length(cfg)
if maybe is not None:
    length = maybe

print(f"{target_short}|{uuid8}|{length}")
PYEOF
)

if [ -z "$HEADERS" ]; then
    echo "ERROR: failed to extract metadata" >&2
    exit 1
fi

TARGET_SHORT=$(echo "$HEADERS" | cut -d'|' -f1)
UUID8=$(echo "$HEADERS" | cut -d'|' -f2)
LENGTH=$(echo "$HEADERS" | cut -d'|' -f3)

PREFIX="${TARGET_SHORT}_HD_LiMPNN_len${LENGTH}_${UUID8}"

echo "================================================================"
echo " HD postprocess"
echo " Run dir:     $RUN_DIR"
echo " Prefix:      $PREFIX"
echo "================================================================"

# Find all MPNN .pdb outputs
mapfile -t MPNN_FILES < <(find "$DESIGNS_DIR" -name "*ligand_mpnn_*.pdb" -type f | sort -V)

if [ ${#MPNN_FILES[@]} -eq 0 ]; then
    echo "WARNING: no *ligand_mpnn_*.pdb files found under $DESIGNS_DIR" >&2
    exit 0
fi

echo "Found ${#MPNN_FILES[@]} MPNN .pdb files to flatten"

n_renamed=0
for src in "${MPNN_FILES[@]}"; do
    fname=$(basename "$src")
    # Pattern: <scaffold>_recycle_<N>ligand_mpnn_<S>_<K>.pdb
    # Extract N, S, K
    cycle=$(echo "$fname" | sed -nE 's/.*_recycle_([0-9]+)ligand_mpnn_[0-9]+_[0-9]+\.pdb$/\1/p')
    sample=$(echo "$fname" | sed -nE 's/.*_recycle_[0-9]+ligand_mpnn_([0-9]+)_[0-9]+\.pdb$/\1/p')
    pack=$(echo "$fname" | sed -nE 's/.*_recycle_[0-9]+ligand_mpnn_[0-9]+_([0-9]+)\.pdb$/\1/p')

    if [ -z "$cycle" ] || [ -z "$sample" ] || [ -z "$pack" ]; then
        echo "  SKIP (unparseable): $fname" >&2
        continue
    fi

    new_name="${PREFIX}_recycle${cycle}_sample${sample}_pack${pack}.pdb"
    dst="$DESIGNS_DIR/$new_name"

    # Use cp -n to be idempotent (don't overwrite)
    if [ -f "$dst" ]; then
        continue
    fi
    cp "$src" "$dst"
    n_renamed=$((n_renamed + 1))
done

echo ""
echo "Canonical-named files added: $n_renamed"
echo "Total canonical at top of $DESIGNS_DIR:"
find "$DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*.pdb" 2>/dev/null | wc -l

# Copy processing_results.csv for easy reference at the canonical-files level
# (it lives in designs/processing_results.csv already, but ok to leave)

echo "================================================================"

# --- canonicalize designs/ -> .cif only; sweep intermediates to _work/ ---
# Retry: a transient apptainer/FS failure under concurrent launch must not
# abort the job (set -e) and orphan designs before the pack step below.
# canonicalize is idempotent, so retrying is safe; give up after 5 tries.
_ctry=0
until apptainer exec --bind /cluster/work/wicky --bind /tmp /cluster/work/wicky/containers/boltzgen.sif python3 /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/canonicalize_designs.py "$DESIGNS_DIR" "$(dirname "$DESIGNS_DIR")/_work"; do
    _ctry=$((_ctry + 1))
    [ "$_ctry" -ge 5 ] && { echo "canonicalize: failed after $_ctry attempts; leaving loose designs for recovery" >&2; exit 1; }
    echo "canonicalize: attempt $_ctry failed; backing off then retrying" >&2
    sleep $(( (RANDOM % 15) + 5 * _ctry ))
done

# --- Phase 1: reclaim inodes — delete disposables after canonicalize ---
source "$(dirname "$0")/archive_pipeline.sh" 2>/dev/null || true
declare -F purge >/dev/null || purge() { rm -rf -- "$@"; }
purge "$RUN_DIR/_work" || true                 # HD: recycle_*/ trees swept here

# Phase 2: pack designs/ for inode reduction (HD designs/ is always per-job).
if [ -n "${SHARED_DESIGNS_DIR:-}" ]; then
    # --num-runs>1: combine into the shared pool as a per-uuid slice (like multif).
    mkdir -p "$SHARED_DESIGNS_DIR"
    mv "$DESIGNS_DIR"/*_${UUID8}_* "$SHARED_DESIGNS_DIR"/ 2>/dev/null || true
    declare -F pack_slice >/dev/null && pack_slice "$SHARED_DESIGNS_DIR" "$UUID8" || true
    # combine done -> relocate this run's logs/meta to the shared pool dir (mirrors the
    # multif stacks' shared slurm_<uuid> logs), then drop the now-empty per-uuid staging
    # dir so it doesn't linger as a phantom run. mv keeps slurm's open output fd valid
    # (same inode) on Linux, so the rest of the job's log still lands in the moved file.
    SHARED_RUN_DIR="$(dirname "$SHARED_DESIGNS_DIR")"
    for _ext in log err; do
        if [ -f "$RUN_DIR/slurm.$_ext" ] && [ ! -e "$SHARED_RUN_DIR/slurm_${UUID8}.$_ext" ]; then
            mv -f "$RUN_DIR/slurm.$_ext" "$SHARED_RUN_DIR/slurm_${UUID8}.$_ext" || true
        fi
    done
    [ -f "$RUN_DIR/run_meta.json" ] && cp -f "$RUN_DIR/run_meta.json" "$SHARED_RUN_DIR/run_meta_${UUID8}.json" || true
    rm -rf "$RUN_DIR" || true
else
    declare -F pack_unit >/dev/null && pack_unit "$DESIGNS_DIR" || true   # single-run: whole dir
fi
