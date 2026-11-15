#!/bin/bash
# boltzgen_multif.sh — BoltzGen design (N backbones) → both inverse folders.
#
# Handles num_designs >= 1. Each backbone gets design<D> index (D = 0..N-1).
# Runs BG IF on all backbones together; loops LigandMPNN per-backbone.
#
# Arguments (positional, required):
#   $1  BACKBONE_DIR    shared backbone home (spec_0.{cif,pdb} ... spec_<N-1>.{cif,pdb})
#   $2  BIF_RUN_DIR     boltzgen__boltzif run dir
#   $3  MPNN_RUN_DIR    boltzgen__ligandmpnn run dir
#   $4  UUID8           shared uuid8
#   $5  LENGTH          protein binder length
#
# Environment variables (set by sbatch wrapper):
#   BG_SIF, BG_HOME, BG_BG_ENVS    — BoltzGen apptainer config
#   MPNN_SIF, MPNN_WEIGHTS         — LigandMPNN apptainer + weights
#   BG_FLAGS, BIF_FLAGS, MPNN_FLAGS — CLI flags for each stage
#   TARGET_RNA_PDB, SPEC_YAML_CONTENT

set -euo pipefail

if [ $# -lt 5 ]; then
    echo "Usage: $0 BACKBONE_DIR BIF_RUN_DIR MPNN_RUN_DIR UUID8 LENGTH [BIF_DESIGNS_DIR] [MPNN_DESIGNS_DIR]" >&2
    exit 2
fi

BACKBONE_DIR="$1"
BIF_RUN_DIR="$2"
MPNN_RUN_DIR="$3"
UUID8="$4"
LENGTH="$5"
# Optional shared designs dirs (--num-runs > 1)
BIF_DESIGNS_DIR="${6:-$BIF_RUN_DIR/designs}"
MPNN_DESIGNS_DIR="${7:-$MPNN_RUN_DIR/designs}"

: "${BG_SIF:?BG_SIF must be set}"
: "${BG_HOME:?BG_HOME must be set}"
: "${MPNN_SIF:?MPNN_SIF must be set}"
: "${MPNN_WEIGHTS:?MPNN_WEIGHTS must be set}"
: "${TARGET_RNA_PDB:?TARGET_RNA_PDB must be set}"
: "${SPEC_YAML_CONTENT:?SPEC_YAML_CONTENT must be set}"

BG_FLAGS="${BG_FLAGS:-}"
BIF_FLAGS="${BIF_FLAGS:-}"
MPNN_FLAGS="${MPNN_FLAGS:-}"
BG_BG_ENVS="${BG_BG_ENVS:-}"

BIF_PREFIX="csy4free_BG_BIF_len${LENGTH}_${UUID8}"
MPNN_PREFIX="csy4free_BG_LiMPNN_len${LENGTH}_${UUID8}"

# --- idempotent restart (safe under #SBATCH --requeue) ---
# BACKBONE_DIR is uuid8-private -> wipe wholesale. The (possibly shared) BIF + LiMPNN
# designs dirs are cleared ONLY of this run's own *_${UUID8}_* files (parallel-safe).
rm -rf "$BACKBONE_DIR"
# uuid-private LigandMPNN scratch root: in --num-runs>1 it lives in the SHARED
# designs dir, so wipe ONLY this run subtree, never a sibling (mpnn_native race fix).
rm -rf "$MPNN_RUN_DIR/mpnn_native_${UUID8}"
for _DD in "$BIF_DESIGNS_DIR" "$MPNN_DESIGNS_DIR"; do
    find "$_DD" -maxdepth 1 -name "*_${UUID8}_*" -delete 2>/dev/null || true
    _WS="$(dirname "$_DD")/_work_${UUID8}"
    [ -d "$_WS" ] && find "$_WS" -maxdepth 1 -name "*_${UUID8}_*" -delete 2>/dev/null || true
done
mkdir -p "$BACKBONE_DIR"
mkdir -p "$BIF_RUN_DIR" "$MPNN_RUN_DIR"   # per-uuid run dirs (never created under --num-runs>1; needed by BIF harvest cp)
mkdir -p "$BIF_DESIGNS_DIR"
mkdir -p "$MPNN_DESIGNS_DIR"

echo "================================================================"
echo " BoltzGen multif worker (N-backbones aware)"
echo " UUID8:         $UUID8"
echo " LENGTH:        $LENGTH"
echo " BACKBONE_DIR:  $BACKBONE_DIR"
echo " BIF_RUN_DIR:   $BIF_RUN_DIR"
echo " MPNN_RUN_DIR:  $MPNN_RUN_DIR"
echo " Started:       $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: BoltzGen design → N backbones in intermediate_designs/
# ============================================================
WORK_BG="$BACKBONE_DIR/bg_workdir"
mkdir -p "$WORK_BG"
cd "$WORK_BG"

cp "$TARGET_RNA_PDB" .
printf '%s\n' "$SPEC_YAML_CONTENT" > spec.yaml

echo ""
echo "=== Stage 1: BoltzGen design ==="
apptainer exec --nv \
    --home "$BG_HOME" \
    --bind /cluster/work/wicky \
    $BG_BG_ENVS \
    "$BG_SIF" \
    boltzgen run spec.yaml \
        --output "$WORK_BG/bg_output" \
        --steps design \
        $BG_FLAGS

# ============================================================
# Discover all backbone files BG produced
# ============================================================
INTERMEDIATE_DIR="$WORK_BG/bg_output/intermediate_designs"
shopt -s nullglob
ALL_BACKBONES=("$INTERMEDIATE_DIR"/spec*.cif)
shopt -u nullglob

if [ ${#ALL_BACKBONES[@]} -eq 0 ]; then
    echo "ERROR: no BG backbones found in $INTERMEDIATE_DIR" >&2
    ls -la "$INTERMEDIATE_DIR" 2>/dev/null || true
    exit 1
fi

# Sort numerically: "spec.cif" → 0, "spec_0.cif" → 1, etc.? No — BG with num_designs=1
# produces "spec.cif" (no index); with num_designs>1 produces "spec_0.cif"..."spec_{N-1}.cif".
# Handle both cases by mapping:
#   spec.cif        → design 0
#   spec_<K>.cif    → design <K>
# When BG produces "spec_0.cif" through "spec_N-1.cif", design indices match the suffix.

declare -A BACKBONE_BY_DESIGN     # design_idx → filename (basename)
declare -A BACKBONE_BY_KEY        # internal key → design index, used by BIF later

for f in "${ALL_BACKBONES[@]}"; do
    base=$(basename "$f" .cif)
    if [ "$base" = "spec" ]; then
        # Single-design case
        idx=0
        key="spec"
    else
        # Multi-design case: spec_<K>
        idx=$(echo "$base" | grep -oP '_\K\d+$')
        key="$base"   # "spec_3"
    fi
    BACKBONE_BY_DESIGN[$idx]="$base.cif"
    BACKBONE_BY_KEY[$key]=$idx
done

N_BACKBONES=${#BACKBONE_BY_DESIGN[@]}
echo "Detected $N_BACKBONES backbones from BG"

# Copy each backbone into the shared _backbones/ dir + convert to PDB
echo ""
echo "=== Copying + converting backbones ==="
for idx in $(echo "${!BACKBONE_BY_DESIGN[@]}" | tr ' ' '\n' | sort -n); do
    src_cif="$INTERMEDIATE_DIR/${BACKBONE_BY_DESIGN[$idx]}"
    dst_cif="$BACKBONE_DIR/spec_${idx}.cif"
    cp "$src_cif" "$dst_cif"
done

# Batch-convert all CIFs to PDBs in one apptainer call
apptainer exec --bind /cluster/work/wicky --bind /tmp \
    "$BG_SIF" \
    python3 -c "
import os
from Bio.PDB import MMCIFParser, PDBIO
bdir = '$BACKBONE_DIR'
p = MMCIFParser(QUIET=True)
io = PDBIO()
for fn in sorted(os.listdir(bdir)):
    if not fn.endswith('.cif'):
        continue
    cif_path = os.path.join(bdir, fn)
    pdb_path = cif_path.replace('.cif', '.pdb')
    if os.path.exists(pdb_path):
        continue
    s = p.get_structure(fn, cif_path)
    io.set_structure(s)
    io.save(pdb_path)
    print(f'  {fn} -> {os.path.basename(pdb_path)}')
"

# Backbone meta
cat > "$BACKBONE_DIR/meta.json" <<EOF
{
  "uuid8": "$UUID8",
  "generator": "boltzgen",
  "length": $LENGTH,
  "num_backbones": $N_BACKBONES,
  "consumed_by_stacks": [
    "boltzgen__boltzif/$(basename $BIF_RUN_DIR)",
    "boltzgen__ligandmpnn/$(basename $MPNN_RUN_DIR)"
  ]
}
EOF

# ============================================================
# Stage 2: BoltzGen --only_inverse_fold (runs IF + folding + analysis + filtering)
# ============================================================
echo ""
echo "=== Stage 2: BoltzGen BIF (inverse fold + analysis) ==="
cd "$WORK_BG"

apptainer exec --nv \
    --home "$BG_HOME" \
    --bind /cluster/work/wicky \
    --bind /cluster/work/wicky/lvalsecchi/rna_binder/boltzgen_patch/writer.py:/app/src/boltzgen/task/predict/writer.py \
    $BG_BG_ENVS \
    "$BG_SIF" \
    boltzgen run spec.yaml \
        --output "$WORK_BG/bg_output" \
        --steps inverse_folding folding \
        $BIF_FLAGS

# Locate ranked BIF outputs
# --steps inverse_folding routes outputs to intermediate_designs_inverse_folded/;
# --only_inverse_fold wrote intermediate_designs/ but collapsed to 1 seq/backbone.
if [ -d "$WORK_BG/bg_output/intermediate_designs_inverse_folded/refold_cif" ]; then
    BIF_NATIVE="$WORK_BG/bg_output/intermediate_designs_inverse_folded/refold_cif"
else
    BIF_NATIVE="$WORK_BG/bg_output/intermediate_designs/refold_cif"
fi

if [ -d "$BIF_NATIVE" ]; then
    cp -r "$BIF_NATIVE" "$BIF_RUN_DIR/bif_native"
    echo "BIF native outputs copied to $BIF_RUN_DIR/bif_native"

    # Filename pattern: rank<R>_spec.cif  OR  rank<R>_spec_<D>.cif
    # We need to extract the design index (which backbone) from the suffix, and
    # assign a sequence index (0, 1, 2, ...) within that backbone in the order
    # we encounter rank-sorted files.

    cd "$BIF_DESIGNS_DIR"

    shopt -s nullglob
    # Process in rank order: rank1, rank2, ...
    for f in $(ls "$BIF_RUN_DIR/bif_native"/spec_*.cif | sort -V); do
        base=$(basename "$f" .cif)   # spec_<bb>_<seq> (N/bb, patched writer) or spec_<bb>
        if [[ "$base" =~ spec_([0-9]+)_([0-9]+)$ ]]; then
            d_idx=$((10#${BASH_REMATCH[1]})); s_idx=$((10#${BASH_REMATCH[2]}))
            new="${BIF_PREFIX}_design${d_idx}_seq${s_idx}.cif"
            npz_new="${BIF_PREFIX}_design${d_idx}_seq${s_idx}.npz"
        elif [[ "$base" =~ spec_([0-9]+)$ ]]; then
            d_idx=$((10#${BASH_REMATCH[1]}))   # base-10 (avoid octal on 08/09)
            new="${BIF_PREFIX}_design${d_idx}.cif"
            npz_new="${BIF_PREFIX}_design${d_idx}.npz"
        else
            echo "Skipping unparseable filename: $base"
            continue
        fi
        cp "$f" "$new"
        if [ -f "${f%.cif}.npz" ]; then
            cp "${f%.cif}.npz" "$npz_new"
        fi
    done
    shopt -u nullglob

    BIF_COUNT=$(find "$BIF_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*.cif" 2>/dev/null | wc -l)
    echo "BIF canonical: $BIF_COUNT files written to $BIF_RUN_DIR/designs/"

    # Copy metrics CSV
    METRICS_CSV="$WORK_BG/bg_output/final_ranked_designs/all_designs_metrics.csv"
    if [ -f "$METRICS_CSV" ]; then
        cp "$METRICS_CSV" "$BIF_RUN_DIR/"
    fi
else
    echo "WARNING: no BIF native output at $BIF_NATIVE" >&2
fi

# ============================================================
# Stage 3: LigandMPNN — one invocation per backbone
# ============================================================
echo ""
echo "=== Stage 3: LigandMPNN over $N_BACKBONES backbones ==="

for idx in $(echo "${!BACKBONE_BY_DESIGN[@]}" | tr ' ' '\n' | sort -n); do
    pdb_in="$BACKBONE_DIR/spec_${idx}.pdb"
    if [ ! -f "$pdb_in" ]; then
        echo "  WARN: missing $pdb_in, skipping design $idx"
        continue
    fi

    mpnn_out="$MPNN_RUN_DIR/mpnn_native_${UUID8}/design${idx}"
    mkdir -p "$mpnn_out"

    echo ""
    echo "  --- LigandMPNN on design $idx (input: $pdb_in) ---"

    apptainer exec --nv --bind /cluster/work/wicky --bind /tmp \
        "$MPNN_SIF" \
        bash -c "cd /app/LigandMPNN && python3 run.py \
            --model_type ligand_mpnn \
            --checkpoint_ligand_mpnn $MPNN_WEIGHTS \
            --pdb_path $pdb_in \
            --out_folder $mpnn_out \
            $MPNN_FLAGS"

    # Rename outputs into canonical. Normalize the backbone index to base-10
    # (BoltzGen zero-pads, e.g. spec_08) so canonical names are unpadded
    # design<N> — matching the BIF path, rfd3, and parse_filename's int round-trip.
    d_idx=$((10#$idx))
    cd "$MPNN_DESIGNS_DIR"

    # FASTA: one per backbone, holds all sequences for this design
    if [ -f "$mpnn_out/seqs/spec_${idx}.fa" ]; then
        cp "$mpnn_out/seqs/spec_${idx}.fa" "${MPNN_PREFIX}_design${d_idx}_sequences.fa"
    elif [ -f "$mpnn_out/seqs/spec.fa" ]; then
        # Fallback: when LiMPNN strips _N from input filename
        cp "$mpnn_out/seqs/spec.fa" "${MPNN_PREFIX}_design${d_idx}_sequences.fa"
    fi

    # Repacked backbones: spec_<idx>_<N>.pdb or spec_<N>.pdb (depending on LiMPNN naming)
    shopt -s nullglob
    for pdb in "$mpnn_out/backbones"/*.pdb; do
        bname=$(basename "$pdb" .pdb)
        # Extract the trailing _N which is the LiMPNN sequence index (1-indexed)
        seqidx=$(echo "$bname" | grep -oP '_\K\d+$' || true)
        [ -z "$seqidx" ] && continue
        s_idx=$((seqidx - 1))
        cp "$pdb" "${MPNN_PREFIX}_design${d_idx}_seq${s_idx}.pdb"
    done
    shopt -u nullglob
done

MPNN_COUNT=$(find "$MPNN_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*.pdb" 2>/dev/null | wc -l)

echo ""
echo "================================================================"
echo " Done at $(date)"
echo " Backbones: $N_BACKBONES in $BACKBONE_DIR/"
echo " BIF designs: $(find "$BIF_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*.cif" 2>/dev/null | wc -l) files"
echo " MPNN designs: $MPNN_COUNT files"
 echo " BIF designs dir:  $BIF_DESIGNS_DIR"
 echo " MPNN designs dir: $MPNN_DESIGNS_DIR"
echo "================================================================"

# --- canonicalize both inverse-folder designs/ dirs ---
# Scope to this job's own files ("_<uuid8>_") so concurrent --num-runs>1 jobs
# sharing these designs/ dirs don't race on each other's outputs.
for DD in "$BIF_DESIGNS_DIR" "$MPNN_DESIGNS_DIR"; do
  # Retry: a transient apptainer/FS failure under concurrent launch must not
  # abort the job (set -e) and orphan designs before the pack step below.
  # canonicalize is idempotent, so retrying is safe; give up after 5 tries.
  _ctry=0
  until apptainer exec --bind /cluster/work/wicky --bind /tmp /cluster/work/wicky/containers/boltzgen.sif python3 /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/canonicalize_designs.py "$DD" "$(dirname "$DD")/_work_${UUID8}" "_${UUID8}_"; do
      _ctry=$((_ctry + 1))
      [ "$_ctry" -ge 5 ] && { echo "canonicalize: failed after $_ctry attempts; leaving loose designs for recovery" >&2; exit 1; }
      echo "canonicalize: attempt $_ctry failed; backing off then retrying" >&2
      sleep $(( (RANDOM % 15) + 5 * _ctry ))
  done
done

# --- Phase 1: reclaim inodes — delete disposables after canonicalize ---
# Two designs dirs (BIF + MPNN) -> two _work dirs. mpnn_native_* is MPNN-side only.
source "$(dirname "$0")/archive_pipeline.sh" 2>/dev/null || true
declare -F purge >/dev/null || purge() { rm -rf -- "$@"; }
purge "$MPNN_RUN_DIR/mpnn_native_${UUID8}" || true
# Raw BoltzIF native specs (already canonicalized into the designs/slices);
# bif_native has no _<uuid>_ so canonicalize/pack_slice never touch it.
purge "$BIF_RUN_DIR/bif_native" "$BIF_DESIGNS_DIR/bif_native" || true
# Two _work sweep dirs (BIF + MPNN), each SHARED when --num-runs>1 (dirname of the
# shared designs/). Purge only THIS job's uuid-scoped swept entries, then rmdir each
# if empty, so concurrent jobs never race. Mirrors the restart-cleanup above.
for _DD in "$MPNN_DESIGNS_DIR" "$BIF_DESIGNS_DIR"; do
    _WS="$(dirname "$_DD")/_work_${UUID8}"
    [ -d "$_WS" ] && find "$_WS" -maxdepth 1 -name "*_${UUID8}_*" -exec rm -rf {} + 2>/dev/null || true
    rmdir "$_WS" 2>/dev/null || true
done
# --num-runs>1: BIF_RUN_DIR / MPNN_RUN_DIR are per-uuid scaffolding (designs +
# logs live in the shared <label>/ dir), so the purge leaves them EMPTY. rmdir
# both — no-op in single-run where they hold designs/ + logs.
rmdir "$MPNN_RUN_DIR" 2>/dev/null || true
rmdir "$BIF_RUN_DIR" 2>/dev/null || true
# Tier B (optional): BoltzGen backbone workdir
purge "$WORK_BG" || true                   # = $BACKBONE_DIR/bg_workdir

# Phase 2: pack designs/ for inode reduction (single-run only; shared
# --num-runs>1 designs/ dirs must be packed once after all N jobs, not here).
if [[ "$MPNN_DESIGNS_DIR" == "$MPNN_RUN_DIR/designs" ]]; then
    declare -F pack_unit >/dev/null && pack_unit "$MPNN_DESIGNS_DIR" || true
else
    # --num-runs>1: shared designs/, so pack only THIS run's *_<uuid8>_* slice
    # (race-safe) -> the summary dir becomes N per-uuid slice archives.
    declare -F pack_slice >/dev/null && pack_slice "$MPNN_DESIGNS_DIR" "$UUID8" || true
fi
if [[ "$BIF_DESIGNS_DIR" == "$BIF_RUN_DIR/designs" ]]; then
    declare -F pack_unit >/dev/null && pack_unit "$BIF_DESIGNS_DIR" || true
else
    declare -F pack_slice >/dev/null && pack_slice "$BIF_DESIGNS_DIR" "$UUID8" || true
fi
