#!/bin/bash
# rfd3_multif.sh — RFD3 diffusion -> external LigandMPNN with mpnn_cli knobs.
#
# Mirrors boltzgen_multif.sh in shape: produce backbones, fan out external
# LiMPNN per backbone, canonicalize outputs.
#
# Arguments (positional, required):
#   $1  BACKBONE_DIR     shared backbone home (spec_0.{cif,pdb} ... spec_<N-1>.{cif,pdb})
#   $2  MPNN_RUN_DIR     rfd3__ligandmpnn run dir
#   $3  UUID8            shared uuid8
#   $4  LENGTH           protein binder length (used only for canonical naming)
#   $5  (optional)       MPNN_DESIGNS_DIR override (shared dir when --num-runs > 1)
#
# Environment variables (set by the launch.py-built sbatch wrapper):
#   RFD3_SIF, RFD3_ACTIVATE, RFD3_CWD, RFD3_INPUTS  -- RFD3 invocation
#   RFD3_N_BATCHES, RFD3_PREVALIDATE, RFD3_SKIP_EXISTING
#   MPNN_SIF, MPNN_WEIGHTS                          -- LigandMPNN apptainer + weights
#   MPNN_FLAGS                                       -- mpnn_cli flags joined as "--k v ..."

set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 BACKBONE_DIR MPNN_RUN_DIR UUID8 LENGTH [MPNN_DESIGNS_DIR]" >&2
    exit 2
fi

BACKBONE_DIR="$1"
MPNN_RUN_DIR="$2"
UUID8="$3"
LENGTH="$4"
MPNN_DESIGNS_DIR="${5:-$MPNN_RUN_DIR/designs}"

: "${RFD3_SIF:?RFD3_SIF must be set}"
: "${RFD3_ACTIVATE:?RFD3_ACTIVATE must be set}"
: "${RFD3_CWD:?RFD3_CWD must be set}"
: "${RFD3_INPUTS:?RFD3_INPUTS must be set}"
: "${MPNN_SIF:?MPNN_SIF must be set}"
: "${MPNN_WEIGHTS:?MPNN_WEIGHTS must be set}"

RFD3_N_BATCHES="${RFD3_N_BATCHES:-1}"
RFD3_DIFFUSION_BATCH_SIZE="${RFD3_DIFFUSION_BATCH_SIZE:-8}"   # backbones = n_batches × this (RFD3 default 8)
RFD3_PREVALIDATE="${RFD3_PREVALIDATE:-true}"
RFD3_SKIP_EXISTING="${RFD3_SKIP_EXISTING:-false}"
MPNN_FLAGS="${MPNN_FLAGS:-}"

MPNN_PREFIX="csy4free_RFD3_LiMPNN_len${LENGTH}_${UUID8}"

# --- idempotent restart (safe under #SBATCH --requeue) ---
# BACKBONE_DIR is uuid8-private -> wipe wholesale. The (possibly shared) designs
# dir is cleared ONLY of this run's own *_${UUID8}_* files, so parallel runs / siblings
# are untouched. (No-op on a fresh run; clears partial output on a requeue.)
rm -rf "$BACKBONE_DIR"
# uuid-private LigandMPNN scratch root - wipe ONLY this run subtree, never a sibling.
rm -rf "$MPNN_RUN_DIR/mpnn_native_${UUID8}"
find "$MPNN_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*" -delete 2>/dev/null || true
_WORK_SWEEP="$(dirname "$MPNN_DESIGNS_DIR")/_work_${UUID8}"
[ -d "$_WORK_SWEEP" ] && find "$_WORK_SWEEP" -maxdepth 1 -name "*_${UUID8}_*" -delete 2>/dev/null || true
mkdir -p "$BACKBONE_DIR"
mkdir -p "$MPNN_DESIGNS_DIR"

echo "================================================================"
echo " RFD3 multif worker"
echo " UUID8:             $UUID8"
echo " LENGTH:            $LENGTH"
echo " BACKBONE_DIR:      $BACKBONE_DIR"
echo " MPNN_RUN_DIR:      $MPNN_RUN_DIR"
echo " MPNN_DESIGNS_DIR:  $MPNN_DESIGNS_DIR"
echo " Started:           $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: RFD3 diffusion (+ bundled LiMPNN, sequence discarded)
# ============================================================
WORK_RFD3="$BACKBONE_DIR/rfd3_workdir"
mkdir -p "$WORK_RFD3"

echo ""
echo "=== Stage 1: RFD3 design ==="
# The activate script sets up frun + the foundry environment for the current shell.
# shellcheck disable=SC1090
source "$RFD3_ACTIVATE"
cd "$RFD3_CWD"

frun rfd3 design \
    out_dir="$WORK_RFD3" \
    inputs="$(basename "$RFD3_INPUTS")" \
    n_batches="$RFD3_N_BATCHES" \
    diffusion_batch_size="$RFD3_DIFFUSION_BATCH_SIZE" \
    prevalidate_inputs="$RFD3_PREVALIDATE" \
    skip_existing="$RFD3_SKIP_EXISTING"

# ============================================================
# Stage 2: collect raw backbones, decompress + CIF->PDB convert
# ============================================================
echo ""
echo "=== Stage 2: collect + convert backbones ==="
shopt -s nullglob
ALL_BACKBONES=("$WORK_RFD3"/*.cif.gz)
shopt -u nullglob

if [ ${#ALL_BACKBONES[@]} -eq 0 ]; then
    echo "ERROR: no RFD3 backbones found in $WORK_RFD3" >&2
    ls -la "$WORK_RFD3" 2>/dev/null || true
    exit 1
fi
N_BACKBONES=${#ALL_BACKBONES[@]}
echo "Detected $N_BACKBONES backbones from RFD3"

# Decompress in sorted order so spec_0..spec_{N-1} indices are stable.
idx=0
for f in $(ls "$WORK_RFD3"/*.cif.gz | sort); do
    dst_cif="$BACKBONE_DIR/spec_${idx}.cif"
    gunzip -c "$f" > "$dst_cif"
    idx=$((idx + 1))
done

# Convert all .cif -> .pdb. Foundry SIF doesn't ship biopython; LigandMPNN SIF
# does (1.79) and is already loaded for stage 3, so we use it for conversion too.
apptainer exec --bind /cluster/work/wicky --bind /tmp \
    "$MPNN_SIF" \
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

# Backbone meta (consumed_by_stacks identifies the only downstream stack)
cat > "$BACKBONE_DIR/meta.json" <<EOF
{
  "uuid8": "$UUID8",
  "generator": "rfd3",
  "length": $LENGTH,
  "num_backbones": $N_BACKBONES,
  "consumed_by_stacks": [
    "rfd3__ligandmpnn/$(basename "$MPNN_RUN_DIR")"
  ]
}
EOF

# ============================================================
# Stage 3: LigandMPNN -- one invocation per backbone
# ============================================================
echo ""
echo "=== Stage 3: LigandMPNN over $N_BACKBONES backbones ==="

# --- RFD3 chain-label fix (2026-06-25) ----------------------------------------
# RFD3 foundry .cif orders the RNA first (chain A) and the protein binder second
# (chain B); the cif->pdb conversion preserves that order, so the configured
# `--chains_to_design A` pointed at the RNA -> LigandMPNN designed 0 residues
# (num_res=0) and just echoed the foundry-bundled sequence identically across all
# batch samples (omit_AA/temperature/seed never applied). Detect the real protein
# chain and design THAT; chain labels are left untouched (RNA stays ligand context).
detect_protein_chain() {
    awk 'BEGIN{split("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL",_a," ");for(_i in _a)aa[_a[_i]]=1}
         (/^ATOM/||/^HETATM/){resn=substr($0,18,3);gsub(/ /,"",resn);ch=substr($0,22,1);if(resn in aa){print ch;exit}}' "$1"
}
# ------------------------------------------------------------------------------

for idx in $(seq 0 $((N_BACKBONES - 1))); do
    pdb_in="$BACKBONE_DIR/spec_${idx}.pdb"
    if [ ! -f "$pdb_in" ]; then
        echo "  WARN: missing $pdb_in, skipping design $idx"
        continue
    fi
    mpnn_out="$MPNN_RUN_DIR/mpnn_native_${UUID8}/design${idx}"
    mkdir -p "$mpnn_out"
    echo ""
    echo "  --- LigandMPNN on design $idx (input: $pdb_in) ---"

    # Pick the protein chain at run time (see detect_protein_chain above) and strip
    # any hardcoded --chains_to_design from MPNN_FLAGS so detection wins.
    prot_chain="$(detect_protein_chain "$pdb_in")"
    if [ -z "$prot_chain" ]; then
        echo "  WARN: no protein chain detected in $pdb_in; using MPNN_FLAGS as-is" >&2
        mpnn_flags_run="$MPNN_FLAGS"
    else
        mpnn_flags_run="--chains_to_design $prot_chain $(printf %s "$MPNN_FLAGS" | sed -E 's/--chains_to_design[[:space:]]+[A-Za-z0-9]+[[:space:]]*//')"
        echo "  protein chain = $prot_chain (designing it; RNA kept as ligand context)"
    fi

    apptainer exec --nv --bind /cluster/work/wicky --bind /tmp \
        "$MPNN_SIF" \
        bash -c "cd /app/LigandMPNN && python3 run.py \
            --model_type ligand_mpnn \
            --checkpoint_ligand_mpnn $MPNN_WEIGHTS \
            --pdb_path $pdb_in \
            --out_folder $mpnn_out \
            $mpnn_flags_run"

    # Canonicalize per-backbone outputs into the shared MPNN_DESIGNS_DIR.
    cd "$MPNN_DESIGNS_DIR"

    # FASTA: one per backbone, holds all sequences for this design.
    if [ -f "$mpnn_out/seqs/spec_${idx}.fa" ]; then
        cp "$mpnn_out/seqs/spec_${idx}.fa" "${MPNN_PREFIX}_design${idx}_sequences.fa"
    elif [ -f "$mpnn_out/seqs/spec.fa" ]; then
        cp "$mpnn_out/seqs/spec.fa" "${MPNN_PREFIX}_design${idx}_sequences.fa"
    fi

    # Repacked backbones: spec_<idx>_<N>.pdb or spec_<N>.pdb.
    shopt -s nullglob
    for pdb in "$mpnn_out/backbones"/*.pdb; do
        bname=$(basename "$pdb" .pdb)
        seqidx=$(echo "$bname" | grep -oP '_\K\d+$' || true)
        [ -z "$seqidx" ] && continue
        s_idx=$((seqidx - 1))
        cp "$pdb" "${MPNN_PREFIX}_design${idx}_seq${s_idx}.pdb"
    done
    shopt -u nullglob
done

MPNN_COUNT=$(find "$MPNN_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID8}_*.pdb" 2>/dev/null | wc -l)

echo ""
echo "================================================================"
echo " Done at $(date)"
echo " Backbones: $N_BACKBONES in $BACKBONE_DIR/"
echo " MPNN designs: $MPNN_COUNT files"
echo " MPNN designs dir: $MPNN_DESIGNS_DIR"
echo "================================================================"

# --- canonicalize designs/ -> .cif only; sweep intermediates to _work/ ---
# Scope to this job's own files ("_<uuid8>_") so concurrent --num-runs>1 jobs
# sharing this designs/ dir don't race on each other's outputs.
# Retry: a transient apptainer/FS failure under concurrent launch must not
# abort the job (set -e) and orphan designs before the pack step below.
# canonicalize is idempotent, so retrying is safe; give up after 5 tries.
_ctry=0
until apptainer exec --bind /cluster/work/wicky --bind /tmp /cluster/work/wicky/containers/boltzgen.sif python3 /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/canonicalize_designs.py "$MPNN_DESIGNS_DIR" "$(dirname "$MPNN_DESIGNS_DIR")/_work_${UUID8}" "_${UUID8}_"; do
    _ctry=$((_ctry + 1))
    [ "$_ctry" -ge 5 ] && { echo "canonicalize: failed after $_ctry attempts; leaving loose designs for recovery" >&2; exit 1; }
    echo "canonicalize: attempt $_ctry failed; backing off then retrying" >&2
    sleep $(( (RANDOM % 15) + 5 * _ctry ))
done

# --- Phase 1: reclaim inodes — delete disposables after canonicalize ---
source "$(dirname "$0")/archive_pipeline.sh" 2>/dev/null || true
declare -F purge >/dev/null || purge() { rm -rf -- "$@"; }
purge "$MPNN_RUN_DIR/mpnn_native_${UUID8}" || true
# The _work sweep dir is SHARED across jobs when --num-runs>1 (it is dirname of
# the shared designs/, not $MPNN_RUN_DIR/_work). Purge only THIS job's uuid-scoped
# swept entries (files + intermediate dirs), then rmdir if now empty (last job out),
# so concurrent jobs never race on it. Mirrors the restart-cleanup scoping above.
[ -d "$_WORK_SWEEP" ] && find "$_WORK_SWEEP" -maxdepth 1 -name "*_${UUID8}_*" -exec rm -rf {} + 2>/dev/null || true
rmdir "$_WORK_SWEEP" 2>/dev/null || true
# --num-runs>1: MPNN_RUN_DIR is per-uuid scaffolding (designs/logs/sbatch live
# in the shared <label>/ dir), so the purge above leaves it EMPTY. rmdir it —
# a no-op in single-run, where it holds designs/ + logs and rmdir fails.
rmdir "$MPNN_RUN_DIR" 2>/dev/null || true
# Tier B (optional): raw diffusion backbones
purge "$WORK_RFD3" || true                 # = $BACKBONE_DIR/rfd3_workdir

# Phase 2: pack designs/ for inode reduction (single-run only; a shared
# --num-runs>1 designs/ must be packed once after all N jobs, not per-job here).
if [[ "$MPNN_DESIGNS_DIR" == "$MPNN_RUN_DIR/designs" ]]; then
    declare -F pack_unit >/dev/null && pack_unit "$MPNN_DESIGNS_DIR" || true
else
    # --num-runs>1: shared designs/, so pack only THIS run's *_<uuid8>_* slice
    # (race-safe) -> the summary dir becomes N per-uuid slice archives.
    declare -F pack_slice >/dev/null && pack_slice "$MPNN_DESIGNS_DIR" "$UUID8" || true
fi
