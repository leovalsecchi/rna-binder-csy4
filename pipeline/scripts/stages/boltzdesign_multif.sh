#!/bin/bash
# boltzdesign_multif.sh — BD1 hallucination only -> external LigandMPNN.
#
# Mirrors rfd3_multif.sh in shape. BD1 runs with --run_ligandmpnn false and
# --run_alphafold false (set in the YAML's generator_params), so the bundled
# MPNN + Boltz1 refold post-stages are skipped. We collect the surviving
# hallucinated backbones from outputs/<target_type>_<target_name>_<suffix>/pdb/
# and fan out external LiMPNN with our mpnn_cli knobs.
#
# Arguments (positional, required):
#   $1  BACKBONE_DIR     shared backbone home (spec_<itr>.pdb after collection)
#   $2  MPNN_RUN_DIR     boltzdesign__ligandmpnn run dir
#   $3  UUID8            shared uuid8
#   $4  LENGTH           protein binder length (used only for canonical naming)
#   $5  (optional)       MPNN_DESIGNS_DIR override (shared dir when --num-runs > 1)
#
# Environment variables (set by the launch.py-built sbatch wrapper):
#   BD1_SIF, BD1_BIND_FLAGS, BD1_ENV_FLAGS    -- BD1 apptainer args
#   BD1_CLI_ARGS                              -- "--target_name X --pdb_path Y ..."
#   BD1_TARGET_TYPE, BD1_TARGET_NAME, BD1_SUFFIX  -- for finding outputs/<...>/pdb/
#   MPNN_SIF, MPNN_WEIGHTS                    -- LigandMPNN apptainer + weights
#   MPNN_FLAGS                                -- mpnn_cli flags joined as "--k v ..."

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

: "${BD1_SIF:?BD1_SIF must be set}"
: "${BD1_CLI_ARGS:?BD1_CLI_ARGS must be set}"
: "${BD1_TARGET_TYPE:?BD1_TARGET_TYPE must be set}"
: "${BD1_TARGET_NAME:?BD1_TARGET_NAME must be set}"
: "${BD1_SUFFIX:?BD1_SUFFIX must be set}"
: "${BOLTZ_CACHE:?BOLTZ_CACHE must be set}"
: "${MPNN_SIF:?MPNN_SIF must be set}"
: "${MPNN_WEIGHTS:?MPNN_WEIGHTS must be set}"

BD1_BIND_FLAGS="${BD1_BIND_FLAGS:-}"
MPNN_FLAGS="${MPNN_FLAGS:-}"

MPNN_PREFIX="csy4free_BD1_LiMPNN_len${LENGTH}_${UUID8}"

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
echo " BoltzDesign1 multif worker"
echo " UUID8:             $UUID8"
echo " LENGTH:            $LENGTH"
echo " BACKBONE_DIR:      $BACKBONE_DIR"
echo " MPNN_RUN_DIR:      $MPNN_RUN_DIR"
echo " MPNN_DESIGNS_DIR:  $MPNN_DESIGNS_DIR"
echo " Started:           $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: BD1 hallucination (bundled MPNN + Boltz refold disabled)
# ============================================================
WORK_BD1="$BACKBONE_DIR/bd1_workdir"
mkdir -p "$WORK_BD1"

# Boltz cache: idempotent setup mirroring the bundled BD1 launcher.
mkdir -p "$BOLTZ_CACHE/.boltz"
ln -sfn /cluster/work/wicky/params/boltz/boltz1_conf.ckpt "$BOLTZ_CACHE/.boltz/boltz1_conf.ckpt"
ln -sfn /cluster/work/wicky/params/boltz/ccd.pkl "$BOLTZ_CACHE/.boltz/ccd.pkl"

# Fake HOME for matplotlib + container writes (must live on a writable path
# inside the container — TMPDIR is per-job and bind-mounted).
FAKE_HOME="${TMPDIR}/fake_home"
mkdir -p "$FAKE_HOME"
export MPLCONFIGDIR="${FAKE_HOME}/.config/matplotlib"
mkdir -p "$MPLCONFIGDIR"

echo ""
echo "=== Stage 1: BD1 hallucination (--run_ligandmpnn false --run_alphafold false) ==="

apptainer exec \
    --nv \
    $BD1_BIND_FLAGS \
    --bind "${TMPDIR}:${TMPDIR}" \
    --bind "${BOLTZ_CACHE}:/cluster/home/${USER}" \
    --env HOME="${FAKE_HOME}" \
    --env MPLCONFIGDIR="${MPLCONFIGDIR}" \
    "$BD1_SIF" \
    bash -c "
set -euo pipefail
cd /app
python boltzdesign.py $BD1_CLI_ARGS --gpu_id 0 --work_dir '$WORK_BD1'
"

# ============================================================
# Stage 2: collect surviving hallucinated backbones
# ============================================================
echo ""
echo "=== Stage 2: collect hallucinated backbones ==="

BD1_OUT_DIR="$WORK_BD1/outputs/${BD1_TARGET_TYPE}_${BD1_TARGET_NAME}_${BD1_SUFFIX}/pdb"
BD1_RESULTS_DIR="$WORK_BD1/outputs/${BD1_TARGET_TYPE}_${BD1_TARGET_NAME}_${BD1_SUFFIX}/results_final"
if [ ! -d "$BD1_OUT_DIR" ]; then
    # i_ptm_cutoff=0 skips the pdb/ stage — extract PDBs from results_final CIFs instead
    echo "  pdb/ not found, falling back to results_final/ CIF extraction" >&2
    if [ ! -d "$BD1_RESULTS_DIR" ]; then
        echo "ERROR: BD1 didn't create $BD1_RESULTS_DIR. Inspect $WORK_BD1/outputs/" >&2
        find "$WORK_BD1/outputs" -maxdepth 3 -type d 2>/dev/null | head -20 >&2 || true
        exit 1
    fi
    mkdir -p "$BD1_OUT_DIR"
    # Convert each result CIF to PDB using python/biotite
    apptainer exec --nv $BD1_BIND_FLAGS         --bind /cluster/work/wicky         /cluster/work/wicky/containers/boltz_design1.sif         bash -c "
            python3 - <<PYEOF
import glob, os
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
io = PDBIO()
results_dir = '${BD1_RESULTS_DIR}'
out_dir = '${BD1_OUT_DIR}'
for cif in sorted(glob.glob(os.path.join(results_dir, '*/predictions/*/*.cif'))):
    base = os.path.basename(cif).replace('.cif', '.pdb')
    # extract itr index from path
    out_pdb = os.path.join(out_dir, base)
    try:
        struct = parser.get_structure('s', cif)
        io.set_structure(struct)
        io.save(out_pdb)
    except Exception as e:
        print(f'  SKIP {base}: {e}')
PYEOF
        "
fi

shopt -s nullglob
ALL_BACKBONES=("$BD1_OUT_DIR"/*_itr*_length*_model_0.pdb)
shopt -u nullglob

if [ ${#ALL_BACKBONES[@]} -eq 0 ]; then
    echo "ERROR: no BD1 backbones found in $BD1_OUT_DIR" >&2
    ls -la "$BD1_OUT_DIR" 2>/dev/null | head -20 >&2 || true
    exit 1
fi
N_BACKBONES=${#ALL_BACKBONES[@]}
echo "Detected $N_BACKBONES surviving hallucinated backbones"

# Use the BD1 iteration number from each source filename as the design index.
# Sort numerically so output order is stable. We map iter<N> -> spec_<N>.pdb
# without renumbering, so MPNN outputs keep BD1's iter index for cross-ref.
declare -A SEEN_IDX
for f in $(ls "$BD1_OUT_DIR"/*_itr*_length*_model_0.pdb | sort -V); do
    base=$(basename "$f")
    if [[ "$base" =~ _itr([0-9]+)_length ]]; then
        idx=$((10#${BASH_REMATCH[1]}))
    else
        echo "  SKIP (unparseable): $base" >&2
        continue
    fi
    if [ -n "${SEEN_IDX[$idx]:-}" ]; then
        echo "  SKIP (duplicate iter $idx): $base" >&2
        continue
    fi
    SEEN_IDX[$idx]=1
    cp "$f" "$BACKBONE_DIR/spec_${idx}.pdb"
done

# Backbone meta
cat > "$BACKBONE_DIR/meta.json" <<EOF
{
  "uuid8": "$UUID8",
  "generator": "boltzdesign",
  "length": $LENGTH,
  "num_backbones": $N_BACKBONES,
  "consumed_by_stacks": [
    "boltzdesign__ligandmpnn/$(basename "$MPNN_RUN_DIR")"
  ]
}
EOF

# ============================================================
# Stage 3: LigandMPNN -- one invocation per backbone
# ============================================================
echo ""
echo "=== Stage 3: LigandMPNN over $N_BACKBONES backbones ==="

shopt -s nullglob
for pdb_in in $(ls "$BACKBONE_DIR"/spec_*.pdb | sort -V); do
    base=$(basename "$pdb_in" .pdb)            # e.g. spec_42
    idx=${base#spec_}                          # 42
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

    cd "$MPNN_DESIGNS_DIR"

    # FASTA: one per backbone
    if [ -f "$mpnn_out/seqs/${base}.fa" ]; then
        cp "$mpnn_out/seqs/${base}.fa" "${MPNN_PREFIX}_design${idx}_sequences.fa"
    elif [ -f "$mpnn_out/seqs/spec.fa" ]; then
        cp "$mpnn_out/seqs/spec.fa" "${MPNN_PREFIX}_design${idx}_sequences.fa"
    fi

    # Repacked backbones
    for pdb in "$mpnn_out/backbones"/*.pdb; do
        bname=$(basename "$pdb" .pdb)
        seqidx=$(echo "$bname" | grep -oP '_\K\d+$' || true)
        [ -z "$seqidx" ] && continue
        s_idx=$((seqidx - 1))
        cp "$pdb" "${MPNN_PREFIX}_design${idx}_seq${s_idx}.pdb"
    done
done
shopt -u nullglob

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
# The _work sweep dir is SHARED when --num-runs>1 (dirname of the shared designs/).
# Purge only THIS job's uuid-scoped swept entries, then rmdir if empty, so concurrent
# jobs never race. Mirrors the restart-cleanup scoping above.
[ -d "$_WORK_SWEEP" ] && find "$_WORK_SWEEP" -maxdepth 1 -name "*_${UUID8}_*" -exec rm -rf {} + 2>/dev/null || true
rmdir "$_WORK_SWEEP" 2>/dev/null || true
# --num-runs>1: MPNN_RUN_DIR is per-uuid scaffolding (designs/logs/sbatch live
# in the shared <label>/ dir), so the purge above leaves it EMPTY. rmdir it —
# a no-op in single-run, where it holds designs/ + logs and rmdir fails.
rmdir "$MPNN_RUN_DIR" 2>/dev/null || true
# Tier B (optional): BD1 hallucination workdir (large)
purge "$WORK_BD1" || true                  # = $BACKBONE_DIR/bd1_workdir

# Phase 2: pack designs/ for inode reduction (single-run only; a shared
# --num-runs>1 designs/ must be packed once after all N jobs, not per-job here).
if [[ "$MPNN_DESIGNS_DIR" == "$MPNN_RUN_DIR/designs" ]]; then
    declare -F pack_unit >/dev/null && pack_unit "$MPNN_DESIGNS_DIR" || true
else
    # --num-runs>1: shared designs/, so pack only THIS run's *_<uuid8>_* slice
    # (race-safe) -> the summary dir becomes N per-uuid slice archives.
    declare -F pack_slice >/dev/null && pack_slice "$MPNN_DESIGNS_DIR" "$UUID8" || true
fi
