#!/bin/bash
# bif_cysregen.sbatch — re-run BoltzIF inverse-folding (Route-B, --inverse_fold_avoid C)
# on the SURVIVING backbones that produced BG+BIF cysteine designs in midrun_v2, to
# regenerate them cysteine-free. One array task per source uuid (manifest line).
#
# Mechanic (validated 2026-07-02, probe job 5541265):
#   BoltzGen `--steps inverse_folding folding` resumes from a hand-staged
#   bg_output/intermediate_designs/ containing the dirty spec_<idx>.cif backbones +
#   the spec-constant design-metadata sidecar spec_<idx>.npz (identical across all
#   backbones of this 80-protein+RNA spec -> one saved reference serves all). The
#   Route-B writer patch gives 10 seqs/backbone; --inverse_fold_avoid C drops Cys.
#   Output design index preserves the source backbone index (spec_<bb>_<seq>).
#
# Idempotent: a uuid whose slice already exists in the cysregen pool is skipped, so a
# requeue/re-run continues where it stopped. Writes ONLY to the new cysregen pool +
# node-local-ish _cache scratch; never touches midrun_v2 or the _backbones archive.
#
# Submit (after choosing --array range): sbatch --array=0-165 bif_cysregen.sbatch
#SBATCH --job-name=bif_cysregen
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --partition=gpuhe.24h
#SBATCH --qos=es_wicky/gpuhe/24
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --gpus=rtx_4090:1
#SBATCH --chdir=/cluster/work/wicky/lvalsecchi/rna_binder/ops/submit/cys_regen_bgbif_2026-07-02
#SBATCH --output=array_logs/%A_%a.log
#SBATCH --error=array_logs/%A_%a.err
set -euo pipefail

RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
CODE=$RB/code
SUBDIR=$RB/ops/submit/cys_regen_bgbif_2026-07-02
MANIFEST="${MANIFEST:-$SUBDIR/bg_bif_dirty_manifest.tsv}"
REF_NPZ="${REF_NPZ:-$SUBDIR/ref_design_meta.npz}"       # spec-constant design metadata
OUT_LABEL="${OUT_LABEL:-cysregen_2026-07-02}"
LENGTH=80

# --- env replicated verbatim from a production BIF run.sbatch -----------------
export BG_SIF=/cluster/work/wicky/containers/boltzgen.sif
export BG_HOME=$RB/_cache/boltzgen_sif/home
export BG_BG_ENVS="--env HF_HOME=$RB/_cache/boltzgen_sif/hf --env HUGGINGFACE_HUB_CACHE=$RB/_cache/boltzgen_sif/hf --env HF_HUB_OFFLINE=1 --env TRANSFORMERS_OFFLINE=1 --env XDG_CACHE_HOME=$RB/_cache/boltzgen_sif/xdg --env MPLCONFIGDIR=$RB/_cache/boltzgen_sif/mpl --env TMPDIR=$RB/_cache/boltzgen_sif/tmp"
export TARGET_RNA_CIF=$RB/targets/csy4_crrna_free/csy4_rna_only.cif   # .pdb cleaned; convert per-task
WRITER_PATCH=$RB/boltzgen_patch/writer.py

# --- pick this task's manifest line (line 1 = header) -------------------------
TASK="${SLURM_ARRAY_TASK_ID:?must run as a job array}"
LINE=$(sed -n "$((TASK + 2))p" "$MANIFEST")
[ -n "$LINE" ] || { echo "no manifest line for task $TASK"; exit 0; }
UUID=$(echo "$LINE" | cut -f1)
SRC_BB=$(echo "$LINE" | cut -f2)
SPEC_CSV=$(echo "$LINE" | cut -f3)
IFS=',' read -r -a SPECS <<< "$SPEC_CSV"

BIF_DESIGNS_DIR="$CAMP/boltzgen__boltzif/$OUT_LABEL/designs"
BIF_PREFIX="csy4free_BG_BIF_len${LENGTH}_${UUID}"
mkdir -p "$BIF_DESIGNS_DIR"

echo "=== bif_cysregen task=$TASK uuid=$UUID backbones=${#SPECS[@]} label=$OUT_LABEL ==="
# Idempotency: slice already built -> nothing to do.
if [ -f "$BIF_DESIGNS_DIR/${UUID}.tar.zst" ]; then
    echo "slice ${UUID}.tar.zst already exists -> skip"; exit 0
fi
nvidia-smi -L || true

WORK=$(mktemp -d "$RB/_cache/bif_cysregen_${UUID}_XXXXXX")
cleanup(){ rm -rf "$WORK" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$WORK/bg_output/intermediate_designs"

# 1. stage dirty backbones + the spec-constant design-metadata npz
for s in "${SPECS[@]}"; do
    [ -f "$SRC_BB/$s" ] || { echo "  WARN missing $SRC_BB/$s"; continue; }
    cp "$SRC_BB/$s" "$WORK/bg_output/intermediate_designs/$s"
    cp "$REF_NPZ"   "$WORK/bg_output/intermediate_designs/${s%.cif}.npz"
done
N_STAGED=$(ls "$WORK/bg_output/intermediate_designs"/*.cif 2>/dev/null | wc -l)
echo "staged $N_STAGED backbones"
[ "$N_STAGED" -gt 0 ] || { echo "ERROR: nothing staged"; exit 1; }

# 2. spec.yaml + target pdb (regenerate .pdb from surviving .cif)
cd "$WORK"
apptainer exec --bind /cluster/work/wicky --bind /tmp "$BG_SIF" python3 -c "
from Bio.PDB import MMCIFParser, PDBIO
s = MMCIFParser(QUIET=True).get_structure('rna', '$TARGET_RNA_CIF')
io = PDBIO(); io.set_structure(s); io.save('$WORK/csy4_rna_only.pdb')
"
cat > spec.yaml <<'YAML'
entities:
- protein:
    id: D
    sequence: 80
- file:
    include:
    - chain:
        id: B
    path: csy4_rna_only.pdb
YAML

# 3. Stage-2 BoltzIF: inverse_folding + folding, avoid Cys, Route-B writer patch
apptainer exec --nv \
    --home "$BG_HOME" \
    --bind /cluster/work/wicky \
    --bind "$WRITER_PATCH:/app/src/boltzgen/task/predict/writer.py" \
    $BG_BG_ENVS \
    "$BG_SIF" \
    boltzgen run spec.yaml \
        --output "$WORK/bg_output" \
        --steps inverse_folding folding \
        --inverse_fold_num_sequences 10 \
        --inverse_fold_avoid C

# 4. harvest ranked BIF outputs -> canonical names in the cysregen designs pool
if [ -d "$WORK/bg_output/intermediate_designs_inverse_folded/refold_cif" ]; then
    BIF_NATIVE="$WORK/bg_output/intermediate_designs_inverse_folded/refold_cif"
else
    BIF_NATIVE="$WORK/bg_output/intermediate_designs/refold_cif"
fi
[ -d "$BIF_NATIVE" ] || { echo "ERROR: no BIF native output at $BIF_NATIVE"; exit 1; }

cd "$BIF_DESIGNS_DIR"
shopt -s nullglob
for f in $(ls "$BIF_NATIVE"/spec_*.cif | sort -V); do
    base=$(basename "$f" .cif)               # spec_<bb>_<seq> (patched writer)
    if [[ "$base" =~ spec_([0-9]+)_([0-9]+)$ ]]; then
        d_idx=$((10#${BASH_REMATCH[1]})); s_idx=$((10#${BASH_REMATCH[2]}))
        new="${BIF_PREFIX}_design${d_idx}_seq${s_idx}.cif"
    elif [[ "$base" =~ spec_([0-9]+)$ ]]; then
        d_idx=$((10#${BASH_REMATCH[1]}))
        new="${BIF_PREFIX}_design${d_idx}.cif"
    else
        echo "  skip unparseable $base"; continue
    fi
    cp "$f" "$new"
    [ -f "${f%.cif}.npz" ] && cp "${f%.cif}.npz" "${new%.cif}.npz"
done
shopt -u nullglob

N_BIF=$(find "$BIF_DESIGNS_DIR" -maxdepth 1 -name "*_${UUID}_*.cif" 2>/dev/null | wc -l)
echo "harvested $N_BIF canonical BIF cifs for uuid $UUID"
[ "$N_BIF" -gt 0 ] || { echo "ERROR: 0 harvested"; exit 1; }

# 5. pack into a per-uuid slice (flock'd, race-safe) — same helper the pipeline uses
source "$CODE/scripts/stages/archive_pipeline.sh"
pack_slice "$BIF_DESIGNS_DIR" "$UUID"

cat > "$BIF_DESIGNS_DIR/${UUID}.cysregen.json" <<EOF
{ "uuid8": "$UUID", "kind": "bif_cysregen_reinverse_fold",
  "src_backbones": "$SRC_BB", "dirty_backbones_used": ${#SPECS[@]},
  "raw_designs": $N_BIF, "inverse_fold_avoid": "C",
  "note": "BoltzIF re-inverse-fold (avoid C) on the surviving dirty midrun_v2 backbones; cysteine-free replacement for the contaminated BG+BIF designs of this uuid." }
EOF

echo "=== DONE uuid=$UUID  slice=$BIF_DESIGNS_DIR/${UUID}.tar.zst ($N_BIF designs) ==="
