#!/bin/bash
# PROBE: can BoltzGen `--steps inverse_folding folding` resume from a HAND-STAGED
# intermediate_designs/ (backbones only, no prior `--steps design` in this output
# dir)? Stages 2 dirty backbones from an existing BG uuid, runs BIF with the
# Route-B writer patch + --inverse_fold_avoid C, and reports output naming + a
# cysteine grep. Writes only to a scratch dir; touches no production pool.
#SBATCH --job-name=bif_probe
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --partition=gpuhe.4h
#SBATCH --qos=es_wicky/gpuhe/4
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --gpus=rtx_4090:1
#SBATCH --chdir=/cluster/work/wicky/lvalsecchi/rna_binder/ops/submit/cys_regen_bgbif_2026-07-02
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
set -euo pipefail

RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
SRC_UUID=0000df81
SRC_BB="$CAMP/_backbones/BG_len80_${SRC_UUID}"
SPECS=(spec_02.cif spec_03.cif)          # 2 known-dirty backbones

# --- env replicated verbatim from a real BIF run.sbatch -----------------------
export BG_SIF=/cluster/work/wicky/containers/boltzgen.sif
export BG_HOME=$RB/_cache/boltzgen_sif/home
export BG_BG_ENVS="--env HF_HOME=$RB/_cache/boltzgen_sif/hf --env HUGGINGFACE_HUB_CACHE=$RB/_cache/boltzgen_sif/hf --env HF_HUB_OFFLINE=1 --env TRANSFORMERS_OFFLINE=1 --env XDG_CACHE_HOME=$RB/_cache/boltzgen_sif/xdg --env MPLCONFIGDIR=$RB/_cache/boltzgen_sif/mpl --env TMPDIR=$RB/_cache/boltzgen_sif/tmp"
export TARGET_RNA_CIF=$RB/targets/csy4_crrna_free/csy4_rna_only.cif   # .pdb was cleaned; convert below
WRITER_PATCH=$RB/boltzgen_patch/writer.py

WORK=$(mktemp -d "$RB/_cache/bif_probe_XXXXXX")
echo "=== PROBE workdir: $WORK ==="
mkdir -p "$WORK/bg_output/intermediate_designs"

# 1. hand-stage the backbones + the spec-constant design-metadata npz that BIF's
#    DataLoader requires (spec_<idx>.npz). The npz is IDENTICAL across all backbones
#    of this spec (design_mask/mol_type/ss_type/... depend only on the 80-protein+RNA
#    token layout), so one saved reference serves every backbone.
REF_NPZ=/cluster/work/wicky/lvalsecchi/rna_binder/ops/submit/cys_regen_bgbif_2026-07-02/ref_design_meta.npz
for s in "${SPECS[@]}"; do
    cp "$SRC_BB/$s" "$WORK/bg_output/intermediate_designs/$s"
    cp "$REF_NPZ"   "$WORK/bg_output/intermediate_designs/${s%.cif}.npz"
done
echo "staged backbones + npz:"; ls -1 "$WORK/bg_output/intermediate_designs"

# 2. spec.yaml + target, exactly as boltzgen_multif.sh Stage 1 sets up.
#    The .pdb target was cleaned, so regenerate it from the surviving .cif with the
#    same MMCIFParser->PDBIO conversion the pipeline uses for backbones.
cd "$WORK"
apptainer exec --bind /cluster/work/wicky --bind /tmp "$BG_SIF" python3 -c "
from Bio.PDB import MMCIFParser, PDBIO
s = MMCIFParser(QUIET=True).get_structure('rna', '$TARGET_RNA_CIF')
io = PDBIO(); io.set_structure(s); io.save('$WORK/csy4_rna_only.pdb')
print('converted target -> csy4_rna_only.pdb')
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

# 3. Stage-2 only: inverse_folding + folding, avoid Cys, Route-B writer patch
echo "=== running BoltzGen --steps inverse_folding folding (avoid C) ==="
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
        --inverse_fold_avoid C || { echo "BIF RUN FAILED (rc=$?)"; }

# 4. report: where did outputs land, how are they named, any cysteine?
echo "=== bg_output tree (2 levels) ==="; find "$WORK/bg_output" -maxdepth 2 -type d | sort
for cand in intermediate_designs_inverse_folded/refold_cif intermediate_designs/refold_cif; do
    d="$WORK/bg_output/$cand"
    if [ -d "$d" ]; then
        echo "=== OUTPUT DIR: $cand ==="
        ls -1 "$d" | head -30
        n=$(ls "$d"/*.cif 2>/dev/null | wc -l)
        cys=$(grep -l -w CYS "$d"/*.cif 2>/dev/null | wc -l)
        echo ">>> cif files: $n   files-containing-CYS: $cys"
    fi
done
echo "=== PROBE DONE (scratch left at $WORK for inspection) ==="
