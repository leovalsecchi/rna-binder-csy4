#!/bin/bash
# refold_b2.sh — refold all designs in a source run with Boltz-2.
#
# Reads designs from <SOURCE_RUN>/designs/, builds per-design YAMLs,
# runs `boltz predict` on the whole directory, renames outputs to canonical.
#
# Arguments (positional):
#   $1  REFOLD_DIR        <source_run>/refold_b2/<label>/
#   $2  SOURCE_RUN        source run dir
#   $3  RNA_SEQUENCE      e.g. ACUGCCGUAUAGGCAG (16nt canonical)
#
# Environment variables (set by sbatch wrapper):
#   B2_SIF                path to boltz2.sif
#   BOLTZ_CACHE           path to /cluster/work/wicky/params/.boltz
#   B2_DIFFUSION_SAMPLES  e.g. 5
#   B2_OUTPUT_FORMAT      e.g. mmcif
#   B2_WRITE_FULL_PAE     "true" or "false"

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 REFOLD_DIR SOURCE_RUN RNA_SEQUENCE" >&2
    exit 2
fi

REFOLD_DIR="$1"
SOURCE_RUN="$2"
RNA_SEQUENCE="$3"
DESIGNS_FILTER="${4:-}"   # per-uuid multif filter ("" = all designs)

: "${B2_SIF:?B2_SIF must be set}"
: "${BOLTZ_CACHE:?BOLTZ_CACHE must be set}"
: "${B2_DIFFUSION_SAMPLES:?B2_DIFFUSION_SAMPLES must be set}"
: "${B2_OUTPUT_FORMAT:?B2_OUTPUT_FORMAT must be set}"
B2_WRITE_FULL_PAE="${B2_WRITE_FULL_PAE:-false}"

# [Fix 2] Already-finalized leaf (designs.tar.zst published): a requeue has nothing
# to do -- skip re-unpacking source / rebuilding inputs / recomputing inference.
if [[ -f "$REFOLD_DIR/designs.tar.zst" ]]; then
    echo "refold already complete (designs.tar.zst present) -- nothing to do"; exit 0
fi

DESIGNS_DIR="$SOURCE_RUN/designs"
# Phase 2: if the source designs are archived (.tar.zst), expand to node-local
# scratch and read from there (host zstd; unpack_unit is flock'd + idempotent).
source "$(dirname "$0")/archive_pipeline.sh" 2>/dev/null || true

# [resume 2026-06-16] (re)started job continues where it left off: clear stale caches
# (avoid boltz/lightning "Missing folder") + restore this leaf's packed partial outputs/.
rm -rf "$REFOLD_DIR"/.cache "$REFOLD_DIR"/.mpl_cache "$REFOLD_DIR"/.numba_cache "$REFOLD_DIR"/.triton_cache "$REFOLD_DIR"/.cueq_triton_cache 2>/dev/null || true
if [[ -f "$REFOLD_DIR/outputs.tar.zst" && ! -d "$REFOLD_DIR/outputs" ]] && declare -F unpack_unit >/dev/null; then
    unpack_unit "$REFOLD_DIR/outputs.tar.zst" || true
fi
if [[ -n "$DESIGNS_FILTER" && -f "$SOURCE_RUN/designs/$DESIGNS_FILTER.tar.zst" ]] && declare -F unpack_slice >/dev/null; then
    # --num-runs>1: this run's designs are a per-uuid slice archive; expand just it.
    _SRC_SCRATCH="$(mktemp -d /tmp/refold_src.XXXXXX)"
    unpack_slice "$SOURCE_RUN/designs/$DESIGNS_FILTER.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH"
elif [[ ! -d "$DESIGNS_DIR" && -f "$SOURCE_RUN/designs.tar.zst" ]] && declare -F unpack_unit >/dev/null; then
    _SRC_SCRATCH="$(mktemp -d /tmp/refold_src.XXXXXX)"   # /tmp is bound into the Boltz container; $TMPDIR (/scratch) is not
    unpack_unit "$SOURCE_RUN/designs.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH/designs"
fi
INPUTS_DIR="$REFOLD_DIR/inputs"
# Resolve REFOLD_DIR to absolute path so subsequent paths are stable across cd
REFOLD_DIR="$(cd "$REFOLD_DIR" && pwd)"
OUTPUTS_DIR="$REFOLD_DIR/outputs"
CANONICAL_DIR="$REFOLD_DIR/designs"

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$CANONICAL_DIR"

echo "================================================================"
echo " Boltz-2 refold worker"
echo " Source run:     $SOURCE_RUN"
echo " Refold dir:     $REFOLD_DIR"
echo " RNA sequence:   $RNA_SEQUENCE"
echo " Diffusion N:    $B2_DIFFUSION_SAMPLES"
echo " Started:        $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: extract sequences + build per-design B2 YAMLs
# ============================================================
echo ""
echo "=== Stage 1: build per-design B2 input YAMLs ==="

apptainer exec --bind /cluster/work/wicky --bind /tmp \
    /cluster/work/wicky/containers/boltzgen.sif \
    python3 << PYEOF
import os
import gzip
from pathlib import Path
import yaml
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

designs_dir = Path("$DESIGNS_DIR")
inputs_dir = Path("$INPUTS_DIR")
rna_seq = "$RNA_SEQUENCE"

design_files = []
filt = "$DESIGNS_FILTER"   # per-uuid multif: keep only this run's designs ("" = all)
for f in sorted(designs_dir.iterdir()):
    if (f.suffix in [".cif", ".pdb"] or f.name.endswith(".cif.gz")) and f.is_file():
        if "sequences.fa" in f.name or f.suffix in [".npz", ".json", ".csv"]:
            continue
        if filt and filt not in f.name:
            continue
        design_files.append(f)

print(f"  found {len(design_files)} design files in {designs_dir}")


def collect_protein_seq(structure):
    best = ""
    for model in structure:
        for chain in model:
            letters = []
            for res in chain:
                resname = res.get_resname()
                if len(resname) == 3:
                    try:
                        one = seq1(resname)
                        if one == "X":
                            continue
                        letters.append(one)
                    except Exception:
                        continue
            seq = "".join(letters)
            if len(seq) > len(best):
                best = seq
    return best


def extract_seq(path):
    if path.name.endswith(".cif.gz"):
        parser = MMCIFParser(QUIET=True)
        with gzip.open(str(path), "rt") as handle:
            structure = parser.get_structure("s", handle)
        return collect_protein_seq(structure)
    if path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    elif path.suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        return None
    structure = parser.get_structure("s", str(path))
    return collect_protein_seq(structure)


count = 0
for path in design_files:
    seq = extract_seq(path)
    if not seq or len(seq) < 5:
        print(f"  SKIP {path.name}: no protein sequence extracted")
        continue
    # Strip .cif.gz / .cif / .pdb to get design stem
    name = path.name
    for suffix in [".cif.gz", ".cif", ".pdb"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    # RNA chain: attach the RNA MSA (a3m) if provided, else single-sequence.
    rna_entry = {"id": "B", "sequence": rna_seq}
    rna_msa_path = os.environ.get("RNA_MSA_PATH", "")
    if rna_msa_path and os.path.exists(rna_msa_path):
        rna_entry["msa"] = rna_msa_path
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": seq, "msa": "empty"}},
            {"rna":     rna_entry},
        ],
    }
    # If a template CIF is provided, add it (applies to chain B by default)
    tpl_path = os.environ.get("RNA_TEMPLATE_PATH", "")
    tpl_chain = os.environ.get("RNA_TEMPLATE_CHAIN", "B")
    if tpl_path and os.path.exists(tpl_path):
        payload["templates"] = [{"cif": tpl_path, "ids": [tpl_chain]}]
    out_path = inputs_dir / f"{name}.yaml"
    with open(out_path, "w") as f:
        yaml.dump(payload, f, default_flow_style=False, sort_keys=False)
    count += 1

print(f"  wrote {count} B2 YAML inputs to {inputs_dir}")
PYEOF

N_INPUTS=$(ls "$INPUTS_DIR"/*.yaml 2>/dev/null | wc -l)
echo "  total YAML inputs: $N_INPUTS"
# [Fix 1] Persist expected source-design count so the aggregator/board can tell a
# leaf preempted/timed-out mid-finalize (partial refold) from a complete one.
printf '{"n_source_designs": %s}\n' "$N_INPUTS" > "$REFOLD_DIR/expected.json"

if [ "$N_INPUTS" -eq 0 ]; then
    echo "ERROR: no B2 inputs were generated. Check $DESIGNS_DIR contents." >&2
    exit 1
fi

# ============================================================
# Stage 2: B2 predict on the inputs directory
# ============================================================
echo ""
echo "=== Stage 2: Boltz-2 inference ==="

WRITE_FULL_PAE_FLAG=""
if [ "$B2_WRITE_FULL_PAE" = "true" ]; then
    WRITE_FULL_PAE_FLAG="--write_full_pae"
fi

apptainer run --nv \
    --bind /cluster/work/wicky:/cluster/work/wicky \
    --env BOLTZ_CACHE="$BOLTZ_CACHE" \
    --env MPLCONFIGDIR="$REFOLD_DIR/.mpl_cache" \
    --env NUMBA_CACHE_DIR="$REFOLD_DIR/.numba_cache" \
    --env TRITON_CACHE_DIR="$REFOLD_DIR/.triton_cache" \
    --env CUEQ_TRITON_CACHE_DIR="$REFOLD_DIR/.cueq_triton_cache" \
    --env XDG_CACHE_HOME="$REFOLD_DIR/.cache" \
    "$B2_SIF" \
    boltz predict "$INPUTS_DIR" \
        --out_dir "$OUTPUTS_DIR" \
        --cache "$BOLTZ_CACHE" \
        --output_format "$B2_OUTPUT_FORMAT" \
        --diffusion_samples "$B2_DIFFUSION_SAMPLES" \
        $WRITE_FULL_PAE_FLAG

# ============================================================
# Stage 3: rename outputs to canonical
# ============================================================
echo ""
echo "=== Stage 3: rename to canonical ==="

# B2 writes outputs to: <OUTPUTS_DIR>/boltz_results_<inputs_dirname>/predictions/<design_name>/
# Each design dir contains: <design_name>_model_0.cif, _model_1.cif, etc. (one per diffusion sample)

cd "$CANONICAL_DIR"
n_canonical=0

shopt -s nullglob
# Find the predictions root
B2_RESULTS_ROOT=$(find "$OUTPUTS_DIR" -maxdepth 1 -type d -name "boltz_results_*" | head -1)
if [ -z "$B2_RESULTS_ROOT" ]; then
    echo "WARN: no boltz_results_* root found under $OUTPUTS_DIR"
    PREDICTIONS_ROOT=""
else
    PREDICTIONS_ROOT="$B2_RESULTS_ROOT/predictions"
fi

if [ -n "$PREDICTIONS_ROOT" ] && [ -d "$PREDICTIONS_ROOT" ]; then
    for design_dir in "$PREDICTIONS_ROOT"/*; do
        [ -d "$design_dir" ] || continue
        design_name=$(basename "$design_dir")
        for cif in "$design_dir"/*.cif; do
            [ -f "$cif" ] || continue
            cif_name=$(basename "$cif" .cif)
            # B2 names like "<design_name>_model_0", "<design_name>_model_1", ...
            sample=$(echo "$cif_name" | grep -oP 'model_\K\d+$' || true)
            [ -z "$sample" ] && sample=0
            new="${design_name}_B2_sample${sample}.cif"
            cp "$cif" "$CANONICAL_DIR/$new"
            n_canonical=$((n_canonical + 1))
        done
        # Copy confidence sidecars if present, renaming to match the canonical _B2_sample<N>
        for json in "$design_dir"/*confidence*.json; do
            [ -f "$json" ] || continue
            json_base=$(basename "$json")
            # Extract sample index from native name like "confidence_<design>_model_<N>.json"
            sample=$(echo "$json_base" | grep -oP 'model_\K\d+' | head -1)
            [ -z "$sample" ] && sample=0
            cp "$json" "$CANONICAL_DIR/${design_name}_B2_sample${sample}_confidence.json"
        done
    done
fi
shopt -u nullglob

# [Fix 2] Pack designs/ FIRST so the leaf is finalized (designs.tar.zst published)
# BEFORE outputs/ is deleted: a requeue after this point hits the top early-exit
# instead of recomputing all inference (the per-design SKIP guard reads outputs/).

# --- IPSAE interface scoring (non-fatal add-on) -----------------------------
# Score each canonical B2 model against Boltz's full PAE (pae_*.npz still under
# outputs/), writing a tiny <stem>_ipsae.json per model. Needs B2_WRITE_FULL_PAE
# (the launcher default = true). Runs on the host python (numpy); never aborts.
IPSAE_DRIVER="$(dirname "$0")/../ipsae_sidecar.py"
if command -v python3 >/dev/null 2>&1 && [ -f "$IPSAE_DRIVER" ]; then
    python3 "$IPSAE_DRIVER" --mode b2 --designs-dir "$CANONICAL_DIR" --raw-dir "$OUTPUTS_DIR" \
        || echo "WARN: IPSAE scoring failed (non-fatal), continuing"
fi

declare -F pack_unit >/dev/null && pack_unit "$REFOLD_DIR/designs" || true

# Raw Boltz outputs/ tree is redundant after flattening to designs/
# (analysis + DockQ read designs/ only); drop it to reclaim disk.
echo "Removing raw outputs/: $OUTPUTS_DIR"
rm -rf "$OUTPUTS_DIR"
declare -F purge >/dev/null && purge "$REFOLD_DIR/inputs" || rm -rf "$REFOLD_DIR/inputs"
[[ -n "${_SRC_SCRATCH:-}" ]] && rm -rf "$_SRC_SCRATCH" || true   # drop the unpacked source-designs scratch

echo ""
echo "================================================================"
echo " Done at $(date)"
echo " Canonical .cif files written: $n_canonical"
echo " Refold dir: $REFOLD_DIR"
echo " Canonical: $CANONICAL_DIR"
echo "================================================================"
