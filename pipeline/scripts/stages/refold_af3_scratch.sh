#!/bin/bash
# refold_af3.sh — refold all designs in a source run with AF3.
#
# Reads designs from <SOURCE_RUN>/designs/, builds per-design AF3 JSONs,
# loops AF3 inference, renames outputs to canonical scheme.
#
# Arguments (positional, required):
#   $1  REFOLD_DIR        the refold pass dir: <source_run>/refold_af3/<label>/
#   $2  SOURCE_RUN        source run dir (parent of designs/)
#   $3  RNA_SEQUENCE      e.g. ACUGCCGUAUAGGCAG (16nt canonical)
#
# Environment variables (set by sbatch wrapper):
#   AF3_SIF               path to alphafold3.sif
#   AF3_PARAMS            path to af3_params dir
#   MODEL_SEEDS_JSON      JSON array string like "[1]" or "[1,2]"

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 REFOLD_DIR SOURCE_RUN RNA_SEQUENCE" >&2
    exit 2
fi

REFOLD_DIR="$1"
SOURCE_RUN="$2"
RNA_SEQUENCE="$3"
DESIGNS_FILTER="${4:-}"   # per-uuid multif filter ("" = all designs)

: "${AF3_SIF:?AF3_SIF must be set}"
: "${AF3_PARAMS:?AF3_PARAMS must be set}"
: "${MODEL_SEEDS_JSON:?MODEL_SEEDS_JSON must be set}"

# [Fix 2] Already-finalized leaf (designs.tar.zst published): a requeue has nothing
# to do -- skip re-unpacking source / rebuilding inputs / recomputing inference.
if [[ -f "$REFOLD_DIR/designs.tar.zst" ]]; then
    echo "refold already complete (designs.tar.zst present) -- nothing to do"; exit 0
fi

# --- node-local scratch refactor: heavy I/O on $TMPDIR (node-local NVMe),
#     only the final packed tarball is copied back to Lustre $REFOLD_DIR. ---
: "${TMPDIR:?TMPDIR unset -- request --tmp=NG so SLURM provisions node-local /scratch}"
WORK="$TMPDIR/rf_${SLURM_JOB_ID:-$$}"; mkdir -p "$WORK"
trap 'rm -rf "$WORK" 2>/dev/null || true' EXIT
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
    _SRC_SCRATCH="$(mktemp -d "$WORK/refold_src.XXXXXX")"
    unpack_slice "$SOURCE_RUN/designs/$DESIGNS_FILTER.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH"
elif [[ ! -d "$DESIGNS_DIR" && -f "$SOURCE_RUN/designs.tar.zst" ]] && declare -F unpack_unit >/dev/null; then
    _SRC_SCRATCH="$(mktemp -d "$WORK/refold_src.XXXXXX")"   # /tmp is bound into the AF3 container; $TMPDIR (/scratch) is not
    unpack_unit "$SOURCE_RUN/designs.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH/designs"
fi
INPUTS_DIR="$WORK/inputs"

# Optionally read RNA MSA file content (if RNA_MSA_PATH is set + file exists)
RNA_MSA_CONTENT=""
if [ -n "${RNA_MSA_PATH:-}" ] && [ -f "$RNA_MSA_PATH" ]; then
    export RNA_MSA_CONTENT=$(cat "$RNA_MSA_PATH")
    echo "Using RNA MSA from: $RNA_MSA_PATH ($(wc -l < $RNA_MSA_PATH) lines)"
else
    echo "No RNA MSA (empty MSA mode)"
fi
OUTPUTS_DIR="$WORK/outputs"
CANONICAL_DIR="$WORK/designs"

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$CANONICAL_DIR"

# ---------------------------------------------------------------------------
# FlashAttention backend selection by GPU compute capability.
# AF3's default Triton FlashAttention kernel needs sm_80+ (Ampere/Ada). On the
# Turing 24G cards (sm_75: quadro_rtx_6000 / nvidia_titan_rtx, eu-g3) it aborts
# in ~1 min with:
#   ValueError: implementation='triton' for FlashAttention is unsupported on
#   this GPU generation. Please use implementation='xla' instead.
# DeepMind's portable fallback is the XLA attention path (the same one used for
# CUDA Capability 7.x devices like V100). The CC-7.x numerical-stability flag
# (XLA_FLAGS=--xla_disable_hlo_passes=custom-kernel-fusion-rewriter) is already
# set unconditionally on the apptainer run below, so the only missing piece is
# the backend flag. Auto-pick by compute_cap so one script Just Works on any
# card: CC < 80 -> xla; otherwise triton (no change for 3090/4090/a100).
# ---------------------------------------------------------------------------
FA_IMPL="triton"
_cc="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')"
if [[ "$_cc" =~ ^[0-9]+$ ]] && (( _cc < 80 )); then
    FA_IMPL="xla"
fi
echo "GPU compute_cap=${_cc:-unknown} -> flash_attention_implementation=$FA_IMPL"

echo "================================================================"
echo " AF3 refold worker"
echo " Source run:     $SOURCE_RUN"
echo " Refold dir:     $REFOLD_DIR"
echo " RNA sequence:   $RNA_SEQUENCE"
echo " Model seeds:    $MODEL_SEEDS_JSON"
echo " Started:        $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: extract sequence from each design + build AF3 JSON
# ============================================================
echo ""
echo "=== Stage 1: build per-design AF3 input JSONs ==="

# Use boltzgen.sif's Python (has biopython) — it's been our reliable conversion env
apptainer exec --bind /cluster/work/wicky --bind /tmp --bind "$TMPDIR" \
    /cluster/work/wicky/containers/boltzgen.sif \
    python3 << PYEOF
import json
import os
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder, CaPPBuilder
from Bio.SeqUtils import seq1

designs_dir = Path("$DESIGNS_DIR")
inputs_dir = Path("$INPUTS_DIR")
rna_seq = "$RNA_SEQUENCE"
model_seeds = $MODEL_SEEDS_JSON
# RNA MSA content is passed via env var (read from file by bash, exported)
rna_msa = os.environ.get("RNA_MSA_CONTENT", "")

# Aim to handle both .cif and .pdb. Walk designs/ and pick out files with sequences.
design_files = []
import os as _os
_allow = _os.environ.get("REFOLD_ALLOWLIST", "")
_allowset = {ln.strip() for ln in open(_allow)} if (_allow and _os.path.exists(_allow)) else set()
def _stem(nm):
    for e in (".cif.gz", ".cif", ".pdb"):
        if nm.endswith(e):
            return nm[:-len(e)]
    return nm
filt = "$DESIGNS_FILTER"   # per-uuid multif: keep only this run's designs ("" = all)
for f in sorted(designs_dir.iterdir()):
    if (f.suffix in [".cif", ".pdb"] or f.name.endswith(".cif.gz")) and f.is_file():
        # Skip metrics CSVs, sidecars, etc.
        if "sequences.fa" in f.name or f.suffix in [".npz", ".json", ".csv"]:
            continue
        if filt and filt not in f.name:
            continue
        if _allowset and _stem(f.name) not in _allowset:
            continue
        design_files.append(f)

print(f"  found {len(design_files)} design files in {designs_dir}")

def _collect_protein_sequence(structure):
    """Walk a Bio.PDB structure, return the longest protein chain as 1-letter seq."""
    best_seq = ""
    for model in structure:
        for chain in model:
            seq_letters = []
            for res in chain:
                resname = res.get_resname()
                if len(resname) == 3:
                    try:
                        one = seq1(resname)
                        if one == "X":
                            continue
                        seq_letters.append(one)
                    except Exception:
                        continue
            seq = "".join(seq_letters)
            if len(seq) > len(best_seq):
                best_seq = seq
    return best_seq


def extract_protein_sequence(path):
    """Return 1-letter sequence of the protein chain in the design file.
    Handles .cif, .pdb, and .cif.gz (transparently decompressed)."""
    import gzip
    if path.name.endswith(".cif.gz"):
        parser = MMCIFParser(QUIET=True)
        handle = gzip.open(str(path), "rt")
        try:
            structure = parser.get_structure("s", handle)
        finally:
            handle.close()
        return _collect_protein_sequence(structure)
    elif path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    elif path.suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        return None
    structure = parser.get_structure("s", str(path))
    return _collect_protein_sequence(structure)

count = 0
for path in design_files:
    seq = extract_protein_sequence(path)
    if not seq or len(seq) < 5:
        print(f"  SKIP {path.name}: no protein sequence extracted")
        continue
    # Build the AF3 JSON — name = the design file stem
    # Strip extension(s) — AF3 treats name as a name, not a filename
    if path.name.endswith(".cif.gz"):
        design_name = path.name[:-len(".cif.gz")]
    else:
        design_name = path.stem  # strips .cif or .pdb
    payload = {
        "name": design_name,
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": model_seeds,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": seq,
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [],
                }
            },
            {
                "rna": {
                    "id": ["B"],
                    "sequence": rna_seq,
                    "unpairedMsa": rna_msa,
                }
            }
        ],
        "bondedAtomPairs": None,
        "userCCD": None,
    }
    out_path = inputs_dir / f"{design_name}.json"
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    count += 1

print(f"  wrote {count} AF3 JSON inputs to {inputs_dir}")
PYEOF

N_INPUTS=$(ls "$INPUTS_DIR"/*.json 2>/dev/null | wc -l)
echo "  total JSON inputs: $N_INPUTS"
# [Fix 1] Persist expected source-design count so the aggregator/board can tell a
# leaf preempted/timed-out mid-finalize (partial refold) from a complete one.
printf '{"n_source_designs": %s}\n' "$N_INPUTS" > "$REFOLD_DIR/expected.json"

if [ "$N_INPUTS" -eq 0 ]; then
    echo "ERROR: no AF3 inputs were generated. Check $DESIGNS_DIR contents." >&2
    exit 1
fi

# ============================================================
# Stage 2: loop AF3 over each JSON
# ============================================================
echo ""
echo "=== Stage 2: AF3 BATCHED inference (one run over not-yet-done designs) ==="

# [batched 2026-06-22] Invoke run_alphafold.py ONCE over --input_dir so the model loads
# + JIT-compiles ONCE per bucket (all our 96-token inputs share bucket 256) instead of
# per-design -> ~10-13x faster (validated baseline-vs-batched, identical structures).
# Both XLA flags: --xla_gpu_enable_triton_gemm=false is DeepMind's compile-time fix that
# the old single-flag override dropped. Resume: only feed designs whose model.cif is
# missing (handles requeue + the old per-design layout); a preempt+requeue just
# re-batches the remainder (per-design model.cif is written incrementally).
TODO_DIR="$WORK/_todo_inputs"; rm -rf "$TODO_DIR"; mkdir -p "$TODO_DIR"
todo=0
for json_path in "$INPUTS_DIR"/*.json; do
    design_name=$(basename "$json_path" .json)
    af3low=$(echo "$design_name" | tr '[:upper:]' '[:lower:]')
    if find "$OUTPUTS_DIR/$design_name" "$OUTPUTS_DIR/$af3low" -name "model.cif" 2>/dev/null | grep -q .; then
        continue
    fi
    ln -sf "$json_path" "$TODO_DIR/$(basename "$json_path")"
    todo=$((todo + 1))
done
echo "  $todo of $N_INPUTS design(s) need folding (rest already done)"
if [ "$todo" -gt 0 ]; then
    apptainer run --nv \
        --bind /cluster/work/wicky:/cluster/work/wicky \
        --bind "$TMPDIR" \
        --env XLA_FLAGS="--xla_gpu_enable_triton_gemm=false --xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \
        "$AF3_SIF" \
        python /app/alphafold/run_alphafold.py \
            --input_dir="$TODO_DIR" \
            --model_dir="$AF3_PARAMS" \
            --output_dir="$OUTPUTS_DIR" \
            --run_data_pipeline=false \
            --flash_attention_implementation="$FA_IMPL" \
            --num_diffusion_samples="${AF3_NUM_DIFFUSION_SAMPLES:-1}"
fi
rm -rf "$TODO_DIR"

# ============================================================
# Stage 3: rename outputs to canonical scheme
# ============================================================
echo ""
echo "=== Stage 3: rename to canonical ==="

cd "$CANONICAL_DIR"
n_canonical=0

# [batched 2026-06-22] Layout-agnostic: model.cif lives at <af3low>/seed-S_sample-N/ in
# the batched layout and <name>/<af3low>/seed-S_sample-N/ in the old per-design layout.
# In both, the seed dir's PARENT is the lowercased AF3 name -> map back to original case
# via the input JSONs so canonical names match B2/ESM (case-sensitive downstream).
declare -A _orig
for _j in "$INPUTS_DIR"/*.json; do
    _n=$(basename "$_j" .json); _l=$(echo "$_n" | tr '[:upper:]' '[:lower:]'); _orig["$_l"]="$_n"
done
while IFS= read -r src_cif; do
    [ -f "$src_cif" ] || continue
    sample_dir=$(dirname "$src_cif")
    af3low=$(basename "$(dirname "$sample_dir")")
    design_name="${_orig[$af3low]:-$af3low}"
    sub=$(basename "$sample_dir")
    seed=$(echo "$sub" | sed -n 's/^seed-\([0-9]*\)_sample-.*/\1/p')
    sample=$(echo "$sub" | sed -n 's/^seed-[0-9]*_sample-\([0-9]*\)/\1/p')
    cp "$src_cif" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}.cif"
    src_side="$sample_dir/summary_confidences.json"
    [ -f "$src_side" ] && cp "$src_side" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}_summary_confidences.json"
    # Retain the FULL confidences.json transiently as *_full_data.json so IPSAE can
    # read the complete PAE matrix in finalize; it is deleted before packing (only
    # the tiny *_ipsae.json sidecar survives).
    src_full="$sample_dir/confidences.json"
    [ -f "$src_full" ] && cp "$src_full" "$CANONICAL_DIR/${design_name}_AF3_seed${seed}_sample${sample}_full_data.json"
    n_canonical=$((n_canonical + 1))
done < <(find "$OUTPUTS_DIR" -name "model.cif" -path "*seed-*_sample-*" 2>/dev/null)
# Top-ranked: <af3low>/<af3low>_model.cif (not under a seed dir)
while IFS= read -r top_cif; do
    [ -f "$top_cif" ] || continue
    af3low=$(basename "$(dirname "$top_cif")")
    design_name="${_orig[$af3low]:-$af3low}"
    cp "$top_cif" "$CANONICAL_DIR/${design_name}_AF3_top.cif"
done < <(find "$OUTPUTS_DIR" -name "*_model.cif" ! -path "*seed-*" 2>/dev/null)

# --- IPSAE interface scoring (non-fatal add-on) -----------------------------
# Score each canonical AF3 model against its retained full PAE, writing a tiny
# <stem>_ipsae.json next to each .cif, then drop the bulky *_full_data.json so
# only the sidecars get packed. Runs on the host python (numpy); never aborts.
IPSAE_DRIVER="$(dirname "$0")/../ipsae_sidecar.py"
if command -v python3 >/dev/null 2>&1 && [ -f "$IPSAE_DRIVER" ]; then
    python3 "$IPSAE_DRIVER" --mode af3 --designs-dir "$CANONICAL_DIR" \
        || echo "WARN: IPSAE scoring failed (non-fatal), continuing"
fi
rm -f "$CANONICAL_DIR"/*_full_data.json

# [Fix 2] Pack designs/ FIRST so the leaf is finalized (designs.tar.zst published)
# BEFORE outputs/ is deleted: a requeue after this point hits the top early-exit
# instead of recomputing all inference (the per-design SKIP guard reads outputs/).
declare -F pack_unit >/dev/null && pack_unit "$CANONICAL_DIR" || true
# node-local scratch auto-cleaned; publish only the single packed tarball to Lustre
cp -f "$WORK"/designs.tar.zst* "$REFOLD_DIR/"
echo "published $(du -h "$REFOLD_DIR/designs.tar.zst" 2>/dev/null | cut -f1) tarball to Lustre (inode cost: 1)"

# Raw AF3 outputs/ tree is redundant after flattening to designs/
# (analysis + DockQ read designs/ only); drop it to reclaim disk.
echo "Removing raw outputs/: $OUTPUTS_DIR"
rm -rf "$OUTPUTS_DIR"
declare -F purge >/dev/null && purge "$REFOLD_DIR/inputs" || rm -rf "$REFOLD_DIR/inputs"
[[ -n "${_SRC_SCRATCH:-}" ]] && rm -rf "$_SRC_SCRATCH" || true   # drop the unpacked source-designs scratch

echo ""
echo "================================================================"
echo " Done at $(date)"
echo " Canonical .cif files written: $n_canonical (+ top-ranked per design)"
echo " Refold dir: $REFOLD_DIR"
echo " Canonical: $CANONICAL_DIR"
echo "================================================================"
