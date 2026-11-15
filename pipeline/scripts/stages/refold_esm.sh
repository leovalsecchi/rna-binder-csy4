#!/bin/bash
# refold_esm.sh — refold all designs in a source run with ESMFold2 (local biohub/ESMFold2).
#
# Mirrors refold_b2.sh: extract the protein sequence from each design, build a
# per-design Boltz-style YAML (protein A + rna B), then run ESMFold2 ONCE over the
# whole input dir (model loaded a single time) writing canonical
# <design>_ESM_sample0.cif. Single-sequence (msaoff).
#
# Arguments (positional):
#   $1  REFOLD_DIR        <source_run>/refold_esm/<label>/
#   $2  SOURCE_RUN        source run dir
#   $3  RNA_SEQUENCE      e.g. ACUGCCGUAUAGGCAG (16nt canonical)
#   $4  DESIGNS_FILTER    per-uuid multif filter ("" = all designs)
#
# Environment variables (set by the sbatch wrapper):
#   ESM_SIF                  path to esmc.sif
#   ESM_HF_CACHE             HF cache holding biohub/ESMFold2 + ESMC-6B weights
#   ESM_NUM_LOOPS            ESMFold2 folding loops (default 20)
#   ESM_NUM_SAMPLING_STEPS   diffusion sampling steps (default 100)
#   ESM_NUM_DIFFUSION_SAMPLES samples per design (default 1)
#   ESM_SEED                 RNG seed (default 0)

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 REFOLD_DIR SOURCE_RUN RNA_SEQUENCE [DESIGNS_FILTER]" >&2
    exit 2
fi

REFOLD_DIR="$1"
SOURCE_RUN="$2"
RNA_SEQUENCE="$3"
DESIGNS_FILTER="${4:-}"

: "${ESM_SIF:?ESM_SIF must be set}"
: "${ESM_HF_CACHE:?ESM_HF_CACHE must be set}"
ESM_NUM_LOOPS="${ESM_NUM_LOOPS:-20}"
ESM_NUM_SAMPLING_STEPS="${ESM_NUM_SAMPLING_STEPS:-100}"
ESM_NUM_DIFFUSION_SAMPLES="${ESM_NUM_DIFFUSION_SAMPLES:-1}"
ESM_SEED="${ESM_SEED:-0}"

# [Fix 2] Already-finalized leaf (designs.tar.zst published): a requeue has nothing
# to do -- skip re-unpacking source / rebuilding inputs / recomputing inference.
if [[ -f "$REFOLD_DIR/designs.tar.zst" ]]; then
    echo "refold already complete (designs.tar.zst present) -- nothing to do"; exit 0
fi

DESIGNS_DIR="$SOURCE_RUN/designs"
# Phase 2: if the source designs are archived (.tar.zst), expand to node-local
# scratch and read from there (host zstd; unpack_unit is flock'd + idempotent).
source "$(dirname "$0")/archive_pipeline.sh" 2>/dev/null || true

# [resume] (re)started job continues where it left off: clear stale caches, and
# the per-design SKIP guard in esmfold2_infer.py keeps already-folded designs.
rm -rf "$REFOLD_DIR"/.cache "$REFOLD_DIR"/.mpl_cache "$REFOLD_DIR"/.numba_cache "$REFOLD_DIR"/.triton_cache 2>/dev/null || true
if [[ -n "$DESIGNS_FILTER" && -f "$SOURCE_RUN/designs/$DESIGNS_FILTER.tar.zst" ]] && declare -F unpack_slice >/dev/null; then
    # --num-runs>1: this run's designs are a per-uuid slice archive; expand just it.
    _SRC_SCRATCH="$(mktemp -d /tmp/refold_src.XXXXXX)"
    unpack_slice "$SOURCE_RUN/designs/$DESIGNS_FILTER.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH"
elif [[ ! -d "$DESIGNS_DIR" && -f "$SOURCE_RUN/designs.tar.zst" ]] && declare -F unpack_unit >/dev/null; then
    _SRC_SCRATCH="$(mktemp -d /tmp/refold_src.XXXXXX)"   # /tmp is bound into the container; $TMPDIR (/scratch) is not
    unpack_unit "$SOURCE_RUN/designs.tar.zst" "$_SRC_SCRATCH"
    DESIGNS_DIR="$_SRC_SCRATCH/designs"
fi
INPUTS_DIR="$REFOLD_DIR/inputs"
# Resolve REFOLD_DIR to absolute path so subsequent paths are stable across cd
REFOLD_DIR="$(cd "$REFOLD_DIR" && pwd)"
CANONICAL_DIR="$REFOLD_DIR/designs"

mkdir -p "$INPUTS_DIR" "$CANONICAL_DIR"

echo "================================================================"
echo " ESMFold2 refold worker"
echo " Source run:     $SOURCE_RUN"
echo " Refold dir:     $REFOLD_DIR"
echo " RNA sequence:   $RNA_SEQUENCE"
echo " Loops/steps:    $ESM_NUM_LOOPS / $ESM_NUM_SAMPLING_STEPS"
echo " Started:        $(date)"
echo "================================================================"
nvidia-smi -L || true

# ============================================================
# Stage 1: extract protein sequences + build per-design YAMLs
# (Boltz-style schema; identical to the Boltz-2 worker. ESMFold2's YAML parser
#  reads only id/sequence here, so the same files drive both refolders.)
# ============================================================
echo ""
echo "=== Stage 1: build per-design ESMFold2 input YAMLs ==="

apptainer exec --bind /cluster/work/wicky --bind /tmp \
    /cluster/work/wicky/containers/boltzgen.sif \
    python3 << PYEOF
import gzip
from pathlib import Path
import yaml
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

designs_dir = Path("$DESIGNS_DIR")
inputs_dir = Path("$INPUTS_DIR")
rna_seq = "$RNA_SEQUENCE"

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
        if "sequences.fa" in f.name or f.suffix in [".npz", ".json", ".csv"]:
            continue
        if filt and filt not in f.name:
            continue
        if _allowset and _stem(f.name) not in _allowset:
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
    name = path.name
    for suffix in [".cif.gz", ".cif", ".pdb"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": seq}},
            {"rna":     {"id": "B", "sequence": rna_seq}},
        ],
    }
    with open(inputs_dir / f"{name}.yaml", "w") as fh:
        yaml.dump(payload, fh, default_flow_style=False, sort_keys=False)
    count += 1

print(f"  wrote {count} ESMFold2 YAML inputs to {inputs_dir}")
PYEOF

N_INPUTS=$(ls "$INPUTS_DIR"/*.yaml 2>/dev/null | wc -l)
echo "  total YAML inputs: $N_INPUTS"
# [Fix 1] Persist expected source-design count so the aggregator/board can tell a
# leaf preempted/timed-out mid-finalize (partial refold) from a complete one.
printf '{"n_source_designs": %s}\n' "$N_INPUTS" > "$REFOLD_DIR/expected.json"

if [ "$N_INPUTS" -eq 0 ]; then
    echo "ERROR: no ESMFold2 inputs were generated. Check $DESIGNS_DIR contents." >&2
    exit 1
fi

# ============================================================
# Stage 2: ESMFold2 inference over the whole input dir (model loaded once).
# Writes canonical <design>_ESM_sample0.cif + _confidence.json into designs/.
# Fully offline: the HF cache is pre-populated, so no proxy/network needed.
# ============================================================
echo ""
echo "=== Stage 2: ESMFold2 inference ==="

WORKER_DIR="$(dirname "$0")"

# ---- Stage the ~25GB weights cache to node-local scratch. ESMFold2 mmaps the
# ESMC-6B safetensors; faulting those pages straight off Lustre is latency-bound
# (~1 MB/s on a busy GPU node -> hours of "startup"), whereas a sequential bulk
# copy runs ~400 MB/s, so staging takes ~1-2 min and the model then mmaps from
# fast local disk. Fall back to Lustre if local scratch is too small. ----
STAGE_BASE="${TMPDIR:-/tmp}"
HF_RUN_CACHE="$ESM_HF_CACHE"
EXTRA_BIND=()
NEED_KB=$(du -sk "$ESM_HF_CACHE" 2>/dev/null | awk '{print $1}')
FREE_KB=$(df -Pk "$STAGE_BASE" 2>/dev/null | awk 'NR==2{print $4}')
if [[ -n "$NEED_KB" && -n "$FREE_KB" && "$FREE_KB" -gt $((NEED_KB + 4000000)) ]]; then
    HF_RUN_CACHE="$STAGE_BASE/esm_hf_cache"
    if [[ ! -d "$HF_RUN_CACHE/models--biohub--ESMFold2" ]]; then
        echo "Staging HF cache -> $HF_RUN_CACHE ($((NEED_KB/1048576))GB; local free $((FREE_KB/1048576))GB)"
        mkdir -p "$HF_RUN_CACHE"
        time cp -a "$ESM_HF_CACHE/." "$HF_RUN_CACHE/"
    else
        echo "HF cache already staged at $HF_RUN_CACHE (resume)"
    fi
    [[ "$STAGE_BASE" != /tmp* ]] && EXTRA_BIND=(--bind "$STAGE_BASE")
else
    echo "WARN: local scratch too small (free ${FREE_KB:-?}KB, need ${NEED_KB:-?}KB) -- loading weights from Lustre"
fi

# ESMFold2's attention is a Triton kernel JIT-compiled at runtime with a C
# compiler. The container ships no gcc, so Triton falls back to the host CC the
# Euler env exports ($CC = the spack gcc under /cluster/software). Bind that tree
# read-only so the compiler resolves inside the container, and give Triton a
# writable per-leaf cache (HOME is not bound here).
apptainer exec --nv \
    --bind /cluster/work/wicky:/cluster/work/wicky \
    --bind /cluster/software:/cluster/software:ro \
    --bind /tmp \
    "${EXTRA_BIND[@]}" \
    --env HF_HOME="$HF_RUN_CACHE" \
    --env TRANSFORMERS_CACHE="$HF_RUN_CACHE" \
    --env HUGGINGFACE_HUB_CACHE="$HF_RUN_CACHE" \
    --env HF_HUB_OFFLINE=1 \
    --env TRANSFORMERS_OFFLINE=1 \
    --env TOKENIZERS_PARALLELISM=false \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env TRITON_CACHE_DIR="$REFOLD_DIR/.triton_cache" \
    --env MPLCONFIGDIR="$REFOLD_DIR/.mpl_cache" \
    --env XDG_CACHE_HOME="$REFOLD_DIR/.cache" \
    "$ESM_SIF" \
    python3 "$WORKER_DIR/esmfold2_infer.py" \
        "$INPUTS_DIR" \
        "$CANONICAL_DIR" \
        --num-loops "$ESM_NUM_LOOPS" \
        --num-sampling-steps "$ESM_NUM_SAMPLING_STEPS" \
        --num-diffusion-samples "$ESM_NUM_DIFFUSION_SAMPLES" \
        --seed "$ESM_SEED" \
        --pae-dir "$REFOLD_DIR/_ipsae_pae"

# ============================================================
# Stage 3: finalize -- pack designs/ then drop scratch.
# ============================================================
echo ""
echo "=== Stage 3: finalize ==="

N_CANON=$(ls "$CANONICAL_DIR"/*_ESM_sample0.cif 2>/dev/null | wc -l)

# --- IPSAE interface scoring (non-fatal add-on) -----------------------------
# esmfold2_infer.py dumped pae_/plddt_ npz into _ipsae_pae; score each model
# (Boltz-npz mode), writing a tiny <stem>_ipsae.json next to its .cif, then drop
# the raw PAE. Runs on the host python (numpy); never aborts the leaf.
IPSAE_DRIVER="$(dirname "$0")/../ipsae_sidecar.py"
if command -v python3 >/dev/null 2>&1 && [ -f "$IPSAE_DRIVER" ]; then
    python3 "$IPSAE_DRIVER" --mode esm --designs-dir "$CANONICAL_DIR" --raw-dir "$REFOLD_DIR/_ipsae_pae" \
        || echo "WARN: IPSAE scoring failed (non-fatal), continuing"
fi
rm -rf "$REFOLD_DIR/_ipsae_pae"

# [Fix 2] Pack designs/ FIRST so the leaf is finalized (designs.tar.zst published);
# a requeue after this hits the top early-exit.
declare -F pack_unit >/dev/null && pack_unit "$CANONICAL_DIR" || true

declare -F purge >/dev/null && purge "$INPUTS_DIR" || rm -rf "$INPUTS_DIR"
[[ -n "${_SRC_SCRATCH:-}" ]] && rm -rf "$_SRC_SCRATCH" || true

echo ""
echo "================================================================"
echo " Done at $(date)"
echo " Canonical .cif files written: $N_CANON"
echo " Refold dir: $REFOLD_DIR"
echo " Canonical: $CANONICAL_DIR"
echo "================================================================"
