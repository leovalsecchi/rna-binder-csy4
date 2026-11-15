# Running ESMFold2 on Euler — the fast recipe

Written 2026-09-08 from Leo Valsecchi's (`lvalsecchi`) production setup, for a labmate
who wants to fold sequences with ESMFold2 on Euler the same way the `rna_binder`
campaigns do it.

Audience: an agent driving this on the labmate's behalf. Everything below is a real,
in-production path — nothing is aspirational.

> **Read §2 before anything else.** There are two ways to load the weights, and the
> obvious one is ~400× slower and will look like a hang, not an error.

---

## 0. What "ESMFold2" is here

`ESMFold2` = EvolutionaryScale/Biohub structure predictor: **ESMC-6B backbone + a
diffusion structure head**. It folds **protein and/or RNA/DNA chains** from sequence
alone (single-sequence — no MSA), and reports pLDDT / pTM / ipTM.

It runs **fully locally and offline** on Euler from a pre-staged HuggingFace cache —
it does *not* call the EvolutionaryScale Forge cloud API that the public Biohub docs
show, and it needs no proxy or network on the compute node.

⚠️ **Naming trap.** In the `rna_binder` project "ESM" is ambiguous. There are two
unrelated models:

| | what | where |
|---|---|---|
| **ESM2 LM** (`esm2_t36_3B_UR50D`) | protein *language* model, embeddings/scoring | inside HalluDesign **generation** |
| **ESMFold2** | *structure* predictor, the thing in this doc | the **refold**/fold backend |

Always say "ESMFold2". This doc is only about ESMFold2.

---

## 1. Fixed assets — paths, owners, permissions

All of these are readable by anyone in `bsse-wicky-euler-storage`. Nothing needs copying.

| Asset | Path | Notes |
|---|---|---|
| **Container** | `/cluster/work/wicky/containers/esmc.sif` | 3.5 GB, built by `mescartin` 2026-06-09. **This is the one production uses.** |
| **Weights (HF cache)** | `/cluster/work/wicky/lvalsecchi/rna_binder/params/esmfold2_hf_cache/huggingface` | **25 GB**. Holds `models--biohub--ESMFold2` + `models--biohub--ESMC-6B` (+ `xet`). Owned by `lvalsecchi`, group-readable. |
| **Inference driver** | `/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/esmfold2_infer.py` | Loads the model **once**, folds a whole directory of YAML inputs. Resume-safe. |
| **Protein-only YAML builder** | `/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/apo_make_yaml.py` | TSV (`id`, `sequence`) → one YAML per entry. Trivial, copy it if you prefer. |
| **Reference job script** | `/cluster/work/wicky/lvalsecchi/rna_binder/apo_fold_array_v9.sh` | The cleanest end-to-end example: array job, node-local staging, packs results. **Best template to copy.** |

Model identifier used inside the container: `biohub/ESMFold2`, loaded as
`transformers.models.esmfold2.modeling_esmfold2.ESMFold2Model.from_pretrained(...)`,
with inputs built via `esm.models.esmfold2` (`ProteinInput` / `RNAInput` / `DNAInput`
→ `StructurePredictionInput` → `ESMFold2InputBuilder().fold(...)`).

**Don't copy the 25 GB cache.** Point `HF_HOME` at Leo's copy and stage it per-job (§2).

---

## 2. THE TWO VERSIONS — which one to use

There is more than one "version" floating around. Here is every axis, and the answer.

### 2a. Weight loading — **this is the optimization that matters**

| | how | startup |
|---|---|---|
| ❌ **Slow (naive)** | `HF_HOME` → the Lustre cache directly | **hours**, often times out |
| ✅ **Fast (use this)** | bulk-copy the cache to node-local `$TMPDIR` first, then `HF_HOME` → the local copy | **~200 s**, then ~4 s/design |

**Why.** ESMFold2 `mmap`s the ESMC-6B safetensors (~25 GB). Faulting those pages
straight off Lustre is a random-read pattern that runs at **~1 MB/s** on a busy GPU
node — measured: `read_bytes` grew ~8 MB in 8 s while the GPU sat at 0%. A
**sequential** bulk copy of the exact same files runs at **~400 MB/s**. So it was never
the disk, only the access pattern. Copy first (~1–2 min for 24 GB), then mmap from
fast local NVMe.

**Failure mode if you get this wrong:** the job looks alive, GPU utilisation is 0%,
nothing is written, and it eventually hits the walltime. It does *not* error out. If
you see that, this is why.

The staging block is ~15 lines and is included verbatim in the template in §4. It
falls back to Lustre if local scratch is too small, and skips the copy if the cache is
already staged (so a requeue on the same node is cheap).

### 2b. Container — two exist, use `esmc.sif`

| Image | What | Use it? |
|---|---|---|
| `/cluster/work/wicky/containers/esmc.sif` | The lab's original ad-hoc image (`mescartin`, 2026-06-09) | ✅ **Yes** — every production refold ran on this |
| `/cluster/work/wicky/lvalsecchi/esmfold2_build/esmfold2.sif`<br>(also symlinked `containers/esmfold2_mlm.sif`) | Leo's *pinned, rebuildable* Apptainer recipe (torch 2.12 / cu130, Biohub `transformers` + `esm` at fixed SHAs). Build recipe: `~/Documents/esmfold2/esmfold2.def` | ⚠️ Only if you specifically need a reproducible/rebuildable image. **Never validated for the refold pipeline.** |

### 2c. Worker script — two variants in the repo

| Script | I/O layout | Status |
|---|---|---|
| `code/scripts/stages/refold_esm.sh` | intermediates on Lustre (`$REFOLD_DIR/inputs`, `designs/`, caches) | ✅ **Production.** This is what `launch.py` calls. |
| `code/scripts/stages/refold_esm_scratch.sh` | *all* heavy I/O on `$TMPDIR`; only the single packed `designs.tar.zst` is copied back to Lustre | ⚠️ Later (2026-07-04) inode-friendly refactor. **Never wired into `launch.py`** — its only caller anywhere is a one-off single-uuid probe, `ops/refold/_esm_scratch_test.sbatch`. |

Both contain the §2a staging fix — they differ only in where the *intermediate* files
live, not in speed of the model load.

**For a fresh standalone run, don't use either.** Use the `apo_fold_array_v9.sh`
pattern (§4): it is the newest, is already all-node-local, and has no `rna_binder`
campaign plumbing to satisfy.

### 2d. Also two, but not versions of each other

`esmfold2_infer.py` has two special-purpose forks that you should **ignore** unless you
know you want them:
- `.../scratch_folds/csy4_dualrepeat_2026-07-06/esmfold2_infer_patched.py` — standalone
  lone-entity folding
- `.../seedtest_2026-07-30/drivers/esmfold2_infer_seedtest.py` — multi-seed
  self-consistency sweeps

The canonical driver is `code/scripts/stages/esmfold2_infer.py`.

---

## 3. Input and output formats

### Input: one Boltz-style YAML per structure, all in one directory

`esmfold2_infer.py` reads only `id` and `sequence` per entity. Supported entity keys:
`protein`, `rna`, `dna`.

Protein alone (apo fold):
```yaml
version: 1
sequences:
- protein:
    id: A
    sequence: ATYKVTYNISGLTTEELEEILKEIKAIAEKHGATVRSRRVSK...
```

Protein + RNA complex (what the campaigns refold):
```yaml
version: 1
sequences:
- protein:
    id: A
    sequence: ATYKVTYNISGLTT...
- rna:
    id: B
    sequence: ACUGCCGUAUAGGCAG
```

The **file stem is the output name** — `foo.yaml` → `foo_ESM_sample0.cif`. Make stems
unique and filesystem-safe.

### Output: per input, into the canonical dir

| File | Content |
|---|---|
| `<stem>_ESM_sample0.cif` | the predicted structure (mmCIF) |
| `<stem>_ESM_sample0_confidence.json` | `{"iptm", "ptm", "plddt", "complex_plddt", "pair_chains_iptm"}` — `plddt` is the **mean** over residues; `pair_chains_iptm` is the chain-pair ipTM matrix (chain A=0, B=1) |

Optional, only if you pass `--pae-dir DIR`:
`pae_<stem>.npz` + `plddt_<stem>.npz` (Boltz npz layout, for downstream ipSAE/pDockQ2
scoring). Costs disk; skip unless you need it.

### Driver CLI

```
python3 esmfold2_infer.py INPUTS_DIR OUTPUT_DIR \
    [--num-loops 20] [--num-sampling-steps 100] \
    [--num-diffusion-samples 1] [--seed 0] \
    [--model biohub/ESMFold2] [--pae-dir DIR]
```

Production values are exactly the defaults above. Behaviour worth knowing:
- **Model loaded once** for the whole directory — so give each job a *big* batch (§6).
- **Resume-safe**: any input whose `_ESM_sample0.cif` already exists is skipped. Writes
  to a temp name then renames, so a preemption can never leave a half-written CIF that
  the skip-guard would mistake for done.
- **One bad input does not kill the batch** — it prints `[FAIL] <name>: <exc>` and
  continues. Grep the log for `[FAIL]` afterwards.

---

## 4. Copy-paste job template (the fast recipe)

Generalised from `apo_fold_array_v9.sh`. Adjust `ACCOUNT`, the paths, and the manifest.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=esmfold2
#SBATCH --account=es_wicky
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --gpus=1
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --tmp=60000                 # ask SLURM for >=60 GB node-local scratch (see §5)
#SBATCH --output=%x_%j.log
# partition / time / nodelist: pass on the sbatch CLI, see §5
set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true

ESM_SIF=/cluster/work/wicky/containers/esmc.sif
ESM_HF_CACHE=/cluster/work/wicky/lvalsecchi/rna_binder/params/esmfold2_hf_cache/huggingface
INFER=/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/esmfold2_infer.py

MANIFEST=$1            # TSV with header: id <TAB> sequence
OUTDIR=$2              # where the final results should land (Lustre)

export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID:-$$}}"; mkdir -p "$TMPDIR"
WORK="$TMPDIR/esm_${SLURM_JOB_ID:-$$}"
INPUTS="$WORK/inputs"; DESIGNS="$WORK/designs"
mkdir -p "$INPUTS" "$DESIGNS"
trap 'rm -rf "$WORK" 2>/dev/null || true' EXIT
nvidia-smi -L || true

# --- 1. build one YAML per sequence (host python) ---
python3 /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/apo_make_yaml.py \
        "$MANIFEST" "$INPUTS"

# --- 2. THE SPEED FIX: stage the 25 GB HF cache to node-local scratch ---
#     mmap off Lustre = ~1 MB/s (hours). Sequential bulk copy = ~400 MB/s (~1-2 min).
STAGE_BASE="$TMPDIR"; HF_RUN_CACHE="$ESM_HF_CACHE"; EXTRA_BIND=()
NEED_KB=$(du -sk "$ESM_HF_CACHE" 2>/dev/null | awk '{print $1}')
FREE_KB=$(df -Pk "$STAGE_BASE" 2>/dev/null | awk 'NR==2{print $4}')
if [[ -n "$NEED_KB" && -n "$FREE_KB" && "$FREE_KB" -gt $((NEED_KB + 4000000)) ]]; then
    HF_RUN_CACHE="$STAGE_BASE/esm_hf_cache"
    if [[ ! -d "$HF_RUN_CACHE/models--biohub--ESMFold2" ]]; then
        echo "staging HF cache -> $HF_RUN_CACHE"
        mkdir -p "$HF_RUN_CACHE"; time cp -a "$ESM_HF_CACHE/." "$HF_RUN_CACHE/"
    else
        echo "HF cache already staged (resume)"
    fi
    [[ "$STAGE_BASE" != /tmp* ]] && EXTRA_BIND=(--bind "$STAGE_BASE")
else
    echo "WARN: local scratch too small -- falling back to Lustre (will be SLOW)"
fi

# --- 3. fold (model loaded once for the whole batch) ---
apptainer exec --nv \
    --bind /cluster/work/wicky:/cluster/work/wicky \
    --bind /cluster/software:/cluster/software:ro \
    --bind /tmp "${EXTRA_BIND[@]}" \
    --env HF_HOME="$HF_RUN_CACHE" \
    --env TRANSFORMERS_CACHE="$HF_RUN_CACHE" \
    --env HUGGINGFACE_HUB_CACHE="$HF_RUN_CACHE" \
    --env HF_HUB_OFFLINE=1 \
    --env TRANSFORMERS_OFFLINE=1 \
    --env TOKENIZERS_PARALLELISM=false \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env TRITON_CACHE_DIR="$WORK/.triton" \
    --env MPLCONFIGDIR="$WORK/.mpl" \
    --env XDG_CACHE_HOME="$WORK/.cache" \
    "$ESM_SIF" \
    python3 "$INFER" "$INPUTS" "$DESIGNS" \
        --num-loops 20 --num-sampling-steps 100 \
        --num-diffusion-samples 1 --seed 0

N=$(ls "$DESIGNS"/*_ESM_sample0.cif 2>/dev/null | wc -l)
echo "folds produced: $N"
[ "$N" -eq 0 ] && { echo "ERROR: no folds"; exit 1; }

# --- 4. publish results: ONE archive on Lustre, not thousands of small files ---
mkdir -p "$OUTDIR"
( cd "$DESIGNS" && tar -cf - . | zstd -q -19 -T0 -o "$TMPDIR/out.tar.zst" )
zstd -tq "$TMPDIR/out.tar.zst"                 # verify before publishing
mv -f "$TMPDIR/out.tar.zst" "$OUTDIR/folds_${SLURM_JOB_ID:-$$}.tar.zst"
```

### Why each non-obvious line is there

- **`--bind /cluster/software:ro`** — ESMFold2's attention includes a **Triton kernel
  JIT-compiled at runtime**, which needs a C compiler. `esmc.sif` ships **no gcc**, so
  Triton falls back to the host `$CC` (the spack gcc under `/cluster/software`). Without
  this bind the compiler doesn't resolve inside the container and the job dies at
  startup. `module load stack/...` is what puts `$CC` in the environment.
- **`TRITON_CACHE_DIR` / `MPLCONFIGDIR` / `XDG_CACHE_HOME`** — `$HOME` is not bound
  inside the container, so these must point at a writable dir or Triton/matplotlib
  crash trying to write their caches.
- **`--bind "$TMPDIR"`** (via `EXTRA_BIND`) — `/tmp` is bound by default but `$TMPDIR`
  (`/scratch/...`) is **not**. If you stage the cache there you must bind it explicitly.
- **`HF_HUB_OFFLINE=1` / `TRANSFORMERS_OFFLINE=1`** — compute nodes have no internet.
  Without these, `from_pretrained` tries to reach huggingface.co and hangs/fails.
- **`--requeue` + `--open-mode=append`** — required on preemptible partitions, so a
  requeued job resumes (the skip-guard picks up where it stopped) instead of truncating
  its own log.
- **verify the archive before `mv`** — a blank-but-correctly-named output is the classic
  silent failure here; `zstd -tq` catches it.

---

## 5. SLURM sizing and routing

### GPU requirements

- **VRAM: ~21 GB** for a small complex (~80 aa protein + 16 nt RNA) → needs a **24 GB
  card**. Bigger complexes **OOM** (validated: antibody–antigen pairs OOM'd 3 of 46).
  Budget VRAM by total token count, and test one large case before a big fan-out.
- **Architecture: not restricted to Ampere.** Older config comments claim ESMFold2 needs
  `sm_80+`; that was **tested and corrected on 2026-06-22** — `esmc.sif` has no mandatory
  fused kernel and logs *"neither xformers nor flash-attn is installed; falling back to
  PyTorch F.scaled_dot_product_attention"*. It runs fine on **Turing (sm_75)** 24 GB
  cards, just slower (no fused attention). Blackwell `rtx_pro_6000` is **untested/known
  bad** for these containers (PyTorch has no `sm_120`).

### Euler node → GPU map

| Nodes | GPU | VRAM | Arch |
|---|---|---|---|
| `eu-lo-g2-*` | RTX 2080 Ti | 11 G | Turing — **too small** |
| `eu-g3-[062-066,083-084]`, `eu-lo-g3-[049-061]` | Quadro RTX 6000 | 24 G | Turing |
| `eu-g3-[067-079]` | Titan RTX | 24 G | Turing |
| `eu-g4-[001-032]` | RTX 3090 | 24 G | Ampere |
| `eu-g6-[001-080]` | RTX 4090 | 24 G | Ada |
| `eu-a65-*` | A100 / A100-80 | 40/80 G | Ampere (`gpupr.*`, reserved — avoid) |
| `eu-g7-*` | RTX PRO 6000 | 96 G | Blackwell — **avoid** |

### The routing rule you must know

**`--gres=gpumem:Ng` is a discrete pool selector, not a "≥ N GB" filter.** It decides
which partition family your job lands in, and the scheduler will **override your
`--partition`**:

| request | lands in |
|---|---|
| `gpumem:8g` | public `gpu.*` (11 GB 2080 Ti — **too small for ESMFold2**) |
| `gpumem:24g` | `gpuhe.*` |
| `gpumem:40g` | `gpupr.*` (reserved — jobs sit `ReqNodeNotAvail` forever) |

So there are exactly **two sane ways** to get a 24 GB card:

1. **`--gres=gpumem:24g`** → routed to `gpuhe`. Simple. This is what the production
   config does.
2. **No `gpumem` at all**, `--gpus=1` plus an explicit
   `--nodelist=eu-g3-[062-077,083-084]` → stays in the **public `gpu.*`** pool on the
   24 GB Turing cards. This is what `apo_fold_array_v9.sh` does, to tap capacity that
   `gpuhe` doesn't have. Slower per design, but far more nodes.

### Which partition

Use the **`.bulk`** tier — `gpu.bulk` or `gpuhe.bulk`:
- `.bulk` is **PriorityTier 250**; the timed tiers (`.4h/.24h/.120h`) are **Tier 1000**,
  the same tier as design generation, so heavy fold jobs there starve everything else.
- `.bulk` is `MaxTime=15-00:00:00`, `AllowQos=ALL`, and the `es_wicky/{gpu,gpuhe}/bulk`
  QOS has **no MaxWall** — so a direct `sbatch --time=06:00:00 -p gpu.bulk` is accepted.
- **Set the walltime at submit.** `scontrol update` can set partition/QOS but **not**
  `TimeLimit`, and a user may only *lower* a running job's limit, never raise it. The
  timed partitions reject `>4h`, so submit *directly* to `.bulk` if you want 6 h.
- `.bulk` is **preemptible** — hence `--requeue`, `--open-mode=append`, and the
  resume-safe skip-guard. This combination is why preemption is cheap here.

A bare `-p bulk` also works at any walltime; `gpumem` then picks the pool and the QOS
is auto-matched.

### Node-local scratch

Request it explicitly with **`--tmp=60000`** (MB) or more: you need ~25 GB for the
staged cache plus room for inputs/outputs. Without a `--tmp` request `$TMPDIR` may be
tiny or unset, the staging block falls back to Lustre, and you're back to the slow path.
A node whose `/scratch` fills up produces **blank-but-correctly-named output files** —
so always verify content, not just file counts.

---

## 6. Time model and batching

Measured on 1× RTX 3090, ~80 aa protein + 16 nt RNA, defaults (20 loops / 100 steps):

```
total ≈ 200 s  +  N × 4 s
```

- stage 24 GB + model load + Triton compile (one-time): **~200 s**
- per design, steady state: **~4 s**

| N designs | wallclock |
|---|---|
| 50 | ~7 min |
| 600 | ~45 min |
| 5000 | ~6 h |

**The startup is a fixed ~200 s tax per job, so batch aggressively.** 10 jobs of 60
designs costs ~33 min of pure startup; 1 job of 600 costs ~3 min. Aim for chunks that
run 2–6 h. Turing cards are slower per design than the numbers above — scale up your
walltime estimate if you pin `eu-g3`.

Sanity values from the validation run (job 3743276, 3 real designs):
`pLDDT 0.882 / pTM 0.876 / ipTM 0.773`, `0.631 / 0.595 / 0.136`, `0.761 / 0.740 / 0.456`.

---

## 7. Failure signatures — diagnose fast

| Symptom | Cause | Fix |
|---|---|---|
| Job alive, GPU 0%, nothing written, hits walltime | Weights being mmap'd off **Lustre** (§2a) | Add the staging block + `--tmp=60000` |
| `RuntimeError: 0 active drivers ([]). There should only be one.` (Triton) | GPU/CUDA driver not visible — the node's GPU is wedged. **Not** a data or code bug | Requeue; it usually lands elsewhere and succeeds. If it repeats on one node, exclude it: `--exclude=<node>` |
| Triton fails to compile / no C compiler | `/cluster/software` not bound, or `$CC` unset | `module load stack/2024-06 ...` on the host **and** `--bind /cluster/software:ro` |
| `from_pretrained` hangs or DNS-errors | offline flags missing | `HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1` |
| CUDA OOM | complex too large for 24 GB | fewer/shorter chains, or a bigger card (test one case first) |
| Job pending forever, `ReqNodeNotAvail` | `gpumem:40g` routed you into reserved `gpupr.*` | drop to `gpumem:24g`, or use the nodelist route |
| Permission-denied writing a cache | `$HOME` isn't bound in the container | set `TRITON_CACHE_DIR` / `MPLCONFIGDIR` / `XDG_CACHE_HOME` to node-local paths |
| Output files exist with the right names but are empty/blank | node `/scratch` filled up mid-run | raise `--tmp`, re-run the chunk, and verify content not counts |

Two logging notes:
- Slurm logs are opened `--open-mode=append`, so **a touched `.err` does not mean the
  job failed** — most jobs hit a transient error, requeue, and finish. Count real
  failures with `sacct --state=FAILED`, never by `.err` mtime.
- Detect a hung job by **StdErr mtime + `nvidia-smi` at 0%**, not by the slurm log.

---

## 8. If you want the full campaign refold pipeline (not just "fold these sequences")

Only relevant if the labmate wants to refold *designs* (extract each design's protein
sequence from a CIF/PDB, pair it with a target RNA, fold, and score) rather than fold a
list of sequences she already has.

Base: `/cluster/work/wicky/lvalsecchi/rna_binder/code/`

| Piece | Path |
|---|---|
| Config | `configs/refold__esm/default.yaml` (`config_kind: refold_esm`) |
| Worker | `scripts/stages/refold_esm.sh` |
| Inference | `scripts/stages/esmfold2_infer.py` |
| Launcher | `scripts/launch.py` → `build_refold_esm_sbatch`, registered as `esm` |
| Write-ups | `ESMFOLD2_REFOLD_NOTES.md`, `ESMFOLD2_REFOLD_CHANGES.md`, `ESMFOLD2_SPEEDFIX.md`, `ESMFOLD2_DOWNSTREAM_PLAN.md`, `GPU_REFOLD_COMPAT.md` |

Three stages: (1) biopython inside `boltzgen.sif` extracts the longest protein chain
from each design and writes a Boltz-style YAML (`protein A` + `rna B`); (2)
`esmfold2_infer.py` runs once over the whole input dir inside `esmc.sif`; (3) results
are packed to `designs.tar.zst` and scratch is purged.

How it's driven:
```bash
# automatic, fans out an ESM leaf per source run
python launch.py configs/rfd3__multif/default.yaml --auto-refold

# tune at submit
--refold-esm-override refolder_params.num_sampling_steps=50

# standalone, one existing run
ESM_SIF=... ESM_HF_CACHE=... \
  bash scripts/stages/refold_esm.sh <refold_dir> <source_run> <rna_seq> [uuid_filter]
```

Config defaults: `num_loops: 20`, `num_sampling_steps: 100`,
`num_diffusion_samples: 1`, `seed: 0`, `use_rna_msa: false` (single-sequence; the
ESMFold2 MSA path is not wired), RNA `ACUGCCGUAUAGGCAG` (16 nt csy4).

**Caveat for a labmate:** this pipeline is wired to `lvalsecchi`'s campaign directory
layout, archive helpers, and aggregator. Reusing it wholesale means inheriting all of
that. For "I have sequences, fold them", **§4 is the right entry point** — copy it and
own it.

---

## 9. Minimal checklist

1. Confirm group access: `ls -l /cluster/work/wicky/containers/esmc.sif` and
   `ls /cluster/work/wicky/lvalsecchi/rna_binder/params/esmfold2_hf_cache/huggingface`
   (both should list — you need `bsse-wicky-euler-storage`).
2. Write a TSV: header `id<TAB>sequence`, unique filesystem-safe ids.
3. Copy the §4 template. Set `--account`, `--tmp=60000`, `-p gpu.bulk`,
   `--time=06:00:00`, and **either** `--gres=gpumem:24g` **or**
   `--nodelist=eu-g3-[062-077,083-084]`.
4. **Smoke test on 3 sequences first.** Expect ~200 s startup then ~4 s each. Confirm
   the log shows `staging HF cache -> ...` and `[OK] <name>: pLDDT ... pTM ... ipTM ...`.
5. Open one output CIF and one `_confidence.json` and check they're real. Do not trust
   file counts.
6. Scale up: chunk into 2–6 h batches, submit as an array, `--requeue`
   `--open-mode=append`.
7. Grep the logs for `[FAIL]` and count terminal failures with
   `sacct -j <jobid> --state=FAILED`.
