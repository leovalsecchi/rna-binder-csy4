# Apo (free-binder) fold + RMSD-vs-ESM metric

**Date:** 2026-07-27
**Status:** proposed — **awaiting user approval before implementation**
**Spans:** `rna_binder` pipeline (new GPU worker + backfill, Euler `$RB/code`) and
`~/Documents/metric_explorer` (drivable metrics + data refresh).

## Goal

Measure whether each designed binder **folds to the same shape on its own as it does bound to
the RNA** — i.e. is it *pre-organized*, or only structured because the RNA scaffolds it
(induced fit)? Today's `ca_rmsd_vs_designed` is computed on the RNA-complex refold and can't
answer this. We add an **apo (free-binder) fold** and score it against the ESM complex refold.

## Why RMSD vs the ESM refold (not the designed backbone)

The apo fold and the existing `esm` refold both come from the **same ESMFold2 model** — the
only difference is the RNA is present in the refold and absent in the apo fold. So
`RMSD(apo protein, esm-refold protein)` isolates **purely the RNA's effect on the binder's
fold**, with no model or generator-intent confound. A clean induced-fit measure.

## Inputs / ground truth (the cleaned CSV)

`analyses/cross_campaign/combined_clean_metrics.csv` — **1,751,683 cleaned designs** (cys-free
+ all-three-refold + seq-deduped). Verified: **100%** have `sequence` (col 4) and
`esm_cif_path` (col 14, the RMSD reference); **every sequence is unique** (no dedup shortcut).
The **100k subset** = `app_metrics_100k.csv` (same columns) — run first.

The CSV is the work-list: `sequence` = what to fold; `esm_cif_path` = the ground-truth
structure to score against.

## What we compute (per design)

Fold the **binder protein alone** (no RNA) with the existing ESMFold2 (`esmc.sif`), then:

| column | meaning | dir |
|---|---|---|
| **`apo_rmsd_vs_esm`** | Cα RMSD of the apo fold vs the ESM refold's protein chain (induced fit) | lower = better (more pre-organized) |
| `apo_plddt` | ESMFold2 pLDDT of the lone fold (how confidently it folds alone) | higher = better |
| `apo_ptm` | ESMFold2 pTM of the lone fold | higher = better |
| `apo_cif_path` | RB-relative ref into the packed archive (`apo_refolds/chunk_NN.tar.zst!<id>_apo.cif`) | — |

Protein-only folding needs **no model change**: `esmfold2_infer.py` builds its input from a
Boltz-style YAML by iterating `sequences` entries — a YAML with only the protein entry (no
`rna`) folds the binder alone and still returns `plddt`/`ptm`.

RMSD reuses `structural.ca_rmsd` (biotite `superimpose`+`rmsd` on Cα). The ESM refold is a
complex → extract its protein chain via `structural.split_chains`. Same sequence ⇒ same residue
count ⇒ Cα superimpose is well-defined.

## Architecture

New worker **`apo_rmsd_infer.py`** (clone of `esmfold2_infer.py`): load ESMFold2 once, then per
design in the chunk →
1. build protein-only input from `sequence`, fold → apo structure + `plddt`/`ptm`;
2. resolve `esm_cif_path`, `split_chains` → ESM protein; `ca_rmsd(apo, esm_protein)`;
3. write the apo CIF to **node-local `$TMPDIR`** (not shared FS — avoids a transient inode
   spike); append `{design_id, apo_rmsd_vs_esm, apo_plddt, apo_ptm, apo_cif_path}` to the chunk CSV.

At chunk end: **`tar --zstd`** the chunk's apo CIFs from `$TMPDIR` into
`$RB/apo_refolds/chunk_NN.tar.zst` (one archive/chunk → ~1 inode/chunk), then the tempdir is
discarded. Merge: **`apo_integrate.py`** (clone of `core_integrate.py`) LEFT-joins the apo
columns onto `app_metrics_100k.csv` (idempotent, `.bak_pre_apo`); the 100k values also flow to
the explorer's `data/metrics_100k.csv`.

### Job submission (heed the GPU + memory + staging rules)
- **Container:** `esmc.sif`, `--nv`. **GPU:** `ampere_24g` (3090 eu-g4 / 4090 eu-g6, sm_80+) —
  monomer folds use *less* than the ~21 G complex, but pin ampere_24g like the ESM refold. Pool
  = `gpu.bulk` (yields to gen).
- **CRITICAL staging (from ESMFOLD2_REFOLD_NOTES):** bulk-copy the 25 GB HF cache to `$TMPDIR`
  first, then mmap from local — mmap'ing off Lustre crawls at ~1 MB/s. Bind `/cluster/software:ro`
  for the Triton JIT gcc; `TRITON_CACHE_DIR` → per-leaf writable.
- **Directives** cloned from the known-good `refold_esm` sbatch (explicit `--mem-per-cpu`, time,
  gres). Chunk size large enough (tens of thousands/worker) to amortize the ~200 s one-time
  staging.

## Scale / cost / footprint

- ~3–4 s/design steady-state (protein-only, likely ≤ the 4 s complex) + ~200 s one-time
  staging per worker.
- **100k subset ≈ 2–3 h**; **full 1.75M ≈ 1,450–1,950 GPU-h ≈ 1.5–2 days on 50 GPUs** (bulk,
  subject to availability).
- **Storage:** ~30 KB/apo CIF packed → ~3 GB (100k) / ~50 GB (1.75M). **Inodes:** ~1 archive +
  1 CSV per chunk (≈ tens total) — flat, watched throughout.

## Testing

- **Unit:** `ca_rmsd` on a structure vs a rotated/translated copy (→ ~0) and vs a perturbed one
  (→ known value); protein-only YAML builds a single-chain `StructurePredictionInput`.
- **Live smoke (few designs):** apo-fold 3–5 real binders, RMSD vs their ESM refold; assert
  sane ranges (well-behaved ~1–3 Å, induced-fit cases higher), `apo_plddt` ∈ [0,1], apo CIF
  packs + re-extracts. Compare one apo vs ESM overlay in ChimeraX by eye.
- **Then** the 100k array; verify 0-error rate + distribution + by-stack signal before the 1.75M run.

## Explorer wiring

Add `apo_rmsd_vs_esm` (lower better) + `apo_plddt` (higher better) to `METRICS_SPEC`; refresh
`data/metrics_100k.csv`. `apo_cif_path` enables opening the apo fold later (opendesign / bridge).

## Decisions (locked 2026-07-27)

1. RMSD reference = **ESM complex refold** (`esm_cif_path`), protein chain.
2. Scope = **100k subset first**, then full 1.75M.
3. Structures = **keep ALL apo CIFs, packed per-chunk `.tar.zst`** + `apo_cif_path` ref.

## Out of scope (YAGNI)

- RMSD vs the designed backbone or vs AF3 (only vs ESM, per decision 1).
- No apo fold of naturals / non-cleaned designs (cleaned CSV only).
