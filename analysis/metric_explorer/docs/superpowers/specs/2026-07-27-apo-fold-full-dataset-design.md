# Apo-fold RMSD — full-dataset scale-up (v2–v7)

**Date:** 2026-07-27
**Status:** proposed — **awaiting user approval before implementation**
**Builds on:** `2026-07-27-apo-fold-rmsd-metric-design.md` (the metric + 100k backfill,
now 98/100 chunks done). This spec scales the **same pipeline** from the 100k subset to the
full cleaned dataset.
**Spans:** `rna_binder` pipeline on Euler (`$RB` = `/cluster/work/wicky/lvalsecchi/rna_binder`).
No metric-explorer changes (the 3.8 GB full table is not loaded into the explorer; the 100k
table keeps its existing apo columns).

## Goal

Run the already-validated apo-fold + `apo_rmsd_vs_esm` metric across **all of
`analyses/cross_campaign/combined_clean_v2v7.csv`** (v2 through v7), not just the 100k subset.
Reuse the exact worker/scoring/packing that produced 98,000 folds with 0 errors.

## Input & scope

- **Source:** `analyses/cross_campaign/combined_clean_v2v7.csv` — **3,913,151 rows**, 3.8 GB.
  Columns present and sufficient: `design_id, stack, uuid8, sequence, esm_cif_path,
  designed_cif_path, …`. Verified **`has_esm` = 100%**, so every row yields a valid
  `apo_rmsd_vs_esm`.
- **Exclusion:** skip the ~100k `design_id`s already folded (union of `design_id` across the
  existing `$RB/apo_chunks/chunk_*_apo.csv`). → **~3.81M new designs to fold.**
- **Namespace-overlap guard (must run first):** compute the set-intersection between the 100k
  apo `design_id`s and `combined_clean_v2v7.design_id`. If the overlap is ~100k → exclusion is
  correct. If the overlap is ~0 → the two files use different ID conventions; **STOP and report**
  rather than silently re-folding 100k or excluding nothing. This is a cheap one-time check.

## Chunking

- **1,000 designs / chunk** (uniform across all four partitions). This is the proven size:
  measured mean **115.6 min** on the 24 G Turing pin → ~124 min of margin under the tightest wall
  (**gpu.4h**, survives a 2× slow node), and comfortably inside bulk's 6 h. 2,000/chunk was
  rejected — ~231 min leaves only ~8 min under the 4 h wall.
- **~3,810 chunks.** Higher chunk count is a non-issue: per-design resume makes any
  TIMEOUT/preempt cheap, and self-consolidating arrays keep the 250-job scheduler window from
  flooding.
- **Manifests:** `apo_full_chunks/chunk_NNNN.tsv` (`design_id, sequence, esm_cif_path, stack`),
  **stratified round-robin by stack** so every chunk is stack-balanced and chunks completing in
  submission order form a representative cross-stack sample.
- Built by new **`apo_build_manifests_full.py`** (streams the 3.8 GB CSV, applies the exclusion
  set, writes 4-digit chunk TSVs). Exact chunk count is fixed at manifest-build time.

## Compute — all four public partitions, bulk-majority

Run concurrently across **`gpu.bulk` + `gpu.4h` + `gpu.24h` + `gpu.120h`** to pull the most GPUs
at once. Each timed partition has its own QOS GPU cap (4h=96, 24h=50, 120h=29) so running all of
them unlocks far more concurrent capacity than bulk alone; **bulk carries the ~80% majority** so
most of the load stays Tier-250 and yields to gen, while the timed partitions add parallel
throughput on top. All four share the same worker, chunk size, and node pin — only the partition/
QOS/walltime differ per family.

| Partition | Share | ~Chunks | Throttle | QOS | Walltime |
|---|---|---|---|---|---|
| **gpu.bulk** | **80%** | ~3,050 | `%40` | `es_wicky/gpu/bulk` | 6 h, `--requeue --open-mode=append` |
| gpu.4h | 11% | ~420 | `%25` | `es_wicky/gpu/4` (cap 96) | ≤4 h |
| gpu.24h | 6% | ~230 | `%20` | `es_wicky/gpu/24` (cap 50) | ≤24 h |
| gpu.120h | 3% | ~115 | `%12` | `es_wicky/gpu/120` (cap 29) | ≤120 h |

- All pinned to **`eu-g3-[062-077,083-084]`** (24 G titan_rtx / quadro_rtx_6000): the OOM-safe
  set — ESM needs ~13.5 G, so 2080Ti (11 G) is excluded **by node pin**, not a fragile `gpumem`
  flag. `--gpus=1`, no `gpumem` gres.
- **1,000/chunk works on all four** — measured ~1.93 h ≪ the 4 h wall, so the same worker is safe
  on the tightest (gpu.4h) partition.
- **Routing:** timed partitions are submitted directly with the matching `--time`; **`gpu.bulk` is
  not directly submittable** (walltime routing) — reach it the same way the refolds do, via
  `scontrol update partition=gpu.bulk` right after submit (keeps the array compact). Chunk-ID
  ranges: bulk 0762–3809, gpu.4h 0000–0419, gpu.24h 0420–0649, gpu.120h 0650–0761.
- Public split weighted by QOS GPU cap (96:50:29). Throttles are initial upper bounds; effective
  concurrency is additionally mem-gated by QOS (`--mem-per-cpu=6G × 8 cpus = 48 G/task`).
- Self-consolidating arrays keep the array-element count out of the 250-window; per-design resume
  (below) makes bulk preemption cheap.

## Per-design resume (new; makes bulk preemption cheap)

`apo_fold_array.sh` gains a resume guard so a requeued/preempted bulk chunk continues instead of
redoing work:

1. On start, if `apo_full_chunks/chunk_NNNN_apo.csv` exists, read its completed `design_id`s
   into a skip-set (append mode).
2. Build protein-only YAMLs (node-local `$TMPDIR`) **only for not-yet-done designs**.
3. Fold → score → **append** rows to the chunk CSV (`--open-mode=append`, dedup-safe by
   design_id).
4. Repack the chunk tar from the full set of apo CIFs at chunk end (idempotent).

With per-design resume, a preempted chunk loses at most the single in-flight fold, so 1,000-design
chunks are safe on bulk despite preemption.

## Score & pack (reused unchanged)

- **`apo_rmsd_score.py`** → `apo_full_chunks/chunk_NNNN_apo.csv`
  (`design_id, stack, apo_rmsd_vs_esm, apo_plddt, apo_ptm, apo_cif_path, apo_error`).
- Pack apo CIFs (node-local) → **`apo_refolds_full/chunk_NNNN.tar.zst`**.
- **Distinct dirs** from the 100k run (`apo_full_chunks/`, `apo_refolds_full/`) so nothing
  collides with the existing `apo_chunks/` / `apo_refolds/`.

## Integration

- Accumulate all results (100k + 3.81M) into one **`apo_master.csv`** keyed on `design_id`
  (`design_id, apo_rmsd_vs_esm, apo_plddt, apo_ptm, apo_cif_path, apo_error`), built by
  concatenating `apo_chunks/chunk_*_apo.csv` + `apo_full_chunks/chunk_*_apo.csv` (dedup on
  design_id, last-write-wins).
- **`apo_integrate_full.py`** LEFT-merges `apo_master` → `combined_clean_v2v7.csv` on
  `design_id` (idempotent; `.bak_pre_apo_<stamp>` backup; adds columns, never drops rows).
- The explorer's `app_metrics_100k.csv` already carries apo columns from the 100k integrate —
  unchanged here.

## Budgets

- **Storage:** 12.3 MB/1,000 designs (measured) → **~48 GB** of apo CIFs total. Trivial.
- **Inodes:** per-chunk pack-and-delete → ~3,810 tars + ~3,810 logs + manifests ≈ **~8k inodes**,
  flat over time (no per-design leak). Same pattern that kept the 100k run flat at ~5.99M.
- **GPU:** ~1.93 GPU-h/1,000 → **~7,350 GPU-h** for 3.81M. Yielding to the campaign on bulk,
  realistically **~1–2 weeks** wall time (faster if the cluster is quiet).

## Monitoring

On-demand (no long-lived background monitor — they get reaped):
`done = ls apo_full_chunks/chunk_*_apo.csv | wc -l`; state via `squeue -r`; aggregate
errors/percentiles via the same python one-liner used for the 100k run. Self-consolidating
arrays keep the array-element count from clogging the 250-window.

## Files

**New (Euler `$RB/code/scripts/`):**
- `apo_build_manifests_full.py` — stream CSV, apply exclusion set, write 4-digit chunk TSVs.
- `apo_integrate_full.py` — build `apo_master.csv`, merge into `combined_clean_v2v7.csv`.

**Modified:**
- `apo_fold_array.sh` — per-design resume guard; 4-digit chunk ids; `apo_full_chunks/` +
  `apo_refolds_full/` paths; tier params on CLI.

**Reused unchanged:** `apo_make_yaml.py`, `apo_rmsd_score.py`, `esmfold2_infer.py`, `esmc.sif`.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| 100k/full `design_id` namespaces differ → wrong exclusion | Overlap guard runs first; STOP-and-report if ~0 overlap. |
| Bulk preemption wastes work | Per-design resume + `--requeue --open-mode=append`. |
| 2080Ti (11 G) OOM | Node pin to `eu-g3-[062-077,083-084]` 24 G only; no `gpumem` flag. |
| Array flood of 250-window | Self-consolidating arrays; bulk `%40`, public throttled. |
| Starving the live campaign | Bulk carries 80% (Tier 250, yields). Timed partitions (20%) are Tier 1000 and DO compete with gen in the top-250 window — accepted tradeoff for throughput; throttled (`%25/%20/%12`) and capped at 20% to bound the pressure. Dial the timed shares down if gen visibly starves. |
| 3.8 GB CSV merge is slow/RAM-heavy | Stream-merge on `design_id` (chunked read), not a full pandas outer-join in RAM. |

## Defaults chosen (tunable)

1,000 designs/chunk; **80/11/6/3 across gpu.bulk / gpu.4h / gpu.24h / gpu.120h** at
`%40/%25/%20/%12`; integrate into `combined_clean_v2v7.csv` only (not the explorer). Change any
before launch.
