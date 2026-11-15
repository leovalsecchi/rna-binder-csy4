# Seed self-consistency test — 2,000 designs × 3 backends × 5 seeds

**Date:** 2026-07-30
**Status:** Approved (design) — ready for implementation plan
**Location when deployed:** Euler `$RB/seedtest_2026-07-30/` (self-contained)

## Question

How reproducible is a refold? Each backend is stochastic — the same sequence folded twice
gives slightly different structures. We want to quantify that spread: **if we fold the same
design 5 times with different seeds, how much do the pose and the confidence scores move?**

This tells us how much of the variation in our existing metrics is real signal versus
sampling noise, and whether "good" designs are also more seed-stable.

## Requirements (locked with the user)

1. **Subsample: random 400 designs per stack = 2,000 designs**, drawn from the 100k set
   (`app_metrics_100k.csv` / local `data/metrics_100k.csv`; 20,000 per stack).
   Uniform random within each stack, fixed RNG seed for reproducibility.
2. **All 3 refold backends**: AF3, Boltz-2, ESMFold2.
3. **5 fresh draws per design per backend.** The existing refold is NOT reused as a draw
   (and could not be for B2 — see below, B2 currently runs unseeded).
4. **Metrics: pairwise DockQ (pose agreement) + confidence spread (iPTM/pTM/pLDDT).**
5. **NO modifications to production drivers.** `refold_{af3,b2,esm}.sh` and
   `esmfold2_infer.py` in `code/scripts/stages/` MUST be left untouched. Every change lives
   in a LOCAL COPY inside `seedtest_2026-07-30/drivers/`.

## Ground truth on seeding (verified 2026-07-30, do not re-litigate)

| backend | mechanism | verified how |
|---|---|---|
| AF3 | `MODEL_SEEDS_JSON` is a JSON **array** → N independent seeds in ONE pass. Canonicaliser already emits `_AF3_seed<S>_sample<N>.cif` for every seed. | read `refold_af3.sh` (`"modelSeeds": model_seeds`, lines ~125/212, canon loop ~316-327) |
| Boltz-2 | `--seed INTEGER` EXISTS ("Seed to use for random number generator. Default is None (no seeding)") + `--diffusion_samples N`. The **production driver never passes `--seed`**, so all existing B2 refolds are unseeded/unreproducible. Canonicaliser already handles `_B2_sample<N>`. | `apptainer exec /cluster/work/wicky/containers/boltz2.sif boltz predict --help` |
| ESMFold2 | `builder.fold(..., seed=int\|None)` = "seed for the **diffusion sampler**" + `num_diffusion_samples`. Driver passes `--seed`, BUT the canonicaliser keeps only `*_ESM_sample0.cif` — samples 1..4 are discarded today. | `esmfold2_infer.py` (~line 118-125); biohub/ESMFold2 API docs |

**Why 5 diffusion samples ≈ 5 seeded runs (and why the 3 backends are comparable):** all three
are a *deterministic trunk* (AF3 Pairformer / Boltz trunk / ESMC-6B embeddings) followed by a
*stochastic diffusion sampler*. Inference dropout is off, so the trunk is fixed given the
sequence; ALL the randomness is in the diffusion sampler. Five draws from one seeded run and
five separately-seeded runs are therefore both i.i.d. samples of the same diffusion
distribution — statistically equivalent, but the multi-sample form is far cheaper. This is why
we use each backend's native multi-sample path rather than 5× separate jobs.

## Design

### 1. Subsample → manifest
`scripts/build_seedtest_manifest.py`
- Read the 100k table; group by `stack`; `sample(n=400, random_state=20260730)` per stack.
- Emit `seedtest_manifest.tsv`: `design_id, stack, sequence, rna_target, uuid8`.
- Sanity: exactly 2,000 rows, 400 per stack, all `sequence` non-null and length 80,
  no duplicate `design_id`.

**No archive extraction is needed** — the binder `sequence` is already a column and the Csy4
RNA target is fixed, so fold inputs are built directly from the manifest. (The 2,000 designs
are scattered over ~2,000 different `designs.tar.zst` archives; extracting them would be
pure waste.)
- RNA target sequence: resolve ONCE from `$RB/targets/` (or an existing refold YAML) and
  record it in the manifest so the run is self-documenting. Verify it is the same 16-nt
  target used by the original refolds before launching.

### 2. Input building
`scripts/build_seedtest_inputs.py` → `inputs/{af3,b2,esm}/`
- One input per design: protein chain **A** (binder sequence) + RNA chain **B** (target).
- Reuse the exact input format the production drivers already consume (Boltz-style YAML for
  B2/ESM; AF3 JSON built the same way `refold_af3.sh` builds it), so the local driver copies
  need no format work.
- MSA: **off** (`msaoff`), matching the existing refolds so results stay comparable.

### 3. Folding — local driver copies in `drivers/`
Copy the three production stage scripts, then make these MINIMAL edits to the copies only:
- `drivers/refold_af3_seedtest.sh` — **no logic change**; invoked with
  `MODEL_SEEDS_JSON="[1,2,3,4,5]"`, `AF3_NUM_DIFFUSION_SAMPLES=1` → 5 structures/design.
- `drivers/refold_b2_seedtest.sh` — add `--seed "$B2_SEED"` to the `boltz predict` call;
  run with `B2_DIFFUSION_SAMPLES=5`, `B2_SEED=1` → 5 structures/design.
- `drivers/esmfold2_infer_seedtest.py` + `drivers/refold_esm_seedtest.sh` — patch the
  writer/canonicaliser to emit **all** `num_diffusion_samples` structures as
  `_ESM_sample{0..4}.cif` (+ per-sample confidence JSON) instead of only `sample0`; run with
  `ESM_NUM_DIFFUSION_SAMPLES=5`, `ESM_SEED=1`.
  - The resume-skip check must key on the LAST expected sample, not `sample0`, or a requeued
    job will think a partially-written design is complete.

### 4. Execution — modelled on `apo_fold_array_v2v7.sh` (proven, inode-safe)
- Chunk the manifest; one SLURM array task per chunk, per backend.
- ALL heavy IO node-local (`$TMPDIR`); each chunk publishes exactly **one packed archive +
  one CSV** to Lustre (keeps the inode footprint flat — see the inode-quota history).
- Copy sbatch directives from a known-good sibling (apo array / `refold_af3_array.sbatch`);
  do not hand-write GPU/partition/QOS lines.
- Isolated from the live v2v7 campaign: own dirs, own job names, throttled arrays (`%N`) to
  respect the 250-job scheduler window.

### 5. Calibration probe BEFORE the full launch
Run ONE small chunk (~20 designs) per backend. Measure wall-time per design at 5 samples and
peak GPU memory (`sacct` TRESUsageInMax). Size chunk count, walltime tier and partition from
the measurement. Do not launch 3 full arrays on estimates.

### 6. Scoring
`scripts/score_seedtest.py` → `seedtest_results.csv`, one row per **design × backend**:
- **Pairwise DockQ** over the 10 unordered pairs of the 5 structures:
  `dockq_pair_mean, dockq_pair_median, dockq_pair_min, dockq_pair_max, dockq_pair_std, n_pairs`.
  Reuse the existing DockQ code path (`dockq_self.py`) rather than a new implementation.
- **Confidence spread** from the per-sample confidence JSONs:
  for each of `iptm, ptm, plddt` → `_mean, _std, _min, _max, _range`.
- Bookkeeping: `n_structures_found` (expect 5), `backend`, `stack`, `design_id`.
- Joinable to the 100k table on `design_id` so seed-stability can be crossed with the
  existing quality metrics.

## Outputs

```
seedtest_2026-07-30/
  seedtest_manifest.tsv          2,000 designs (400/stack)
  inputs/{af3,b2,esm}/           fold inputs
  drivers/                       LOCAL copies (production untouched)
  chunks/                        per-chunk manifests + logs
  outputs/                       one packed archive + one csv per chunk
  seedtest_results.csv           design × backend metrics
  REPORT.md                      findings
```

## Analysis questions the result should answer
1. How large is seed-to-seed spread per backend (pose and confidence)?
2. Which backend is most reproducible?
3. Does seed-stability correlate with quality (iPTM/DockQ from the 100k table)? — i.e. are
   good designs also stable, or is high iPTM sometimes a lucky draw?
4. How much of the existing consensus metrics' spread is sampling noise rather than signal?

## Non-goals (YAGNI)
- No re-scoring of the full interface suite per seed (user chose DockQ + confidence only).
- No change to production drivers or to the existing 100k refolds.
- Not a new campaign — this is a measurement on existing designs, no generation.

## Risks / gotchas
- **ESM multi-sample patch** is the only real code work; if the model's `fold()` returns a
  single structure regardless of `num_diffusion_samples`, fall back to 5 separate seeded runs
  (ESM is ~4 s/design, so this is cheap). Verify in the probe.
- **AF3 output volume**: 5× the CIFs + confidence JSONs per design. Packing per chunk is what
  keeps this inode-safe; do not leave loose files.
- **Walltime routing**: Euler routes by walltime to a tier and can override `--partition`;
  gpu.bulk submits get remapped. Use the documented submit recipe.
- **Do not disturb the running v2v7 apo campaign** — throttle and keep to the public pool.
