# Union Fold Map Explorer — Design

Date: 2026-09-28

## Purpose

Show the campaign sample, the filtered pool and the ordered library **in one shared
fold space**, so the narrowing measured by `compare/` on 2026-09-28 can be seen
rather than read off a table. Three independent UMAPs would not support this: their
coordinate systems are unrelated, so a point landing in the same place in two of them
would be coincidence. The three node sets are therefore **co-embedded** into a single
map, and the page draws each as a toggleable layer.

This is a visualization. It produces no selection, no gate, and no new comparison
statistic.

## The node sets — measured 2026-09-28, not estimated

| set | source tree | nodes |
|---|---|---|
| `sample500k` | `results_full` | 499,999 |
| `passing60k` | `results_sw_k30_f50_mut` | 31,763 |
| `library12k` | `results_lib12k` | 9,759 |
| **union** | new `results_union` | **526,467** |

Pairwise, by `coord_hash`:

| intersection | nodes |
|---|---|
| sample ∩ passing | 5,295 |
| sample ∩ library | 1,553 |
| passing ∩ library | **9,759 — the library is entirely inside the passing pool** |
| all three | 1,553 |

Two cross-checks passed: 5,295 equals the `n_passing` count stage 9 derived
independently, and `sample ∩ library` equals `overlap` in stage 0's `sample.json`.
The union arithmetic closes exactly:
499,999 + 31,763 + 9,759 − 5,295 − 1,553 − 9,759 + 1,553 = 526,467.

**"Conserved backbones" means the intersections above** — a backbone physically
present in more than one set, detected by coord-hash identity when stage 2 collapses
the union. It does not mean recurring folds, cluster-level overlap, or funnel
survival inferred from flags.

## Decisions taken during design

1. **Co-embed the union**, rather than projecting into the existing v6 layout or
   toggling two unrelated maps. Overlay is then spatially meaningful.
   **Consequence to accept:** the familiar v6 picture gets new coordinates. The
   existing `dist/binder_map_v6.html` is not modified and remains the old view.
2. **Fold axis only.** Revised 2026-09-28 — the original decision was both axes.
   The fold axis is what the comparison asks about, and the UMAP at 526k nodes is the
   riskiest step in the build, so epitope is dropped to halve it. The 48-dim contact
   fingerprints are still assembled into `results_union/design/fp/`, so the second axis
   can be added later without redoing the search.
3. **Points plus click-to-ChimeraX.** The page carries coordinates, cluster ids and
   layer flags only; clicking calls the existing local REST server on `:45000`.
   Embedding 3D traces as v6 does would put the page near 400 MB at this scale.

## Pipeline — no new extraction

Every structure needed is already on disk: `results_shared/design/shards` holds the
campaign sample and the library, `results_v6/design/shards` holds the passing pool.

| stage | what | notes |
|---|---|---|
| U0 | assemble `results_union` shards + hashes from both source trees | copy or link; no archive reads |
| 2 | dedup + `createdb` | 526,467 nodes; the coord-hash collapse detects conservation |
| 3 | search, 512 chunks | target DB staged to node-local scratch (the 2026-09-25 fix) |
| 4 | graph, fold axis | frozen parameters |
| 5 | Leiden, fold axis, r05…r80 | `cluster_r80` operative |
| 6 | UMAP, fold axis | the expensive and only genuinely risky step |
| U8 | build the page | new builder, modelled on `stage8_build_html.py` |

**Frozen parameters carry over unchanged:** `--alignment-type 1`,
`-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a`, `TM = min(qTM,tTM)`,
`--tm-floor 0.5`, `--epi-floor 0.6`, `--k 30`, `--mutual`, Leiden fixed seed.
Not because the union must match map 2 numerically — it cannot, being a different
node set — but so the picture is built by the same rule as everything else here.

## Layer flags

Each union node carries `in_sample500k`, `in_passing60k`, `in_12k`, assigned by
joining each source tree's `dedup/nodes.tsv` on `coord_hash`. The counts in the
tables above are the expected totals and must be asserted, not trusted.

## The page

- **Layer toggles:** campaign / passing / library, independently on and off.
- **Conserved marking:** nodes in more than one active layer drawn distinctly, with
  the count of what is currently shown.
- **Rendering:** canvas or WebGL. 526,467 points must pan and zoom smoothly; SVG
  will not. Coordinates packed binary (base64 `Float32Array`), not JSON per point.
- **Hover:** design id, stack, union cluster id, layer membership.
- **Click:** POST to the local ChimeraX server on `:45000`; degrade quietly to
  hover-only when it is not running.
- **Self-contained** single file under `dist/`, target ≤ 25 MB.

## What this is not

**The union clustering is a fourth, new clustering.** The comparison in `compare/`
used three independent per-map clusterings by design; cluster ids here are union
cluster ids and are not comparable to it. The page must say so on its face, so no
one quotes an explorer cluster count as a comparison result.

## Risks

- **UMAP at 526k nodes is the failure point.** The kNN block size already had to
  shrink once for memory at this scale (4096 → 1024, an 8.3 GB allocation). Size the
  embed job from a measured run, not an estimate, exactly as stage 3 was sized.
- **Page weight.** 526k points × 2 axes must stay inside the size target; if it does
  not, the campaign layer is decimated for display and the page says by how much.

## Out of scope

- Rebuilding or modifying `dist/binder_map_v6.html`.
- The af3 and esm layers. Design layer only.
- Any change to `compare/`, the 12k library, or the three existing maps.
- Medoids, traces and the complex viewer that v6 carries.

## Validation

1. Union node count is exactly 526,467, and the three layer-flag totals match the
   measured table. Hard asserts.
2. Every node carries at least one layer flag; none carries zero.
3. Conserved counts reproduce: 5,295 / 1,553 / 9,759 / 1,553.
4. The fold embedding has one row per node — no silent drops between stage 5 and 6.
5. The page opens and renders at 526k points without exceeding the size target.

## Known hazard carried in from the source trees

`results_shared/design/hashes` holds **508,207 rows**, but its shards hold fewer
structures: stage 1 records a hash row for every design and writes only one
structure per coord_hash (`stage1_extract_nested.py:200` vs `:214-216`). The
per-map trees were reconciled on 2026-09-25 by filtering their hashes to the
structures actually written — but **`results_shared` itself was never filtered**,
because the fix lives in `stage1b_split.py`, which the union assembly does not run.

Assembling the union straight from `results_shared/design/hashes` therefore
reintroduces exactly the bug that killed the 2026-09-24 run: stage 2 elects a
canonical design whose structure is not in the DB, the name join comes up one
short, and `stage2_builddb.sh` aborts with
`FATAL: N canonical structures are not in the Foldseek DB`.

**Stage U0 must filter the assembled hashes to design_ids present in the assembled
shards**, the same way `stage1b_split.py` does, and record what it dropped. The same
applies to `results_v6`'s hashes, which predate that fix entirely.

This is the one place where the union build can fail late and expensively, so it is
checked at assembly time, not left to stage 2 to discover.
