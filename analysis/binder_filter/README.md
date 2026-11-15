# binder_filter

Interactive filter bench for the Csy4 binder campaign. Move a threshold, watch the yield and
its composition move with it, then shift-click into ChimeraX to look at what you selected.

```bash
~/miniconda3/bin/python ~/Documents/binder_filter/app.py     # http://127.0.0.1:8010
openrefolds                                                   # ChimeraX bridge, for clicking
```

## Boltz-2 is excluded — and that is a data trap, not just a policy

Only AF3 and ESMFold2 values drive a filter. This matters more than it sounds: **every bare
consensus column in the master CSV folds B2 in** — `consensus_iptm`, `iptm_min3`,
`iptm_best2`, `dockq_self`, `dockq_LRMSD`, `n_h_bonds_interface`, `pi_stack_score`, `red_*`,
`dsasa_*`. Reading any of them would silently contaminate the tool.

Three layers enforce the exclusion:

1. `build_filter_store.py` refuses to map a column that is a known bare consensus name, or
   whose backend token does not match its tool key (`_assert_no_consensus`).
2. The store has no b2 column, so nothing downstream can reach one.
3. `filters.BACKENDS` is `("af3", "esm")`; asking for `b2` raises.

`iptm_min` and `iptm_mean` are recomputed in `store.py` from the two backends.

## The two plots

| plot | x | y |
|---|---|---|
| sc_RMSD | `dockq_LRMSD_af3` | `dockq_LRMSD_esm` |
| sc_DockQ | `dockq_self_af3` | `dockq_self_esm` |

Both are in the **RNA-aligned** (target-fixed) convention that became the project default on
2026-08-04. `dockq_std_*` holds the stock-convention values and is not used here.

⚠️ **sc_RMSD here is pose RMSD, which is a DockQ input.** The two plots therefore move
together — they are not independent evidence. Fold RMSD (Cα binder-on-binder vs the designed
backbone) would be orthogonal, but has never been computed campaign-wide per refolder.

## Refolder agreement — `pose_disagree` / `fold_disagree`

Both plots measure each refolder against the **design**. These two tier-3 metrics measure the
refolders against **each other**, which is a different question and not derivable from the
plots — the triangle inequality only bounds it:

```
| rmsd_af3 − rmsd_esm |  ≤  pose_disagree  ≤  rmsd_af3 + rmsd_esm
```

| metric | source | asks |
|---|---|---|
| `pose_disagree` | `rmsd_pose_af3_esm` | superimpose both refolds on the RNA, apply that transform to the binder, Cα RMSD with **no second fit** — do they agree **where** it binds? |
| `fold_disagree` | `rmsd_prot_af3_esm` | binder-on-binder Kabsch fit — do they agree **what it looks like**? |

Neither has a backend split — they *are* the comparison — so they take no `_af3`/`_esm` suffix
and sit out the agreement panel. They are orthogonal to everything else in the tool
(|ρ| ≤ 0.31) and to each other (ρ = +0.23).

⚠️ **How much they add depends entirely on how you gate**, because that upper bound does real
work. Measured on the full store:

| gate | passing | + `pose_disagree ≤ 10 Å` | kept |
|---|---|---|---|
| sc_DockQ ≥ 0.5 **AND** (both refolders) | 145,140 | 141,097 | 97.2% |
| sc_DockQ ≥ 0.5 **OR** (either refolder) | 635,877 | 263,022 | **41.4%** |

Requiring both refolders near the design already forces them near each other, so on a strict
AND shortlist this is nearly a no-op. Under OR gating it removes 59% — designs where one
refolder liked it and the two disagree completely about where it sits. That is what makes
per-metric OR in tier 1 safe to use.

Campaign-wide the refolders rarely agree at all: median `pose_disagree` **22.2 Å**, 85% > 10 Å,
0.35% within 2 Å.

## Filter model

Atoms are `(metric, backend, cmp, value)`. Atoms combine inside a **group**, groups combine
inside a **tier**, tiers combine with a top-level AND/OR — an op at every level:

```
Tier 1   ┌ sc_RMSD    [AND]  af3 ≤ 2.0 · esm ≤ 2.0
         │  ── AND ──                                  <- tier op, combines the groups
         └ sc_DockQ   [OR]   af3 ≥ 0.50 · esm ≥ 0.50
```

Tier 1 puts **each metric in its own group**, so sc_RMSD can demand both refolders agree while
sc_DockQ accepts either. One op for the whole tier could not express that: over those four
atoms it can only reach 2,269 designs (all AND) or 635,898 (all OR), while the mixed settings
land anywhere in between. Tier 3 has no groups and still reads as a single tier-wide op —
`tier_groups()` normalises it to one group, so both shapes evaluate identically.

Three rules that are easy to get wrong:

- **NaN fails closed.** A missing value fails any filter on that metric — never silently passes.
- **An empty group or tier is neutral**, under both AND and OR. Unticking the last box in an OR
  tier must not zero the result.
- **Neutral absorbs under OR, at every level.** "No constraint OR anything" is no constraint, so
  an empty group inside an OR neutralises its whole tier. Correct, but it reads as a broken
  filter, so `evaluate` returns `neutralized` / `unfiltered_or` and the UI names the cause.

The ipTM gate sits between the tiers and runs in one of two modes: a plain threshold, or a
rank cut (top N / top X%) scoped either to the filtered set (yield is exactly N) or globally.

## Counting discipline

**Every count and percentage comes from the full store; the display subsample is only there so
the browser can draw something.** The subsample is drawn once with a fixed seed, so moving a
threshold recolours points rather than redrawing them.

## Data gaps

| metric | state |
|---|---|
| globularity (`rog_af3` / `rog_esm`) | all-NaN placeholder — the filter greys out. Needs a structural pass over the campaign. |
| refolder agreement (`pose_disagree` / `fold_disagree`) | complete, 5,343,523 / 5,343,524. v9 will need `pairrmsd_chunk.sbatch` rerun on its chunks. |
| π-stacking | weighted **score**, not a count. The per-nucleotide count exists for 100k designs only and is not recoverable from the score. |
| v9 (~1.19M designs) | not in the master yet. Rebuilding with a new `--csv` is all it takes. |

Backfilling a placeholder is a pure data swap: rebuild the store, no schema or UI change.

## Layout

| file | role |
|---|---|
| `store.py` | columnar load + CIF path reconstruction |
| `filters.py` | the filter-tree evaluator (the unit worth testing hardest) |
| `app.py` | Flask API; also proxies the ChimeraX bridge so the page never goes cross-origin |
| `export.py` | CSV + CIF export; borrows the hardened ssh extraction from `metric_explorer` |
| `templates/index.html` | the UI |
| `test_filters.py`, `test_store.py` | plain-assert tests — `~/miniconda3/bin/python test_*.py` |

Built on Euler by `code/scripts/build_filter_store.py` (via `ops/submit/build_filter_store.sbatch`).
Cross-checked by `code/scripts/verify_filter_counts.py`, which re-implements the filter
semantics independently over the raw CSV.

## Path compression

The store keeps no path strings. Each row carries a template index; the path is rebuilt by
substituting `<ID>` and `<U>`. The builder asserts a byte-exact round-trip on **every** row —
worth keeping, because the campaign contains both packed (`…tar.zst!designs/…`) and loose
(`…/designs/…`) layouts, and a prefix assumption would have silently corrupted the loose ones.

The design_id *contains* the uuid8, so substitution order matters: design_id first, then uuid8.

## ChimeraX

Clicking opens **design + AF3 + ESM**, superposed on the RNA (chain B), with the designed
backbone as model 1.

- shift-click a point → that design
- shift-click empty space → 20 random designs from that region, pre-extracted, then step with
  `nextd` / `backd` (`whered`, `gotod`, `markd` also available)

`DesignCycler.py` lives in `~/Desktop/MA/Data/ChimeraX/CustomPreset/` next to `PDBcyclerTools.py`
and keeps its own session state, so it coexists with `next` / `back` / `skip` rather than
fighting them.
