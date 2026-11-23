# repredictor_comparison

Why Boltz-2 was dropped. Two bare panels, cited from Section 2.

```bash
~/.venvs/binder_map/bin/python3 code/supp/repredictor_comparison/make_figure.py
```

Emits `output/repredictor_comparison/panel_{a,b}_*.{pdf,png,svg}` — (a) 93 × 96.5 mm
and (b) 52 × 36 mm, **side by side**: 93 + 6 + 52 = 151 mm inside the book's 182 mm
text block, so the two share one row with 31 mm to spare and (b) is top-aligned to (a).
Those are the EXACT emitted boxes: both go out with `style.save(..., tight=False)`,
because a tight bounding box quietly resizes what assembly places. No panel letters, no
titles; both are assembly's job.

## ⚠️ The gate panel was removed, 2026-11-20

Leo: *"Remove panel b. Try to get c next to a and rename it."* The gate cascade — 93 ×
42 mm, the campaign gate applied one term at a time over all 250,000 designs, with the
campaign's own either-backend rule dashed and the three single-re-predictor
counterfactuals coloured — is gone from `output/`. Its last render, the table it drew,
and **an audit of which of its readings survive in other figures and which do not**,
are in `archive/repredictor_comparison_panel_b_gate_cascade_2026-11-20/`.

The short version of that audit:

- **ligand_RMSD does almost all the work** (100% → 1–2% at the first term) survives in the
  main text — `Fig3a_filter_funnel` draws the same cascade campaign-wide with a count
  on every step (6,992,515 → 171,843 → 133,863 → 93,022 → 80,983 → 55,706), and
  `Fig3b_filter_venn` carries the four metric terms as sets with their region counts.
- **Boltz-2's permissiveness is confined to ipTM** — this figure's title claim. **Settled
  2026-11-20, Leo: *"I think we can see it in the 3x3 refolder comparison."* Panel (a) is
  its home, and it does carry it** — checked rather than assumed, on the per-cell
  quadrant counts `report()` already prints:

  | ipTM row, Boltz-2 vs ESMFold2 | only Boltz-2 passes | only ESMFold2 passes | same call |
  |---|--:|--:|--:|
  | **ipTM ≥ 0.6** | **151,668 (60.7%)** | 4,112 (1.6%) | 37.7% |
  | sc_DockQ ≥ 0.7 | 941 | **4,500** | 97.8% |
  | ligand_RMSD ≤ 3 Å | 1,137 | **3,290** | 98.2% |

  Read **down that third column**: the one-sided lobe is in the ipTM cell and nowhere
  else, and in the two rows below it the asymmetry **reverses** — Boltz-2 is the
  *stricter* of the pair on both metrics that read the designed backbone. That is
  "permissive on ipTM, and only on ipTM", drawn, in a single column of the panel. The
  caption now names the cell and the direction.

  ⚠️ **What (a) cannot show, and the caption says so, is the 99.2% / 59.8% itself.**
  Those are cascade quantities — the share of a subset that has *already* cleared
  ligand_RMSD, sc_DockQ and the seq-run term (666 and 2,123 designs) — and an unconditional
  density over all 250,000 cannot express a conditional rate. They stay in the caption
  as **stated, not drawn**, with `report()` deriving them on every run. Everything in the
  caption's new (a) sentence, by contrast, **is** on the panel and is printed by
  `report()`'s per-pair `[both pass … only … neither …]` line, so it cannot drift.
- The **OR/AND quantifier per term**, the **applied order** ending on ipTM, and the
  **three counterfactual chains** are drawn nowhere else either.

**The arithmetic did not go with the drawing.** `GATE`, `chain()`, `campaign_chain()`
and `report()`'s gate block all stay in `make_figure.py`, so every number the panel
carried is still derived and printed on every run and cannot drift out of this README —
the same treatment `bif_effect` gave its round-uniformity χ² when that lost its panel.

**(b) is the old (c), renamed.** The book letters panels by **file-name order**, so
dropping the gate panel moves the seed-ICC box to (b) and the stem has to say so:
`panel_c_seed_icc` → **`panel_b_seed_icc`**. The descriptive half is unchanged, because
what it draws is unchanged — it is still the ICC across the five refolding seeds, per
re-predictor. `make_supp_book.py` lists `panel_b_gate_pass_rate` under `SUPERSEDED`, so
restoring those files into `output/` cannot take the (b) letter back off the panel that
now holds it.

> The supplementary NUMBER is not settled, so it appears nowhere in this folder or in
> the output path. Renumbering is a caption edit, not a rename.

## Sources

| file | mtime | what is taken |
|---|---|---|
| `~/rna_binder_plots/_data/refolder_agreement_designs.parquet` | 2026-11-03 10:47 | 250,000 designs, 50,000 per route, all nine metric columns (`{af3,b2,esm}_{iptm,dockq,sc_rmsd}`) plus the gate flags `g1,g2,g3,g5` and `passes_all` |
| `~/rna_binder_plots/_data/refolder_agreement_stats.csv` | 2026-11-03 10:40 | cross-check only; every value is recomputed from the parquet |
| `~/rna_binder_plots/_data/README.txt` | 2026-11-03 10:48 | §4–5: the gate definition, and per-route generated and Boltz-2-eligible counts for the two reweightings |
| `~/rna_binder_plots/seed_metric_reliability/seed_metric_reliability.csv` | 2026-08-06 14:30 | panel (b): ICC per metric per re-predictor |

**BG+BIF is on the CORRECTED values in every panel.** The parquet's main `dockq` and
`sc_rmsd` columns already hold the true-backbone rescore (`dockq_self_*_true`,
`dockq_LRMSD_*_true`, 2026-10-21); the stored ones survive only as `*_stored` and are
not read here. This matches `Fig3a_backend_agreement`'s `fix_bif()`. It also means
**the chain `report()` prints is not the campaign's funnel**: the campaign gated on the stored BG+BIF
values, so its own chain is 6,992,515 → 55,706 and the corrected chain is → 56,179
(`cascade_by_stack.csv`, variants `old_bif` / `new_bif`). The caption says so.

**Why the sample and not the campaign.** Boltz-2 is in no other local table: the
canonical filter store and `inputs/large/campaign_metrics_full.csv.gz` are AF3 +
ESMFold2 only. It does exist campaign-wide on the Euler master (`b2_iptm`,
`dockq_self_b2`, `dockq_LRMSD_b2`); 6,729,199 of 6,992,515 designs carry a complete set
of all nine metrics, which is the export's eligibility rule, so a campaign-wide version
is possible at the cost of another Euler pass. This figure uses the ready-made
route-balanced export, which means its rates are **deliberately not campaign rates**:
every route contributes 50,000 designs regardless of how many it generated. The two
reweightings printed by the script bracket the campaign-wide value.

## Panels

- **(a) 3 × 3.** Rows are the gate's three metrics, columns the three re-predictor pairs
  in pipeline order. Both axes of a cell carry the row's metric, so the metric is
  written once per row in the left margin and each cell's own axis labels are just the
  two re-predictors — six short labels a row instead of six long ones. Hexbin on a log
  count scale, scaled per cell; dashed = the campaign threshold on both axes, dotted =
  `y = x`. Two numbers per cell, both load-bearing: **R²** (the share of variance one
  re-predictor's score explains in the other) and `same call %` (how often the pair reaches
  the same verdict at that threshold). They disagree, and that disagreement is the
  figure's argument. **Since 2026-11-20 it is also the home of the Boltz-2 claim**
  (Leo: *"I think we can see it in the 3x3 refolder comparison"*): read down the
  Boltz-2 vs ESMFold2 column and the one-sided lobe is in the ipTM cell only, reversing
  in the two rows below — see the audit above for the counts.
- **(b)** Seed reproducibility per re-predictor: one **box** over the ICC(1,1) of 16
  interface metrics across 5 refolding seeds (box = quartiles, rule = median, whiskers =
  full range, no fliers — there are 16 values and no reason to hide any). This is the
  `metric_reliability` figure's panel (a) collapsed to one box per re-predictor; **its n is
  different** (2,000 designs per re-predictor, a different set from (a)). A box hides
  that the ordering is unanimous, so the two counts that carry it are annotated. Was
  `panel_c_seed_icc` until 2026-11-20.

## Encoding

**Palette v0.6, 2026-11-19: re-predictors have their own colours and no hatch.** AF3
`#1b7f79`, Boltz-2 `#2e7d32`, ESMFold2 `#8c6d4f`, from `palette.REFOLDERS` in
`palette.REFOLDER_ORDER` (pipeline order, never by rank) — COLORS.md rule 7. The theme
step + `///` + `...` texture this figure used until today is **withdrawn**: Leo asked
twice for the dots and the stripes to go. No mark in any panel carries a hatch.
**Corrected 2026-11-19:** this paragraph used to end by saying `code/supp/metric_reliability`
(then `seed_reliability`) and `code/supp/interface_chemistry` still keyed re-predictor by
hatch and needed converting. They were converted the same day, so that was false by the
time anyone read it. `grep -rn "hatch=" code/supp/ --include="*.py"` now returns nothing;
the only `hatch=` calls left in the repo are two main-text ones that key a design-mode
class and a "removed" bar, not a re-predictor — see `supp/metric_panel`'s Open list.

**Palette v0.9 (`config/palette.json`, 2026-11-20) changes nothing here — verified by
re-reading the file from disk, not from memory.** v0.9's move is `scaffolds.mod1`–`mod4`
onto theme 600 `#0678a6`; this figure draws no scaffold and never reads
`palette.SCAFFOLDS` (`grep -n SCAFFOLD make_figure.py` is empty). The three re-predictor
hexes in the file are unchanged at AF3 `#1b7f79`, Boltz-2 `#2e7d32`, ESMFold2 `#8c6d4f`,
and `palette.py` parses the JSON at import, so there is no cached copy to go stale.

Panel (a)'s ramp and count scale are `style.density_cmap()` / `style.density_norm()`,
**shared with `Fig3a_backend_agreement`** rather than copied — see below.

## Reconciliation with `Fig3a_backend_agreement`

That main-text figure is AF3 vs ESMFold2 across the same three metrics, so it *is* this
matrix's middle column and the two must not be allowed to disagree.

- **Mechanically**: the ramp and the count scale moved into `code/common/style.py`
  (`density_cmap`, `density_norm`) and both scripts call them, so neither can drift.
  Fig3a's stdout is unchanged by that move, checked line by line.
- **Conventions**: thresholds, axis ranges, tick sets, clipping into the edge bins,
  hexbin gridsize, the `y = x` rule and the annotated statistic are all copied from it.
- **Statistic**: `R²`, not `r`. Fig3a was switched to R² on 2026-11-19 at Leo's
  instruction ("the share of variance one backend's score explains in the other"); this
  panel now follows, so the book annotates one statistic and not two. Literal U+00B2,
  never mathtext (2026-11-16 font audit). `r` and Spearman ρ are still printed by the
  script.
- **Numbers**: they differ only in population, and by very little.

| AF3 vs ESMFold2 | ipTM | sc_DockQ | ligand_RMSD |
|---|---|---|---|
| campaign-wide, 6,992,515 designs (Fig3a) — r | 0.251 | 0.495 | 0.568 |
| campaign-wide — **R²** | **0.063** | **0.245** | **0.323** |
| route-balanced 250,000 (this figure) — r | 0.2545 | 0.4979 | 0.5707 |
| route-balanced — **R²** | **0.065** | **0.248** | **0.326** |

## Verified — every number re-derived by the script and printed on every run

Pairwise, n = 250,000 in every cell:

| metric | pair | R² | Pearson r | Spearman ρ | same call |
|---|---|---|---|---|---|
| ipTM ≥ 0.6 | AF3–Boltz-2 | 0.097 | 0.3121 | 0.2849 | 130,431/250,000 (52.17%) |
| | AF3–ESMFold2 | 0.065 | **0.2545** | **0.2470** | 155,559/250,000 (62.22%) |
| | Boltz-2–ESMFold2 | 0.103 | 0.3211 | 0.3469 | **94,220/250,000 (37.69%)** |
| sc_DockQ ≥ 0.7 | AF3–Boltz-2 | 0.222 | 0.4708 | 0.4956 | 246,308/250,000 (98.52%) |
| | AF3–ESMFold2 | 0.248 | 0.4979 | 0.5437 | 243,748/250,000 (97.50%) |
| | Boltz-2–ESMFold2 | 0.234 | 0.4835 | 0.5354 | 244,560/250,000 (97.82%) |
| ligand_RMSD ≤ 3 Å | AF3–Boltz-2 | 0.275 | 0.5245 | 0.5324 | 246,498/250,000 (98.60%) |
| | AF3–ESMFold2 | 0.326 | 0.5707 | 0.5948 | 244,933/250,000 (97.97%) |
| | Boltz-2–ESMFold2 | 0.355 | 0.5956 | 0.6110 | 245,573/250,000 (98.23%) |

ligand_RMSD on log10 instead: R² 0.258 / 0.295 / 0.300. The sc_DockQ and ligand_RMSD rows sit at
97–99% same-call almost entirely because both re-predictors **fail**; read the quadrant
counts the script prints, not the agreement share, for those two.

Per-predictor, each term and then all three, n = 250,000:

| | ipTM ≥ 0.6 | sc_DockQ ≥ 0.7 | ligand_RMSD ≤ 3 Å | **all three** |
|---|---|---|---|---|
| AF3 (MSA off) | 102,154 (40.862%) | 3,014 (1.206%) | 2,550 (1.020%) | **1,765 (0.7060%)** |
| Boltz-2 | 194,843 (77.937%) | 1,262 (0.505%) | 1,454 (0.582%) | **848 (0.3392%)** |
| ESMFold2 | 47,287 (18.915%) | 4,821 (1.928%) | 3,607 (1.443%) | **1,751 (0.7004%)** |

⚠️ **Those three middle columns are each term on its own.** They do not chain into the
last column and were never meant to: a design can clear sc_DockQ and fail ligand_RMSD.

### The gate AS APPLIED — what the retired panel (b) drew, now `report()` output

The applied order is `g1` ligand_RMSD → `g2` sc_DockQ → `g3` seq run → `g4` foot → `g5` ipTM,
from the export's README section 4, where it is validated against the documented campaign
chain on all five steps. **`g4` is not a column in this export**, so the panel's chain has
four steps and ends at 2,894, which is `passes_all` and is asserted against `g1 & g2 & g3
& g5` on all 250,000 rows.

| chain, applied order | start | ligand_RMSD ≤ 3 Å | + sc_DockQ ≥ 0.7 | + seq run ≤ 5 | + ipTM ≥ 0.6 |
|---|---|---|---|---|---|
| AF3 (MSA off) | 250,000 | 2,550 (1.0200%) | 1,933 (0.7732%) | 1,456 (0.5824%) | **1,339 (0.5356%)** |
| Boltz-2 | 250,000 | 1,454 (0.5816%) | 853 (0.3412%) | 666 (0.2664%) | **661 (0.2644%)** |
| ESMFold2 | 250,000 | 3,607 (1.4428%) | 3,027 (1.2108%) | 2,123 (0.8492%) | **1,269 (0.5076%)** |
| campaign, either backend | 250,000 | 5,612 (2.2448%) | 4,600 (1.8400%) | 3,288 (1.3152%) | **2,894 (1.1576%)** |

All of /250,000. Quantifier per term: **either / either / no backend / either**, with the
missing `g4` being the one **both**-backends term — that asymmetry is on the panel.

**⭐ The ipTM step is the whole argument, and the applied order isolates it.** ipTM is
applied *last*, and it is the only step at which the three re-predictors behave
differently:

| at the ipTM step | before | after | keeps |
|---|---|---|---|
| Boltz-2 | 666 | 661 | **99.2%** |
| AF3 | 1,456 | 1,339 | 92.0% |
| campaign | 3,288 | 2,894 | 88.0% |
| ESMFold2 | 2,123 | 1,269 | **59.8%** |

On ipTM alone Boltz-2 passes **4.12×** ESMFold2 (194,843 against 47,287); the gate's
ipTM term then removes **0.8%** of what Boltz-2 has left and **40.2%** of what ESMFold2
has left. Its permissiveness is confined to the one metric that never reads the designed
backbone, which is why its line is flat at exactly that step. `report()` prints all of it.

### ⚠️ That chain is NOT the campaign's funnel

Three separate reasons, all of which the caption carries:

1. **`g4` (foot ≤ 5, both backends) is not in this export** — footnb was joined only for
   `cascade_by_stack.csv`. The panel's chain is four of five terms.
2. **These are route-balanced sample rates**, 50,000 per route, not campaign rates.
3. **BG+BIF is on corrected true-backbone values** here; the campaign gated on the stored
   ones.

The campaign's own five-step chain, corpus-wide, from `cascade_by_stack.csv` (`ALL`
block) — `new_bif` is the corrected one, `old_bif` is what the campaign actually ran:

| variant | generated | g1 ligand_RMSD | g2 sc_DockQ | g3 seq run | g4 foot | g5 ipTM |
|---|---|---|---|---|---|---|
| `old_bif` (as run) | 6,992,515 | 171,843 | 133,863 | 93,022 | 64,185 | **55,706** |
| `new_bif` (corrected) | 6,992,515 | 166,019 | 136,020 | 95,261 | 64,575 | **56,179** |

The applied gate's JSON lives in the header of `passing_v6_final.csv`. ⚠️ **The gate in
`all_methods.py` gives 66,253 and is NOT the applied one** — do not use it.

Campaign gate on the same sample (each term on *either* gated backend): ligand_RMSD
5,612/250,000 (2.245%), sc_DockQ 7,045/250,000 (2.818%), ipTM 121,941/250,000 (48.776%),
all three **3,961/250,000 (1.584%)** — what the dashed rule in (b) drew. The script asserts that
the export's own `g1/g2/g5` flags reproduce exactly from the corrected columns on all
250,000 rows, which is also what proves the flags are on the corrected BG+BIF values.

Reweighting the all-three rate moves it by under a tenth of a percentage point: by
designs generated AF3 0.761%, Boltz-2 0.365%, ESMFold2 0.758%; by eligible pool
0.78091%, 0.375%, 0.78087%. **Per-route variation is far larger** — AF3 2.120%
(RFdiffusion3+LigandMPNN) to 0.058% (HalluDesign+LigandMPNN), ESMFold2 2.350% to
0.020%. (AF3 and ESMFold2 landing on the same eligible-pool figure is a coincidence of
two different route profiles, not a copy-paste: 0.78091% against 0.78087%.)

Median ipTM: Boltz-2 0.7531, AF3 0.5600, ESMFold2 0.4294.

Seed ICC, 2,000 designs × 5 seeds per re-predictor, 16 metrics:

| | min | Q1 | median | Q3 | max | highest of 3 | lowest of 3 |
|---|---|---|---|---|---|---|---|
| AF3 | 0.1764 | 0.2979 | 0.5929 | 0.7193 | 0.8745 | 0/16 | 3/16 |
| Boltz-2 | 0.1507 | 0.2339 | 0.4717 | 0.7181 | 0.8237 | 0/16 | **13/16** |
| ESMFold2 | 0.2293 | 0.4752 | 0.6470 | 0.8088 | 0.9278 | **16/16** | 0/16 |

Boltz-2's three exceptions are BUNS (optH), contacts backbone and ΔUnsat (optH), where
it beats AF3 but not ESMFold2.

## ⚠️ Corrections this figure has to keep carrying

1. **"Boltz-2 agreed least with either of the other two" is not what the data say.**
   **Say which agreement you mean.** On RANK agreement the weakest pair is **AF3 vs
   ESMFold2** — Pearson r 0.2545, Spearman ρ **0.2470** — and those are the two
   re-predictors the gate actually used. Boltz-2 correlates *better* with each of them
   on both statistics: r 0.3121 / 0.3211, ρ 0.2849 / 0.3469. The Boltz-2 claim is true
   only of VERDICT agreement at the 0.6 threshold, never of ranking. What *is* true, and what panel (a)
   shows, is that Boltz-2's two pairs have the **worst verdict agreement** at ipTM ≥ 0.6
   (52.2% and 37.7% vs 62.2% for AF3–ESMFold2), because its ipTM sits systematically
   higher (median 0.753 vs 0.560 and 0.429), not because it ranks designs differently.
   The three-way shape survives; the sentence that carried it does not.
2. **The ipTM-only story reverses under the whole gate.** On ipTM alone Boltz-2 passes
   4.1× as many designs as ESMFold2 (77.9% vs 18.9%); on all three terms it passes
   **half** as many (848 vs 1,751). Its permissiveness is confined to the one metric
   that never reads the designed backbone. In the gate's APPLIED order ipTM comes last,
   which isolates the effect at a single step: that step keeps **99.2%** of what Boltz-2
   still has and **59.8%** of what ESMFold2 still has. Panel (b) drew it as a line that
   did not fall where the others do; since 2026-11-20 **no figure draws it** — see the
   archive README.
3. **Boltz-2 never faced the campaign gate**, which used AF3 and ESMFold2 only, and its
   coverage of HalluDesign+LigandMPNN is 846,924 of 1,023,007 designs because the v8
   top-up ran AF3 and ESMFold2 without it. Its line in (b) was a comparison constructed
   here — as were AF3's and ESMFold2's, since the campaign gated on the two together and
   never on one alone. Only the dashed line was a rule the campaign applied.

## Open

- **Panel (b) is interface-metric reproducibility, not ipTM reproducibility.** No
  ipTM-per-seed table exists: `seed_metrics_merged_v2.csv` carries the interface metrics
  only, with no ipTM column, although `metric_std_vs_quality.csv` has `rho_cv_iptm`, so
  the per-design ipTM existed when that analysis was run. **NEED** a per-seed ipTM table
  if the panel is to be about ipTM directly.
- **COLORS.md rule 7's validation figures do not reproduce here.** Rule 7 claims the
  re-predictor set sits at ΔE ≥ 45 from `MAIN` and from theme 400 and that its own worst
  internal pair is ΔE 29.3. Measured under plain CIE76: AF3 `#1b7f79` is **ΔE 32 from
  `MAIN`** (theme 700) and 36 from theme 400 — it is the closest re-predictor to the house
  blue — while the worst internal pair is AF3–Boltz-2 at 39.8. The same implementation
  reads the design-route palette (`palette.STACKS`) worst pair as 23.6 against the
  documented 16.5, so the two
  scales differ by ~1.4×; AF3 is under 45 either way. Practically this is cross-panel
  only — (a)'s blue means *count*, (b)'s colours mean *re-predictor* — but the palette
  is meant to be the single source of truth, so **this is for Leo to settle centrally**,
  not for this script to work around.

## Caption draft

**Boltz-2 differs from the other two re-predictors in level, not in ranking: most
permissive of the three on ipTM, least under the gate as a whole.** Both panels use a
route-balanced sample of 250,000 designs, 50,000 per design route, so these are sample
rates, not campaign rates; BG+BIF carries its corrected true-backbone sc_DockQ and
ligand_RMSD, whereas the campaign's funnel was gated on the stored values. (a) Rows are the
gate's metrics, columns the re-predictor pairs; both axes of a cell carry the row's metric.
Density on a log count scale, scaled per cell; dashed = the campaign thresholds
(ipTM ≥ 0.6, sc_DockQ ≥ 0.7, ligand_RMSD ≤ 3 Å), dotted = y = x; R² and the share of designs on
which the pair reaches the same verdict. The pairs rank designs almost equally weakly on
ipTM (R² 0.06–0.10) and differ in level, which drives Boltz-2's two pairs to the worst
verdict agreement, 37.7% and 52.2% against 62.2% for AF3–ESMFold2. **The permissiveness
is in one cell and the restraint in the two below it:** in the ipTM row's Boltz-2 vs
ESMFold2 panel, 151,668/250,000 (60.7%) designs clear 0.6 on Boltz-2 and fail it on
ESMFold2 — the dense lobe right of the vertical threshold and below the horizontal one —
against 4,112 (1.6%) the other way; in the sc_DockQ and ligand_RMSD cells of that same
column the asymmetry **reverses** (941 against 4,500, and 1,137 against 3,290) and the
pair agrees on 97.8% and 98.2%. Boltz-2's level difference is confined to the one gate
metric that never reads the designed backbone. The middle column is
the route-balanced version of Fig. 3a. (b) ICC over 5 refolding
seeds for 16 interface metrics, on a separate 2,000 designs (box = quartiles, whiskers =
full range): ESMFold2 most reproducible on 16 of 16 metrics, Boltz-2 least on 13 of 16.
**What (a) shows and what is stated.** The panel carries the mechanism — the ipTM
level difference, its one-cell location and its reversal in the two rows below. The
consequence under the gate is **stated here and drawn nowhere**, because it is a
cascade quantity conditional on a design having already cleared ligand_RMSD, sc_DockQ and
the sequence-run term, which an unconditional density over all 250,000 cannot show:
applied one term at a time, ipTM comes **last** and is the only step at which the three
differ, keeping **99.2%** of what Boltz-2 still has against **59.8%** of ESMFold2's, so
all three finish within a factor of 2.6 of each other (AF3 1,339/250,000, ESMFold2
1,269, Boltz-2 661) although Boltz-2 passes 4.12× as many designs as ESMFold2 on ipTM
alone. `report()` prints the whole chain on every run; the panel that used to draw it is
in `archive/repredictor_comparison_panel_b_gate_cascade_2026-11-20/`.

## Naming, settled 2026-11-20

**`sc_DockQ` and `ligand_RMSD (Å)`**, matching the main text (Fig 2c), wherever these two
quantities are an axis label or appear in prose. The bare `DockQ` and the hyphenated
`ligand-RMSD` are not used as display forms anywhere in `code/supp`.

Renamed **`sc_RMSD` -> `ligand_RMSD` on 2026-11-21**: the gate's first term is the
RNA-aligned *ligand* RMSD (`dockq_LRMSD_af3` / `_esm`), a pose metric, not the Ca fold
RMSD it was once written up as. Same measurement throughout, new name.

⚠️ **CSV column names are untouched and stay untouched** — `sc_dockq`, `sc_rmsd`,
`dockq_af3`, `dockq_esm`, `rmsd_af3`, `rmsd_esm`, `dockq_self`, `dockq_LRMSD`,
`dockq_fnat`, `dockq_iRMSD` and the `_true` variants all keep their lower-case spelling,
because they are what is on disk. The sweeps only ever matched capitalised display forms,
never a column name.

## Vocabulary, settled 2026-11-19, swept 2026-11-20

**"re-predictor"** for AF3 / Boltz-2 / ESMFold2, never "refolder", "predictor" or
"repredictor". **"refolding seeds"** for the five seeds. **"design route"**, never
"stack". In all prose, labels and legends.

⚠️ **Identifiers are deliberately exempt. Do not "fix" them — each one breaks a read or
an import:** the `refolder` and `stack` columns in the source CSVs, `palette.REFOLDERS`,
`palette.REFOLDER_ORDER`, `palette.STACK_ORDER`, `palette.SHORT`, the source file names
that contain `_by_stack_refolder`, and the directory name `Fig5fg_h1null_refolder`.
