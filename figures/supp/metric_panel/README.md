# metric_panel

The campaign-scale metrics the main text does **not** plot, per design route, plus the
one stored-vs-corrected pair that exists. One bare artefact — a 2 × 3 grid of six
unshown metrics over a two-panel strip — cited from Section 2. Every body is a **split
violin** in main-text Fig 2c's idiom: AF3 (MSA off) on the left, ESMFold2 on the right.
The one exception is apo RMSD, which has a single re-predictor and so a single whole
body — see **Split violins, 2026-11-20**.

```bash
~/.venvs/binder_map/bin/python3 code/supp/metric_panel/make_figure.py
# ~9 min, reads a 259 MB, a 233 MB, a 22 MB and a 20 MB gz. While iterating:
METRIC_PANEL_CACHE=/tmp/mp2.parquet ~/.venvs/binder_map/bin/python3 code/supp/metric_panel/make_figure.py
```

Emits `output/metric_panel/metric_panel.{pdf,png,svg}` at the authored 178 × 145 mm
(`style.save(..., tight=False)` — every axes is placed in mm, so a tight bounding box
would crop to the ink and quietly resize what assembly places). No letters, no titles.
Grid panels are 44 × 33 mm on an 8 mm row gap; the strip's two are **72 × 33 mm**, which
puts its right edge at 13 + 2×72 + 12 = 169 mm, the same as the grid's 13 + 3×44 + 2×12,
so the two blocks share a right margin while the strip's six positions get 12 mm each
against the grid's 8.8. `main()` asserts that agreement rather than trusting the
arithmetic.

**Both blocks of rotated x labels are measured, never estimated**
(`figure-measure-dont-estimate-text`). The grid's longest, "RFD3+LiMPNN", is 13.11 mm at
the house 5.5 pt and needs 15.05 mm with its tick and pad; it now hangs into the 18 mm
gap above the strip. The strip's longest is **"BG+BIF corrected" at 15.30 mm, needing
17.24 mm**, which is why the bottom margin is **18.5 mm and not the 16.5 that held the
grid** — the pair's labels are longer than any route name, and with `tight=False` a
margin a fraction short clips the label instead of being absorbed by a tight bounding
box. `main()` asserts each margin against `style.text_w` of every label that lands in
it, and asserts the three rows fit the canvas.

Terminology, repo-wide: **re-predictor**, never "refolder" — AF3, Boltz-2 and ESMFold2
predict a structure from sequence rather than refolding a given one (COLORS.md rule 7) —
and **design route**, never "stack", for a generator + inverse-folder combination. CSV
column names keep their historical spellings (`stack`, `bif_true_per_backend`).

> The supplementary NUMBER is not settled and appears nowhere in this folder or in the
> output path.

---

## Rebuilt 2026-11-19, on Leo's comment

> "most of these metrics we already have in the main panel. can you show the metrics we
> have not shown? Im not sure which ones are left. I think we have a lot like apo_rmsd
> or delta sasa etc. buns etc. Show the old bg+bif wherever available and different from
> the correction"

He was right, and the audit below is the answer to "which ones are left". The previous
cut drew ipTM, sc_DockQ, ligand_RMSD, buried ΔSASA total, Rg and foot clearance. **Five of
those six are in the main text** — only `foot_tp_min_dist` was not, and even that is the
distance twin of a gate term the main text draws as an atom count. The figure is now
built from what the audit left.

## Split violins, 2026-11-20

> "can you do the split violins for metric_panel the same way as its done in the main in
> fig2"

Done, in **`code/Fig2cd_benchmark/make_violins.py`'s idiom**: every body with two gate
re-predictors behind it is now one violin split down its centre line — **AF3 (MSA off)
the left half at full strength, ESMFold2 the right at alpha 0.5**, both in the design
route's own colour, with one key in the grid's row gutter. Fig 2c splits top/bottom
because its violins lie on their side; these stand up, so the same split runs left/right
and AF3 still reads first.

**What that changes is the quantity on the canvas.** The panels used to draw this figure's
worst-of-two **consensus**; they now draw its two **inputs**. The consensus rule is
unchanged, still asserted against a column somebody else computed before anything is
drawn, and still what `report()` tabulates and what the verified medians below are — so it
left the canvas and nothing else. `report()` now prints the **halves' medians beside it**,
because a reader checking a median off the panel needs the number that is drawn.

**Apo RMSD is NOT split, and that is Fig 2c's own first-panel decision.** There, "Designs"
is a bar because it is a count with no distribution behind it and a violin would invent a
spread. Here, apo RMSD is **one ESMFold2 measurement** — the binder folded alone, against
the ESM complex prediction — with no AF3 counterpart in existence, so a split would invent
a re-predictor. It keeps one whole symmetric body, and its n note now reads
`ESMFold2 only · n = …` so the canvas says why.

**The strip's pair is untouched in kind.** `BG+BIF stored` and `BG+BIF corrected` are still
two positions in the route colour, **stored open and corrected filled**, side by side in
BoltzGen+BoltzIF's own slot — and both of them split, because both sides have
per-re-predictor columns (stored in the campaign's own `dockq_/rmsd_{af3,esm}`, corrected
in the rescore's `*_true` twins). ⚠️ **A tone cannot show on an unfilled body**: alpha on an
outline leaves the ESMFold2 half a half-visible hairline, so `body()` deliberately does not
apply it there and the **half itself** carries the re-predictor on the stored position. That
is the cue the grid leans on too — in Fig 2c as here, the alpha is the redundant one.

⚠️ **An open body is drawn as its density CURVE, not as a closed half** — changed after
looking at the first render. `fill_betweenx` closes each half along the centre line, which
on a filled body disappears under the fill and on the strip's open stored body came out as
the loudest mark in the panel: a full-height stroke bisecting the one body whose two halves
are nearly identical, crossed by its median rules into a plus sign. Two curves meeting at
the trim ends read as the hollow violin the position had before the split, with each side
now traced on its own values. At `TRIM` the curve is 0.5% of the half width off the centre,
so the ends close themselves.

**Fig 2c's `_half()` is not imported, though the idiom is its**, and that is a considered
no rather than an oversight: it draws horizontally (`vert=False`), writes to module-global
collectors for its own alignment check, and its density is matplotlib's KDE over a
**60,000-design subsample per route**. This figure's density is a hand-rolled histogram
over **every value** — 1.4 M per route, which is what makes the KDE unaffordable and the
subsample necessary there. Importing it would have changed what six panels draw in order
to share nine lines. What **is** shared is the convention, and it is declared **once**, in
`BACKENDS`: which half, which tone, which name, used by the grid, the strip and the key
alike, so a panel cannot draw the halves in one order and the key explain them in another.
(Correction while checking this: the module docstring said matplotlib's `violinplot` was
ruled out because it calls `scipy.stats.gaussian_kde`. It does not — it uses its own
`mlab.GaussianKDE` and works fine in this venv. The histogram stays for the reasons above,
which are the real ones; scipy is still unusable here, which is why every statistic is
hand-rolled.)

**The key.** Two touching swatches in neutral `palette.INK` at the two alphas, each
re-predictor's name beside its own side, left-aligned in the **8 mm gutter between the
grid's two rows** — the only empty full-width band on a canvas whose three rows and both
blocks of rotated labels are measured to the millimetre. Measured at **35.10 mm** wide
against the 100 mm of two panels, and asserted. It is neutral and never a
`palette.REFOLDERS` hue for the same reason Fig 2c's is: every body here is a design route
in its route colour, and a re-predictor colour beside them would read as a sixth route
(COLORS.md rule 7). The gutter's **right** end is not free — the apo panel's n note hangs
there, at axes y = 1.04 — which is why the key is left-aligned.

## THE AUDIT — every metric the main text plots, and where

> ⚠️ **Two rows of this table are drawn here anyway**, and that is not an oversight:
> sc_DockQ and ligand_RMSD are re-admitted as the strip, because they are the only two metrics
> with a BoltzGen+BoltzIF corrected counterpart in existence. See **THE STRIP** below;
> the caption says it on the page. The audit still governs the six grid panels, and
> nothing else is re-admitted.

Read on 2026-11-19 from `code/Fig2*`–`code/Fig5*` and `captions/Fig1.md`–`Fig5.md`.
**Promoted** = on the page in `selection/FigN.svg`. **Built** = the script and its panel
exist and `FIGURES.md` lists it as a main-text panel, but the current assembly does not
place it — for Fig 5 d–f that is a 2026-11-18 change (`code/Fig5_assembly`: "with d-f
dropped") that `FIGURES.md` has not caught up with. Treat **built** as claimed by the
main text: it may come back.

| quantity | column | where | status | form | n |
|---|---|---|---|---|---|
| ipTM | `af3_msaoff_iptm`, `esm_iptm`, `consensus_iptm` | Fig 2c, Fig 3a, Fig 3b | promoted | violin per re-predictor; hexbin; Venn term | 6,992,515 |
| sc_DockQ | `dockq_self` | Fig 2c, Fig 3a, Fig 3b | promoted | same | 6,992,515 |
| ligand_RMSD | `dockq_LRMSD` | Fig 2c, Fig 3a, Fig 3b | promoted | same | 6,992,515 |
| Rg | `radius_of_gyration_ca` | Fig 2c | promoted | violin per re-predictor | 6,992,515 |
| Rg, elongation, helix fraction, longest helix, interface residues, interface atom pairs, nucleotides touched, base contact fraction, contact entropy, K+R at interface | AF3 re-prediction columns | Fig 4h | built | percentile strip + densities | 32 vs 12,000 |
| Designs generated; designs per GPU-hour | — | Fig 2c, Fig 3f | promoted | bars; scatter | 5 routes |
| Gate pass rate | — | Fig 3f | promoted | scatter | 5 routes |
| Longest same-residue run | gate `seq_run` | Fig 3b | promoted | Venn set | 6,992,515 |
| Foot clash, **atom count** | gate `footnb_010` | Fig 3b | promoted | coda bar | 80,983 → 55,706 |
| SSE α / β / loop fraction | 100k `sse_*` | Fig 3c | promoted | ternary density | 31,763 backbones |
| Per-nucleotide contacts, split phosphate / sugar / base | measured on design models | Fig 2d, Fig 4e | promoted | bars | 614–800 per route; 497,007 nodes |
| Binder occupancy density | — | Fig 2d | promoted | volume renders | 800 per route |
| π-stacking distance, groove split, moiety split | measured per call-out | Fig 2b | promoted | annotation on 8 designs | 8 |
| Backbone-cluster coverage; pool composition; hit rate | — | Fig 3d, Fig 3e, Fig 4d | promoted | cells; pies; bars | 917 clusters; 12,000; 32 |
| Splicing ratio, log₂ fold change; ON/OFF | — | Fig 4b, Fig 5b | promoted | scatter; ranked dots | 16,954; 32 |
| **ipTM, sc_DockQ, ligand_RMSD, interface ΔSASA total, heavy contacts, interface H-bonds, shape complementarity, packstat, buried unsatisfied (as `buns_vsasa`, the WHOLE-POSE column — see the BUNS box), π-stacking** | `dsasa_total`, `n_heavy_contacts_4A`, `n_h_bonds_interface`, `sc_value`, `packstat`, `buns_vsasa`, `pi_stack_score` | Fig 5d (`Fig5def_asn_interface/metrics`; the 10-metric `Fig5fg_h1null_refolder` variant is the same list per re-predictor) | built | P(binder > inactive) dots with bootstrap intervals | 15 vs 17 designs |
| 19 amino acids, sequence and interface share; interface Asn | `if_*` | Fig 5e, Fig 5f | built | volcano; densities | 15 / 17 / 600 / 59 |

**This is why BUNS and ΔSASA-total count as "shown" but are still worth drawing.** Fig 5d
shows them as a **separation statistic on 32 designs** — one dot per metric, is it
predictive. Nothing anywhere shows their **campaign distribution per design route**.
Different quantity, same metric name. The caption says so.

### Drawn in a sibling supplementary, so not re-drawn here

| quantity | figure | n |
|---|---|---|
| Base / backbone / phosphate contacts and H-bonds against 58 natural RBP–RNA complexes | `code/supp/natural_rbp_benchmark` | 58 naturals vs designs |
| The five gate terms, gated pool vs ordered panel, as ECDFs | `code/supp/metric_significance` | 55,706 vs 12,000 |
| Boltz-2 against the two gate re-predictors | `code/supp/repredictor_comparison` | — |
| The BoltzIF defect, stored vs corrected | `code/supp/bif_effect`, `code/supp/bif_quasicontrol` | 1.48 M; 2,397 |
| Precision–recall of ten of the Fig 5d metrics | `code/Fig5pr_precision_recall` (built, no page letter) | 15 vs 8,849 |

⚠️ **`code/supp/interface_chemistry` left this table on 2026-11-20.** Archived whole
(Leo: *"Archive for now"*) to `archive/interface_chemistry_retired_2026-11-20/`, so it is
no longer a sibling that draws anything. It carried **interface H-bonds per design, base
share of H-bonds, and π-stacking per design route and per nucleotide**, 20,000 per route ×
3 re-predictors — and by the time it went, only the per-nucleotide π-stack panel was still
drawn there; (a), (b) and (c) had been deleted earlier the same day. **None of those
quantities is drawn anywhere now.** The two cut rows below that leaned on it
(`n_h_bonds_interface`, `pi_stack_score`) are rewritten accordingly: both verdicts stand
on their own measured grounds, which is why neither metric comes back here — the sibling
was a second reason, not the reason.

### What is LEFT — unplotted anywhere, at campaign scale

All on the Euler master `analyses/cross_campaign/csy4_master_foot.csv`
(14,553,905,627 bytes, 2026-09-01 13:35, 6,992,515 rows × 218 columns), per re-predictor
`_af3` / `_b2` / `_esm` plus a stored 3-re-predictor consensus.

| family | columns | drawn here |
|---|---|---|
| Rosetta packing | `sc_value`, `packstat` | **packstat** |
| Buried unsatisfied polars | `delta_unsat_hbonds`, `buns_vsasa`, `delta_unsat_hbonds_opth`, `buns_vsasa_opth`, `n_polar_h_moved` | **delta_unsat_hbonds** |
| ΔSASA, protein side | `dsasa_polar`, `dsasa_hphobic` | no |
| ΔSASA, RNA side | `dsasa_rna_total`, `_backbone`, `_phosphate`, `_sugar`, `_o2prime`, `_base`, `_base_major`, `_base_minor`, `_base_wc`, `_base_unassigned`, `frac_dsasa_base`, `frac_dsasa_major` | **frac_dsasa_base**, **frac_dsasa_major** |
| Contact counts | `n_heavy_contacts_4A`, `_4p5A`, `n_contact_res_prot`, `n_contact_res_rna`, `n_contacts_base`, `_backbone`, `_phosphate`, `_sugar`, `frac_base_contacts` | no |
| Unfavourable contacts | `red_total_per_ifaceres`, `red_clashes_`, `red_acid_phos_`, `red_hphob_polar_`, `n_interface_clashes`, `clash_overlap_sum`, `worst_clash` | **red_total_per_ifaceres** |
| sc_DockQ components | `dockq_fnat`, `dockq_iRMSD`, `dockq_std_self`, `dockq_std_LRMSD` | no — see below |
| ipTM variants | `iptm_min3`, `iptm_best2`, `mean_iptm3`, `b2_iptm`, `af3_msaoff_ptm` | no |
| Size | `n_rna_res`, `fnat_factor`, `length` | no |

**Not at campaign scale, and the figure says so:**

- **apo RMSD is bigger than the 100k table suggests.** `apo_rmsd_vs_esm` is on the
  100,000-row `~/Documents/metric_explorer/data/metrics_100k.csv` (247 columns) and
  **not** on the master — but the 2026-07-29 re-launch finished, and the chunk tables are
  still on Euler: `rna_binder/apo_chunks/` (100 chunks, the pilot) and
  `rna_binder/apo_chunks_v2v7/` (3,835 chunks). Together they cover the **v2–v7 corpus**,
  so apo is drawn on **3,934,826 designs with 0 missing values**, 56.3% of the campaign,
  not on 100,000 — the 100k table's `apo_rmsd_vs_esm` is a 100,000-row slice of the same
  run. Coverage per route is 52.6% to 69.1%, so it is **mildly uneven rather than skewed**,
  and in the opposite direction to the obvious guess: HalluDesign is the best covered
  (706,862/1,023,007 = 69.1%) and BoltzGen+BoltzIF the least (778,141/1,480,140 = 52.6%).
  `report()` prints the per-route share and the panel carries the n.
- **100k only**, so not drawn: `sse_*`, `rna_contact_nt_*`, `aa_base_*` / `aa_bb_*`,
  `tip_*`, `interface_fa_rep`, `interface_bsa_approx`, `helix_fraction` /
  `sheet_fraction` / `loop_fraction`, `n_interface_contacts`.
- **Nowhere at all**, unchanged from the previous cut and re-checked 2026-11-19 by
  grepping both headers (0 matches): **ipSAE** (`ipsae*` / `pdockq*` / `lis*`) — wired
  into the re-prediction pipeline 2026-07-02 but forward-only, needing the full L×L PAE the
  pipeline discarded, so retroactive scoring means re-predicting the corpus; and a
  **Rosetta interface energy** (`dG_separated`, `interface_energy`, total score). The
  nearest stored thing is `interface_fa_rep`, the repulsive term alone, on the 100k subset
  only. **NEED:** re-predict a defined shortlist if either is wanted.

## Which five campaign metrics, and why — measured, not asserted

`report()` re-derives both columns on every run for the five drawn metrics: **route
separation** as (max − min of the five route medians) ÷ the campaign IQR, and **Pearson r
against `dsasa_total`**, which the main text already plots. A member with a high |r| is
the interface-size panel restated under another name. The **cut** rows come from the
one-off candidate sweep over all eighteen members (2026-11-19, same consensus rules, all
6,992,515 designs) and are not re-derived — the drawn rows are.

| candidate | spread ÷ IQR | r vs ΔSASA total | verdict |
|---|---|---|---|
| `delta_unsat_hbonds` | 1.00 | **+0.628** | **drawn** — BUNS, named in the comment. The one drawn metric that FAILS the independence test, and it is stated rather than dodged: see the box below |
| `packstat` | **0.76** | **−0.328** | **drawn** — the clearest diffusion-vs-hallucination split of the five, 0.62–0.64 against 0.56. Its 45,313 exact zeros are dropped (see below), which is why its r is −0.328 here and −0.136 in the candidate sweep that kept them |
| `frac_dsasa_base` | **0.66** | **−0.305** | **drawn** — the ΔSASA member that is not interface size |
| `frac_dsasa_major` | **0.63** | +0.382 | **drawn** — which base edge the area comes off |
| `red_total_per_ifaceres` | 0.42 | **+0.022** | **drawn** — the one candidate fully independent of interface size, and monotone across routes |
| `buns_vsasa` | **1.15** | +0.368 | **cut** — the better-looking numbers, the wrong quantity. See the box below |
| `dsasa_polar` | 0.98 | +0.975 | cut — ΔSASA total again |
| `dsasa_hphobic` | 0.93 | +0.927 | cut — same |
| `dsasa_rna_phosphate` | 0.98 | +0.838 | cut — same |
| `dsasa_rna_base_major` | 0.61 | +0.783 | cut — same |
| `dsasa_rna_base` | 0.68 | +0.746 | cut — same |
| `dsasa_rna_sugar` | 0.69 | +0.694 | cut — same |
| `dsasa_rna_base_minor` | 0.55 | +0.238 | cut — independent, but a 25–45 Å² area with a weak route signal; the base-share fraction carries it normalised |
| `dsasa_rna_o2prime` | 0.37 | +0.388 | cut — same, 20–30 Å² |
| `n_h_bonds_interface` | 0.50 | +0.507 | cut — the medians are 1–2 H-bonds. ⚠️ Until 2026-11-20 this row also said "already in `supp/interface_chemistry`"; that figure is archived, so the metric is now drawn nowhere |
| `sc_value` | **0.20** | +0.078 | cut — 0.364 to 0.381 across all five routes. Nothing to see at campaign scale |
| `pi_stack_score` | — | +0.021 | cut — the worst-case consensus is **0 at the median of every route with zero IQR**. A violin of it is a line. ⚠️ `supp/interface_chemistry` drew it per route until it was archived 2026-11-20; it is now drawn nowhere |

### ⚠️ Which BUNS — the one place this figure draws the metric with the WORSE numbers

The master carries two BUNS columns and they are **not two estimators of one quantity**:

- `delta_unsat_hbonds` is `InterfaceAnalyzerMover.get_all_data().delta_unsat_hbonds` —
  a **bound-minus-unbound interface delta**, which is what BUNS means.
- `buns_vsasa` is a **whole-pose absolute count** (VSASA burial, probe 1.1, all heavy-atom
  unsatisfieds). It is absolute and not a delta because
  `BuriedUnsatHbondFilter.set_use_ddG_style(True)` **returns 0.0 for every protein–RNA
  pose** in this PyRosetta build — six configurations were tried (`rna-binder-buns-metric`).

`buns_vsasa` has the prettier numbers on both of this figure's criteria (1.15 × IQR
against 1.00, r +0.368 against +0.628). **It is still cut**, because its apparent
independence from interface size is a consequence of not being an interface measurement
at all, so taking it would have bought a better panel by quietly swapping the quantity.
`report()` prints its spread, r and route medians on every run so the comparison stays
visible. `delta_unsat_hbonds` is drawn, its +0.628 is stated in the caption, and the
reader is told that most of the route difference in interface BUNS is interface size.

⚠️ **The main text's BUNS row uses `buns_vsasa`.** `Fig5def_asn_interface` scores
`buns_vsasa_sel` and `Fig5fg_h1null_refolder` the same, both labelled "Buried unsat.".
Not changed here — but if those panels return to the page, the label needs the
whole-pose caveat or the column needs swapping. `Fig5pr_precision_recall` has it right:
`delta_unsat_hbonds` is its primary and `buns_vsasa` a `primary = 0` variant.

**That is the answer to "ΔSASA is a family, not one column".** Every absolute member
correlates +0.69 to +0.98 with the total the main text already shows; the two
**fractions** do not (−0.305, +0.382). So the two fractions are drawn and the ten
absolute members are not.

`dockq_fnat` and `dockq_iRMSD` are unplotted anywhere and are **deliberately not drawn**:
they are measured against the **designed backbone**, so they are wrong for BoltzGen+BoltzIF
in exactly the way the stored sc_DockQ and ligand_RMSD were — and the true-backbone rescore never
re-emitted them (next section). Drawing them would put a known-wrong 1.48 M-design series
on the page with no correction in existence.

## THE STRIP — BoltzGen+BoltzIF stored vs corrected, restored 2026-11-19

> "why did you remove the pre-fix bg+bif. I thought it was informative"

He is right, and the rebuild had dropped it by accident rather than on purpose: the pair
rode on sc_DockQ and ligand_RMSD, the audit cut those two as already-in-the-main-text, and the
comparison went with them. It is back, as a **two-panel strip under the grid**, because
it is the only place a reader sees the size of the defect **against the campaign's other
metrics and the four design routes it never touched**.

### Why only two metrics carry a pair — a file-level fact, not a judgement

The true-backbone rescore
`euler:rna_binder/analyses/bif_true_bb_rescore/full/bif_true_bb_scores.csv.gz` has
**16 columns**: `design_id, uuid8, campaign, D, paired, unpaired_reason`, then
`dockq_LRMSD_*_true`, `dockq_self_*_true`, `true_status_*` for af3 / b2 / esm, and
`true_bb_flag`. It emits **two quantities per re-predictor and nothing else**. So for
every one of the six metrics in the grid above there is **no `*_true` counterpart in
existence** — not merely no difference, nothing to draw. The reason is mechanical: the
correction only touches quantities measured **against the designed backbone**, and all
six grid metrics are read off the re-predicted complex (apo, off two ESMFold2
predictions), which the defect never touched. The defect was in which structure the sc_DockQ
reference came from. **The caption says this**, because it is the answer to the first
question a reader will ask of the strip.

### Why these two are re-admitted when the figure's rule is "not already shown"

sc_DockQ and ligand_RMSD **are** in the main text — Fig 2c, Fig 3a, Fig 3b — so drawing them
breaks this figure's own rule. They are re-admitted deliberately and the **caption flags
them as such**, for one reason the main text cannot carry: **the main text plots the
corrected values only**. Fig 2c substitutes them, Fig 3a does too, and
`code/supp/bif_effect` shows stored against corrected on the 2,397 **ordered-panel**
members — a set selected for scoring well under the wrong reference, which is why that
figure carries its own "not a route-wide failure rate" caution. Nowhere else does a
reader see **stored and corrected at campaign scale with four unaffected design routes
beside them for scale**. The other six panels keep the rule untouched.

Do not let this drift from `code/supp/bif_effect`: that figure owns the **consequence**
(81.9% of the ordered panel flips pass → fail), this strip owns the **campaign-scale
distribution**. Neither restates the other.

### Verified on every run, re-derived not quoted

Stored n = **1,478,680** of 1,480,140, corrected n = **1,436,724**, **1,435,617**
carrying both. Median sc_DockQ falls **0.0581 → 0.0202**, median ligand_RMSD rises
**30.4566 → 40.2242 Å**. On the **common subset** — the designs carrying both values,
which is what answers the objection that the shift is a coverage artefact — **0.0593 →
0.0202** and **30.1879 → 40.2211 Å**. The same size either way, so it is not coverage.

The strip applies **the same consensus rule as the rest of the figure** to the `*_true`
columns (`min` for sc_DockQ, `max` for ligand_RMSD, over AF3-MSA-off and ESMFold2), so the two
sides are compared like with like, and `strip_frame()` asserts the join is one-to-one
over all 1,480,140 designs. It also asserts the master's own **`bif_rescored` flag marks
exactly** the designs the rescore could score — a cross-file check, and the only thing
here that would catch one of the two files being swapped for another vintage of itself.

⚠️ **On BoltzGen+BoltzIF only, `campaign_metrics_full.csv.gz`'s `sc_dockq` / `sc_rmsd`
are the CORRECTED consensus where one exists**, while its per-re-predictor `dockq_af3` /
`dockq_esm` / `rmsd_af3` / `rmsd_esm` are the STORED values. Measured 2026-11-19: on the
43,416 rows with `bif_rescored = 0`, `sc_dockq` equals `min(dockq_af3, dockq_esm)` with
0 mismatches; on the 1,436,724 with `bif_rescored = 1` it does not, and its median
matches the true-backbone file's to 4 dp. **This is why the consensus assertion runs on
RFdiffusion3+LigandMPNN** and must keep doing so — running it on BoltzGen+BoltzIF would
fail for a reason that is not a bug. The strip takes stored from the per-re-predictor
columns and corrected from the `*_true` file by name, so it never depends on this.

## Consensus — one rule for the whole figure

⚠️ **Since 2026-11-20 the consensus is REPORTED, not DRAWN.** The split violins draw its
two inputs, one per half. Everything below still holds: it is the rule the assertion
checks, the rule `report()` tabulates, and the rule every verified median in this README
is taken under. Read it as the figure's definition of a worse-case number, not as a
description of a body on the canvas.

The **worst of the two gate re-predictors**, AF3 (MSA off) and ESMFold2: `min` where high is
good (`packstat`), `max` where low is good (`delta_unsat_hbonds`, `red_total_per_ifaceres`), and
**`mean` for the two composition fractions**, which have no better or worse end — the
same exception Rg took in the previous cut. Boltz-2 is excluded, as it was from the gate.
apo RMSD takes **no consensus**: it is a single ESMFold2 measurement, the binder folded
alone against the ESM complex refold, so there is no second re-predictor to be worse than.

**The rule is asserted, not described.** `assert_consensus_rule()` reads
`campaign_metrics_full.csv.gz` and checks that `min(dockq_af3, dockq_esm)` reproduces the
stored `sc_dockq` and `max(rmsd_af3, rmsd_esm)` the stored `sc_rmsd`, to 1e-4, on all
1,524,796 RFdiffusion3+LigandMPNN designs — a route the BoltzIF defect never touched. It
runs **before** anything is drawn. This is the only place in the repo that checks the
worst-of-two convention against a column computed by someone else, which is why it stays
even though neither sc_DockQ nor ligand_RMSD is now on the page.

## Sources

| file | mtime | what is taken |
|---|---|---|
| `inputs/large/metric_panel_v2.csv.gz` | **new, 2026-11-19** | `stack` + `packstat`, `delta_unsat_hbonds`, `frac_dsasa_base`, `frac_dsasa_major`, `red_total_per_ifaceres`, plus `dsasa_total` and `buns_vsasa` (reported, not drawn), each `_af3` and `_esm`, all 6,992,515 rows |
| `inputs/large/apo_rmsd_all.csv.gz` | **new, 2026-11-19** | `stack, apo_rmsd_vs_esm` from both apo runs |
| `inputs/large/campaign_metrics_full.csv.gz` | 2026-11-05 11:48 | read **once**, for `assert_consensus_rule()` and for the strip's **stored** sc_DockQ / ligand_RMSD: `design_id, stack, sc_dockq, sc_rmsd, dockq_{af3,esm}, rmsd_{af3,esm}, bif_rescored` |
| `inputs/large/bif_true_per_backend.csv.gz` | 2026-11-17 15:18 | the strip's **corrected** values: `dockq_{af3,esm}_true`, `rmsd_{af3,esm}_true`, 1,480,140 rows. Shared with `code/supp/bif_effect`, which reads the same file |

The two 2026-11-19 files were extracted on Euler in one `mawk` pass each, **columns
resolved by header name, never by index**, and carry `stack` so **nothing is joined** for
the grid — three frames, each keyed by nothing, each counted independently. The **strip
is the one place this figure joins**, on `design_id`, one-to-one, asserted, and only
because stored and corrected live in two files; the common-subset medians it prints need
the designs matched. `metric_panel_v2.csv.gz` comes from
the 218-column master; `apo_rmsd_all.csv.gz` is `cat` over 3,935 `chunk_*_apo.csv` files
in `rna_binder/apo_chunks/` and `rna_binder/apo_chunks_v2v7/`. Added to
`inputs/MANIFEST.md`. The previous cut's `inputs/large/metric_panel_dsasa_foot.csv.gz` is
**no longer read by this figure** — left in place, since it is the only local copy of
`foot_tp_min_dist`, and its MANIFEST row keeps the warning that `foot_*` is a distance and
not the gate's `footnb_010` atom count.

## Metric definitions, as drawn

The `min` / `max` / `mean` in each note is the **consensus rule** — reported, and the rule
the verified medians are taken under. Since the split, each panel **draws the two columns
themselves**, AF3 as the left half and ESMFold2 as the right. Apo RMSD has one column and
so one whole body.

| panel | column | note |
|---|---|---|
| Packstat | `packstat_{af3,esm}` | min. Rosetta packing-quality score, 0–1. **Exact zeros dropped**, see the coverage note |
| Interface BUNS ↓ | `delta_unsat_hbonds_{af3,esm}` | max. Rosetta InterfaceAnalyzer `delta_unsat_hbonds`: buried unsatisfied polar **atoms** at the interface, bound minus unbound, `set_pack_separated(False)` so no repack. It counts ATOMS; `n_h_bonds_interface` counts BONDS, and the two must never be ratioed without saying so. The `_opth` variants (after polar-hydrogen optimisation) and the whole-pose `buns_vsasa` are on the master and not drawn |
| Unfav. contacts ↓ | `red_total_per_ifaceres_{af3,esm}` | max. Unfavourable contacts per interface residue: clashes + acid–phosphate + hydrophobic–polar, summed |
| ΔSASA base share | `frac_dsasa_base_{af3,esm}` | mean. RNA-side buried area falling on bases rather than on phosphate, sugar or 2′-OH |
| ΔSASA major share | `frac_dsasa_major_{af3,esm}` | mean. Of that base area, the share off the major-groove edge |
| Apo RMSD (Å) ↓ | `apo_rmsd_vs_esm` | no consensus — one ESMFold2 measurement. Binder folded **alone**, Cα-RMSD against the ESM complex prediction: does the binder hold its fold without the RNA. Log axis; route medians run 1.01 Å (BG+BIF) to 9.94 Å (HD) |
| **strip** sc_DockQ | `dockq_{af3,esm}` (stored) · `dockq_{af3,esm}_true` (corrected) | min, the same rule on both sides. **Log axis** — see the drawing choices |
| **strip** ligand_RMSD (Å) ↓ | `rmsd_{af3,esm}` (stored) · `rmsd_{af3,esm}_true` (corrected) | max, the same rule on both sides. Linear |

**Two things about `red_total_per_ifaceres` the panel does not say.** The acid–phosphate
term inside it was tested as a predictor of binding and **withdrawn** (it did not survive
the best-re-prediction re-run); nothing here claims it predicts anything. And it is a **rate per
interface residue**, which is why it is nearly independent of interface size — the
numerator and the denominator both scale.

## Drawing choices

- **Split violins, in Fig 2c's idiom** (2026-11-20): AF3 (MSA off) is the left half at
  full strength, ESMFold2 the right at alpha 0.5, both in the design route's colour, one
  key in the row gutter. The **half** carries the re-predictor; the alpha is the
  redundant cue, as it is in Fig 2c, which is why its absence on the strip's one open
  body costs nothing. Each half keeps the whole body's 0.40 reach, so the split changed
  what is inside the envelope and not the envelope: row pitch, limits and both label
  blocks are untouched by it. **Apo RMSD stays one whole body** — one ESMFold2
  measurement, no second half in existence — the way Fig 2c's "Designs" panel stays a
  bar. The convention lives in `BACKENDS` and nowhere else.
- The density is a histogram (160 bins) smoothed with an 11-bin Hann window, computed in
  the panel's own axis space (log where the axis is log), over **every value** — 1.4 M
  per route. That, not scipy, is why Fig 2c's `_half()` is not imported: it KDEs a
  60,000-design subsample per route, which is what a KDE costs at this size.
  ⚠️ **Corrected 2026-11-20:** this bullet used to say matplotlib's `violinplot` was
  ruled out because it calls `scipy.stats.gaussian_kde`. It does not — it uses its own
  `mlab.GaussianKDE` and runs fine here. scipy IS unusable in `~/.venvs/binder_map`
  (`_propack` dlopen failure), which is why every statistic in `report()` — median, IQR,
  Pearson r, the spread ratios — is hand-rolled; it was never why the density is.
- Each body is drawn only where its density is ≥ 0.5% of its own peak. Without that trim
  every violin ends in a hairline needle at the axis limit — the clipped tail drawn at
  zero width — which reads as a stray rule rather than as data.
- **Two blocks of x labels, one per block of panels.** Inside the grid only its bottom
  row is labelled — all six grid panels have the same five positions in the same order,
  so repeating five names six times is ink, not information (STYLE_NOTES §0). The strip
  carries its own, because its positions are **six, not five**, and two of them are not
  route names.
- **No hatch anywhere** (palette v0.6, 2026-11-19). The restored stored-vs-corrected pair
  is carried by **fill versus outline plus adjacency** — corrected filled, stored open,
  side by side in BoltzGen+BoltzIF's own slot in the fixed route order — which is what
  the previous cut used too, so the hatch withdrawal was a no-op here:
  `grep -n "hatch=" code/supp/metric_panel/make_figure.py` returns nothing, before the
  change and after. **Both pair members are the same colour**, the route's own: a route
  colour means route identity and nothing else (COLORS.md rule 2), and the pair is one
  route measured two ways. Nothing in this figure distinguishes a re-predictor, so no
  mark takes a `palette.REFOLDERS` colour (rule 7): every body is a design route in its
  route colour, in `palette.STACKS_SHORT` order. **No hue is hand-picked anywhere** — all
  six colours come from `palette.stack_color()`.
- The rule inside each body is the median, drawn in `palette.SURFACE` on a filled body,
  matching Fig 2c's white median marker. **On the strip's one open body it is the route
  colour instead** — white on white is invisible.
- Apo RMSD is on a log axis because its route medians span 1.01 Å to 9.94 Å. **The
  strip's sc_DockQ is on one for the same reason and a sharper one**: its route medians run
  0.021 to 0.081 under a limit of 0.75, so on a linear axis the comparison the strip
  exists for — stored 0.058 against corrected 0.020 — is a 5% sliver against the frame,
  where on the log axis it is about a fifth of the panel. A violin is not a bar, so
  STYLE_NOTES §8 does not apply: the metric axis carries position and the density is
  encoded by width, not by length from zero. ligand_RMSD needs no log — its medians sit
  mid-axis and the shift is 30.5 → 40.2 Å.
- **Axis limits are measured, not chosen.** Each panel runs from about the 0.5th to the
  99.5th percentile. A limit set inside that range leaves the clipped mass piled against
  the frame as a flat-topped body, which `TRIM` cannot remove — it only removes thin
  tails. The first cut of this panel did exactly that on unfavourable contacts (cut at 4,
  where the 99.5th percentile is 4.38 and the 99.9th 6.31) and on both ΔSASA shares.
  ⚠️ **The percentile is taken PER SERIES, not pooled**, and the first render of the
  strip got this wrong. Its six series are far more unlike each other than the grid's
  five routes, so a pooled limit is set by the crowd: pooled p99.5 is 0.665 for sc_DockQ
  and 66.1 Å for ligand_RMSD, and at those limits **1.02% of RFdiffusion3's sc_DockQ and
  1.83% of HalluDesign's ligand_RMSD** sat above the frame as flat tops — visible in the
  render, which is how it was caught.
  **Since the split the series are the HALVES**, and `check_limits()` covers the grid as
  well as the strip: 26 series instead of 6, every one of them something actually drawn,
  printed as a table on every run and asserted in one place. The halves reach further
  than the consensus they replaced in whichever direction the min/max rule was not
  looking — a `min` consensus hides how high AF3 and ESMFold2 each go, a `mean` hides
  both ends — which is why three limits moved on 2026-11-20: **ΔSASA base share
  0.1–0.75 → 0.05–0.77** (RFdiffusion3's AF3 half needs 0.0552, BoltzGen+BoltzIF's
  ESMFold2 half 0.758), **ΔSASA major share 0.1–1 → 0–1** (see the zero mass below) and
  the strip's **sc_DockQ 0.75 → 0.85** (RFdiffusion3's ESMFold2 half needs 0.8308 where
  the consensus series needed 0.7358). Packstat, interface BUNS, unfavourable contacts,
  apo RMSD and ligand_RMSD were already wide enough for their halves.
- A series that fails the check is drawn as a body running off the frame at full width,
  not as a flat top: `density()` has not clamped since 2026-11-19. Either way the limit
  is wrong, and the assertion names the series and its percentile range rather than
  leaving the next person to find it in the render.
- Direction arrows only where the good end is not self-evident (STYLE_NOTES §2): `↓` on
  buried unsatisfied, unfavourable contacts and apo RMSD. Packstat is a score and the two
  ΔSASA shares are compositions with no good end, so neither takes one.

## Verified — re-derived and printed on every run

Designs per route: RFdiffusion3+LigandMPNN 1,524,796 · BoltzGen+LigandMPNN 1,519,958 ·
BoltzGen+BoltzIF 1,480,140 · BoltzDesign1+LigandMPNN 1,444,614 ·
HalluDesign+LigandMPNN 1,023,007 — 6,992,515 in all, the whole campaign, no sampling.

Medians under the figure's own rule, route order fixed (2026-11-19). ⚠️ **These are the
consensus, which since 2026-11-20 is reported and not drawn** — the panels draw the two
halves, whose medians `report()` prints beside each row. The split is worth most where
those halves disagree, and interface BUNS is the clearest: RFdiffusion3 is **20 (AF3)
against 29 (ESMFold2)** under a `max` consensus of 30, BoltzDesign1 **23 against 32**
under 34. ΔSASA major share runs the other way, AF3 above ESMFold2 in every route
(RFdiffusion3 0.8046 against 0.7508), and unfavourable contacts barely splits at all
(RFdiffusion3 0.7826 against 0.7857), which is the honest result for that panel.

| metric | RFD3+LiMPNN | BG+LiMPNN | BG+BIF | BD1+LiMPNN | HD+LiMPNN |
|---|---|---|---|---|---|
| Packstat | 0.6277 | 0.6357 | 0.6228 | 0.5676 | 0.5628 |
| Interface BUNS (atoms) | 30 | 20 | 18 | 34 | 28 |
| Unfav. contacts / iface res | 1.067 | 1.227 | 1.200 | 1.385 | 1.412 |
| ΔSASA base share | 0.2992 | 0.3662 | 0.3668 | 0.3601 | 0.3371 |
| ΔSASA major share | 0.7721 | 0.6817 | 0.6569 | 0.6869 | 0.7652 |
| Apo RMSD (Å), n = 3,934,826 | 3.91 | 1.16 | 1.01 | 8.30 | 9.94 |

The strip, same rule, six positions (BoltzGen+BoltzIF twice). `n with a value` is the
column that matters here — stored and corrected cover different numbers of designs:

| metric | RFD3+LiMPNN | BG+LiMPNN | **BG+BIF stored** | **BG+BIF corrected** | BD1+LiMPNN | HD+LiMPNN |
|---|---|---|---|---|---|---|
| sc_DockQ | 0.0613 | 0.0210 | **0.0581** | **0.0202** | 0.0418 | 0.0814 |
| ligand_RMSD (Å) | 30.19 | 39.84 | **30.46** | **40.22** | 29.42 | 26.20 |
| n with a value | 1,524,766 | 1,476,308 | **1,478,680** | **1,436,724** | 1,444,524 | 1,016,258 |

Of 1,480,140 BoltzGen+BoltzIF designs, **1,435,617 carry both** a stored and a corrected
value. On that common subset the shift is the same size as on the full series —
**0.0593 → 0.0202** and **30.1879 → 40.2211 Å** against **0.0581 → 0.0202** and
**30.4566 → 40.2242 Å** — so it is the correction, not a coverage artefact. Both pairs
are printed on every run.

Read the corrected body against **BG+LiMPNN**, the route that shares its backbones: the
corrected medians (0.0202, 40.22 Å) land essentially on BoltzGen+LigandMPNN's (0.0210,
39.84 Å), where the stored ones (0.0581, 30.46 Å) did not. That is the strip's second
reading, and it is the one the grid cannot give.

⚠️ **ΔSASA major share has a zero mass of its own, and it is drawn, not dropped.**
Also AF3's: **39,277 (0.567%) of `frac_dsasa_major_af3` are exactly 0 against 575 of
`frac_dsasa_major_esm`**, running 0.18% (HalluDesign) to 0.96% (RFdiffusion3) per route.
The `mean` consensus hid it — a 0 averaged against a 0.4 reports 0.2 — so it appeared
only when the halves went on the canvas. It is kept because, unlike packstat's zeros, it
is a **real composition value**: base area with no major-groove-edge share, not Rosetta
reporting that there is nothing to measure. There is no empty gap above it either, so it
draws as a thickened lower tail rather than as a hairline, and the panel's axis now
honestly starts at 0.

Coverage, campaign-scale panels: `delta_unsat_hbonds` 6,992,515/6,992,515 (100.0%),
`red_total_per_ifaceres` 6,992,508/6,992,515 (100.0%, 7 missing), `packstat`
6,947,202/6,992,515 (99.4%, see the zeros below), `frac_dsasa_base`
6,927,921/6,992,515 (99.1%) and `frac_dsasa_major` 6,921,342/6,992,515 (99.0%) — the two
fractions are undefined where the RNA-side base area is zero.

⚠️ **The zeros are AF3's, which only the split could show.** Of the 45,313, **45,310
are in `packstat_af3` and 5 in `packstat_esm`** — so Rosetta finds no packable interface
in 1.86% of BoltzDesign1's and 1.49% of HalluDesign's **AF3** predictions and in
essentially none of anybody's ESMFold2 ones. The `min` consensus inherited every one of
them, which is why the count below is unchanged. `derive()` now drops zeros from **the
consensus column and both halves**, after the consensus is computed, so no reported
number moves and the panels cannot quietly put the point mass back; `check_limits()`
caught exactly that, as a 0.5th percentile of 0 on two AF3 halves.

⚠️ **45,313/6,992,515 (0.65%) of designs have `packstat` exactly 0 and are dropped from
that panel alone** — BoltzDesign1+LigandMPNN 26,808/1,444,614 (1.86%) and
HalluDesign+LigandMPNN 15,244/1,023,007 (1.49%) carry almost all of them. A packstat of
exactly 0 is Rosetta saying there is no packable interface, not that the interface packs
badly, and between 0 and 0.35 there are only 1,310 designs campaign-wide, so it is a point
mass with an empty gap above it. A violin cannot draw a point mass: clipped into the
bottom bin it is a bump against the frame, drawn from zero it is a hairline across the
gap. `report()`'s "n with a value" states the drop on every run. apo RMSD 3,934,826/3,934,826 (100.0%) within its corpus, which is
3,934,826/6,992,515 (56.3%) of the campaign: RFdiffusion3+LigandMPNN
815,447/1,524,796 (53.5%), BoltzGen+LigandMPNN 817,960/1,519,958 (53.8%),
BoltzGen+BoltzIF 778,141/1,480,140 (52.6%), BoltzDesign1+LigandMPNN
816,416/1,444,614 (56.5%), HalluDesign+LigandMPNN 706,862/1,023,007 (69.1%).

Not drawn but printed every run, so the BUNS choice stays checkable: `buns_vsasa`
spread 1.15 × IQR, r vs ΔSASA total +0.368, route medians 37 / 37 / 37 / 52 / 46.

`report()` also prints the apo per-route coverage, which is **the number to read before
comparing apo across routes**: the apo runs covered v2–v7 only.

## Caption draft

**The metrics the campaign measured and the main text does not show.** Split violins over
all 6,992,515 designs, one per design route in fixed order: each body is split down its
centre line into the two gate re-predictors, **AF3 with MSA off on the left, ESMFold2 on
the right**, with the rule inside each half its own median and the density a smoothed
histogram trimmed at 0.5% of its own peak. **Apo RMSD is one whole body**, being a single
ESMFold2 measurement with no second re-predictor to split. Where this figure quotes one
number per design route — the medians in the text, and the independence tests behind the
choice of panels — it is under the campaign's own worse-of-the-two rule (`min` for
packstat, `max` for interface BUNS and unfavourable contacts, the `mean` for the two
ΔSASA shares, which have no better end). ΔSASA is a family, and the two **shares** are drawn rather than
its ten absolute members because each of those correlates +0.69 to +0.98 with the buried
ΔSASA total the main text already plots, while the shares do not (−0.305, +0.382).
Interface BUNS is the InterfaceAnalyzer bound-minus-unbound delta, and most of its
difference between routes is interface size (r = +0.628). Apo RMSD — the binder folded
alone by ESMFold2 against its complex prediction — covers 3,934,826 of the 6,992,515 designs,
stated on the panel; packstat drops the 45,313 at exactly zero, where there is no packable
interface. ipSAE and a Rosetta interface energy are on no campaign-scale file.

**The two lower panels, sc_DockQ and ligand_RMSD, are the one exception to "not already
shown", and are repeated from the main text (Fig 2c, Fig 3a, Fig 3b) on purpose.**
BoltzGen + BoltzIF's stored scores were measured against a BoltzIF re-prediction rather
than against the generated backbone; the true-backbone rescore corrects them, and it
emitted **sc_DockQ and ligand_RMSD and nothing else**, so these are the only two metrics in the
campaign with a corrected counterpart in existence — which is why no panel above carries
a pair, and why the comparison can be drawn here and nowhere else in this figure. Open,
the 1,478,680 values the campaign stored and gated on; filled, the 1,436,724 the rescore
produced; 1,435,617 designs carry both. Both are split by re-predictor like every other
body, and on the open one the half alone says which, a tone being invisible on an
outline. Median sc_DockQ falls 0.058 → 0.020 and median
ligand_RMSD rises 30.5 → 40.2 Å, the same shift on the common subset (0.0593 → 0.0202,
30.19 → 40.22 Å), landing the route on BoltzGen + LigandMPNN, which shares its backbones.
The main text plots the corrected values only; `code/supp/bif_effect` shows what the
defect cost the ordered panel. The other four design routes are unaffected and carry one
body each.

## Open

- **Fig 5 d–f are off the assembled page but still in `FIGURES.md`.** The audit treats
  them as main-text, which is why this figure does not lead on BUNS-as-a-predictor. If
  they stay off, the twelve-metric list becomes genuinely unshown and a seventh panel
  (shape complementarity, heavy contacts) becomes defensible — though `sc_value`'s route
  spread of 0.20 × IQR argues against it either way. Settle Fig 5 first.
- **Two `hatch=` calls survive in main-text scripts**, found while applying palette v0.6:
  `code/Fig2cd_benchmark/make_figure.py:149` (design-mode class key) and
  `code/Fig3b_filter_venn/make_figure.py:273` (the "removed" foot-clash bar). Neither
  keys a re-predictor, so COLORS.md rule 7 does not reach them, but Leo's instruction was
  "remove the dots and the stripes" from all figures. Not touched here. The three supp
  figures that did key re-predictors by hatch — `interface_chemistry`,
  `metric_reliability`, `repredictor_comparison` — are converted in code **and**
  (corrected 2026-11-19) their READMEs now describe the colour encoding too, so these two
  main-text calls are the only `hatch=` left in the repo.
- **apo coverage is uneven across routes (52.6–69.1%)** and the panel states only the
  total n. It is mild and does not run with design mode, so it probably needs no further
  caveat — but each apo violin is a different *sample* of its route than the five
  campaign-scale violins above it, and nothing on the canvas says so. Decide whether one
  clause in the caption closes it.
- **Overlap with `code/supp/bif_effect` is now real and has to stay deliberate.** Both
  figures draw stored against corrected, from the same two files. The split: **this strip
  is campaign-scale** (1.48 M designs, four unaffected design routes beside it, the
  distribution) and **`bif_effect` is consequence-scale** (the 2,397 ordered-panel
  members, paired per design, 81.9% flipping pass → fail, with its own "not a route-wide
  failure rate" caution). Each states the other's role; the medians agree because they
  are derived from the same columns under rules each figure names. If one is re-cut,
  check the other — a silent disagreement between them is the failure mode here.
- **The stored body is the one split a reader cannot see.** `BG+BIF stored` is an
  outline, so it takes no tone, and its two halves differ too little (sc_DockQ 0.0878
  against 0.0814, ligand_RMSD 26.07 against 25.69) for the shape asymmetry to carry them
  either. The pair it belongs to is unaffected — hollow against filled two-tone is if
  anything a sharper contrast than hollow against filled was — but on that one position
  the split is decorative. Leaving it split, because a position that opted out of the
  figure's own idiom would need explaining more than this does. `NEED:` Leo's call if it
  should instead be drawn whole.
- ~~**The strip's y labels repeat the main text's.**~~ **CLOSED 2026-11-20.** The strip
  now reads `sc_DockQ` and `ligand_RMSD (Å)`, which is exactly what Fig 2c reads, so a reader
  moving between them sees one quantity under one name. If either is ever relabelled
  again, move both together.

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
