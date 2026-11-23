# metric_reliability

*Was `seed_reliability` until 2026-11-19, and `S2` in the 2026-10-02 set. Built under its
slug: supplementary numbers are assigned in one pass once the main text stops moving, so no
S-number appears in any file name.*

How precisely each interface metric can be measured at all. Two bare panels.
Source: `~/rna_binder_plots/seed_metric_reliability/` (2026-08-06) — 2,000 designs per
re-predictor, 5 refolding seeds each.
Run: `~/.venvs/binder_map/bin/python3 code/supp/metric_reliability/make_figure.py`.
Emits `output/metric_reliability/panel_{a,b}_*.{pdf,png,svg}`. **(a) is
178 × 68.5 mm** (double column, flipped 2026-11-20); (b) is 85 × 46 mm.

🚨 **`NEED:` a citing sentence in Section 3.** This figure is currently uncited, and it is
what justifies using H-bond and stacking counts descriptively rather than as gates —
without it that choice looks arbitrary and the figure is an orphan.

## Why the slug is `metric_reliability` (renamed 2026-11-19)

Leo: *"I think the name is bad. For me its more a precision / reliabiltiy test. The seed was
a the way we got there just a tool."* He is right — `seed_reliability` named the
**instrument**. The five seeds are only how the measurement error is exposed; the quantity
on both panels is a property of the **metric**.

Of his two words, `reliability` is the one taken, for two reasons:

1. **It is the name of the statistic the figure plots.** ICC(1,1) *is* a reliability
   coefficient — the axis reads `ICC (5 seeds)` — and the source tables are already
   `seed_metric_reliability.csv` / `pistack_presence_reliability.csv`, with a `refolder`
   column. The slug now matches the data and the axis.
2. **`metric_precision` would collide.** The neighbouring `supp/metric_significance` draws
   a **precision**-recall curve, where precision means positive predictive value. Two
   adjacent supplementary figures called "precision" in incompatible senses is a worse
   name than the one it replaces.

The residual risk is that *reliability* can be read as "is this metric any good" — i.e.
validity. It is not that: nothing here says whether a metric measures the right thing, only
how tightly it can be measured. The caption's first sentence is written to say so.

Not reverted to `FigS2_metric_reliability`: the S-number drop on 2026-11-11 was right, it
was the subject half of that rename that was wrong.

## The re-predictor encoding (hatch withdrawn 2026-11-19)

Re-predictors take their own colours from `palette.REFOLDERS` — AF3 `#1b7f79`, Boltz-2
`#2e7d32`, ESMFold2 `#8c6d4f` — iterated in `palette.REFOLDER_ORDER`, the pipeline's order
and never the result's. See COLORS.md rule 7.

**What changed and why.** This figure previously keyed re-predictor by *texture*: one theme hue
at steps 900 / 600 / 300, with `///` on Boltz-2 and `...` on ESMFold2, on the reasoning
that a re-predictor is a measurement instrument rather than a design route, so it must not
borrow a route colour (rule 2), and `palette.json` then carried no categorical triple for
it. The reasoning about route colours still holds; the texture answer does not. Leo asked
on 2026-11-19 for the dots and stripes to come off every figure, and palette v0.6 settled
the question properly by giving re-predictors their own validated triple (ΔE ≥ 45 from
every route colour, so the rule-2 worry is met by colour rather than by hatch). Every
`hatch=` in this script is gone.

Two consequences worth recording:

- **Bar edges dropped.** The hairline ink edge existed because theme 300 is under 3:1 on
  white. The three re-predictor colours are 4.8 / 5.1 / 4.8 : 1, so the edge is now ink for
  nothing.
- **Greyscale identity moved from texture to position.** The three colours sit at
  relative luminance 0.168 / 0.155 / 0.171 — a greyscale print cannot separate them, and
  the hatch is what used to. Identity therefore rests on **order within the group**, which
  is fixed and identical in both panels (checklist item 4). This is the real cost of
  the withdrawal; it is paid, not denied.

The frameless key is drawn once, at the bottom of panel (a), and serves the figure.

## Panels

- **(a)** ICC per metric, 16 metrics in 5 families, 3 re-predictors, as **vertical bars
  across the full double-column width** (2026-11-20; it was 16 horizontal rows on 85 mm
  until then). Families separated by a vertical hairline and named **above** their block.
  Metric order: by family (`SASA, contacts, BUNS, H-bonds, pi-stack`), then by mean ICC
  descending within the family.
- **(b)** In how many of the five seeds a π-stack is present, over designs that show one in
  at least one seed. Five grouped bars (1…5 seeds) per re-predictor; the all-5 bars are
  labelled because 1 design in 2,000 is a bar you cannot see.

### (a) was flipped to the full page width, 2026-11-20

Leo: *"Flip a to be horizontal full width of page. Give the sets of three a bit more
space."* Both done, and the second is the reason the first is worth it.

**The space.** The panel is now 178 mm (`style.CFG["width_double_mm"]`) over 16 metrics,
so a metric gets a **10.25 mm** slot. Three bars at `BAR_W = 0.22` data units span 0.66 of
it — **2.25 mm each** — and the 0.34 left over is a **3.48 mm gap** between neighbouring
sets of three. On the old 85 mm rows the same three bars sat at `h = 0.26`, leaving 0.22 of
the row pitch, and the sets ran together. The bars themselves are unchanged in meaning;
only the axis they stand on turned.

**Where the family separators and names went.** They could not stay where they were. The
hairlines were *horizontal* rules between row blocks with the family named at the **right
edge**, and in the flipped panel the right edge is one metric (`π-stack count`), not one
family — there is no "end of the row" to name any more. So the rules are now **vertical**,
between the column blocks, and each family is named **above** its own block, horizontally.
The top strip is the only clear one left once the bottom carries rotated metric labels, and
horizontal is the only orientation that does not add a second reading angle to a panel that
already has one.

**The labels, measured not estimated.** 16 metric names along an x axis cannot be set
horizontally: the longest, `base-contact fraction`, is **17.75 mm** against a 10.25 mm
slot. They are rotated **45°**, right-anchored, and three things are asserted on every run
from `style.text_w`, never from a per-character guess:

1. **Nothing drops off the canvas.** The longest label's diagonal footprint is
   17.75/√2 = **12.55 mm**; with one line of type it drops **14.8 mm** into a **19.5 mm**
   bottom margin, which also has to hold the key.
2. **No two labels collide.** At 45° adjacent labels are parallel lines offset
   perpendicular by `pitch/√2` = **7.25 mm**, against **2.23 mm** of type. It is the
   perpendicular gap, not the label width, that has to clear — getting that test wrong is
   how 16 rotated labels hide a collision.
3. **No family name is wider than the block it names.** The tight ones are the two-metric
   blocks: `π-stacking` is 8.93 mm against 20.5 mm, `BUNS` 5.15 mm against the same.

**The key had to be rebuilt for the width.** `refolder_key()` laid its three entries out on
a fixed 0.30-of-figure pitch, tuned for 85 mm; at 178 mm that is 53 mm per entry and the
three labels would have scattered over half the panel. It now **measures** each entry with
`style.text_w` and centres the block, in mm off `fig.get_size_inches()`, the same way
`backbone_landscape`'s key does. Panel (b) does not draw a key and is untouched.

**Panel (b) is unchanged** — still 85 × 46 mm, still the same five grouped bars. The
figure is now mixed-width, (a) double and (b) single; assembly places them, and nothing in
`make_supp_book.py` needed changing, since it already scales anything wider than the
178 mm text block and (a) is exactly that.

### The CV-vs-ipTM panel was cut, 2026-11-19

Leo: *"remove b, a&c are enough."* The old panel (b) drew the Spearman ρ between a design's
across-seed CV and its ipTM for 6 representative metrics. It is gone — the drawing code
with it, not just the call — and what is left of the figure renumbers, so the π-stack panel
is now **(b)** and its file is `panel_b_pistack_reseed`.

The three reasons, recorded because they are the reasons and not the instruction:

1. **It argued a different thing.** This figure's claim is about *metrics* — the chemically
   discriminating ones cannot be measured precisely. That panel was about *designs*:
   unstable ones score lower. A side observation that changes nothing about how the metrics
   are used.
2. **It was close to a sanity check.** ipTM is itself a per-prediction confidence output,
   so predictions the model is unsure about being also the ones that move between seeds is
   near-tautological. Its absence would have been alarming; its presence is not a finding.
3. **One caption clause replaces it, and carries more.** All 45 values rather than the 6
   drawn: ρ between −0.540 and −0.090 across 15 metrics × 3 re-predictors.

`metric_std_vs_quality.csv` is therefore still read — by `report()`, not by a panel — so
that clause's range is re-derived and printed on every run rather than copied from here.

### (b): the 0-seed column is not drawn

Designs with no π-stack in any seed are 1,535 / 1,716 / 1,579 of 2,000 — five times the
tallest bar on the panel, and plotting them flattens the shape the panel exists to show.
The column is exactly `2,000 − ever` and the caption gives all three counts, so nothing is
hidden. Bars stay linear from zero (STYLE_NOTES §8); no log axis was considered.

## Numbers the caption carries (all re-derived by the script, printed on every run)

- Family ICC ranges: burial/SASA 0.634–0.928, contacts 0.459–0.888, **BUNS 0.676–0.788**,
  H-bonds 0.163–0.514, π-stacking 0.151–0.519.
- **ESMFold2 has the highest ICC on every one of the 16 metrics**, not merely on every
  family mean.
- Cut panel's clause: all 45 `rho_cv_iptm` values negative, −0.540 to −0.090
  (15 metrics × 3 re-predictors).
- π-stack reseed, per re-predictor, out of n = 2,000:

| seeds with a π-stack | AF3 | Boltz-2 | ESMFold2 |
|---|---|---|---|
| 0 (not drawn) | 1,535/2,000 (76.8%) | 1,716/2,000 (85.8%) | 1,579/2,000 (79.0%) |
| 1 | 301/2,000 (15.0%) | 218/2,000 (10.9%) | 258/2,000 (12.9%) |
| 2 | 101/2,000 (5.1%) | 41/2,000 (2.1%) | 76/2,000 (3.8%) |
| 3 | 30/2,000 (1.5%) | 17/2,000 (0.85%) | 25/2,000 (1.2%) |
| 4 | 15/2,000 (0.75%) | 7/2,000 (0.35%) | 22/2,000 (1.1%) |
| 5 | 18/2,000 (0.90%) | 1/2,000 (0.05%) | 40/2,000 (2.0%) |
| ever (1–5) | 465/2,000 (23.2%) | 284/2,000 (14.2%) | 421/2,000 (21.1%) |

- Of designs with a π-stack in any seed, 90.5% (ESMFold2) to 99.6% (Boltz-2) lose it in at
  least one other seed — and **exactly one seed of five is the commonest outcome**:
  64.7% (AF3), 76.8% (Boltz-2), 61.3% (ESMFold2) of those designs.
- The decay is **not** monotone: AF3 and ESMFold2 both tick back up at all-5 (15 → 18 and
  22 → 40), so their all-5 column is not their smallest. Only Boltz-2 falls the whole way
  (7 → 1). A small genuinely always-stacking population exists and Boltz-2 does not find it.

Three corrections against the 2026-10-02 `CAPTIONS.md`, files win:

1. **"0.1–2.0% stack reproducibly" is wrong at the lower bound.** Boltz-2 is 1 design in
   2,000 = 0.05%. Given as counts here, so the floor cannot be misread.
2. **"90–99.6%" → 90.5–99.6%.**
3. **BUNS is a fifth family that the old caption omitted while the figure drew it.** It
   belongs in the caption: a burial-like metric that replicates (0.676–0.788) while the
   chemistry metrics do not *strengthens* the argument.

## Estimator and provenance, for Methods

One-way random-effects **ICC(1,1)**, k = 5 seeds, n = 2,000 designs per re-predictor —
established by reproducing the stored values to 10 decimals. The merged table carries
**six** AF3 seeds (1–6, seed 6 a later top-up) against five for Boltz-2 and ESMFold2; the
stored ICCs use **AF3 seeds 2–6**, so k = 5 throughout and AF3 seed 1 is unused.

**Panel (b) is counted from the per-seed rows, not from the summary.**
`pistack_presence_reliability.csv` carries only `ever` and `always`, so the 1–5
distribution comes from `seed_metrics_merged_v2.csv` (32,000 rows = 3 re-predictors × seeds ×
2,000 designs), counting designs whose `pis_count_nodangle > 0`, on AF3 seeds 2–6 as above.
`pistack_seed_counts()` asserts that both stored endpoints come back out of the per-seed
counts — all six reproduce exactly (465/18, 284/1, 421/40). That assert is the panel's
check: it is what proves the column and the seed subset are the ones the stored summary
used. `pis_count_all` does **not** reproduce them (478/20, 301/1, 467/43), confirming the
no-dangle variant.

No scipy anywhere — it cannot import in this venv (`_propack` dlopen) — and none is needed:
the ICCs and the Spearman ρ are read from the stored tables, and everything the script
computes is plain counting and a min/max, hand-rolled over those tables.

`seed_metric_reliability.csv`'s **`cv` column is not plotted and must not be**: it is
pinned at exactly √5 for all six π-stack rows, the algebraic maximum when one of five
seeds is non-zero, so it is degenerate there.

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the n. 197 words.*

**How precisely each interface metric can be measured.** Five refolding seeds put an error
bar on every interface metric; the metrics that discriminate designs chemically carry the
widest one. (**a**) Intraclass correlation across the five seeds, for 16 interface metrics
in five families: burial and buried-unsatisfied counts are measured reproducibly (ΔSASA
0.634–0.928, BUNS 0.676–0.788), hydrogen-bond (0.163–0.514) and π-stacking (0.151–0.519)
metrics are not. ESMFold2 gives the highest ICC on every one of the 16. Designs whose
values scatter most across seeds also score lowest on confidence (Spearman ρ between a
design's across-seed coefficient of variation and its ipTM, −0.540 to −0.090 over 15
metrics × 3 re-predictors). (**b**) In how many of the five seeds a π-stack is present, over
the designs that show one in at least one seed — AF3 465/2,000 (23.2%), Boltz-2 284/2,000
(14.2%), ESMFold2 421/2,000 (21.1%); designs showing none in any seed are not drawn.
Exactly one seed is the commonest outcome; all five is reached by 18/2,000 (0.90%),
1/2,000 (0.05%) and 40/2,000 (2.0%), a floor AF3 and ESMFold2 both tick back up to.
π-stacking is the no-dangle variant, excluding the 5′ A1. ICC(1,1), k = 5 seeds, n = 2,000
designs per re-predictor.

## Vocabulary (settled 2026-11-19)

**"re-predictor", never "refolder"**, and **"design route", never "stack"**, in all prose,
labels and legends — supplement-wide, not just here. This figure moved first. Only
human-facing text changes: the CSV column is still `refolder` and the palette API is still
`palette.REFOLDERS` / `palette.REFOLDER_ORDER`, so nothing downstream of a column name or
an import has to move with it.

**Swept supplement-wide on 2026-11-20**, and the decision is now unanimous: "re-predictor"
for the three tools and **"refolding seeds"** for the five seeds, with no competing form
left in any README or script prose under `code/supp`.

⚠️ **Identifiers are deliberately exempt. Do not "fix" them into breakage:** the
`refolder` and `stack` CSV columns, `palette.REFOLDERS`, `palette.REFOLDER_ORDER`,
`palette.STACK_ORDER`, `palette.SHORT`, the source file names containing
`_by_stack_refolder`, the local `refolder_key` helper, and the directory name
`Fig5fg_h1null_refolder`. The one remaining competing word in `code/supp` is **"backend"**
(≈90 uses, nearly all "either gated backend" for the two gate re-predictors); it was left
alone deliberately in the 2026-11-20 sweep because it was not part of the decision, and it
is flagged here rather than changed silently.

## Open

- **ICC-by-family panel ownership.** `fig_F7_triage.py` (2026-10-05) draws an ICC-by-family
  panel that is this figure's panel (a) at family resolution. Decide which figure owns it;
  do not ship both. Panel (a) here is the metric-resolution version and is the stronger of
  the two.
- The 0.5 / 0.75 "poor / good" guide lines from the old version are **not** drawn
  (STYLE_NOTES §0: the caption can carry a threshold). Restore only if a reader must read
  a band off the panel.
