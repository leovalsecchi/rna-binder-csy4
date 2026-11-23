# arrayed_flow

*Built under its slug: supplementary numbers are assigned in one pass once the main text
stops moving, so no S-number appears in any file name.*

The arrayed re-test of the 32 ordered candidates (LV2.33, transfected 2026-11-04; read on
an LSRFortessa 2026-11-06 as LV2.34.1). Three bare panels over one variant set and one
variant order: **(a) runs horizontally**, (b) and (c) run down, and a variant is at the
same place along (a) as it is down (b) and (c).

Run: `~/.venvs/binder_map/bin/python3 code/supp/arrayed_flow/make_figure.py`
Emits `output/arrayed_flow/panel_{a,b,c}_*.{pdf,png,svg}` (`style.save` writes the SVG
too, and that is what page assembly should consume). Authored at 152 × 38 mm for (a) and
86 × 103.6 mm for (b) and (c) (2.4 mm per row; 5.5 pt type is ~1.94 mm, so that pitch is
the floor); `bbox_inches="tight"` crops them to **152.8 × 37.2**, **86.8 × 104.3** and
**86.8 × 104.4 mm**.

Those are not free numbers. `make_supp_book.py` packs panels left to right, in filename
order, by their **emitted** size into rows no wider than the A4 text block (210 − 2 × 14 =
182 mm, 6 mm gutter), so the sizes are what decide the page: **(a) takes its own row and
(b) and (c) sit side by side under it** (179.5 of 182 mm), the page body running 177.7 of
283 mm at scale 1.0 — checked by calling the book's own `rows_for`, not inferred.
`main()` asserts all of it on the files it has just written, so a panel that grows past
its budget fails the run instead of silently re-packing the page.
⚠️ (b) and (c) still differ in height by 0.1 mm, and (a) is a different object entirely,
so **align the 34-row block of (b)/(c) at assembly, not the panel edges**; their rows are
at identical y positions.

**`scipy` is broken in this venv** (`_propack` dlopen), so both tests are hand-rolled and
exact, not normal approximations: the two-sided binomial by enumerating the pmf, the
Wilcoxon signed-rank by DP over the rank set (32 ranks, 528 states). The script asserts
there are no ties in |d|, which is the condition the exact version needs.

## Sources

| file | mtime | what it gives |
|---|---|---|
| `inputs/lv2_34_1_flow_on_off.csv` | 2026-11-09 18:47 | the 34 kept plate columns × 4 targets, `design_id`, route |
| `inputs/candidates32_index.tsv` | 2026-11-09 18:31 | ordered binder ID (`twist_id`), generator / inverse-folder, `hts_comparison`, `flow_well` |
| `/Volumes/wicky/expdata/flow/lvalsecchi/LV2_34_1/LV2_34_1_Analysis.xlsx` | 2026-11-09 14:27 | the FlowJo export the CSV comes from (140 rows = 35 columns × 4) |
| `~/LV2_34_1_ON_OFF_plots/plot_on_off.py` | 2026-11-09 14:36 | the drop rule and the averaging convention |
| `LV_FlowCytometry/analysis/lv2_33_screen.py` | 2026-11-07 10:39 | plate map, including which column is which control |
| `LV_FlowCytometry/analysis/lv2_33_stats.py` | 2026-11-07 10:51 | the empirical-null idea this figure's band is built on |
| `LV_Transfection/LV2.33_transfection_AS-PERFORMED_2026-11-04.html` | 2026-11-05 10:43 | 34 variants × 4 targets × 1 replicate, 136 wells, 84 ng binder |
| `LV_Cloning/pLV_NNN/pLV0{13,14,15,16,17}` + `pLV16_asbuilt_LV2-26-004c_annotated.gb` | 2026-10-22 16:20 / 2026-11-16 19:59 | the four targets and the unused negative |
| `~/Downloads/design29_742db4d3/an32/m32.pkl` | 2026-10-21 11:43 | 12,000 × 281; `pos` marks the 32 (POS_A 7, POS_B 25) — **not read by the script**, see below |

`m32.pkl` is not a dependency. Its `pos` column (POS_A 7, POS_B 25) and the
`design_id`/route of the 32 are already in `candidates32_index.tsv`, which carries the
`flow_well` the join needs; both were checked to agree before the pkl was dropped. Loading
a 35 MB / 281-column table to recover two columns that are already in the repo would be
the only slow step in the build.

## Panels

- **(a)** Every well, **horizontally** since 2026-11-20 (below). Normalized mScarlet per
  variant × target, as log2 of that variant's
  own mean over its four targets. The absolute level spans 0.00294–0.44215 across the 128
  candidate wells (150×) and is a transfection property, so the raw number would show
  nothing else. Scale clipped at ±0.6 log2; the script **asserts** that the only candidate
  wells outside it are the two dead positives (binder 4 posB 0.0061, binder 6 posA
  0.0029), so a third one cannot appear silently. Clipped cells are outlined — which
  includes all eight control cells, off-scale by construction, with their ON/OFF printed
  **above their column** (`11x`, `97x`). Variants run left to right in ON/OFF order,
  strongest at the right, and the two Csy4-H29A doses sit past the rule where they used
  to sit above it; the four targets are the four rows.
- **(b)** ON/OFF ranked, `avg(posA, posB) / avg(NegA, no-scaffold)`, dot coloured by design
  route. **Built like the main text's `Fig5a_PosAvg_vs_NegAvg`** (`code/Fig5a_wetlab_on_off/`)
  — same 32 designs from the same CSV, same ranking, same scale rule, same tick set, same
  clamp-on-the-edge idiom, same colour by stack. Two things it carries that the main-text
  panel does not, and they are the reason the supplementary version exists: the line
  through each dot = range over the four positive × negative pairings, and the band =
  **±2 SD of log2(no-scaffold/NegA)**, the same ratio measured on two wells where no
  binding difference is expected — the only honest spread a single replicate affords. Its
  *centre* is a plasmid offset (−0.115 log2) and is deliberately **not** used; the band is
  centred on the no-effect line, after `lv2_33_stats.py`, which does this on the per-cell
  mScarlet-positive fraction instead.
- **(c)** The four pairings per variant, on **(b)**'s axis, so a row reads straight across
  the block. Colour = which positive (theme steps, since these are conditions and never a
  route colour, COLORS.md rules 2–3), fill = which negative. The key names the one matched
  pairing. This panel exists because the four ratios are **not four measurements of one
  quantity** (below).

### The x axis is the 32 designs' — changed 2026-11-19

Leo: *"B and c need to be built like fig 5b from the main figures. You should cut the
controls and super large deviations on the x axis."* Two things were setting the extent
and neither is a candidate effect:

| | was on the axis | now |
|---|---|---|
| Csy4-H29A 100 ng / 10 ng | dots at 10.81× and 96.68× | **printed as their number in their own row** (`11x` / `97x` in (b), `10–12x` / `79–117x` in (c)), as panel (a) already does |
| binder 4 posB, binder 6 posA (failed wells) | pulled a pairing ratio to 0.059 | **clamped on the left edge with a caret**, (b) labelling binder 4's averaged `0.48x` |

The axis ran **0.05–200** and the 32 designs — 0.478 to 1.563, a factor of 3 — occupied
about a fifth of it as a vertical smear. It now runs **0.492–2.780**, Fig5a's rule
(candidate ratios ≥ 0.5, padded ×/÷ 1.35) computed once and shared by (b) and (c). **No
statistic moved**; every number this README quotes is printed unchanged by `report()`,
which now also prints the limits and what landed outside them.

Three deliberate departures from Fig5a, all of them because this panel carries more:
1. **No gridlines.** Fig5a draws vertical ones; STYLE_NOTES §7 allows them only where a
   value must be read off precisely, which a ranked dotplot never needs. The supplementary
   panels follow the house rule, not the main-text panel.
2. **The 1.2 cut is a dotted rule, not Fig5a's two-tone background.** The background fill
   is spent on the ±2 SD null band, and two pale fills in one frame read as one thing.
   The band is the evidence this figure exists for, so it keeps the slot.
3. **The band is NEUTRAL 100 with a hairline at each edge**, not the NEUTRAL 200 fill it
   was. On the clipped axis the band covers two thirds of the frame; at 200 it read as the
   panel background and took the contrast out of every dot on it.

Also in (c): the carets carry **no** value label. A label long enough for both off-scale
numbers runs into the same row's in-range dots (binder 4's nearest sits at 0.97), and
these are dead-well readings, which the caption names.

### The layout — changed 2026-11-20

Leo: *"the arrayed flow, you can make the A horizontal, and then B and C, you can really
expand the x-axis so that they're almost like a square size plot, each of them. Because
then it also allows to read a bit better. So feel free to expand it horizontally."*

**No data, no statistic and no scale moved** — `report()` prints the same numbers it
printed before the change, line for line. What changed is orientation and width.

| | was | now |
|---|---|---|
| (a) | 46.7 × 102.7 mm, a 34-row column block beside (b)/(c) | **152.8 × 37.2 mm**, 34 columns × 4 target rows, on its own book row |
| (b) | 66.7 mm wide, plot box 60.9 × 81.8 (w/h **0.63**) | **86.8 mm**, plot box **71.2 × 81.8** (w/h **0.87**) |
| (c) | 62.9 mm wide, plot box 60.9 × 80.5 (w/h **0.57**) | **86.8 mm**, plot box **71.2 × 80.5** (w/h **0.88**) |

The height of (b)/(c) is fixed by the 34 rows at the 2.4 mm pitch that 5.5 pt type floors,
so the only lever on the aspect is width, and the only width available is the book's:
(182 − 6) / 2 = 88 mm each, side by side. 0.87 is as square as this panel can be made
without dropping rows or type size.

**What the width was actually being spent on.** Not margins in general — one word. The two
control rows were labelled `Csy4-H29A 100 ng` / `Csy4-H29A 10 ng`, **16.1 mm** of y tick
label on an 86 mm panel, and the same string was the single largest item in the left
margin of both (b) and (c). Split into a dose on the row (`100 ng`, 6.1 mm) and the name
once beside the pair, the plot box went from 60.9 to **71.2 mm** — **two thirds of the
total gain came from that**, not from the extra 20 mm of canvas. Same two labels, two
tiers, nothing dropped. (a) uses the same split under its 3.9 mm columns.

**The labels are measured, never estimated** (`style.text_w`; the usual 0.52 em/char
stand-in is out by 10–33% and by a different amount per string). `check_columns()` runs on
the laid-out figure, because the column pitch is not knowable before tight_layout has seen
the labels:

| | measured | fits? |
|---|--:|---|
| column pitch, (a) | **3.92 mm** | — |
| widest binder number (`25`) | 2.02 mm | upright, 1.90 mm of clearance |
| narrowest (`1`) | 0.51 mm | upright |
| `100 ng` | 5.59 mm | **over the pitch → rotated 90°** |
| `Csy4-H29A` alone | 9.84 mm | over the pitch at any size this panel could have |
| `Csy4-H29A 100 ng` whole | 16.06 mm | 16 mm of hanging label under a 3.9 mm column |
| tightest rendered gap, all 34 | **1.75 mm** | asserted > 0 on the real boxes |

So the 32 binder numbers stay **upright**, which is what a reader wants on a number, and
only the two control doses turn. The last check is the catch-all: the rendered
`get_window_extent` of all 34 labels, in order, asserted not to overlap — which is what
caught the one collision this revision did have (below).

**Three traps hit and fixed, all visible only in the render:**

1. **An inset colourbar + `tight_layout` is a feedback loop.** The bar is positioned in
   *axes fractions* above the axes, so shrinking the axes to make room moves it down and
   the layout shrinks again. It flattened the grid to ~4 mm and stacked the four target
   labels on top of each other. Fixed by laying out into `rect=[0, 0, 1, 0.78]` *before*
   the bar exists, which reserves the top band outright.
2. **The rule had vanished.** (b)/(c) draw it on `split` because nothing else is at that
   y; in (a), `split` is the boundary the first control *cell* is drawn on, and all eight
   control cells are outlined (off-scale by construction), so the rule sat under an INK
   edge. Moved to the middle of the blank column.
3. **A rotated label is aligned after it is rotated.** `ha="center"` on the rotated
   `Csy4-H29A` put half its ink height back over the gap — the `4` of Csy4 sat on the `1`
   of `100 ng`, in both (b) and (c). `ha="right"`.

A rotated colourbar label was clipped in the first attempt too: at 17.7 mm it is longer
than a bar bounded by the 17 mm grid height, so the bar is horizontal, above the grid at
the left end, where the panel's free space actually is.

## Verified — every number re-derived by the script, printed on every run

- **34 = 32 candidates + 2 Csy4 controls.** The export holds **35** plate columns; `P6_C3`
  is dropped by both of `plot_on_off.py`'s tests — BFP+ 0.093 / 0.00367 / 0.025 / 0.028 %
  against a kept-column median minimum of 30.8 %, *and* no `Normalized mScarlet` in its
  posB and NegA wells (`561(D)-A+ | Freq. of Parent` = 0, mean blank). 35 − 1 = 34 = 32
  rows carrying a `design_id` + `P6_C4` / `P6_C5`, which are blank. **P6_C3 was never a
  33rd candidate**: `candidates32_index.tsv` maps the 32 to `P1_C1…P6_C2`, and
  `lv2_33_screen.py:18` says "Plate 6 col 3 untransfected". So Methods V2's 34 and the
  main text's 32 are both right and count different things; the 34 is 32 designs + one
  plasmid at two doses.
- `ONOFF_avg` re-derives from the four raw condition columns to `rtol` 1e-3 (the stored
  columns are rounded to 5 dp and the smallest wells are ~0.003); against the unrounded
  xlsx the agreement is 4.6e-6.
- **26/32 (81%) above 1.0**, max **1.5627** (binder 23 = `P4_C5`), median 1.1859, IQR
  1.0432–1.3398. The six below are 0.478, 0.700, 0.842, 0.856, 0.919, 0.992 — none at
  exactly 1.0, so `>` vs `≥` does not move it. **15/32 reach the 1.2 cut** (proposed
  binders); the 15 are strictly inside the 26, so the two numbers are one column at two
  thresholds, not two findings.
- Exact two-sided binomial, 26/32 vs 0.5: **p = 5.35e-4**. Exact Wilcoxon signed-rank,
  `ON_avg` vs `OFF_avg`, n = 32, W+ = 455: **p = 1.62e-4**. Both reproduce, to the one significant figure they are
  quoted at, the prose values in `LV_writing/style/thesis-handoff.md:183`
  ("p=0.0002", "sign test p=0.0005"), which no script had ever computed.
- **The empirical null is the headline caution.** SD of log2(no-scaffold/NegA) = 0.302,
  i.e. **1.23×**, so ±2 SD is 1.52× — wider than every candidate effect but two.
  **2/32 sit above +2 SD and 1/32 below −2 SD.** The median candidate (1.186) is inside
  1 SD of a contrast with no binding in it. The panel-level shift is significant; no
  individual design is.
- **The two negatives disagree by median 1.083× (range 0.742–1.585).** Choosing the
  negative moves a variant's ON/OFF by about as much as the biology does.
- Controls: `P6_C4` = **Csy4-H29A 100 ng, 10.81×**; `P6_C5` = **10 ng, 96.68×**
  (`lv2_33_screen.py:45`, `VARIANT[(6,4)] = "Csy4 100 ng"`, `[(6,5)] = "Csy4 10 ng"`). The
  lower dose gives the larger ratio, through a 12× lower negative baseline (0.00465 vs
  0.05370), not a higher positive. ⚠️ The dose names are the plan's: the same 2.40 → 2.00 µL
  deviation was applied to both control wells "so the positive control carries the same
  dose", so the delivered amounts are **84 ng and 8.4 ng**. Panel labels keep the plan's
  names; the caption should not quote 100/10 ng as delivered.
- Routes among the 32: BG+LiMPNN 11, BG+BIF 9, BD1+LiMPNN 7, RFD3+LiMPNN 4, HD+LiMPNN 1.
- Plasmids, counted from sequence (circular-aware motif count of the 16 nt hairpin core
  `GTTCACTGCCGTATAG` and of the 19 nt replacement `ATTTCGATCAATGGATTAT`, not from the
  annotations): pLV13 posA 6,479 bp / **16** sites / 8 F30 three-way-junction units;
  pLV14 posB 5,951 bp / **8** sites / **no F30 feature at all** (one de novo scaffold);
  pLV15 5,879 bp / 0 sites / **8** replacement loops, no F30; pLV16 NegA 6,335 bp / 0 sites
  / **16** replacement loops / 8 F30 units; pLV17 5,551 bp, no scaffold.

## The three things to settle

### 1. 34 against 32 — closed, see above

34 = the 32 ordered candidates (each with a `design_id`) + **pLV18, the catalytically
inactive Csy4-H29A binder, at 100 ng and 10 ng**: one plasmid, two doses. A 35th column
was seeded and acquired untransfected and is not part of the 136 transfected wells.

### 2. The pLV15 / pLV16 substitution — what it does to panel (c)

Established from the maps, by length arithmetic and by motif count, not from labels:

- **pLV15 is posB's matched negative.** pLV14 (posB) − pLV15 = **72 bp = 8 × 9**, which is
  exactly 8 sites each losing 9 nt (28 nt repeat → 19 nt replacement); pLV15 carries 8
  replacement loops, 0 hairpins, and no F30 unit — i.e. pLV14's single de novo scaffold
  with its 8 hairpins knocked out.
- **pLV16 is posA's matched negative.** pLV13 (posA) − pLV16 = **144 bp = 16 × 9**; pLV16
  carries 16 replacement loops on 8 F30 three-way-junction units.
- pLV15 was never transfected — and **it was never planned**. `pLV15` occurs **0 times**
  in `LV2.33_transfection_plan.html` and 0 times in
  `LV2.33_transfection_AS-PERFORMED_2026-11-04.html` (same for the `~/Downloads` copies).
  The planned four master mixes are MM1–MM4 = pLV13 / pLV14 / pLV16 / pLV17, and every
  variant's row C is NegA/pLV16 whichever positive it is paired with. The as-performed
  deviations table lists only pipetting volumes; no plasmid swap. **So this is a design
  decision, not a bench slip** — "pLV16 went in by accident" is not what the record shows,
  and the figure should say the experiment was designed with one negative scaffold for two
  different positive scaffolds.

So **posB was read against a negative built on the other target's scaffold**. `PosB/NegA`
divides an 8-site single-scaffold transcript by a 16-loop F30-array transcript: numerator
and denominator differ in scaffold architecture, in length and in loop number, not only in
whether the Csy4 hairpin is there. That ratio is therefore **not** an ON/OFF in the sense
`PosA/NegA` is, and the two must not be averaged into one number without saying so.
`ONOFF_avg`, which the main ranking uses, does exactly that average — it is kept here
because it is the published statistic, and panel (c) is what keeps it honest.

What it costs, measured:

| pairing | > 1.0 | ≥ 1.2 | median |
|---|--:|--:|--:|
| posA / NegA **(matched)** | 27/32 | 15/32 | 1.188 |
| posA / no scaffold | 27/32 | 20/32 | 1.239 |
| posB / NegA **(unmatched)** | 20/32 | 7/32 | 1.065 |
| posB / no scaffold | 23/32 | 13/32 | 1.129 |

The count of proposed binders runs 7 to 20 out of 32 depending on which pairing is read.
An upper bound on the damage: the two negatives that *were* measured differ per variant by
a median 1.083× (0.742–1.585), so swapping one inert scaffold for another moves the
denominator by up to 1.6× — the size of the entire candidate signal range.

⚠️ **And the negative that was used is itself half missing.** The as-built pLV16 (clone
`pLV2.25.4_c`, Microsynth Full PlasmidSeq order 4909870, 27-OCT-2026) is **5,943 bp with 8
replacement loops in 4 F30 units** against a design of 6,335 bp with 16 loops in 8 units —
a 392 bp deletion, counted here independently of the annotation (the 19 nt replacement
occurs 8× in the as-built and 16× in the design). So NegA is not architecture-matched to
posA's 8 units either; it matches the *family*, at half the array. pLV13 and pLV14 have no
as-built record at all, so their arrays are the designed ones by assumption.

### 3. `NEED:` the as-built loop count of pLV16 — **closed, and the brief's caution is the
wrong one**

The brief says the map file carries unannotated deletions, so no copy count should be
quoted from it. That applies to the **design** map read naively; it does not apply here,
because the count above comes from counting the replacement motif in the ORIGIN sequence
of both files, circular-aware, which needs no annotation. The as-built record also rules
out assembler repeat collapse on its own numbers (mean read depth over the array 642.7×
against 669.1× elsewhere, ratio 0.961; no insertion signal above 1.09%).

**One caveat stands:** 11 of 29 uncertain positions fall inside the array, 8 of them one
per loop at the same internal offset. The loop *count* and array period are reliable;
**an individual loop sequence from this assembly is not**, and must not be quoted as
verified without Sanger across it.

## `NEED:` and what did not reproduce

- 🚨 **"roughly 6% of the low-dose Csy4 control" does not reproduce, in either direction.**
  The low dose is 10 ng = `P6_C5` = 96.68× (`lv2_33_screen.py:45`), and the largest
  candidate, 1.5627×, is **1.6%** of it (0.6% of its lift above 1). Against the 100 ng
  control, 10.81×, it is **14.4%**, and **5.7% of that control's lift** — the only reading
  anywhere near 6%, and it is the *high*-dose control on a lift basis. Panel (b) therefore
  puts both doses on the axis and quotes no percentage; the caption gives 1.56× against
  10.81× and 96.68×. `NEED:` which quantity "6%" was meant to be.
- 🚨 **The two positives do not agree, so the context-independence argument is not
  available.** Across the 32, log2(posA/NegA) vs log2(posB/NegA) gives Pearson **+0.015**,
  permutation p = **0.89** (20,000 draws, hand-rolled). On the no-scaffold negative it is
  +0.279, p = 0.044 — not robust to the choice of denominator. `lv2_33_stats.py` finds the
  same thing on the per-cell readout (p = 0.203). A design responding to both geometries
  would be the strongest specificity evidence here; **the data do not support that any
  design does**, and the figure must not imply it.
- ⚠️ **posA and posB do not "differ only in geometry".** Same 16 nt Csy4 hairpin core, but
  **16 sites on 8 F30 three-way junctions** against **8 sites on one de novo scaffold** —
  valency differs 2× and the scaffold is a different molecule (`MAINTEXT_DATA_REQUESTS`
  A/B: pLV14's insert is pART252, pLV13's is pART158). Two presentations, differing in
  valency *and* architecture.
- ⚠️ The FCS file names say `85ng`; the as-performed record says **84 ng** (2.0 µL instead
  of the planned 2.4 µL), uniform across all binder and control wells. 85 is the nominal
  label on the run, 84 the delivered amount. Uniform, so it does not bias between
  variants.
- ⚠️ **An independent analysis of the same wells reaches a weaker conclusion.**
  `lv2_33_stats.py` works from raw FCS with BFP-matched per-cell statistics (not the
  FlowJo per-well normalized mean this figure uses) and finds **1 of 19 variants at
  BH q < 0.05 — V18, in the wrong direction**. It covers only the 19 variants whose four
  wells survived that export. Both analyses agree that no candidate is individually
  distinguishable; they disagree on nothing. Worth a Methods sentence saying which
  readout the figure uses.
- ✅ One correction in `SI_BUILD_STATUS_2026-11-18.md` is already spent: it asks for
  `inputs/MANIFEST.md` to stop calling `P6_C3` a failed transfection. The string "failed
  transfection" does not occur anywhere in `MANIFEST.md`; its entry for this CSV states the
  drop rule and the four BFP+ readings without a cause. Nothing to fix.
- The pooled-screen comparison a design was selected on does not predict which arrayed
  target it responds to: the 25 picked on POS_B give a median posB/NegA of 1.067 against
  the 7 picked on POS_A at 1.059, while their posA/NegA medians are 1.219 and 1.108.
- `NEED:` the **plate-column → binder-identity mapping** is not re-derivable from any file.
  `inputs/MANIFEST.md` says so outright — it is the as-performed plate map plus assertion.
  Everything in this figure that names a design rests on it.
- `NEED:` an as-built sequencing record for **pLV13** and **pLV14**. Only pLV16 has one,
  and it turned out to be half-deleted.

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the n. 211 words.*

**Arrayed re-test of the 32 ordered candidates.** The panel sits above its no-effect line;
no individual design clears the experiment's own noise floor. (**a**) Normalized mScarlet
(mScarlet⁺ mean × %mScarlet⁺ ÷ BFP⁺ mean × %BFP⁺) per well, as log2 of each variant's
four-target mean; variants ranked by ON/OFF left to right, Csy4-H29A past the rule. Outlined cells exceed
the ±0.6 scale; binder 4 posB and binder 6 posA are failed wells. (**b**) ON/OFF,
avg(posA, posB)/avg(NegA, no-scaffold); line = range over the four positive × negative
pairings, band = ±2 SD of log2(no-scaffold/NegA), where no binding differs.
26/32 (81%) exceed 1.0 (exact binomial p = 5.4 × 10⁻⁴; Wilcoxon p = 1.6 × 10⁻⁴) 15/32 (47%)
reach 1.2; only 2/32 (6%) leave the band. (**c**) Each pairing alone:
only posA/NegA pairs a target with a negative on its own scaffold; the count above 1.2
runs 7/32 (22%) to 20/32 (63%). (**b**, **c**) are scaled to the designs: Csy4-H29A is
printed in its rows (10.8×, 96.7×); ◀ marks an off-scale value (binder 4, 0.48×; failed
wells, 0.06–0.16×). mScarlet reports receptor dimerisation on a multivalent
scaffold — an avidity proxy, not affinity (no Kd). One biological replicate, one well per
variant × target, 84 ng binder plasmid. The 32 are a ranked selection, not a hit rate.

## Open

- 🚨 **Overlap with `Fig5a_wetlab_on_off` — now close enough to need a decision.** Making
  (b) look like `Fig5a_PosAvg_vs_NegAvg` is what Leo asked for, and it worked: same 32
  dots, same order, same colours, same axis, **the two panels are now near-identical in
  everything except what (b) adds**. What (b) adds is real and is not in the main-text
  panel: the ±2 SD empirical null band (the "no individual design clears its own noise
  floor" claim lives there and nowhere else), the per-design range over the four pairings,
  and the two controls as a scale reference. The main-text panel adds the ≥ 1.2 shading
  and the two band labels, and is 100 mm wide against 66 mm here.
  **Recommendation: the main text owns the ranked view** — it is the result, and Fig 5
  is where a reader meets these 32 — **and the supplementary (b) is read as the noise-floor
  version of it**, which both captions must say in one clause ("the same ranking against
  the experiment's empirical null, Supplementary Fig. X b" / "the main-text ranking,
  Fig. 5 b"). If only one may survive, keep this one: the ranking without the null band
  overstates what the experiment can resolve. Do not ship both silently.
- **(c) and the pass/fail matrix.** `Fig5a_pass_matrix_gate*` and panel (c) are two views
  of the same four pairings — the matrix thresholds them, (c) keeps the values. Decide
  which figure owns it; (c) is the more informative of the two.
- **`ONOFF_avg` averages a matched and an unmatched pairing.** It is kept because it is the
  statistic the selection and Fig 5 are built on. If the text ever wants one number per
  design, posA/NegA is the defensible one, and the 15/32 at the 1.2 cut happens to be the
  same count either way.
