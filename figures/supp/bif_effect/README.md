# bif_effect

What the BoltzIF scoring defect did to the delivered panel, and to the design route
behind it. Three bare panels, cited from Section 3.

```bash
~/.venvs/binder_map/bin/python3 code/supp/bif_effect/make_figure.py
```

Emits `output/bif_effect/panel_{a,b,c}_*.{pdf,png,svg}`. Authored boxes (a) 85 × 46 mm,
(b) 85 × 38 mm, (c) 177.8 × 86.4 mm; **emitted** 87.6 × 48.1, 88.3 × 43.0 and 177.8 ×
86.4, because (a) and (b) save with the default `tight=True` and their top annotations
overflow the authored box — only (c) goes out with `tight=False`. Every number below is
re-derived and printed on every run, and `check_emitted()` re-measures all three
**emitted PDFs** at the end of every build: no text off any MediaBox, and the (a)|(b)
pair still fitting one row of the book.

## ⚠️ (b) and (c) swapped letters, 2026-11-20

Leo: *"move c up next to a and relabel it."* The ipTM-decile curve moves up beside (a)
and is now **(b)**; the four-quadrant route density is now **(c)** and stays, unchanged
in every respect. The book letters panels by **file-name order**, so the letter lives in
the emitted stem and nowhere else — `panel_c_flip_by_iptm_decile` →
**`panel_b_flip_by_iptm_decile`**, `panel_b_stored_vs_corrected` →
**`panel_c_stored_vs_corrected`**. The two stale stems were deleted from `output/` in the
same pass, because `make_supp_book.py` globs `output/<slug>/*.pdf` and leaving them
would have put five panels on the page.

🚨 **The side-by-side row has 0.1 mm of slack — now asserted, not just noted.** (a) and
(b) are emitted at 87.6 and 88.3 mm, and the book's text block is 182 mm with a 6 mm
gutter: 87.6 + 6 + 88.3 = **181.9 mm**. They share a row, but anything that widens either
tight bounding box — one longer tick label is enough — would drop (b) onto its own row
and undo the ask **silently**, during some later edit, and get blamed on something else.

`check_emitted()` closes that off. It runs at the end of every build and asserts
`(a) + 6 + (b) ≤ 182`, so the reflow cannot happen without the build stopping and saying
why. The threshold is the **real limit**, not `182 − ROW_SLACK_MM`: the pair emits at
181.9 mm today, so a slack-based assert would fail every run. The thin margin is
**warned** about instead, on every run:

```
  ⚠️  (a) + 6 mm + (b) = 181.9 mm leaves 0.1 mm in the 182 mm text block, under the
      1.5 mm we want. One longer tick label and the row breaks - see README.
```

Buying the slack means changing panel geometry: the authored boxes (85 + 6 + 85 =
176 mm) have 6 mm, and `tight=False` would take them — and **would also clip** the
`226/240` and `8.7` annotations that sit outside those boxes, which is why both panels
are tight in the first place. `NEED:` Leo's call on whether to spend that.

**Populations differ by panel, deliberately.** (a) and (b) are the 2,397 rescored
members of the ordered panel. (c) is the **whole BoltzGen + BoltzIF design route**,
1,435,617 designs, drawn bare — a four-quadrant density with nothing overlaid on it.
See "Why (c) is the route".

## The defect

BoltzGen + BoltzIF's stored `dockq_*` and `sc_rmsd_*` were scored against BoltzIF's own
Boltz **re-prediction**, not against the generated backbone. The `*_true` columns are
the correction. **STORED** is campaign history — what the gate actually saw, and
therefore what was ordered — and **CORRECTED** is what those backbones really are.
They are never mixed inside one number here, both panel (c) axes are labelled with
which one they carry, and the panels keep those two words.

## Sources

| file | mtime | what it gives |
|---|---|---|
| `inputs/large/bif_true_per_backend.csv.gz` | 2026-11-17 15:18 | `rmsd_{af3,esm}_true`, `dockq_{af3,esm}_true`, 1,480,140 rows |
| `inputs/large/bif_route_stored_per_backend.csv.gz` | 2026-11-19 | the route's **STORED** `rmsd_{af3,esm}`, `dockq_{af3,esm}`, 1,480,140 rows |
| `inputs/csy4_binder_library_12k_v6_k30.csv` | 2026-11-12 14:07 | the delivered 12,000-design panel; 2,400 are BoltzGen + BoltzIF |

Joined on `design_id`, one-to-one, asserted in the script. All mtimes confirmed on
disk. `inputs/MANIFEST.md` carries each file's own provenance.

## Why (c) is the route and not the ordered panel

Leo, 2026-11-19: *"B i want to see all quadrants."* **Labelling the old panel could
not have produced them.** It drew the ordered panel, and the ordered panel *is* the
gated set, so every point passed on the stored axis by construction: stored
best-of-two ligand_RMSD topped out at exactly **3.0000** and stored sc_DockQ bottomed out at
exactly **0.7000**. The gate wall was the edge of the data; two of the four quadrants
held no points at all, and no annotation fixes an empty quadrant. The fix is the
population.

**Finding the stored half of the route.** `bif_true_per_backend.csv.gz` carries
corrected scores for all 1,480,140 route designs but no stored ones, and the library
CSV carries stored scores for the 2,400 ordered designs only — which is precisely why
the panel had no stored-fail point to draw. The route's stored per-backend scores are
in `large/campaign_metrics_full.csv.gz`: its **consensus** `sc_dockq` / `sc_rmsd` are
the rescore, but its **per-backend** `dockq_{af3,esm}` / `rmsd_{af3,esm}` were left
alone on purpose — `export_campaign_metrics.py`, "the raw per-backend columns are kept
as-is so the correction stays visible". Verified rather than trusted: against the
library's own stored columns on all 2,400 ordered designs, max |d| is **one 4-dp
rounding unit** (0.0005 Å, 0.0001 sc_DockQ). One `awk -F,` pass over it produced
`inputs/large/bif_route_stored_per_backend.csv.gz`.

**Population 1,435,617.** 1,480,140 route designs; 1,436,724 carry a true-backbone
score (the `bif_rescored` flag); 1,435,617 of those also carry a complete stored
quadruple. All 2,397 rescored ordered designs survive, so **(c)'s population contains
(a)'s exactly** rather than merely overlapping it. Asserted in `load_route()`.

**No overlay, and no per-quadrant "ordered" line (2026-11-19).** Leo: *"there's just
too much going on now. So please remove all the dots. Just make it a clean hex plot
with the 4 quadrants, and then the labels are enough."* The `INK` layer that marked
the 2,397 ordered designs is gone, and the count line that went with it
(`1,896 ordered` / `none ordered`) went too, in all eight blocks: a label naming
marks that are no longer drawn is worse than no label, and keeping it in one block
only would read as an exception rather than a result. Each block is now two lines,
its quadrant name and its `n/N (%)`.

**The "none ordered" fact is not lost, it moved to the caption**, which is where Leo
put it. Every ordered design is a stored pass by construction, so neither stored-fail
quadrant can hold one — **`Falsely rejected` is 11,986/1,435,617 (0.835%) and none of
it was bought**. `report_route()` still counts the ordered members of all eight
quadrants on every run and prints them, so the number cannot drift out of the caption
unnoticed; the tables below are that output.

## Verified against the files

- **Denominator 2,397 for (a) and (b).** All 2,400 panel members are present in the
  true-backbone table; three carry empty values. They are **not one thing**, and
  calling all three "no rescore" understates the result:
  - `csy4free_BG_BIF_len80_00946892_design6_seq3` (v2) — `unpaired_no_limpnn_uuid`,
    genuinely unmeasurable: no LigandMPNN sibling, so the true backbone does not exist.
  - `csy4free_BG_BIF_len80_5532a9d9_design23_seq1` (v4) and
    `csy4free_BG_BIF_len80_787d9574_design7_seq5` (v7) — `no_native_interface`, which
    would **fail the gate's sc_DockQ term outright**.

  Reason codes read from the Euler source,
  `analyses/bif_true_bb_rescore/full/bif_true_bb_scores.csv.gz` (`true_bb_flag`).
  The denominator does not change the answer: 1,963/2,397 = 81.89%, 1,965/2,399 =
  81.91% (the two `no_native_interface` rows counted as flips), 1,965/2,400 = 81.88%.
  All round to **81.9%**.
- **1,963/2,397 (81.9%) of the ordered panel flip pass → fail.** Predicate, both sides
  identical: `(sc_rmsd_af3 ≤ 3 OR sc_rmsd_esm ≤ 3) AND (dockq_af3 ≥ 0.7 OR dockq_esm ≥
  0.7)`, on the stored columns and then on `*_true`; flip = stored pass AND NOT true
  pass. `stored_pass` is 2,397/2,397 (that panel is drawn from the gated set,
  asserted); `true_pass` is 434.
- **Medians are the POOLED AF3 + ESMFold2 convention** and the caption says so:
  ligand_RMSD **2.744 → 8.663 Å**, sc_DockQ **0.709 → 0.406**. Best-of-two — what the gate
  actually uses — is 2.181 → 5.799 Å and 0.742 → 0.553. Panel (a) plots the pooled
  convention, so it draws 4,794 paired segments (2,397 designs × 2 re-predictors).
- **Panel (c)'s four quadrants, per gate term, on the gate's own best-of-two
  convention** — `min(sc_rmsd_af3, sc_rmsd_esm)` and `max(dockq_af3, dockq_esm)`, one
  point per design, N = **1,435,617**. ⚠️ The `ordered` columns are **printed, not
  drawn** — they are `report_route()` output, and since 2026-11-19 nothing on the
  panel refers to them:

  | quadrant | ligand_RMSD | ordered | sc_DockQ | ordered |
  |---|--:|--:|--:|--:|
  | Both rejected | 1,396,698/N (97.3%) | 0 | 1,404,151/N (97.8%) | 0 |
  | Falsely accepted | 20,342/N (1.4%) | 1,896 | 11,248/N (0.8%) | 1,793 |
  | Both accepted | 3,981/N (0.3%) | 501 | 3,517/N (0.2%) | 604 |
  | Falsely rejected | 14,596/N (1.0%) | **0** | 16,701/N (1.2%) | **0** |

  🚨 **These are PER TERM and the two terms overlap, so no one of them is the combined
  flip.** The ordered columns reproduce the ordered panel's own per-term numbers
  exactly: 1,896/2,397 (79.1%) fail corrected ligand_RMSD, 1,793/2,397 (74.8%) fail
  corrected sc_DockQ, overlapping on 1,726.
- **Both structural terms together, same population** — the numbers the caption quotes:

  | quadrant | n/N | (%) | ordered |
  |---|--:|--:|--:|
  | Both rejected | 1,411,499/1,435,617 | 98.320% | 0 |
  | Falsely accepted | 9,803/1,435,617 | 0.683% | 1,963 |
  | Both accepted | 2,329/1,435,617 | 0.162% | 434 |
  | **Falsely rejected** | **11,986/1,435,617** | **0.835%** | **0** |

- 🚨 **`Falsely rejected` is a number nobody had produced.** Designs the
  campaign discarded on stored scores that the corrected scores accept:
  **11,986/1,435,617 (0.835%)** on both terms, **none of them ordered**; per term
  14,596/1,435,617 (1.0%) on ligand_RMSD and 16,701/1,435,617 (1.2%) on sc_DockQ. The
  corrected gate passes **14,315** designs against the stored gate's **12,132**, so the
  defect cost this route more than it gained it — the same direction as the corrected
  funnel below, where BoltzGen + BoltzIF goes 9,042 → 9,515.
- **Cross-checks, asserted on every run** rather than compared hopefully, each against
  the population it was counted on. All come from files built by other scripts:
  - stored best-of-two ligand_RMSD pass = **24,401** over all 1,480,140, exactly
    `venn_regions_by_stack.csv` for BG+BIF (full-precision store floats, a different
    code path). 24,323 of those are drawn; the other **78** carry no complete pair of
    scores. The assert is taken in `load_route()` *before* the population is narrowed,
    so it compares like with like.
  - corrected best-of-two ligand_RMSD pass = **18,577**, exactly the corrected cascade's
    `g1_sc_rmsd` below. It lives in `report_route()` because the cascade counts only
    designs that *have* a corrected score, which is this population exactly.
  - the cascade's step ratios fall out too — stored 12,162/24,401 = **49.8%**,
    corrected 14,315/18,577 = **77.1%**, both as already published.
  - stored both terms is **12,162** over the full route against venn's **12,158**, and
    **12,132** drawn. It reconciles exactly: **+4** from the campaign master's 4-dp
    sc_DockQ rounding (23 route designs land on exactly 0.7000) and **−30** stored-gate
    passers with no true-backbone score to plot.
- **Uniform across rounds: χ² = 6.2672, df = 6, p = 0.394.** Seven rounds, no v3.
  survive/flip/total — v2 39/195/234, v4 19/118/137, v5 35/123/158, v6 26/107/133,
  v7 92/449/541, v8 126/503/629, v9 97/468/565. Smallest expected count 24.08, so no
  exact test is needed.
  ⚠️ **This result has no panel.** It is **stated in the caption with its statistic**,
  and `main()` still derives it on every run and `report()` still prints the table, so
  it cannot drift. It is what shows the defect is one shared mechanism rather than a
  bad batch, so it must not be dropped from the caption.
- **ipTM.** The "gap 0.027" holds **only** as a pooled AF3 + ESM mean difference
  (+0.0271). Per re-predictor: `iptm_af3` +0.0189, `iptm_esm` +0.0354, `iptm_max`
  +0.0204. The caption quotes `iptm_max` **0.020**, the gate's own quantity, and names
  it.
- **Top decile 74–77%, against 81.9% overall.** "74%" is AF3-specific: `iptm_af3`
  74.92% (227/303 at ≥ 0.8200 — n = 303 from ties), `iptm_esm` 76.67% (184/240),
  `iptm_max` 76.65% (197/257). Panel (b) draws `iptm_max` deciles; the flip rate falls
  from 94.2% (226/240) in the lowest decile to 76.7% (197/257) in the highest, so ipTM
  carries a little signal and buys almost no protection — it never reads the designed
  backbone.
- **`VERIFY:` the corrected funnel is 56,179, and it checks out.**
  `/Users/leo/rna_binder_plots/_data/cascade_by_stack.csv` (2026-11-03 10:40) carries
  both variants, and only `old_bif` sums to the quoted 55,706:
  `old_bif` ALL g5_iptm 55,706 = 27,881 + 14,752 + 9,042 + 3,229 + 802;
  `new_bif` ALL g5_iptm 56,179 = the same four routes with BoltzGen + BoltzIF at 9,515.
  The whole 473-design difference is that one route, 9,042 → 9,515. It goes **up**
  because the cascade steps are cumulative over different subsets: on the true
  backbones the route loses more at `g1_sc_rmsd` (24,401 → 18,577) and then keeps far
  more at `g2_dockq` (49.8% → 77.1% of the previous step).

## Not drawn, deliberately

- **The 3 Å and 0.7 hairlines in panels (a) and (c) ARE drawn**, against STYLE_NOTES §0's
  default that a threshold is caption work. Crossing them is the flip: without them (a)
  shows two distributions moving and (c) has no quadrants at all. Leo asked for them on
  both axes of (c) by name.
- **No colorbar on (c).** The ramp only has to show where designs are; the quantities
  are the eight quadrant counts, and the caption states that darker is more designs on
  a log count scale, per panel. A colorbar would repeat what the labels already say.
- **(c) is the house DOUBLE-column width, not 85 mm.** Eight labelled quadrants do not
  fit in 85 mm, and that was **measured**, not guessed: the widest count line is 20.9 mm
  and sc_DockQ's stored-pass strip is 13.7 mm at 85 mm against 23.3 mm at 178 mm. Equal
  aspect is kept — y = x stays a true 45° and the two hairlines cross on it — so the
  axes are square and the figure **height is derived from the width** (86.4 mm, read
  off the emitted PDF's MediaBox; the README said 85.9 until 2026-11-19) rather
  than chosen. `_fits()` re-measures all eight blocks with `style.text_w` on every run,
  in the weight each is actually drawn in, and asserts, so a relabelled quadrant cannot
  reach the PDF clipped.
  sc_DockQ's axis runs to **1.0**, not the old 0.95: its natural ceiling, and 0.95 left the
  stored-pass strip at 26% of the axis where a 20.9 mm line needs 30%. ligand_RMSD runs to
  **70**, not 65 — one route design sits at 66.1.
- **"Falsely rejected" is set bold, and nothing else is.** It is the figure's result
  and Leo asked for it to be unmissable. Weight, not colour, so no palette decision is
  smuggled in and it survives a greyscale print. The label is Leo's own 2026-11-19
  wording, replacing "should have been ordered"; the rename keeps that label's logic
  intact — it is a counterfactual about designs that were *not* ordered, so it stays
  true of every design in the quadrant however large the population grows. The other
  three are unchanged: "Both rejected", "Falsely accepted", "Both accepted".
- **The quadrant blocks sit on a soft `SURFACE` patch.** A quadrant's outer corner is
  not always empty — sc_DockQ's "both rejected" corner is the densest part of the figure —
  and without it two of the eight blocks are unreadable on near-black hexes. Legibility,
  not decoration. A patch and **not** a path-effects halo, which was tried first and
  rejected: a stroked glyph has to be rendered as an outline, which silently turned
  these eight blocks into paths in the SVG. `style.save()` writes SVG precisely so
  assembly gets live text it can restyle, so the elements that most need to stay text
  cannot be the ones outlined. Verified on the emitted SVG: **38** `<text>` elements,
  16 of them the eight blocks' two lines each (46 and 24 when the blocks were three
  lines).
- **No "ordered" clause and no "ordered" count on any of the four labels.** Leo's
  original wording had two of them ending "and ordered", which did not survive widening
  the population to the route: 7,840 of the 9,803 designs the stored gate falsely
  accepted were never ordered, so the clause would be false of most of that quadrant.
  Those two were carried instead by each block's third line (`1,896 ordered` /
  `none ordered`) and by the marked points — and when the marks went, 2026-11-19, the
  third line had nothing left to point at and went with them. The one quadrant whose
  label is itself about not being ordered keeps saying so in words ("Falsely rejected"),
  and its `none ordered` is in the caption.
- **A scatter, for (c)'s 1.4M points.** Rejected: the three quadrants that carry the
  result hold ~1% each and vanish at any alpha low enough to stop the 97% corner going
  solid. `style.density_cmap()` / `style.density_norm()` are the shared house ramp and
  count scale, so this panel cannot drift from the main text's 2-D densities.
- **The ordered-design overlay itself, since 2026-11-19** — `ax.plot(a[o], b[o], …)`
  in `palette.INK`, deleted on Leo's *"remove all the dots"*. While it was drawn it took
  no new colour: a *selection* of the same series is not a second entity, so it was
  carried in `INK` exactly the way panel (a) carries its median (COLORS.md rule 1).
  Restoring it is two lines and the mask `o` is still computed — but the panel reads
  better without it, and the fact it carried is in the caption.
- **Panel (c) is best-of-two, not pooled like panel (a), and that is the whole reason a
  gate can be drawn on it.** On the pooled convention **1,969/4,794 (41.1%)** stored
  ligand_RMSD values and **2,058/4,794 (42.9%)** stored sc_DockQ values sit on the fail side of
  a line that every one of those designs passed, because the gate takes either
  re-predictor; drawing the gate there would show a rejection that never happened.
  Best-of-two is the quantity the gate actually read, so it is the only convention a
  drawn threshold is honest on. Panel (a) keeps pooled because its number is a pooled
  median, and the caption names which convention each panel uses.
- **The two re-predictors are not separated in panel (a) or (c).** Pooled and
  best-of-two are both conventions that by definition do not distinguish them, so no
  re-predictor colour (COLORS.md rule 7) applies here. Per-re-predictor values are
  above.
- **No hatch anywhere in this figure, and there never was.** The dot/stripe texture
  encoding for re-predictors was withdrawn repo-wide on 2026-11-19 (palette v0.6,
  COLORS.md rule 7); `grep -n 'hatch=' code/supp/bif_effect/make_figure.py` returns
  nothing. (Grep for the **kwarg**, not the word — these paragraphs match on the word
  and would answer for themselves.)
- **No y = x diagonal in panel (a)**, which has no x quantity to have one against;
  panel (c) carries it, as the reference the displacement is read from.
- **Panel (b) is a decile curve, not two ipTM densities.** The distributions overlap so
  heavily (means 0.020 apart) that two densities would show nothing; the decile curve
  answers the actual question — does high confidence protect a design — and carries the
  top-decile rate natively.
- No route colour anywhere: one design route, nothing for a route colour to
  distinguish, so every panel is a single series in `MAIN` (COLORS.md rule 1).

## 🚨 The caution that travels with the 81.9%

**1,963/2,397 (81.9%) is not a route-wide failure rate.** Those 2,400 were selected
precisely for scoring well under the wrong reference, which maximises expected
regression.

**The route-wide rate is 9,803/12,132 (80.8%)** — stored-gate passers across the whole
route that the corrected gate rejects, against the ordered panel's 1,963/2,397 (81.9%).
Widening panel (c)'s population is what finally allowed it to be stated. The two must
be kept apart because they are different quantities over different populations, but the
comparison is now possible and the selection effect turns out to be **1.1 points**. The
funnel number, 9,042 stored-gate survivors against 9,515 true, is the five-term version
of the same thing.

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the n. 193 words. The per-term "falsely rejected" counts are ON the panel, so the caption carries only the both-terms one - and since 2026-11-19 it also carries "none of it ordered", which the panel no longer marks.*

**Figure Sx. Scoring BoltzGen + BoltzIF against a re-prediction rather than the
generated backbone failed both ways: it wrongly bought four of every five ordered
designs, and discarded as many again that were sound.** (**a**) Stored and corrected
ligand_RMSD and sc_DockQ for the ordered designs, one line per design per re-predictor,
AF3 and ESMFold2 pooled; heavy line, the pooled median — 2.74 → 8.66 Å and
0.709 → 0.406; grey, the 3 Å and 0.7 gates. (**b**) Flip
rate by `max(AF3, ESM)` ipTM decile: confidence barely protects, ipTM never reading the designed backbone (gap
0.020). Uniform across seven rounds (χ² = 6.27, df = 6, p = 0.39). (**c**) Stored against corrected across
the design route, best-of-two as the gate read it, threshold on both axes; dashed,
y = x; darker, more designs (log counts, per panel). Quadrants name what each score
decided, counted per term. Every ordered design passes on the stored axis by
construction, so the stored-fail quadrants are the campaign's discard pile: falsely
rejected is 11,986/1,435,617 (0.8%) on both terms, none of it ordered. n = 1,435,617 of
1,480,140 route designs in (**c**), 2,397 ordered; (**a**) and (**b**) are those.
81.9% is that panel's rate, not the route's: route-wide, 9,803/12,132 (80.8%).

## ⚠️ The x label was clipped out of the PDF, 2026-11-20 — fixed

**(c)'s left x label printed as `Stored ligand RMSD (Å) ↓`.** The underscore hangs about
0.1 em below the baseline, matplotlib put the label's baseline low enough that it fell
**0.78 pt off the bottom of the MediaBox**, and the glyph was clipped out of the PDF.
The **PNG was always correct**, which is exactly how it survived an earlier renaming
sweep — and matplotlib's own renderer bbox put the label 0.5 pt *inside* the canvas,
because it measures the laid-out box and not the glyphs' descender metrics. The right
half's last x tick had the same disease: `1.00` is centred on the axes' right edge and
ran **2.81 pt** past it, so the PDF lost its final digit.

Fixed with the two levers that **move no geometry**, because this panel's height is
derived from its width under equal aspect and every margin lever would have moved the
documented 177.8 × 86.4 mm box:

- **`labelpad=2.5`** on both x labels (the default is 4.0), which lifts the label inside
  the canvas. Set on both halves so the two stay on one line.
- **`{:g}` tick labels** on the linear half: `0 / 0.25 / 0.5 / 0.75 / 1`, not
  `0.00 / … / 1.00`. `1` is 3.1 pt of type against `1.00`'s 10.7 and sits inside. It also
  makes this half agree with the log half, which has always written its ticks bare
  (`1 / 3 / 10 / 30`), with `repredictor_comparison` panel (a), and with STYLE_NOTES §7's
  "round numbers elsewhere". Applied to both axes of that half so they match.

Measured after: the emitted page is **unchanged at 504.00 × 244.80 pt = 177.8 × 86.4 mm**
and the worst text margin is **−0.72 pt** (inside), against **+2.81 pt** (outside) before.

**It cannot come back.** `emitted_box()` re-opens each panel's **PDF** after it is written
and asserts that every text span lies inside the MediaBox. The PNG is not the artefact of
record and is not what is checked. `repredictor_comparison`'s two panels pass the same
test (worst −1.09 and −0.20 pt), which is how I know this was local to (c).

## Palette

**Palette v0.9 (`config/palette.json`, 2026-11-20) changes nothing here — verified by
re-reading the file from disk, not from memory.** v0.9's move is `scaffolds.mod1`–`mod4`
onto theme 600 `#0678a6`; this figure draws no scaffold and never reads
`palette.SCAFFOLDS` (`grep -n SCAFFOLD make_figure.py` is empty). It is a single-route
figure, so every panel is one series in `MAIN` with the house density ramp, all read
through `palette.py`, which parses the JSON at import — there is no cached copy to go
stale.

## Open

- The caption's figure number is a placeholder. Supplementary numbering is not settled,
  so it appears in no file or directory name here.

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
