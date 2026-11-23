# pooled_readout

*Built under its slug. Supplementary numbers are assigned in one pass once the main
text stops moving, so no S-number appears in any file name.*

How deeply the 12,000-member library was actually read, and how hard the 32 candidates
depend on where the read floor was set. Three bare panels.

```bash
~/.venvs/binder_map/bin/python3 code/supp/pooled_readout/make_figure.py
```

Emits `output/pooled_readout/panel_{a,b,c}_*.{pdf,png,svg}`. Takes about a minute: it
reads eight ~1 MB count tables and a 9 MB manifest off the NAS on every run.

### The panel widths ARE the page layout

> "pull B next to A and C can be below them. That's fine. That figure is good."
> — Leo, 2026-11-19

`code/supp/make_supp_book.py` packs panels greedily into rows no wider than its text
block — `A4_W` 210 mm less 2 × `MARGIN` 14 = **182 mm** (that file's own comment says
"178 mm double-column"; it is wrong, derive it) — in **filename order**, and draws the
panel letters itself at assembly. So
"(a) and (b) side by side, (c) spanning underneath" is not expressed anywhere in this
script except in the **widths**, which is why they are named constants (`W_A`, `W_B`,
`W_C`) rather than literals in three `plt.subplots` calls.

The width the book measures is the **emitted** one, not the authored `figsize`:
`style.save` crops to the ink box and pads it, which *adds* a few mm to a panel whose
labels overhang and *removes* a lot from one with generous fractional side margins.

| panel | authored | emitted | row |
|---|---|---|---|
| (a) library recovery | 80.5 × 48 mm | **84.4 × 51.0 mm** | top, left |
| (b) reads per member | 78.0 × 52 mm | **85.2 × 55.7 mm** | top, right |
| (c) candidates vs floor | 188.0 × 48 mm | **175.9 × 51.1 mm** | spans below |

84.4 + 6 mm gutter + 85.2 = 175.6 mm, inside the block, so (a) and (b) share a row;
(c) at 175.9 mm cannot join them and fills the row underneath, matching the pair's width.
Both rows clear the 182 mm block by ~6 mm, so a small `MARGIN` change cannot reflow them. Nothing about the
panel contents changed — only these three numbers. **If you edit a `subplots_adjust`
margin here, re-measure: the emitted width moves with it and the page can silently
reflow to three rows again.**

## ⚠️ This figure was listed as blocked on screen data. It was not blocked.

`SM_BUILD_STATUS_2026-11-19.md` has `pooled_readout` as `PENDING` on screen data for
panels (a), (b) and (d), with "(c) exists". Both halves of that are wrong:

* **The pooled ONT data has been on the NAS since 2026-10-19**, at
  `LV_BinderDesign/LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment/counts__*.csv` —
  eight count tables, 12,000 rows each. Panels (a) and (b) are computed straight
  off them.
* **Panel (c) did not exist as a plot.** What exists is the gated hit list at one
  floor (`12_analysis/gated_Q1_75Reads_OrderedHits.xlsx`) and an interactive
  gating widget inside the 48 MB `splicing_report.html`. Neither is a curve of
  candidate count against read depth. It is built here.

Only panel (d) is genuinely blocked, and on wet-lab numbers rather than screen data
(see **Panel (d)** below).

## Sources, with mtimes

All primary, all read in place. `inputs/screen_points.csv` is **deliberately not
used**: it is reconstructed out of `splicing_report.html` and holds only the members
with at least one POS read (10,999 in context A, 10,515 in B), so it cannot answer
panel (b), which is about the members that were never seen.

| file | mtime | what is taken |
|---|---|---|
| `…/11_alignment/counts__{POS,NEG}_{A,B}_rep{1,2}.csv` | 2026-10-19 10:50 – 16:37 | `intact_count`, `removed_count`, `total` per library member, 8 × 12,000 rows |
| `…/11_alignment/reference_manifest.csv` | 2026-09-30 16:58 | `reference_id` → `binder_design` per context; 48,000 rows = 12,000 × 4 contexts |
| `…/12_analysis/gated_Q1_75Reads_OrderedHits.xlsx` | 2026-10-19 17:47 | the 32 rows the gate reproduction is checked against |
| `…/12_analysis/report_render_provenance.json` | 2026-10-19 16:57 | the screen's own default gate, and the pipeline/renderer commits |

The counts are the **combined run 1 + run 2** tables (`scope` in the provenance JSON:
*"Combined Run 1 + Run 2; presentation-only regeneration from validated count
tables"*). `rep1` and `rep2` are **sequencing replicates of one library**, summed as
the screen's own report sums them, not biological replicates.

## The four samples

Each context is a backbone pair, from `reference_manifest.csv`:

| sample | role | backbone | assigned reads |
|---|---|---|---|
| POS_A | positive | `pLV005.gb` | 15,408,861 |
| NEG_A | negative | `pLV008.gb` | 32,890,281 |
| POS_B | positive | `pLV006.gb` | 40,821,823 |
| NEG_B | negative | `pLV007.gb` | 35,850,445 |

Every context carries its own 12,000 references and `context_equivalence_size` is 1
on all 48,000 rows, so inside a context no two library members are indistinguishable
and the reference → design mapping is one-to-one. The script asserts this.

## The gate is reproduced, not trusted

`reads ≥ 75` in POS **and** in the matched NEG, `splicing ratio ≥ 0.07`,
`fold change ≥ 10`, where fold change is POS splicing over matched-NEG splicing with
the report's own pseudocount (NEG splicing = 1/NEG reads when no NEG read was
spliced). Off the count tables this returns **7 (context A) + 25 (context B) = 32**,
which is `gated_Q1_75Reads_OrderedHits.xlsx` exactly. `load()` asserts it, so every
number below fails loudly if a source changes.

For contrast, the screen's own *default* gate (`report_render_provenance.json`:
`min_reads` 100, `min_splicing` 0.15, `fc_threshold` 5, `fc_cap` 50) returns 6. The
32 are a **selection on ranked splicing ratio and fold change, not a hit rate**.

## Panels

* **(a)** Library recovery per sample, out of 12,000: filled = seen at all, open =
  at or above the floor. The count is labelled on each bar; the denominator is the
  same 12,000 for all eight bars and is drawn once as the reference line.
* **(b)** Reads per library member, each sample ranked independently, log y. The
  75-read floor is drawn across the axis. Members with zero reads cannot sit on a log
  axis, so each curve simply ends at that sample's recovered count — the gap to
  12,000 on the x axis is the unrecovered library, quantified in (a).
* **(c)** Candidate count against the read floor, 1 to 500 reads, holding the
  splicing and fold-change cuts fixed. Both contexts and each separately.

### (c): the area under each curve is filled in that curve's own colour

> "In c: Can you color the volume under the curve the color of the hits that it
> represents? Otherwise good" — Leo, 2026-11-20

**How the overlap is handled, which had to be decided rather than left to the
drawing order.** `total = A + B` exactly, so the total area *contains* both context
areas, and below a floor of about 120 the two context areas also overlap each other,
because both contexts still hold candidates there. Three options, and what was chosen:

| | |
|---|---|
| **stacked** | would partition the total exactly — A's share, then B's — but the upper boundary of the B band would not be the drawn B line, so the band labelled B would not be the curve labelled B. **Rejected.** |
| **separately translucent, fixed order** | each curve fills its own area, α 0.10, painted **total → B → A** and back-to-front by construction, never by value. **Chosen.** |
| one fill only | loses exactly what Leo asked for — the per-context colour. **Rejected.** |

So the panel reads as three nested backdrops, each under its own line: the palest
region is the total's headroom over both contexts, the middle one is where context B
contributes, and the darkest wedge at the left is where **both** contexts still hold
candidates. Two layers reach α 0.19, which is why the 75-read rule and the annotated
point are lifted above every fill by `zorder` (5 and 6) rather than merely drawn
afterwards — a later fill would have washed the floor out, and the floor is the one
mark the panel is read against. The rule was also darkened from `NEUTRAL["400"]` to
`NEUTRAL["600"]` so it still reads over the darkest wedge.

Nothing else about the panel moved: the sweep, the gate, the key and the annotated
point are as approved.

Sample encoding: a sequencing sample is not a design route and takes no route colour
(COLORS.md rule 2). Context carries the theme step (A = 900, B = 500), role carries
the line style (POS solid, NEG dashed), and the tick or key label names both.

**Palette v0.9 (2026-11-20) does not touch this figure.** The only change is
`scaffolds.mod1`–`mod4` moving to theme 600 `#0678a6`; this figure uses `THEME`
steps 900/500, `MAIN`, `NEUTRAL` and `SURFACE`, and never reads `palette.SCAFFOLDS`.
Re-rendered against v0.9: (a) and (b) are unchanged, down to their emitted sizes
(84.4 and 85.2 mm, as in the table above), and (c) differs only by its new fills —
it still emits at **175.9 mm**, so the page does not reflow.

## Numbers the caption carries (all printed on every run)

**Recovery, out of 12,000 ordered members**

| sample | seen at all | ≥ 75 reads |
|---|---|---|
| POS_A | 10,999/12,000 (91.7%) | 8,185/12,000 (68.2%) |
| NEG_A | 11,793/12,000 (98.3%) | 11,111/12,000 (92.6%) |
| POS_B | 10,515/12,000 (87.6%) | 9,272/12,000 (77.3%) |
| NEG_B | 11,137/12,000 (92.8%) | 9,255/12,000 (77.1%) |

124,971,410 reads assigned across the four samples. **The POS samples are the
limiting ones**: a design is only testable where both its POS and its matched NEG
clear the floor, and POS_A clears it for just over two thirds of the library.

**🚨 `CHECK:` discharged — every scale figure here is computed at 12,000 members.**
Nothing is inherited from the lab protocol's 24,682-member worked example; no quantity
in this folder was derived from it. At 12,000 members a 100-fold read target is
**1,200,000 assigned reads per sample**, which all four samples exceed by 13× to 34×.
That is read depth only, not library coverage through the cloning and transfection
bottlenecks — see panel (d).

**Candidates against the read floor**

| floor | context A | context B | total |
|---|---|---|---|
| 1 | 29 | 28 | 57 |
| 50 | 12 | 27 | 39 |
| **75** | **7** | **25** | **32** |
| 100 | 2 | 22 | 24 |
| 200 | 0 | 15 | 15 |
| 300 | 0 | 5 | 5 |
| 500 | 0 | 0 | 0 |

## ⚠️ What panel (c) actually shows, which is not what the build plan hoped

The build plan calls (c) *"a genuine robustness result… whether the 32 are a stable
call or an artefact of where sequencing stopped"*. **The files answer: the count is
not stable, but the set is nested.** Two separate facts, and the caption needs both
or it misleads:

1. **The count falls monotonically with the floor**, 57 at 1 read to 0 at 500, with
   no plateau anywhere. At 50 it would have been 39, at 100 it would have been 24.
   32 is the value at 75 and nothing distinguishes 75 from its neighbours in the
   data — the floor is a judgement call and the candidate count follows it directly.
   The script asserts monotonicity.
2. **No floor above 75 ever admits a design the 32 does not already contain.** The
   script asserts this for 100, 150, 200 and 300. Raising the floor only removes:
   24, 19, 15 and 5 of the 32 survive those floors. So the 32 is a strict superset
   of every stricter call, and the instability is entirely on the permissive side —
   a lower floor would add low-read designs, not reshuffle the ranking.
3. **The two contexts are not equally supported.** Context A collapses from 29
   candidates at floor 1 to 7 at 75 and 0 by 200; context B holds 15 at 200. POS_B was
   sequenced **2.65× deeper** than POS_A — 40,821,823 against 15,408,861 assigned
   reads — so the A/B split in the ordered set (7 against 25) is partly a **depth
   artefact**, not biology. **This belongs in the caption**: any per-context
   comparison downstream inherits it.

## Panel (d) — optional, and genuinely blocked

`NEED:` **the cloning and transfection bottleneck counts.** Coverage against a
100-fold target needs the number of library copies surviving each step —
transformation colonies, cells transfected, cells sorted or harvested. The
sequencing-depth bottleneck is measurable here (above); the others are not on disk.
The colony counts are already a recorded `NEED:` in `code/supp/library_pcr_cycle_gel/`
(journal LV1.40 records the intervention and its reason but no number; LV2.3 blames
the pLV008 failure on a low count without giving it), and the pooled screen is
LV2.19, which is **past the end of the Notion journal page** *LV Wet-Lab Journal
(LV01.01 – LV2.18)* — so there is no journal entry for this run to read them from.

**Who supplies it:** Leo, from the LV1.40/LV2.3 plate records and the LV2.19
transfection sheet; or Matteo Lampis, who ran the alignment
(`gobi_work_directory: /home/mlampis/ont_matched_reports/LV2_19_combined_20261019`)
and may hold the cell counts for that run. Panel (d) is marked optional in the build
plan and is **not drawn**; drawing one measurable bottleneck out of four would read
as a coverage audit that was never done.

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the n. 191 words.*

**Pooled readout depth, and the read floor behind the 32.** Two thirds of the library
is testable in the shallower context, and the size of the candidate set follows the
read floor directly. (**a**) Members recovered out of the 12,000 ordered,
filled for any read and open for at or above the floor; the positive samples are
limiting, at 8,185/12,000 (68.2%) in context A and 9,272/12,000 (77.3%) in context B.
(**b**) Reads per library member in each sample, ranked independently; the dashed line
is the 75-read floor applied to both the positive sample and its matched negative.
Curves end where a sample's recovered members end. (**c**) Candidates against the read
floor at fixed splicing (≥ 0.07) and fold-change (≥ 10) cuts. The count falls
monotonically, from 57 at one read to 32 at the applied floor of 75 and 0 by 500, with
no plateau; no floor above 75 admits a design the 32 does not contain, so the ordered
set is a strict superset of every stricter call. Context B is sequenced 2.65× deeper
than context A (40,821,823 against 15,408,861 assigned reads), so the 7/25 context
split of the ordered set is partly a depth artefact. 12,000 library members, four
samples, two sequencing replicates each, one biological replicate.

## Open

* `NEED:` **a citing sentence in Section 4**, at the ONT readout. The build plan
  assigns one but none is written.
* The A/B depth imbalance in (c) — POS_B 2.65× deeper than POS_A — argues the ordered
  set's 7/25 context split is partly a sequencing artefact. Whether that changes anything downstream — `bif_quasicontrol`
  and `arrayed_flow` both split by context — is **not** settled here.
* `inputs/screen_points.csv` is now reproducible from primary files rather than from
  the HTML report. If `12_analysis/matched_comparisons.csv` ever becomes readable
  (it is chmod 700), regenerate it from there and check all three against each other.
