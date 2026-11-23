# negative_controls — negative control design and construction

**Six** bare panels in `output/negative_controls/`, all **authored 85 mm wide**
and all **emitted at 84.8 mm** (a 27 mm tall, b 30, c 26, d 52, e 46, f 45;
matplotlib rounds the SVG canvas to whole pixels at the save dpi, so 85 mm comes
out as 84.8 — the point is that the six agree, not the exact figure).
`style.save()` writes PDF, PNG and SVG; assembly should consume the SVG, whose
text is live. Panels are placed in
**file-name order** by `make_supp_book.py`, which is why the panel letter lives
in the file name and why the old `panel_d_failed_scaffolds` had to be
**renamed**, not just reordered — see the 2026-11-20 revisions below.

Every panel goes out through `emit()`, which saves with `tight=False`, so the
**emitted** box IS the authored 85 mm canvas. That is what puts the page in **three rows of
two** — see the layout revision immediately below. Content runs
**189 mm** of the 283 mm page body, which clears the 299 mm overflow this
README used to record as pre-existing.

```bash
~/.venvs/binder_map/bin/python3 code/supp/negative_controls/make_inputs.py
~/.venvs/binder_map/bin/python3 code/supp/negative_controls/make_figure.py
```

`make_inputs.py --skip-fold` leaves out the ViennaRNA step, which runs on
Euler because ViennaRNA is not in `~/.venvs/binder_map`.

`make_figure.py` runs `_selfcheck()` before it draws anything: the drawn fold
has to be the fold the table reports (pair count, length), the stems have to
carry the chemistry claimed below, the naview layout has to come out stem-
horizontal, the replacement's best Watson–Crick stem has to be the 3 bp `AUU`,
pART256's insert has to still fold over 100 pairs, panel (d)'s shortfall has to
close (`_check_missing_repeat`: two pairs short a repeat and by how much, each
shortfall equal to the units lost plus 9 nt per hairpin position kept, and equal
to the array-span difference), and the two failed assemblies have to still be
the LONGEST at their copy number, since panel (f) prints the word "longer"
beside a measured number. Corrupting one base in the
layout CSV aborts the build rather than quietly drawing a different molecule —
checked, not assumed.

Three are new on 2026-11-20. `_check_missing_repeat` no longer expects
**one** pair short a repeat but **two, by one unit and by four** — pART257 and
pLV16 — and it re-closes pLV16's as-built shortfall against its own ordered map
(392 nt = 4 × 98) on top of the shortfall against its positive (536 nt). And
`check_emitted()` re-derives the **page packing** from the SVGs the run just
wrote, using `make_supp_book`'s own `rows_for` and `svg_size`: one emitted width
across the six, a pair plus the gutter inside the text block, every row exactly
two, content inside the page body. And `emit()`, which every panel is saved
through, compares each figure's tight bounding box with its canvas **before**
saving, because `tight=False` CLIPS what `tight` padded around and a clipped
label is silent — it just loses its last word.

## Revision 2026-11-20 (Leo, later) — two per row, palette v0.9, and the arrayed screen's negative

> "Align the panels. Should be two per horizontal level always."
> "for d add the negative with less copies in the 1 by 1 resulting in pLV15"

### 1. The page is three rows of two, and the 299 mm overflow is gone

**The panels were never 85 mm on disk.** They were authored at 85 and emitted at
**89.9 to 92.0**, because `style.save()` defaults to `bbox_inches="tight"`, which
adds `pad_inches=0.1` — 2.54 mm each side — and expands around any ink that left
the canvas. `make_supp_book.py` packs greedily into the 182 mm text block with a
6 mm gutter:

| | two per row needs | fits 182 mm? |
|---|--:|---|
| authored 85 mm | 85 + 6 + 85 = 176 | yes |
| as emitted ~90 mm | 90 + 6 + 90 = 186 | **no** |

So all six landed **one per row** and the page ran **299 mm** of content on a
283 mm body. That warning was not pre-existing in the sense of being
unavoidable — it was this.

The fix is the house idiom five sibling figures already use (`amplicon_design`,
`metric_panel`, `repredictor_comparison`, `library_pcr_cycle_gel`,
`missing_gel_records`): emit the authored canvas with **`tight=False`**. Every
panel here emits exactly **85.0 mm** now, so the pairing is deterministic rather
than a function of whichever label happens to be longest.

`tight=False` **clips** ink outside the canvas where `tight` padded around it, so
the four places that overflowed were brought inside. All four were **measured**
with `fig.get_tightbbox` against the canvas, not guessed:

| panel | overflowed by | fix |
|---|---|---|
| a | left edge, **2.11 mm** — the right-anchored `pART256, negative` label | `x0` is now `max(w_mm(label)) + 2.3`, not the round `16.0` |
| c | bottom, **0.22 mm** — the axis caption | canvas 25 → **26 mm** |
| f | bottom, **2.52 mm** — the caption row | canvas 42 → **45 mm** |
| e | right, **1.84 mm** — the axes' own right spine | `subplots_adjust(right=0.99 → 0.955)` |

**Two checks close it off.** `emit()` wraps every save: it draws the figure,
compares `fig.get_tightbbox` with the canvas, and **fails the build naming the
side and the millimetres** if anything hangs off, because `tight=False` clips
silently and the layout fix would otherwise have traded a visible page overflow
for an invisible missing label. It caught panel (e) still **0.03 mm** outside at
`right=0.96`, which is why that is now `0.955` — a hairline, and still a clip.

**`check_emitted()`** then runs at the end of the build, reads the
SVGs the run just wrote, and re-derives the packing with
**`make_supp_book.rows_for` and `make_supp_book.svg_size` — the book's own two
functions, imported, not reimplemented.** It asserts the six emit **one**
width — not `== 85.0`, which fails at the 84.8 mm matplotlib actually writes —
that a pair plus the 6 mm gutter clears the text block (**175.6 of 182 mm**),
that every row holds exactly two, and that the content fits the body. It prints
the rows. One longer tick label, or a `tight=True` slipping back in, now stops
the build instead of silently reflowing the page.

⚠️ **Rows pair by file name, not by height**, so (c) at 26 mm sits beside (d) at
68 mm. That is the cost of "two per horizontal level always" on a figure whose
panels are not the same height, and geometry cannot fix it: the book places
panels in file-name order and draws the letters from that order, so balancing
the rows would move the letters.

### 2. Palette v0.9: the four architectures move from theme 700 to 600

`scaffolds.mod1`–`mod4` are all **`#0678a6`, theme step 600**, where v0.8 had
them all at `palette.MAIN` `#015d82`, step 700. Lighter, on Leo's instruction.
Still **one colour for all four**. `palette.F30` `#bfe0f5` (step 200) is
unchanged, and so is `DIVERGING[-1]` `#ae4e1a` for the Csy4 site.

Nothing in `make_figure.py` names the colour — it reads
`palette.SCAFFOLDS[arch]` — so the change landed as a re-render. The only code
that moved was the comment, plus `FAMILY`, which was a hardcoded
`{"pART252": MOD4, "pART256": NEUT, …}` and is now `family(design)`, built from
the table, because panel (d) gained two constructs the dict did not hold.

**Remeasured, with the repo's own converter** (`common/style._lab`, plain
CIELAB ΔE*ab):

| pair | v0.8 | v0.9 |
|---|--:|--:|
| scaffold blue vs `palette.F30` | 52.4 | **44.5** |
| scaffold blue vs `NEUTRAL["300"]`, the negatives | 54.5 | **48.8** |
| scaffold blue vs the Csy4 site `#ae4e1a` | 95.7 | **92.1** |

Both separations shrink and both stay large. F30 still reads as a different kind
of object from the four designed architectures.

⚠️ **This README's old "ΔE 54.6 from MAIN" was against the wrong pair.** 54.5 is
`MAIN` ↔ `NEUTRAL["300"]`; `MAIN` ↔ `palette.F30` is **52.4**. The number was
very nearly right and attached to the wrong comparison, which is why it is
replaced with a measured table rather than carried forward.

⚠️ **One contrast is weak, and v0.9 did not cause it.** `palette.F30` `#bfe0f5`
and `NEUTRAL["300"]` `#cfd5d9` are **ΔE 12.5** apart, and panel (d) stacks an F30
positive bar directly above its grey negative in two of its three pairs. That
pair reads by **label and by length** — which is the rule this whole figure runs
on since v0.8 — and the hairline unit divisions inside the F30 bars carry the
count. It is recorded here so nobody later asks those two fills to carry more.

### 3. Panel (d) gains a third pair: the arrayed screen's negative

> Leo: "for d add the negative with less copies in the 1 by 1 resulting in pLV15"

⚠️ **It is pLV16, not pLV15 — and pLV15 is a real construct, which is why the
name had to be settled before anything was drawn.** Put to Leo, who answered
"ah yeah pLV16". Both maps are on disk and they are different plasmids:

| | file | bp | Csy4 sites | replacement loops | F30 units | from |
|---|---|--:|--:|--:|--:|---|
| pLV13, posA | `pLV013/pLV013.gb` | 6,479 | **16** | 0 | **8** | pART158 |
| pLV14, posB | `pLV014/pLV014.gb` | 5,951 | **8** | 0 | 0 | pART252 |
| **pLV15** | `pLV015/pLV015.gb` | 5,879 | 0 | **8** | **0** | pART256 |
| pLV16, NegA, as **ordered** | `pLV016/pLV16.gb` | 6,335 | 0 | **16** | **8** | pART257 |
| **pLV16, NegA, as built** | `pLV016/pLV16_asbuilt_LV2-26-004c_annotated.gb` | **5,943** | 0 | **8** | **4** | pART257 |

**pLV15 is posB's matched negative and was never transfected.** Counted here from
its own sequence: 5,951 − 5,879 = **72 = 8 × 9**, the same per-position
arithmetic as pART252 → pART256, which is what identifies it — the plain tandem
scaffold with its eight hairpins knocked out, no F30 unit at all. It appears
**0 times** in the LV2.33 transfection plan and 0 times in the as-performed
record, where `pLV16` appears **41** times, and
`inputs/lv2_34_1_flow_on_off.csv` carries four conditions — `PosA`, `PosB`,
`NegA`, `FullNeg` — with no column for a second scaffold negative. So the
one-by-one screen used **pLV16 as the scaffold negative for both positives**,
which is a design decision visible in the plan, not a bench slip. (Two other
files in the repo call it an accident; see `Open` at the end.)

**The bar drawn is pLV13 → pLV16 as built**, because that is the pair the arrayed
screen ran, and the as-built sequence is what was in the well.

#### The closure, and it closes twice like the other two

Measured, not transcribed, under the same two-sources rule as the rest of the
panel. The as-built map is annotated by **pLannotate**, which numbers every
feature (`NegativeControlLoop 3 of 8 (as-built)`, `F30 unit 3 of 4 (as-built)`)
where SnapGene repeats one label, so `_labels()` now counts a **prefix** rather
than an exact string: 8 loop features and 4 unit features, against 8 and 4 from
the sequence.

| | nt | hairpin-array span | period |
|---|--:|--:|--:|
| pLV13 | 6,479 | **888** (16 Csy4 sites) | 116 |
| pLV16 as built | 5,943 | **352** (8 loops, 4 units) | 98 |

* **6,479 − 5,943 = 536**, and **888 − 352 = 536** — the whole difference is
  inside the hairpin array, as in the other two pairs.
* **536 = 4 × 116 + 8 × 9** — four whole F30 units at pLV13's own period, plus
  9 nt at each of the eight positions that remain. This is what the panel labels
  and draws.
* **536 = 4 × 98 + 16 × 9** — four of the *negative's* own units plus 9 nt at all
  sixteen positions the design had. Equally true; the same pair of readings this
  README records for pART257's 242, and mixing them is the error that gave 224
  there.
* It also splits by **stage**: 6,479 − 6,335 = **144 = 16 × 9** (the ordered map
  was valency-matched to pLV13 — eight units, sixteen positions, shorter only by
  the replacement), and 6,335 − 5,943 = **392 = 4 × 98** (the as-built deletion).
  144 + 392 = 536.

All of those are `assert`ed in `make_inputs.design()`, and the first two plus the
392 are re-asserted against the CSV in `_check_missing_repeat`.

#### The two F30 shortfalls are different events, and the panel says which

The one thing here that could mislead, so it is labelled on the panel rather than
left to the caption:

| | shortfall | what happened |
|---|---|---|
| pART257 | 1 unit of 8 | **never placed.** pART257 *is* the seven-unit clone the spike-in screen adopted, so its own map is its as-built record. Panel (e) is that screen: **0/56 (0%)** colonies at eight units. |
| pLV16 | 4 units of 8 | **absent from the sequenced plasmid.** Ordered at eight units and sixteen loops; whole-plasmid ONT gave four and eight. |

`negative_design_units` is what separates them — the negative's own **ordered**
unit count, which for pART256 and pART257 equals their as-built count and for
pLV16 does not. The panel reads its label off that column, and `_selfcheck`
asserts **pLV16 is the only row** permitted to print an "ordered at …" line.

**No cause is claimed for either.** One array was built by spike-in and came back
one unit short; the other was ordered whole and came back half. The panel reports
each shortfall's size and the construct it is measured against, and nothing else.

⚠️ **Two numbers in the as-built file's own `COMMENT` are wrong and are not
used.** It gives the design as "6336 bp" and the deletion as "393 bp"; the
ordered map measures **6,335** and the deletion is **392** — which the comment's
own later line, `4 x 98 = 392`, agrees with. Every number here is read from the
sequence, so the file's prose being off by one never reaches the panel. Worth
correcting in the file.

### 4. The prose came off panel (d) — it was the caption's job

> Leo: "dont add to much text. thats for the caption to do."

Six explanatory sentences had accumulated under the bars. Four were moved to the
caption draft below, nothing dropped:

| moved to the caption | why it is caption material |
|---|---|
| `72 nt shorter than pART252: 9 nt at each of 8 positions, no unit lost` | arithmetic |
| `242 nt shorter than pART158: 1 × 116 nt F30 unit short, plus 9 nt at each of 14 positions` | arithmetic |
| `536 nt shorter than pLV13: 4 × 116 nt F30 units short, plus 9 nt at each of 8 positions` | arithmetic |
| `ordered at 8 units and 16 loops; sequenced at 4 and 8` | provenance |

The test applied: **a label earns the canvas only if the drawing cannot say it.**
None of those four is needed to read the bars.

**What stayed, and why.** The row names, the value labels (`888 nt   8 units,
16 Csy4 sites`), the pair names — and a **two-or-three-word marker on each
unfilled group**, because the distinction between the two shortfalls is the
thing that was asked to be *visualised* ("Exactly, this was known. Thats why i
want it visualized") and a caption cannot do it: the reader has to tell the two
unfilled groups apart while looking at them.

| | marker | was |
|---|---|---|
| pART257 | **`never placed`** | `one F30 unit, 116 nt, never placed` |
| pLV16 | **`designed, not built`** | `4 F30 units, 116 nt each, absent from the sequenced plasmid` |

Three words against the previous eight, and the contrast survives: one unit was
never built, four were designed in and are missing from the plasmid. The unit
count and the 116 nt are **not** repeated in the marker — the blocks are drawn
one per unit and the row's own value label gives the nt, so the old marker said
three times over what the drawing says once. If `designed, not built` ever reads
as ambiguous beside `never placed`, the longer `not in the sequenced plasmid` is
a one-string change and the caption carries the full provenance either way.

**The canvas shrank rather than the bars spreading**: 68 → **52 mm**. Bar
height, the mm-per-nt scale and the row spacing are untouched; the canvas lost
the band the text occupied. The lowest ink now sits at 50.8 mm, measured.

### 5. The same test, applied to (a) and (e)

The instruction is a principle, not a note about one panel, so (a) and (e) went
through it too. **(b), (c) and (f) were already clean** — every label there is a
value measured off the mark it sits beside.

**Panel (a), 33 → 27 mm.**

| line | verdict |
|---|---|
| `pART252 6,194 nt, pART256 6,122 nt, 72 nt apart` | **gone** — pure arithmetic, and the caption already carried it |
| `8 positions replaced, 0 Csy4 sites left` | **cut to `the same at all 8 positions`** |
| `46% GC across the 28 nt replaced, against 26%` | **kept, shortened to `28 nt replaced: 46% GC, against 26%`** |
| `shared flank, mod4 scaffold` | **kept** — a row name, and under v0.9 the colour cannot say it |

⚠️ **`8 positions replaced` did not go wholesale, and the reason matters.** The
test is whether the *drawing* can say it, and **this panel draws ONE position**.
A reader cannot count eight off the marks, so without that line the swap reads
as happening at a single site — which is the one misreading (a) exists to
prevent. What went is the `0 Csy4 sites left` half, which *is* arithmetic and is
in the caption. The scope stayed, at six words instead of eleven.

**The GC line stayed and now names its window first.** `28 nt replaced: 46% GC,
against 26%` — 32 mm against 41. It is the one place the figure says *which*
window it means, the 28 nt replaced and not the 16 nt hairpin core, and both
figures are right under their own rule (see "Which GC window the figure means").
A caption can say it; on the panel it stops the misreading where it would
happen.

**Panel (e): `4 clones at seven units; one of them was adopted` →
`one became pART257`.** The first half was the bar's own value label said
twice — the bar is already annotated with its count — so only the provenance is
left, and it now names the construct instead of saying "one of them". That
clause stays rather than going to the caption because removing it leaves the
panel ambiguous about what happened to those clones: without it (e) is a
histogram, and nothing on the page says which clone became the negative that
panel (d) draws. Shorter *and* more informative than what it replaced, which is
the right trade; the panel's height is unchanged, the annotation being inside
the axes.

#### What the panel draws, now

**Six bars in three pairs** at **one** mm-per-nt scale, in ascending shortfall —
mod4 (none), pART158/pART257 (one unit), pLV13/pLV16 (four) — each pair under its
own name from the `pair` column. That name is not decoration: two of the three
pairs are F30 arrays and under palette v0.9 all four architectures are one
colour, so without it pLV13/pLV16 read as two more library pARTs instead of the
arrayed screen's constructs.

* Missing units are **one unfilled block per unit**, not one long gap, so four
  short reads as a count. Absence of fill is absence; no hatch; it survives a
  greyscale print.
* The blocks are drawn at the **positive's** 116 nt period, so they come out
  wider than the negatives' own 98 nt blocks — the decomposition the panel
  labels. pLV16's four close 464 of its 536, and the 72 nt left over is the
  8 × 9.
* The unfilled group's marker is anchored **left at the blocks, or right at the
  canvas edge** where the string would otherwise run off it — chosen on a
  `style.text_w` measurement, since `tight=False` clips rather than pads. Its
  leader sits at the **first missing block's left edge**, not the group's
  centre, where it came down in the middle of the string and read as a stray
  mark.
* Each pair's shortfall arithmetic, and pLV16's ordered-against-sequenced line,
  are **in the caption, not on the panel** — see the sub-revision above.

## Revision 2026-11-20 (Leo) — the right failure, and palette v0.8

> "with the failed scaffolds, I meant the negatives that failed and had missing
> repeats because of the spike in fail."

**Panel (d) was drawing the wrong failure.** It showed pART250 and pART251, the
two mod3 scaffolds whose assembly was never obtained. That result is correct and
stays — see **panel (f)** below — but it is not what the note asked for. The
failure Leo meant is the **spike-in** one: the negative that came back a whole
F30 repeat short.

| | panel | |
|---|---|---|
| **d** | `panel_d_missing_repeat` | **new** — each negative against its own positive, and the F30 unit the spike-in never placed |
| **e** | `panel_e_colony_screen` | unchanged; it is the screen that explains (d) |
| **f** | `panel_f_failed_scaffolds` | the old panel (d), **renamed**, kept whole |

**The mod3 result is kept as its own panel (f), not folded into the caption.**
It is measured, it is checked by `_check_failed`, and it is a *different failure
mode* — a construct that was never assembled, against a construct that was
assembled one unit short. Putting them in one panel would say they are the same
kind of event. (f) is last so it cannot be read as part of the spike-in story.

### Palette v0.8: colour no longer names an architecture

⚠️ **Superseded by v0.9, above: the one colour is now theme 600 `#0678a6`, not
`palette.MAIN` `#015d82`.** Everything else in this subsection still holds — one
colour for all four, so every panel names its architecture.

`scaffolds.mod1`–`mod4` are all `palette.MAIN` `#015d82` at v0.8 (COLORS.md
rule 9, superseded 2026-11-20); the indigo ramp is withdrawn. `palette.F30` `#bfe0f5` is
unchanged. The substantive consequence for this figure is that **the mod4 ring
and the F30 array can no longer be told apart by hue alone**, so every panel that
draws a construct now names its architecture:

| panel | before | now |
|---|---|---|
| a | "shared flank" (colour said mod4) | "shared flank, **mod4** scaffold" |
| c | `pART252` | `pART252  **mod4**`, `pART158  **F30**` on every bar |
| d | — (new) | `pART252  **mod4**` … `pART257  **F30**` on every bar |
| f | `pART250  mod3 2×` | unchanged; it already named every row |

The architecture is carried as **data** — a new `arch` column in
`supp_negative_controls_design.csv`, from `PAIRS` in `make_inputs.py` — and read
by `arch_of()`, so a name a panel prints is not a string the panel invented.

Panel (f) loses nothing to the change: the distinction it encodes is *obtained
vs not obtained*, which was always fill against no-fill, never hue. All eight of
its bars are now one blue and the panel reads exactly as before.

## Revision 2026-11-19 (Leo)

> "Negatives is good rearrange the pannels (maybe also split b into two panels)
> also would be good to see two failed scaffolds and how much shorter they were"

* **Panel (b) is split**, into the element scale and the construct scale. It
  carried both, and they are the two different arguments this figure exists to
  keep apart. Splitting them also let two drawn headers go: "the replaced
  element on its own" and "the whole scaffold insert" were scale dividers that
  only existed because one frame held two scales, and as standalone panels they
  would be titles, which STYLE_NOTES bans.
* **Panels rearranged**, five now:

  | | | scale |
  |---|---|---|
  | a | the swap at one hairpin position | sequence |
  | b | the two elements, ViennaRNA MFE fold | element |
  | c | the whole insert, folding energy | construct |
  | d | the two scaffolds never obtained | construction, design level |
  | e | the colony screen | construction, clone level |

  **Why this order.** a–b–c is one ascending scale, which is the distinction
  the figure exists to protect. The claim is made at (b) — *the replacement
  cannot form **the** hairpin* — and bounded immediately at (c), where the
  negative insert still folds 101 pairs over 319 nt, so no reader carries "the
  negative is unstructured" past the panel that disproves it. Putting the
  construct folds anywhere later would let that error travel. d–e is then what
  was obtained, design level before clone level, ending on the construct that
  was actually used, which is the figure's conclusion.
* **Panel (d) is new**: the two failed scaffolds, below.
* **Palette v0.7**: the two construct families now read as different colours.

### Which colours, under palette v0.9 (COLORS rule 9, superseded)

This figure draws two families and they are different kinds of object, so they
no longer share a colour:

| | colour | |
|---|---|---|
| pART252 (mod4 ring), pLV13/pLV16 | `palette.SCAFFOLDS[arch]` = `THEME["600"]` `#0678a6` | one of the four designed architectures — **all four share this one colour since v0.8, at step 600 since v0.9** |
| pART158 (F30 tandem array) | `palette.F30` = `THEME["200"]` | **not** one of them; the Fig 1 construct-map blue |
| Csy4 site | `DIVERGING[-1]` `#ae4e1a` | rule 4, as in the Fig 1 maps |
| both negatives | `NEUTRAL["300"]` | a reference, not a series (rule 3) |

Blue scaffold, red sites, exactly as in the Fig 1 construct maps — which is
why panel (a)'s shared flank is mod4 blue rather than grey: the flank *is* the
mod4 scaffold, and the panel **says so in text**, since from v0.8 on the colour
no longer can.

⚠️ **Do not read architecture off colour anywhere in this figure.** All four
designed architectures are one blue. **Every mark carries its label**
(STYLE_NOTES §9): every bar in (c), (d) and (f) is named with its construct
number *and its architecture*, (d) adds its pair's name and (f) the copy number.
F30 stays apart by colour as well as by name (ΔE **44.5** from the scaffold blue
at v0.9, remeasured — see the v0.9 revision for why the old 54.6 was the wrong
pair), because it is not one of the four.

⚠️ **One tension, recorded rather than hidden.** Rule 1 reserves theme steps
100–300 for fills and backgrounds and 400–900 for data marks, and `palette.F30`
is step 200. Rule 9 nevertheless assigns step 200 to the F30 array by name, so
the more specific and newer rule is followed here, with the hairline `INK` edge
and the label on every bar that STYLE_NOTES §9 requires of a pale fill. Where
the mark is **not** a construct the tension does not arise and rule 1 wins:
panel (e)'s bars are a count of colonies, a single series, so they stay `MAIN`.

## Sources

| source | what it gives |
|---|---|
| `LV_Cloning/pLV_NNN/pLV00{1..4}/*.gb` | pART158, pART252, pART256, pART257 |
| `LV_Cloning/pLV_NNN/pLV01{3,4,5,6}/*.gb` | the arrayed screen's plasmids: posA pLV13, posB pLV14, pLV15, NegA pLV16 (ordered **and** as built) |
| `LV_Transfection/LV2.33_transfection_{plan,AS-PERFORMED_2026-11-04}.html` | which of them were actually transfected — pLV15 in neither |
| `LV_Cloning/Projects/F30_mod/Assembled_Final/*.dna` | the **12 designed maps**, the three never obtained included |
| `LV_Sequencing/LV1.16_SangerSeq/*.dna` | which of the 12 were assembled |
| `LV_Sequencing/LV1.24_NightSeq/Sequence Report {1,2}_2 [V4858402]/*.fasta` | plate **49761**, colonies 1 to 24 |
| `LV_Sequencing/LV1.30_NightSeq/Sequencing_Results/*.fasta` | plate **49780**, colonies 25 to 56 |
| Euler `pyfurnace_road/venv` (ViennaRNA 2.7.2) | MFE folds **and the naview layouts panel (b) draws** |
| journal LV1.22 to LV1.34 (Notion, LV Wet-Lab Journal) | the build narrative |

Derived: `inputs/supp_negative_controls_{design,scaffolds,fold,foldxy,colonies}.csv`.

**Plate 49780 is LV1.28's plate, and its reads are filed under LV1.30, not
LV1.28.** `LV1.28_NightSeq/` holds only the submission sheet.

**`make_inputs.py` no longer hardcodes `/Volumes/wicky`.** The share mounts
under either that root or `/Volumes/bsse/wicky` depending on how it was
connected, and the hardcoded form fails outright under the other one. It now
takes whichever exists. (Several other scripts in this repo still hardcode it.)

## Panel (d): the F30 repeats the arrays came back without

> Leo: "with the failed scaffolds, I meant the negatives that failed and had
> missing repeats because of the spike in fail."
> Leo, later the same day: "for d add the negative with less copies in the 1 by 1
> resulting in pLV15" — which is **pLV16**; see the v0.9 revision above for why,
> and for what pLV15 actually is.

**Counted from the `.gb` annotations, not transcribed.** Every number below was
read twice — once as a motif count in the sequence, once as a `/label=` feature
count in the SnapGene-exported GenBank file — and `make_inputs.design()` asserts
the two agree, so a count here is two independent records agreeing.

| construct | file | bp | Csy4 extended | NegativeControlLoop | F30 junction features | hairpin-array span |
|---|---|--:|--:|--:|--:|--:|
| pART158 | `pLV001/pLV001_part158.gb` | 6,722 | **16** | 0 | **8** | 888 nt |
| pART257 | `pLV004/pLV004_part257.gb` | 6,480 | 0 | **14** | **7** | 646 nt |
| pART252 | `pLV002/pLV002_part252.gb` | 6,194 | **8** | 0 | — | 341 nt |
| pART256 | `pLV003/pLV003_part256.gb` | 6,122 | 0 | **8** | — | 269 nt |
| pLV13 | `pLV013/pLV013.gb` | 6,479 | **16** | 0 | **8** | 888 nt |
| pLV16, ordered | `pLV016/pLV16.gb` | 6,335 | 0 | **16** | **8** | 744 nt |
| pLV16, **as built** | `pLV016/pLV16_asbuilt_…annotated.gb` | **5,943** | 0 | **8** | **4** | 352 nt |
| pLV15, never transfected | `pLV015/pLV015.gb` | 5,879 | 0 | **8** | — | 269 nt |

The last four rows are new on 2026-11-20. pLV13's array is the **same 888 nt
span** as pART158's and at the same 116 nt period, which is what makes the two
F30 pairs comparable at one scale. The as-built map is annotated by
**pLannotate**, which numbers each feature individually — `NegativeControlLoop
3 of 8 (as-built)`, `F30 unit 3 of 4 (as-built)` — so `_labels()` counts a
**prefix** rather than an exact string; its 8 and 4 agree with the 8 and 4 the
sequence gives. pLV15 is drawn nowhere and is tabulated only so the construct
Leo named is on the record.

`F30_Stem1`, `F30_Stem2` and `F30_Stem3` are annotated 8× in pART158 and 7× in
pART257 as well, and the per-unit linkers (`SpacerSequence`, `TAGA-overhang`)
8−1 = 7× and 7−1 = 6×, all consistent with eight units against seven.

**pART257 is one whole F30 unit short of pART158**: seven units and fourteen
replacement loops against eight and sixteen. That is a **spike-in** failure, not
an assembly failure — panel (e) is the screen it came from, where **0/56 (0%)
colonies carried eight units** and a seven-unit clone was adopted.

**pLV16 is four whole F30 units short of pLV13**: four units and eight loops
against eight and sixteen. That one is a different event — it was *ordered* at
eight and sixteen, and the sequenced plasmid is 392 nt shorter than its own
design, 4 × its own 98 nt unit. The panel labels the two differently and the
v0.9 revision above says why.

### The period is 116 nt, and the arithmetic closes twice

Measured from the hairpin-position gaps in each file, not assumed:

| | gaps between successive positions | period |
|---|---|--:|
| pART158, 16 Csy4 sites | 48, 68, 48, 68, … | **116 nt** |
| pART257, 14 replacements | 39, 59, 39, 59, … | **98 nt** = 116 − 2 × 9 |
| pLV13, 16 Csy4 sites | 48, 68, 48, 68, … | **116 nt** |
| pLV16 as built, 8 replacements | 39, 59, 39, 59, 39, 59, 39 | **98 nt** |
| pART252, 8 Csy4 sites | 43 to 47, irregular | **none** |
| pART256, 8 replacements | 34 to 38, irregular | **none** |
| pLV15, 8 replacements | 34 to 38, irregular | **none** |

The mod4 ring has **no repeat unit at all** — its eight positions sit at
irregular gaps — which is precisely why it could not lose one, and why it is the
right control to draw beside the F30 pair.

Each replacement is 9 nt shorter than the Csy4 repeat it stands in for
(19 against 28), so:

* **pART158 → pART257**: 6,722 − 6,480 = **242 nt** = one 116 nt F30 unit + 14 × 9
* **pART252 → pART256**: 6,194 − 6,122 = **72 nt** = 8 × 9, no unit lost
* **pLV13 → pLV16 as built**: 6,479 − 5,943 = **536 nt** = four 116 nt F30 units
  + 8 × 9 — and it splits by stage as 144 (= 16 × 9, the ordered design) + 392
  (= 4 × 98, the as-built deletion)
* **pLV14 → pLV15**: 5,951 − 5,879 = **72 nt** = 8 × 9, the arithmetic that
  identifies pLV15 as posB's matched negative. Not drawn — it was never
  transfected.

The 8-site pair fixes the 9 nt per position independently of the F30 pair.
Both decompositions of the 242 close (116 + 126, or the negative's own 98 nt
unit + 16 × 9 = 98 + 144); mixing them is what gave the wrong 224, below.

**The length difference lies entirely inside the hairpin array** — 888 − 646 =
242, 341 − 269 = 72 and 888 − 352 = 536, the same numbers as the whole-plasmid
differences — so the panel's one drawn quantity carries each comparison.
`make_inputs` asserts that for all three pairs.

### What the panel draws

**Six** bars at **one** mm-per-nt scale, in three pairs, each pair's positive
above its negative and under its own name from the `pair` column, each bar
labelled with its construct **and its architecture** (palette v0.9 — colour no
longer says which it is). The F30 bars are divided into their units by
hairlines, so eight blocks against seven, or against four, is readable as a
count and not only as a length, and the negatives' blocks come out visibly
narrower because their unit is 98 nt against 116.

The units the array came back without are **unfilled** blocks appended where
they would have sat, **one block per unit** — absence of fill is absence, the
same convention panel (f) uses, and it survives a greyscale print. A hairline
carries the **positive's** array end down across its negative, so each shortfall
is a bounded gap: pART257's one block closes 116 of the 242 and the 126 nt left
are the 14 × 9; pLV16's four close 464 of the 536 and the 72 nt left are the
8 × 9.

Three things to know when reading it:

* The **last block of each F30 bar is short.** The span is measured first
  position to last position, so it stops *inside* the final unit (888 = 7 × 116
  + 76). Drawing eight equal 116 nt blocks instead would break the arithmetic.
* The unfilled blocks are drawn at **116 nt, the positive's unit**, so they are
  wider than the negatives' own 98 nt blocks. That is the decomposition the
  panel labels (116 + 14 × 9, and 4 × 116 + 8 × 9); the other one (98 + 16 × 9,
  4 × 98 + 16 × 9) is equally true and the label would then have to count
  sixteen positions neither negative has.
* **Two of the three pairs are F30 arrays** and all four designed architectures
  are one colour, so the pair names are load-bearing: pLV13/pLV16 are the
  arrayed screen's plasmids, pART158/pART257 the library constructs, and nothing
  but the printed name distinguishes them.

**What the panel does not claim.** Nothing here says why the eighth unit was not
placed, or why half of pLV16's array is missing. One is a repeat array built by
spike-in, where 0/56 (0%) colonies came back at eight units; the other was
ordered whole and sequenced at half. The panel reports the size of each
shortfall and the construct it is measured against.

## Panel (f): the two failed scaffolds, and which way the length went

> Leo: "would be good to see two failed scaffolds and how much shorter they were"

**Which two, defined from the files rather than assumed.** Of the twelve
designed F30_mod combinations — four architectures at one, two and three
copies — the outcome of each is derivable from file names alone:

* `pART249_notcloned.dna` — the **fragment was never delivered**, so this one
  was never attempted.
* present among `LV_Sequencing/LV1.16_SangerSeq/*.dna` — **assembled** (9 of
  them, plus pART158, which is not one of the 12).
* neither — **the assembly was not obtained**: `pART250` and `pART251`.

So 12 designed, 11 delivered, **9/12 (75%) assembled**, and the two failed
scaffolds are **pART250 and pART251**, the mod3 two- and three-copy
constructs. They are the two that were attempted and failed; pART249 is a
different failure mode — a synthesis that never arrived — and is left off the
panel, which is also why panel (f) shows only 2× and 3×: at one copy there is
no assembly outcome to compare with.

**pART249's architecture is proved from sequence, not from the numbering.** Its
file is the bare `pART249_notcloned.dna` and carries no architecture in its
name, unlike the other eleven. `make_inputs.py` asserts that its Csy4 array is
contained in pART250's, which is what makes it the one-copy mod3 member.

### 🚨 They are not shorter. They are longer, and by how much is measured

Designed Csy4-array length (first site to last, from the SnapGene map):

| copies | mod1 | mod2 | mod3 | mod4 |
|---|--:|--:|--:|--:|
| 1× | 363 | 357 | **405** | 341 |
| 2× | 755 | 769 | **843** | 746 |
| 3× | 1,147 | 1,181 | **1,281** | 1,184 |

Bold is mod3, the architecture that was lost. Against the three siblings that
*did* assemble at the same copy number:

| | designed | against the three obtained | |
|---|--:|---|---|
| pART250, mod3 2× | 843 nt | 746 / 755 / 769 | **+74 to +97 nt longer** |
| pART251, mod3 3× | 1,281 nt | 1,147 / 1,181 / 1,184 | **+97 to +134 nt longer** |

The reason is the repeat unit itself. Per architecture, from the 1× → 2× step:
**mod1 392, mod2 412, mod3 438, mod4 405 nt** — mod3's unit is the **longest of
the four**, not the shortest. (The same arithmetic re-confirms the sibling
README's note that pART254's third unit is a mod3 unit and not a mod4 one: its
2× → 3× step is 438, mod3's, not mod4's 405.)

⚠️ **Which "shorter" this is.** The comparison drawn is **designed length
against the siblings that did assemble**, because the other reading has no
measurement behind it: nothing was ever built for pART250 or pART251, so there
is **no as-built length** to compare with their design. There is no shorter
product on disk for either. The panel therefore reports a designed size and its
direction, and says nothing about cause — n is one architecture, and a 74 to
134 nt excess over three siblings is not evidence that length is why the
assembly failed.

`NEED:` whether Leo meant **pART253** instead. That is the one construct in
this series with a *measured* shortfall: it assembled, gave the highest
first-run splicing ratio (47.8%), and was then dropped on full-plasmid
sequencing for extensive deletions (journal LV1.17, LV1.18). The LV1.18 ONT
consensus `LV1-18-001` is **2,497 bp against pART253's 3,079 bp design, 582 bp
shorter, with 0 of its 16 Csy4 sites left**; the second sample `LV1-18-002` is
6,043 bp with 3 sites, which is not a deletion product at all. ⚠️ **Neither can
be attributed yet**: `LV1.18_ONT/` holds no submission sheet, both `.Addendum`
files are empty, and the only reference map in the folder is pART253's, so
"LV1-18-001 is a pART253 clone" is an inference from the folder, not a record.
Confirm against the LV1.18 submission sheet before any of those numbers is
drawn or quoted. Until then panel (f) stays on the two assembly failures, where
every number comes from a file.

## Verified against the construct files

* The replacement is `ATTTCGATCAATGGATTAT`, **19 nt, 26.3% GC**, the same
  sequence at every position in both negatives. Found by walking the 18 nt
  flank either side of each Csy4 site from the positive into the negative, and
  asserted, not transcribed. It is **19 nt against the 28 nt it replaces — not
  the same length**, which is where every 9 nt in the arithmetic below comes
  from.
* **pART256** carries 8 replacements and 0 Csy4 sites; it is 72 nt shorter
  than pART252, which is exactly 8 × (28 − 19).
* **pART257** carries 14 replacements and 0 Csy4 sites, so seven F30 units at
  two positions each. It is 242 nt shorter than pART158.

### ⚠️ The F30 repeat period is 116 nt, not 98 (corrected 2026-11-19)

Panel (d) draws this. The decomposition it labels is the first of the two below.

An earlier version of this file decomposed pART257's 242 nt as "14 × 9 = 126
for the loops plus one whole 98 nt unit". That sums to 224, not 242; it mixes
two different accountings. Measured from the site positions instead:

| | site starts | gaps | period |
|---|---|---|--:|
| pART158, 16 Csy4 sites | 2298, 2346, 2414, … | 48, 68, 48, 68, … | **116 nt** |
| pART257, 14 replacements | 2298, 2337, 2396, … | 39, 59, 39, 59, … | **98 nt** |

So **98 nt is real but it is the *negative's* unit** (116 − 2 × 9), not the
positive's. Either decomposition of the 242 then closes exactly:

* 14 × 9 = 126 at the fourteen positions that remain, **plus one whole 116 nt
  positive unit** = 242 ✓
* 16 × 9 = 144 had all eight units been replaced, **plus one whole 98 nt
  negative unit** = 242 ✓

Both are right; mixing them is what gave 224. Journal LV1.34's working note
"3058 − 98 bp" is the negative unit, so it agrees with the 116 nt period rather
than contradicting it.

* 24 colonies at LV1.24 and 32 at LV1.28, **56 in total**, from the two
  submission sheets and from the read files, which agree.

### Which GC window the figure means

| window | length | GC |
|---|--:|--:|
| replacement | 19 nt | **26.3%** |
| Csy4 hairpin core `ACTGCCGTATAGGCAG` | 16 nt | **56.2%** |
| Csy4 extended repeat | 28 nt | **46.4%** |

"26% against 56%" is **the replacement against the 16 nt hairpin core**, and it
is right. An earlier note here called 56% wrong on the basis that the 28 nt
repeat is 46.4%; that compared the wrong window. Panel (a) is about the swap,
so it uses the 28 nt that actually comes out and **says so on the panel** —
"46% GC across the 28 nt replaced, against 26%" — rather than leaving a reader
to meet two different GC figures for the same element with nothing to tell
them apart.

### The colony screen reproduces all three journal counts

Unit count per colony is the replacement loop count divided by two, read off
that colony's own Sanger file. Two counts are kept: exact matches, and matches
allowing one miscalled base (a 700 nt near-identical repeat read turns one
copy odd easily). A colony is `clean` when the exact count is even, `rescued`
when the 1-mismatch count is even, `ambiguous` otherwise.

| | journal | re-derived here |
|---|---|---|
| LV1.24 colony 21 | "5× repeats total" | 10 loops, **5 units** |
| LV1.24 colony 23 | "9× repeats total" | 17 exact / 18 at ≤ 1 mismatch, **9 units** |
| LV1.28 colony 33 | "28.33 (7× as neg.)" → pART257 | 14 loops, **7 units** |

Distribution, clean plus rescued: plate 49761 gives 2 units ×3, 3 ×9, 4 ×3,
5 ×2, 6 ×5, 9 ×1, with one colony ambiguous; plate 49780 gives 1 ×1, 2 ×8,
3 ×12, 4 ×1, 5 ×6, 7 ×4, none ambiguous. **0/56 colonies (0%) carried eight
units**, which is what journal LV1.28's "No suitable 8× Neg" records, and why
a seven-unit clone was adopted instead.

`VERIFY:` resolved, and the brief's worry does not hold. "LV1.24 records 5×
and 9× for clones 21 and 23 but LV1.28 adopts 28.33 as the 7× negative" is not
a contradiction: **28.33 is the sample name LV1.28.33**, a colony on the second
plate, not a repeat count. Different colonies on different plates, all three
counts confirmed from the reads.

## Panel (b): the folds are drawn, not asserted

Leo: *"if you want you can also show the viennfold for the pos and neg
hairpin."* Done, and kept after looking at the render. What it replaced was
two literal dot-bracket strings, `(((((.....)))))` against `((((....))))`,
which at a glance say the two elements fold alike — the opposite of the truth
and of the MFE numbers printed beside them.

**ViennaRNA is reachable**, on Euler, at
`/cluster/work/wicky/lvalsecchi/rna_binder/pyfurnace_road/venv/bin/python`
(2.7.2, `import RNA` works). It is not in `~/.venvs/binder_map`, which is why
`make_inputs.py` has always shelled out to Euler for the folds. It also pulls
`RNA.naview_xy_coordinates()` for the two **elements** into
`inputs/supp_negative_controls_foldxy.csv`; `make_figure.py` draws from that.
ViennaRNA does the layout, matplotlib does the drawing in house colours, so
`RNA.svg_rna_plot()`'s own styling never enters the book. Inserts get no
layout — a naview plot of 938 nt is not a drawable object at 85 mm.

Both folds are drawn at **one** mm-per-naview-unit scale, so the 19 nt element
comes out visibly smaller than the 28 nt it replaces. Fitting each to its own
box would have hidden the 9 nt per site that identifies the crossed pair.
Rung weight is the pair chemistry — G–C heavy, A–U light, G·U dotted — which
is texture and not colour, and the rungs are neutral, so the encoding survives
a greyscale print.

### What the drawing shows

| | stem | pairs | MFE |
|---|---|---|--:|
| Csy4 repeat, 28 nt | 5 bp | 4 G–C, 1 A–U | −8.6 |
| replacement, 19 nt | 4 bp | 1 G–C, 2 A–U, **1 G·U wobble** | −0.9 |

**The replacement's MFE fold is a hairpin.** A marginal one — a single real
G–C, a wobble, −0.9 kcal/mol against −8.6 — but it is not a flat line, and the
panel draws it as what it is. This is the honest form of the claim the figure
exists to protect, and it is the stronger form: *the replacement cannot form
**the** hairpin*, and *the negative is **not** unstructured*. Both are visible
instead of asserted, and since 2026-11-19 they are visible in two panels rather
than one, which is what stops the first being read as the second.

### ⚠️ "3 bp, all A–T" and "4 pairs" are both right, under different rules

The longest self-complementary stem anywhere in the 19 nt is **3 bp and all
A–U** (`AUU`/`AAU`, `ATT`/`AAT` as DNA) **if only Watson–Crick pairs count** —
re-derived here two ways and asserted in `_selfcheck`. ViennaRNA's MFE reaches
**4** because it allows the **G·U wobble** at the stem foot. Neither number is
wrong; they answer different questions. Wherever the figure, the caption or the
Methods quotes one, **say which rule is meant** — a reader who meets "3 bp" in
the text and counts four rungs in the figure will otherwise conclude one of
them is a mistake. The panel names the wobble in the stem line beside the
drawing (`4 bp stem: 1 G–C, 2 A–U, 1 G·U wobble`), which is the rung a
Watson–Crick count drops, and the caption gives the 3 bp figure with its rule.

## The scaffold is a different answer from the element

| | nt | predicted pairs | MFE (kcal/mol) | per nt |
|---|---|---|---|---|
| Csy4 repeat alone | 28 | 5 | −8.6 | −0.31 |
| replacement alone | 19 | 4 | **−0.9** | −0.05 |
| pART252 insert | 391 | 110 | −138.3 | −0.35 |
| pART256 insert | 319 | **101** | −88.5 | −0.28 |
| pART158 insert | 938 | 337 | −476.7 | −0.51 |
| pART257 insert | 696 | **275** | −331.3 | −0.48 |

The 19 nt element on its own is barely folded at −0.9 kcal/mol. The negative
*scaffold* is not: the mod4 negative still folds 101 pairs over 319 nt, and the
F30 negative folds at −0.48 kcal/mol per nucleotide against its positive's
−0.51. **The negative removes the Csy4 site; it does not remove the fold.**

`NEED:` whether "essentially unstructured" in the thesis was ever meant about
the 19 nt element rather than the construct. If it was meant about the
construct, it has to go.

Folds are MFE only, on the scaffold insert taken from the first to the last
hairpin position plus 25 nt either side. No ensemble, no co-transcriptional
folding, no protein.

## The caution this figure must state once

**The replacement removes secondary structure as well as sequence at every
hairpin position.** A binder that scores on the target and not on the negative
has been shown to need the Csy4 site; it has not been shown to be reading the
site's sequence rather than its shape. That is one caution, in the caption and
in the Methods, not a paragraph. Drawing the folds makes it easier to see, not
less true: both elements fold, and the one that is removed is the one that
folds into the Csy4 hairpin.

## Open

* `NEED:` whether the 19 nt replacement was chosen for a defined property or
  is arbitrary filler. Nothing in the journal entries LV1.22 to LV1.34, in the
  construct files, or in the fold tables states a design criterion. Its 26.3%
  GC and −0.9 kcal/mol are consistent with "low GC, unstructured", but that is
  a description of what it is, not evidence of why it was picked. **Either
  answer is reportable; unstated is not.**
* `NEED:` the LV1.18 sample attribution, above — without it pART253's measured
  582 bp shortfall cannot be used.
* `VERIFY:` the pairing. pART256 is the negative for pART252 and pART257 for
  pART158, so the matched pairs in the screen are **pLV009 ↔ pLV012 and
  pLV010 ↔ pLV011**, crossed against the numbering. Confirmed three ways in
  `code/supp/scaffold_design/README.md`, and by arithmetic here: 9 nt shorter
  per hairpin position, 72 = 8 × 9.
* `VERIFY:` **two files in this repo call the pLV15 → pLV16 substitution an
  accident, and the transfection plan says otherwise.**
  `LV_writing/SI_DATA_BUILD_SPEC_2026-11-18.md` ("pLV16 was transfected in its
  place by accident, so posB was read against the wrong negative") and
  `code/supp/tables/README.md` repeat it. But `pLV15` occurs **0 times** in
  `LV2.33_transfection_plan.html` as well as in the as-performed record, where
  `pLV16` occurs **41** times and is the planned MM3; the as-performed
  deviations table lists pipetting volumes only, no plasmid swap. So the record
  shows a **design decision** — one scaffold negative for two different positive
  scaffolds — not a bench slip. Those two files should be corrected. It does not
  change panel (d), which draws pLV13 against pLV16 and says what each is.
* The consequence is worth carrying wherever `ONOFF_avg` is used: posB was read
  against an F30-array negative, so `PosB/NegA` divides an 8-site single-scaffold
  transcript by a 4-unit/8-loop F30-array one — different architecture, length
  and valency, not just ± hairpin. `code/supp/arrayed_flow/README.md` quantifies
  what that costs.
* `NEED:` an as-built sequencing record for **pLV13**. Only pLV16 has one in the
  013–017 range, and it turned out half-deleted, so pLV13's eight units are its
  **designed** count. Panel (d) draws pLV13 as designed and the bar would move
  if a record existed.
* Journal LV1.29 records that **ART255, the pART243 counterpart, was delivered
  but never used** (it has no tick, where ART256 does). pART243 therefore has
  no matched negative, which is one reason only pART252 and pART158 went
  forward.

## Caption draft

**The negative controls replace the Csy4 site, not the scaffold.** (a) One
hairpin position in pART252 and its negative pART256: a 28 nt Csy4 extended
repeat, 46% GC across the window replaced, becomes a 19 nt loop at 26% GC,
flanks unchanged. **The drawing shows one position; the swap is the same at all
eight**, which is what the panel's own line says. pART252 is 6,194 nt and
pART256 6,122 nt, 72 nt apart; pART256 **keeps no Csy4 site at all**, and the
72 nt is 8 × 9 — nine nucleotides at each of the eight positions. The shared flank is the mod4
scaffold. All four designed architectures share one blue (palette v0.9, theme step 600),
so **every construct is named with its architecture** rather than coloured by
it; the F30 array keeps the lighter Fig 1 blue, not being one of the four. (b) ViennaRNA MFE fold of
the two elements at one scale: the repeat forms a 5 bp stem, four pairs G–C, at
−8.6 kcal/mol; the replacement reaches 4 bp — one G–C, two A–U and one G·U
wobble — at −0.9, and its best **Watson–Crick** stem anywhere in the 19 nt is
3 bp, all A–U, the wobble being the fourth rung drawn. So the replacement
cannot form *the* Csy4 hairpin. (c) The construct does not follow: pART256's
insert still folds 101 pairs over 319 nt, −0.28 kcal/mol per nucleotide against
pART252's −0.35, and pART257 −0.48 against pART158's −0.51 — **the negative is
not unstructured.** (d) Each negative against the positive it was built from, as the
hairpin array, first position to the end of the last, at one scale, three pairs
in ascending shortfall. The mod4 pair kept every position: pART256 carries 8
replacement loops against pART252's 8 Csy4 sites and is 72 nt shorter than
pART252, 9 nt at each position. The F30 library pair did not: **pART257 came
back one whole F30 unit short of pART158**, 7 units and 14 loops against 8 and
16, and is 242 nt shorter than pART158 — one 116 nt F30 unit plus 14 × 9 nt at
the positions it kept; 0/56 (0%) colonies in the spike-in screen of (e) reached
eight units, so that unit was never placed. The arrayed screen's pair did worse:
**pLV16, the negative that screen used, carries 4 F30 units and 8 loops against
pLV13's 8 and 16**, and is 536 nt shorter than pLV13 — four 116 nt F30 units
plus 9 nt at each of its 8 remaining positions. That one is a different event,
which is what the two markers on the panel separate: pART257's unfilled unit is
marked **never placed**, pLV16's four **designed, not built**. pLV16 was
*ordered* valency-matched to pLV13 at 8 units and 16 loops (6,335 bp, itself
only 16 × 9 = 144 nt shorter than pLV13), and whole-plasmid sequencing of the
clone used gives 5,943 bp with half the array absent, 392 nt = 4 × its own 98 nt
unit; 144 + 392 = 536. Each shortfall decomposes two ways and both close: 242 =
one 116 nt positive unit + 14 × 9 nt, or one 98 nt negative unit + 16 × 9; and
536 = four 116 nt positive units + 8 × 9, or four 98 nt negative units + 16 × 9.
The panel draws and labels the first of each. Unfilled blocks are the units the
array came back without, one block per unit, drawn at the positive's 116 nt
period measured from its own site spacing; the hairline is where the positive's
array ends. The whole length difference of every pair lies inside the drawn
array span — 341 − 269 = 72, 888 − 646 = 242, 888 − 352 = 536, the same numbers
as the whole-plasmid differences. No cause is claimed for any shortfall. Counts
are from the `.gb` annotations — pLannotate's per-feature labels for the
as-built map — and agree with the motif counts in the sequence. (e) The spike-in
screen that explains (d): F30 units per clone, 56 colonies on two plates.
**0/56 (0%) carried eight units**. Four clones carried seven, and **one of those
four, LV1.28.33, became pART257** — the panel marks that bar, since it is the
link between this screen and the negative drawn in (d). (f) A different failure mode — two scaffolds whose assembly was never
obtained, against the three that assembled at the same copy number. Of twelve
designed combinations 9/12 (75%) were assembled, one fragment was never
delivered (mod3 1×, pART249) and the mod3 two- and three-copy assemblies were
not obtained. Length is the **designed** Csy4-array length — nothing was built
for either, so no as-built length exists — and it runs **opposite** to
expectation: mod3 carries the longest repeat unit of the four architectures
(438 nt against mod1's 392, mod2's 412 and mod4's 405), so pART250 and pART251
are 74 to 134 nt *longer* than the siblings that assembled at the same copy
number, not shorter. No cause is claimed in (d) or (f). The replacement removes
structure as well as sequence, so the negative does not separate sequence
readout from structure readout.
