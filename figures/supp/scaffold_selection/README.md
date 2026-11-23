# scaffold_selection — scaffold construction and selection

Three bare panels in `output/scaffold_selection/`, authored 85 mm wide
(a 38 mm tall, b 44, c 52). `style.save()` writes PDF, PNG and SVG; assembly
should consume the SVG, whose text is live. The tight bounding box emits them
at 245.6 / 248.3 / 248.3 pt wide (85 mm = 240.9 pt), so scaling to an exact
85 mm column shrinks type by ~3%, the same as every other panel in the repo.

```bash
~/.venvs/binder_map/bin/python code/supp/scaffold_selection/make_inputs.py
~/.venvs/binder_map/bin/python code/supp/scaffold_selection/make_figure.py
```

## Revision 2026-11-20 (Leo): palette v0.9, and panel (b) grouped and counted

> "Adjust colors. You could do a light background behind three bars of the same
> mod in b to bunch them together more visually. Please add to b how many csy4
> hairpins each column includes."

* **Colour.** Nothing in this script picks a colour; the change is
  `config/palette.json` moving `scaffolds` from `MAIN` `#015d82` to theme 600
  `#0678a6`. See the colour section below for the measured numbers. Looked at
  all three rendered panels afterwards: the lighter blue is still clearly
  separated from F30's pale fill in (b) and (c) and from the neutral-300
  negative controls in (c), the white dividers in (a)'s stacked bars still
  read, and the in-bar white `mod` ids in (a) still sit on a fill dark enough
  to carry them (4.95:1, down from 7.27 — the lightest the ids can take before
  they need an outline, so **do not take this blue lighter again** without
  re-checking panel (a)).
* **An architecture grouping band in (b).** `NEUTRAL["100"]` `#f0f6fa` behind
  each architecture's bars, `axvspan` at `zorder=0` with the bars at 2.
  Deliberately a **tint and not a hue**: the band groups, it does not encode,
  and a fifth colour in a panel whose colour already means "architecture"
  would read as a fifth series (COLORS.md rule 3). Matplotlib draws spines at
  2.5, so the axis is untouched, and the bars keep their 0.3 pt INK edge over
  it — F30's pale fill is 1.27:1 against the band (1.38:1 against white), so
  that edge is what holds it, exactly as STYLE_NOTES section 9 requires. The
  band stops 0.5 either side of the end bar, leaving the full 0.6 `GAP` of
  white between blocks. F30 gets a band of its own: it is an architecture
  block like the others, just a block of one.
* **A Csy4 hairpin count under every column**, as a third tick row, with the
  row named `Csy4 / hairpins` in the left gutter. **Derived, not typed**:
  `csy4_sites()` multiplies the copy number by `SITES_PER_UNIT`, 8 for a
  designed mod unit and 2 for an F30 repeat, which gives 8 / 16 / 24 across a
  mod block and 16 for pART158's eight units. `_check_sites()` then asserts
  every one of the ten against the counts taken off the plasmids themselves —
  `supp_scaffold_design_space.csv` (the nine mods, counted in the SnapGene
  maps) and `supp_negative_controls_design.csv` (`csy4_in_positive` = 16 for
  pART158, counted in the GenBank). All ten agree.
* ⚠️ **The gutter label is on two lines for a measured reason.** "Csy4
  hairpins" on one line is 11.6 mm at the house 5.5 pt against a 7.7 mm gutter,
  and the tight bounding box took the panel from 248.3 pt to **262.3 pt** —
  8.9% over the 85 mm column against the 2–3% every sibling panel carries, so
  its type would have come out visibly smaller than theirs in the book. Broken
  in two the widest line is "hairpins", 6.6 mm, and the panel is back at
  **248.3 pt**, the same width as panel (c). Measured with `style.text_w`.

## Revision 2026-11-19, second round (Leo)

> "a is good but rename sequence checked to something like signal strength
> measured / b is good but just put the ratios in 0.xx format / c same for this
> also why is there a range for the ratios? should only be one fixed number"

* **Panel (a)'s fourth stage is now "Splicing measured".** Checked before
  renaming, and **sequence-checking and measuring are two different records**:
  `LV_Sequencing/LV1.16_SangerSeq/` holds one `.dna` per assembled construct,
  and the LV1.15 ONT run holds read counts. They cover the **identical nine**
  (pART243–248, 252–254, with pART158 alongside in both), which is why the bar
  was 9/12 under either name. Drawing one bar therefore loses no construct, and
  the stage now counts the readout the rest of the figure is built on. The
  caption states that the same nine were Sanger-verified, so the QC step is
  recorded rather than deleted. "Splicing measured" rather than Leo's literal
  "signal strength measured": the measured quantity is the splicing ratio, and
  he wrote "something like".
* **The ratio is printed as `0.xx`** above every bar in (b) and (c), with the
  raw counts under it (STYLE_NOTES.md section 5 wants the inputs behind a derived
  value). Sizes are **measured, not estimated** (`style.text_w`): "408/1009" is
  5.80 mm at 4.0 pt against panel (b)'s 6.03 mm bar pitch and 6.53 mm at
  4.5 pt, which is why (b)'s counts sit at 4.0 pt and (c)'s, with 12.1 mm of
  pitch, keep 4.5. Panel (b)'s `ylim` went 0.56 → 0.62 for the second line and
  its left spine is now bounded at 0.4, as (c)'s already was.
* **The fold-change brackets print one number.** See *The fold changes on panel
  (c)* below — the parenthesis was a confidence interval, not a range, and it
  moves to the caption.
* **Colours are the `scaffolds` block, and F30 is out of it.** See the table
  below. (The step moved to theme 600 on 2026-11-20,
  palette v0.9.)

## Revision 2026-11-19 (Leo)

> "good but group the mod repeats with colors. A is good b needs refinement. No
> need to show the correct pairing of ont an dpart in the plot (no lines) just
> put the correct part nr. Also show the fc ratios between plots in the c"

* **Colour is the architecture, in all three panels.** mod1 / mod2 / mod3 /
  mod4 and the pre-existing F30 each hold one step of the theme ramp, the same
  step everywhere, so a reader sees the three-copy collapse as a property of
  copy number *within* an architecture rather than as twelve unrelated bars.
* **Panel (b)'s bow-tie strip is gone.** Each bar now carries its corrected
  pART number under its copy number. The correction itself is asserted in
  `_check_inversion()` and stated in the caption, not drawn.
* **Panel (c) carries the fold changes**, with confidence intervals, and
  handles the two zero denominators as a bound rather than an infinity.
* **Panel (a) is structurally unchanged** — same five stages, same bar lengths,
  same `n/12` labels. Its bars are now segmented by architecture, which is the
  only way the colour key means anything in (a), and it is what makes the loss
  of the whole mod3 series visible. Reverting is one `barh` call.

### Which colour block, and why not the stack colours

COLORS.md rule 2 reserves the five `stacks` colours for **design-route identity
only**; an RNA scaffold is not a design route, so `palette.stack_color()` is
deliberately not called here. The `synbio` block is molecular identity, the
`refolders` block is the re-predictors — measurement instruments — and
`thesis_div` is for signed data; none of them fit either. **Palette v0.7
(2026-11-19) added the block this figure actually wants**, `scaffolds`
(COLORS.md rule 9), and it replaced the theme-petrol steps used until then.
**Current values, palette v0.9 (2026-11-20):**

| | source | hex |
|---|---|---|
| mod1, tandem compact | `palette.SCAFFOLDS["mod1"]` | `#0678a6` |
| mod2, branched | `palette.SCAFFOLDS["mod2"]` | `#0678a6` |
| mod3, tandem splayed | `palette.SCAFFOLDS["mod3"]` | `#0678a6` |
| mod4, ring | `palette.SCAFFOLDS["mod4"]` | `#0678a6` |
| F30 array (pART158) | `palette.F30` (= `THEME["200"]`) | `#bfe0f5` |
| architecture grouping band, panel (b) | `palette.NEUTRAL["100"]` | `#f0f6fa` |
| negative controls, panel (c) | `palette.NEUTRAL["300"]` | `#cfd5d9` |

⚠️ **The four are ONE colour, not four steps — and v0.9 moved that colour.**
v0.7's four-step indigo ramp was withdrawn on 2026-11-20 (adjacent ΔE was only
9.8, and all four sat too dark to tell apart); v0.8 made all four `MAIN`
`#015d82`; **v0.9 moved the single step to theme 600 `#0678a6`** on Leo's
instruction — *"adjust the rna presentation modules blue to a lighter one"*.
Measured on the new step: contrast on white **4.95** (7.27 before), ΔE **44.5**
from `palette.F30`, **48.8** from the negative controls' neutral 300 and
**58.3** from the grouping band.

**Any ΔE 9.8 still quoted below is the v0.7 ramp and no longer exists** — the
four are now ΔE 0 from each other, which makes the label-on-every-mark rule
below stronger, not weaker, and makes the white divider in panel (a) load-
bearing rather than merely helpful.

**F30 is deliberately not one of them.** The four architectures share one
colour because they are four variants of the same kind of object; the
pre-existing F30 array is not a variant of them, and it keeps the blue the main
figures' construct maps give it — ΔE 44.5 from the four. The old scheme gave it
theme 900, the mods' own hue family, which is exactly what Leo asked to change.
Negative controls stay neutral 300 (rule 3, a reference rather than a series).

The four architectures are now the SAME colour, so COLORS.md rule 9 and
STYLE_NOTES.md section 9 both demand **a label on every mark**. Supplied:
architecture blocks are separated by a gap on
the axis and named *in full* underneath, every bar carries its pART number, and
panel (a) keys all four architectures over its top bar at 5.0 pt (measured: the
widest name is 12.1 mm against a 13.3 mm segment, and 13.4 mm at the house
5.5 pt). The pale F30 fill keeps its hairline INK edge, which section 9
requires of any fill under 3:1 against white. Identity never rests on colour
alone (checklist item 4).

Panel (a)'s stacked segments are divided by a **white** hairline, not INK:
four identically coloured blocks with a near-black edge between them read as
one block, and the architecture split — the whole reason (a) is segmented —
disappears. Panels (b) and (c) keep the INK edge; their bars stand alone (on
the grouping band in (b), on white in (c)).

**Every segment in (a) carries its own `mod` id**, in white inside the bar at
4.0 pt (measured: "mod4" is 3.37 mm against the 4.43 mm one-unit segment in
"Carried forward"; 3.79 mm at 4.5 pt). This was added 2026-11-19 because the
key alone is wrong for three of the five rows: as soon as one architecture goes
to zero the segments after it slide left, so in "Assembled" and "Splicing
measured" the mod4 segment sits under the "tandem, splayed" key, and in
"Carried forward" it sits under "tandem, compact". The geometry is unchanged and
inherent to a stacked funnel; with all four the same colour a reader cannot
recover the mapping from colour at all, so the label is repeated rather than
inferred. The full names stay
on the key above the top row, where each sits directly over its own id.

## Sources

| source | mtime | what it gives |
|---|---|---|
| `LV_Sequencing/LV1.15_ONT/RawReadsAnalysis/LV1_15_AllResults.csv` | 2026-08-27 11:20 | first run, raw counts |
| `…/LV1_15_AllResults_CorrectedpART.csv` | 2026-08-31 13:38 | the corrected sample → construct assignment |
| `LV_Sequencing/LV1.27.1_ONT/RawReadsAnalysis/LV1_27_1_AllResults.csv` | 2026-09-06 15:10 | repeat run, raw counts |
| `…/LV1.27.1_ONT/RawReadsAnalysis/run_alignment.sh` | — | the counter and its settings |
| `LV_Cloning/Projects/F30_mod/Assembled_Final/` | — | the 12 SnapGene maps; their names carry the architecture and copy number |
| journal LV1.17, LV1.24, LV1.28 (Notion, LV Wet-Lab Journal) | 2026-10-14 | the inversion, the repeat's transfection list |

Derived: `inputs/supp_scaffold_selection_ont.csv`.

**The ONT outputs are the source, not the chart pasted into the journal.** Both
NAS plot folders also hold a `*_ByFamily_thesis.svg`; those were read only to
recover the repeat run's construct labels and are not an input here.

## Verified

Every ratio below is `removed / (intact + removed)`, recomputed from the counts
and asserted against the counter's own `total` column, which excludes
`ambiguous_splicing` (LV1.15.001 is 426 + 299 = 725 with one ambiguous read set
aside). `ambiguous_transcript` is excluded by the counter.

**First run (LV1.15), corrected assignment, 3.1% to 47.8%:**

| construct | | removed/evaluable | ratio |
|---|---|---|---|
| pART243 | mod1 1× | 408/1,009 | 40.4% |
| pART244 | mod1 2× | 65/804 | 8.1% |
| pART245 | mod1 3× | 18/582 | **3.1%** |
| pART246 | mod2 1× | 214/738 | 29.0% |
| pART247 | mod2 2× | 322/1,021 | 31.5% |
| pART248 | mod2 3× | 50/758 | 6.6% |
| pART252 | mod4 1× | 567/1,555 | 36.5% |
| pART253 | mod4 2× | 351/734 | **47.8%** |
| pART254 | mod4 3× | 37/747 | 5.0% |
| pART158 | F30, 8 units | 299/725 | 41.2% |

Three-copy constructs collapse in every architecture measured (3.1%, 6.6%,
5.0%), which is the one pattern that survives both runs.

**Repeat (LV1.27.1), on 267 to 329 evaluable reads — all five confirmed:**

| construct | removed/evaluable | ratio | uncertainty |
|---|---|---|---|
| pART158 | 117/300 | 39.0% | ± 2.8 pp |
| pART252 | 99/279 | 35.5% | ± 2.9 pp |
| pART243 | 82/287 | 28.6% | ± 2.7 pp |
| pART215, negative control, 4 hairpin positions | 0/329 | 0.0% | < 0.9% |
| pART242, negative control, 12 hairpin positions | 0/267 | 0.0% | < 1.1% |

Binomial standard error, `sqrt(p(1−p)/n)`, hand-rolled: scipy is broken in this
venv (`_propack` dlopen). For the two zero counts an SE is meaningless, so the
panel gives the one-sided 95% rule-of-three bound `3/n` instead. **At ~300 reads
the error on a ratio near 0.4 is 2.8 to 2.9 percentage points**, so the
39.0% against 35.5% gap is inside it: pART158 and pART252 are not separated by
this experiment, which is consistent with the journal carrying both forward
rather than picking one.

## The fold changes on panel (c)

Ratios of two proportions, with a two-sided 95% CI from the delta method on
`log(p1/p2)` — `var = (1−p1)/(n1·p1) + (1−p2)/(n2·p2)`, z = 1.96 written out
because there is no `scipy.stats.norm` to ask (`fc_ci()`, with a `_selfcheck()`
that fails if it breaks).

| comparison | fold change | 95% CI | |
|---|---|---|---|
| pART158 / pART252 | 1.10× | 0.89 – 1.36 | contains 1 |
| pART252 / pART243 | 1.24× | 0.98 – 1.58 | contains 1 |
| pART158 / pART243 | 1.37× | 1.08 – 1.72 | 1 outside — *not drawn*, see below |

**The bracket prints the point estimate only** (changed 2026-11-19). It used to
read `1.10× (0.89–1.36)` and Leo read the parenthesis as a range of values —
"why is there a range for the ratios? should only be one fixed number". It was
never a range: there is **one measurement per construct per run**, so there is
no spread across replicates to collapse, and the parenthesis was an interval
*estimate* of that one number. So the number stays single on the panel and the
interval moves to the caption, where a word of prose can say what it is. The
uncertainty is not dropped — the binomial standard errors are still drawn on
every bar, `report()` prints every interval on every run, and `panel_c` now
**asserts** that both drawn comparisons still contain 1, so the caption cannot
quietly go stale. Both do, so this experiment does not rank the target
scaffolds. The 158/243 comparison is
the one that clears 1, and it is left to the caption rather than drawn as a
third nested bracket: it is the weakest and the strongest of three conditions,
it is not the comparison the selection turned on, and a third bracket costs
more ink than it returns.

**The two negative controls spliced 0/329 and 0/267, so a ratio against them is
undefined.** The panel prints a *lower bound* instead of an infinity: divide by
the **weaker** of the two one-sided 95% bounds (3/267 = 1.12%, the larger
number, hence the most conservative denominator) and take the **weakest**
target, so the single figure printed holds for every target against both
controls.

| target | vs 3/267 = 1.12% |
|---|---|
| pART158 | > 34.7× |
| pART252 | > 31.6× |
| pART243 | > 25.4× |

The panel therefore reads **"> 25×, vs 95% control bound"**, which is true of
all six target × control pairs.

### An alternative reading of "the fc ratios between plots"

Read literally, "between plots" could mean between panels (b) and (c) — the
three constructs measured twice. Those ratios are **not** drawn; if that is what
was wanted, say so and they are two lines of code:

| construct | run 1 (LV1.15) | run 2 (LV1.27.1) | run2 / run1 |
|---|---|---|---|
| pART158 | 41.2% | 39.0% | 0.95× |
| pART252 | 36.5% | 35.5% | 0.97× |
| pART243 | 40.4% | 28.6% | 0.71× |

`report()` prints them on every run so the comparison is never lost.

## The sample inversion, and what the files say about it

Journal LV1.17 records it as *"the samples were inverted 1 ↔ 10"*. The
corrected assignment file reverses the **whole** series: sample k is construct
11 − k. Three independent spot checks in LV1.17's own table agree with the
reversal, not with a single swap (LV1.16.001 → pART158, the 10th construct;
.004 → pART252, the 7th; .010 → pART243, the 1st). **Files win: the correction
is a full order reversal.**

It is no longer drawn. `_check_inversion()` asserts it instead — `RUN1_ORDER[i]`
must be sample `10 − i` — so if the input CSV ever reverted to the uncorrected
file, every bar in panel (b) would be mislabelled and the script would fail
rather than render.

`NEED:` resolved. The corrected splicing-ratio graph LV1.17 refers to **does
exist** — `LV1.15_ONT/RawReadsAnalysis/Plots/
LV1_15_SplicingRatio_CorrectedpART_HighQuality.{png,svg}`, written 2026-08-31,
the same day as the corrected CSV. It was simply never pasted into the journal
frame. Nothing is missing; the journal is.

`VERIFY:` the repeat run's five construct labels. Every sample in that run is
aligned to one merged reference (`MLS12_18_1_AllSamples`), so the counter's
output cannot say which plasmid a sample was, and the mapping in
`make_inputs.py` is journal LV1.17's transfection list in its own order. It
agrees with the labels on `LV1_27_1_SplicingRatio_ByFamily_thesis.svg`, but
both trace to the same hand. Check against the LV1.21 platemap before this
reaches the thesis.

## Panel (a) stage definitions

Each stage counts over the same 12 designed combinations — four architectures
at one, two and three copies — and each is traceable to files. The
per-architecture split read by the segments comes from the 12 SnapGene map
names in `LV_Cloning/Projects/F30_mod/Assembled_Final/`: **243–245 mod1,
246–248 mod2, 249–251 mod3, 252–254 mod4**.

* **Designed, 12/12** — 12 maps, three per architecture.
* **Fragment delivered, 11/12** — `pART249_notcloned.dna` is the exception, so
  mod3 loses its one-copy member.
* **Assembled, 9/12** — the 9 that appear both in the `LV1.16_SangerSeq`
  submission and in the LV1.15 ONT run. The whole mod3 series (pART250,
  pART251) is absent from both, independently.
* **Splicing measured, 9/12** — the nine designed constructs in the LV1.15 ONT
  run, which is panel (b). Renamed from "Sequence checked" 2026-11-19.
  **The two are not the same record**: `LV_Sequencing/LV1.16_SangerSeq/` holds
  one `.dna` per assembled construct (ten files — the nine plus pART158), and
  the ONT run holds read counts for the same ten. They cover the identical
  nine, so the bar is 9/12 either way and nothing is lost by drawing the one
  the figure is built on; the caption records the Sanger verification.
* **Carried forward, 1/12** — pART252 (mod4), plus pART158, which already
  existed and is therefore not one of the 12. The panel says so directly.

`VERIFY:` "assembled" is inferred from two downstream records, not from a build
log. The written record (context snapshot, 2026-11-05) says assemblies 7 and 8,
pART250 and pART251, failed; no assembly log was found on the NAS to confirm
it.

`NOTE:` **pART254's third unit is a mod3 repeat, not a mod4 one.** Its map is
named `(F30_mod4…)x2_(F30_mod3…)`, so the mod4 "3×" bar is three repeat units
but not three identical ones. It is still a three-copy construct and it still
collapses (5.0%), so the claim the caption makes holds; but if the panel ever
needs to say "three copies of the same module", this one does not qualify.
Found 2026-11-19 while reading the map names for the panel (a) split.

## Deliberate

* pART253 gave the highest first-run ratio (47.8%) and is **not** carried
  forward: LV1.17 and LV1.18 record extensive deletions on full-plasmid
  sequencing. The panel plots it at its measured value and the caption says
  why it stops there.
* The negatives are grey, not a hue: they are a reference, not a series
  (COLORS.md rule 3).
* No hatching anywhere. COLORS.md rule 7 records the texture encoding being
  withdrawn 2026-11-19 — "remove the dots and the stripes" — so F30 is
  separated from the mods by a ramp step, not by a fill pattern.
* Panel (c)'s left spine is bounded at 0.4, where the ticks stop, so the axis
  does not run up through the space the fold-change brackets use.
* A splicing ratio is a reporter proxy for **receptor dimerisation** on a
  multivalent scaffold, and so for **avidity**. The binder is Csy4-**H29A**,
  catalytically dead (H29 is the catalytic-dyad residue contacting the hairpin
  in 4AL5), fused through a (GGS-GGT)x2 linker to **MgIRE1** — so the readout
  cannot be cleavage of the scaffold. Binding the hairpin array brings the
  MgIRE1 halves together and their reconstitution is what splices. It is not an
  affinity, and nothing in this figure is a Kd.
  ⚠️ Corrected 2026-11-19: this README and the script docstring previously said
  "a reporter proxy for Csy4 cleavage of the scaffold", which attributes the
  signal to an enzymatic activity the construct deliberately abolishes.

## Caption draft

**Nine of twelve designed scaffolds were obtained, and two were carried
forward.** (a) Attrition across the twelve designed combinations — four
architectures (mod1–mod4) at one, two or three copies. Segments are
architectures — mod1 a compact tandem array, mod2 branched, mod3 a splayed
tandem array, mod4 a ring — keyed above the top bar; one fragment was never delivered
and the whole mod3 series was lost, leaving nine. Those nine were both
Sanger-verified (LV1.16) and measured by ONT (LV1.15); the stage counts the
measurement. pART252 was carried forward
alongside the pre-existing pART158. (b) First ONT readout (LV1.15): splicing
ratio per construct, printed above each bar with the removed/evaluable reads it
came from beneath it, and copy number, corrected pART identifier and the number
of Csy4 hairpins the construct presents below the axis. Colour is the
architecture and a light band groups each architecture's bars. Three-copy
constructs collapse in all three architectures measured. Labels carry the
LV1.17 correction: the ten ONT samples had been assigned to constructs in
reverse order, so the first analysis was wrong throughout and (c) exists
because of it. (c) Repeat readout (LV1.27.1), five conditions, 267–329
evaluable reads. Bars are labelled with the splicing ratio and, beneath it, the
removed/evaluable reads it was computed from. Error bars are binomial standard
errors. Brackets give fold changes as **point estimates**; their **95%
confidence intervals** are 0.89–1.36 for pART158/pART252 (1.10×) and 0.98–1.58
for pART252/pART243 (1.24×). **Both include 1**, so this experiment does not
separate the two target scaffolds and both were carried forward. The controls
spliced 0/329 (0.0%) and 0/267 (0.0%), so a ratio against them is undefined;
each control bar carries its one-sided 95% rule-of-three bound and the bracket
gives the lower bound that holds for every target against both controls. Splicing ratio is a reporter proxy for receptor dimerisation on a
multivalent scaffold, and so for avidity, not an affinity.
