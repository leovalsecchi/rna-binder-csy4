# scaffold_design — RNA scaffold design

Four bare panels in `output/scaffold_design/`, authored 66 / 178 / 178 / 85 mm
wide; (c) is 88 mm tall, being two rows of folds at the figure's one `SCALE`.

The tight bounding box adds `style.save`'s fixed pad, so the emitted files are
200.9 / 518.4 / 519.3 / **254.9** pt against 187.1 / 504.6 / 504.6 / 240.9 pt
for those widths — 7.4 / 2.7 / 2.9 / **5.8%** over. The fixed pad costs a
narrow panel proportionally more, which is why (a) and (d) are the two over
5%; (d) grew on 2026-11-20 when panel (d)'s right-hand key widened, and its
content now reaches to 82.9 mm of the 85 mm canvas (asserted in `panel_d`).

```bash
~/.venvs/binder_map/bin/python3 code/supp/scaffold_design/make_inputs.py
~/.venvs/binder_map/bin/python3 code/supp/scaffold_design/make_figure.py
```

`make_inputs.py --skip-folds` leaves the two Euler steps out; the fold check
takes about three minutes because ViennaRNA folds four 1,200 nt sequences.

> ✅ **Which panel draws which design generation — SETTLED, verified on disk
> 2026-11-20.** **All four panels are now on the ORDERED second-generation
> D-series** (22-JUL-2026; units 396 / 409 / 442 / 402 nt). (b) prints those
> four lengths, from `supp_scaffold_design_ordered.csv`, and its mod2 tree is
> parsed from the ordered target; (c) is a single panel,
> `panel_c_ordered_fragments`, drawing each ordered 1× fragment twice — design
> target above, RNAfold prediction below; (d) folds the twelve ordered
> constructs. **`panel_c_target_structures` is no longer emitted** and the
> first-generation drawing is in `archive/`. The earlier state of this note —
> "(b) is still the superseded first generation", "two (c) panels, both
> emitted, Leo picks" — is **stale and was corrected here**, not re-done.
> Only mod4 is the same in both generations, its geometry being unchanged.

## Sources

| source | what it gives | mtime |
|---|---|---|
| `LV_Cloning/Projects/F30_mod/Assembled_Final/pART2{43..54}*.dna` | the 12 designed constructs, as SnapGene maps | 2026 (unchanged since the build) |
| `LV_Cloning/pLV_NNN/pLV009/pLV009_POS_A_example.fa` | the ordered eight-unit F30 array, panel (a)'s sequence | 2026 (unchanged since the build) |
| `LV_Cloning/Projects/F30_mod/WIP/F30_final_assembly/F30_D{1..4}-(N x Csy4)-f{1..3}.gb` | **the 24 ORDERED fragments** - panel (d) | 22-JUL-2026 |
| `inputs/dseries_targets/mod{1..4}.dbn` | their design targets, vendored | 2026-11-20 |
| `inputs/csy4_dotbracket/F30_mod{1..4}_pos.dbn` | the FIRST-GENERATION targets, superseded; still what (b) and (c) draw | 2026-10-05 09:10-09:11 |
| Euler `rna_binder/pyfurnace_road/runs/*/{best_seq,target}.txt` | every PyFurnace/ROAD run, re-folded here | 2026-07 to 2026-10 |
| Euler `pyfurnace_road/venv` (ViennaRNA 2.7.2) | the MFE fold and the naview layouts | — |
| `Fig1f_scaffold_splicing/README.md` | the D-series is a different target set from the mod series | 2026-11-17 |

Derived tables, all written by `make_inputs.py`:
`inputs/supp_scaffold_design_space.csv`, `..._targets.csv`, `..._folds.csv`,
`..._layouts.csv`, `..._ordered.csv`.

## The architectures are NAMED, not just numbered — renamed 2026-11-20

| token | name | design file |
|---|---|---|
| mod1 | **tandem, compact** | `D1_tandem_lean` |
| mod2 | **branched** | `D2_branched_lean` |
| mod3 | **tandem, splayed** | `D3_splayed_tandem` |
| mod4 | **ring** | `F30_mod4` |

Note the order: mod2 and mod3 are not what the numbering suggests - **mod2 is
the branched one and mod3 is the second tandem array**, so mod1 and mod3 are
the pair. The name is printed on **(b), (c) and (d)**, and the `modN` token
stays beside it everywhere, because the files, the SnapGene maps and the
journal all use the token. `make_figure.py` asserts that every architecture has
a name - under palette v0.8 and v0.9 alike all four share one colour, so the
name is the only thing telling them apart.

### Why they were renamed (do not re-derive this in six months)

They were first named from Leo's own words on 2026-11-19: *"Mod1 is f30 long
arms and Mod3 is f30 short arms and 2 is branched and 4 is ring ro circular"*.
Those words describe the **superseded first-generation** designs. When the
panels moved to the designs that were actually ordered, the names broke twice
over:

* **The difference they name no longer exists.** Measured off the layout tables
  - base pairs in the helix carrying each Csy4 site - the first generation had
  mod1 at 8.5 bp against mod3 at 5.5. **In the ordered designs all four
  architectures are at 5.0 bp.** The low-structure redesign spent the arms to
  meet a budget of at most 30 bp in any 120 nt window
  (`csy4_redesigns_lowstructure.txt`: *"Tandem 4 modules, stems 5/5/5"*).
* **And they invert.** mod3 is now the LARGER fold (442 nt, 18.5 mm of reach)
  and mod1 the smaller (396 nt, 13.3 mm) - the reverse of the first generation,
  so a reader taking "short arms" at face value reads the panel backwards.

What separates them now is the linker between modules - per-site spacings
`41 57 …` for mod1 against `38 75 …` for mod3 - and the size of the exterior
loop it opens: mod1's core radius 6.7 mm against mod3's 10.5 mm. "Compact" and
"splayed" name that, and they are the design files' own words. COLORS.md rule 9
carries the same correction.

## Colour: palette v0.9, one colour for the four, and F30 keeps its own

`palette.SCAFFOLDS` is **one colour for all four architectures**, and since
**palette v0.9 (2026-11-20)** that colour is **theme 600 `#0678a6`** rather
than `MAIN` theme 700 `#015d82`. Leo: *"adjust the rna presentation modules
blue to a lighter one."* **Only the step moved** — the one-colour decision of
v0.8 (and the withdrawal of v0.7's four-step indigo ramp) stands, and the four
keys in the `scaffolds` dict still hold four identical values on purpose.
Measured for the amendment: contrast on white 7.27 → **4.95**, ΔE **11.8** from
the old step (so the change is visible), ΔE **44.5** from `palette.F30` and
**21.2** from `synbio.rna`; theme 500 was rejected at ΔE 11.1 from `synbio.rna`
and theme 400 *is* `synbio.rna`. COLORS.md carries the amendment under rule 9
and `config/palette.json` is the single source of truth.

**Colour therefore separates nothing in this figure**, and the drawn fold plus
the printed NAME carry the architecture on (b), (c) and (d) alike.
`make_figure.py` asserts `palette.SCAFFOLD_ORDER == ARCH` and that every
architecture has a name, so the one thing now telling them apart cannot go
missing.

Looked at all four rendered panels under v0.9: the lighter blue still reads
against white at every weight the figure uses — the 1.0 pt glyph strokes in
(b), the 0.5 pt fold backbones in (c) — and it stays clearly apart from the
warm `#ae4e1a` Csy4 repeats drawn over it and from F30's pale `#bfe0f5` in (a).
Nothing in the figure needed re-tuning for the step.

The F30 array keeps `palette.F30` (`THEME["200"]`), the blue the Fig 1 construct
maps give it, with its Csy4 sites in the diverging warm end - so panel (a) reads
as F30 and not as a designed architecture.

### ✅ DECIDED (Leo, 2026-11-20): F30's pale blue stays. Do not re-raise.

Leo: *"I want you to keep the F30 the color that it is in the main."*

`palette.F30` is `#bfe0f5`, which measures **1.38:1** against white, and panel
(a) draws it as a bare 0.75 pt stroke with no edge. That is below what
STYLE_NOTES section 9 asks for (a hairline edge and a visible label under 3:1)
and outside what COLORS.md rule 1 reserves steps 100-300 for ("tints for fills
and backgrounds only; 400 to 900 for data marks") - in `constructs.py` F30 is a
filled bar *with* a `synbio.outline` edge, which is why the same colour works
there. Under v0.8 it sat beside architectures drawn in the much darker `MAIN`;
under v0.9 those are theme 600 `#0678a6`, still 44.5 ΔE away, so the array
still reads as a different object from them.

**It stays anyway, deliberately.** F30's colour is not chosen in this figure; it
is **inherited from the main text, where F30 is the target module**, and a
reader moving from Fig 1 to this supplement has to see the same object in the
same colour. Cross-figure identity wins over local contrast here. `THEME["300"]`
(1.78:1) and `THEME["400"]` (2.44:1) were the palette-legal alternatives and
were **rejected** - do not switch to them.

## All five scaffolds are drawn as folds (2026-11-19)

Leo: *"I want to have all 5 scaffolds (including the f30 with csy4 repredicted
in viennafold to use in the figure) use the dotbrace representation to get the
nice folds. A is unusable. b is good. c is not usable and d is good but has to
be refined."*

Panel (a) was a hand-drawn three-way junction and panel (c) was an arc diagram.
Both are now **secondary-structure drawings laid out by ViennaRNA's naview**,
the layout `RNAplot` itself uses — stems as ladders, loops opened out. Panel (b)
is untouched.

`inputs/supp_scaffold_design_layouts.csv` holds one row per nucleotide
(`scaffold, i, base, x, y, pair, csy4`) for all five, so the figure redraws
without Euler. ViennaRNA is **not** importable in `~/.venvs/binder_map`, so the
layout runs on Euler; the data travels inside the script on stdin, because the
login banner breaks a separate stdin stream (`euler-ssh-transfer-quirks`).

Of the three engines the 2.7.2 Python module exposes, `simple` overlaps its own
arms on these targets and `turtle` emits one polyline per element rather than
one per chain. naview is both the nicest and the least code.

**Panel (a) is a real fold of the ordered sequence, not of a target.** The
eight-unit F30 array of pLV009 is folded whole (`RNA.fold`), and the drawn unit
is the **outermost closed domain holding exactly two Csy4 sites**, taken from
the repeat nearest the middle of the array so no dangling end is in it. That is
array positions 432–569, **138 nt**, 52 bp, 2 sites — a three-way junction with
a closing stem and two Csy4 arms, which is the two-sites-per-unit claim drawn
rather than asserted.

**The repeat period is 116 nt, and 138 is not a second value for it.** The
sixteen Csy4 sites sit at alternating 48 / 68 nt spacings, so one unit is 116 nt
and the array is 8 × 116 = 928 nt. The folded domain is 138 nt because its
closing helix reaches 11 nt into each neighbouring period. `make_inputs.py`
measures the spacings and asserts them rather than assuming.

**This corrects this README's own "98 nt F30 unit".** 6,427 − 6,185 = 242, and
98 + 14 × 9 = 224, which is 18 short. 116 + 126 = 242 exactly. pLV009 carries 16
Csy4 repeats in 6,427 nt and pLV012 carries 14 negative loops in 6,185 nt (both
counted from the `_example.fa` files), so the A pair differ by one whole 116 nt
unit plus fourteen 9 nt loop swaps. The crossed-pairing conclusion below is
unaffected — only the unit length was wrong.

## The VERIFY that other figures depend on: which scaffolds were ordered, and how many sites

**Settled from the DNA files, three independent ways. Leo's "both carried-forward
constructs carry 16 Csy4 sites" is wrong, and Fig. 1F's "8" for the circular
scaffold is right.**

| plasmid | role | Csy4 sites | length |
|---|---|---|---|
| pLV009 | POS_A, from pART158 (F30, 8 tandem units) | **16** | 6,427 nt |
| pLV010 | POS_B, from pART252 (mod4, 1×, the ring) | **8** | 5,899 nt |
| pLV011 | NEG_B, from pART256 | 0 | 5,827 nt |
| pLV012 | NEG_A, from pART257 | 0 | 6,185 nt |

Counted by searching the Csy4 extended repeat `GTTCACTGCCGTATAGGCAGCTAAGAAA` in
`LV_Cloning/pLV_NNN/pLV0{09..12}/*_example.fa`. The same count on the parent
plasmids gives pART158 16/6,722 nt and pART252 8/6,194 nt, which are the two
figures Methods 2.12.2 already reports, so the two records agree.

Copy number and site count are one axis, and the SnapGene maps say so exactly:
1× = 8 sites, 2× = 16, 3× = 24, in all four architectures (panel d). pART158 is
not a mod construct at all: it is eight tandem F30 units at two sites each.

**The crossed pairing is confirmed by sequence arithmetic, not only by the file
names.** The negative replaces a 28 nt Csy4 repeat with a 19 nt loop at every
hairpin position, so a negative is 9 nt shorter per position than its own
positive. 6,194 − 6,122 = **72 = 8 × 9** (pART252/pART256, the B pair), and the
same 72 separates pLV010 from pLV011. The A pair differ by 242, which is
14 × 9 = 126 for the replaced loops plus one whole 116 nt F30 unit, pART257
being a seven-unit clone against pART158's eight. So pLV009 pairs with pLV012
and pLV010 with pLV011, and pairing by adjacent number compares each scaffold
against the wrong negative. (The A-pair arithmetic is 126 + **116**, one whole
F30 unit; see the unit-length correction above.)

## Panel (d): the twelve ORDERED constructs, folded at their real copy number (2026-11-20)

Leo, on the grid of "not run" cells: *"not sure this is true we have the data
for at least nine"*. He is right, and the data is better than "at least nine" -
**every one of the twelve cells is a real designed sequence that was actually
ordered.**

**The ordered material is the D series, not `F30_mod*_pos`.** Twenty-four
GenBank files of 22-JUL-2026 in
`LV_Cloning/Projects/F30_mod/WIP/F30_final_assembly/`:
`F30_D{1,2,3,4}-(8xCsy4)-f1.gb`, `-(16xCsy4)-f{1,2}.gb`, `-(24xCsy4)-f{1,2,3}.gb`.
**D1 = mod1, D2 = mod2, D3 = mod3, D4 = mod4**, which `make_inputs.py` confirms
rather than assumes: the per-unit Csy4 site-spacing fingerprint of each D
fragment matches its pART map exactly (`41 57 41 57 41 57 41` for D1/mod1,
`41 55 41 57 41 53 41` for D2/mod2, `38 75 38 75 38 75 38` for D3/mod3,
`44 45 47 43 43 46 45` for D4/mod4).

⚠️ **A 24x construct is three DISTINCT fragments, not one unit tripled.** f1, f2
and f3 are different sequences - `design_fragments_seq.py` designed each in the
locked context of the previous ones, to break the repeat symmetry that makes
tandem Csy4 arrays misregister. A construct at k copies is therefore the first k
fragments assembled **head to tail**: 8x = f1, 16x = f1+f2, 24x = f1+f2+f3.
`assemble()` asserts that the 8x / 16x / 24x file sets really are the same
fragments (they are, by md5) and that each prefix carries exactly 8k Csy4 sites.

| arch | fragments (nt) | assembled 1x / 2x / 3x |
|---|---|---|
| mod1 | 396, 396, 396 | 396 / 792 / 1188 |
| mod2 | 409, 408, 409 | 409 / 817 / 1226 |
| mod3 | 442, 442, 442 | 442 / 884 / 1326 |
| mod4 | 402, 401, 402 | 402 / 803 / 1205 |

**The design targets, vendored into `inputs/dseries_targets/`** so the figure
rebuilds without hunting them down again. Two kinds, because the two pairs of
architectures were designed differently:

* **`whole_3x` (mod2, mod4).** These converged as one molecule, so the 3x design
  IS the target and the ordered fragments are its three synthesis cuts.
  mod2 = `D2_branched_lean_3x`, 1226 nt / 396 bp - exactly the assembled length.
  mod4 = `F30_mod4_pos`, 1205 nt / 336 bp - also exact, because **mod4's
  geometry never changed between the two design generations**; only its spacers
  were re-rolled, which is what `inputs/F30_mod4_f1_mask280.txt` is for.
* **`flanked_unit` (mod1, mod3).** These never converged whole, so each ordered
  fragment is **5 nt 5' flank + the 1x design + 5 nt 3' flank** and the 3x
  target is that unit three times. The GenBank annotates both flanks
  (`misc_feature 1..5` and `392..396` for D1, `438..442` for D3), the 1x designs
  are `D1_tandem_lean` 386 nt / 120 bp and `D3_splayed_tandem` 432 nt / 120 bp,
  and 5 + 386 + 5 = 396 and 5 + 432 + 5 = 442 exactly. The +5 offset is
  **confirmed, not assumed**: every one of the eight Csy4 sites in each fragment
  sits exactly 5 nt later than in the 1x design, and `make_inputs.py` asserts
  the whole 24-site registration between target and ordered sequence.

Each cell is then `RNA.fold` of the assembled construct against the matching
prefix of its target, with `RNA.bp_distance`; the denominator is twice the
target's pair count, the maximum the distance can reach. Table:
`inputs/supp_scaffold_design_ordered.csv`.

| | 1x (8 sites) | 2x (16 sites) | 3x (24 sites) |
|---|---|---|---|
| mod1 — tandem, compact | **2** / 240 | **71** / 480 | **53** / 720 |
| mod2 — branched | **0** / 264 | **0** / 528 | **0** / 792 |
| mod3 — tandem, splayed | **11** / 240 | **99** / 480 | **165** / 720 |
| mod4 — ring | **0** / 224 | **0** / 448 | **0** / 672 |

### The panel has to explain itself — rebuilt 2026-11-20

Leo: *"Also explain what panel d does."* He is right that it did not. A reader
met a 4 × 3 grid of boxes each holding a bare number and a bare `2 / 240`, with
the only key two grey words — "Csy4 sites", "bp off target" — floating beside
the first row, and neither of those said what the **denominator** was or what
the grid was indexed by. STYLE_NOTES section 0 sends explanatory notes to the
caption, but a reader meets the panel first, and an **axis or row heading is
not a note**. So the explanation went on as headings and keys, never as prose:

* **Both grid axes are named.** `architecture` sits over the row labels and
  `fragments assembled, 8 Csy4 sites each` over the columns, so what the twelve
  cells are indexed by is stated rather than inferred from `mod1` and `1×`.
* **Both cell lines are keyed in numerator / denominator form**, beside the
  first row: `Csy4 hairpins formed / designed` and `bp off target / 2 × target
  bp`. The second one is the thing that was actually missing — the denominator
  is **twice the target's pair count**, the maximum a `RNA.bp_distance` can
  reach, and nothing on the panel had said so.
* **The top line now draws `hairpins_ok / n_csy4`** — 8/8, 16/16, 24/24 — where
  it drew the site count alone. ⭐ This is the point of the rebuild: it puts the
  figure's actual result **on the panel as data** instead of leaving it to the
  caption. Every Csy4 hairpin folds in all twelve constructs, so where mod1 and
  mod3 depart from target it is the **scaffold between the sites and never a
  site**, which is the part that matters for a Csy4 display scaffold. `main()`
  asserts `hairpins_ok == n_csy4` over all twelve, so a regenerated table
  cannot leave the key claiming something the data stopped saying.

`cw` went **14.0 → 13.0 mm** to pay for the wider right-hand key: the gutter
needs 25.1 mm for `Csy4 hairpins formed / designed` at 5.0 pt (`style.text_w`,
and `panel_d` asserts each key against `W`), and at `cw` 14.0 only 24.2 mm was
left. The cells still hold `165 / 720` comfortably at 12 mm inner width.

**What the panel now says.** The two architectures that converged as whole
molecules - mod2 branched and mod4 the ring - hit their target **exactly at
every copy number**. The two tandem arrays drift as copies are added, and mod3
(the splayed array) drifts worst, reaching 165 bp off at three copies. That is the
misregistration of tandem repeats of identical Csy4 sites showing up as a
number, and it is the reason the fragments were designed heterogeneous in the
first place.

**And in all twelve, every Csy4 hairpin folds.** `hairpins_ok` counts the Csy4
stem-loops the MFE actually forms: 8 / 8, 16 / 16, 24 / 24 in every cell. Where
mod1 and mod3 depart from target it is the scaffold between the sites, never a
site itself - which is the part that matters for a Csy4 display scaffold.

## Panel (c), rebuilt from the ORDERED fragments (2026-11-20) — the exchange HAPPENED

Leo: *"Can you recreate the c panel (dot brace based) with the real ordered
designs? If yes we can exchange them."* The exchange was made, and then the
panel gained a second row (*"we also have the actual predictions of the
sequences by viennafold. Can you add them to figure."*). **Verified on disk
2026-11-20: `panel_c` writes one file.**

| file | what it draws |
|---|---|
| `panel_c_ordered_fragments.*` | the **ordered 1x fragments**, 396 / 409 / 442 / 402 nt — the same objects panel (d)'s 1x column scores — each drawn twice, design target above and RNAfold prediction below, on a 178 x 88 mm canvas |
| ~~`panel_c_target_structures.*`~~ | the **first generation**, no longer emitted; the drawing is in `archive/` |

Both rows come from the one `RNA.fold` call whose distance panel (d) reports,
and that distance is asserted against `supp_scaffold_design_ordered.csv` when
the layout is built, so (c) and (d) cannot disagree about what a sequence folds
into. All eight folds are at **the same `SCALE`** as the rest of the figure,
which is why the panel is 88 mm tall rather than redrawn smaller: the
comparison being asked for is structural, and a target and a prediction at
different sizes cannot be read against each other. mod2 and mod4 are at
distance 0, so their two rows are the SAME drawing — drawn twice anyway,
because two empty cells would read as missing data rather than as agreement.

### 🚨 Area is not magnitude — the caption has to say so

**The drawn difference between a target and its prediction is NOT proportional
to the structural difference**, and on this panel the gap is large enough to
mislead. naview places an exterior loop from its pairing context, so a handful
of pairs opening or closing near a multiloop reorganises the entire drawing.
mod3's two rows look dramatically different — an open ring above, a tight core
with splayed arms below — while mod1's look near-identical; the measured
distances are **11 / 240 (4.6%) against 2 / 240**. A reader given the drawings
alone will take mod3 as catastrophically wrong.

The panel shows **where** a design departs from its target; the `bp off` figure
printed under each column is what says **how much**. The caption draft above
separates the two readings explicitly, and it must keep doing so — without that
sentence the figure argues something stronger than the data.

This is also what sharpens the panel's point. mod3's eleven pairs are in the
**linker between modules**, not at a site: all eight Csy4 hairpins form in the
prediction (asserted in `main()`, which requires 8 Csy4 runs in every `target`
and `mfe` layout), and panel (d) asserts the same across all twelve constructs.
The honest reading is **"the linker region reorganises, the sites do not"** —
which the drawing supports only once the reader is told not to read area as
magnitude.
The layouts live in their own table,
`inputs/supp_scaffold_design_layouts_ordered.csv`, so the table panel (a) uses
is untouched.

**The flanks are drawn.** mod1 and mod3 are `flanked_unit`, so their ordered
fragment is 5 nt + design + 5 nt; drawing the fragment *as ordered* keeps one
definition of "1x" across all four architectures and across (c) and (d), and
makes the length printed under each fold the number (d) scores. They show as two
short unpaired tails. mod2 and mod4 have no flanks — ⚠️ **their 1x is the first
synthesis CUT of a 3x design that converged whole**, not an independently
converged 1x design, so for those two the panel shows a piece of a molecule
rather than a molecule that was optimised at that size. The caption says so.

### ⚠️ The architecture names no longer match what the ordered folds show

Measured off both layouts — the base pairs in the helix carrying each Csy4 site:

| | first generation | **ordered** |
|---|--:|--:|
| mod1 (now "tandem, compact") | 8.5 bp | **5.0 bp** |
| mod2 "branched" | 8.5 bp | **5.0 bp** |
| mod3 (now "tandem, splayed") | 5.5 bp | **5.0 bp** |
| mod4 "ring" | 5.1 bp | **5.1 bp** |

**In the ordered designs every arm is 5 bp. The long-arm / short-arm
distinction does not exist any more** — the low-structure redesign removed it on
purpose, to meet a budget of at most 30 bp in any 120 nt window
(`csy4_redesigns_lowstructure.txt`: *"Tandem 4 modules, stems 5/5/5"*), and
shortening every arm is how it got there.

What separates mod1 from mod3 now is the **linker between modules**, not the
arm: per-site spacings `41 57 41 57 …` for mod1 against `38 75 38 75 …` for
mod3, so mod3's four modules sit much further apart. In the drawing that is a
bigger exterior loop — core radius 10.5 mm against mod1's 6.7 mm — and it makes
**mod3 the LARGER fold (442 nt, reach 18.5 mm) where mod1 is the smaller
(396 nt, 13.3 mm), the reverse of the first generation.** A reader who takes
"long arms" and "short arms" at face value would read the new panel backwards.
The design files' own words for the two are `D1_tandem_lean` and
`D3_splayed_tandem`.

**Leo's answer (2026-11-20): rename them.** Done - see the naming section
above. The new words are the design files' own, and COLORS.md rule 9 carries the
same correction.

**The topology is unharmed, which is what the panel is for.** Counting
top-level helices and the Csy4 sites under each: mod1 and mod3 are `[2,2,2,2]`
in both generations — four F30 modules of two sites each — while mod2 and mod4
are a single helix carrying all eight, drawn as a branched trunk and as a ring.
So the three architecture *classes* still read apart at a glance; it is only
mod1-against-mod3 that has lost its visual cue.

## ⚠️ Why the two sources disagreed, and what is superseded

`inputs/F30_mod{1,2,3,4}_pos.txt` on Euler, and the
`inputs/csy4_dotbracket/F30_mod*_pos.dbn` copies in this repo, are the
**first-generation targets of 20-JUL-2026**, superseded two days later. Their
designs sit in the sibling folder `.../F30_mod/WIP/F30_mod(8xCsy4)/` and **none
of them occurs in any ordered map**. The per-unit site spacings say so directly:
first generation `48 68 ...` for mod1 against the ordered `41 57 ...`.

⚠️ **The trap:** the Twist-style `(seq_XXXXXXXXXX)` IDs and the
`~/Downloads/D{1,2,3,4}_8xCsy4.txt` files (which say "ORIGINAL design", 458 /
521 / 410 / 402 nt) sit on the **superseded** set. The ordered fragments are
396 / 409 / 442 / 402.

**mod4 is the exception that explains the confusion.** It was never a different
architecture: between the generations its eight Csy4 site starts are identical
(29, 73, 118, 165, 208, 251, 297, 342) and all seven spacers keep their length,
every one simply re-rolled. mod1-mod3 got new architectures and their geometry
moved; mod4 kept its geometry, which is why it alone matched in both sources.

**This replaces the VERIFY this README carried on 2026-11-19** ("the ordered
pART constructs are not the `.dbn` designs for mod1-mod3"). That observation was
right about the mismatch and wrong about its cause: the pART maps are not a
different design, they are the *second* generation, and the `.dbn` files are the
superseded first. The `NEED:` for "the dot-bracket targets the pART inserts were
designed against" is **answered** - they are `inputs/dseries_targets/`.

The two sources are the same design family: the D fragments and the pART maps
share the per-unit spacing fingerprint exactly, and the 1x ordered fragment is
literally present in the mod2 and mod4 1x maps. At higher copy number the
plasmid carries Golden Gate assembly junctions, so a plain concatenation is not
a substring of it - which is expected, not a discrepancy.

## What else the files corrected

* **pART254 is the mod4 row's 3-copy construct, not a 2-copy one.** Its file
  name reads `(mod4)x2_(mod3)`, which is composition; it carries 24 sites. The
  per-site spacing fingerprint along the insert is
  `44 45 47 43 43 46 45 | 78 | 38 75 38 75 38 75 38 | 75 | 44 45 47 43 43 46 45`:
  mod4, then **mod3 in the middle position**, then mod4. The
  `LV1_15_SplicingRatio_ByFamily_thesis.svg` footnote ("a mod3 unit in its
  middle position") is right and the file name is the misleading one.
* **`pART249_notcloned.dna` is mod3 1×.** Its name carries no architecture;
  `make_inputs.py` assigns it from the same spacing fingerprint, which matches
  mod3 and nothing else.
* The four `.dbn` files carry **24** Csy4 sites each, so they are the 3-copy
  forms; one repeat is the 8-site construct. Unit lengths 463 / 520 / 415 /
  401 nt.
* Panel (c)'s folds reproduce the topology the targets encode: mod1 and mod3
  open into a rosette of four two-armed units around one exterior loop, mod2
  branches, mod4 closes into a ring with its eight hairpins on the outside.

## Deliberate

* Panel (c)'s four folds are at **one scale**, set by mod2, the widest. mod1's
  compact core against mod3's splayed one is therefore something the reader
  measures off the panel. Panel (a) is at **twice** that scale, stated in the
  caption: at 138 nt against their 400-odd it would otherwise be a smudge.
* **No per-nucleotide dots.** At 400 nt in 40 mm they are a texture, not
  information. The fold is three strokes: backbone, neutral base-pair rungs,
  and the 28 nt Csy4 repeats drawn over the backbone in the warm colour.
* Panel (c) still carries its **per-repeat lengths**, and panel (a) its own, so
  the common scale can be checked against a number.
* **(b) and (c) do not disagree about mod1 and mod3.** naview draws an exterior
  loop as a circle, so an array of four independent units comes out as a
  rosette; (b) draws the same four units strung along the transcript. The
  pairing is identical, the layout convention is not, and the caption calls (c)
  a layout rather than a shape.
* Panel (d) shows the **design** space only. Which assemblies worked is
  `scaffold_selection` panel (a), so nothing here implies pART250 or pART251
  were obtained.
* Under palette v0.9 colour carries **nothing** here: all four architectures
  are one blue (theme 600 `#0678a6` since 2026-11-20, `MAIN` before it), so the
  drawn fold and the printed name are the whole encoding. That is why
  `make_figure.py` asserts every architecture has a name.
* Glyphs come from `common/scaffold_glyphs.py`, the module Fig 1f and Fig 4a
  share, so a target scaffold is the same drawing in all three figures. Those
  glyphs draw a helix as one line, not the two rails Fig 1f's own README asks
  for; changing that is a change to the shared module, not to this figure.

## Caption draft

**Four scaffold architectures, each at one, two or three copies.** Secondary
structures are ViennaRNA naview layouts; base pairs grey, the 28 nt Csy4
repeats orange. (a) One unit of the F30 array the designs start from, as the
eight-unit array of pLV009 folds (138 nt, two Csy4 sites on the two arms of a
three-way junction); drawn in F30's own blue and at twice the scale of (c).
(b) The four designed architectures - mod1 a compact tandem array, mod2
branched, mod3 a splayed tandem array, mod4 a ring - with the length of the
ordered 1x fragment; the four
share one colour because they are four variants of one module, and
each is named on the panel. (c) The 1× fragment of each architecture as ordered
— the same molecules scored in (d) — drawn twice: its dot-bracket design target
above and the ViennaRNA prediction of the ordered sequence below, all eight
folds at one scale; mod1 and mod3 carry their 5 nt synthesis flanks as two
short unpaired tails, and the mod2 and mod4 fragments are synthesis cuts of a
3× design that converged whole rather than designs optimised at this size, so
their two rows are the same drawing. **The pair of drawings shows *where* a
design departs from its target; the `bp off` figure beside it is what says *how
much*, and the two must not be read as one.** A naview layout places an
exterior loop from its pairing context, so a handful of pairs opening or
closing near a multiloop reorganises the whole drawing: mod3's prediction looks
far more different from its target than mod1's does, but it is 11 base pairs
out of 240 — 4.6% — against mod1's 2. Those eleven sit in the scaffold between
the sites, and all eight Csy4 hairpins form in the prediction, as they do in
every construct in (d). (d) The design space: architecture against the
number of fragments assembled. Each cell gives, above, the Csy4 hairpins the
ViennaRNA MFE fold of that construct actually forms out of the sites it was
designed with (1× = 8, 2× = 16, 3× = 24) and, below, the base pairs by which
that fold departs from the design target, out of twice the target's pair count,
which is the most the distance can reach. Each construct is the first one, two or
three of its three distinct synthesis fragments assembled head to tail. mod2 and
mod4, which converged as whole molecules, meet their target exactly at every
copy number; the two tandem arrays drift as copies are added, mod3 worst, which
is the misregistration of repeated Csy4 sites. Every Csy4 hairpin folds in all
twelve (8/8, 16/16, 24/24) — the departures are in the scaffold between the
sites, never at a site. A dashed cell is a map that was designed but never
cloned; pART254, the mod4 row's three-copy construct, carries a mod3 repeat in
its middle position. n = 12 constructs, 4 architectures.

## Open

* Both earlier `NEED:`s are **closed**. There is no 1× ROAD run for mod1-mod3
  and the panel no longer wants one; and the design targets for the ordered
  inserts are found, vendored in `inputs/dseries_targets/`.
* F30's pale blue in panel (a) is **decided, not open** - see the Colour section
  above. Leo keeps it so F30 matches the main text; the 1.38:1 contrast is an
  accepted exception and `THEME["300"]` / `["400"]` were rejected. Do not
  re-raise it.
* ✅ **(c) is CLOSED, and so is (b).** Verified on disk 2026-11-20: `panel_c`
  emits one file, `panel_c_ordered_fragments`, and it is now the two-row panel
  (design target above, RNAfold prediction below) rather than the one-row
  target drawing this README used to describe; `panel_c_target_structures` is
  no longer written and the first-generation drawing is in `archive/`.
  **(b) followed**: it takes its four lengths from
  `supp_scaffold_design_ordered.csv` (396 / 409 / 442 / 402 nt, confirmed in
  the render) and parses mod2's branching from the ordered target. The "two
  panels, Leo picks" and "(b) is still first-generation" bullets that stood
  here are **stale** and were corrected rather than re-done. The rename landed
  too: `NAME` is `tandem, compact` / `branched` / `tandem, splayed` / `ring`
  and is printed on (b), (c) and (d).
* ⚠️ **(b) and (c) still draw the same glyph for mod1 and mod3.** Both arms are
  5 bp in the ordered designs, so there is nothing true to draw differently;
  what separates them - how far apart the modules sit - has no parameter in
  `common/scaffold_glyphs.py`, the module Fig 1f and Fig 4a share. The names
  carry it. Adding a spacing parameter is a change to the shared module, so it
  stays flagged rather than made here.
* The engine question is **closed**. The table was first computed with
  ViennaRNA 2.7.0 while Euler was unreachable, then regenerated on Euler with
  **2.7.2** through `make_inputs.ordered_folds()` as the documented command
  does. All twelve rows are **byte-identical** between the two versions, MFE
  included, so the numbers do not depend on which of the two folded them. The
  `engine` column records 2.7.2, the Euler run.
* The old `VERIFY:` on panel (a) is **closed**: the panel no longer draws the
  junction schematically. It folds the ordered pLV009 array and draws the
  domain that fold puts the two sites in, so two-sites-per-unit now rests on a
  structure, not only on the 16 sites counted over eight units.
* `VERIFY:` panel (a) is the **MFE of the array**, panel (c) the **design
  targets**. The two are not the same kind of object, and the caption says so.
  If the F30 unit should instead be shown as its own isolated 116 nt fold, note
  that the isolated period folds with one Csy4 hairpin dangling outside the
  junction — the array fold is the one that reproduces the F30 three-way
  junction.
