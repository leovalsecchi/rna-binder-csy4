# Fig1f_scaffold_splicing — target-presentation optimisation workflow

Output: `output/Fig1f_scaffold_splicing/scaffold_workflow.{svg,pdf,png}`.
Rebuild with `python code/Fig1f_scaffold_splicing/make_figure.py`.
Placed as **panel d** of Figure 1 (the letters are the page order, not the
folder names), at 174 x 35.5 mm.

## What it shows

One wide strip, seven steps, in the order Methods 2.12.2-2.12.3 performs them:
module design -> sequence design (PyFurnace 1.2.2, against AlphaFold3 and
ViennaRNA) -> structure re-prediction -> cloning into the ccdB entry plasmid
(PaqCI, pTWIST-Kan, *E. coli* Top10) -> the landing-pad-compatible backbone
(BbsI with the ccdB-weak part pART177, pTWIST-Amp, DB3.1) -> Bxb1 integration
into the HEK293-H C8 Spartan landing pad (integrase in trans, pMLS596) ->
splicing readout (RNeasy Plus, SuperScript IV, KAPA HiFi ~850 bp, ONT;
ratio = removed / (intact + removed)).

The revise arrow returns to the DESIGN step, not to the sketch: Methods frames
design and re-prediction as one iteration, "iterating until each fragment
satisfied the synthesis constraints of the supplier".

House workflow idiom (as `readout_workflow` sets out): line glyphs, one short
label each, thin chevrons, **no boxes**. The tinted boxes the earlier version
used are gone with it.

## History

- **2026-11-17.** Split: this panel used to carry the splicing-ratio bars on its
  right. Those moved into `Fig1e_assay_bridge`, which removed a real duplication
  — its two ONT bars were the same two samples (pART242 0/267, pART158 117/300,
  one LV1.27.1 run) this panel drew again.
- **2026-11-17, later.** Rebuilt from three nodes at 100 x 62 mm to seven at
  178 x 38, on Leo's "more horizontal less vertical" and "we also had more
  steps". The assembly gives this panel a fixed height and derives its width
  from the ink aspect, so **trimming ink height is what widens it on the page**.

## The four architectures

All four are redrawn from Leo's own schematics at
`~/rna_binder_plots/thesis_figures/RNA_Scaff_Mod/schematics/schematics.py`
(which asserts eight hairpin rings per design), with the topology checked
against the dot-bracket targets.

**Which dot-bracket files are current.** `~/Downloads/csy4_dotbracket/F30_mod{1..4}_pos.dbn`,
not `~/rna_binder_plots/thesis_figures/targets/D*.db`. Two reasons:
`scaffold_target_structures.md` says so outright ("The current targets are the
`.dbn` files"), and the D-numbering is a DIFFERENT scheme from the mod-numbering
— D1 is the pART158 family, not mod1, and only D4 matches a construct map
(pART252) exactly. The `.dbn` files carry the mod names the Methods itself uses
(`F30_mod1_variablespac`, `F30_mod2_stemwith_2x3`, `F30_mod4_variablespac`).
The structures agree where they overlap: `D1.db` is byte-identical to the first
458 nt of `F30_mod1_pos.dbn`. The `.dbn` files are the 3-copy, 24-site forms;
the panel draws one repeat, which is the 8-site construct.

| | topology in the .dbn | drawn |
|---|---|---|
| mod1 | 4 x 3WJ array, branch stems 9 and 8 bp | tall-armed linear array |
| mod2 | balanced binary tree, trunk -> 2 -> 4 -> 8 | branched |
| mod3 | the SAME 3WJ array as mod1, stems 6 and 5 bp | short-armed linear array |
| mod4 | long closing stem whose loop carries all 8 hairpins | ring |

mod3 never physically existed (ART249 undelivered, pART250/251 failed to
assemble) but is still drawn: step 1 is what was DESIGNED, and four
architectures were designed. What survived assembly belongs in the caption.

## Two drawing rules that are not cosmetic

**One scale for all four.** Fitting each geometry to its own box normalises away
arm length — which is the ONLY thing separating mod1 from mod3 — and drew them
as near-identical. The largest extent sets k and the rest follow.

**Every helix is two rails.** Single lines read as unpaired strand; these are
duplexes. `DUPLEX = 9.0` source units: at the common scale (~0.067 mm/unit)
7.0 left only 0.77 pt of white between 0.55 pt strokes, which fuses in print.
The Csy4 loop stays a ring, wider than the duplex, as a loop closing a stem is.

## Open points

- The two geometries in step 1 are not individually labelled; the caption has to
  name them. There is no room for labels at 24 mm per step.
- Delta O1 of the Methods contradicts §2.12.3 on the entry vector for the
  scaffold layer (pMLS633 against pMLS543). The panel names neither, so it is
  unaffected, but the caption should not name one until that is resolved.
- No selection or sorting of integrated cells is described anywhere in the
  Methods; the strip shows integration and then readout, with nothing between.

## The four architectures are drawn 2 x 2 (2026-11-17)

Across then down, in NUMERIC order: mod1 mod2 / mod3 mod4. Two reasons it beats
the single row of four it replaced:

- **It suits the shapes.** mod4 is square and mod2 nearly so. In a wide, short
  cell those two were what bound the shared scale, so all four were drawn
  smaller than they needed to be. A squarer cell lets every one sit larger.
- **It gives the strip 18 mm back.** Step 1 went from 54 mm wide to 36, and the
  other six steps took the difference; they were the tightest part of the row.

Numeric order also removes a real ambiguity: the row grouped the two linear
designs (mod1, mod3) together, which put the branched design third and made it
read as mod3. It is mod2.

## The four shapes are verified against the dot-brackets (2026-11-17)

**The Methods never states any scaffold's shape.** It gives the process only -
a sketch converted to dot-bracket, four architectures, PyFurnace against
AlphaFold3 and ViennaRNA. So Methods cannot confirm these glyphs; that is a gap
in Methods, not a tick.

The dot-bracket targets can, and do. Parsing one repeat of each
`~/Downloads/csy4_dotbracket/F30_mod*_pos.dbn` and walking the structure:

| | the dot-bracket | drawn as |
|---|---|---|
| mod1 | 8 hairpins, **4 top-level helices**, each opening one 2-way junction | 4 Y units on the strand |
| mod3 | the SAME topology, unit 415 nt against mod1's 463 | the same array, shorter arms |
| mod2 | 8 hairpins, **1 top-level helix**, **7 two-way junctions** | trunk -> 2 -> 4 -> 8 |
| mod4 | 8 hairpins, 1 top-level helix, **one junction opening 8 branches** | a ring carrying 8 arms |

Seven binary junctions producing eight leaves is exactly a trunk splitting
twice, so mod2's tree is right; mod4 being a SINGLE 8-way junction is why it is
a ring and not a tree.

**Caveat.** mod2's seven junctions sit at uneven sequence depths (7, 13, 14, 19,
20, 32, 32), so the molecule is not the balanced tree the glyph draws. The glyph
is right about TOPOLOGY - what branches into what - and schematic about
proportion. If the chapter ever states the shapes, write that paragraph from the
dot-brackets, which are the ground truth here.
