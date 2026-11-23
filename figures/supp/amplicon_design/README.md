# amplicon_design — amplicon design and primer test

Two bare panels in `output/amplicon_design/`, 85 mm wide (a **26.5** mm tall, b 61).

```bash
~/.venvs/binder_map/bin/python3 code/supp/amplicon_design/make_figure.py
```

No derived table. The lengths are in a worksheet, the module spans in the
construct files, the gel is an image. All are named below.

## Sources

| source | what it gives |
|---|---|
| `LV_PCR/LV1.4.1.xlsx` (2026) | the three reactions: one forward primer, three reverse, three expected sizes, three lanes |
| `LV_PCR/LV1.27.xlsx` | the same pair PR10148 + PR10137, expected ~850 bp, on the same template |
| `LV_Sequencing/LV1.27.1_ONT/MLS12_18_1_AllSamples.fa` | that amplicon, sequenced: **865 bp** |
| `LV_Cloning/pLV_NNN/pLV018/pLV018.gb` | the ORF the amplicon runs into: Kozak, Csy4-H29A, linker, MgIRE1 |
| `LV_Cloning/pLV_NNN/pLV001/pLV001_part158.gb` | the AT cassette, the F30 scaffold, and `PR249_BindingSite` / `PR250_BindingSite` |
| `working/Fig1d_integrated_construct/Fig1d_integrated_construct.svg` | **the glyphs**: main figure 1's construct map, every proportion and colour in panel (a) |
| `LV_Gels/gel_2026-07-15_primer-test_inverted.png` | the gel |
| `LV_Protocols/gel_protocol.md` | 1.2% agarose, 1x TAE, 135 V, "1 kb DNA ladder" |
| journal LV01.04, 14.07.26 (Notion) | the same experiment, same pairs, on `MLS12.18.1.1_cDNA` |

## Verified

`LV1.4.1.xlsx` is journal LV01.04. It records one master mix: KAPA HiFi,
50 µL, template `MLS12.18.1.1_cDNA`, forward primer **PR10148**, and:

| reaction | reverse primer | expected | lane |
|---|---|---|---|
| LV1.4.1 | **PR10137** | ~850 bp | 3 |
| LV1.4.2 | **PR249** | ~2,650 bp | 5 |
| LV1.4.3 | **PR250** | ~2,650 bp | 7 |

Cycling: 95 °C 3 min, then 98 °C 20 s / 65 °C 15 s / 72 °C 3 min, **24
cycles**, 72 °C 5 min. Methods 2.12.5 gives the same pair at 0.5 µM each over
24 cycles, so the two records agree.

**Both PR249 and PR250 are the target-capturing options.** They are the two
~2,650 bp designs that would have spanned sensor and target in one amplicon;
PR10137 is the ~850 bp design that was adopted.

**The gel alone does not show "too long to amplify".** It shows two empty
lanes. The lengths in panel (a) carry that claim, which is why both designs
are drawn and not only the one that worked.

The gel date, `2026-07-15`, is one day after the journal's `14.07.26` for
LV01.04. Treat them as the same experiment run over two days or as a filing
date; nothing else distinguishes them.

## Where each module span comes from

Panel (a) is drawn in **amplicon coordinates**: base 1 is the 5′ end of
PR10148. Everything below is read off files, nothing is sketched.

The chain that fixes the origin: `LV1.27.xlsx` records PR10148 + PR10137,
~850 bp, on this template; `LV1.27.1_ONT/Universal_Reference_fasta/` holds the
reference that run aligned to, `MLS12_18_1_AllSamples.fa`, **865 bp**. So the
adopted amplicon exists on disk as sequence. Aligning it to `pLV018.gb`:

* amplicon **1..62** is identical to pLV018 1180..1241 — the Twist 5′ UTR.
* amplicon **63..221** is a **159 bp insert absent from pLV018**, sitting
  immediately before the Kozak. Its 3′ 30-mer is found in `pLV001_part158.gb`
  at position 130, inside `At_hairpinATGremoved_Overhangs` (44..160). That is
  the **recorder**.
* amplicon **222..865** is identical to pLV018 1240..1883, which carries
  Kozak 223..232 (ATG at 229), **Csy4-H29A** 241..801, the (GGS-GGT)x2 linker
  808..843, and the first 19 bp of **MgIRE1** from 847. That is the
  **sensor**: one ORF, binder slot fused to SMS.

Carrying pLV018's own annotation past the end of the amplicon puts the MgIRE1
stop at amplicon **2328**, so the sensor is 229..2328 (2,100 bp). The
**target** starts there, at amplicon 2329.

**Where it ends is the 2026-11-19 correction.** It used to run to ~2,650 —
"the remainder up to the reverse primers" — which is the *product's* end, not
the module's, so the block swallowed the 3′ UTR and the primer binding site
with it. That is what Leo caught: *"The reverse primers should be after so
more right than the target module."* Both primers anneal **3′ of** the target
module, and the block was drawn straight through them.

It now ends at **2580**, derived: a product ends at its reverse primer's 5′
end, so PR249's 20 bp footprint is amplicon 2631..2650, and
`pLV001_part158.gb` puts **50 bp** between the F30 array's 3′ end (3202) and
PR249's site (3253) — Twist 3′ UTR 40, X6-GTCA overhang 4, 3′ UTR 6. So the
target is **2329..2580 (252 bp)**, and **2581..2650 is bare strand**: the 3′
UTR and the primer site, drawn as Fig 1 draws its 3′ lead. That is the strand
the reverse primer marks sit on.

`PR249_BindingSite` (3253..3272) and `PR250_BindingSite` (3273..3292) are
annotated in `pLV001_part158.gb`, **3′ of** the eight-copy
`F30_three_way_junction_scaffold` array (2281..3202, each copy 110 bp) and 3′
of the Twist 3′ UTR (3203..3242). That is the file evidence that the two long
designs capture the target module; it is not their position on this template
(see Still open).

**The two sites ABUT**, so PR250 anneals exactly **20 bp 3′ of PR249** and
gives the longer product. Re-read off the file 2026-11-19. That is the one
thing about the pair this panel can state, and it is now drawn: PR250's bar
ends 20 bp further right than PR249's, which is 0.44 mm at this scale — small,
because 20 bp out of 2,650 *is* small. Both bars keep the worksheet's ~2,650
label, because the **order carries across templates and the positions do
not**: see Still open for why the donor's own coordinates imply ~3,360 bp.
`_selfcheck` asserts the 20 bp, so the panel cannot quietly stop answering it.

## The glyphs are main figure 1's

Changed **2026-11-19** (Leo: *"take the elements from main figure 1 and align
everything properly"*). The panel used `common/constructs.py`'s SBOL bars —
CDS pentagons, module names on a baseline with leader lines. Fig 1's construct
map does not look like that, and `constructs.py`'s bars are on Leo's excluded
list (2026-11-12: *"use Matteo-derived elements only, never the ones we drew
from scratch"*). So every shape, proportion and offset now comes from
`working/Fig1d_integrated_construct/Fig1d_integrated_construct.svg`, which is
what reaches the page as **Fig 1 panel e** (`code/Fig1_assembly/make_figure.py`,
`PANELS["e"]`).

Carried over verbatim, in that file's own units — one unit is one point, the
`6 / 5.5` block-to-caption ratio `Fig1d/make_figure.py:_check` asserts:

| | Fig 1 | here |
|---|---|---|
| block | rounded rect, height 6, corner radius 0.8 | same, 2.117 mm tall |
| outline | `synbio.outline` at 0.4 | same |
| strand | INK hairline at 0.4, `5′` caption | same |
| linker | **thin**: 3 of the 6 units (28.78..31.78), neutral 300 | same |
| SMS | kinase `enzyme` + **RNase lobe** `enzyme_alt`, one outline | same |
| brackets | spine at −10.0, feet at −7.4, caption −12.78 | same |
| label ticks | +4.0..+6.1, caption +11.92 | same |
| primers | below the strand, neutral 700, schematic length | same, one tier lower |

What changed in the drawing, and why:

* the fused glyph's first part is labelled **Binder**, not Csy4-H29A. It is
  Fig 1's binder *slot*; on `MLS12.18.1.1` that slot holds Csy4-H29A, which
  the caption says and the canvas does not.
* the **grey linker and the amber RNase lobe are back**. The pentagon version
  had neither, so the sensor read as two flat blocks where Fig 1 shows four.
* the **Twist 5′ UTR (1..62)** is drawn, as a neutral 200 block.
* module names are **brackets above** — Recorder / Sensor module / Target
  module, Fig 1's own three — instead of a shared baseline with leaders. That
  is the "align everything properly": a bracket states its own extent, where a
  leader line only points.
* the intron red stays **`synbio.intron` `#a8232a`**, *not* Fig 1's `#ae4e1a`.
  That file predates palette v0.5; COLORS.md rule 4 and `constructs.py` both
  moved on. `_selfcheck` asserts the fills still equal `constructs.py`'s, so a
  palette edit cannot move one without the other.
* no **ATG tick**. Fig 1 has one on an elbow leader; here ATG (229) and the
  Binder block (241) are 0.26 mm apart, so it would need the elbow and buy
  nothing (STYLE_NOTES §0).
* no **`3′` caption**. Fig 1 has one because it draws the whole mRNA. This
  panel stops at the reverse primers, and the transcript does not.

⚠️ **One boundary on this panel is not base pairs.** The RNase lobe's share of
the SMS block, `RNASE_FRAC`, is Fig 1's own proportion read off its paths
(19.904 of 51.3685, so 38.7%). No kinase/RNase boundary for MgIRE1 was read off
this template. Every *other* boundary is the alignment's. Said in the caption.

## Where the primers sit

Leo, 2026-11-19: *"Primers have to be after or before the elements."* The old
panel drew the forward primer as an arrow on the strand at base 1, and since
the recorder starts at base 63 — 1.2 mm along — the arrowhead landed **on the
recorder block**.

PR10137 anneals at ~850, which is **inside the sensor**, so "before or after"
is not something the x axis can deliver for every primer. Fig 1 settled this
already: a primer is a short neutral-700 arrow **below the strand**, position =
where it anneals, direction = extension (`Fig1d/make_pcr.py`, and
[[thesis-fig1d-panels]] records it as the web-checked convention). Its
`PRIMER_Y = 35.3` is commented *"in the lane the label ticks used"*, and
`_bare()` strips those ticks to free it.

This panel keeps the ticks, so the primer lane drops one tier further, below
their captions. The forward primer therefore starts at base 1 — before every
element — and no primer mark can touch a module bar at any x. `_selfcheck`
asserts the lane is below the captions. The arrow is **schematic in length**
(2.2 mm): a 20-nt primer is 0.44 mm here, and Fig 1 draws its primers
schematic for the same reason. ⚠️ Pass `shrinkA=0, shrinkB=0`: matplotlib
trims an annotate arrow by 2 pt at **each** end by default, which is 1.4 mm
here — it ate most of a 2.2 mm shaft and left a bare triangle, and it also
moved the drawn head away from the coordinate the assertion measures.

The three reverse primers stay on their own amplicon bars, one per row, where
each bar's 3′ end already is its annealing site on the shared axis. Their mark
is the **bar's own terminus** — a leftward head *at* that end, the ordinary
way an amplicon is drawn under a map — not a second copy of the primer glyph.

### How big the arrowheads can be

Leo, 2026-11-19: *"these primer arrows could be a little larger. I can barely
see it. The reverse primer arrow for the lane 3 piece you made."*

The head was `mutation_scale = 4`, which matplotlib draws as 0.4 × 4 pt =
**0.56 mm long and 0.56 mm wide** — exactly the amplicon bar's own 1.6 pt
thickness. On lanes 5 and 7 that triangle at least sits on a neutral-400 bar
in near-black ink; on **lane 3 the bar is neutral-700 and the head is
`INK` #1d2428**, two darks 0.2 apart in lightness, and with no silhouette
sticking out past the bar there was nothing left to read it by. That is why
lane 3's was the one Leo could not see.

It is now `HEAD_SCALE = 8`: **1.13 mm**, twice the bar's thickness, so the head
reads by *shape* and no longer depends on the ink-against-neutral-700
contrast. Verified by zooming the render — the triangle now stands proud of
the dark bar top and bottom.

🚨 **8 is the ceiling, not a preference**, and it is a tight one. Two
assertions bound it and both are near their limit:

| bound | value |
|---|---|
| bar end past the target block | 69 bp = **1.52 mm** |
| so `REV_HEAD` must be < | **1.32 mm** (tip must clear the block by 0.2) |
| and < the bare 3′ strand | **1.54 mm** (70 bp) |
| `REV_HEAD` as drawn | **1.20 mm** |
| the head matplotlib draws | **1.13 mm** (0.4 × 8 pt) |

As rendered, PR249's tip lands **+0.32 mm** past the target module's right
edge and PR250's **+0.76 mm** (they were +0.72 and +1.16 at the old head
length). A `mutation_scale` of 9 would make the head 1.27 mm, `REV_HEAD` would
have to go to ~1.3, and the tip would clear the block by 0.02 mm — inside the
width of the rule itself. **So the heads cannot be made larger than this
without moving the target module's 3′ boundary**, which is a sequence claim
and not a drawing choice. All four marks are the one glyph at the one size:
making lane 3's bigger on its own is the other thing that is not available.

A third assertion now ties the two together: `HEAD_LEN <= REV_HEAD`.
`mutation_scale` sets the drawn triangle independently of the segment the
other two assertions measure, so without it the head could be enlarged
straight back over the target block while both of those still passed — they
measure the tip, not the triangle.

The forward primer's arrow went **2.2 → 2.8 mm** so that the bigger head
leaves a visible shaft rather than becoming all triangle; `_selfcheck` asserts
the shaft stays over 1 mm (it is 1.67). Both arrows keep `shrinkA=0,
shrinkB=0`.

🚨 **That head used to sit 2.6 mm (118 bp) inside the bar**, which put it deep
inside the target block, and it is the fault Leo found by zooming in. It is
now 1.2 mm, so the tip lands at the product's 3′ end. Two alternatives were
tried first and are recorded so they are not retried: a full primer arrow
shrunk to fit the 70 bp of bare 3′ strand is ~1.2 mm, all head and no shaft,
reading as a stray tick; and moving that arrow off the bar to its own lane
below reads as a smudge hanging under the bar end. The original glyph was
right — only its position was wrong, and the block's end with it.

**Both lanes were always on one x axis.** The blocks and the bars go through
the same `bp()` with the same `x0` and scale, and `_axis()` spans the longest
product drawn. The misplacement was never a coordinate-space problem; it was
a shaft length and a block that ran past its module.

## How panel (b) is calibrated and aligned

**The ladder is GeneRuler 1 kb DNA Ladder (Thermo).** The protocol says only
"1 kb DNA ladder", so this is identified from the lane: twelve bands, eight
bunched at the well end, then 2000/1500/1000/750 pulling apart — and the three
darkest *within their own neighbourhoods* are 1000, 3000 and 6000 (darkness
0.485, 0.316, 0.278 against 0.20–0.29 for their neighbours), which are exactly
GeneRuler 1 kb's three reference bands. Its 500 and 250 bands are not visible
on this run.

Sizes on the panel are interpolated **linearly in log bp between adjacent
ladder bands**, measured as darkness centroids in lane 1. No mobility model is
fitted; the axis is the ladder and nothing else.

**Alignment.** Well centres were re-measured off the raw image and fit
`x = 239.7 + 109.39·(lane − 1)` to within 4 px over lanes 1–7. The previous
version used `215 + 111·(lane − 1)`, which put every label about 20 px left of
its lane. Lane 3's product band sits at x = 458.7 against the fit's 458.5. The
well row is level (−0.016°, −0.3 px over 1000), so the crop is **not** rotated.
The panel is now emitted with `tight=False` at its authored 85 mm: with
`tight` the on-canvas note was wider than the gel, so the bounding box grew
around the text and the gel sat left of centre in a panel that was no longer
85 mm. That note is deleted — it was caption material (STYLE_NOTES §0).

**What the marks mean.** Short grey ticks: every ladder band. Longer ink
ticks, labelled 850 and 2,650: the expected products. A **dashed box** inside
lanes 3, 5 and 7: the region where that lane's expected product would have
run, captioned once on the canvas as **"expected band position"**.

**One label for the set** (added 2026-11-19, Leo: *"maybe you can label the
right ones like expected location of a band"*). It sits on lane 3's box — the
one box with a band under it, so it shows what the mark means — left-aligned
to that box's own left edge and 0.5 mm above its top rule, which labels it by
position and needs no leader. The run to its right is empty lane-4 gel. One
label, not three: the other two boxes are the same mark. `_selfcheck` measures
the string with `style.text_w` (19.75 mm) against the gel that is free to its
right (20.91 mm) and refuses if it would reach lane 5's box.

**Boxes, not rules** (changed 2026-11-19, Leo: *"could you put a box instead of
the lines? box doesnt have to be precise and can be bigger than the exact
location of the band"*). A hairline claims a precision this panel does not
have — ~850 and ~2,650 are the worksheet's nominal figures, and the one product
that did appear reads **812**. A box says "a product of about this size would
appear in this region", which is what the evidence supports.

The box spans the expected size **±25% in bp** (`BAND_TOL`), symmetric in log
bp because that is how the size axis interpolates: 680–1,062 and 2,120–3,312.
The window is a drawing choice, stated in the caption, not a measurement.
`_selfcheck` asserts lane 3's box **contains the band that actually ran** — if
the window were ever narrowed below the panel's own 812-against-850, the mark
would be claiming more than the panel can show, and the build stops. The
ladder ticks are unchanged.

### The box did not go round the band — three separate causes

Leo, 2026-11-19: *"make them literally a bit taller and wider. As you can see
with the leftmost and lower one, it's not even going around the band fully."*
The leftmost and lower box is lane 3's, the only one with a band in it. It was
cutting that band on **two** sides, for **three** reasons:

**1. 🚨 `bp_to_y` clamped, so the box could not be made taller at all.** It
interpolated the ladder with `np.interp`, which **clamps** outside its range.
The ladder's smallest visible band is **750 bp at y 194.2**; the box's lower
edge is at `850 / BAND_TOL`, which is below 750 for any tolerance past 1.133.
So **every** tolerance from 1.14 upward drew the same lower edge, 194.2, two
pixels inside the band's own lower edge at 196 — and raising `BAND_TOL`, the
obvious fix, would have changed nothing. `bp_to_y` now extrapolates on the
terminal interval's slope instead. It is 0.04 of a decade past the last band
(680 against 750) and is used only to place a drawn box edge. **`y_to_bp`
deliberately keeps the clamp**: reading a *size* off this gel outside the
ladder is not something the ladder supports.

**2. The box was narrower than the band.** 76 px against a band that measures
**90 px** (`BAND3_BOX`, the 15%-of-depth contour against a local background —
the first contour that is stable: x is 409..499 at 15, 20 and 25% and 410..498
at 35%, while below 15% it runs away into the lane's smear, 127 px at 10% and
160 at 5%). The box cut 11 px into the band's left end.

**3. Widening it symmetrically ran into the other wall.** Lane 7's centre is
at 896 and the crop **ends at 948**, because lane 8's well rim starts there
and the panel shows lanes 1–7. A 106 px box was tried first: its right edge
landed at 949 and the render showed it hanging off the edge of the gel into
the panel margin. The crop cannot grow — re-measured here, the eighth well's
bright interior starts at ~950.

So the width is **100 px** and that is a ceiling in both directions. A 100 px
box centred on the fitted well centre satisfies both bounds arithmetically —
and the render *still* showed the rule sitting on the band, because its left
clearance is 0.5 px = 0.05 mm while the dashed rule is 0.5 pt = **0.18 mm**
wide. Hence `box_dx = -4` px.

**`box_dx` is a drawing offset, not a correction to the lane model.** The band
is not centred on its lane *by extent*: it reaches 49.5 px left of the fitted
centre and 40.5 px right. Its darkness centroid is 458.7 against the fit's
458.5, so the *ink* is centred and the left **smear** is not. The well fit was
re-measured for this change and is good — well 3's bright interior centres on
**460.5** against the model's 458.5 — and 4 px is inside the ±4 px residual
the fit already carries over lanes 1–7, so the box is placed at one end of
that residual rather than against the fit. The lane labels still sit on the
fitted centre.

As drawn, lane 3's box is **9.56 × 7.68 mm** (was 7.26 × 4.56) and clears the
band by 4.5 px left, 5.5 px right, 19 px below and 36 px above — the top
margin is larger because the product ran *small*, 812 against the 850 the box
is centred on. Lane 7's box ends 0.57 mm inside the gel's right edge.

## Still open

* `NEED:` **where PR249 and PR250 anneal on this template.** Their binding
  sites are annotated in the library donors `pLV001`–`pLV004`, but those are
  not `MLS12.18.1.1`: the donor's AT cassette is 117 bp against this
  template's 159, and its binder slot holds the ccdB placeholder, not
  Csy4-H29A. Carrying the donor's spacing over gives ~3,360 bp, not the
  worksheet's ~2,650, so it is **not** carried over. Panel (a) therefore draws
  the two long products at the worksheet's length and puts the target block
  between the sensor stop and that length. No construct file for
  `MLS12.18.1.1` was found on the NAS.
  **Their relative order is a different question and it is settled**: the two
  sites abut in the donor, so the 20 bp between them is a property of the pair
  rather than of the backbone, and that offset is drawn. The absolute
  positions remain open; this `NEED:` stands.
* `VERIFY:` **the sensor's 3′ end and the target's extent.** Both rest on
  pLV018's MgIRE1 annotation continuing past the end of the amplicon. pLV018
  is the same Csy4-H29A–linker–MgIRE1 part and matches this template base for
  base over the 644 bp they share, but the stop codon itself was not read off
  `MLS12.18.1.1`. The target block (**252 bp as drawn**, 2329..2580) now ends
  70 bp before the product does, so the 3′ UTR and the primer site are bare
  strand rather than part of the module — see "Where each module span comes
  from". That 70 bp is **20 bp of primer footprint plus the donor's own 50 bp
  of 3′ UTR and cloning overhang**, and the 50 is the part that is still open:
  pLV018's Twist 3′ UTR is 64 bp, not 50, so the two constructs' 3′ ends are
  not identical and the boundary is an estimate. It is a better-founded one
  than "the remainder up to the primers", and it is the one that puts the
  primers where the files say they anneal. The block is still **not** the
  922 bp `part158` scaffold, whose length would not fit a ~2,650 bp product.
* `VERIFY:` **the adopted amplicon has three lengths on record.** ~850 bp
  (worksheet, `MLS12.18.1.1` cDNA), **865 bp** (the sequenced LV1.27.1
  reference, same template), and **509 bp** uniformly across ~48,000 pipeline
  rows on backbones pLV005–008. The first two agree to 2%; the third is a
  different template, so not a contradiction, but any sentence quoting a
  length has to name its template. The 509 bp figure was not re-derived here.
* `VERIFY:` **PR8919 against PR8916.** Both appear in the journal, at
  different entries: LV1.16 writes *"PR 8919 dillulute at 10 µM"* and LV1.42
  writes *"Used PR 8916 for both sequencings"*. A separate report finds 30
  sequencing files naming PR8916 and one other hit for PR8919. **The files
  support PR8916, and the journal supports both**, at different steps, so this
  is flagged rather than silently corrected.

## Deliberate

* The gel is flipped vertically. The raw capture has the wells at the bottom
  and migration running upward; the flip is the only change to the pixels, and
  the caption says so. The grey levels are the crop's own range, linearly
  stretched — no gamma, no contrast enhancement.
* Lanes 1–7 are cropped in and the wells are kept. The 18 mm of near-empty gel
  between the wells and the 10,000 bp band is where a >10 kb product would
  have run, and the wells are what show that lanes 5 and 7 were loaded at all.
* **Adopted against abandoned is lightness, not hue.** The amplicon bars are
  neutral 700 and neutral 400; the hues on panel (a) are molecules, and MAIN
  already means *binder* inside the sensor glyph on the same canvas.
* **The brackets are inset 0.35 mm each side.** The recorder ends at 221 and
  the sensor starts at 229 — 0.18 mm apart here — so without it the two feet
  merge into a single rule with a notch in it, which is exactly what the first
  render of this revision did. Fig 1 insets its own brackets for the same
  reason (sensor 42.2..120.4 over blocks 42..120.6). `_selfcheck` measures the
  gap and refuses below 0.6 mm.
* **The expected-band boxes are offset 4 px left of the fitted well centre**
  (`box_dx`), which is inside the ±4 px residual the well fit already carries
  and is what lets a 100 px box clear the band's left smear while lane 7's
  stays on the gel. The lane LABELS are not offset. See "The box did not go
  round the band".
* **Both long bars keep the "~2,650 bp" label** although one is drawn 20 bp
  longer. The worksheet gives one figure for both, nominal to ±50; writing
  ~2,670 on one of them would be inventing a length the record does not carry.
* The size axis labels only 10,000 / 2,000 / 1,000. Adding 3,000 would put a
  label 1.44 mm from the 2,650 one; `_selfcheck` measures that gap and refuses
  to draw.

## Checks

`_selfcheck()` runs before anything is drawn and asserts the ladder is
monotone, that its spacing is smooth in log size (no adjacent interval more
than 2x another, which is what a band assigned to the wrong size looks like),
that the calibration puts lane 3's band at 800–830 bp, that the longest size
label fits the gutter, that the blocks run 5′ to 3′ without overlapping and
reach 2,650, that every bracket spans a block, and that no two size-axis
labels land within 1.95 mm. Proved by mutation: swapping the 1500/2000
assignment → *"a smaller fragment must sit further from the well"*; moving the
1000 band 20 px → *"lane 3 product should read ~815 bp"*; adding 3,000 to the
labels → *"size-axis labels overlap: 1.44 mm"*. All three abort before a file
is written.

Added with the 2026-11-19 revision, each guarding something a reader cannot
check by eye:

* **lane 3's box must contain the 812 bp band** — the window cannot be
  narrowed below the panel's own evidence.
* **PR250 must be drawn 20 bp 3′ of PR249** — zero that and the panel silently
  stops answering the question Leo asked for.
* **adjacent brackets must clear each other by 0.6 mm** — this one caught a
  real defect, see Deliberate.
* **the primer lane must sit below the label captions** — the rule that keeps
  every primer mark off the module bars.
* **the block fills must equal `constructs.py`'s** — the shapes are
  deliberately Fig 1's and not that module's, so only the colours are shared,
  and this is what stops the two drifting apart.
* **`RNASE_FRAC` must stay in 0.30–0.45** — a reminder that it is Fig 1's
  proportion and has to be re-read from that file, not nudged.
* 🚨 **each reverse primer's arrowhead must sit past the target module's
  right edge**, measured in millimetres on the panel's own axis. This is the
  one thing that was silently wrong, so it is the one that is now measured
  rather than intended. PR10137 is exempt and says why: it anneals at ~850,
  inside the sensor, so its head belongs under a block.
  As drawn: PR249 **+0.72 mm** past the edge, PR250 **+1.16 mm**.
* **the bar's head must be shorter than the bare 3′ strand it sits on** —
  the constraint that the 2.6 mm shaft broke.

Added with the 2026-11-19 **second** revision (bigger heads, bigger boxes):

* **the drawn arrowhead must fit the segment the other two measure**
  (`HEAD_LEN <= REV_HEAD`). `mutation_scale` sizes the triangle independently
  of that segment, so without this the head could be enlarged back over the
  target block with both of the above still passing. Proved by mutation:
  `HEAD_SCALE = 10` → *"the drawn arrowhead is 1.41 mm (mutation_scale 10) and
  the segment it sits in is 1.20 mm"*.
* **the forward primer must keep a shaft** — over 1 mm of it, so the arrow
  does not become a triangle on a stick.
* 🚨 **lane 3's box must go AROUND the band**, in raw-image pixels, on all
  four sides with 2 px to spare, against `BAND3_BOX`. This is the one Leo
  reported and so the one that is now measured rather than eyeballed. Proved
  by mutation three ways: `lane_w = 76` → *"lane 3's box spans x 416..492 and
  the band spans 409..499"*; restoring `np.interp` in `bp_to_y` → *"lane 3's
  box spans y 194.2..256.9 and the band spans 196..221"*; and `box_dx = 0` →
  *"spans x 408..508"*, which is the 2 px margin catching the version that
  contained the band arithmetically but drew its rule on it.
  The width half of it carries to lanes 5 and 7, which share `lane_w`;
  **the height half cannot** — neither lane produced a band, so there is
  nothing there to contain and no honest check to write for them.
* **no box may run off the cropped gel** — this caught the 106 px box whose
  right edge landed 1 px past the crop, hanging in the panel's margin.
  `lane_w = 112` → *"lane 7's box runs off the cropped gel (836..948 against
  the crop's 182..948)"*.
* **the one on-canvas label must fit beside the box it labels**, measured with
  `style.text_w` (STYLE_NOTES: measure, never estimate). *"expected location
  of a band"*, Leo's own phrasing, is **23.5 mm** against the 20.9 mm of empty
  gel that runs to lane 5's box, and the assertion names both figures — which
  is how the shorter *"expected band position"* (19.75 mm) was chosen.

The arrowhead rule proved by mutation against the **exact state Leo found**:
putting the target block back to 2,650 → *"PR249's arrowhead is +0.82 mm from
the target module's right edge — it must sit AFTER the block, not on it"*, and
restoring the old 2.6 mm shaft → the same assertion at +1.08 mm.

All six of the earlier ones proved the same way, each aborting before a file
is written:
zeroing `REV_OFFSET["PR250"]` → *"PR250 anneals 20 bp 3′ of PR249 (3273..3292
against 3253..3272)"*; `BRK_INSET = 0` → *"Recorder and Sensor module brackets
are 0.18 mm apart"*; `BAND_TOL = 1.02` → *"lane 3's box must contain the band
that actually ran"*; raising the primer lane into the caption row → *"the
primer lane overlaps the label captions"*; hand-picking `#015d83` for the
binder → *"binder has drifted from constructs.py"*; overlapping two blocks →
*"linker overlaps the block before it"*.

## Caption draft

**Of two amplicon designs, only the short one amplified.** (a) The
recorder–sensor–target transcript drawn to scale from the 5′ end of the shared
forward primer PR10148, in the module glyphs of figure 1, with the three
candidate amplicons below it on the same axis. The binder slot holds Csy4-H29A
on this template. PR10137 gives a ~850 bp product covering the recorder and the
binder half of the sensor; PR249 and PR250 anneal 3′ of the target module and
give ~2,650 bp products spanning sensor and target in one read; both anneal in
the 3′ UTR beyond the target module, which is why their marks sit past that
block and the strand runs on after it. Their binding
sites abut in the library donor `pLV001_part158`, so PR250 anneals 20 bp 3′ of
PR249 and is drawn that way — the **order** comes from that annotation, the
**lengths** from the PCR worksheet, which gives ~2,650 for both. Module
boundaries come from the sequenced LV1.27.1 amplicon (865 bp, same template)
aligned to pLV018; the sensor's 3′ end and the target's extent come from that
part and from the expected product size, not from this template; the target
module's 3′ end is placed 70 bp before the product's, from the donor's own 3′
architecture, and the RNase lobe's boundary within SMS is figure 1's schematic
proportion rather than a sequence coordinate. (b) The LV01.04 test on MLS12.18.1.1 cDNA — one master
mix, KAPA HiFi, 50 µL, 24 cycles, three reverse primers. Sizes are read from
the GeneRuler 1 kb ladder in lane 1; the dashed box in each lane spans ±25% of
that lane's expected size, the region a product of about that length would run
in. A single product appears in the PR10137 lane, at 812 bp against ~850
expected, and nothing at ~2,650 in either long-design lane. One reaction per
design. The gel is flipped vertically so the wells sit at the top; nothing else
is changed.
