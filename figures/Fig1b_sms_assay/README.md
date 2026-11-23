# Fig1 panel B - SMS reporter assay

Output: `output/Fig1b_sms_assay/`. Rebuild with
`python code/Fig1b_sms_assay/make_panel.py`.

**The figure number is a guess.** It is one string,
`paths.fig_output("Fig1b_sms_assay")`, on the assumption this is the assay
belonging to the Fig1 pipeline strip.

## What it shows
- **Sensor, left.** The RNA target in its F30 scaffold, held by two SMS
  receptors, with mNeonGreen translated from the same transcript.
- **Constructs, right.** The three transfected pieces, ruled rather than boxed.
- **States, below.** The reporter mRNA before and after cleavage.

## It draws none of its own parts
| part | comes from |
|---|---|
| construct bars | `common.constructs` - the same `BARS` entries panel d uses |
| scaffold | `cartoon_library.f30`, **imported, not copied** |
| hairpins | `v19.18.svg` through `common.cartoon` |
| binder, SMS, GFP | Matteo's 1x components through `common.cartoon` |

Importing `f30` rather than embedding a copy means the panel follows any change
to that model without being edited. The receptors are placed on the scaffold's
**actual arm tips**, computed from `f30`'s own constants, not at guessed pixels;
if those names change, `_arm_tips()` raises rather than letting the receptors
drift off the arms quietly.

## Two things that were wrong in the sketch
- The **stop codon sat after mScarlet**. It has to come first - that is what
  holds the reporter out of frame until the hairpin is cleaved.
- The hairpin was drawn *and* boxed as `AT` on the same strand. The hairpin is
  the AT module, so in the states it is drawn once, as a hairpin.

No title, no standfirst, no panel letter (STYLE_NOTES.md items 1 and 8).

## Open points
- The scissors are gone: cleavage is the arrow between the states plus the
  hairpin drawn cleaved, and the dashed line from the activated receptor.
- `AT` is the spliced hairpin; expand the name if Methods calls it something
  else.
- The scaffold is drawn with two target sites. The real ones carry more
  (see Fig1f_scaffold_splicing); two is schematic, not a count.

---

## 2026-11-12 — rebuilt as `make_panel_b.py`

`make_panel.py` above is **superseded**. It drew its parts here from scratch;
Leo asked for his own and Matteo's. The new builder is
`make_panel_b.py` → `output/Fig1b_sms_assay/Fig1b_panel_b.{svg,pdf,png}`.

**Leo's arrangement is the point of the rebuild**: mechanism down the left,
the uncut reporter beneath it, a Frameshift arrow across to the spliced
reporter making mScarlet, and a grey box of components along the top.

| slot | element | where from |
|---|---|---|
| mechanism | `Fig1d_mechanism_hairpin_V2.svg` — the **full four-state** version, not the compact one Fig 1d uses | Leo's Illustrator arrangement of Matteo's VHH v7.9, SMS v6.5, linker, target, hairpin |
| hairpins in the rows | `HP v19.18 Unspliced` / `Spliced` | Matteo's `v19.18.svg` via `common.cartoon` |
| mScarlet protein | his GFP barrel, green swapped for `synbio.reporter_alt` | Matteo's `Selected_GFP_v1_1x.svg` |
| module blocks | his idiom: rounded box 6 deep, hairline in `synbio.outline` | same as the retheme pipeline gives the cassette |

**Nothing from `f30.py`, `receptor.py`, `mechanism.py` or `constructs.py`.**
Those are ours. If a part seems missing, take it from `inputs/cartoons/`, not
from the drawn library.

### Two traps
- **Keep the mechanism's `<defs>`.** It carries the `<style>` block that fills
  his shapes; stripping it renders the whole mechanism as solid black
  silhouettes.
- **The excised intron imports neutral.** The standard map in
  `cartoon_colors.json` treats it as structure, but COLORS.md rule 4 gives the
  warm end of the diverging ramp to the piece the hairpin loses, and the
  mechanism already draws it that way. `hairpin()` recolours that one path
  locally rather than changing the global map.

### Open
- The glow around the barrel is the one mark with **no source in Matteo's
  set** — he has no lightbulb. A few rays in `synbio.enzyme_alt`.
- The two rows carry no stop codon glyph yet; "uncut = stop in frame" is
  currently told by the grey block alone.
- Space under the component box is empty; the box or the rows may want moving
  once the panel is placed in the deck.

## `fix_fluorescence.py` — Leo's Illustrator template, corrected

`working/Fig1d_integrated_construct/Fluorescence Splicing V1.svg` is Leo's
layout. The script only corrects it; it moves no geometry except the grey
container, which it redraws as a clean rounded rect.

Out: `working/Fig1d_integrated_construct/Fig1b_fluorescence_splicing.{svg,pdf,png}`

| was | now | why |
|---|---|---|
| two stroke weights: none (so 1) on the art, 1.76 on the blocks | **0.73 everywhere** | his blocks are 10.95 deep where Matteo's are 6, so the file is at 1.825x and Matteo's 0.4 hairline becomes 0.73 |
| `#d4145a`, `#22b573`, `#e6e6e6`, `#000000` | reporter_alt, reporter, neutral 100, ink | off-palette placeholders |
| 3 sizes, 2 families, one text x-squashed | **Arial, 10.0 on the page** (8.38 before the 1.193 matrix); in-block captions 8.0 | — |
| coloured squares | **Matteo's GFP barrel**, green as mNeonGreen and recoloured as mScarlet | they were placeholders |
| unlabelled blocks | mNeonGreen / RNA scaffold / Binder / SMS / AT / mScarlet | — |
| "mScarlett" | "mScarlet" | — |

**It matches `Fig1b_panel_b.svg`**: that draws 6-deep blocks with a 0.4 hairline
and 5.5 type; this is the same drawing at 1.825x, so 0.73 and 10.0 are those
numbers multiplied through. Caption height over block depth is 0.913 here and
0.917 there. Place the two at the same block depth and they agree.

### The trap that cost a render
The quoted font family is written `&apos;Source Sans 3&apos;`, and **that
entity ends in a semicolon**. A `font-family:[^;}]+` stops inside it and leaves
an unclosed quote, which makes every CSS rule *after that point* unparseable -
so the component box rendered solid black while the art above it looked fine.
Match the quoted form explicitly.

### Open
- "SMS" is centred on the kinase lobe, not the whole outline, or it straddles
  the kinase/RNase boundary. The RNase lobe is unlabelled.
- Bar 1 has no 5' UTR block, unlike the other panel's component list.
- The barrels are at `BARREL_K` 1.5; his placeholder squares were smaller.

## Audit against the old panel b (2026-11-12)

Cross-referenced `Fig1b_fluorescence_splicing.svg` against the slide-deck panel
b, and every element name and fill against the three figures that draw the same
blocks.

### Names: his ids were the ids of whatever he duplicated

`rename()` fixes six. Three of them are **keys in
`retheme_integration.py`'s `FILL_BY_ID`, which colours by id** - left alone, a
re-theme of this file paints the mScarlet CDS reporter-green and the grey
container theme 200.

| was | is | why |
|---|---|---|
| `RNA_scaffold` (the grey container) | `Component_box` | copied from the scaffold block |
| `RNA_scaffold-3` (green, bar 1) | `mNeonGreen_replacement_cassette` | Matteo's name for that block |
| `RNA_scaffold-2` (theme 200, bar 1) | `RNA_scaffold` | the real one |
| `mNeonGreen_replacement_cassette` (red, bar 3) | `mScarlet_encoding` | copied from the landing pad |
| `SMS_kinase_lobe` / `SMS_RNase_lobe` | `SMS_kinase` / `SMS_rnase` | the ids `_fuse_receptor` ships |
| the two strand CDSs (unnamed) | `mScarlet_encoding_{unspliced,spliced}` | — |

Renaming is keyed on the **id**, and the data-name is rewritten with it: his
two disagree (`id="RNA_scaffold-3" data-name="RNA scaffold"` on the green
block), so renaming by text alone misses half of them.

`Fig1d_integrated_construct.svg` was the one file still using Matteo's
pre-fusion `Kinase fill` / `RNase fill`; it now says `SMS kinase lobe` /
`SMS RNase lobe` like the cassette and Fig 4a, and `Fig1d_full` was rebuilt.
Binder, Linker, SMS kinase lobe, SMS RNase lobe, SMS outline, RNA scaffold and
the three hairpin segments now carry the same name **and the same fill** in
Fig 1b, Fig 1d, `cassette_integration` and Fig 4a.

### Alignment: nothing in the component box sat on its own line

His three lead lines are an even grid (27.277 / 54.071 / 81.133); the blocks
had drifted off it, so the line ran through each bar at a different height.
`snap_rows()` moves the blocks, not the lines - the lines are the grid.

| | was | moved |
|---|---|---|
| bar 1 | its two blocks 0.720 apart, line between them | +0.090 / −0.629 |
| bar 2 | line 1.096 below the blocks | +1.096 |
| bar 3 | line 1.126 below, and the AT block 0.892 right of the other bars | +1.126, −0.892 |
| 1 / 2 / 3 | 2.87 / 2.57 / 1.40 above their rows, the "2" 0.185 off the column | centred |

### What the old panel had and this one did not

* **AUG and STOP** - added, as a tick off the strand with the codon under it,
  Fig1d's idiom and Matteo's own (`ATG site` + `Start codon`). Without them
  the Frameshift arrow acts on nothing: the point is that the stop is *before*
  the CDS unspliced and *after* it spliced, which is now what the panel shows.
  The viewBox grew 197.5 -> 207.5 to hold the codon row; no existing element
  moved.
* **The callouts 1 2 3** - added. The numbers in the box were orphans, with
  nothing in the panel saying what they named. Ring radius is 5.5 because
  bigger does not fit the 13.2 between his cut-site arrowhead and the mScarlet
  bracket.
* **"RNA Target in Scaffold" / "SMS Receptor" / "XBP1 Hairpin and Reporter"**
  - not added. The rings carry the same information and the box spells each
  one out; the bold captions are a title row (STYLE_NOTES section 1).
* **A 5' UTR block in bar 1** - not added, and this settles the question left
  open last pass. Matteo's own name for the AT block is
  `Hairpin encoding | 5 prime UTR`: in his set the hairpin **is** the 5' UTR.
  A second, generic UTR block would make "5' UTR" mean two things in one box.
* **Scissors at the cut sites** - not added; his file marks them with the two
  black arrowheads, which is less ink and already in place.
* **A lamp beside each barrel, and mNG / mScar labels** - not added. There is
  no lamp in Matteo's set, and the barrels are the only green and red objects
  in the panel, so the box's own blocks key them. There is also no room: the
  green barrel is 5 units from the RNA target.

### The one residual difference from Fig 1d

Fig 1d's blocks are 6 deep over a 78.6 receptor; his are 10.95 over 115.84.
Depth scales by 1.825 but width only by 1.474, so his blocks are 24% chunkier
than Fig1d's. Everything keyed to depth matches (stroke 0.73 = 0.4 x 1.825,
type 10.0 = 5.5 x 1.825, caption-over-depth 0.913 against 0.917). Fixing the
aspect means re-cutting his bar widths, which is his layout, so it is left.

Weights are Fig1d's two, not one: **0.73** for blocks, strand and leads,
**0.64** for brackets and ticks (Fig1d's 0.4 / 0.35). Matteo's barrel keeps
its own 0.72 outline and 0.45 inner cap - his hierarchy inside the glyph.

## Second pass: widths, names from the key, smaller annotation

**Every bar now spans the receptor.** Bar 2 is one polypeptide at Fig1d's
proportions, so its span (245.640 -> 361.482) is fixed and the other two are
cut to it. Only the **gaps** absorb the difference - 16.966 in bars 1 and 3 -
so no block's width changes except the AT one.

**The AT block was redrawn, 21.892 -> 38.646 across.** His 8.756 intron
segment was narrower than the "AT" sitting on it, so the caption ran onto the
blue exons either side; the intron is now 15.46 and holds the caption with
2.2 clear. His 0.30 / 0.40 / 0.30 exon-intron split is kept. The corners were
also 0.876 horizontally against 1.460 everywhere else - the block had been
squashed in x - so it is rebuilt with Binder's circular 1.460.

**Captions now follow their blocks.** They were absolute coordinates, so the
first pass's row snap left all six behind by 0.09 to 1.13. Each is placed from
its block's id after the block has moved, which also puts "AT" on the intron
segment rather than across the whole block (the same reason "SMS" sits on the
kinase lobe).

**Names from `construct_key.svg`,** which is the key for these: `at` -> "AT"
(so that one was right), `aug` -> "AUG", `stop` -> **"\*"**, and
`mscarlet_off` -> **"mScarlet\*"**. "STOP" is not a name used anywhere in the
set, so the stop ticks are asterisks and the unspliced reporter is
`mScarlet*`. The same asterisk in both places is what makes the lone glyph
readable: the stop marked `*` is what puts the reporter out of frame.
`panel_b_sms_rna_sensor.svg` already read "AUG * mScarlet".

**Smaller annotation.** Codon captions and ring digits moved from the body
10.0 to the in-block 8.0; ticks from 6.8 to 3.6 long (Fig1d's 2.1 x 1.825);
rings from r 5.5 to 4.4, stroked at the lighter 0.64. The asterisk is drawn
at 12.0 and dropped 2.78, because an asterisk is a small glyph in its em and
at 8.0 it reads as a speck. Page height 207.5 -> 202.5.

**Callouts moved** to (73.5, 17.5), (74.0, 95.0) and (91.6, 178.0); a
bounding-box check against every path in the file reports no overlap for any
of the three. 2 had to clear the linker curve, which reaches x 72.2 down to
y 89.2; 3 threads the 13.2 between his cut-site arrowhead (ends 84.9) and the
mScarlet bracket (starts 98.09), level with that bracket.

## Third pass: the component box, and the reporters

### The box is laid out from its content

His was 230.4 across holding 152 of content, with the numbers 26 from their
bars and a 22.0 lead stub on the left against 5.6 on the right - the numbers
read as a detached list rather than row labels. Now every number is derived:

    MARGIN  15.85   one margin all round, = the white between two rows
    STUB     9.0    lead either side of the blocks, equal
    NUM_GAP  9.0    ring to the start of its lead

from a fixed right edge of 391.3, which is where he had it and is level with
the strand below. That gives the box 183.3 x 96.5 with all four margins at
15.85 exactly, and the blocks shifted 4.968 right. The ring digits dropped
from the caption 8.0 to their own 6.6, so a digit sits in its ring instead of
filling it.

**It leaves a wider gap between the mechanism and the box** - the content is
only 152 wide and a tighter box cannot cover 230. Re-cutting the bars to the
full 1.825x of Fig1d (span 115.84 -> 143.4) would close most of it and fix
the 24% aspect difference in the same move; it needs his say-so, since the
bar widths are his.

### The reporters

`Selected_GFP_v1_1x.svg` is the only fluorescent-protein glyph in Matteo's
set - his cytosolic sheets draw GFP as the same cylinder - so it stays. What
made it read as a tin can was ours, not his:

* **scale 1.5 -> 1.2.** His own drawings use it near 1x.
* **the rim.** `cartoon_colors.json` maps his `#6C6970` to neutral 700; put
  back to neutral 600, which is his own value to within a hair.
* **a lit cap.** The lens between his dome and his rim arc, filled white at
  0.3. A cylinder with no lit top is a flat shape.
* **it shines.** A faint halo and eight rays off the rim, in the protein's
  own colour, at the tick weight. There is no glow in his set and the slide
  used a clip-art lamp; rays are an annotation, and they say the protein
  itself is the light rather than putting a lamp beside it.

The red one is lifted 4.5 off his placeholder centre so his translation arrow
stops short of the rays instead of running into them. Barrel strokes are
0.58 and 0.36 against the panel's 0.73, which is his own ratio inside the
glyph carried through at 1.2 - deliberate, and of a piece with the lighter
rim colour.
