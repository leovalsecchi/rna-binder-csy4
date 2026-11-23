# Figure 5 panel c: the four best hits as molecular surfaces, X-rayed

The four best hits by flow fold change, drawn as molecular surfaces, with a
square window cut out of each one and the cartoon showing through it — the
envelope outside the window, the fold and the bases inside it. See **The X-ray
window** below.

**This panel was d and is now c** (Leo, 2026-11-18: "replace the renderings in c
by the ones in d"). It began as a second row under `code/Fig5c_best_hits`, the
cartoon close-ups of the same four designs, and took that row's place: one row
of four, labelled, not two rows of the same four. `Fig5c_best_hits` still builds
and is still the source of the four cameras, imported here, but it is no longer
promoted.

The asparagine panels that used to be d/e/f came off the page on the same day
("probably d-f will be removed. You can remove them for now"); the rows are
kept, commented, in `code/Fig5_assembly/make_figure.py`, and
`code/Fig5def_asn_interface` still builds untouched.

## What is shown

Ranking is `flow_FC` in `LV_32Candidates/index.tsv`, the LV2.34.1 `ONOFF_avg`
column, mean(PosA, PosB) / mean(FullNeg, NegA).

| panel | twist | FC | stack | design | model used |
| --- | --- | --- | --- | --- | --- |
| 1 | 23 | 1.56 | BG+LiMPNN | e8cb617b_design26_seq2 | AF3 |
| 2 | 15 | 1.54 | BG+BIF | 25a69495_design18_seq0 | ESM |
| 3 | 16 | 1.47 | RFD3+LiMPNN | 73d23e78_design49_seq4 | ESM |
| 4 | 07 | 1.41 | BD1+LiMPNN | c5072321_design5_seq1 | ESM |

Each candidate PDB holds three superposed models (1 = design, 2 = AF3 refold,
3 = ESMFold2 refold, per its own REMARK). The model with the lowest
self-consistency RMSD is used, never AF3 by default — AF3 is the worse model for
three of the four. Same choice as `Fig5c_best_hits`, verified there against
`inputs/csy4_binder_library_12k_v6_k30.csv` (`sc_rmsd_af3`, `sc_rmsd_esm`).

RNA is `thesis_grey_400`, each binder its stack colour, and **nothing else is
coloured**. `--hilite` paints the RNA-contacting residues `thesis_petrol_800`;
it was tried on the panel and dropped — a second colour over the stack colour
reads as a chemistry claim the panel cannot support, and the shape is the point.

Soft lighting, depth cueing off, `thesis_grey_700` silhouettes at width 1.5,
transparent background. The panel carries no labels of its own: it sits under c,
which already labels the same four in the same order, so repeating them would be
the same four labels twice two centimetres apart. `LABELS = True` in
`make_figure.py` puts them back for standalone use.

## The X-ray window

Leo, 2026-11-18: *"keep the exact positioning and sizing of figure D, and then
render it in the same style as figure C was rendered ... you use the surface
one as the background, the larger background, and then you just crop the figure
C style one, have it on top, crop it, and then it looks kind of like an X-ray
into the figure."*

So each panel is **two renderings of the same design through the same camera**:

| file | what | drawn as |
| --- | --- | --- |
| `w<twist>.png` | the background | solid surfaces, this panel's own recipe |
| `x<twist>.png` | the window | panel c's cartoon, turned up (see below) |

`render.py --xray` does not re-frame anything. It builds, colours and frames on
the surfaces exactly as the normal run does, and only then swaps the
representation for c's and saves. The camera is therefore the surface run's
camera to the last digit and the two layers land on the same pixels - which
re-framing on the cartoon would NOT give, since `view` fits whatever geometry
is displayed and a surface stands ~3 A proud of it. It also turns depth cueing
back on, because c has it on and this layer has to look like it was rendered
there.

The cartoon started as c's recipe setting for setting and has moved off it
twice, both on Leo's call, 2026-11-18:

- *"remove the dark aminoacids. They should be the same color."* c paints the
  contacted bases a darker grey and leaves them opaque while the rest of the
  RNA is 60 % transparent. At 45 mm that marks them; inside a window a sixth
  that size it is a dark blob. One grey across the whole hairpin now. (They are
  nucleotides, not amino acids - no side chain has been drawn here since the
  sticks came off - but the call is the same either way: nothing in the window
  takes a second colour.)
- *"make the hairpin and protein in the image a bit more color intensive"* -
  every transparency roughly halved: protein 65 -> 40 %, RNA 25 -> 10 %,
  nucleotides 60 -> 35 %. c's values are tuned for a 45 mm panel and wash out
  at 13 mm.

`make_figure.py` cuts the window: a square, hard-edged, with a hairline rule
around it (Leo, 2026-11-18: *"a square with a clear cut (no blur needed (just
thin lines))"*). Where it sits is **measured, not chosen**: `render.py` aims every panel's
camera at that design's own interface - the contact residues and the bases they
read - so the interface is at the centre of every saved frame, and the window
rides along with whatever shift `lay_out` gives that panel. Four constants:

    XR_FRAC = 0.24    HALF the window's side / the panel's shorter side
    XR_LW = 0.5       pt, the rule around it: the book's thinnest
    XR_EDGE = "600"   its neutral step
    BG_IN = 0.45      how much surface is left inside the window

`BG_IN` is the only real judgement call. 0 is a clean hole; a ghost of the
envelope over the cartoon is what makes it read as an X-ray rather than as a
bite taken out of the surface. 0.65 is where the envelope starts swallowing the
cartoon, worst on twist 23, whose binder wraps the hairpin. 0.45 holds only
because the cartoon was strengthened in the same pass; at c's transparencies it
was mud much earlier.

Sizes rendered: a disc of radius 0.42 first (*"the zoom in is too large"*),
then squares at half-side 0.30, 0.25, 0.20 and 0.24. The ceiling is where the
window eats so much of the surface that what is left falls into disconnected
islands, which read as a rendering fault; 0.24 keeps every surface one
connected mass, so the square reads as a window into it.

`XRAY = False` in `make_figure.py` gives the plain surface row back.

## The layout (`make_figure.py`)

One box SIZE for the four, so they stay at one scale and one baseline, but each
structure centred **sideways** in its own copy of it. The four cells are
equally spaced, so the structures standing in them have to be too; a single
shared crop box holds each one where its camera left it, which put the widest
binder against its cell's edge and left a hole beside the next one. The box is
as wide as the widest ink, so nothing is cut and the narrower panels simply get
transparent margins. Vertically the box stays shared, which keeps the row on
one line.

The structures are still not all the same width on the page, and should not be:
the hairpin is the same size in all four (one camera, see below), and what
differs is the binders, which really are that different - twist 23's is compact
and wraps the hairpin, twist 16's is long.

## The frame (`--frame hits`, the default)

The four per-panel cameras are **imported** from `Fig5c_best_hits/render.py`,
not copied, so the two panels cannot drift apart: its `PANELS` table supplies
each design's `turn y` angle (60 / 60 / 90 / 30), the residues and nucleotides
the view is framed on, and which model to use. The alignment target is twist
23's own chosen model, as there — the C1' fits reproduce that script's own
numbers exactly (0.000 / 1.486 / 1.111 / 1.374 Å).

**One scale across the row: one camera, moved, never re-framed.** The camera is
global — `view` on the reference hairpin, then `HITS_ZOOM` — and a panel's only
private setting is a `move x/y` that slides its own interface to the middle of
the frame, where `make_figure.py` cuts the window. Translation cannot change
scale, so the hairpin, the same molecule in all four, comes out the same size in
all four. The move is undone after the save, so the next panel's
`turn y center 0,0,0` still turns about the origin.

This replaced a per-panel `view` on each design's framing selection, corrected
by that selection's radius (15.3 / 9.2 / 15.5 / 8.5 Å — c frames twist 23 on
three residues and two nucleotides, twist 15 on one and one). The correction was
close but not exact, because `view` fits a bounding BOX while the correction
divided out a centroid RADIUS: measured on the renders, the hairpin came out a
fifth smaller in twist 7's panel than in twist 23's (Leo, 2026-11-18: "check
that they are visually all about the same size and spaced correctly").
`scale_check.py` measures it.

Each panel draws its **own** RNA copy. Fig 2b and Fig 3b draw only the reference
RNA because their binders are cartoons laid over it; here both are surfaces, and
pairing binder *i* against a hairpin it was never modelled with would invent
clashes and gaps at exactly the interface the panel is read for.

## The other two frames

`--frame fig4` is the Fig 2b shared RNA frame plus `turn y 120`, the camera Fig
3b and the 32-hit grid of Fig 4 use, referenced to the book's own hairpin
(`working/structures/fig02b_examples/cifs/csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0.cif`).
Use it if these four should read as four of those thirty-two rather than as
companions to c.

`--cut --frame axis` is the cut-away this folder was started for: the camera
looks straight down the hairpin's helical axis from the loop end and a 10 Å slab
is clipped out, so the reader looks into a section. The slab is centred on the
**measured** interface — the mean depth of every contacting protein heavy atom,
pooled over the four — not on the hairpin head. That is a result in itself: the
head centroid is at z +10.97 Å and the interface centroid at +0.62 Å, so **these
four grip the stem, not the loop**, and a cut at the head comes back nearly
empty for three of them (`--probe` renders twist 23 at four depths and shows it).

All write into the same folder: `w*.png` whole, `x*.png` the cartoon X-ray
layer, `c*.png` cut, `z*.png` zoomed-and-recentred. A run clears only its own
prefix, so one variant never wipes another. `make_figure.py`'s `VARIANT` picks which the panel
assembles; it is `"w"`.

## Run

    ~/.venvs/binder_map/bin/python render.py ../../output/Fig5c_cutaway_surface --width 2000
    ~/.venvs/binder_map/bin/python render.py ../../output/Fig5c_cutaway_surface --width 2000 --xray
    ~/.venvs/binder_map/bin/python make_figure.py ../../output/Fig5c_cutaway_surface ../../output/Fig5c_cutaway_surface
    ~/.venvs/binder_map/bin/python ../Fig5_assembly/make_figure.py
    ~/.venvs/binder_map/bin/python ../Fig_book/make_book.py --promote 5:Fig5_V2_adjusted

The other two variants:

    ... render.py ../../output/Fig5c_cutaway_surface --width 2000 --frame fig4
    ... render.py ../../output/Fig5c_cutaway_surface --width 2000 --cut --frame axis

## VERIFY markers

- `VERIFY:` the caption must say that **twist 15 fails the fidelity gate**. Its
  best-refold DockQ is 0.020, the worst of the 32, so neither refolder
  reproduces the complex drawn for it. It is in because it is flow-FC rank 2.
  The other three are above the 0.23 gate.
- `VERIFY:` the caption must name the model used per panel. Three of the four
  are ESMFold2 re-predictions, not designed coordinates, and a surface panel
  reads as if it were a structure.
- `VERIFY:` the caption must say d is c's four in c's order, since d carries no
  labels of its own.
- No claim about shape complementarity, buried area or contact quality is
  supported by this panel on its own; it shows poses, not energetics.

## Traps hit while building this, worth not re-learning

- **`view <spec>` on a spec with no displayed ATOMS kills ChimeraX silently.**
  The error ("No displayed objects specified.") goes to the GUI log, never to
  stdout; the command file aborts, the process still exits 0, and not one PNG is
  written. A molecular surface over the atoms does not satisfy it. This is why
  the `--surface` mode of `code/Fig5c_best_hits/render.py` reported
  `wrote 0 png`. Show the atoms, frame, hide them again.
- **Dump the ChimeraX log to see any error at all**: a second
  `--cmd "log save <file>"` runs even when the first `--cmd` aborted.
  `render.py` does this and greps the result.
- Depth cueing is on by default and, on a white page, fades the far half of each
  binder into a pale smear that reads as a rendering fault. `lighting depthCue
  false`.
- Transparent surfaces were tried at 40/55/70/82 % and dropped: past ~55 % the
  stack colour stops being readable and the grey RNA disappears into the page,
  which costs the panel both of the things it encodes. `--alpha` keeps the knob.
- `import render` from this folder would import **this** file — the panel script
  next door is also `render.py`. It is loaded by explicit file path instead.
- Verify a promote by CONTENT. `make_book.py` reports success whether or not the
  upstream build worked. The check here: `selection/Fig5.svg` hashes identical
  to the page just built, panel letters a/b/c/d present, and zero occurrences of
  "Asn" inside `<text>` elements (grepping the raw file matches base64 image
  data by chance — 127 false hits).
