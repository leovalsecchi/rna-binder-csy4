# Fig4b_cell_populations

The four cell populations the library lands in. Follows Fig 4a.

Run: `python code/Fig4b_cell_populations/make_figure.py`
Output: `output/Fig4b_cell_populations/cell_populations.{svg,pdf,png}`

## Both drawings are reused whole

* The cluster is `output/cartoon_library/svg/cell_cluster.svg`, used whole
  and unedited - three cells, one large and two small, nuclei and all.
* The target is the Fig1f glyph for that construct - the same four files Fig
  4a's call-out uses, so a reader can lay the panels side by side and map a
  population onto its target.

**The target is not drawn loose in the cell.** A small circle marks a patch of
the large cell's cytoplasm and a cone of two outer tangents opens from it into
a zoom circle alongside, where the scaffold is drawn. The magnification is
then explicit - molecules in the cytoplasm, seen closer - rather than one
giant scaffold floating in a cell it could not fit inside. An earlier draft
put the glyph straight into the large cell in place of its nucleus; this
replaced it.

The cytoplasm ring is only about 27 wide between the nucleus (r 19) and the
membrane (r ~46), so `SRC` was chosen by searching that ring: it clears both
by 7, and sits on the side the zoom is on so the cone stays short.

**The label hangs under the cells, not under the whole unit.** The zoom
circle reaches well below the cluster, so measuring the label from the unit's
bounding box left it stranded far beneath - and closing that gap only pushed
it into the circle's lower arc. Under the cluster there is nothing to its
right at that height, so it sits close.

Laid out **2x2** in Fig 4a's arrangement: the two geometries across,
positive and negative down. **A is the Y array (mod1), B the circular
scaffold (mod4)**, set by Leo. 4x1 was tried and dropped - it came out about
8:1, too wide for anything but a full-width row.

A pale cytoplasm wash for positive against negative was tried and dropped
too. If it ever comes back: it has to be mixed to an opaque hex, not applied
as fill-opacity, because the three cells overlap and a translucent fill
doubles up along every seam.

**Each zoom reads top to bottom as a reaction**: the target, a plus, and
below it three binder variants with "or" between them. A cell takes a single
insert from the oligo pool, so it holds ONE of the three - which is what the
"or" says, and which is exactly what is unknown for any given cell. The
helices come from `common.glyphs.ribbon_helix` (a two-tone ribbon off real
helix geometry).

They are three **steps of the binder's own hue** (theme 400 / 600 / 800), not
three unrelated colours: these are three members of one library, so they have
to read as the same kind of thing (COLORS.md rule 2, colour follows the
entity). The spread is wide enough that they are plainly three.

The binders sit beside the target rather than docked on it: the negative
scaffolds have no site to dock to, and one drawing serves all four panels.

**The glyphs are fitted by their ink, not their viewBox.** A matplotlib SVG
is mostly margin: the Y array's drawing is 62 x 21 inside an 88 x 59 viewBox,
so three quarters of what was being scaled was empty padding, and the targets
stayed small through three rounds of enlarging. `glyph()` measures the path
coordinates instead. Fig 4a's call-out had the same bug and the same fix, with
an explicit inset there since its cells want a margin.

Sizing inside the zoom is also cramped by the circle: the content is three
stacked rows and a circle has no corners to put them in, so the scaffold's box
is the largest whose top corners stay inside it at that offset.

The grid is measured from the drawing's own bounding box, not the cluster
SVG's viewBox: about a quarter of that viewBox is empty, which left the labels
stranded and the clusters off-centre.

## Open points

* **No landing pad is drawn inside the cells.** Leo asked for it, then said
  not to: the integration is carried by Fig 4a and the caption, and a locus
  glyph at this size would be a smudge.
* Which construct is A and which is B is a guess from Fig 4a's reading order.
* Site counts in the glyphs are schematic, not the real 8 and 16.

Draft caption: The pooled library is integrated into four cell populations,
one per target construct - the circular scaffold (mod4) and the Y array
(mod1), each with its Csy4 sites present (positive) or abolished (negative).
Each population is drawn as a cell cluster carrying its target.
