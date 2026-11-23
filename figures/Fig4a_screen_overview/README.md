# Fig4a_screen_overview

Figure 4a assembled: the screen from cassette to reads to hits, in one
horizontal panel.

Run: `python code/Fig4a_screen_overview/make_figure.py`
Output: `output/Fig4a_screen_overview/Fig4a_screen_overview.{svg,pdf,png}`

## Nothing is redrawn here

This script reads three SVGs this repo already emits and only rescales and
places them:

| piece | source |
|---|---|
| donor cassette, oligo pool, target-scaffold call-out | `output/Fig4a_library_cassette/library_cassette.svg` |
| the four populations | `output/Fig4b_cell_populations/cell_populations.svg` |
| RNA extraction → min. reads filtering | `output/readout_workflow/readout_workflow.svg` |

The fourth piece, the hit gate, is drawn here — see below.

Change any of them and rerun this; the panel cannot drift from its parts.

**Run order.** `run_all.sh` walks `Fig*/` alphabetically, so
`Fig4b_cell_populations` runs *after* this script. A change to Fig 4b needs
its own script run first, or a second pass of `run_all.sh`.

## One label size, computed not eyeballed

Each source is brought in at `5.5 / its own label size`, so the panel has one
type size throughout and a source can be re-emitted at a different internal
scale without this script changing. The cassette and the populations already
draw at 5.5. The workflow is a single-column 85 mm panel in its own right and
draws its 7 pt labels at 13.65 of its units, so it comes in at **0.403**.

That one number is also why the workflow's circles land at about the size of
a population's zoom circle — the match follows from the type, it was not set
by eye.

The one deliberately smaller mark is the populations' "or", at 4.5. `_check()`
knows all three raw sizes by name and fails on a fourth, which is what would
happen if a source changed its type and its scale went unrechecked.

## The two links this panel adds

Plasmid → cells ("Bxb1 landing pad transfection and integration") and cells →
RNA extraction. Both are **the workflow's own connector**, scaled to the
panel, because both are process steps exactly like the ones inside the
workflow. The cells → RNA extraction link is deliberately NOT a funnel or
zoom cone: the cone already means "seen closer" twice inside the very box it
would be leaving, and the box's outline is what says all four populations
feed it. A bracket was the other option and would have duplicated that
outline.

The integration label sits beside its arrow, not under it: the gap is 20
units and two stacked lines would leave the arrow no room to be an arrow in.

## Layout

**The left column's two containers share a left edge.** The cassette is
placed by its *ink*, not by its viewBox: that viewBox is padded about 13 units
left of the plasmid, which left the plasmid visibly inset from the box below
it. `ink_x0()` measures the paths; `_check()` re-tests that the two edges
still agree.

**One gap between blocks, everywhere.** Plasmid down to the populations box,
and that box across to the workflow, are the same `GAP`. Leo spaced those by
hand at 16.2 and 20.0 and asked for them even; both are now the larger, since
the smaller would not hold the integration arrow and its label. `_check()`
fails if either drifts — including if something in the left column grows
wider than the populations box, which would silently widen the second gap.
`CORNER_PAD` is left out of this: it is clearance from a label to a frame
tucked inside the L, and at 20 it would eat 13 units of the plot.

**The workflow's column pitch (`COL_DY`) is the pitch Leo set by hand**, made
even: he spaced it 44.0 / 39.6 / 38.4, and 101 readout units is that average
applied to all three.

**The workflow is cropped to its own ink** (`RO_INK`, derived from its step
positions). It carries a panel's margin, which inside this figure is just a
hole between it and the populations box.

**The plot frame sits in the L's own empty corner.** Where that corner starts
— right of the column's longest label, below the row's deepest one — is read
off the workflow script by importing it, not typed in here. Move a step or
re-word a label there and the frame still clears it; `_check()` fails if it
does not, and its right edge is flush with the workflow's, so the panel has
no dead margin on that side.

The workflow is centred vertically against the left column, which is the
taller of the two.

## The hit gate is Fig 4c, redrawn in panel units

`gate_plot()` reads `inputs/screen_points.csv` — Fig 4c's own input — and
applies Fig 4c's own gate (reads ≥ 75 in POS **and** matched NEG, splicing
ratio ≥ 0.07, fold change ≥ 10). So the cloud's shape, and where the hits sit
in it, are the screen's rather than a drawing of what a screen looks like.
It also keeps Fig 4c's encoding — POS_A / NEG_A on theme 700, POS_B / NEG_B
on theme 400, everything out of gate in one neutral grey — its frame, ticks,
axis wording and its legend, so a reader meets one plot twice rather than
two. `_check()` fails if the gate stops yielding Fig 4c's **7 + 25 = 32**.

**Why it is redrawn rather than embedded.** Dropping
`output/Fig4c_fc_vs_splicing/fc_vs_splicing.svg` in whole would bring its own 9, 9.5, 10.5
and 15 pt type into a panel whose premise is one size, and it is square where
this corner is about 2:1. Redrawing it in panel units is what keeps both.

Two departures from Fig 4c, both forced by the size: the legend's third entry
("outside gate") is dropped, since one grey against two coloured series reads
without a swatch and the row would otherwise overrun; and the legend starts at
the box's left edge rather than the frame's, since it is about 12 units wider
than the frame. The alternative was to lengthen the workflow's arrows until
the plot grew to fit a legend — the layout driving the wrong thing.

The two gate values are named on the lines themselves (`0.07`, `FC 10`), in
the empty bottom-right quadrant. Fig 4c leaves them to its caption; a panel
has no caption to leave them to.

Only the out-of-gate cloud is thinned, to `BG_N` points, with a fixed seed so
the panel is reproducible. All 32 hits are drawn.

**It fills the L's corner exactly**: left edge clear of the workflow's column
labels, top clear of its row labels, right edge flush with the workflow's own.
Lengthening the row's arrows is therefore what widens the plot; they are one
move, not two numbers to keep in step.

## Open points

* No panel letter is drawn (STYLE_NOTES section 11); **a** goes on at
  assembly.
