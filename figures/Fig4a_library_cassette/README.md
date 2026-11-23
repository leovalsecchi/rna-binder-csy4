# Fig4a_library_cassette

The library cassette before integration. Rebuilt 2026-11-12.

Run: `python code/Fig4a_library_cassette/make_figure.py`
Output: `output/Fig4a_library_cassette/library_cassette.{svg,pdf,png}`

## This is Matteo's donor plasmid, moved — not a redrawing

The script reads
`output/cartoon_library/integration/Promoterless-donor-plasmid.svg`, the first
row of the cassette-integration figure, and does exactly two things to it:

1. **The plasmid outline is left untouched**, byte for byte: his rounded
   rectangle, 0.2 to 150.329 across, 11.78 to 43.78 down.
2. **Everything else is translated up by 32**, which is that rectangle's own
   height — the distance from its bottom edge to its top edge. The cassette,
   which sat on the bottom edge, now sits on the top edge and the loop hangs
   below it.

Nothing is scaled, mirrored or re-spaced. Every block keeps its geometry,
width, colour and spacing; the brackets and label ticks travel with the blocks
they belong to, so Recorder / Sensor module / Target module stay above the
strand and the ATG, Binder, SMS and att captions end up inside the loop.

Rerun this after any change to `retheme_integration.py` — the blocks come from
its output, so the two figures cannot drift apart.

## The two additions

* **Oligo pool** into the binder, the one variable block. The sensor-module
  bracket runs straight over the binder, so the arrow has to cross it; it
  crosses the way Matteo's own exchange curves do, with a white casing under
  the ink so it reads as passing in front.
* **RNA presenting module call-out**, the scaffold blown up into the four screened
  constructs: mod1 (Y array, "A") and mod4 (circular, "B"), each with its
  Csy4 sites live or abolished. Those four glyphs are the files Fig1f already
  emits (`output/Fig1f_scaffold_splicing/elements/`), reused rather than
  redrawn, in the same order as Fig 4b's 2x2 read row-major.

  **The leaders land on the target-module BRACKET's ends, not on the block.**
  Running them to the block means crossing the bracket row, and once the
  call-out went wide the left one cut straight through the sensor-module
  bracket. Stopping at the bracket also says the right thing: what the
  call-out expands is the module that bracket names.

  They start on the straight part of the box's bottom edge, inset by the
  corner radius: from the rect's corner point a leader springs out of the
  rounded curve and reads as misaligned.

  **The left leader has to thread a 14-unit corridor** between the "Sensor
  module" and "Target module" captions, which is what fixes `BOX`'s left
  edge. `check_leader()` samples each leader against the real caption boxes
  and stops the build if one crosses. That check exists because the clash was
  "fixed" twice by testing the caption's BASELINE only - the leader cleared
  at the baseline and cut through the text's top edge both times.

  It is **one flat row**, not a 2x2 block, and its **title sits inside the
  box** on its own band rather than floating above it. Both are for height:
  the call-out was what made the panel tall, and the gap over its title was
  the largest piece of nothing in it. The on/off row labels went with the
  2x2 - live and dead sites are already warm against grey, and Fig 4b names
  the four.

The only deletion is his "Promoterless donor plasmid" heading: it is a title,
and titles belong in the caption (STYLE_NOTES section 1).

## Open points

* **The loop is reshaped** (`reshape()`), the one change to his rectangle:
  widened 4 on each side (`LOOP_WIDEN`) and its floor raised 14
  (`LOOP_LIFT`), from 32 deep to 18. His rectangle was 32 deep because it sat
  ABOVE a cassette; below one it is mostly empty and was taking a third of
  the panel. Each wall's points move by the same amount rather than the shape
  being scaled, so every corner arc stays the circle he drew, and the top
  wall - the strand - never moves. The widening is there because `attB-GT`
  and `attB-GA` are wider than the gap between their site and his corner, so
  once the cassette moved up and those captions landed inside the loop they
  grazed the walls. Clearance is now 3.9 left, 3.2 right, 4.6 below the
  captions' descenders.
* **Not to scale**, deliberately: the binder is drawn small and the scaffold
  large so the later scaffold panels follow on. In base pairs it is the other
  way round. Say so in the caption.
* Site counts in the glyphs are schematic, not the real 8 and 16.
* The follow-on panel is `Fig4b_cell_populations`.

Draft caption: Library cassette before integration. A pooled oligo library
supplies the binder, the one variable block; the recorder (AT hairpin), the
sensor (binder–SMS) and the target scaffold are otherwise fixed. The scaffold
position carries one of four target constructs — mod4 (circular) and mod1
(Y array), each with Csy4 sites present or abolished. Flanking attB-GT and
attB-GA sites recombine into the genomic landing pad. Not to scale.
