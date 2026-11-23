# Figure 5: the arrayed screen and what separates the binders

Assembled on an A4 page by `make_figure.py`. Two versions are written:
`Fig5_V1_as_sketched` (panels at the size they were drawn) and
`Fig5_V2_adjusted` (every panel redrawn to the page block at house type).
**V2 is the figure.** V1 exists to check the sketch, and its panel c still
carries 7.5 pt type because it places the standalone SVG rather than
re-requesting it.

Rebuild:

    python3 code/Fig5_assembly/make_figure.py

Run it with no argument. The output path must be absolute, because panel
scripts are invoked with their own cwd and a relative path lands under
`code/` instead.

Current geometry, V2: block 178 mm, scale to page 1.000, figure ends 270.3 mm
down a 297 mm page, 21 mm left below.

## The argument

The page moves from the assay to the hits. Panels a and b establish that 32
ordered binders were screened and how they rank. Panels c and d show the best
four of them as structures, the same four in the same order: c the contact in
detail as cartoon close-ups, d the shape that contact sits in as molecular
surfaces, in c's own four cameras.

**2026-11-18: d, e and f were removed** (Leo: "probably d-f will be removed.
You can remove them for now"), and the surface row took the letter d. They were
the negative point of the page - that none of the twelve metrics the field
ranks designed binders on separates binder from non-binder here, and that the
only compositional feature surviving multiple-testing correction is asparagine.
Their rows are kept, commented, in `make_figure.py`, and
`code/Fig5def_asn_interface` is untouched and still builds, so restoring them
is uncommenting two rows and re-lettering. Everything below about the old d/e/f
is kept for that.

## Panels

| panel | size, V2 | source |
| --- | --- | --- |
| a | 99.3 x 88.0 mm | `working/screen_matrix/screen_matrix.svg`, hand-built |
| b | 72.7 x 88.0 mm | `output/Fig5a_wetlab_on_off/Fig5a_PosAvg_vs_NegAvg.svg` |
| c | 177.9 x 44.5 mm | `code/Fig5c_best_hits` |
| d | 178.0 x 34.2 mm | `code/Fig5c_cutaway_surface` |
| ~~d~~ | 86.0 x 52.0 mm | removed: `code/Fig5def_asn_interface`, `metrics.svg` |
| ~~e~~ | 86.0 x 52.0 mm | removed: `code/Fig5def_asn_interface`, `aminoacid_volcano.svg` |
| ~~f~~ | 178.0 x 46.0 mm | removed: `code/Fig5def_asn_interface`, `asndist.svg` |

V2 now ends 200.7 mm down the page with 81 mm to spare, so there is room if
d/e/f come back or c and d want to be larger.

**a** The transfection matrix: RNA target scaffold and control plasmid against
the well grid, four rows (Target 1, Target 2, Neg. 1, Neg. 2) by the candidate
columns.

**b** The 32 ordered binders ranked by ON/OFF, avg Pos over avg Neg, coloured
by design stack, with the 0.48x marker at the foot and the 1.2 threshold
dividing inactive designs from proposed binders. `Fig5a_wetlab_on_off` emits
eighteen candidate plots; this is the one matched on the 0.48x marker.

**c** The four best hits by flow fold change, each on the refold model with the
lowest self-consistency RMSD. Twist 23, 15 and 16 are ranks 1 to 3; twist 7 is
rank 4, present for its interface asparagine rather than its rank. Dashes mark
contacts passing the ChimeraX `hbonds` test, with the measured distance. Only
twist 23 (TYR24 to the A11 phosphate, 2.59 A) and twist 16 (LYS42 to the U10
phosphate, 3.11 A) carry one. Full provenance in
`code/Fig5c_best_hits/README.md`.

**d** The same four designs as c, as molecular surfaces: RNA in neutral grey,
each binder in its stack colour, nothing else coloured. Rendered in c's own
per-panel cameras - the `turn y` angles, view specs and model choice are
imported from `code/Fig5c_best_hits/render.py`, not copied - so each surface
sits under the close-up of the same design. It carries no labels of its own,
because c labels the same four in the same order directly above it; the caption
must therefore say d is c's four in c's order. Provenance and the verify
markers are in `code/Fig5c_cutaway_surface/README.md`.

### Removed 2026-11-18, kept for restoring

**d (old)** The twelve metrics the field ranks designed binders on, as
P(binder > non) with bootstrap intervals. All twelve intervals cover 0.5.
ligand-RMSD and DockQ appear per backend, never as the better of the two, because
the refold used elsewhere is picked by lowest ligand-RMSD and testing the selected
value would test the variable that did the selecting.

**e (old)** All nineteen amino acids present in these designs (Cys is excluded by the
design protocol), same statistic, for sequence and for interface composition.
Asparagine is the only one that survives BH correction, q = 0.046 for sequence
and 0.064 for interface. Histidine is nominally significant and does not
survive, q = 0.20.

**f (old)** Asparagine across the 12k library, the flow inactive designs, the flow
binders and 59 natural RBPs, as a percentage of sequence and of interface.

## The caption has to carry this

From `Fig5def_asn_interface`, stated there as the honest shape of the result:
the binders are **not** unusually asparagine-rich. At 5.7% of interface
residues they sit close to the 600-design library sample at 4.6%. The
**inactive designs are the anomaly**, at 1.5%, with 71% of them having no
interface asparagine at all against 41% of the library. The result predicts
failure from a missing contact, not success from a clever one. A caption that
says "binders are asparagine-rich" inverts it.

## Open markers

- `VERIFY:` panel c's legend must state that a dash marks a contact passing the
  ChimeraX `hbonds` test. Without it the two dashless panels read as an
  oversight rather than as the result they are.
- `VERIFY:` that hbonds test ran on structures with no hydrogens, so
  constraints were relaxed by 0.4 A and 20 degrees and donor-H..acceptor angles
  were never evaluated. It is a donor-acceptor distance and geometry screen.
  Per `LV_writing/style/claim-boundaries.md`, do not word it as binding.
- `NEED:` panel letters a to f are placed by the assembly, but the figure has no
  caption text yet. Leo writes it.

## Provenance of this file

Panel c was re-verified against primary files on 2026-11-16: ranking and FC
from `LV_32Candidates/index.tsv`, model choice from
`inputs/csy4_binder_library_12k_v6_k30.csv`, contacts from a fresh ChimeraX run
kept at `code/Fig5c_best_hits/verify_hbonds.log`.

Panels a, b, d, e and f are described from their own scripts' docstrings and
the assembly's reported geometry. Their underlying numbers were **not**
independently re-checked in that pass. Treat the d/e/f statistics above as
quoted from `Fig5def_asn_interface`, not as verified here.
