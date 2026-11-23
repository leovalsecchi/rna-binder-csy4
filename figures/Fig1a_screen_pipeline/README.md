# Fig1 panel A - screening pipeline

The route from designed binders to validated hits, as one row of seven steps.
Methods/overview figure. Output: `output/Fig1a_screen_pipeline/`.

## What changed on 2026-11-13
**The two phases now sit on grounds.** A neutral-100 band behind IN SILICO and
a theme-100 band behind EXPERIMENTAL, 12 units of white between them, no
stroke and no label of their own. The phase rules stay on top of them. Driven
by `PHASES`, so regrouping the phases regroups the bands; `BAND_FILLS` needs
one fill per phase and `_check` fails if it does not.

Why, after "no filled bands" was the 2026-11-10 decision: that call was made
on this panel standing alone, when a pale box here was the only one in the
figure. Panels b, d and f now all carry pale containers, so a tint per phase
is consistent with the page rather than against it. It is also the thing the
one-row layout could not say otherwise - the split was two words on a rule.

Chosen from six drafts in `working/Fig1a_drafts/` (gap, ground, divider,
bands, bands_soft, bands_named), built by `make_drafts2.py`. The richer
molecule glyphs drafted alongside them in `rich_glyphs.py` were **not**
adopted; the panel keeps the shipped one-colour line glyphs.

**The measured caveat, kept in the script.** Neutral 100 is relative luminance
0.913 against theme 100's 0.863 - a contrast ratio of 1.055, about 5%. The
split is carried by hue and will not survive a greyscale print. Neutral 200
is 1.066, barely better and the other way round. For greyscale, neutral 300 or
leave the dry half white.

Panel is 178 x 33.8 mm (was 178 x 32.8): the bands need air under the labels.

## What changed on 2026-11-11
Re-proportioned. It was a 178 x 21.8 mm letterbox slit whose labels printed at
4.3 pt, because the sizes were hand-picked viewBox units against `W=1600` and the
SVG was only ever looked at zoomed in.

- **Sizes are written in points now**, converted by `U` (viewBox units per point,
  set by the panel width). `LABEL_PT` comes from `style.yaml`, so the labels are
  the house 7 pt by construction and cannot drift again.
- `W` 1600 -> 1000, so a unit is 0.178 mm. Box is now 178 x 32.8 mm.
- Glyphs ~11.5 mm printed (`ICON_K` 1.45 -> 1.34 against the smaller W), chevrons
  and hairlines scaled with them. The row is dense instead of airy.
- `python3 make_figure.py --check` asserts the label size, that no label line is
  wider than its slot, and that the box is not a strip again. Nothing is written.

A two-row 4+3 variant (178 x 71 mm, serpentine with a wrap arrow) was built and
rejected on 2026-11-11: it reads as two sweeps rather than one pipeline. The
third option discussed was a quantity caption under each step (designs made,
designs passing, oligos ordered, hits); not done, the numbers were not confirmed.

## What changed on 2026-11-10
Reworked from the earlier `panel_a.svg` (kept in
`archive/fig01_pipeline_prerework_2026-11-10/`):

- **No filled bands.** The four pale rounded boxes were the slide-deck tell and
  carried nothing the step labels did not. Phases are now a hairline above the
  row with the phase name set into a gap in it.
- **Two phases, not four.** `IN SILICO` (steps 1-2) and `EXPERIMENTAL` (3-7).
  `PHASES` near the top of the script is a list of `(name, step indices)`; put
  the four back by editing that one list, nothing else depends on it.
- **One accent.** Every glyph is `palette.MAIN`; the stray orange is gone (rule 1).
- **Shared glyphs.** The seven placeholder glyphs defined inside the script were
  replaced by `code/common/glyphs.py` entries, so they are reusable and recoloured
  by the parent. Five are new: `binder-rna`, `funnel`, `integration`,
  `long-reads`, `cytometry`; `oligo-pool` and `tube` already existed.
- **Chevrons, not arrows** between steps, in neutral 400.
- `BxB1` corrected to `Bxb1`.
- **No panel letter** in the SVG. The `a` the old version drew is now added at
  assembly (STYLE_NOTES.md section 11), so the panel can be placed anywhere.

## Steps
| # | Glyph | Label |
|---|---|---|
| 1 | binder-rna | RNA-protein binder generation |
| 2 | funnel | Selection of favourable designs |
| 3 | oligo-pool | Synthetic DNA pool ordering |
| 4 | integration | Bxb1 integration |
| 5 | tube | RNA extraction |
| 6 | long-reads | Long-read NGS and analysis |
| 7 | cytometry | Flow cytometry validation |

Step 6 is drawn as reads tiled under the construct rather than as a nanopore: a
pore glyph did not survive the size, and the alignment also covers "and analysis".

## Inputs
None. Schematic only, no data read from `inputs/`.

## Icons
Written as standalone SVGs to `output/Fig1a_screen_pipeline/icons/` for
reuse in slides. They use `currentColor`, so the parent sets the colour.
`write_icons` deletes icons whose glyph no longer exists, so the folder never
offers a stale name.

Matteo's `0_Illustrator_Common_Cartoons` set was checked first: it is molecular
cartoons (hairpin, binder monomer/dimer, plasmid and landing-pad maps), with no
lab-workflow glyphs, so these were drawn fresh. `integration` borrows its
donor-above-locus arrangement from the landing-pad panels there.
