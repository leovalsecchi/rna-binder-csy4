# Fig02 - design pipeline

Schematic of the route from the RNA target to the ordered oligo library.
Methods figure (Notion figure plan row F2).

## Caption material
The figure carries no title or standfirst of its own (STYLE_NOTES.md, checklist
item 1). The line the artefact used to draw, kept here for the caption:

> Four backbone generators in two mechanistic classes, two inverse folders,
> three re-predictors, one gate.

## What it shows
Target, then the four backbone generators split into their two mechanistic
classes, then inverse folding, re-prediction, the filter cascade, and the
ordered library.

The generator split is the claim worth checking: RFdiffusion3 and BoltzGen
sample from a learned denoising prior (diffusion), while BoltzDesign1 and
HalluDesign optimise through a pretrained structure predictor (hallucination).
BoltzDesign1 inverts Boltz-1 but applies a stop-gradient to its diffusion
module, so it is not diffusion-based design.

## Panels
Single assembled figure, no manual assembly step, so the same artefact is
copied to `final/ by hand if it is ever promoted` as well as living in `output/`.

## Inputs
None. The figure reads no data from `inputs/`.

The three counts are hardcoded from Methods 2.9:

| Number | Meaning |
|---|---|
| 6,992,515 | designs at generation freeze |
| 55,706 | surviving the five-criterion gate |
| 12,000 | ordered panel after clustering |

If those change in Methods, change them here too; they are not derived.

## Icons
The glyphs come from `code/common/glyphs.py` and are also written as
standalone SVGs to `output/Fig2a_pipeline/icons/` for reuse in slides and other
figures. They use `currentColor`, so colour is set by the parent.

Ribbon glyphs (protein, diffusion, hallucination, inverse-folding, refolding)
are two-tone fills and need roughly 80 px to read. The hairpin is stroke-based
and survives down to 24 px.

## Colour
The semantic palette (teal for diffusion, amber for hallucination, indigo for
sequence design and selection) is defined in `make_figure.py`, not in
`config/style.yaml`, because that file still carries the placeholder matplotlib
colours. Worth promoting once the thesis palette is fixed.

The library step uses `oligo-pool` (staggered single strands, brightness varies
per strand). The older `tube` glyph is still in `glyphs.py` but unused here.
All three funnel clouds use one dot radius (`DOT_R`); only the count shrinks.

## The hairpin
The target slot holds a **real render of the 16-nt Csy4 crRNA hairpin**, not a
glyph: `render_hairpin.py` draws the campaign target itself (chain RNA of the
Fig2b shared-frame reference `csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0`)
in the same canonical pose as the Fig2b panels, loop up and open end down.

It is styled to sit among the glyphs rather than beside them: `lighting flat`
(ambient only, no specular, no shadow), `nucleotides ladder`, backbone
`thesis_petrol_700` which is the colour the old glyph used, bases
`thesis_petrol_300`, silhouettes `thesis_petrol_900` at width 3, transparent
background.

```
python code/Fig2a_pipeline/render_hairpin.py      # variants into output/.../hairpin/
python code/Fig2a_pipeline/make_figure.py         # embeds the pick
```

### Getting the bases horizontal
The camera is orthographic (`camera ortho`), so the rungs stay parallel instead
of converging: a glyph has no vanishing point.

The pose was chosen by measurement, not by eye. Sweeping tilt and azimuth and
scoring each by the mean angle of the six Watson-Crick C1'-C1' vectors off
screen-horizontal, the flattest view of this hairpin is **x -5, y 30**, where
the rungs still sit about **21.5 degrees** off horizontal. That residue is the
molecule, the A-form twist plus the hairpin's own bend, not a camera fault: no
single camera makes all six pairs horizontal at once.

`--roll 21.5` spends that 21.5 degrees as an in-plane roll, which levels the
rungs and tips the stem by the same amount. The figure uses the rolled view,
because horizontal rungs are what make the hairpin read as a hairpin at glyph
size.

`--roll auto` computes that levelling roll per azimuth, so a whole rotation
sweep comes out with level rungs and only the face changes. The 24-step sweep
is kept in `output/Fig2a_pipeline/hairpin/options/` with its contact sheet:

```
python code/Fig2a_pipeline/render_hairpin.py --styles ladder --tilts=-5 \
    --roll auto --angles 0,15,30,...,345
```

`render_hairpin.py` writes `ladder` and `slab` (the house tube/slab) and takes
`--tilts`, `--angles`, `--roll` and `--spins` (note `--tilts=-5`, since a leading minus
needs the equals form). `--spins` turns the hairpin about the rung axis, so it rotates THROUGH the
base-pair plane and the rungs stay level while the faces come round. Sweeps
kept beside the pick: `options/` (rotation about the hairpin axis),
`spins/` (rotation through the base-pair plane), `v345/` (the two crossed,
around s345).

The pick is `hairpin_ladder_x-005_y000_r+14.8_s345.png`, ladder bases, which
hold up better than slab at glyph size, named in `make_figure.py` as `HP_PNG`; change that
constant to switch. The render is
embedded base64 so the SVG stays one file, and `make_figure.py` falls back to
the old `rna-hairpin` glyph if the PNG is missing.

### Angle labels before 2026-11-14
Turns in a `.cxc` are relative, and the frame loop applied only the step
between angles, so the FIRST azimuth of a run was never turned to. A run whose
angle list did not start at 0 came out one step low while being named for the
angle asked for: `--angles 15,30,45` rendered 0, 15, 30. Fixed, and the files
already in `v345/` and `spins/` were renamed to the angle they actually show.
The `options/` sweep started at 0 and was always right.

This is why the figure's pick is `y000` with a roll of 14.8: the roll was
computed for the requested 15 and the frame came out at 0. The pose is the one
chosen off the sheet, so it is kept as it is rather than "corrected" into a
different picture.
