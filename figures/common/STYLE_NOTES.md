# Figure style notes

House conventions for thesis/paper figures in this repo, derived from
**RFD4-Proteina** (Baker Lab / NVIDIA, Oct 2026) as the reference for what a good
figure in this field looks like. Read this before building a new panel.

The mechanical settings - font, sizes, DPI, column widths - live in
`config/style.yaml` and `code/common/style.py`. This file covers the decisions
those cannot encode. Colours: follow `COLORS.md` at the repo root (thesis palette v0.1).

> The reference PDF is **not** in this repo (55 MB). Ask Leo for
> `RFD4-Proteina.pdf`; the figures cited below are on pages 5 (Fig. 1), 7 (Fig. 2)
> and 13 (Fig. 4). Captions sit on the page *after* each figure.

---

## 0. Paper-like means the least ink that still reads
Every mark has to earn its place. A panel is read inside a figure with a caption,
not alone on a slide, so anything the caption can carry is **not drawn**. Before
adding any element, ask what breaks if it is absent; if the answer is nothing,
leave it out.

In practice, the things that creep in and should not:

| Creeps in | Where it belongs instead |
|---|---|
| Title, subtitle, standfirst | caption (section 1) |
| Gridlines | nowhere, almost always (section 7) |
| A legend box | direct labels (section 6) |
| A number on every mark | selective annotation (section 5) |
| "lower is better" as text | an arrow in the axis label (section 2) |
| A long axis label | short label + caption (section 3) |
| Explanatory notes on the canvas | caption |

This is the single rule the others are instances of. When a panel looks busy,
the fix is removal, not rearrangement.

## 1. No chart titles. Ever.
A panel in a paper has no title and no descriptive subtitle. What a title would
have said goes in the **caption**.

- The quantity and its unit are the **axis label**: `Clustered pass rate`,
  `Designs per GPU-hour`, `Experimental hit rate (%)`, `Motif RMSD (Å)`.
- Panels get a **bold lowercase letter** (`a`, `b`, `c`) in the left margin,
  outside the plot content, aligned to the panel top. No parentheses, no `a)`.
- At most a short bold topic label beside the letter, two or three words:
  `a  Protein-protein binding`. Never a sentence.
- Panel letters are added at **assembly** (Inkscape/Illustrator), not by the
  plotting script. Scripts emit bare panels.

## 2. Direction belongs in the axis label, as an arrow — but only where it is news
For a metric whose good direction is **not self-evident**, append an arrow:
`H-bond energy (TMol) ↓`, `Motif RMSD (Å) ↓`. Compact, and it survives being read
out of context. Do not write "lower better" as a subtitle.

**Do not arrow a quantity that explains itself.** `Pass rate ↑`,
`Designs per GPU-hour ↑`, `Yield ↑` tell the reader nothing they did not already
know and cost a glyph on every axis — exactly the kind of creep section 0 is
about. The test: would a reader hesitate about which end is good? Energies,
RMSDs and error rates yes; rates, counts, throughputs and scores no.

## 3. Keep axis labels short
Rotated y labels are bounded by the panel height. At the house 5.5 pt
(`config/style.yaml`; 7 until 2026-11-13) a label beyond ~17 characters is
clipped on a ~19 mm tall panel. Abbreviate on the panel (`Designs
(M)`, `Rg (Å)`) and expand in the caption. Check the rendered panel; do not assume.

## 4. Encode class by shape or texture, not by colour
When marks carry both an identity and a category, spend colour on the identity and
something else on the category. The reference keys method type by **marker
shape** (`● Generative  ■ Hallucination  ▲ Refinement`) in a single small key
shared by all panels, leaving colour free for the individual method. In bar
panels, hatch is the equivalent lever.

Series in line plots are distinguished by **colour and line style together**
(solid / dashed / dotted), so they survive greyscale printing and colourblind
readers.

## 5. Label the raw counts, not the plotted number
Where a bar shows a derived quantity, annotate it with the inputs:
`158/188`, `18/60` for a hit rate; `22.8k / 1.52M` for GPU-hours per design. The
axis already carries the derived value - repeating it on the bar wastes the slot
that could show what it came from.

Never a number on every point of a dense series; annotate selectively.

## 6. Direct labels beat legend boxes
Where there is room, name the series at the end of its line or beside its mark, in
the series colour, rather than in a legend box (reference Fig. 1e, Fig. 2a-d).
Italic qualifiers go on a second line: *RFD4-Proteina (high thinking)*.

Where a legend is unavoidable it is **frameless**, small, and placed outside the
plot area - right of the panel or along the bottom, shared across panels.

## 7. Marks and axes are recessive
- Hairline axes, ticks outward, sparse: decade ticks on log axes
  (`10 / 100 / 1,000 / 10,000`), round numbers elsewhere (`0% 25% 50% 75% 100%`).
- No top or right spine.
- **No gridlines.** The default is off and `code/common/style.py` pins
  `axes.grid: False`. They are allowed only when a reader must read a precise
  value off the panel, which for a scatter, a ranked bar or anything with direct
  labels is never. A comparison between marks is carried by position; a number
  that must be exact belongs in the caption or a table. If a panel does turn them
  on, it is one axis only, hairline, in `palette.GRID`, and behind the data
  (`ax.set_axisbelow(True)`).
- Markers ~3-4 pt, not 5+. Lines ~0.7-1 pt.
- Error bars: thin, capped, and the caption states what they are (mean ± 1 SD,
  n per group).
- Pale grey rounded containers group sub-panels; a pale vertical band can tie a
  stage across stacked panels (Fig. 1d-f).

## 8. Bars start at zero; log axes need point marks
A bar encodes magnitude by **length from zero**, so a bar chart cannot go on a log
axis - zero is at minus infinity and bar length stops encoding the ratio. The
reference follows this strictly: log axes appear only on its scatter panels
(Fig. 2a-d), while its bar panel (Fig. 4f) is linear from zero.

If a quantity spans decades and must be logged, draw it as **dots** - position,
not length, carries the value. See `Fig2cd_benchmark/make_figure.py`
(`gpu_per_design` vs `gpu_per_design_log`) for both forms of the same data.

## 9. Colour
- One restrained palette across the whole figure set. The reference uses roughly
  purple / teal / navy / gold plus greys, and reuses them consistently so a colour
  means the same thing in every panel.
- **Stack order is fixed: RFD3, BG+LiMPNN, BG+BIF, BD1, HD** — in every figure,
  left to right and top to bottom. Iterate `palette.STACK_ORDER`; never sort by
  value. See COLORS.md rule 2.
- **Colour follows the entity, not its rank.** A stack, method or series keeps its
  colour in every panel it appears in; it must not be reassigned because the
  sort order changed.
- Validate any categorical palette rather than eyeballing it. Lightness ramps
  within one hue will never pass a categorical check - what matters there is
  adjacent-pair separation plus a label on every mark.
- A fill below 3:1 against white needs a hairline edge and a visible label.

## 10. Two places the reference diverges from general dataviz advice
Noted so nobody copies them unaware:

- **Dual axes.** Fig. 1e carries a second y-axis (sequence recovery, %), coloured
  to match its series. General guidance says never; this field tolerates it when
  the axis is colour-matched to the one series that uses it. Prefer two panels.
- **Series-coloured label text.** The reference writes series names in the series
  colour. General guidance keeps text in neutral ink with a coloured mark beside
  it. Either is defensible; be consistent within a figure.

---

## 11. No panel letters in the artefact
A script never draws the a / b / c marker. Letters are placed at assembly, where
the arrangement is decided, so the same panel can be reused in another position
or another figure without being re-rendered. The panel's letter lives in its file
name (`Fig<N>_Pan_<L>_<description>`) and in the assembly document, never in the
SVG or PNG.

---

## Checklist before calling a panel done
0. Every element on the canvas is load-bearing; anything the caption could carry
   has been removed (section 0). No gridlines (section 7).
1. No title and no standfirst drawn in the artefact - schematics included, not
   just plots. Both belong to the caption, drafted in the figure's README. A
   number or test statistic annotating one panel is not a title and may stay.
2. Axis labels carry the unit, and a direction arrow ONLY where the good
   direction is not self-evident (section 2).
3. Nothing clipped - look at the rendered PNG, do not infer from the code.
4. Identity never rests on colour alone.
5. Bars from zero; log only with point marks.
6. Panel emitted bare, at its authored width, into `output/figNN_*/`.
7. Nothing written to `final/` by a script - assembly exports there.
8. No panel letter drawn in the artefact (section 11).
