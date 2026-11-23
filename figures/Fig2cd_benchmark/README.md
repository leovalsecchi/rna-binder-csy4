# Fig03 - campaign benchmark (quality and compute cost)

Per-stack results of the Csy4 hairpin binder campaign: the quality metrics as a
row of small multiples, the compute cost in a panel of its own. The claim is that
RFdiffusion3 leads on quality at a fraction of HalluDesign's compute, so cost is
deliberately NOT plotted against quality on a shared axis.

## Style
Follows the house conventions in `code/common/STYLE_NOTES.md`, derived from
the RFD4-Proteina paper (not in this repo - ask Leo for the PDF):

- **No chart titles.** The metric name and its unit are the y-axis label; panel
  letters (a, b, c) are added at assembly, not by the script.
- **Class by hatch, not colour.** Colour identifies the stack; diffusion is solid
  and hallucination is hatched, with a small frameless key. This frees colour from
  having to carry two things at once (the paper marks method type by marker shape
  for the same reason).
- **Bars carry the raw counts** (`22.8k / 1.52M` = GPU-hours / designs), not the
  plotted ratio, which the axis already gives.
- Sparse ticks, hairline axes, sentence case, no frames.

Axis labels must stay **short**: the metric panels are ~19 mm tall, so a rotated y
label beyond ~13 characters at 7 pt is clipped. `M` and `Rg` are expanded in the
caption instead.

## Panels

| File | What | Size |
|---|---|---|
| `metrics.*` | Five bar charts in one row: designs, ipTM, sc_DockQ, ligand_RMSD, Rg | 178 x 25 mm (double column) |
| `gpu_per_design.*` | GPU-hours per 1,000 designs, per stack | 85 x 40 mm |
| `gpu_per_design_log.*` | The same cost on a log axis, as dots - see below | 85 x 40 mm |
| `gpu_hours_bar.*` | Raw GPU-hours per stack. Alternative | 85 x 38 mm |

Panels of one figure, so the script writes to `output/` only and touches nothing
in `final/`. Assemble in Inkscape/Illustrator (source in `working/vector_sources/`)
and export the assembled figure to `final/Fig03_*`.

## Why the log version is dots and not bars
A bar encodes magnitude by its **length from zero**. On a log axis zero is at
minus infinity, so the baseline is arbitrary and bar length no longer encodes the
ratio - a bar chart on a log axis misstates its own data. The reference paper uses
log only on its **scatter** panels (a-d, designs per GPU-hour); its bar panel (f,
experimental hit rate) is linear from zero.

So `gpu_per_design.*` is bars on a linear axis and `gpu_per_design_log.*` is the
same numbers as dots, where position rather than length carries the value and log
is legitimate. The 5.5x spread compresses visibly on the log version - that is the
cost of the log scale, not a drawing error. Use whichever suits the argument; do
not combine them.

## Inputs
`inputs/campaign_summary.csv` - means over all designs, consensus = worst of the
AF3 and ESMFold2 refolds. **Transcribed from the campaign summary slide, not
recomputed from raw scores.** If Methods changes, change the CSV, not the script.

`designs_m` is rounded as shown on that slide, so the per-design figures inherit
that and are good to about +-1%. Exact per-stack design counts are not recorded;
the campaign total is 6,992,515.

## Charging the shared BoltzGen job
The 6,671 h BoltzGen generation job is shared by BG + LiMPNN and BG + BIF and
belongs to no single stack, so it stays in the CSV as a `kind = overhead` row and
`SHARED_MODE` (top of `make_figure.py`) decides how the compute panels charge it:

| `SHARED_MODE` | BG stacks | Reads as |
|---|---|---|
| `"each"` (default) | full 6,671 h to each | what each stack cost to run on its own |
| `"split"` | 3,335.5 h to each | shares of one real total, but neither BG cost is real |

`"each"` counts 6,671 h twice, so per-stack costs do not sum to a real total - fine
for a per-design comparison, not fine for anything presented as a share of the
campaign. `gpu_hours_bar.*` plots the hours as accounted and does NOT apply
`SHARED_MODE`. Neither mode equals the 197,797 h campaign total, which also
contains 13,749 h of self-consistency and other work.

## Colour
Palette as supplied: Shadow Grey `#161925`, Twilight Indigo `#23395B`, Rich
Cerulean `#406E8E`, Powder Blue `#8EA8C3`, Frozen Water `#CBF7ED`, assigned per
stack in `STACK_COLOR` and used in every compute panel, so a stack keeps its
colour throughout.

| Stack | Hex | Class (hatch) |
|---|---|---|
| RFD3 + LiMPNN | `#406E8E` | diffusion (solid) |
| BoltzGen + LiMPNN | `#23395B` | diffusion (solid) |
| BoltzGen + BIF | `#8EA8C3` | diffusion (solid) |
| BoltzDesign1 + LiMPNN | `#CBF7ED` | hallucination (hatched) |
| HalluDesign + LiMPNN | `#161925` | hallucination (hatched) |

It is a lightness ramp, not a set of distinct hues, so it does not meet the
categorical chroma check and never will - the steps share a hue by design. What
was validated and must hold is **separation**: worst adjacent pair dE 17.6 in the
bar order the script draws. Frozen Water is 1.13:1 against white, so every bar
also carries a label and a 0.4 pt hairline edge. Re-run
`scripts/validate_palette.js` before reordering or reassigning a stack.

The metrics panel is deliberately NOT stack-coloured: it keeps the Fig01 accent
(`#1F6F7A`) on RFD3 against neutral grey, since there the comparison is within a
metric, not across stacks.

## Caption material (belongs in the text, not drawn on the figure)
Suggested wording for the compute panel:

> Compute cost per design generated, in GPU-hours per 1,000 designs. Bar labels
> give GPU-hours and designs generated. BoltzGen's shared 6,671 h generation job
> is charged in full to both BoltzGen stacks, so their per-design figures are what
> each stack would cost run on its own. Solid bars, diffusion generators; hatched,
> hallucination generators.

For the metrics panel: **M**, millions; **Rg**, radius of gyration. ipTM and
sc_DockQ are higher-better, ligand_RMSD lower-better. Means over all designs;
consensus is the worst of the AF3 and ESMFold2 refolds.

Also for the caption or Methods:

- **BG + BIF** is scored against the matching BG + LiMPNN backbone (BoltzGen's
  true design). Values stored against BIF's own re-prediction were inflated
  (ligand_RMSD 29.6 -> 39.0 A, sc_DockQ 0.098 -> 0.042); 42,447 designs without a
  designed interface excluded.
- The campaign total is 197,797 GPU-hours including self-consistency runs.
- HalluDesign returned the fewest designs at the highest cost per design (82.9
  GPU-h per 1,000 designs against RFdiffusion3's 15.0, a 5.5-fold difference).
- BG + LiMPNN and BG + BIF share one BoltzGen generation job and the same
  backbones, which is why they track each other on every metric.

## Notes
Figure number is a placeholder: `code/Fig4c_fc_vs_splicing/` also claims `fig03`. Rename
this folder and the `Fig03_` stems to match the final thesis numbering.

`config/style.yaml` asks for Arial. The script falls back to the first available
of Arial / Helvetica / Liberation Sans / Nimbus Sans / DejaVu Sans rather than
letting matplotlib substitute silently, so rebuild on a machine with Arial before
submission if the typeface matters.
