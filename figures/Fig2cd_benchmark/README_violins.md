# Fig03 metrics — violin version

`make_violins.py` → `output/Fig2cd_benchmark/metrics.{svg,pdf,png}`

Replaces the bar version, archived at
`archive/fig03_campaign_benchmark/metrics_bars_2026-11-05.*`. The bars showed one
transcribed mean per stack; these show the full per-design distribution behind
each of those means.

## Data
`inputs/large/campaign_metrics_full.csv.gz` — all 6,992,515 designs, built by
`export_campaign_metrics.py` from the Euler store (also mirrored on the NAS at
`users/lvalsecchi/euler_work/rna_binder/analyses/`). Consensus is the worst of
the two refolds, verified to reproduce every per-stack mean on the campaign
summary slide:

| metric | definition |
|---|---|
| ipTM | min(af3, esm) |
| sc_DockQ | min(af3, esm) |
| ligand_RMSD | max(af3, esm) |
| Rg | mean(af3, esm) |

BG+BIF's sc_DockQ/ligand_RMSD are the true-backbone rescore, as on the slide; see
the note in `inputs/MANIFEST.md`.

## Layout (changed 2026-11-10)
Horizontal. Five stacks across five panels meant 25 tick labels set at 45
degrees; on their side the stack names are written once, horizontally, in the
first panel. `STACK_ORDER` reads top to bottom (rule 2), so the y axis is
inverted - row 0 is the top, not the bottom. Metric names are left-aligned panel
titles rather than rotated y labels, violins carry no outline, and the summary
is a white median dot on an interquartile bar instead of a median rule. The
count panel labels its bars directly and has no x axis.

## Four things to know
1. **Panel 1 is a bar, not a violin.** "Designs (M)" is a count — there is no
   distribution behind it, and a violin would invent a spread that does not exist.
2. **sc_DockQ is strongly zero-inflated**: 61% of all designs sit below 0.05 and
   81% below 0.1, while only 0.33% reach the 0.7 gate. The violin is therefore a
   narrow body with a long thin tail — that shape is the finding, not a plotting
   artefact, and it is why DockQ is the step that empties the funnel (fig07).
3. **Axes are clipped to the 0.5–99.5 percentile window.** ligand_RMSD reaches 155 Å
   and the KDE tails run further; without the window every body flattens to a
   line. Medians are drawn from **all** rows, not the clipped window.
4. **Quartiles and medians come from all rows**, only the KDE is sampled.
5. **The KDE uses a 60,000-design sample per stack** (`NSUB`, seeded); medians and
   every number quoted use all 6,992,515 rows. Raising `NSUB` changes nothing
   visible and costs time.

## Colour
Strictly `config/palette.json` via `common.palette` — no hand-picked hexes
(COLORS.md rule 1). Each stack carries its fixed colour in the canonical
`STACK_ORDER` (rule 2). This **changes the figure from the bar version**, which
highlighted RFD3 in a teal accent and greyed the other four: that predates the
palette and conflicts with stack colours meaning stack identity. The x order
changes with it, to RFD3, BD1, HD, BG+BIF, BG+LiM.

x tick labels double as the direct labels rule 5 asks for, since gold, rose and
light purple converge in greyscale.

## ⚠️ The rest of the repo is not yet on the palette
Only this script and `_template` import `common.palette`. Every other figure
still hard-codes hexes and predates `palette.json` (v0.1, 2026-11-05):

| script | hard-coded colours |
|---|--:|
| Fig1a_screen_pipeline, Fig2a_pipeline | 8 each |
| Fig2cd_benchmark/make_figure.py | 7 |
| Fig4c_fc_vs_splicing | 7 |
| Fig1f_scaffold_splicing | 3 |
| fig05_fold_clustering | 6 |
| fig06_fold_diversity | 1 |
| fig07_filter_funnel | 3 |
| fig07b_funnel_shape | 18 |

So `metrics.*` and `gpu_*.*` in this same folder are now in **different
palettes**. They should not be shown side by side until the others are migrated.
EOF
echo written