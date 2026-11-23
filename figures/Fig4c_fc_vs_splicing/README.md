# Fig4c - enrichment vs splicing (PoolCraft screen)

Overlaid scatter of the RNA-binder splicing screen `LV2.1.1_RNAbinders_PoolCraft_V1`:
both comparisons in one panel, coloured by pair — **POS_A / NEG_A** (dark theme blue) and
**POS_B / NEG_B** (light theme blue). Clean framed style, no gridlines, no per-point labels.
Only designs with **reads >= 75 in both POS and its matched NEG** are plotted
(the NEG cut is what makes the ordered set 32); 

## Not a volcano
Called one until 2026-11-16, wrongly: a volcano plots effect size against
SIGNIFICANCE, and there is no p-value on this panel - it is log2 fold change
against the splicing ratio, two measurements, with the gate drawn on. Fig 5e
(`aminoacid_volcano`) is a real volcano and keeps the name.

## Axes
- **y = log2 fold change** of ON/OFF splicing, UNCAPPED (`fold_change_raw`): each
  design's splicing ratio in the POS (ON) condition divided by its matched negative
  control (OFF); when no NEG read was spliced, NEG = 1/NEG reads (the report's own
  pseudocount). Axis runs 0 to 8 (data max 7.5). The report's cap at 50
  (`matching_fold_change`) is display/export only and is no longer applied. The
  bottom is floored at log2 FC = 0: every design with FC <= 1 (no enrichment,
  including FC = 0 = no ON splicing) stacks on the bottom spine.
- **x = splicing ratio** = `removed_fraction`, the fraction of reads in which the
  hairpin was spliced out. This is the quantity the 0.07 gate is applied to.
  Clipped at 0.25; the 12 points beyond stack on the right spine.

## Gate (highlighting)
reads >= 75 (POS and NEG), splicing >= 0.07, FC >= 10. In-gate points are solid and coloured by
comparison; all other points are the same size but semi-transparent. Dotted lines mark
the FC and splicing thresholds; the reads cut is not an axis (noted on the
figure). Hits: 7 (A), 25 (B) = the 32 ordered (`candidates32_mapped.csv`). With the cut on
POS only it was 7 + 27: `RFD3_LiMPNN_a5dd135d_design27_seq0` (NEG 64 reads) and
`BD1_LiMPNN_a8e3107f_design57_seq4` (NEG 71) were never ordered.

## Inputs
`inputs/screen_points.csv` - one row per POS design (21,514), columns
`transcript, comparison, reads, splicing_ratio, matching_fold_change, log2fc,
fold_change_raw, neg_reads` (last two added 2026-11-12, other columns unchanged).

## Provenance (important)
The per-measurement source `matched_comparisons.csv` in the screen's
`12_analysis/` folder is **owner-only (chmod 700)** and was not readable. The same
data is embedded in `splicing_report.html` (objects `ALL_POINTS`, `NEG_CTRL_MAP`,
`TRANSCRIPT_META`); `extract_from_report.py` reconstructs the points using the
report's own `_computeFC` / `FC_CAP = 50` logic. **Validated**: the 32 rows of
`gated_Q1_75Reads_OrderedHits.xlsx` were re-derived and matched exactly on reads,
splicing ratio and fold change. If/when `matched_comparisons.csv` becomes
readable, regenerate `screen_points.csv` from it and confirm it is identical.

## Notes
- Figure number `fig03` is a placeholder; renumber to match the thesis.
- Font/raster: the committed PNG/PDF may use a fallback font if Arial is absent
  where they were rendered; `python3 make_figure.py` (or `run_all.sh`) regenerates
  them in Arial on a machine that has it. The SVG keeps text as text.
