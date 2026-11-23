# Fig3c_sse_ternary

Alpha / Beta / Loop ternary plots of backbone secondary-structure content, one column per
generator (RFD3, BoltzGen, BD1, HD), top row all backbones, bottom row backbones that passed the gate.
Replaces Fig3b_sse_cluster_map (left untouched). Input: `inputs/sse_stage.parquet`. Run: `~/.venvs/binder_map/bin/python code/Fig3c_sse_ternary/make_figure.py`.
Emits `output/Fig3c_sse_ternary/panel_c.{pdf,png}`, 177.8 x 62.0 mm.

- `stage` is exclusive (0 fail, 1 gate, 2 library): before = all rows, after = `stage >= 1`.
  The script asserts the eight counts (RFD3 144,015/14,312; BoltzGen 149,933/9,895; BD1 113,175/1,661; HD 91,941/444).
- loop = 1 - helix - strand. BG_BIF + BG_LiMPNN merged as BoltzGen (colour: BoltzGen + LigandMPNN).
- Density: 2-D histogram in triangle coordinates (140 bins across), hand-rolled separable Gaussian, filled at the 90% and 50%
  highest-density regions, clipped to the triangle. Sigma (bins) is the SAME for both rows of a generator and varies between
  generators by their smaller set: RFD3 2.5, BoltzGen 2.5, BD1 3.5, HD 6.
- **HD's passed-gate region is estimated from only 444 backbones and is the least reliable triangle in the figure.** It is drawn
  as a filled area for consistency (Leo's request) with the widest blur, and with the 444 raw backbones faintly overlaid so the
  sparsity is visible. n is printed on every triangle.
- Three nested fills: 50% (dark), 90%, 100% (lightest). "100%" = level of the least-dense DATA POINT on the smoothed grid, so every
  backbone lies inside its own 100% region. Dashed grey (after row) = the 100% region of that generator's unfiltered set, same sigma.
  Because passed-gate is a subset, nesting is guaranteed; verified: passed backbones outside the dashed line 0.00% and after-100% area
  outside 0.00% for all four generators. The 100% regions are large and somewhat ragged (set by the most isolated backbone; RFD3
  unfiltered has stray islands, BD1 fills most of the triangle) - the cost of an all-inclusive outline.
- Ticks (0.25, 0.5, 0.75) with numbers on ONE triangle only (RFD3, after row); every triangle carries Alpha/Beta/Loop vertex names.
  Alpha is read on the left edge, Loop on the right edge, Beta on the bottom edge.

Caption draft: Secondary-structure composition (helix = Alpha, strand = Beta, remainder = Loop) of designed backbones per
generator before and after gate filtering; shading encloses 50% and 90% of backbones; smoothing width increases as n falls
(HD after filtering, n = 444, is the least reliable estimate); dashed line, 90% region before filtering.


## Key (frameless, one line below the grid, no sentences)
Swatches "50%", "90%", "100%" and a dashed sample tagged "unfiltered". Meaning lives in the caption (captions/Fig3.md): each fill is the
smallest region holding that share of the generator's backbones; the dashed line is the same construction (100%) over the unfiltered set;
blur widens as n falls. Sigma (bins of 1/140): RFD3 2.5, BoltzGen 2.5, BD1 3.5, HD 6 - identical for a generator's two rows, so sharpness
is not comparable between columns. HD passed-gate (n = 444, drawn as a smoothed area with no points) is the least reliable estimate.
