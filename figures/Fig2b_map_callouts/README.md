# Fig2b_map_callouts (Figure 2b)

Campaign fold space as **one map**, dots coloured by design stack, in the same
shared UMAP as `fig05_fold_clustering`, so a position means the same backbone in
both figures.

Output is **two exactly registered layers**, 178 mm square, PNG at 1200 dpi
(8409 x 8409 px):

| File | Holds |
|---|---|
| `dots.*` | the point cloud and nothing else (zero text elements) |
| `labels.*` | axis cues, swatch key and caption, transparent background |

Stack them in Illustrator and you get the complete figure back; verified by
compositing the two PNGs, which reproduces the single-panel version pixel for
pixel. Keeping them apart means the dots can be moved, scaled or recoloured
without dragging the text along, and the text can be restyled without
re-rendering half a million marks.

**Registration is load-bearing.** The axes rect is fixed with `add_axes(RECT)`
and both layers save with `bbox_inches=None`. Do NOT reintroduce
`tight_layout()` or `bbox_inches="tight"`: those crop each layer to its own
content, and since the two layers have different content they would come out
different sizes and no longer line up.

| Colour | Group | Plotted backbones |
|---|---|---|
| `#a83902` | RFD3 + LigandMPNN | 144,015 |
| `#a376e9` | BoltzGen (LigandMPNN + BoltzIF) | 149,933 |
| `#bd8703` | BoltzDesign + LigandMPNN | 113,175 |
| `#d35182` | HalluDesign + LigandMPNN | 91,941 |

Nodes are unique backbones (coordinate-hash collapsed), **not** designs: the
6,992,515 campaign designs collapse to 551,516 distinct backbones because
LigandMPNN draws many sequences per backbone.

## Draft caption
> **Fold-space occupancy by design stack.** Each point is one unique protein
> backbone, positioned by UMAP over the Foldseek TMalign graph used for
> clustering (mutual k = 30, TM >= 0.5), coloured by the stack that produced it.
> The two BoltzGen stacks are shown as one colour: BoltzIF is an inverse folding
> step on BoltzGen backbones, so the two differ in sequence design, not backbone
> source. The 52,452 backbones with no surviving mutual k-NN edge are omitted,
> leaving 499,064 of 551,516. The stacks occupy largely distinct regions of fold
> space.

## Four things to know before using this
1. **The two BoltzGen stacks are merged by provenance, not by identity.** They
   are disjoint by coordinate hash: the 87 plotted BG+BIF backbones are not
   byte-identical to any BG+LiMPNN backbone, because repacking moves the
   coordinates. Merging says they come from the same generator, not that they
   are the same structures. The merged group keeps the BG+LiMPNN hex, the
   palette colour of 149,846 of its 149,933 members (99.94%); no new hex was
   invented. If a figure ever needs them apart again, the `stack` column in the
   parquet still distinguishes them.
2. **Every backbone belongs to exactly one stack.** Verified over all 7.0M
   design-to-hash rows in `results_bbmap/design/hashes/`: zero coord_hashes are
   produced by more than one stack, so the colour is unambiguous and no
   election rule is hiding a tie.
3. **The 52,452 degree-0 isolates are dropped** (9.5% of 551,516). They have no
   mutual k-NN edge surviving the TM 0.5 floor, so UMAP's random init scatters
   them around the perimeter. Coloured by stack at this size they form a
   near-solid ring that crowds the real structure inward and carries no angular
   meaning. Dropping them changes the per-stack totals, so the panel quotes
   plotted counts. Campaign totals including isolates: RFD3 152,512, BoltzGen
   151,996, BD1 144,468, HD 102,540.
4. **Draw order is randomised** (seeded, `default_rng(0)`). With half a million
   overlapping marks, whichever stack is drawn last otherwise dominates the
   visible colour. That is a plotting artifact, not a result. Dense overlap
   regions should not be read as one stack excluding another.

## Structure panels (the callouts around the map)

Three ChimeraX panels sit beside the fold map, all built from the AF3 refolds
or the design CIFs under `working/structures/fig02b_examples/`.

| Output | Shows | Design |
|---|---|---|
| `minorgroove.png` | binder tracking the backbone, reading the MINOR groove | RFD3 `8b13398e_design50_seq0` |
| `majorgroove.png` | binder reading the MAJOR groove of the stem | RFD3 `383063df_design6_seq0` |
| `pipi_spotlight.png` | TRP40 stacked on the flipped-out U10 | RFD3 `e2a4d6b9_design57_seq0` |

Built by `render_shared_frame.py` (groove panels) and `render_pipi_panel.py`
plus `make_spotlight.py` (the stacking panel). Head-binding was dropped on
request; `find_loop_stacker.py` and `analyse_examples.py` produce the numbers
quoted below.

### Things that are deliberate
1. **One shared RNA frame.** Both groove panels superpose onto the same
   reference hairpin and use one camera, so the RNA is in the identical
   position and only the binder moves. Comparing binding sites is the whole
   point; per-panel poses would make that impossible.
2. **Hairpin upright, open end down, side on.** The orientation is computed
   from the geometry (helix axis = first principal direction of the six
   Watson-Crick pair midpoints, signed so the loop is up), not eyeballed, and
   `turn y` then spins about that axis so the bottom stays at the bottom.
3. **The groove tint is chemistry, not a painted patch.** The surface is
   coloured from the groove-edge atoms themselves (major: N7/O6/N6/N4/O4/C5/C6;
   minor: N3/N2/C2/O2/O2'), the same definition `groove.py` uses to classify
   these binders. The tint is patchy because only parts of each groove edge
   reach the surface, which is a property of the surface, not a glitch.
4. **No re-prediction overlay.** Dropped on request, so these panels do NOT
   show design-vs-refold agreement the way the Fig03a filter panels do.
5. **Sidechains come from the AF3 refold, never the design.** Design CIFs are
   backbone-only (N, CA, C, O), so no sidechain, H-bond or stacking can be
   drawn from them. See `code/common/CHIMERAX_RENDERING.md`.
6. **The pi-stack is on the TOP flipped-out U10**, measured: U10 sits 13.0 A off
   the helix axis at the top, A1 13.6 A off it at the bottom, every stem base
   2-5 A. TRP40 stacks at 4.07 A centroid separation, 8.5 deg interplanar.
   An earlier candidate (`998adbc2`) stacked on G16/A1 at the BOTTOM instead.

## Rebuild
```
~/.venvs/binder_map/bin/python3 code/Fig2b_map_callouts/make_figure.py
```
Reads `inputs/large/fold_map_scatter.parquet` (see `inputs/MANIFEST.md`). The
point cloud is saved `rasterized=True`, so it embeds as a bitmap; a fully vector
SVG of half a million circles will not open in Illustrator. The labels layer has
no raster content at all, so it stays true vector and fully editable (4 KB SVG).
Nothing is written to `final/`.

## Callout layer (`callouts.*`, `map_callouts.*`)

`make_callouts.py` marks the eight highlighted binders from
`Fig2b_stack_binders` on this map, two per generator. All eight pass every term
of `passing_v6_final` and all eight are in the ordered 12k.

| # | Panel | Generator | Fold | UMAP x, y | ipTM | DockQ | ligand_RMSD |
|---|---|---|---|---|---|---|---|
| 1 | `callout_hd_a` | HD+LiMPNN | 7 | 1.2, 1.8 | 0.86 | 0.83 | 1.96 |
| 2 | `callout_bd1_a` | BD1+LiMPNN | 8 | 4.1, 2.6 | 0.82 | 0.87 | 1.51 |
| 3 | `callout_rfd3_a` | RFD3+LiMPNN | 5 | 7.5, -0.3 | 0.84 | 0.85 | 1.42 |
| 4 | `callout_hd_b` | HD+LiMPNN | 3 | 7.5, 3.8 | 0.77 | 0.86 | 1.80 |
| 5 | `callout_bg_a` | BoltzGen | 1 | 2.7, 4.0 | 0.92 | 0.90 | 1.64 |
| 6 | `callout_rfd3_b` | RFD3+LiMPNN | 4 | 4.7, 5.5 | 0.83 | 0.84 | 0.96 |
| 7 | `callout_bd1_b` | BD1+LiMPNN | 4 | 4.1, 7.1 | 0.68 | 0.86 | 1.66 |
| 8 | `callout_bg_b` | BoltzGen | 14 | 9.4, 9.6 | 0.85 | 0.87 | 1.11 |

Numbers match the structure panels in `output/Fig2b_stack_binders/callout_*.png`,
all rendered at y120 in the shared RNA frame.

`callouts.*` is a transparent layer in the **same frame** as `dots.*` and
`labels.*` (same figsize, same `RECT`, same limits, saved with
`bbox_inches=None`), so the three stack in Illustrator with no alignment work.
`map_callouts.*` is the same content flattened over the dot cloud at
alpha 0.40, for reading on screen.

### Three things that are deliberate
1. **Numbers sit against their marker, not on the perimeter.** Placement is
   solved, not typed: 24 directions at five radii, scored on how many backbones
   fall under the number, nearest-empty wins, with a minimum separation between
   numbers. A leader line is drawn only when a badge had to be pushed more than
   1.7 UMAP units clear. Hand-set offsets were tried first and needed re-tuning
   whenever a pick changed.
2. **No cluster hulls.** The eight sit in clusters of 8,813 to 34,804 backbones
   with x/y spread of 0.6 to 3.3 UMAP units, two of them in the same cluster
   (6 and 7 are both fold 4). A hull around those is a large blob that says
   nothing about where the binder is. The cluster id is in the key instead.
3. **Markers carry a faint keyline on their white halo**, so a mark looks the
   same over a dense patch as over white paper.

### Three of the eight sit in another generator's neighbourhood
Within 0.5 UMAP units of each marker, by generator:

| # | Own generator | Local composition |
|---|---|---|
| 2 | BD1+LiMPNN | RFD3 79%, BD1 13%, BoltzGen 4% |
| 5 | BG+LiMPNN | BoltzGen 50%, BD1 37%, HD 12% |
| 6 | RFD3+LiMPNN | BD1 81%, RFD3 18% |

The other five sit in neighbourhoods 97-100% their own generator. This is fold
convergence, not a plotting error: the marker is on its own backbone's
coordinate in every case.

### Rebuild
```
python code/Fig2b_map_callouts/make_callouts.py
```
Reads `inputs/large/fold_map_scatter.parquet` and `inputs/map8_callouts.csv`,
writes `callouts.{svg,pdf,png}` and `map_callouts.{svg,pdf,png}` at 1200 dpi.
The dot cloud is rasterized for the same reason as `dots.*`.

## `render_groove_ladder.py` - the two groove ZOOM panels (2026-11-17)

Not this figure's own three structure panels above: it re-renders two of the
EIGHT callout zooms that `output/Fig2b_map_composite/symmetric/zoom/` holds,
`zoom_hd_b.png` ("major groove") and `zoom_bd1_b.png` ("minor groove"), on
request - as `nucleotides tube/slab` ladders like `zoom_rfd3_b.png`, instead of
molecular surfaces.

With no surface to tint, the groove moves onto the slabs: each contacted
nucleotide is painted by the face the binder reads (petrol_800 major,
grey_800 minor, grey where it ties), classified per base type from
`analyse_examples.MAJOR/MINOR`. Both faces are drawn in both panels, so the
caption is checkable: hd_b is 11 petrol slabs to 1 grey, bd1_b 10 grey to 3.

That branch belongs in `zoom/make_zooms.py`, which builds all eight; this was
written as a separate script only because that file was out of scope. Re-run
`make_zooms.py hd_b bd1_b` and the surfaces come back (it also needs its
`sys.path` fixed first - it still points at `code/fig02b_fold_stacks`, the old
name of this directory, so it raises ImportError as it stands).
```
~/.venvs/binder_map/bin/python code/Fig2b_map_callouts/render_groove_ladder.py
```
