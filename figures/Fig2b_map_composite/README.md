# Fig2b_map_composite

Assembly components for the campaign-map figure. The map is drawn here; the
zoom-ins are placed by hand.

| File | What it is |
|---|---|
| `map_only.{svg,pdf,png}` | The map and the eight numbered marks |
| `map_only_territories.{svg,pdf,png}` | Same frame and same marks, with each generator's density shaded and outlined |
| `../Fig2b_stack_binders/callout_*.png` | The eight renderings, one per file, transparent background |

170 x 170 mm. The two map variants share one frame, so they drop onto each
other. **Neither carries the stack key**: it is placed at assembly time with
the renderings. The only text is the one-line note on thinning at the foot,
which is a single text object if it needs to move or go.

## Numbers on the map to rendering files

| # | File | UMAP x, y | Fold |
|---|---|---|---|
| 1 | `callout_hd_a.png` | 1.18, 1.82 | 7 |
| 2 | `callout_bd1_a.png` | 4.14, 2.60 | 8 |
| 3 | `callout_rfd3_a.png` | 7.48, -0.27 | 5 |
| 4 | `callout_hd_b.png` | 7.46, 3.75 | 3 |
| 5 | `callout_bg_a.png` | 2.67, 4.04 | 1 |
| 6 | `callout_rfd3_b.png` | 4.66, 5.45 | 4 |
| 7 | `callout_bd1_b.png` | 4.14, 7.08 | 4 |
| 8 | `callout_bg_b.png` | 9.41, 9.55 | 14 |

Metrics for the eight are in `code/Fig2b_map_callouts/README.md` and
`inputs/map8_callouts.csv`; they are deliberately not drawn on anything.

All eight renderings share ONE crop, taken over the eight alpha bounding boxes
together, so relative size is real: a binder that is physically larger is drawn
larger. Scale them as a group when placing them, not individually.

## Three decisions

### 1. The dot cloud is thinned
499,064 marks in a 128 mm square is a solid field in which an 8-point callout
cannot be seen. The map is binned at 0.12 UMAP units and at most 6 backbones are
kept per bin, so every sparse arm keeps all of its points and only the saturated
masses lose any. **41,072 marks survive, 8.2%**, and the silhouette is unchanged
because the parts that were thinned were already opaque. The figure says so on
its face.

Drawing decision only. No count in this or any other figure comes from the
thinned set.

### 2. The dots go out as true vector
`rasterized=False`, which is only affordable because of the thinning: 41,072
circles is a 6.6 MB SVG that Illustrator opens, 499,064 is about 80 MB and it
does not. Every mark stays selectable and recolourable.

### 3. Territories are a separate file, not the default
At the settings used here (0.10-unit grid, gaussian sigma 3.0, contour at 5% of
each generator's own maximum) the outlined area is **not** proportional to the
number of backbones:

| Generator | Outlined area | Backbones | Area per 1,000 |
|---|---|---|---|
| BoltzGen | 129.1 | 149,933 | 0.86 |
| RFD3+LiMPNN | 40.6 | 144,015 | 0.28 |
| BD1+LiMPNN | 23.6 | 113,175 | 0.21 |
| HD+LiMPNN | 17.1 | 91,941 | 0.19 |

BoltzGen and RFD3 hold almost the same number of backbones and BoltzGen's
territory comes out three times larger per backbone. BoltzGen's footprint is
dozens of small, widely separated islands; the other three sit in compact
masses, and they agree within 15% of each other, which is the check that the
effect is specific to fragmented topology rather than general. Any smoothing
wide enough to outline a scatter of islands as one territory inflates it.

So the shading is a reading aid and **the dots underneath are the honest
record**. If `map_only_territories` is used, say so in the caption.

(Earlier numbers of 240.2 / 54.0 / 38.0 / 29.7, about five times, came from
sigma 4.0 at a 2% level, which is not what these files use.)

## Rebuild
```
python code/Fig2b_map_composite/make_figure.py
```
Reads `inputs/large/fold_map_scatter.parquet` and `inputs/map8_callouts.csv`.
The eight renderings are rebuilt separately:
```
python code/Fig2b_stack_binders/make_figure.py OUT --set map8 --width 3000 --angle 120
```

## Open
The map is at full extent, so the dense core occupies roughly the middle third:
about 9% of backbones sit in far-flung islands that take up most of the area.
Cropping to x -4..16, y -5..16 would enlarge the core by about a quarter and
drop 22,323 backbones (4.5%); not done, since `fig02b` and `fig05` draw the same
embedding at full extent.

`map_composite.{svg,pdf,png}` in the output folder is the earlier assembled
version, with the renderings and leaders baked in. This script no longer
produces it.

`map_composite.{svg,pdf,png}` in the output folder is an earlier assembled
version with the renderings and leaders baked in. This script no longer produces
it, so it will not update; delete it once the hand-assembled figure exists.
