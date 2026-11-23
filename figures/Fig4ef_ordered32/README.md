# Fig 3c: the 32 ordered candidates

Two outputs, both about the same 32 designs: the Echo order LV2.27
(`32_Binders_ON`), screened arrayed in LV2.33.

## Source of truth

`working/structures/hits32/`, copied from
`LV_BinderDesign/LV2.1.1_RNAbinders_PoolCraft_V1_old/12_analysis/hits32_browse/`
(32 PDB complexes plus `index.tsv`, which carries the comparison and the ON/OFF
ratio). The same 32 rows are the ones validated against
`gated_Q1_75Reads_OrderedHits.xlsx` when `inputs/screen_points.csv` was built.

The browse PDBs hold **three MODEL records each**. ChimeraX opens a multi-model
PDB as submodels, which makes `align #N/B:...` select nothing and the whole
`.cxc` abort with no error on stdout. The copies in `working/structures/hits32`
are therefore truncated to MODEL 1. These files have full sidechains, unlike
the backbone-only design CIFs used elsewhere.

## Stack distribution

    python3 code/Fig4ef_ordered32/make_pie.py     -> output/Fig4ef_ordered32/Fig4e_stack_pie.{png,pdf}

| Stack | Ordered |
| --- | --- |
| RFD3+LiMPNN | 4 |
| BG+LiMPNN | 11 |
| BG+BIF | 9 |
| BD1+LiMPNN | 7 |
| HD+LiMPNN | 1 |

Slice order and colour are `palette.STACK_ORDER`, labels are direct with leader
lines rather than a legend.

## 32-panel grid

    python3 code/Fig2b_stack_binders/make_figure.py working/structures/hits32_renders \
        --set hits32 --width 1200 --frames 12 --angle 120
    python3 code/Fig4ef_ordered32/make_grid.py PNGDIR [--rows 5 --cols 7 --cell 600]

The 32 renders live in `working/structures/hits32_renders` (re-made 2026-11-14;
the original PNGDIR was a scratch directory and had been lost, so the sheet
could not be rebuilt without re-rendering).

Panels use the Fig 2b shared RNA frame and the same y120 rotation as the Fig 3b
panels and the map callouts, so the hairpin is identical in all 32 cells and
only the binder moves. Binder colour is the stack colour.

Cells are grouped by stack in the fixed stack order, then by ON/OFF ratio, and
stack identity is carried by the binder colour, the caption colour and, with
`--boxes`, the tint behind the cells.

`--boxes` additionally tints each stack's cells as one connected shape (its
row-runs painted through a shared mask and bridged over their horizontal
overlap, so a stack that wraps rows is one shape, not a box per row) and
writes to `..._boxed.png`. Off by default, but the BOXED sheet is the one
the Fig 4 assembly places (since 2026-11-14).

Caption per cell is ONE line since 2026-11-14: the ON/OFF ratio, centred under
the cell in the stack colour. Its pixel size is solved from where it lands -
`--label-pt` (5.5) printed once the sheet is placed at `--page-mm` (178) - not
taken as a fraction of the cell, which is how the old caption ended up at
about 3 pt. `--names` restores the design name above it; it is off because at
eight columns it was unreadable, and the identifiers belong in a supplementary
table.

The stack key goes in the slots the 32 panels leave free - three of them at
7 x 5 - left-aligned on the first free slot's cell so the swatches sit on the
same column grid as the renderings. Frameless and at the ON/OFF label size
(STYLE_NOTES section 6); the swatches are the same 0.80 tint as the blocks
plus a hairline edge, which an 0.80 tint needs to be legible at swatch size
(COLORS rule 9). It is drawn only with `--boxes` and only when a slot is
free, so a layout with no gap simply has no key.

The default grid is 7 x 5, not 8 x 4: 25.4 mm cells instead of 22.2 at 178 mm,
and a 123 mm sheet instead of 86, which is most of the white Fig 4 had left
under it. 35 slots for 32, so the last row runs three short.

The two BG stacks are adjacent in the stack order and both palette colours are
purple, which is why `--boxes` is hard to read for those two.

All 32 panels are cropped to a window of ONE size - the union of all 32 ink
boxes, grown by `--margin` (default 6% of its larger side) - placed on each
panel's own ink centre. One size keeps the scale common, so binder sizes are
comparable across the grid and nothing is clipped (asserted per panel at
build time); centring per panel is what puts each rendering in the middle of
its stack tint.

Until 2026-11-14 the window was also POSITIONED in common, at the union box,
and slid 10% left by a `--shift` option. Both are gone. The hairpin sits at a
fixed place in every render and the binder extends to one side only, so a
union-positioned window is as wide as the widest binder: the median panel sat
11.7% of the window left of its cell centre, the compact ones further, and
the left column's hairpins were clipped by the shift. The cost of centring is
that the hairpin no longer lands at the same place in every cell - the shared
camera, rotation and scale, which is what the comparison rests on, are
untouched.

Layout is strict slot geometry. Every rendering owns one slot of identical
size, a slot being the rendering plus its caption band plus half a gutter on
each side, and the stack shapes are slots inset by the same amount on all four
sides. So every shape boundary falls exactly midway between two renderings,
every inner margin is equal, and the outer margin matches. The caption band
height is computed from the two font sizes rather than guessed, so nothing
clips at the bottom row.

Output: `output/Fig4ef_ordered32/Fig4f_grid_5x7_y120.png`.

## Open

The request said a "4x16" grid, which cannot hold exactly 32 cells. Built as
4 rows by 8 columns; `--rows` and `--cols` change it.
