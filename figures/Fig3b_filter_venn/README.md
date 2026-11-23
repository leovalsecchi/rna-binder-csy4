# Fig3b_filter_venn

Four-set Venn of the filter criteria (ligand_RMSD, sc_DockQ, AA run, ipTM) with the
foot-clash step shown separately in a column to the right. Replaces `Fig3a_filter_funnel`.
178 x 60 mm. Run from repo root: `~/.venvs/binder_map/bin/python code/Fig3b_filter_venn/make_figure.py`.
Output: `output/Fig3b_filter_venn/Fig3b_filter_venn.{svg,pdf,png}`.

Data: `inputs/venn_regions.csv` (16 cells, sum asserted = 6,992,515; all-four
asserted = 80,983). Renders: `selection/Fig3a_box_{scrmsd,dockq,aarun,iptm,footclash}.png`.

## Caption notes
- **Areas are schematic, not proportional**: the 16 regions span 1,646 to
  2,901,047 designs, so counts are printed in each region.
- Layout is the standard four-ellipse Venn (two pairs of rotated ellipses).
- Foot clash is not a fifth set. The vertical bar at right is linear: of the 80,983 designs
  passing all four criteria, foot clash removes 25,277 (hatched), leaving
  55,706 (the campaign's passing set).
- The ellipse layout is stretched wider than tall (affine; topology unchanged).
- Renders are example designs; brackets mark the set each represents.
- Bottom-left number (617,398) = designs failing all four criteria.
