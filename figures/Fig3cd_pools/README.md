# Fig3cd_pools

Replaces panels c and d from `Fig3d_selection_difference` (that directory is untouched).
Emits `output/Fig3cd_pools/panel_c.{svg,pdf,png}` (108 x 45 mm) and `panel_d.{...}` (90 x 38 mm).
Text is live (`svg.fonttype: none`). No titles or letters; add at assembly.

## Panel c: 50 largest backbone clusters (178 x 44 mm, full width)
Column = backbones of that backbone cluster cluster in the gate-passing pool. Two narrow bars inside it, side
by side: left (cool) = taken by the ipTM-ranked selection, right (warm) = in the diversity-optimised
12,000 that was actually ordered. Largest cluster left; x axis is the cluster id, every column
labelled (checked with `text_w`: 1.25 mm clear between labels, no thinning needed). Insets are
1.19 mm wide (slot 3.31 mm). Class is carried by colour and by position; hatch is not practical at
that width. y axis 0-800.

Coverage matrix (upper right): one cell per backbone cluster, all 917, largest pool first, read left to right
then top to bottom (46 x 20 cells, last 3 empty); the outline marks the 50 drawn. Grey = reached by
both selections (326 backbone clusters), warm = reached by ordered only (591). Nothing is ipTM-only (asserted).
Coverage, asserted from `fold_cluster_arms_k30.csv`: ipTM touches 326/917 backbone clusters (206/256 multi-member,
120/661 singletons) and misses 591, of which 541 are singletons; ordered touches 917/917.
The 591 missed backbone clusters are all SMALL: the largest has 5 backbones and together they hold only 657
backbones. The coverage gap is about BACKBONE-CLUSTER DIVERSITY, not backbone count.
Defensible statement: ipTM ranking reaches 326 of 917 backbone clusters; the backbone clusters it misses are almost all
singletons and none has more than 5 backbones. Diversity ordering reaches every backbone cluster, by
construction, taking the whole backbone cluster where the backbone cluster is smaller than its per-backbone cluster cap
(ordered = min(cluster size, ~50)). Ordered reaching 917/917 is partly BY CONSTRUCTION (a backbone-cluster
round-robin touches every backbone cluster given enough picks): the claim is about ipTM's misses, not ordered's
perfection. The 50 drawn columns total 17,867 pool, 4,104 ipTM and 3,548 ordered backbones.
Trendlines and a cumulative-share inset were tried and removed.

Source: `inputs/fold_cluster_arms_k30.csv` (917 rows), the single k30 clustering the library
selection acted on (mutual k = 30, TM >= 0.5; Leiden on the Foldseek TM-align graph). Totals are the
headline ones: pool 31,763, ipTM 7,813, ordered 9,759, all asserted, as is ipTM and ordered <= pool
in every cluster. The ordered arm is a backbone cluster round-robin, so it sits near 50-60 in most clusters
while ipTM varies.

This file CANNOT be joined by cluster id to `fold_map_scatter.parquet` or `sse_arms.csv`: those
carry different partitions and embeddings. Panel c uses neither.

## Panel d: stack mix of three pools (90 x 40 mm)
Pies from `inputs/stackmix_k30.csv`: Passed the gate (55,706), ipTM-ranked (12,000),
diversity-optimised (ordered) (12,000; this is the `library` column). Totals asserted. Radii: pool
10.0 mm, the two selections 8.3 mm each (equal); the pool pie is not area-scaled to its total. The
selection labels sit to the right of their pies, left-aligned on a common x ("Diversity-optimised"
wraps to two lines); the pool label is beneath its pie because the arrows leave from its right
side. Slices clockwise from 12 o'clock in `palette.STACK_ORDER`; the five-stack key is one shared
horizontal row at the foot (78.4 mm of 90 mm). Emitted at 90 mm wide on purpose: the page must
place it at about scale 1.

## panel_picks (85 x 45 mm): picks per backbone cluster against cluster size
One point per backbone cluster per arm, all 917: x = pool backbones in the backbone cluster (log), y = backbones that arm
selected. Circles = ipTM-ranked, triangles = ordered, alpha 0.4 because integer sizes overplot.
ZEROS: y = 0 cannot sit on a log axis, so the 591 ipTM zeros are drawn in a separate lane below a
dotted break, labelled "0" on the axis, with a small vertical jitter inside the lane for visibility
only; none is clipped onto the axis as a small positive value. Ordered has no zeros (asserted).
The thin grey diagonal is y = x, labelled "whole cluster taken" along the line.
What it shows: ordered = min(cluster size, ~50). That one sentence explains both the y = x segment
(small backbone clusters, up to about 40 backbones, taken whole) and the flat band of roughly 20-60 per backbone cluster
on large backbone clusters; it never reaches zero. ipTM's zeros occur only in small backbone clusters (pool <= 5, mass at
1-2); on large backbone clusters it takes 5-230 per backbone cluster and overlaps ordered. Same honesty note as above:
ordered never reaching zero is partly by construction.
