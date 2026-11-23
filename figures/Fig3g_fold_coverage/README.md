# Fig3g_fold_coverage

**Terminology:** the panel's "backbone clusters" are the k30 Leiden communities on a Foldseek TM-align graph (mutual k = 30, TM >= 0.5),
column `cluster` in `inputs/fold_cluster_arms_k30.csv`. They are NOT folds in the structural-biology sense (nothing here tests whether a cluster
is a recognised fold), and they are a different partition from the `cluster_r10` clusters of the Figure 2 fold map: the two cannot be joined by id
(that mix-up has already produced one wrong figure). File/variable names keep "fold"; the reader sees "backbone clusters".

Standalone 178 x 45 mm full-width panel (above the pies and Fig3e); emits `output/Fig3g_fold_coverage/panel.{svg,pdf,png}` (live text). Does not touch `Fig3cd_pools`
(which keeps the small inset matrix). No title or letter; add at assembly.

Source: `inputs/fold_cluster_arms_k30.csv` only (917 backbone clusters). Do not mix with `fold_map_scatter.parquet` / `sse_arms.csv` (different clustering).

## What is drawn
**Matrix:** equal cells, one per backbone cluster (917), 71 columns x 13 rows, largest pool first, left to right then top to bottom (last 6 slots
empty). Pitch 2.49 mm (2.24 mm drawn, 0.25 mm gap). Grey = reached by both ipTM and ordered (326); warm = ordered only (591). No outline (the 50-cluster bar panel is no longer on the page). A hairline arrow in the top-right of the key row reads
"Largest -> smallest" (the row-by-row wrap is NOT drawn; say it in the caption). **Reading order, for the caption in one clause: cells run from the largest cluster
(top-left) to the smallest (bottom-right), left to right along each row and then down row by row.** Because pool sizes are ties among
singletons, order within the 661 singletons is arbitrary (stable sort of the CSV). The size gradient is visible from the pattern alone: the
first ~2.1 rows (163 clusters, all of size >= 6) are solid grey, then warm cells appear and take over (this IS the ipTM-reaches-every-cluster-of-6+ result).
The grey band shows no sense of scale by itself; the cue gives the direction, not magnitudes.

**Small print removed (Leo, 2026-11-17):** the line under the matrix ("657 pool backbones (2% of 31,763) lie in clusters only ordered reaches; each such cluster has
5 backbones or fewer") and the cue's parenthetical "(read along each row, then down)" are no longer on the panel. The 657 fact is the only statement that the 591
missed clusters hold almost no material, so it MUST move into the figure caption (pool backbones, not selected ones: the arms select only 7,813 / 9,759 unique
backbones). Earlier strip A (31,106 / 657 pool backbones) was also removed.

**Strip B (bottom):** the 917 backbone clusters as one bar: 326 grey (reached by both selections), 591 warm (reached by ordered only, 64%),
"+181% more clusters".

**Not drawn, for the caption:** ipTM reaches EVERY cluster of 6 or more backbones; the whole gap is clusters of 5 or fewer backbones.
By cluster-size class (clusters / ipTM reaches / ordered reaches / % more): 1: 661 / 120 / 661 / +451%; 2-5: 93 / 43 / 93 / +116%;
6-20: 9 / 9 / 9 / 0%; 21-100: 40 / 40 / 40 / 0%; >100: 114 / 114 / 114 / 0%; all: 917 / 326 / 917 / +181% (asserted earlier from the CSV;
a by-class bar version was tried and dropped as not what was asked for).

**Area-weighted treemap tried and REJECTED:** with cell area = backbones, the 591 missed clusters shrank to sub-0.5 mm cells
(0.46 x 0.42 mm) packed into one corner, below what prints, and the panel read as "the gap is negligible". Do not re-try it.

Asserted: pool 31,763 / ipTM 7,813 / ordered 9,759; ipTM touches 326 clusters (206 of 256 multi-member, 120 of 661 singletons);
ordered touches 917; ipTM misses 591 (541 singletons) holding 657 backbones (2.1%); no cluster is ipTM-only.

## Caption caveat (do not overclaim)
Ordered reaching 917/917 is partly **by construction**: a per-cluster round-robin touches every cluster given enough picks, so the rule is
effectively `ordered = min(cluster size, ~50)`. The defensible claim is about ipTM: it misses 591 clusters (64% of clusters), but every one
is small (largest = 5 backbones) and together they hold only 657 backbones (2% of the pool).

## Caption caveat: not a matched draw
Both arms select the same 12,000 designs but a different number of unique BACKBONES: ipTM 7,813, ordered 9,759 (about 25% more
distinct backbones for ordered to spread). The coverage comparison is therefore not a matched draw. The gap is far too large for that to
explain, and the misses are specifically small clusters, so the finding stands, but the caption must say it.
