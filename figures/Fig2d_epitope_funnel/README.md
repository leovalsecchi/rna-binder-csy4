# Fig2d - what the funnel did to the epitope

    PY=~/.venvs/binder_map/bin/python3
    $PY code/Fig2d_epitope_funnel/make_columns.py   -> output/.../epitope_funnel_columns.{svg,pdf,png}
    $PY code/Fig2d_epitope_funnel/make_summary.py   -> output/.../epitope_funnel_summary.{svg,pdf,png}

Panels **e** and **f** of the assembled page (`code/Fig2_assembly`, bottom
row: the columns take the width they need, the square U10 panel rides beside
them at the same height). Letters are by position; see the clash note below.

    $PY code/Fig2d_epitope_funnel/make_u10_rise.py  -> output/.../u10_rise.{svg,pdf,png}

**f** answers the one question e makes you hunt for: U10's share against the
other fifteen nucleotides, 0.20 -> 0.31 -> 0.36 -> 0.39. Two bands stacked to 1,
because U10 and the rest are the whole; free-floating lines would invite the
reader to wonder what the gap between them means. Per nucleotide the gap is
wider than the bands look - the other fifteen average 0.053 -> 0.041 each, so
U10 runs from **3.9x to 9.5x the average nucleotide**. Put that in the
caption; the panel does not draw it.

Two forms of one result. `make_columns.py` is the panel Leo picked for 2d:
per-nucleotide contact share, campaign / gate / library side by side, each
column split phosphate / sugar / base. `make_summary.py` is the same result
compressed to five lines on three stages, for when 2d has to be narrow.

`--layers campaign` draws the campaign on its own
(`epitope_columns_campaign.*`); the default is all four, including the 32
ordered candidates. Four groups is where the layer labels start to overlap -
Library and Hits differ by 0.025 at U10, under one line height on a 45 mm
panel - so the label y is nudged apart while each leader still lands on its
true bar top.

## ⚠️ Panel letter clashes with Fig2cd_benchmark

`Fig2cd_benchmark` already emits **c** (quality metrics) and **d** (GPU cost).
Putting this at 2d means the benchmark's cost panel has to become 2e, or this
one takes another letter. Nothing is renamed yet - `Fig2cd_benchmark` is
untouched, and its README still describes itself as c and d.

## Geometry is set by the assembled page

`make_columns.py` authors at 178 x **45** mm, not 52: with the stem/loop
brackets below the axis the rendered ink comes out about 1.2x the authored
height, and at 52 the Fig 2 page overflowed by 3 mm. At 45 the page closes at
276.4 of 282 mm with every label still at T = 5.5 pt and no scale-to-fit.
If anything above e on that page grows, shrink this number rather than
letting the assembler scale the page, which would take the type off T.

Layer labels sit in a right-aligned column left of the U10 group, each joined
to its own bar by a hairline leader. Stacked above the bars (the obvious
placement) the campaign and gate labels run straight through the taller bars
beside them. The leaders are drawn by `annotate`'s own `arrowprops`, which
starts them at the text's bounding box, so the gap is right without measuring
the string (see [[figure-measure-dont-estimate-text]]).

## The result

| Layer | Nodes | Mean contacts | U10 | Loop (nt7-11) | phos / sugar / base |
| --- | --: | --: | --: | --: | --- |
| Campaign | 497,007 | 199 | 0.205 | 0.397 | 0.364 / 0.267 / 0.369 |
| Gate | 31,760 | 117 | 0.308 | 0.455 | 0.371 / 0.180 / 0.449 |
| Library | 9,757 | 123 | 0.362 | 0.505 | 0.372 / 0.175 / 0.453 |
| Hits | 32 | 134 | 0.387 | 0.553 | 0.369 / 0.165 / 0.466 |

### ⚠️ The 32 hits are NOT a fourth step of the funnel
The hits column continues the trend in DIRECTION only. At n = 32 the bootstrap
CI is +/- 0.024 per nucleotide on average, and the library mean falls outside
it at **2 of 16** nucleotides (G13, A15 - both near-zero shares), against 0.8
expected by chance. U10 0.362 -> 0.387 sits comfortably inside. So the panel
says *the screen hits read the library's epitope*, not *the screen sharpened
it further*. Say exactly that in the caption; the wide whiskers on the fourth
bar are the point, not a blemish. This matches the independent AF3-refold
measurement in `Fig4g_hit_epitope` (also 2 of 16, different two).

All 32 resolve to 32 distinct coordinate-hash nodes, all flagged `in_12k` and
`in_passing60k`, so the four layers nest properly. Resolution is
design_id -> coord_hash, falling back to the backbone representative
`struct_name` only where the design_id is absent from `hashes/*.tsv` (9 of 32);
the two differ for 2 designs that resolve both ways, and design_id wins because
the fingerprint measures the ordered design's own coordinates.

U10 nearly doubles its share of contacts, the loop goes from 40% to 51%, base
contacts go from 37% to 45% - and **all of that comes out of the sugar**
(0.267 -> 0.175) while **phosphate does not move** (0.364 -> 0.372). Designs
that gripped the 3' end are removed outright: G16 0.048 -> 0.003,
A15 0.025 -> 0.007.

**None of this was selected for.** The gate is ipTM / ligand_RMSD / DockQ /
seq_run / foot-clash ([[rna-binder-real-campaign-gate]]) with no epitope term
anywhere in it, so the epitope sharpening is a by-product - most likely of
DockQ, which is the filter that actually bites (it keeps 5.0% of what reaches
it).

In `make_summary.py` the two families are kept apart by colour temperature and
line style: warm solid for *where* on the hairpin (U10, loop), petrol dashed
for *what* is touched (base, sugar, phosphate), reusing the same three tones
as the stacked columns. The two families share an axis because both are
shares of one denominator, but they are not mutually exclusive - U10 is inside
the loop, and base/sugar/phosphate sum to 1 over the whole hairpin. Say so in
the caption.

## Source and its two caveats

`~/binder_filter_store/` (all 6,992,515 designs) has **no per-nucleotide
column** - only `fbase` and the `foot*` summaries. These numbers come from the
binder_map union run on Euler,
`foldseek_map/results_union/design/fp/canonical.parquet`: 526,467 nodes x a
48-dim fingerprint (16 nt x 3 moieties, heavy atoms within 4.5 A), joined to
`meta/layers.parquet` for the population flags. Aggregated by
`/cluster/work/wicky/lvalsecchi/epi_agg2.py` -> `epitope_layers_hits.csv`,
copied here as `inputs/epitope_layers_union.csv`. (`epi_agg.py` is the earlier
three-layer version; epi_agg2 reproduces its three rows bit-for-bit and adds
the hits row, verified on copy-down.) The 32 design ids go up as `h32.tsv`.
Run it with `/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python3` -
Euler's module python has no pyarrow. The `nt<NN>_lo` / `nt<NN>_hi` columns
are a bootstrap 95% interval over nodes. **They are no longer drawn** - Leo
asked on 2026-11-18 for "just the mean" - but they are still in the CSV, and
`make_columns.py` says in a comment what restoring them costs (two lines).

1. **Nodes are backbones, not designs**, and the campaign layer is the
   499,999-node sample, not all 6,992,515. Unbiased at node level, NOT at fold
   level - never scale fold counts off it.
2. **Design layer, not a refold.** Design-model coordinates, so the library
   reads U10 = 0.362 here against 0.231 on AF3 refolds in
   `Fig4g_hit_epitope/make_columns.py`. Different structures, not a
   contradiction - but never put the two on one axis.

Join trap: `node_id` in the fingerprint table is a row index, NOT the
coord_hash, and the tables are not in the same order (checked). Go
design_id -> coord_hash (`hashes/*.tsv`) -> `layers.parquet`; the script
asserts all 526,467 land. Euler login nodes have node-local `/tmp`, so write
anything you mean to `scp` into the work directory instead.

Hairpin topology (5-bp stem, GUAUA loop at nt7-11, unpaired A1) is measured
from 4AL5 and shared via `code/common/csy4.py`.

Caption draft: *Epitope of the design campaign, its computational gate and the
ordered library. Mean share of each design's protein-RNA heavy-atom contacts
(4.5 A) per nucleotide of the Csy4 repeat hairpin, split by moiety; shares sum
to 1 per design (n = 497,007 / 31,760 / 9,757 / 32 coordinate-hash nodes,
campaign being a 499,999-node sample). Bars are means with no uncertainty
shown; the hits column rests on 32 backbones, where the bootstrap 95% interval
is about +/-0.024 per nucleotide and overlaps the library at 14 of 16
nucleotides. Contacts are measured on the design models. The gate
contains no epitope term, yet selection concentrates contacts onto the GUAUA
loop and onto bases at the expense of the sugar, and removes designs binding
the 3' end. The 32 ordered candidates (right-hand bar) fall within the
library's distribution at 14 of 16 nucleotides and are not a further
sharpening step.*
