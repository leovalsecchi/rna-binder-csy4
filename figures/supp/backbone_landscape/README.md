# backbone_landscape

The backbone fold landscape through selection, drawn **on main-text Fig 2b's own
map**. **One bare panel, three maps in a row**, cited from Section 3.

> **2026-11-20, Leo:** *"remove the iptm counterfactual and all the text below in
> grey."* Both done. The ipTM counterfactual (map 4) and the grey foot note are
> gone; see "Map 4 is CUT" and "Where the foot note went" below for what left
> with them and where it landed.

```bash
# once, on the cluster (everything it reads lives only there)
scp code/supp/backbone_landscape/make_inputs.py euler:/cluster/work/wicky/lvalsecchi/
ssh euler '/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python3 \
    /cluster/work/wicky/lvalsecchi/make_inputs.py --fig2b'
# -> _landscape_export/fig2b_layers.parquet, copied into inputs/large/
# (the no-flag run rebuilds the union-partition tables the figure only quotes)

~/.venvs/binder_map/bin/python3 code/supp/backbone_landscape/make_figure.py
```

Emits `output/backbone_landscape/landscape_stages.{pdf,png,svg}`, **178 × 72.8 mm**
(double column), three 57.0 × 60.8 mm maps. Dropping map 4 and the 8 mm note strip
gave the remaining three 57.0 mm each instead of 42.2, which is the one unlooked-for
gain: at 42.2 mm the campaign map's four dense masses went to flat plate, and at
57.0 mm their internal texture is back without the marker having to change. Every
number below is re-derived and
printed on every run, and the script asserts the layer sizes rather than trusting
them.

`output/backbone_landscape/panel_{a,b,c,d}_*.{pdf,png}` are the superseded
four-panel version. **The script no longer produces them**, so they will not
update. They are deliberately still on disk (an Illustrator assembly may link
them) and `code/supp/make_supp_book.py`'s `SUPERSEDED` dict already keeps them out
of the book, so the `*.pdf` glob cannot place a cut panel beside its replacement.
The book picks up `landscape_stages` unchanged — the stem did not move.

## The three maps

All three sit on **main-text Fig 2b's embedding** — the `results_bbmap` Foldseek
fold layout, 551,516 backbones, of which 499,064 are drawn — so a position means
the same backbone here and in Fig 2b, and a map differs from its neighbour only by
which backbones are **still there**. Nothing is greyed: each stage *removes* what
it filtered out, which is what Leo asked for. One marker size across all three, and
it is **Fig 2b's own**; the thinning is the message and must not be flattered by
growing the dots on the maps that have fewer of them.

| | map | designs | drawn backbones | outline area |
|---|---|--:|---|--:|
| 1 | The campaign | 6,992,515 | 499,064 of 551,516 drawn | 96.1 |
| 2 | Five applied gates | 55,706 | 26,312/499,064 (5.3%) | 53.2 |
| 3 | Ordered panel | 12,000 | 8,312/499,064 (1.7%) | 52.7 |
| — | ~~ipTM counterfactual~~ **cut 2026-11-20** | 12,000 | 6,627/499,064 (1.3%) | **36.7, measured not drawn** |

Before the 52,452/551,516 (9.5%) isolates are dropped those layers are 551,516 /
26,658 / 8,552 / 6,698.

All four areas are still **measured on every run**, on the same grid: `outline()`
takes `ax=None` and returns the area without drawing the contour, so the cut arm's
36.7 and the 69.0% it gives are re-derived by `report()` rather than frozen into
this table. Nothing on the canvas draws them.

### 🚨 The maps count BACKBONES; the question was asked in DESIGNS

2026-11-19, Leo, looking at the four maps when they were captioned with backbone
counts alone:

> "I'm not sure why there's 3 maps after the first one. The first subsample should
> be the 55,000 we filter to. And then the second subsample should be the one with
> 12,000 designs."

**The figure already did exactly that.** Map 2 *is* the 55,706 gated designs; map 3
*is* the 12,000-design library. They collapse onto 26,312 and 8,312 **drawn
backbones**, because LigandMPNN draws many sequences per backbone, and only the
labels never said which space the number was in. **Every map now states both**: the
design count it represents and the backbone count actually drawn. The design counts
come from `funnel_real_gate.csv` (campaign total and the gate's survivors) and from
`backbone_landscape_coverage.csv` (the two 12,000 arms); `designs()` asserts the two
files agree on the gate, and nothing is typed into the figure script.

### Map 1 is the whole campaign, and the only omission is Fig 2b's isolates

🚨 **This question cost a round of confusion once, and the text that answered it
was in the grey note that is now cut.** Where it went: **map 1's own label**, which
reads `The campaign / 6,992,515 designs / 499,064 of 551,516 backbones drawn` — the
design count is the campaign total and the backbone fraction names the only
omission, both on the map that provokes the question. The isolate drop and
"Nothing else is filtered out of it" are in the **caption draft** below, in full.
Neither is typed in: the label is built in `figure()` from `designs()` and the
drawn table.

> "these 499,000, these are the full campaign, right? There's no other backbones
> that weren't included?"

Yes. **6,992,515 designs → 551,516 backbones → 52,452 (9.5%) unclustered isolates
dropped → 499,064 drawn.** Nothing is filtered, sampled or thinned out of map 1; the
isolate drop is **inherited from main-text Fig 2b**, not a choice made here, and it
is there because a degree-0 node has no mutual k-NN edge surviving the TM 0.5 floor,
so UMAP's random init throws it onto the perimeter where its position means nothing.
That sentence is now in the panel's own foot note, in full, so the question cannot be
asked of the figure again.

**The design-route key ships with the panel**, in a strip under the row, for the
same reason Fig2b_map_composite's does (Leo, 2026-11-16): without it a standalone
export shows four unexplained colours, and the colour is the one thing on these
maps a reader cannot decode from the maps themselves. Labels are Fig 2b's own
short forms, so the two panels read as one system; no counts in the key, because
they differ per map and each map already carries its own. Entry widths are
measured with `style.text_w`, not estimated per character, and the script asserts
the key is narrower than the panel.

**Colour is Fig 2b's and means design-route identity only** (COLORS.md rule 2),
including Fig 2b's BoltzGen merge — BoltzIF is an inverse folding step on BoltzGen
backbones, so on this partition the two routes are the same backbone, and the
merged group keeps the BG+LiMPNN hex. No new hex is invented. The one non-route
mark is the outline, in `MAIN`, a single series (rule 1). No hatch: withdrawn
repo-wide on 2026-11-19 (rule 7); this figure never carried one.

## 🚨 Map 4 is CUT — and this is what left with it

2026-11-20, Leo: *"remove the iptm counterfactual and all the text below in grey."*
Done. The row is three maps: the campaign, the five applied gates, the ordered
panel.

It was the right cut on its own terms — Leo never asked for a fourth map, nothing
passed through it, and until 2026-11-19 the panel gave no sign of that, so a reader
walking the row left to right would have taken it for a fourth filtering step. The
grey label it then got ("ipTM counterfactual (not a stage)") was a patch on a
mislabelled object, not a reason to keep it.

### ⚠️ The one result that leaves the supplement with it

> the ipTM arm contracts to **69.0%** of the gated region, while the ordered panel
> does not contract at all — **52.7 against 53.2 units², 99.2%**.

**That claim now appears in NO figure and NO caption in this repo.** Checked
2026-11-20 across `captions/`, `supp_text/`, every `code/**/README.md` and every
figure script: the only hits for 69.0 / 52.7 / 53.2 / 36.7 / "contracts to" are
this figure's own README and script, and `code/RUNLOG.md`. The nearest live
relative is **main-text Fig 3d**, which compares the same two arms — but it
compares them by **backbone and fold-cluster counts against a design-matched random
null**, not by the territory they occupy, so it carries no area claim and the
contraction is not implied by anything it draws.

**What is NOT lost**, and it is the larger half of the comparison: the k = 30
fold-coverage argument — the ordered panel reaches 917/917 (100.0%) fold clusters
against the ipTM arm's 326/917 (35.6%), with 591 the arm never sees — is a caption
number on a different partition that the maps never drew anyway. It is unaffected
and still in the caption draft below.

**What IS lost is the visual.** The contraction was *visible* in map 4 in a way no
number is: the arm lost almost every peripheral island the other three keep. A
sentence cannot do that.

✅ **Settled 2026-11-20: the sentence is IN the caption** — see "The contraction
sentence" under the caption draft, which records that this was the **coordinating
agent's call applied on Leo's behalf, not Leo's**, and how to reverse it. It is
marked **`stated, not drawn`** (the `bif_effect` / `repredictor_comparison` idiom),
and it reports the areas without implying the reader can see the region. It was
flagged rather than inserted silently in the first place because a caption sentence
describing a region the reader cannot see is a different kind of claim from one
that labels what is drawn — that distinction survives into how the sentence is
worded.

**The measurement itself survives the cut.** `outline(ax=None, ...)` computes the
region without drawing it, so `figure()` still measures the ipTM arm on the same
grid as the three drawn maps and `report()` prints `36.7 UMAP units²` and the
69.0% / 99.2% on every run. The number cannot drift into a stale README value —
the same treatment `metric_reliability` gives the caption clause of *its* cut
panel.

## Where the foot note went

The three-line grey note under the row is cut with the map, same instruction. It
carried three things; here is each one's new address.

| the note said | where it is now |
|---|---|
| **the campaign-scope chain** — 6,992,515 designs → 551,516 backbones → 52,452 (9.5%) unclustered isolates omitted → 499,064 drawn, *"Nothing else is filtered out of it"* | **map 1's own label** states `6,992,515 designs` and `499,064 of 551,516 backbones drawn`; the isolate drop and the "nothing else" sentence are in the **caption** (draft below) |
| the outline holds 90% of that stage | **caption** (it already said so) |
| the fourth map is a counterfactual, not a stage | **gone with the map** |

🚨 The first row is the one that matters: *"these 499,000, these are the full
campaign, right? There's no other backbones that weren't included?"* was asked
once and cost a round. It is now answered **on the panel**, by the map that
provokes it, and again in the caption — not only in a note strip that STYLE_NOTES
§0 was against in the first place (explanatory notes on the canvas belong to the
caption).

## Why map 1 exists, and why the main text's low-N idiom was not copied

Worth recording, because the two instructions pull against each other. Leo asked for
filtered-out dots to be **removed**, not greyed. With no grey field under maps 2–4
there is nothing standing in for the pool they were drawn from, which is exactly
why the campaign map has to be map 1 — it *is* the grey field, unrolled into its own
panel. Cutting map 4 does not touch that: map 1 stands in for the pool under maps 2
and 3 exactly as it did under maps 2, 3 and 4.

The main text's own low-N treatment is `Fig2b_map_composite`'s eight callouts: a big
ringed mark, halo plus route colour, over the **full coloured cloud**. That is a
**per-point** idiom for **eight** points. It does not transpose to 26,312, and
adopting it would mean putting the filtered-out backbones back under every map —
the greying Leo rejected. **So the low-N fix taken here is Fig 2b's marker, not Fig
2b's callout**, and nothing is greyed or thinned. If Leo later prefers a grey pool
under each map, map 1 becomes redundant and the panel can go to two.

## 2026-11-19: the embedding swap, and the previous README's error

> "I want to use exactly the campaing umap, like it is in main fig2b, then remove
> the dots that go filtered out and so on. I want to keep the clustering we have
> in the main conserved. Also with the colors."

Done, on all four counts (then four maps; three since 2026-11-20). 🚨 **The previous version of this figure did none of it,
and the reasoning it gave was wrong on its decisive point. That argument is
withdrawn.** It claimed the two tables were "NOT joinable by id" and that Fig 2b's
cascade was "the wrong gate — the one that produced the retired 5,264". Against the
files:

- **Three of the four maps need no join at all.**
  `inputs/large/fold_map_scatter.parquet` already carries a `stage` column. The
  previous README quoted it (524,858 / 18,106 / 8,552) and did not notice that it
  *is* the cascade. ⚠️ `stage` is the **furthest** stage a backbone reached, so the
  gated set is `stage >= 1` = **26,658**, not `stage == 1` = 18,106.
- **The fourth map joins too**, by `design_id`, through
  `build_bbmap/meta/design_to_node.parquet` — which carries every one of the
  6,992,515 designs. The "not joinable" line in `inputs/MANIFEST.md` is true of the
  *exported* parquet, which drops `node_id` and `coord_hash`; it was never true of
  the tables upstream on Euler, and `counterfactual3.csv` had already been placed
  on this same map by this same bridge.
- **The embedding and the clustering are reused, not regenerated.** Nothing is
  re-embedded and nothing is re-clustered: the same 551,516 nodes with the same
  Leiden `cluster_r10`. This figure's own k = 30 r8.0 / 917-fold clustering no
  longer draws anything — it is only quoted, for the coverage numbers.

`make_figure.py` **proves** the reuse rather than asserting it: `fig2b_layers.parquet`
is Fig 2b's own table plus one flag, and `load()` checks that its `x`, `y`,
`stage`, `cluster_r10`, `stack` and `deg0` reproduce `fold_map_scatter.parquet` row
for row. If the embedding ever drifts, the figure stops building.

**Nothing is left that is not Fig 2b's**, as of the 2026-11-19 second pass. The
marker used to be the one exception — a shrunken local `DOT = 0.22`, invented here
on the argument that four 42.2 mm maps cannot carry a 178 mm map's dots. It cost
the thin maps their legibility: at 26,312, 8,312 and 6,627 marks they read as dust.
Leo: *"for the umap just take the umap from the main figure and use its code and
display setting for low number binders."* So `fig2b_marks()` now **reads `s` and
`alpha` out of `Fig2b_map_callouts/make_figure.py`'s own `ax.scatter` call** —
`s = 0.5`, `alpha = 0.55` — rather than keeping a copy that can drift. That file is
a script, not a module, so importing it would rebuild the main-text map; reading its
source is the cheap way to get the same guarantee `load()` gives for the table.
⚠️ **The cost, and it is real:** at `s = 0.5` the campaign map's four dense masses
went to flat saturated plate and lost their internal texture. The silhouette, the
filaments and every peripheral island still read, and main-text Fig 2b saturates the
same way at 178 mm, so the character matched — but if that texture were wanted back,
the only honest fix would be a smaller marker **on every map**, not on map 1 alone.
✅ **Largely moot since 2026-11-20:** with map 4 and the note strip gone, each map is
57.0 mm instead of 42.2, and at that size the dense masses hold their texture at the
same `s = 0.5`. Nothing was retuned to get it.
Everything else is carried over: the
window (full extent of the non-isolate cloud, 3% pad, aspect equal — *not* forced
square, because Fig 2b's is not), the bare frame with no spines, ticks or axis
labels, the isolate drop, the randomised draw order on seed 0, the colours, and the
outline's
`GB` = 0.10, `SIGMA` = 3.0, `HDR_MASS` = 0.90 at Fig 2b's own 0.5 pt / alpha 0.6,
beneath the marks.

Not thinned. Fig 2b's composite bins at 0.12 and caps 6 per bin because it needed
41,072 true-vector marks Illustrator would open and 8-point callouts that could be
seen over them. Neither applies here: these maps are rasterised at 600 dpi and
carry no callouts, and thinning the campaign map but not the thin layers would make
26,312 backbones look as dense as the 499,064.

## 🚨 Settled: Fig 2b's cascade IS the applied gate

This mattered — Fig 2b is a main-text figure, and the previous README called its
gate wrong. It is not. Re-derived on Euler, from the design ids rather than from
the stored flags:

- `build_bbmap/meta/design_to_node.parquet` carries per-design `in_passing` =
  **55,706** and `in_12k` = **12,000** — the applied gate's design count and the
  delivered library's, exactly.
- Recomputing the node set from `passing_v6_final.csv`'s 55,706 design ids gives
  26,664 nodes, **set-identical** to the stored `in_passing`; from the 12k library's
  12,000 ids, 8,557 nodes, **set-identical** to `in_12k`. 26,658 and 8,552 of those
  survive into the drawn table (6 and 5 nodes carry no embedding or cluster row).
- The rival `all_methods.py` gate (66,253) appears nowhere on this path. **5,264
  came from `counterfactual3.csv`** — the `all_methods.py` pool *placed on* this
  embedding. Same map, different pool. The previous README read that as the map's
  own cascade.

**So the 26,658 against 31,763 difference is a partition effect, not a gate
difference** — and the cause is now named. The bbmap backbone key is
`<gen>_<uuid8>_<designN|recycleN>` with BG+BIF and BG+LiMPNN merged
(`build_bbmap/meta/build.json`), **not a coordinate hash**, so it is coarser than
the k30 partition's. Per design route, same applied gate:

| design route | k30 gated | bbmap gated | k30 ordered | bbmap ordered |
|---|--:|--:|--:|--:|
| RFD3 + LigandMPNN | 14,455 | **14,455** | 2,782 | **2,782** |
| BoltzGen (+ LiMPNN, + BoltzIF) | 14,732 | 9,982 | 4,778 | 3,926 |
| BoltzDesign1 + LigandMPNN | 1,774 | **1,774** | 1,397 | **1,397** |
| HalluDesign + LigandMPNN | 802 | 447 | 802 | 447 |
| total | 31,763 | 26,658 | 9,759 | 8,552 |

**Two design routes agree to the unit, on both layers.** A different gate cannot do
that. The whole deficit is BoltzGen's BoltzIF merge (−4,750 gated, −852 ordered;
BG+BIF survives as 107/551,516 (0.02%) nodes on this partition) plus HalluDesign's
name collapse (−355 on both). −4,750 − 355 = −5,105 = 31,763 − 26,658, and
−852 − 355 = −1,207 = 9,759 − 8,552. The arithmetic closes exactly.

## 🚨 The ipTM counterfactual CAN be placed here, and is

**12,000/12,000 (100.0%)** of the `method == "top ipTM"` designs in
`hier_k30_max_selections.parquet` place on a bbmap node by `design_id`, giving
6,700 backbones, 6,698 in the drawn table and **6,627 after the isolate drop**.
Nothing is lost to coverage.

The bridge is **proved by a control, not asserted**: the ordered arm of the same
file, placed the same way, lands on **8,552/8,552 (100.0%)** stage-2 backbones and
nowhere else. `make_inputs.py` asserts that equality, so a wrong bridge cannot ship.

On this partition the arm holds 6,698 backbones against the panel's 8,552:
4,151/6,698 (62.0%) are shared, 4,401 are the panel's alone and 2,547 are the arm's
alone. ⚠️ Those are **bbmap** numbers and are not the union-partition 4,615 / 5,144
/ 3,198; the same coarser key is why (7,813 k30 backbones collapse to 6,700 here).
The caption quotes the union/k30 figures, which is where the coverage argument was
measured.

## ⭐ The outline-area finding, re-measured — and half of it does not carry over

The 90% density outline is recomputed per map, so it contracts across the row, and
its enclosed area is a result rather than decoration. On the union partition the
areas were gated 39.5, ordered 47.7, ipTM 27.8 UMAP units², and the striking half
of that was that **the ordered panel's region was *wider* than the gated pool's it
was drawn from**. On Fig 2b's embedding:

| | area | against the gated pool |
|---|--:|--:|
| the campaign | 96.1 | — |
| after the five applied gates | 53.2 | 100% |
| the ordered panel | 52.7 | **99.2%** |
| ipTM-ranked top 12,000 | 36.7 | **69.0%** |

⚠️ 99.2%, not the 99.1% this table carried before: the script now prints the ratio
from the **unrounded** areas, where dividing the rounded 52.7 by the rounded 53.2
gives 99.06%. Same measurement, one more digit of it. The 69.0% is unchanged.

🚨 **The "wider" half does NOT survive the change of embedding.** Here the ordered
panel's region is the *same size* as the gated pool's — 52.7 against 53.2, a 0.8%
difference — not wider. Do not repeat the 47.7-against-39.5 claim for this figure;
those numbers were measured on a different layout and are not comparable to these.
What does survive, and strongly, is the other half: the **ipTM arm contracts to
69.0% of the gated region** while the ordered panel does not contract at all. That
is the comparison the figure is for, and it is visible in the panel — the ipTM map
loses almost every peripheral island the other three keep.

The areas are not directly comparable between the two embeddings in any case: they
are in each layout's own UMAP units, and the campaign cloud here spans
24.6 × 26.2 units.

## Sources

| file | mtime | what it gives |
|---|---|---|
| `inputs/large/fig2b_layers.parquet` | 2026-11-19 20:19 | 551,516 × (x, y, stage, cluster_r10, stack, deg0, in_counterfactual) — **the drawn table** |
| `inputs/large/fold_map_scatter.parquet` | 2026-11-05 14:46 | main-text Fig 2b's own table; the figure asserts the two agree |
| `inputs/large/backbone_landscape_layers.parquet` | 2026-11-19 11:08 | 526,467 × (node_id, x, y, 6 flags) — union partition, **no longer drawn**, still asserted |
| `inputs/backbone_landscape_coverage.csv` | 2026-11-19 11:10 | per-arm designs / backbones / folds / epitopes / pairs / mean ipTM |
| `inputs/backbone_landscape_quantities.json` | 2026-11-19 11:11 | 41 set quantities over pool, ordered and counterfactual |

`fig2b_layers.parquet` is built by `make_inputs.py --fig2b` on Euler from
`foldseek_map/results_bbmap{,_pass,_lib}/design/` (the embedding, the clusters, the
fold graph and the two subset trees), `build_bbmap/meta/{nodes,design_to_node}.parquet`
(the bridge and the gate flags), `passing_v6_final.csv`, the 12k library csv and
`library_selection/hier_k30_max_selections.parquet`. Construction of the six shared
columns is `fig05_fold_clustering/export_fig.py`'s, unchanged — that is what
produced `fold_map_scatter.parquet` in the first place.

### The four partitions, which do not share ids

`union` (526,467 nodes), `k30` (31,763), `lib12k` (9,759) and **`bbmap`
(551,516, the one drawn here)** are different Foldseek trees. The first three key
on `coord_hash`; **bbmap keys on a name**, which is why it is coarser. All four are
reachable from `design_id`, which is the bridge this figure uses and the one the
previous README said did not exist.

## Verified against the files

- **The applied gate**, from the header of `foldseek_map/passing_v6_final.csv`
  (2026-09-01 17:12), which states `n=55706 store_n=6992515`: `fmax(iptm_af3,
  iptm_esm) ≥ 0.6`; `rmsd_af3 ≤ 3 OR rmsd_esm ≤ 3`; `dockq_af3 ≥ 0.7 OR dockq_esm
  ≥ 0.7`; `seq_run ≤ 5`; `footnb_010_af3 ≤ 5 AND footnb_010_esm ≤ 5`. **The foot
  clash term is AND; the other three are OR.** Every other atom in that JSON is
  `"on": false`. `v6k30_pool.parquet` is exactly 55,706 rows.
- **Layer totals on the drawn partition**, re-read 2026-11-19: 551,516 total,
  gated 26,658, ordered 8,552, ipTM arm 6,698, isolates 52,452/551,516 (9.5%).
  Furthest-stage buckets 524,858 / 18,106 / 8,552.
- **Layer totals on the union partition**, unchanged and still asserted:
  sample500k 499,999, passing 31,763, in12k 9,759, counterfactual 7,813. Overlaps:
  sample ∩ passing 5,295; sample ∩ 12k 1,553; passing ∩ 12k 9,759 (the panel is a
  strict subset of the gated set).
- **The coverage map is k = 30, r8.0** — `results_sw_k30_f50_mut`, **917 fold and
  437 epitope clusters**, 661 of the fold clusters singletons. 🚨 **The Methods text
  saying k = 15 / 2,166 clusters is wrong**, and the file it describes exists:
  `results_v6` ( = `results_sw_k15_f50_mut`) gives `cluster_r80` = 2,166 on the same
  31,763 backbones.
- **The diversity argument**, re-read 2026-11-19. Hierarchical fold → epitope
  selection reaches 917/917 (100.0%) folds, 341/437 (78.0%) epitopes and
  7,784/8,819 (88.3%) (fold, epitope) pairs at mean ipTM 0.7647; the pure top-ipTM
  counterfactual reaches 326/917 (35.6%), 196/437 (44.9%) and 3,980/8,819 (45.1%)
  at 0.8182. **The cost is 0.054 mean ipTM for 917 folds instead of 326.**
- ⚠️ **The coverage numbers are k30, the drawn counts are bbmap**, and the figure
  keeps them apart on purpose: a k30 fold count beside a bbmap backbone count on
  the same map would mix partitions. The maps carry only what you can see on them;
  the coverage goes in the caption, which names the partition once.

## 🚨 The suspect numbers, re-derived

| number | verdict |
|---|---|
| **591** | **Reproduces.** Fold clusters the ordered panel reaches and the counterfactual does not = 917 − 326. The counterfactual's folds are a **strict subset** of the panel's (counterfactual-only = 0), which is why the subtraction is also the set difference. |
| **5,144** | **Reproduces**, on the union partition. Backbones in the ordered panel and not in the counterfactual (both = 4,615, counterfactual-only = 3,198). On the drawn bbmap partition the same three numbers are 4,401 / 4,151 / 2,547. |
| **1,946** | **Reproduces.** Net backbone difference, 9,759 − 7,813 (union). |
| **5,264** | **Retired, and its origin is confirmed.** The library-only backbone count inside `inputs/counterfactual3.csv` — the `all_methods.py` pool, i.e. the wrong gate, *placed on this embedding*. It is not this map's cascade, which is the applied gate. |
| **528** | **Does not exist.** No quantity over pool / ordered / counterfactual backbones, folds, epitopes or (fold, epitope) pairs equals it. "Both" is **326** folds, not 528. The script asserts this. |
| **1,515** | **Does not exist**, likewise asserted. The nearest real number in the neighbourhood is 1,553 — backbones in both the campaign sample and the panel — but nothing establishes that 1,515 came from it, and it is not reported as if it did. |
| **39.5 / 47.7** | **Withdrawn for this figure.** Union-partition outline areas; on Fig 2b's embedding the gated and ordered regions are 53.2 and 52.7, so the "ordered is wider than the pool" reading does not hold here. |

## Caption draft

⚠️ **Rewritten 2026-11-20 for three maps.** Three changes beyond dropping the
fourth: the **campaign-scope chain** moved in here from the cut grey note (second
sentence, in full, including *"nothing else is filtered out of the first map"*);
the ipTM arm is named as a selection that was **not** made and is **not drawn**;
and the **outline-area contraction** is now carried as prose — see "The
contraction sentence" below for whose call that was.

**Figure Sx. Selecting by fold and epitope keeps the whole backbone landscape;
selecting by confidence collapses it onto one lobe.** Main-text Fig. 2b's
Foldseek/UMAP fold embedding of the campaign's design backbones, drawn three times
with the same clustering and the same design-route colours; each map keeps only the
backbones still present at that stage, so the map thins from left to right. The
campaign's 6,992,515 designs collapse onto 551,516 backbones, of which 499,064 are
drawn: the 52,452/551,516 (9.5%) degree-0 isolates are omitted from all three maps,
their position in the layout carrying no meaning, and **nothing else is filtered
out of the first map**. From the left: the whole campaign; the 26,312/499,064
(5.3%) backbones of the 55,706 designs surviving the five applied gates; and the
8,312/499,064 (1.7%) backbones of the 12,000-design ordered panel. The outline is
the smallest smoothed-density region holding 90% of that stage; it encloses 96.1
UMAP units² for the campaign, 53.2 after the gates and 52.7 for the ordered panel,
so the ordered panel holds 99.2% of the gated pool's footprint on 1.7% of its
backbones. **Not drawn here:** the same twelve thousand designs ranked by ipTM
alone instead occupy 36.7 units², 69.0% of the gated region. Measured on the
k = 30 Leiden r8.0 fold map of the same gated pool, the ordered panel reaches
917/917 (100.0%) fold clusters, 341/437 (78.0%) epitope clusters and 7,784/8,819
(88.3%) (fold, epitope) pairs at mean ipTM 0.7647, against 326/917 (35.6%),
196/437 (44.9%) and 3,980/8,819 (45.1%) at 0.8182 for that ipTM-ranked selection:
a cost of 0.054 mean ipTM for 917 fold clusters instead of 326, with 591 fold
clusters the ipTM arm never sees. n = 12,000 designs in the third map, from a
gated pool of 55,706.

### The contraction sentence — added 2026-11-20, and whose call it was

🚨 **This was NOT Leo's call. It was the coordinating agent's, applied on his
behalf**, and it is recorded that way so it can be reversed cleanly: delete the
sentence beginning *"Not drawn here:"* and the `96.1 ... 99.2%` clause before it,
and nothing else in the figure changes. Leo's instruction on 2026-11-20 was only
*"remove the iptm counterfactual and all the text below in grey"* — he left the
question of the sentence open across several turns after being offered a default.

The reasoning on the record: the claim appears in **no other figure and no other
caption** in this repo (checked across `captions/`, `supp_text/`, every
`code/**/README.md` and every figure script), and the measurement still runs on
every build, so the only live question was whether a reader ever meets it.
Dropping the map was about page economy, not about retiring the finding.

**It is marked `stated, not drawn`**, the same way `bif_effect` marks its χ² and
`repredictor_comparison` marks its 99.2% / 59.8% ipTM-step rates: the caption says
in words that the ipTM arm is not on the panel, so no reader can take the sentence
for a description of something in front of them.

⚠️ **What the sentence deliberately does NOT claim.** The contraction was *visible*
in map 4 — the arm lost almost every peripheral island the other three keep — and a
number cannot carry that. The caption therefore reports the areas and stops; it
does not say the arm "collapses onto one lobe" of anything the reader can see, and
it must not be edited into implying the evidence is on the page. (The figure's
title sentence makes that claim about *selection*, which the fold-coverage numbers
support independently of any map.)

**It cannot drift.** Every number in it — 96.1, 53.2, 52.7, 36.7, 99.2%, 69.0% —
is derived on every run: the three drawn outlines by `figure()`, the ipTM arm's by
`outline(ax=None, ...)` on the same grid, and all of them printed by `report()`.
None is typed into this README as a source.

## Open

- The caption's figure number is a placeholder; supplementary numbering is not
  settled, so it appears in no file or directory name here.
- ✅ ~~`NEED:` the 69.0% contraction sentence~~ — **settled 2026-11-20, it is in
  the caption**, marked `stated, not drawn`. Decided by the coordinating agent on
  Leo's behalf, not by Leo; reversible in one edit (see "The contraction
  sentence"). Worth a glance from him, since it is the only caption claim here
  about a region no map shows.
- **A per-route outline** — Fig 2b's four-colour background line — is now available
  in principle, since `fig2b_layers.parquet` carries the design route. Not done:
  four outlines × three maps is twelve lines on a 178 mm row, and the maps are
  about stages, not routes. A single series takes `MAIN` (COLORS.md rule 1).
- The maps keep Fig 2b's full-extent window, so each is about 55% occupied
  horizontally, and the row is 72.8 mm tall because of it. That is deliberate — a tighter crop would stop the panel
  registering with the main-text map — but it is why the row looks airy.
- `inputs/large/backbone_landscape_layers.parquet` and the union partition are no
  longer drawn. They are still read and asserted, because the coverage numbers and
  the set quantities in the caption come from that side. Do not mix the counts.
