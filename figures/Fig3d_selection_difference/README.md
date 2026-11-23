# Fig3d_selection_difference

The difference between the two selection rules, drawn directly, instead of two
absolute maps the reader has to subtract by eye.

**State on 2026-11-12: sections b and c, both in the right column.** The
canvas is still 178 x 100 mm and the left column is still free for whatever
goes there.

| Section | What |
|---|---|
| b | what each rule reduces: unique backbones and fold clusters kept |
| c | what each rule does to the stack mix, as grouped columns |

## The two arms
Both are selections of **12,000 designs from the same candidate pool**, made
by different rules. They are pre-selections in the same sense: neither is an
experimental result, both are a choice of what to order.

| in the figure | what it is |
|---|---|
| **ipTM-ranked 12,000** | the counterfactual: take the 12,000 highest ipTM designs and nothing else |
| **diversity-optimised 12,000 (ordered)** | the rule that was actually used, maximising fold-cluster coverage; this is the library that was synthesised |

Named this way, changed from "top 12,000 by ipTM" and "ordered library" on
2026-11-12, because "ordered library" named the *outcome* and said nothing
about the rule, which made the panel look like a comparison of a method
against a product rather than of two rules. "(ordered)" keeps the fact that
this arm is the real library and not a hypothetical. Both names are set once
in `ARMS` in the script and used by both panels.

Warm is the diversity-optimised arm, cool is the ipTM-ranked one.

## What it shows

### c, what each rule reduces

Rebuilt onto the k30 partition on 2026-11-16. See "Rebuilt onto k30" below for
why, and for the numbers this replaced. Pool: 55,706 gated designs on 31,763
unique backbones in 917 fold clusters (Foldseek TM-align, mutual k = 30 at
TM >= 0.50, Leiden r8.0). Both arms are exactly 12,000 designs from that pool.

The 917 clusters are 256 holding two or more backbones and 661 holding one, and
the panel draws them as separate rows because the two rules separate on exactly
that line.

| | backbones / 31,763 | | populated folds / 256 | | one-member folds / 661 | |
|---|--:|--:|--:|--:|--:|--:|
| ipTM-ranked 12,000 | 7,813 | 25% | 206 | 80% | 120 | 18% |
| diversity-optimised 12,000 | 9,759 | 31% | 256 | 100% | 661 | 100% |
| random draw, 300 x 12,000 designs | 9,994 +/- 38 | 31% | 210 +/- 5 | 82% | 194 +/- 11 | 29% |

**What the panel argues.** Not that the library covers more, which is true by
construction. That ipTM ranking is AT CHANCE on populated folds (206 against
210 +/- 5, 0.8 SD) and far below it on isolated backbones (120 against
194 +/- 11, 6.6 SD). Score ranking does not narrow the populated fold
repertoire; it rejects structural orphans.

**The null is design-matched**, 300 draws of 12,000 designs from the pool's
55,706, seed 0, computed in `figures/k30_fig3c_inputs.py` on Euler and drawn as
a third bar with a +/- 1 SD whisker. It is design-matched and not
backbone-matched so that the backbone row says something: a draw of N backbones
holds N backbones.

#### Three things to be careful about
1. **Several clusterings in this repo use the same column names.** SUPERSEDED
   2026-11-16: this warning described the map partition the panel used until
   then (`fold_map_scatter.parquet`, `counterfactual3.csv`, 139 folds). The
   panel is now on k30 (`counterfactual_k30.csv`, 917 folds), which is the
   partition the selection itself ran on, so the mixing hazard is gone for this
   panel. It is not gone for the repo: `inputs/fold_coverage_pool.csv` is a
   Leiden r1.0 partition of the passing set at 732 folds, and the map partition
   is still used by the Fig 2b panels. A coverage percentage is meaningless
   without naming its partition. Never compare across them.
2. **The library's side is partly circular.** The 12k selection was built to
   maximise cluster coverage. It ran on the k30 fold *and* epitope map while
   this measures fold-only coverage in a different partition, so it is not
   pure tautology, but do not caption it as a discovery. The ipTM number is
   the one that is not circular.
3. **346 isolate backbones are excluded** (`deg0`, degree 0 in the thresholded
   clustering graph), so the pool is 26,312 here and 26,658 where isolates are
   counted, as in the archived `fig06_fold_diversity`.

### Which pool the selections actually came from - RESOLVED 2026-11-16

Resolved by the k30 rebuild, not by tracing the old files. The panel no longer
reads `counterfactual3.csv`, so the gate mismatch below no longer touches it:
both arms are now taken from `hier_k30_max_selections.parquet`, which selects
12,000 from exactly the 55,706 that `funnel_real_gate_v2.csv` counts. The grey
bar and the two coloured bars are therefore the same gate by construction, and
no arm can hold a fold outside the pool. The account below is kept because it
documents a real trap in the older files, which the Fig 2b panels still use.

`counterfactual3.csv` was built from the **`all_methods.py` candidate pool**
over `cross_campaign/filter_store_v2v9` (MANIFEST line for that file).
`funnel_real_gate_v2.csv`, the grey bar in c, is the **`passing_v6_final`
gate** over `filter_store_v2v9_hdtop`, which reproduces 55,706. Those are a
different gate definition over a different store, and the campaign notes
record that the `all_methods.py` gate yields about 66k designs rather than
55,706.

Direct evidence in the data: **2 of the 78 folds the ipTM arm holds, and 8 of
cluster_flat's 94, do not exist in the passing partition at all.** A selection
drawn from the 55,706 could not do that.

So, until this is traced back on Euler:
- c's grey bar is the 55,706 gate and is a legitimate reference point, but it
  should NOT be described as the set the two selections were drawn from. An
  earlier draft of this README said exactly that; it is corrected here.
- b's pool (26,312 backbones, 139 folds) is the passing partition from
  `fold_map_scatter.parquet`, which is the pool the *map* defines. The arms
  sit almost entirely inside it, but not perfectly, for the same reason.

### The partition question, settled 2026-11-16

`139` is not a cut on a finer partition. It is the campaign `results_bbmap`
partition restricted to gated non-isolate backbones: 26,658 gated minus 346
degree-zero isolates gives 26,312, and 485 clusters minus those same 346 gives
139, because every degree-0 node is its own singleton community. There is no
minimum cluster size anywhere in the script, and 44 of the 139 clusters hold
exactly one backbone.

Methods 2.9 describes a DIFFERENT partition: `results_v6`, which is the sweep
variant `results_sw_k15_f50_mut` (mutual k-NN k=15, TM floor 0.50), 31,763
backbones in 2,166 clusters at resolution 8.0. It holds exactly the same 55,706
passing designs. The whole difference is the dedup key:

| | backbones |
|---|---|
| v6 structural dedup | 31,763 |
| less coordinate splitting (1,817 backbones split 2-10 ways) | -3,311 |
| bbmap key `<gen>_<uuid8>_<designN>`, BG unmerged | 28,452 |
| less the BG+BIF / BG+LiMPNN merge | -1,788 |
| bbmap key, merged (`build_bbmap/meta/build.json`) | 26,664 |
| less 6 not in the exported parquet | -6 |
| this panel's gated set | 26,658 |

**Resolution is not the knob.** Across r0.5 to r8.0 every partition moves by
under 10%. `k` is what matters: at TM 0.50 mutual, k=8 gives 4,538 clusters,
k=15 gives 2,166, k=30 gives 917, k=50 gives 526. Methods 2.9 cites the
resolution and omits k.

**Coverage recomputed in the finer partitions** (isolates dropped as here, null
size-matched to each arm, `output/Fig3d_selection_difference/fig3c_sizematched.csv`):

| partition | folds | library | its null | lib delta | ipTM | its null | ipTM delta |
|---|---|---|---|---|---|---|---|
| bbmap (this panel) | 139 | 82.0% | 74.8% | +7.2 | 56.1% | 71.2% | -15.1 |
| v6 r1.0 (k15) | 423 | 80.9% | 67.1% | +13.8 | 46.6% | 53.7% | -7.2 |
| v6 r8.0 (Methods 2.9) | 600 | 86.2% | 76.5% | +9.7 | 61.5% | 66.8% | -5.3 |
| k30 r1.0 | 147 | 94.6% | 71.8% | +22.8 | 53.7% | 60.2% | -6.5 |
| k30 r8.0 | 256 | 96.9% | 83.6% | +13.2 | 72.3% | 76.6% | -4.3 |

Both of the panel's claims survive in every partition: the library beats its
size-matched null everywhere (+9.7 to +22.8 points) and ipTM is below its own
null everywhere (-4.3 to -7.2). The DIRECTION is robust. The MAGNITUDES are
not: this panel's coarse partition understates the library's margin by half and
overstates ipTM's shortfall by two to three fold.

**Why the panel nevertheless stays on bbmap.** The v6 family is built on the
passing set only, so it holds 100% of the library arm but just 71% of the ipTM
arm (3,479 of 11,895 ipTM designs are not in the 55,706). Scoring a full arm
against a three-quarter arm is worse than scoring both coarsely. `results_bbmap`
is campaign-wide and is the only partition that carries both selections whole,
exactly as the partition warning at the top of `make_figure.py` says.

VERIFY CLOSED 2026-11-16, the k15/k30 question. The library was selected on
**k30**, not k15. Evidence, all primary: `library_selection/k30_library_stats.txt`
reports the shipped 12,000 as "hierarchical fold→epitope, k=30, max-ipTM
ranking" over a pool of 55,706, covering 917/917 fold clusters and 661/661 fold
singletons on 9,759 backbones; `hier_k30_max_selections.parquet` reproduces that
arm exactly, including its per-stack composition (3,132 RFD3 / 3,266 BG+LiMPNN /
2,400 BG+BIF / 2,400 BD1 / 802 HD). So Methods 2.9, which documents k15 and
quotes 2,166 clusters, describes a partition that selected nothing. It must be
rewritten to k=30 at r8.0, 917 clusters. This is a Methods error, not a figure
error.

VERIFY: 2,166 at r8.0 includes 1,566 singleton isolates; only 600 are
multi-member communities. The same holds for k30, 917 of which 661 are
isolates. If Methods 2.9 quotes a fold count, say what fraction of them are
single backbones. Still open only as a wording task in Methods; the figure now
states the split itself.

### Rebuilt onto k30, 2026-11-16

Panel c no longer uses `counterfactual3.csv` or the 139-cluster map partition.
It reads `inputs/counterfactual_k30.csv` and `inputs/stackmix_k30.csv`, built by
`figures/k30_fig3c_inputs.py` on Euler from `hier_k30_max_selections.parquet`
and `v6k30_pool.parquet`. Both arms are now exactly 12,000 designs from the same
55,706, on the partition that selected the library, which removes the 71%
problem described above: that argument was for keeping the map partition when
the only alternative was mapping the OLD ipTM arm into a passing-set clustering.
It does not apply here, because the ipTM arm is re-derived inside the correct
pool rather than imported.

What the panel now shows, per 12,000 designs:

| | backbones / 31,763 | populated folds / 256 | one-member folds / 661 | mean ipTM |
| --- | --- | --- | --- | --- |
| ipTM-ranked | 7,813 (25%) | 206 (80%) | 120 (18%) | 0.818 |
| library (ordered) | 9,759 (31%) | 256 (100%) | 661 (100%) | 0.765 |
| random, 300 draws | 9,994 ± 38 | 210 ± 5 | 194 ± 11 | 0.726 |

**The finding changed.** On the old partition the claim was "ranking by ipTM is
worse than drawing at random". Split by cluster size, that is true only for
one-member clusters (120 against 194 ± 11, 6.6 SD below) and false for
populated ones (206 against 210 ± 5, 0.8 SD, i.e. at chance). The aggregate
326/917 against 404 ± 12 is driven entirely by the isolates. The defensible
claim is narrower and more specific: **score ranking does not collapse the
populated fold repertoire, it rejects structurally isolated backbones.**

The library's 100% on both rows is **by construction** and must never be
reported as a result: the rule is a fold round-robin and 12,000 picks over 917
clusters guarantee full coverage. What is measurable is the price, 0.053 mean
ipTM against the ipTM arm, paid on slightly FEWER backbones than chance would
give (9,759 against 9,994), so the coverage came from where the picks landed.

Superseded files, kept because older drafts quote them:
`output/Fig3d_selection_difference/fig3c_sizematched.csv` and
`fig3c_partition_sweep.csv` (the five-partition sweep that led here),
`k30_threeway.csv` and `k30_fold_split.csv` (all six selection rules, of which
the panel draws two). `make_figure.py.bak_pre_k30` is the 139-cluster version.

### The two attempts this replaced, and why they failed

Both tried to answer fold diversity through secondary structure, and both died
on the same fact: **74% of the library and 84% of ipTM have no beta-strand at
all.**

- **Helix content as a split violin.** Readable, but helix content cannot tell
  a beta-rich backbone from a floppy one, because the non-helical part is
  mostly loop either way (below 60% helix: mean 0.22 strand against 0.37 loop).
  Its left end needed a footnote to say what it meant.
- **Alpha/beta balance, (helix - strand) / (helix + strand).** Makes the left
  end beta by construction, but every strand-free backbone lands on exactly
  +1. Drawn over all designs, the smoother turns that point mass into a wall at
  the right edge, the real distribution flattens to a sliver and both medians
  read +1.00. Drawn over only the strand-bearing 26% / 16%, it is honest but
  silently about a quarter of the data.

The composition question itself, if it is ever wanted, belongs in **class bars**
rather than a density, because a point mass is what a bar is for:

| % of arm | all-alpha | alpha+beta | all-beta |
|---|--:|--:|--:|
| diversity-optimised 12,000 | 73.7 | 25.3 | 1.0 |
| ipTM-ranked 12,000 | 83.8 | 16.0 | 0.2 |

Still printed on every rebuild, not drawn.

### References behind c
- Proteina, *Scaling Flow-based Protein Structure Generative Models*
  (arXiv 2503.00710), appendix F.1: Foldseek `easy-cluster` at TM 0.5 as the
  diversity metric for generated backbone sets. Borrowed: cluster count as the
  diversity statistic.
- *Protein FID: Improved Evaluation of Protein Structure Generative Models*
  (arXiv 2505.08041): argues designability/novelty/diversity miss how well a
  set covers the space, and works on Foldseek cluster representatives.
  Borrowed: coverage of an available set, rather than a bare count.
- *Is Retrieval All You Need?* (arXiv 2608.10598): reports coverage counts at
  topology and superfamily level across generators. Borrowed: reporting
  coverage against a fixed reference partition.
- `archive/unplaced_2026-11-11/code/fig06_fold_diversity` (this repo):
  dots-not-curves, the equal-size random draw as the null, error bar on the
  null only, and the evenness caveat.

Not borrowed: rarefaction/accumulation curves, which imply a continuum of
measurements that does not exist here (one library, one ranking), and ternary
helix/strand/loop plots, which would put 74-84% of each arm on a single edge.

### d, what each rule does to the stack mix

Share of each set that comes from each generator, three columns per stack: the
campaign gate in grey, the ipTM-ranked 12,000 in the cool state colour, the
diversity-optimised 12,000 in the warm one.

| % of the set | RFD3 | BG+LiM | BG+BIF | BD1 | HD |
|---|---|---|---|---|---|
| all 6,992,515 designs, in the note | 21.8 | 21.7 | 21.2 | 20.7 | 14.6 |
| **passed the campaign gate, 55,706** | 50.1 | 26.5 | 16.2 | 5.8 | 1.4 |
| **ipTM-ranked 12,000, n=11,895** | 46.4 | 33.8 | 15.1 | 3.2 | 1.5 |
| **diversity-optimised 12,000, n=11,623** | 26.1 | 27.7 | 20.3 | 19.1 | 6.9 |

**Shares, not counts.** The sets are 55,706 against about 11,600 designs, so
counts would only show which is bigger.

**The grey bar is the gate, not the campaign** (changed 2026-11-12). See the
unresolved-pool note above: the gate is a legitimate reference point but is
not established to be the exact set the two selections were drawn from.
Against it:

- Ranking on ipTM barely moves the mix at all. 50% RFD3 becomes 46%, BD1 goes
  6% to 3%, HD stays at about 1%. It mostly keeps what the gate already chose.
- The diversity-optimised arm is the step that moves it: BD1 6% to 19%, HD 1% to 7%,
  RFD3 50% down to 26%. It is the only one of the three that pushes the mix
  back towards the campaign's own spread.

The campaign before the gate is near flat - the five stacks were run to roughly
equal size - and is not in the panel at all, because against a flat baseline
both selections look like the same big collapse and the difference between the
two rules, which is what this figure is about, disappears. The collapse is the
gate's doing, and the gate is one step upstream of this comparison.

**The panel was stripped on 2026-11-12** to the bars, their numbers, the three
names and the stack labels. No section title, no campaign note under the axis,
no `n=` in the key. **So the caption has to carry what the panel no longer
says:** the three set sizes (55,706 / 11,895 / 11,623) and, if the campaign
baseline is wanted, its flat mix. All of it is in the table above and still
printed to the console on every rebuild.

## What changed
### 2026-11-12
- **The helix-content violin was replaced by this one.** Helix content cannot
  separate a beta-rich backbone from a floppy one - the non-helical part is
  mostly loop either way (below 60% helix: mean 0.22 strand against 0.37 loop,
  and strand beats loop in only 14% of them), so its left end had to be
  explained in a footnote. This axis makes the left end beta by construction.
- **Section c added**, the grouped column chart, in the lower-right slot the
  old per-fold ratio panel used to have.
- The helix numbers the old section carried are still printed on a rebuild,
  marked `[helix, not drawn]`, so the claims above can be rechecked: medians
  0.800 library and 0.838 ipTM, and 19.4% of the library against 9.0% of ipTM
  below 60% helix.

### 2026-11-11
- **a removed** - the map of which rule holds each backbone, over generator
  territory outlines. With it went the `counterfactual3.csv` and
  `fold_map_scatter.parquet` reads, the union/state classification, the
  territory contouring and the local `blur`.
- **c removed** - the per-fold log2 abundance ratio by dominant generator.
- **The difference strip under the violin removed.** It carried the figure's x
  axis, so the violin draws its own ticks and axis label now.
- **The two halves swapped.** ipTM is the upper half, library the lower.

All of it is in git history if it is wanted back; the last version with a map
is the one committed under the old folder name `fig06d_selection_difference`.

## Numbers the cut sections carried
Kept here because they are cited elsewhere and the code that produced them is
gone from the script.

| | Backbones |
|---|---|
| selected by either rule | 11,559 |
| in both | 3,048 (26%) |
| ipTM-ranked only | 3,247 |
| diversity-optimised only | 5,264 |

Median log2 library:ipTM ratio per fold cluster, over the 64 of 139 pool folds
holding at least 20 backbones: RFD3+LiMPNN -0.33 (16 folds), BoltzGen +1.07
(43), BD1+LiMPNN +1.49 (2), HD+LiMPNN +1.56 (3).

## Deliberate choices
1. **No p-values anywhere.** Both arms are deterministic selections from one
   pool. There is no sampling error to test, so a p-value would be a category
   error. Everything here is descriptive.
2. **Both halves of the violin share one bandwidth**, which is the thing that
   makes the halves comparable. Density is normalised to sum to 1 per side and
   both are scaled by the same factor.
3. **The two violin labels sit on their own half**, not above and below the
   whole axes, so colour is not the only thing tying a label to a density.

## Rebuild
```
python code/Fig3d_selection_difference/make_figure.py
```
Reads `inputs/sse_arms.csv`, `inputs/funnel_real_gate_v2.csv` and
`inputs/counterfactual3.csv`. No cluster job. Writes
`output/Fig3d_selection_difference/selection_difference.{svg,pdf,png}`.

## Related
- `Fig3b_sse_cluster_map` is the source of the per-backbone secondary structure.
- `archive/unplaced_2026-11-11/code/fig06b_selection_counterfactual` is the
  absolute version: the pool and the two selections side by side, plus rare-fold
  retention and per-generator enrichment.
