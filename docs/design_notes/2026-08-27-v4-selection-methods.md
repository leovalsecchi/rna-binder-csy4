# Library selection on `passing_v4` — seven methods compared

**Date:** 2026-08-27 · **Pool:** 43,771 gated designs · **Target:** K = 12,000

---

## 1. The pool

`passing_v4.csv`, exported from the binder_filter tool against the **LRMSD-corrected**
store (6,821,957 designs).

**Gate:**

| tier | criterion |
|---|---|
| tier 1 | `sc_RMSD ≤ 3 Å` on **AF3 OR ESM** |
| tier 1 | `sc_DockQ ≥ 0.7` on **AF3 OR ESM** |
| tier 3 | `seq_run ≤ 6` (longest same-AA run) |
| tier 3 | `footnb_1.0Å ≤ 5 atoms` on **AF3 AND ESM** |
| top-level | `iptm_mean ≥ 0.6` |

Tier ops are AND throughout; the OR sits *inside* the rmsd and dockq groups, across
refolders. Boltz-2 is excluded everywhere — only AF3 and ESMFold2 are read.

**This differs from `passing_v1` in exactly one place:** the ipTM gate is `mean` rather
than `min`, which is what takes the pool from 29,992 to 43,771.

| stack | designs | share of pool |
|---|---:|---:|
| RFD3+LiMPNN | 22,443 | 51.3% |
| BG+LiMPNN | 11,779 | 26.9% |
| BG+BIF | 6,893 | 15.7% |
| BD1+LiMPNN | 2,087 | 4.8% |
| HD+LiMPNN | 569 | 1.3% |

---

## 2. Methodology

### 2.1 Quality score

`iptm_mean = mean(ipTM_AF3, ipTM_ESM)`. Mean, not min — minima are extrema and amplify
seed noise. ipTM is the only quality-bearing column **not** in the gate, so it keeps full
range inside the pool; DockQ and sc_RMSD are both truncated by the gate and cannot rank.

### 2.2 The three metric axes

| axis | definition | source |
|---|---|---|
| `rog` | Cα radius of gyration of the **designed** backbone | Foldseek CA traces |
| `dsasa` | interface ΔSASA, total | refold, sc_RMSD-winning backend |
| `fbase` | fraction of contacts to RNA base | refold, sc_RMSD-winning backend |

Where a metric exists per refolder, the value is taken from the backend with the better
`sc_RMSD` (= `dockq_LRMSD`, a pose measure) — a right-fold/wrong-place refold gives a
meaningless interface. `rog` is design-side, so it has no backend and no draw noise.

**All three are rank-normalised WITHIN STACK before bucketing.** This cancels per-stack
offsets (generator Rg range, sidechain provenance, backend-win skew) for free.

### 2.3 Stack allocation (used by every "nested" method)

Three phases, in order:

1. **Protected floor** — `min(1200, available)` per stack. Never redistributed away.
2. **Equal pass** — top up toward the base quota of 12,000 / 5 = 2,400.
3. **Quality pass** — leftovers by best ipTM as a true greedy on the merged sorted
   tails, capped at 1.5 × base.

Resulting quotas: BD1 **2,087** · BG+BIF **2,400** · BG+LiMPNN **3,420** ·
HD **569** · RFD3 **3,524**.

BD1 and HD receive their *entire* gated pools — there is no selection inside them at
all. HD at 569 sits below the 1,200 protected floor and near the ~500 "present but
mute" line.

### 2.4 Fold clustering

Foldseek all-vs-all on the **designed** backbones, `--alignment-type 1` (TMalign),
`-c 0.8`, `--max-seqs 300`. Similarity is `TM = min(qTM, tTM)`. Leiden
(RBConfiguration) over a k-NN graph, k = 15.

43,771 designs → **25,359 unique backbones** (Cα dedup) → 6,203,930 raw hits →
4,665,632 undirected pairs.

| granularity | TM floor | edges | Leiden res | clusters | singletons | largest |
|---|---:|---:|---:|---:|---:|---:|
| **fine (r80)** | 0.50, mutual | 92,640 | 8.0 | **1,688** | 1,196 | 461 bb |
| **coarse (~20)** | 0.30, non-mutual | 286,007 | 0.677 | **20** + 13 isolates | 13 | 3,151 bb |

Singleton count is identical at every fine resolution (r05→r80): those 1,196 backbones
have no edge above the floor, and resolution cannot split what has no neighbours. The
coarse partition therefore needs a *different graph*, not a lower resolution — dropping
the floor to 0.30 collapses isolates from 1,196 to 13.

---

## 3. The seven methods

| # | method | what it does |
|---|---|---|
| 1 | **random 12k** | uniform draw, seed 42. The null. |
| 2 | **top ipTM 12k** | pure quality ranking, no structure at all. |
| 3 | **grid 27 flat (proportional)** | 3 axes × 3 terciles = 27 joint cells over the whole pool; largest-remainder apportionment by cell size; fill each cell by ipTM. No stack logic. |
| 4 | **grid 27 flat (equal)** | same 27 cells, but an **equal** quota per cell. No stack logic. |
| 5 | **grid 27 (nested)** | method 3 run *inside* each stack's §2.3 quota. |
| 6 | **cluster r80 RR (nested)** | round-robin over the 1,688 fine fold clusters: rank within cluster by ipTM, then sort by rank, then ipTM. Nested in the stack quota. |
| 7 | **cluster20 RR (nested)** | same round-robin over the 20 coarse clusters. |

---

## 4. Scoring — the four numbers

| metric | definition | direction |
|---|---|---|
| **ipTM mean** | mean `iptm_mean` over the 12,000 | higher better |
| **stack CV** | coefficient of variation of the five per-stack counts | lower = more even |
| **tercile drift** | deviation from the **43,771 pool's** 33/33/33 split, in percentage points, over all 9 axis × tercile cells — reported as the mean and the worst cell | lower = split better preserved |
| KS | max deviation between the library's CDF and the pool's deciles, averaged over the 3 axes | lower = more representative |
| **designs/stack** | raw counts | — |

KS is a descriptive statistic here, not a test — no p-value is computed, and at n = 12,000
a real test would reject everything.

**It is reported both ways, and they disagree.** §5.1 is **within stack**: each stack's
library is compared to that stack's own pool deciles, then averaged over 5 stacks × 3 axes.
§5.1b is **global**: the whole 12,000 against the whole 43,771's deciles, no grouping.
Stack fairness and global axis uniformity are in direct conflict — the nested methods
re-weight the stacks on purpose — so the two measures rank the methods in opposite orders.
Quote whichever matches the question being asked, and say which one it is.

⚠️ **Caveat on the within-stack average:** it weights all five stacks equally, and every
nested method takes BD1 and HD *whole*, which forces KS = 0 on those two cells by
construction. Restricted to the three stacks where selection actually happens, the nested
grid is 0.0151 and random 0.0152 — a tie, not the win §5.1 appears to show.

---

## 5. RESULTS

### 5.1 Headline

**Everything here is measured against the whole 43,771-design pool.** The headline
representativeness number is **tercile drift**: each axis is cut at the pool's terciles, so
an unbiased 12,000 would put 33.3% in each; drift is how far the library actually lands
from that, in percentage points.

| method | ipTM mean | stack CV | tercile drift (mean pp) | worst cell | clusters (r80) | singletons | backbones |
|---|---:|---:|---:|---:|---:|---:|---:|
| **random 12k (base)** | 0.6875 | 1.019 | **0.3** | **0.8** | 853 / 1688 | 457 | **9,654** |
| **top ipTM 12k** | **0.7633** | 0.853 | 5.2 | 12.2 | 708 / 1688 | 324 | 7,797 |
| grid 27 flat (prop) | 0.7626 | 0.864 | 4.1 | 9.5 | 707 / 1688 | 329 | 7,807 |
| grid 27 flat (equal) | 0.7621 | 0.855 | 4.2 | 9.8 | 726 / 1688 | 337 | 7,829 |
| grid 27 (nested) | 0.7513 | **0.500** | 5.2 | 12.9 | 1034 / 1688 | 611 | 8,091 |
| **cluster r80 RR (nested)** | 0.7338 | **0.500** | 5.8 | 14.9 | **1688 / 1688** | **1,196** | 8,226 |
| cluster20 RR (nested) | 0.7399 | **0.500** | 5.5 | 12.4 | 1090 / 1688 | 670 | 8,122 |

**Only a random draw preserves the split.** 0.3 pp mean drift, worst cell 0.8 pp — the
12,000 is a faithful miniature of the 43,771 on all three axes. Every ipTM-driven method
lands at 4–6 pp mean and 9–15 pp in its worst cell.

**The damage is concentrated in ΔSASA HIGH.** It should hold 33.3% and instead holds
18.5–23.8% depending on method — roughly half. The mirror image is ΔSASA LOW at 38–46%.
Rg HIGH is the second casualty at 22.8–26.9%. The base-contact fraction survives almost
intact (30–37%), consistent with it being the near-inert axis.

⚠️ **None of the grid variants fix this, and that is a design flaw, not a surprise.** They
bucket on axes rank-normalised **within stack** (§2.2), so they equalise each stack's own
terciles — which is what §5.1b measures — while leaving the pool-level split untouched. To
preserve the 33/33/33 of the 43,771 you have to bucket on the **pool's** terciles. That
method has not been run; see §11.

### 5.1b KS by axis, and the within-stack alternative

Per-axis breakdown of the §5.1 global KS:

| method | rog | dsasa | fbase | **mean** |
|---|---:|---:|---:|---:|
| **random 12k (base)** | 0.0088 | 0.0050 | 0.0074 | **0.0071** |
| grid 27 flat (prop) | 0.0961 | 0.0996 | 0.0214 | 0.0724 |
| grid 27 flat (equal) | 0.1023 | 0.1022 | 0.0226 | 0.0757 |
| top ipTM 12k | 0.1128 | 0.1236 | 0.0318 | 0.0894 |
| cluster20 RR (nested) | 0.1369 | 0.1229 | 0.0226 | 0.0941 |
| grid 27 (nested) | 0.1264 | 0.1494 | 0.0304 | 0.1021 |
| cluster r80 RR (nested) | 0.1277 | 0.1567 | 0.0398 | 0.1081 |

**The within-stack variant ranks the methods in the opposite order** — random 0.015,
top-ipTM 0.046, nested grid 0.009 — because it asks a different question: is each stack's
slice representative *of that stack*, rather than of the pool. Two effects separate them,
and only one is a flaw.

**The shape effect** — every quality-driven method pulls toward low Rg and low ΔSASA.
Global terciles for top-ipTM: Rg LOW 42.0% / HIGH 26.4%, ΔSASA LOW 40.6% / HIGH 21.1%.
That is a real distortion, worth ~0.09 on its own.

**The stack-mix effect** — the nested methods move RFD3 from 51.3% to 29.4% and BD1 from
4.8% to 17.4% *by design*. Rg has a between-stack variance share of 0.589, so changing
the stack mix mechanically moves the global Rg distribution. The nested methods are not
failing to represent the pool; they represent a deliberately different population.

**Which number to quote depends on the question.** "Does the library look like the pool I
filtered?" → global, and only random does. "Within each stack, is the library
representative of that stack?" → within-stack, and the nested grid wins. They are not
competing answers to one question.

### 5.2 Designs per stack

Equal share = 2,400 = 20.0%.

| method | BD1+LiMPNN | BG+BIF | BG+LiMPNN | HD+LiMPNN | RFD3+LiMPNN |
|---|---:|---:|---:|---:|---:|
| POOL (43,771) | 2,087 (4.8%) | 6,893 (15.7%) | 11,779 (26.9%) | 569 (1.3%) | 22,443 (51.3%) |
| random 12k | 568 (4.7%) | 1,885 (15.7%) | 3,136 (26.1%) | 167 (1.4%) | 6,244 (52.0%) |
| top ipTM 12k | 560 (4.7%) | 2,260 (18.8%) | 4,180 (34.8%) | 247 (2.1%) | 4,753 (39.6%) |
| grid 27 flat (prop) | 589 (4.9%) | 2,190 (18.2%) | 3,999 (33.3%) | 241 (2.0%) | 4,981 (41.5%) |
| grid 27 flat (equal) | 600 (5.0%) | 2,213 (18.4%) | 4,041 (33.7%) | 247 (2.1%) | 4,899 (40.8%) |
| grid 27 (nested) | 2,087 (17.4%) | 2,400 (20.0%) | 3,420 (28.5%) | 569 (4.7%) | 3,524 (29.4%) |
| cluster r80 RR (nested) | 2,087 (17.4%) | 2,400 (20.0%) | 3,420 (28.5%) | 569 (4.7%) | 3,524 (29.4%) |
| cluster20 RR (nested) | 2,087 (17.4%) | 2,400 (20.0%) | 3,420 (28.5%) | 569 (4.7%) | 3,524 (29.4%) |

All three nested methods land on identical counts — the §2.3 allocation fixes them
before any within-stack selection runs.

### 5.3 Metric distribution — LOW / MID / HIGH vs the pool

Percent of the 12,000 in each tercile **of the 43,771-design pool**. An unbiased library
sits at 33.3 / 33.3 / 33.3. Figure: `fig3_terciles.png`.

| method | axis | LOW | MID | HIGH |
|---|---|---:|---:|---:|
| **random 12k (base)** | Rg | 32.6% | 34.0% | 33.4% |
| | ΔSASA | 33.3% | 33.0% | 33.8% |
| | base frac | 33.1% | 33.3% | 33.5% |
| top ipTM 12k | Rg | 42.0% | 31.6% | **26.4%** |
| | ΔSASA | 40.6% | 38.3% | **21.1%** |
| | base frac | 30.7% | 33.4% | 35.9% |
| grid 27 flat (prop) | Rg | 40.5% | 32.6% | 26.9% |
| | ΔSASA | 37.9% | 38.3% | 23.8% |
| | base frac | 31.7% | 33.7% | 34.6% |
| grid 27 flat (equal) | Rg | 40.8% | 32.7% | 26.5% |
| | ΔSASA | 38.1% | 38.3% | 23.6% |
| | base frac | 31.5% | 34.0% | 34.5% |
| grid 27 (nested) | Rg | 41.3% | 33.5% | 25.2% |
| | ΔSASA | **45.9%** | 33.7% | **20.4%** |
| | base frac | 30.9% | 33.4% | 35.8% |
| cluster r80 RR (nested) | Rg | 40.0% | 34.4% | 25.7% |
| | ΔSASA | 43.6% | 37.9% | **18.5%** |
| | base frac | 30.0% | 33.2% | 36.8% |
| cluster20 RR (nested) | Rg | 40.2% | 37.0% | **22.8%** |
| | ΔSASA | 41.0% | 38.0% | 21.0% |
| | base frac | 31.4% | 33.5% | 35.0% |

**Random reproduces the pool to within 1 pp on every one of the nine cells.** Nothing else
comes close.

**Every other method tilts the same way:** low Rg and low ΔSASA over-represented, high
under-represented. ΔSASA HIGH is the worst cell throughout, bottoming out at **18.5%** for
the r80 round-robin — barely half its share. This is the pool's own latent
ipTM ↔ ΔSASA correlation of −0.171 (§10) amplified by ranking.

⚠️ **This matters for the learning objective, not just for tidiness.** A library holding
18–24% high-ΔSASA designs instead of 33% is a weak instrument for testing whether large
interfaces bind better — the cells you would need for that comparison have been thinned by
a third to a half before the assay starts.

### 5.4 Pairwise overlap (/12,000)

| | random | top ipTM | grid flat prop | grid flat eq | grid nested | clust r80 | clust20 |
|---|---:|---:|---:|---:|---:|---:|---:|
| random | 100.0% | 27.5% | 27.6% | 27.6% | 27.4% | 27.4% | 27.2% |
| top ipTM | 27.5% | 100.0% | **94.0%** | **91.8%** | 81.1% | 62.3% | 65.6% |
| grid flat prop | 27.6% | 94.0% | 100.0% | 94.3% | 80.8% | 61.7% | 65.0% |
| grid flat eq | 27.6% | 91.8% | 94.3% | 100.0% | 81.1% | 61.8% | 65.1% |
| grid nested | 27.4% | 81.1% | 80.8% | 81.1% | 100.0% | 73.9% | 76.4% |
| cluster r80 | 27.4% | 62.3% | 61.7% | 61.8% | 73.9% | 100.0% | 70.2% |
| cluster20 | 27.2% | 65.6% | 65.0% | 65.1% | 76.4% | 70.2% | 100.0% |

Decision rule from the spec: **> 85% overlap ⇒ the method choice is immaterial.**

---

## 6. Conclusions

**The flat metric grid is top-ipTM wearing a hat.** 94.0% overlap with pure ipTM ranking,
the same KS (0.048 vs 0.046), the same ipTM to three decimals, the same stack CV. Whether
cells get proportional or equal quotas makes no difference either (94.3% overlap with
each other). **Metric bucketing without stack nesting buys nothing.** This reproduces the
earlier result on the v1 pool, where three separate metric methods collapsed onto the
baseline and one returned the *identical* 12,000.

**Nesting is where the grid earns its keep** — within-stack KS **0.009** against
top-ipTM's 0.046, a fivefold improvement, for 0.0120 ipTM, stack CV 0.853 → 0.500 and
+326 fold clusters.

⚠️ **It does NOT beat random, and an earlier draft of this document said it did.** Two
corrections, both material:

1. The within-stack average weights all five stacks equally, and every nested method takes
   **BD1 and HD whole** (2,087/2,087 and 569/569). When the library *is* the pool, KS is 0
   by construction — so two of five cells are free zeros. On the three stacks where
   selection actually happens, the nested grid scores **0.0151 against random's 0.0152**.
   A tie, not a win.
2. Measured **globally** (§5.1b), the ranking inverts outright: random 0.0071, nested grid
   0.1021. The nested methods re-weight the stacks on purpose, and Rg is 59% stack
   identity, so changing the mix mechanically moves the global Rg distribution.

**Random is a much stronger diversity baseline than it looks** — the only method that
reproduces the pool globally (KS 0.007, an order of magnitude ahead of everything else),
the most unique backbones (9,654), and more fold clusters than top-ipTM (853 vs 708). It
pays **0.0758 ipTM**, by far the largest quality drop here, and it inherits the pool's
stack imbalance wholesale (CV 1.019, RFD3 at 52.0%).

**Fold-cluster coverage is the axis nothing else sees.** The quality-ranked control
reaches only **708 of 1,688 clusters (42%)** and 324 of 1,196 structural singletons — a
top-ipTM library silently omits **58% of the fold space**. Only the round-robin reaches
100%, at a cost of 0.0295 ipTM. Unique-backbone count barely separates the nested methods
(8,091 / 8,122 / 8,226), so what the clustering buys is **fold-space coverage, not
backbone count**. No combination of Rg / ΔSASA / base-fraction can substitute: those are
continuous scalars, and fold identity is not.

**The coarse 20-cluster partition is not a substitute for the fine one.** It reaches 1,090
of 1,688 fine clusters (65%) — better than the grid but well short of the round-robin —
and it is the *only* method with a visibly uneven tercile profile (ΔSASA LOW 29.5% /
MID 36.0%). With the largest coarse cluster holding 3,151 backbones, one round of
round-robin over 20 clusters is simply too blunt to act as a coverage constraint.

**Stack CV 0.500 is a floor, not a result.** It is identical for all three nested methods
because BD1 (2,087) and HD (569) receive their entire pools; no allocation scheme can
equalise further. **HD is effectively mute** at 569 designs — below the 1,200 protected
floor and near the ~500 line where a stack is present but statistically silent.

---

## 7. Caveats

- **`fbase` is a near-inert axis.** It is `frac_base_contacts`, measured at cross-backend
  r = 0.16 (AF3 vs ESM on the same design) — close to uncorrelated, so it functions as a
  near-random third coordinate. Every method including random lands at KS 0.006–0.025 on
  it. The real ΔSASA base fraction exists but covers only ~1.3% of this pool.
- **Best-of-2 backend selection is optimistically biased** (winner's curse). Named, not
  corrected.
- **No per-uuid or per-sequence-cluster cap** — deliberate. Multiple LiMPNN sequences on
  one backbone is the only way to learn whether sequence choice matters independently of
  backbone.
- **KS is descriptive.** No hypothesis test is implied.
- Cluster counts are over **backbones** (25,359), not designs (43,771). A per-cluster
  quota buys backbones; extra sequences on one backbone are sequence diversity at zero
  fold diversity.

---

## 8. Provenance

| artefact | path (Euler) |
|---|---|
| pool | `analyses/library_selection/passing_v4.csv` |
| method comparison | `analyses/library_selection/v4_methods.py` |
| coarse sweep | `analyses/library_selection/v4_coarse.py` |
| fold map | `analyses/foldseek_map/results_v4/design/` |
| fine clusters | `.../clusters/fold/clusters.parquet` (`cluster_r05…r80`) |
| coarse clusters | `.../clusters/fold_coarse/clusters.parquet` (`cluster_coarse`) |

Python: `/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python` (default python3 has
no pyarrow). Foldseek v8-ef4e960 via
`module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960`.

---

## 9. Figures

PNG, 200 dpi, Okabe-Ito CVD-safe palette matching `analyses/plots_tools/config.py`.
Local: `~/Documents/binder_filter/exports/v4_figs/` · Euler:
`analyses/library_selection/v4_figs/` · generator: `v4_plots.py`.

| file | figure | what it shows | goes with |
|---|---|---|---|
| `fig1_stack_share.png` | **Designs per stack** | 100% stacked bars, one row per method plus the pool, counts printed in-bar, dotted line at the 20% equal share | §5.2 |
| `fig2_tradeoff.png` | **Two trade-off scatters** | ipTM mean vs KS (left) and vs fold-cluster coverage (right); blue = nested in the stack allocation, orange = no stack logic | §5.1 |
| `fig3_terciles.png` | **LOW / MID / HIGH** | three panels (Rg, ΔSASA, base fraction), stacked tercile bars per method, dotted lines at 33.3 / 66.7% | §5.3 |
| `fig4_cluster_coverage.png` | **Fold-space coverage** | paired bars — clusters reached (of 1,688) and structural singletons (of 1,196) | §6 |
| `fig5_overlap.png` | **Pairwise overlap** | 7×7 heatmap, % of the 12,000 shared between every pair | §5.4 |

**Suggested slide order:** fig1 (the fairness problem) → fig3 (what the grid actually
controls) → fig2 (the trade-off, both panels) → fig4 (what only clustering buys) →
fig5 (why the three unnested methods are one method).

Two figures carry the argument on their own. **fig2 left** shows the three unnested
methods collapsed into a single orange blob at the bottom-right — near-identical ipTM,
near-identical KS. **fig4** shows the coverage gap: 708 clusters for top-ipTM against
1,688 for the round-robin.

---

## 10. Slide: "Basecase" — the starting pool

*Figure: `fig0_pool_profile.png`. Everything below is the 43,771-design pool before any
selection runs, on the same dimensions used to score the methods.*

### Headline

**43,771 designs · 25,359 unique backbones · 1,688 fold clusters · 5 stacks · campaigns v2–v9**

### Composition — stack CV **1.008**

| stack | designs | share | backbones | designs/backbone | clusters touched | ipTM med | Rg med | ΔSASA med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| RFD3+LiMPNN | 22,443 | 51.3% | 11,970 | 1.87 | 573 | 0.670 | 17.26 | 2,377 |
| BG+LiMPNN | 11,779 | 26.9% | 4,786 | 2.46 | 382 | 0.694 | 11.54 | 1,708 |
| BG+BIF | 6,893 | 15.7% | 6,893 | 1.00 | 442 | 0.691 | 11.86 | 1,647 |
| BD1+LiMPNN | 2,087 | 4.8% | 1,141 | 1.83 | 517 | 0.678 | 16.05 | 1,678 |
| HD+LiMPNN | 569 | 1.3% | 569 | 1.00 | 161 | 0.707 | 12.38 | 1,624 |

Cluster counts sum to 2,075 against 1,688 distinct clusters, so some folds are reached by
more than one stack — but most are stack-specific.

### Metric ranges

| axis | min | Q1 | median | Q3 | max | tercile cuts |
|---|---:|---:|---:|---:|---:|---|
| ipTM (mean AF3/ESM) | 0.600 | 0.641 | **0.679** | 0.726 | 0.896 | 0.653 / 0.709 |
| Rg, designed backbone (Å) | 10.19 | 11.88 | **15.43** | 17.47 | 28.33 | 12.80 / 16.72 |
| interface ΔSASA (Å²) | 755 | 1,662 | **1,973** | 2,400 | 3,824 | 1,765 / 2,249 |
| base-contact fraction | 0.011 | 0.291 | **0.360** | 0.433 | 0.901 | 0.316 / 0.407 |

### Structure of the pool

| quantity | value |
|---|---|
| designs per backbone | mean 1.73 · **68.1% of backbones carry exactly one** · max 11 |
| fold clusters (r80) | 1,688 · **1,196 singletons (71%)** · largest 461 backbones |
| fold concentration | the 5 largest clusters hold only **7.2%** of backbones |
| between-stack variance share | Rg **0.589** · ΔSASA 0.387 · base fraction 0.029 |
| correlation with ipTM | Rg **−0.135** · ΔSASA **−0.171** · base fraction +0.047 |

### Talking points

1. **One stack owns half the pool.** RFD3 51.3% against HD 1.3% — a **39× spread**. Every
   fairness decision downstream starts from this imbalance, and it is the imbalance, not the
   selection algorithm, that dominates the outcome.
2. **HD is already below the statistical floor before selection begins.** Its entire gated
   pool is 569 designs, under the 1,200 protected floor and near the ~500 line where a stack
   is present in the library but statistically mute.
3. **Rg is bimodal, and the two modes are the stacks.** BG sits at ~11.5 Å, RFD3 at ~17.3 Å;
   between-stack variance share **0.589**. "Spread the library over Rg" and "balance the
   stacks" are very nearly the same instruction. The base-contact fraction, at 0.029, carries
   essentially no stack signal at all.
4. **ipTM is already truncated by the gate.** Floor 0.6 by construction, IQR only 0.085 wide.
   It is still the widest quality column available — DockQ and sc_RMSD are gated harder — but
   its ranking power is limited.
5. **The skew that appears in every selected library is already latent here.** ipTM correlates
   **negatively** with both Rg (−0.135) and ΔSASA (−0.171). Weak, but signed: rank on ipTM and
   you pull toward smaller binders that bury less surface. That is the mechanism behind the
   global KS results.
6. **Fold space is flat and long-tailed, not dominated.** 71% of clusters are singletons and
   the five largest hold just 7.2% of backbones. There is no big redundant fold family to
   prune — coverage is entirely a question of reaching the tail.
7. **Sequence redundancy is modest and uneven.** 1.73 designs per backbone overall, but
   BG+BIF and HD sit at exactly 1.00 (BoltzIF is all-atom generative, so each design is its
   own backbone; HD's gate is tight enough that only one sequence per backbone survives),
   against 2.46 for BG+LiMPNN.

---

## 11. The method that is missing

Every bucketing method here rank-normalises its axes **within stack** before cutting
terciles (§2.2). That was the right call for stack fairness — it cancels per-stack offsets
in Rg range and backend-win skew — but it means the 27 cells equalise *each stack's* terciles
and say nothing about the pool's.

**To preserve the pool's 33/33/33, cut the terciles on the raw pooled values instead**, and
apportion the 27 cells equally, either flat or nested inside the §2.3 stack quotas. Expected
behaviour: tercile drift falls toward random's 0.3 pp while ipTM stays well above random's
0.6875, because within each cell designs are still taken best-ipTM-first.

Not yet run. It is a small edit to `v4_methods.py` — `cell_key()` switches from `r_<ax>`
(within-stack rank) to a `pd.qcut` on the pooled raw axis — and about ten minutes.
