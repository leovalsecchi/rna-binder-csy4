# 12k library selection — method and findings

Handover for slide-building. Csy4 RNA-binder campaign, one-shot wet-lab library.
Everything below is measured, not projected. Where a number is uncertain or a caveat
applies it says so explicitly — please keep those caveats on the slides.

---

## 1. The problem

Pick **exactly 12,000 designs** from a gated pool for a **one-shot** wet-lab library
(synthesised and assayed in one batch, no feedback loop). Two objectives at once:

1. **Hits** — maximise the chance that something binds.
2. **Learning** — find out which metrics and which generator stacks predict binding.

These pull in opposite directions. Maximising hits means taking the highest-scoring
designs; learning means deliberately spending wells on designs you are less sure about,
so the assay can tell you what matters.

Five generator stacks: `RFD3+LiMPNN`, `BG+LiMPNN`, `BG+BIF`, `BD1+LiMPNN`, `HD+LiMPNN`.

---

## 2. The pool

Export `passing_v1.csv` — **29,992 designs**, gate:

```
(LRMSD_af3 <= 3  OR  LRMSD_esm <= 3)
AND (DockQ_af3 >= 0.7 OR DockQ_esm >= 0.7)
AND seq_run <= 6
AND footnb_010 <= 5 on BOTH backends        # <=5 foot atoms within 1.0 A
```

Plus a top-level **ipTM gate: `min(AF3, ESMFold2) >= 0.6`** — both backends must clear
0.6, not their average. Note this differs from the design spec, which argued for the *mean*
on the grounds that a minimum is an extremum and amplifies noise. Ranking inside both
libraries uses `iptm_mean`, so the gate and the ranker currently use different statistics.
Boltz-2 excluded throughout — it is the least reproducible backend while reporting the
highest ipTM.

346 designs (1.2%) are not in the Foldseek map and were **excluded from both libraries**
so the two are directly comparable. Working pool: **29,646 designs on 18,138 backbones**.

**K/N = 40.5%** — we keep two designs in five. This matters more than any algorithm choice
(see Finding 1).

| stack | pool | share |
|---|---|---|
| RFD3+LiMPNN | 13,881 | 46.8% |
| BG+LiMPNN | 8,833 | 29.8% |
| BG+BIF | 4,797 | 16.2% |
| BD1+LiMPNN | 1,651 | 5.6% |
| HD+LiMPNN | 830 | 2.8% |

---

## 3. The fold clustering

Foldseek all-vs-all on the **designed backbone** (the generator's own output, *not* a
refold), `TM = min(qTM, tTM)`, Leiden at fixed seed. Designs are first deduplicated by a
hash of their Cα coordinates, so the several LiMPNN sequences that share one backbone
collapse to a single node — 29,646 designs become 18,138 backbones.

Two runs exist. The **fine** run (`TM >= 0.50`, mutual kNN) is the usable one:

| run | clusters | singletons | largest cluster | top-10 hold |
|---|---|---|---|---|
| coarse TM0.40 / r05 | 104 | 73 | 4,162 designs | 87.7% |
| coarse TM0.40 / r10 | 117 | 73 | 1,472 | 63.5% |
| fine TM0.50 / r10 | 1,236 | 940 | 921 | 30.8% |
| **fine TM0.50 / r80 (chosen)** | **1,435** | **941** | **293** | **13.4%** |

The coarse run is unusable for selection: one cluster holds 14% of the pool, so a
per-cluster quota would be meaningless. **r80** was chosen for the flattest profile.

**Why the clusters matter more than the metric axes:** cluster identity is a property of
the *design*, carrying zero refold noise, and it tracks how a design was built. See
Finding 3.

---

## 4. Methods evaluated

Fourteen variants were run. They fall into three groups.

**Group A — cluster-based (all worked)**
- *cluster round-robin* — round 1 takes the best-ipTM design from every cluster, round 2
  each cluster's second-best, and so on to 12,000
- *cluster floor + grid* — reserve one design per cluster, then fill the rest by a 27-cell
  joint metric grid
- *cluster floor + terciles* — same floor, filled by 3 axes × 3 terciles

**Group B — metric axes nested inside a stack allocation (worked, but the allocation did
the work)**
- terciles, 27-cell grid, and kNN density deciles, each run within per-stack quotas

**Group C — metric axes with no stack layer (all failed)**
- flat 27-cell grid, equal-per-cell and proportional
- pure kNN density thinning
- thirds-per-axis decile union

Group C is the informative one — see Finding 2.

---

## 5. Findings

### Finding 1 — At 40% retention, metric diversity comes for free

Tercile coverage on every axis is even **without any metric bucketing at all**. Measured
within stack, where even the naive top-12k-by-ipTM lands near 33/33/33:

| method | Rg L/M/H | ΔSASA L/M/H | base-fraction L/M/H |
|---|---|---|---|
| base (top-12k ipTM) | 34.8 / 32.9 / 32.2 | 31.7 / 35.6 / 32.8 | 32.0 / 34.9 / 33.1 |
| cluster round-robin | 32.3 / 31.2 / 36.6 | 30.0 / 33.8 / 36.1 | 32.0 / 34.6 / 33.3 |

**The elaborate bucketing machinery was solving a problem that does not occur at this
scale.** You cannot concentrate hard enough on 12,000 of 29,646 to abandon a tercile.
Bucketing still buys *refinement* — it tracks the pool's shape more closely inside each
band — but not *coverage*.

### Finding 2 — Metric axes cannot balance the stacks; fold clusters can

Every method that partitioned only on metric axes, with no stack layer, reproduced the
naive baseline almost exactly:

| method | stack CV | clusters reached | verdict |
|---|---|---|---|
| pure kNN density | 0.849 | 623 | ≈ base case |
| thirds-per-axis decile | 0.847 | 624 | **identical to base case** |
| flat 27-cell grid | 0.872 | 691 | ≈ base case |
| cluster round-robin | **0.337** | **2,311** | works |

*(this table is from the earlier, wider pool of 71,866 where all 14 methods were compared)*

The reason: **the metric axes do not separate generators.** RFD3 sits at 44–58% of every
density decile and is spread across most of the 27 grid cells. 1,435 fold clusters do
separate them, because fold identity tracks how a design was built.

The thirds-per-axis method is the cleanest demonstration — it produced *exactly* the same
12,000 designs as ranking by ipTM. Any method that allocates less than K/N per stratum and
then tops up the remainder by global ipTM collapses to the baseline, because the top-up
undoes the diversity work.

### Finding 3 — Only design-derived metrics are reliable

Cross-backend agreement, AF3 vs ESMFold2, **same metric on the same design** (n ≈ 100k):

| metric | global r | within-stack r |
|---|---|---|
| ΔSASA total | 0.590 | 0.36 – 0.60 |
| ΔSASA base fraction | 0.293 | 0.17 – 0.28 |
| base-contact fraction | **0.159** | **0.08 – 0.17** |
| H-bonds at interface | 0.219 | 0.13 – 0.22 |

Interface metrics measured on refolds are **weak descriptors of the design** — the value
depends heavily on which model folded it. Within stack, where the selection actually
operates, agreement is worse still.

**Rg and fold-cluster identity are the only exceptions**, because both are computed on the
designed backbone and have no backend to disagree with.

This is a harsher test than the seed-to-seed reliability previously measured (ICC
0.78/0.63/0.87), which compared two draws from the *same* model. Caveat, do not overstate:
low cross-backend agreement may partly reflect genuinely different interfaces rather than
measurement error — but either way the value fails as a stable per-design attribute.

### Finding 4 — Designs grip the backbone, not the bases

From the 100k reference set, RNA-side buried surface split by moiety:

**~30% base / ~70% backbone** — designs are roughly 2:1 backbone-grippers. Of the base
portion, ~70% is major-groove edge.

| stack | backbone % | base % | ΔSASA_rna | phosphate |
|---|---|---|---|---|
| BG+LiMPNN | 65.1 | 34.9 | 661 | 276 |
| BG+BIF | 65.2 | 34.8 | 632 | 261 |
| BD1+LiMPNN | 65.4 | 34.6 | 839 | 323 |
| HD+LiMPNN | 67.7 | 32.3 | 754 | 337 |
| RFD3+LiMPNN | **71.4** | **28.6** | **869** | **401** |

**RFD3 builds the biggest interface and reads the least base with it** — the excess is
almost entirely phosphate (401 Å²), i.e. large charge-driven backbone contacts. The rank
holds independently in all three backends, so it is not a refold artefact.

### Finding 5 — Both libraries under-sample large, extended binders

Measured *globally* (not within stack), against tercile cuts of the pool:

| axis | tercile cuts | base 12k L/M/H | cluster 12k L/M/H |
|---|---|---|---|
| Rg | < 12.44 Å / 12.44–16.66 / > 16.66 Å | 41.1 / 33.2 / **25.7** | 41.0 / 35.6 / **23.4** |
| ΔSASA | < 1,673 Å² / 1,673–2,069 / > 2,069 Å² | 35.8 / 37.4 / **26.8** | 36.7 / 38.3 / **25.0** |
| base fraction | < 0.292 / 0.292–0.388 / > 0.388 | 30.9 / 35.0 / 34.2 | 29.6 / 34.8 / 35.6 |

ipTM correlates negatively with Rg (−0.136) and ΔSASA (−0.173), so ranking by ipTM inside
any stratum favours compact binders with smaller interfaces. Both libraries inherit it.
The within-stack view hides this because Rg is more than half stack-determined.

**Consequence for the analysis:** comparing stacks is fine. Asking "does interface size
predict binding" is weakened — the large-interface end is tested with about two-thirds the
designs you would want.

---

## 6. The two libraries delivered

Both are 12,000 designs from the same 29,646-design pool. Neither imposes a stack quota.

| | **A — base** | **B — cluster round-robin** |
|---|---|---|
| rule | top 12,000 by ipTM | round *r* takes each cluster's *r*-th best design |
| BD1+LiMPNN | 593 (4.9%) | 921 (7.7%) |
| BG+BIF | 2,115 (17.6%) | 2,635 (22.0%) |
| BG+LiMPNN | 4,469 (37.2%) | 4,306 (35.9%) |
| HD+LiMPNN | 447 (3.7%) | 726 (6.0%) |
| RFD3+LiMPNN | 4,376 (36.5%) | 3,412 (28.4%) |
| stack CV | 0.816 | **0.649** |
| ipTM mean | **0.7663** | 0.7486 |
| fold clusters | 808 / 1,435 | **1,435 / 1,435** |
| structural singletons | 366 / 941 | **941 / 941** |
| unique backbones | 7,680 | **8,295** |

**Overlap 69.0%** — 8,280 designs shared, ~3,700 different.

The trade in one line: **B costs 0.018 mean ipTM and buys the 627 fold clusters and 575
structural singletons that A never touches.** B runs 39 rounds; round 1 places 1,435
designs (one per cluster, 12% of the library).

**Neither reaches an even stack split.** HD at 447/726 and BD1 at 593/921 are both well
below a 2,400 equal share. Adding a per-stack quota layer would fix this and was measured
to work, but was not applied here.

---

## 7. Caveats to keep on the slides

1. **A residue-register bug inflated sc_RMSD for four of five stacks. Now fixed.**

   The generators number the RNA target differently — RFD3 renumbers it 1-16, BoltzGen keeps
   2-17, BoltzDesign and HD emit a complete 18-mer as 1-18 — while every refold is folded
   from a bare 16-mer and comes out 1-16. The scoring code paired atoms on
   `(res_id, atom_name)`, so for the four non-RFD3 stacks it matched nucleotide *i* to
   nucleotide *i* across two different numbering origins. On a quasi-periodic A-form helix
   that **still converges** (~4.9 Å rather than ~1.3 Å) and silently returns a superposition
   frame off by roughly one base-pair screw.

   **sc_RMSD is measured in that frame, so it inherits the error.** DockQ moved only as a
   downstream consequence — it is a composite of fnat, iRMSD and LRMSD, and only LRMSD
   changed. **RFD3 numbers its RNA 1-16, the same as the refolds, so it was never affected**
   — its count changes by exactly zero, which is the cleanest available check that the fix
   does what it claims.

   The per-design effect is **bidirectional**, not uniform inflation: correcting it makes
   54.7% of AF3 values *smaller* and 45.3% larger (mean 28.02 → 27.75 Å). A frame off by one
   screw does not degrade every pose equally, and for 3,169 designs it happened to score
   them *below* the 3 Å gate when they should not have been.

   **Status: fixed 2026-08-27.** The store had been patched for DockQ but not sc_RMSD; both
   `rmsd_af3`/`rmsd_esm` are now written from the corrected values (backups at
   `*.bak_pre_lrfix_20260827`) and verified post-fix at 100.00% by three-way comparison.
   Zero designs became NaN — the code's `MIN_REGISTER` guard never fires on this data.
   Fix in `code/scripts/lrmsd_rna_aligned.py` (`.bak_pre_register_fix_20260826`).

   **Measured impact on the pool** — the gate reproduces exactly at 29,992, so these are
   reliable:

   | stack | before | corrected | delta | lost |
   |---|---|---|---|---|
   | BD1+LiMPNN | 1,651 | 1,342 | −309 | 506 |
   | BG+BIF | 4,797 | 4,779 | −18 | 580 |
   | BG+LiMPNN | 8,833 | 7,883 | −950 | 1,630 |
   | **HD+LiMPNN** | **830** | **430** | **−400** | 453 |
   | RFD3+LiMPNN | 13,881 | 13,881 | **+0** | **0** |
   | total | 29,992 | 28,315 | −5.6% | 3,169 |

   Every one of the 3,169 losses is a real value that moved above 3 Å, not a dropped
   measurement. K/N moves 40.5% → 42.4%; the method comparison and all three findings are
   unaffected.

2. **`HD+LiMPNN` may not be viable as an experimental arm.** It falls from 830 to **430**
   under the corrected sc_RMSD — below the ~500 threshold where a stack is present but
   cannot support a hit-rate comparison. It is unresolved whether this is real or an
   artefact of how 18-nt designs (BD1 and HD, the two stacks that lose most) are handled
   by the fix. **Hedge any claim about the HD arm.**
3. **The base-fraction axis currently in use is the contacts version**, which the design
   spec rejected on reliability grounds. The ΔSASA version exists but covers only 1.3% of
   this pool; extending it costs ~43 CPU-h.
4. **Wanting "high base contacts" conflicts with diversity on that axis.** Diversity means
   covering low/mid/high; preferring high means skewing to one end. Currently neither
   library enriches it meaningfully (high tercile 34.2% and 35.6% against 33.3%).
   Selecting the upper tail of a metric with r = 0.16 across backends would mostly select
   refold luck.
5. **Statistical power.** At a 1% hit rate, 2,400 designs per stack detects a ~2× hit-rate
   difference and 1,200 detects ~2.5×. At a 0.1% hit rate, nothing is detectable at this
   library size. Decide consciously before ordering.

---

## 8. Files

| file | what |
|---|---|
| `binder_filter/exports/library_12k_base_iptm.csv` | library A, 12,000 × 93 cols |
| `binder_filter/exports/library_12k_cluster_rr_r80.csv` | library B, 12,000 × 93 cols |
| `binder_filter/exports/library_map.html` | interactive fold map, click-to-ChimeraX |
| `binder_filter/build_library_map.py` | rebuilds the map |
| `rna_binder/analyses/library_selection/*.py` | every method run, on Euler |

Both CSVs carry `cluster`, `node` (backbone id), `sel_round` and `sel_reason`, so every
design traces back to why it was selected.

**The map** (`library_map.html`) plots one dot per design on the Foldseek design-layer
UMAP, coloured by stack / library membership / ipTM / cluster. Clicking a dot posts its
fold cluster to a local bridge on port 45001; ChimeraX opens the design plus its AF3 and
ESM refolds superposed on the RNA, and cycling steps through the rest of that cluster.
Requires `python3 ~/Documents/rna_binder_gui/open_refolds_server.py` running.

---

## Suggested slide narrative

1. **The problem** — 12,000 wells, one shot, two competing objectives
2. **What we thought we needed** — elaborate metric-space bucketing
3. **Finding 1** — at 40% retention, metric coverage is free; the machinery was unnecessary
4. **Finding 2** — the negative result: metric axes cannot balance stacks, and one method
   returned *exactly* the baseline
5. **Finding 3** — only design-derived metrics are reliable; refold interface metrics
   agree at r = 0.16–0.59 across backends
6. **What we built instead** — fold clusters from the designed backbone
7. **The two libraries** — the 69% overlap, and what the 0.018 ipTM buys
8. **Finding 4** — designs grip backbone not bases; RFD3 is the extreme case
9. **Caveats** — the register fix, HD's viability, and the power limit
