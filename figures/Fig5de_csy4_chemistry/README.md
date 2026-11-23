# Fig5d,e — interface chemistry runs toward natural Csy4

Placed by Leo 2026-11-13 as **Fig 5 d and e**: `love_oriented` left (d), `curves`
right (e). Panels a-c of Fig 5 are not built yet, so the assembled page currently
carries only this row.

Bare panels from `make_figure.py`:

| file | what it is |
|---|---|
| `love_oriented.{svg,pdf,png}` | **the main panel.** Standardised-difference dumbbell, one row per feature, every feature flipped so **+x = toward Csy4**; non-binder and binder dots with 95% bootstrap intervals, Csy4 as a diamond. 105 × 50 mm |
| `features.{svg,pdf,png}` | the five score features one at a time, in raw units: library step histogram, the 32 as two dot rows with means, Csy4 as a line. 178 × 42 mm |
| `score.{svg,pdf,png}` | the Csy4-direction score: library density with the 32 riding **on** the curve, binders filled / non-binders open, Csy4's position. 105 × 46 mm |
| `curves.{svg,pdf,png}` | the same score as three overlapping densities (library, non-binders, binders), each mean a thin tick riding its own curve. 105 × 46 mm |

`score` and `curves` are two views of one quantity — pick one for the figure, not both.

**Removed 2026-11-13, do not resurrect without a reason:** `ladder` (the four-column
slopegraph) and `love` (the unoriented dumbbell). `love_oriented` supersedes both,
the script no longer emits them, and the stale output files are deleted. Both are in
git history if ever needed.

## Why the five panels are step histograms, not densities
Every one of these features is a count over a denominator, so they sit on a grid and
carry point masses. `aa_A` and `arom` are counts over the 80-residue design — a
0.0125 grid, 44 and 22 distinct values. The three interface shares are counts over
the interface and spike at exactly zero: **46.9% of the library has no Asp/Glu at the
interface, 51.6% no aromatic, 7.8% no Lys.** A KDE smears those spikes into humps and
invents density between grid points, and for two of these features **the bar at zero
is the finding**. Grid features get one bar per level; the ratio features get 26 bins
from 0 to the library's 99.5th percentile, so the longest tail is not drawn (nothing
of interest is out there — Csy4 is at 0 on two of them).

The starkest single panel is **Lys at interface**: Csy4 is at zero, **0 of 15 binders
and 1 of 17 non-binders** reach it, and it is the only feature where the binder mean
sits on the *far* side of the library from Csy4 (+0.20 SD). Csy4 is at the extreme
edge of the design range, matched by exactly one of the 32 — *not* outside it, which
an earlier draft of this README claimed.

⚠️ **Two rows are half-zero, which weakens their dumbbell summary.** `af3_if_neg` is
exactly 0 for 46.9% of the library (and 12 of the 32), `af3_if_arom` for 51.6% (14 of
the 32). A mean and an SD are poor summaries of a spike-at-zero variable, and the
dumbbell cannot show that the variable is one. For those two rows a proportion test
(fraction with none at all) is the better statistic, and the counts are already in
the main write-up: acid-free interface 8/15 binders vs 4/17 non-binders. This is the
main reason `features` is worth keeping somewhere even though it overlaps
`love_oriented`.

## What "Csy4-direction score" means

One number per design, in library standard deviations, saying how far it sits along
the direction natural Csy4 defines in interface chemistry. Built in four steps:

1. **Five features**, all interface chemistry, no fold and no epitope:
   `aa_A` (Ala fraction of the sequence), `arom` (F+W+Y fraction of the sequence),
   `af3_if_neg` (share of interface residues that are Asp/Glu), `af3_if_K` (share
   that are Lys), `af3_if_arom` (share that are aromatic).
2. **z-score each feature against the 12k library**, so the library is centred on 0
   with SD 1 on every feature and the five become comparable.
3. **Orient each feature by Csy4's own measured value** — +1 if Csy4 sits above the
   library mean, −1 if below:

   | feature | Csy4 | library mean | sign |
   |---|---|---|---|
   | `aa_A` | 0.069 | 0.181 | −1 |
   | `arom` | 0.085 | 0.029 | +1 |
   | `af3_if_neg` | 0.000 | 0.036 | −1 |
   | `af3_if_K` | 0.000 | 0.157 | −1 |
   | `af3_if_arom` | 0.100 | 0.037 | +1 |

4. **Score = the mean of the five signed z-scores.**

So **0 is an average library design and +1 means one library SD toward Csy4**,
averaged over the five. Pushing Csy4's own values through the same formula gives
**+1.342**. Measured: binders **+0.373**, non-binders **−0.031**, library 0 by
construction. The binders are therefore about **28%** of the way from the library
to the natural protein, and the non-binders are sitting on the library.

**What it is not.** Not fitted to the flow data: the five signs come from Csy4's
measured values and nothing else, which is exactly what makes the P = 0.010 test
legitimate instead of circular. Not a similarity score to Csy4 — fold, sequence and
epitope show no resemblance at all, and none of them is in the formula. And not
fully independent of the earlier exploratory sweep: aromatics and interface Asp/Glu
had already surfaced there as near-significant, so the feature *list* is partly
data-informed even though the *direction* is not. Say that in the caption.

## Why `love` replaces `ladder`
The slopegraph gave equal billing to a column that is 50 by construction (the
library, zero variance) and one that is a single protein, which implies a four-step
trajectory among comparable groups. In the love plot the library is the **origin**,
Csy4 is **one mark**, and the only things drawn as samples-with-uncertainty are the
two that are samples. Percentile rank also flattened effect size; SD units keep it.
It is the covariate-balance convention, so readers decode it without a key.

**Two features the ladder had to leave out are back**, because extra rows are free
in this form and both earn their place:

- **Asp/Glu at interface** — non-binders +0.74 SD, binders −0.23, Csy4 −0.87. The
  dumbbell straddles zero pointing at Csy4: non-binders move away from the natural
  protein, binders toward it. This is the acid–phosphate kill flag as a picture.
- **Lys at interface** — both groups sit **above** the library (+0.17, +0.20) while
  Csy4 is at **−1.56**. The one row where the designs move *away* from the natural
  protein, and the clearest statement of the Arg-not-Lys gap.

Base-contact fraction stays off this panel for the same reason as before (reverses
on ESM).

⚠️ **Most intervals cross zero.** At n = 15 and 17 only Ala fraction and mean helix
length have a binder interval clear of the library, and no single row is significant
after correction. The panel's claim is the pattern across rows plus the composite
P = 0.010 in `score` / `curves`, never an individual row. Caption accordingly.

### The bandwidth rule in `curves`
All three densities use **one bandwidth**, the Silverman width for n = 15 (the
smallest group) computed from the **library's** spread (0.219). Per-group Silverman
would hand the n = 15 curve a bandwidth several times the library's purely because
it is small, so the binder curve would look broad for a reason that has nothing to
do with binders. With a common width a difference in shape is a difference in the
data. The library curve is mildly over-smoothed in exchange; its true shape is in
`score`.

⚠️ **The tail bumps are single designs.** At n = 15 and 17 the little hills in the
right tail of both small-group curves (binders near 2.05, non-binders near 1.5 and
0.85) are one or two designs each, not modes. That is what the rug under the axis is
for, and the caption must say it. If a reader is likely to over-read them, use
`score` instead — it cannot manufacture a mode because it plots the designs
themselves.

Input: `inputs/csy4_ladder_12k.csv.gz` (12,000 designs + one Csy4 row).
Built by `~/hits32_browse_flowFC/` → `export_ladder.py`; the underlying analysis
is `FLOW_FC_ANALYSIS.md` Appendix 2 in that folder.

## What the panels claim, and what they do not

**Binder = flow-cytometry fold change (mean POS / mean NEG) ≥ 1.2**, giving 15
binders and 17 non-binders out of the 32 ordered designs. This is a *different*
ground truth from the HTS enrichment the 32 were selected on — the two are
uncorrelated (Spearman −0.16, P = 0.40).

Csy4 = **4AL5 chain A** with chain B nt 5–20, which is the exact design target
(1.16 Å C1′ RMSD onto the design RNA frame), run through the **same `geom.py`
pipeline** as the 12k so every number is comparable.

The 32 share **no** fold similarity with Csy4 (TM 0.275 binders vs 0.276
non-binders vs 0.271 for 300 random library designs — all unrelated) and **no**
sequence similarity (Smith–Waterman no better than their own shuffled
sequences). They also read a different epitope: Csy4 puts 23.7% of its contacts
on nt16 and the designs put 0.7%. **Only interface chemistry lines up**, and
that is all these two panels show.

### Statistic in panel a
Mean **within-library percentile rank** of a group's members — not the
percentile of the group's mean. The library then sits at exactly 50 by
construction, ranks survive the skew and the ties (`af3_if_arom` is 0 for 51.6%
of the library, so anything median-based is unusable), and it is the same
statistic the permutation tests use. Error bars are 95% bootstrap intervals,
2,000 resamples.

### Two things the caption must not overstate
1. **Not five monotone ladders.** On this rank statistic the two aromatic
   features dip slightly below the library at the non-binder step (46.8, 49.3)
   before rising to 69.4 and 62.0. On group means that step reads as a small
   rise instead. It is noise either way; what the panel shows is the **binders
   displaced from both the library and the non-binders toward Csy4**.
2. **Base-contact fraction is deliberately absent** although it makes the same
   ladder on AF3 (0.365 / 0.384 / 0.449 / 0.561). It is the only feature of the
   set that reverses on the second refolder (ESM: 0.389 binders vs 0.403
   non-binders), so it is an AF3 artefact. Helix fraction and mean helix length
   do replicate on ESM; Ala and aromatic content are sequence features and
   model-independent.

The score in panel b is also partly AF3-dependent: its interface terms
(`af3_if_neg`, `af3_if_K`, `af3_if_arom`) are individually null on ESM. The
sequence half (Ala, aromatics) carries P = 0.038 on its own, the AF3-interface
half P = 0.032, the combination P = 0.010.

## Caption draft

> **Interface chemistry of the flow-validated binders runs toward natural Csy4,
> while fold, sequence and epitope do not.**
> **(a)** Five interface-chemistry features, each plotted as the mean percentile
> rank of a group's members within the 12,000-design ordered library (grey, 50 by
> construction): the 17 designs with flow fold change < 1.2 (dark grey), the 15
> with fold change ≥ 1.2 (blue), and Csy4 itself (navy diamonds, 4AL5 chain A put
> through the same pipeline). Bars are 95% bootstrap intervals (2,000 resamples;
> n = 17 and 15). The binders sit between the library and the natural protein on
> all five; the non-binder step is within noise.
> **(b)** Csy4-direction score — the mean of five library z-scores oriented by
> Csy4's own values (low Ala, high aromatic, low interface Asp/Glu, low interface
> Lys, high interface aromatic). Grey, the library distribution; points, the 32
> ordered designs with their group means as ticks; navy line, Csy4. The direction
> was fixed by the natural protein before the flow data were consulted.
> *P* from 20,000 label permutations within design stack.
> The same 32 designs show no fold similarity to Csy4 (mean TM-score 0.27, the
> value for random library members), no sequence similarity, and a different
> epitope.

## Rebuild

```
cd LV_ThesisFigures
~/.venvs/binder_map/bin/python code/Fig4i_csy4_ladder/make_figure.py
```

`run_all.sh` globs `Fig*/`, so it is picked up automatically.
Nothing is written to `final/`; no panel letters are drawn (STYLE_NOTES §11).
