# Library selection: 12k from 40k, diversity-maxxed across 4 metric axes

Date: 2026-08-20
Status: design agreed; numeric parameters pending the gated CSV
Related: [[rna-binder-interface-metrics]], [[rna-binder-seed-selfconsistency]],
[[rna-binder-hairpin-foot-clash]], [[rna-binder-consensus-iptm-stale-min2]]

## 1. Problem

Select K = 12,000 designs from N ≈ 40,000 that pass the hard gate, for a
**one-shot wet-lab library**. No feedback loop: the library is synthesised and
assayed in a single batch.

Two objectives, both real:
- **Hits** — as many working binders as possible.
- **Learning** — the hit/no-hit data should tell us which metrics and which
  generative stack are worth trusting next campaign.

K/N ≈ 30%. This is a high retention fraction: the task is not sparse sampling
but **density equalisation** — thin the crowded regions, keep the sparse ones.

### 1.1 What the gate has already done

The input CSV arrives **already filtered** on **RMSD, DockQ, aa-repeats, and
hairpin-foot clash**. Those four are settled questions by the time this pipeline
runs. **They are therefore not diversification axes** — spreading across a
quantity that has already been truncated wastes budget on a coordinate whose
range no longer means anything.

This has a second consequence, developed in §3.1b: it removes the obvious
candidates for the within-cell quality score, because everything discriminating
is now either gated or is itself an axis.

### 1.2 K = 12,000 is EXACT

Confirmed by the user 2026-08-20: 12,000 is a hard synthesis constraint, not a
target. Four consequences, all load-bearing:

1. **Integer apportionment is mandatory.** 12,000 splits cleanly across stacks
   (2,400 at five, 3,000 at four), but a stack quota does not divide evenly into
   cells — 2,400 / 81 = 29.63. Per-cell quotas must use **largest-remainder
   apportionment**: floor each cell's share, then allocate leftovers to the
   largest fractional remainders, with a **deterministic tiebreak on cell index**
   so the run reproduces. Naive rounding drifts and the drift accumulates.
2. **The redistribution ledger becomes an ASSERTION.** Slots out must equal slots
   in, and the final count must equal 12,000 exactly. A hard failure, not a
   diagnostic noticed after the fact.
3. **The primary method needs a search to land exactly.** Union-then-deduplicate
   (§5.5) produces an unpredictable count, so N per bucket is binary-searched to
   the stack quota and topped up with the highest-ipTM unselected designs. Both
   steps are deterministic; the tiebreak on `design_id` is what makes the run
   reproducible.
4. ⚠️ **The 1.5× cap yields to the count** (but never to a protected floor —
   see §5.4). §5.4 caps any stack at 1.5× its base
   quota. If several stacks come in short, capped headroom on the remainder may
   not reach 12,000 — the cap and the exact count then contradict. **The count
   wins:** relax the cap **uniformly** across eligible stacks until 12,000 is
   reachable, and **log the relaxation loudly**. The cap is a guardrail; the
   count is a synthesis constraint. A relaxation must never happen invisibly.

## 2. Priority order (set by the user)

1. **Fair distribution across stacks** — the most important distribution.
2. **Interpretability** — every selection decision must be explainable.
3. Coverage of the metric axes below.

## 3. Axes

**Three** bucketing axes, plus ipTM as the quality score (§3.1b).
All rank-normalised within stack before use (§5.2).

| axis | definition | why |
|---|---|---|
| `rg_backbone` | Cα radius of gyration of the **designed** backbone | binder shape / compactness |
| `dsasa_total` | total interface ΔSASA, **sc_RMSD-selected backend** (§3.1c) | interface size — the master confounder |
| `dsasa_base_fraction` | base ΔSASA ÷ total ΔSASA, same selected backend | base-reading vs backbone-reading **mode**, size divided out |

### 3.1 Decisions and rationale

**D1 — ipTM is mean(AF3, ESMFold2); Boltz-2 excluded.**
B2 is the least reproducible backend (seed-seed DockQ 0.363, std 0.153) while
reporting the highest ipTM (0.737) — least reliable, most confident. `mean` not
`min`: a minimum is an extremum, the same statistic that gave
`foot_tp_min_dist` its r = 0.166 cross-backend agreement. A mean of two draws is
strictly more stable than the lesser of two.

⚠️ **Do not use the stored `consensus_iptm`.** In `app_metrics_100k` it is a
stale min-of-*two* that excludes ESM entirely (matches `min(af3, b2)` 100% of the
time). Derive `iptm_mean = (iptm_af3 + iptm_esm) / 2` at load; verify the column
names against the real v2–v9 header first.

Backend ipTM scales differ, but stage ②'s rank transform removes scale entirely.
This would matter for a threshold; it does not matter for a coordinate.

**D2 — Rg is computed on the designed backbone, Cα only.**
Two properties follow. It carries **zero draw noise** (it is a property of the
committed pose, not of a refold), and it is **stack-unbiased**: all five stacks
have N/CA/C/O in their design files, whereas any sidechain-dependent quantity
would silently favour the three LiMPNN stacks, which are backbone-only
(RFD3+LiMPNN, BG+LiMPNN, BD1+LiMPNN). Binder length is 80 across every stack, so
Rg measures shape, not size.

Note Rg is strongly stack-associated globally (BG+BIF 12.0 → RFD3+LiMPNN 18.9).
Used flat it would largely re-encode stack. The nested design (§5.4) confines it
to within-stack variation, which is where it does real work.

**D3 — the mode axis is `dsasa_base_fraction` (base ÷ total), not raw base.**
Raw base and total are nested — base is ≈32% of total and both scale with
interface size, so they would correlate ≈0.8 and stage ③ would drop one,
collapsing the two ΔSASA axes into one. The fraction divides size out, leaving
`dsasa_total` to carry size and `dsasa_base_fraction` to carry binding mode as
genuinely separate coordinates.

⚠️ The fraction's own reliability is **unmeasured**. ICC(1) is known for raw
ΔSASA base (.78/.63/.87) and total (.87/.81/.93), but a ratio of two correlated
quantities can be less reliable than either part. This is computable from
`seed_metrics_merged_v2.csv` with **no new compute** — it already carries SASA
base and backbone across all 30k seed structures. See §7 item 8.

**D3b — Both ΔSASA axes are read from the refolder that best agrees with the
designed backbone**, per design, measured by `sc_RMSD` (§3.1c).

**D4 — One SASA family throughout.**
Biotite and Rosetta correlate at r = 0.80 with a 1.26 ratio. Mixing them
silently rescales the axis. Pick one, record which, never join across.

### 3.1b The quality score

**DECIDED 2026-08-20 (supervisor's design, user-confirmed): the objective is to
MAXIMISE `iptm_mean` = mean(AF3, ESMFold2) ipTM.**

ipTM is **not** a diversification axis. It is the quantity every bucket is ranked
by, descending. See §5.5 for the selection method built on it.

**What the method actually optimises**, stated plainly so the write-up is honest:
*maximise ipTM subject to a floor of N designs in every tercile of every axis.*
This is **not** diversity maximisation. The bucketing is a representation floor
that prevents ipTM-ranking from abandoning whole regions of the three axes.

**⚠️ CORRECTION 2026-08-20 — ipTM IS gated.** An earlier draft argued for ipTM on
the grounds that it was the only quality column the gate had *not* truncated.
**That was wrong.** The binder_filter export applies `iptm mean >= 0.6` in
threshold mode, so ipTM is range-restricted like DockQ and `sc_RMSD`.

Measured on `passing_footincomplete.csv` (n=41,890): `iptm_mean` min = 0.600
exactly, IQR = 0.086, full range 0.281. Against a mean-of-two seed noise of
≈0.025, the IQR spans only **≈3.4 noise units** — the middle half of the pool is
separated by about three times the amount a re-run would move a design. The
ranking carries real but limited signal.

It remains the *least* truncated of the candidates (DockQ is gated at ≥0.7,
`sc_RMSD` at ≤3), so the choice still stands — but on narrower grounds than
originally claimed.

**The accepted tradeoff.** Full range, weaker validity: ipTM correlates with
RNA-aligned DockQ at only ρ ≈ 0.21–0.27, does not certify pose (81% of
fold-agreeing designs still dock >10 Å apart), and B2 calls 27% of wrong folds
≥0.8 — which is part of why B2 is excluded (D1). The alternative was truncated
range with better validity. This is a deliberate choice, not an oversight.

⚠️ **Use a freshly computed `mean(iptm_af3, iptm_esm)`.** The stored
`consensus_iptm` is a stale min-of-*two* excluding ESM. See D1.

**Superseded, recorded so they are not re-introduced:**
`mean(dockq_self_af3, dockq_self_esm)` (gated); `sc_RMSD` ascending (gated);
epitope spread (rejected by the user). `sc_RMSD` retains **only** its §3.1c role
of choosing which backend's ΔSASA to read.

**Stack-level redistribution (§5.4) uses ipTM**, so "quality" has one meaning
throughout this document.

**ipTM is a model confidence score, not an affinity prediction.** No ranking in
this pipeline may be read as a strength ranking.

### 3.1c ΔSASA backend selection

Per design, take both ΔSASA values from whichever refolder best reproduces the
designed backbone, ranked by **`sc_RMSD` = `dockq_LRMSD_af3` / `_esm`**, the
RNA-aligned ligand RMSD.

**Why a refold at all.** The design files *cannot* supply these axes for three
of the five stacks: RFD3+LiMPNN, BG+LiMPNN and BD1+LiMPNN are backbone-only
(N/CA/C/O), because LigandMPNN is inverse folding and its side-chain packer is a
separate opt-in module the multif stage scripts never invoke. ΔSASA is
sidechain-dependent, so it must be measured on a full-atom structure — i.e. on a
refold. The only remaining question is *which* refold.

**Why sc_RMSD picks the right one.** It aligns the RNA and measures the binder
with no second fit — "did this refold put the binder where the design put it?".
Since the designed backbone is the pose actually committed to, the refold that
best reproduces it is the closest available stand-in for the structure that
cannot be measured directly. Fold RMSD (`ca_rmsd_vs_designed`) would be the wrong
criterion, and in any case exists only as a single column on the 100k app table,
not campaign-wide.

**Implementation requirements:**

1. **Rank-normalise within backend FIRST, then select.** Selecting the backend
   and rank-normalising afterwards mixes raw AF3 and ESMFold2 values that may
   sit on different scales, and the rank transform cannot undo that after the
   fact. Normalise each backend's ΔSASA to quantiles separately, then take the
   selected backend's **percentile**. Percentiles mix; raw values do not.
2. **Candidates are AF3 and ESMFold2 only.** Consistent with D1, and
   independently justified: B2 has the worst ΔSASA ICC (.81 vs .87/.93) and
   every production B2 refold is unseeded — one uncontrolled draw. Best-of-2
   also halves the selection-on-noise exposure versus best-of-3.
3. **Report which backend won, per stack** (deliverable 8). If ESMFold2
   systematically wins one stack and AF3 another, the ΔSASA axes partially
   encode stack identity — which matters because stack is the outer stratifier.

**This survives the gate (§1.1).** Selecting a backend by `sc_RMSD` is not
diversifying on RMSD: it is a *relative* comparison between two backends for one
design, and gating the absolute value does not flatten the difference between
them. ⚠️ But the gate's boolean structure matters — if RMSD was applied as
`[OR]` across af3/esm, some designs have only one backend passing and the
selector is forced. Record which op was used and how many designs were forced.

**Recorded property, not a defect:** best-of-2 selection makes the resulting
ΔSASA **optimistically biased** relative to a random draw — the same
winner's-curse structure measured on the seed replicates (−0.013 af3 / −0.049 b2
/ −0.009 esm). This is the bias chosen deliberately; it is named so it is not
rediscovered as a surprise.

⚠️ **Coupling with the tiebreak:** `sc_RMSD` is a **DockQ input** (LRMSD is one
of DockQ's three terms), so backend selection is mildly correlated with the
within-cell ranking (§3.1b) — the same quantity now fills both roles. This is a
deliberate single point of dependence rather than a hidden one: it governs
*which* ΔSASA measurement is read and *which* design wins a slot, but not what
the axes themselves measure. `sc_RMSD` is also gated (§1.1) and therefore
range-restricted in both roles.

### 3.1d Side-chain packing — the alternative that would retire §3.1c

§3.1c exists only because three of five stacks lack sidechains in their design
files. **That is fixable**: LigandMPNN's `--pack_side_chains` places sidechains
on an existing backbone+sequence with no re-design, and
`postprocess_halludesign.sh` already runs it for HD via
`number_of_packs_per_design`.

**Cost.** The per-design figure is already measured inside the HD postprocess
step — pull it from there rather than estimating. Upper bound from campaign data:
the seed run cost 22.19 GPU-h for 30,000 *full structure predictions*
(2.7 s/structure), and packing is strictly cheaper than a prediction. **40k packs
is single-digit GPU-h.** The lever is scope, not speed — this is the gated 40k
(or just the selected 12k), not the 5.34M campaign.

**Pack ONLY the three that are missing** — RFD3+LiMPNN, BG+LiMPNN, BD1+LiMPNN.
(An earlier draft argued for re-packing all five for uniform provenance. That was
wrong, for three reasons:)

1. **HD is already LigandMPNN-packed**, so "all five" was really "the three plus
   a re-pack of BG+BIF". HD needs nothing.
2. **BG+BIF's sidechains are part of what BIF designed.** BoltzIF is generative
   all-atom and places sidechains jointly with the backbone. For the three
   unpacked stacks, packing *adds* information that never existed; for BG+BIF it
   would *replace* information the generator committed to — measuring a different
   molecule than the one BoltzIF produced.
3. **Within-stack rank normalisation (§5.2) already cancels per-stack
   provenance offsets** for selection purposes, at zero cost.

What remains of the provenance concern is narrow: cross-stack *interpretation*.
Per-stack distributions appear side by side in the §8 deliverables, and a
provenance offset would make BG+BIF's ΔSASA not commensurable with the
LiMPNN-packed stacks'. Mitigation is a diagnostic, not a re-pack: **record
sidechain provenance as a column** and check once whether BG+BIF's ΔSASA
distribution is offset.

**What it buys:** §3.1c retires entirely — no backend selection, no
winner's-curse bias, ΔSASA becomes a property of the committed pose with **zero
draw noise**, uniform across stacks.

**What it does not fix:** packed sidechains are *predicted*, not observed — the
packer's guess replaces the refolder's guess, though conditioned on the exact
committed backbone, which is the stronger position. And it **breaks
comparability**: every production metric in the 5.34M was computed on refolds, so
design-based ΔSASA would not sit on the same footing as any campaign-wide
distribution.

### 3.1d.1 The packing run — SPECIFIED, to execute when the CSV lands

**Decided 2026-08-20:** pack the three unpacked stacks, write the output to a
**separate archived folder**, and do **not** integrate it into the pipeline. It
is a standing fallback, kept off the critical path.

**Scope:** gated 40k ∩ {RFD3+LiMPNN, BG+LiMPNN, BD1+LiMPNN}. Cannot be launched
before the gated CSV exists, since the design list is the gate's output.
Confirm whether BD1 is present at all — it was disabled for csy4.

**Inode discipline — the reason for the separate archive.** Designs already live
inside `<uuid>.tar.zst`. Staging 40k loose CIFs would cost 40k inodes against a
project cap that has already caused one crisis. The pattern is **archive in →
pack → archive out**: one input archive per task, output written as chunked
`.tar.zst`, so inode cost is one file per output chunk rather than one per
design. Never stage loose.

**Carried-over gotchas** (all previously paid for):
- `archive_read.resolve` requires **`cwd=$RB`** or every extract fails *silently*
- cap `ARCH_READ_CACHE_MAX` or `/local/scratch` leaks unbounded
- **container bind trap** — if the packer runs under apptainer, per-task work
  dirs must be on bound Lustre, not node-local `$TMPDIR`; verify what the
  packer's container actually binds
- copy the invocation from `postprocess_halludesign.sh` (the known-good sibling
  already running this packer) rather than composing flags fresh

**Backbone invariance — verify, do not assume.** LigandMPNN's packer predicts
side-chain torsions on a *fixed* backbone, so N/CA/C/O should be untouched and
the packed files remain valid DockQ/RMSD ground truth. Confirm with the
**coordinate-md5 check** already used to prove `design12_seq0/seq4/seq9` share a
Cα hash. Two benign failure modes: precision round-trip through a different
format path (compare rounded coordinates, not raw bytes), and any spliced-in
relax step — which *would* move the backbone and must be kept out of this path.

⚠️ **Determinism.** If `number_of_packs_per_design > 1` or the packer samples,
design-based ΔSASA carries **its own draw noise**, and the zero-draw-noise
argument for switching away from §3.1c partly collapses. Pin to a single
deterministic pack, or measure the across-pack spread before relying on it.

**Recommendation: run it as a parallel VALIDATION, not a swap.** Pack the 40k,
compute design-based ΔSASA, correlate against the §3.1c refold-selected values.
High agreement validates the spec as written and ships unchanged; divergence is
worth knowing before ordering. The pipeline does not block on it — §3.1c is the
standing fallback, which is what makes this safe to attempt.

**Independent of this library:** sidechain-dependent analysis currently returns
*silently empty* on three of five stacks — "no clashes found" means "nothing was
tested". Packing at generation time closes that gap permanently.

### 3.2 Axes considered and rejected

From the seed reliability run (`seed_metrics_merged_v2.csv`, 30k structures),
ICC(1) as af3/b2/esm:

| candidate | ICC | verdict |
|---|---|---|
| ΔSASA total | .87 / .81 / .93 | included |
| ΔSASA base | .78 / .63 / .87 | included |
| BUNS (optH) | .70 / .77 / .79 | rejected — 77% correlated with ΔSASA total, redundant |
| contacts to base | .61 / .48 / .68 | rejected |
| **contacts**-to-base fraction | .57 / .46 / .61 | rejected — the original candidate axis; **not** the same as the included `dsasa_base_fraction`, which is built from ΔSASA (.78/.63/.87 and .87/.81/.93) |
| H-bonds total | .31 / .29 / .50 | rejected |
| π-stack count | .27 / .15 / .39 | rejected |
| H-bond base fraction | .18 / .16 / .23 | rejected — noise exceeds signal |

The rejected metrics are properties of **which refold draw you got**, not of the
design. The prior conclusion from that run stands: never use them as a
per-design attribute or a selection filter.

### 3.3 Metric inventory — what is measured on what

Reference for anyone extending this pipeline. Confirm the load-bearing rows
against the actual scripts before relying on them; two are flagged as unverified.

**Refold-based, SIDECHAIN-DEPENDENT** — the interface suite. All safe, because
refolds are full-atom for every stack:
ΔSASA (total / base / backbone / groove splits) · contacts (backbone/base,
pairwise AA↔nt) · H-bonds (`n_h_bonds_interface`, optH) · BUNS
(`delta_unsat_hbonds`) · π-stacking · red/unfavourable contacts · foot clash
(`foot_tp_min_dist`, explicitly "protein backbone AND sidechains") ·
tip-weighted contact score *(atom set unverified)*.

**Refold-based, BACKBONE-ONLY:** DockQ (`dockq_self`) · sc_RMSD
(`dockq_LRMSD`) · `ca_rmsd_vs_designed` · pairwise fold/pose/rna RMSDs · SSE
fractions. `lrmsd_rna_aligned.py` uses `PROT_BB=(CA,C,N,O)`; RMSD is Cα.
*Unverified:* DockQ's third term `fnat` is a contact fraction and its atom set
is not documented in the notes.

**Neither — model scores:** ipTM, pTM, pLDDT.

**Design-based:** only two, both backbone-derived and therefore safe — the
designed backbone as DockQ/RMSD ground truth, and `rg_backbone` (§3.1 D2).

**The consequence.** Every sidechain-dependent metric in the campaign is already
computed on refolds, which is why the missing design sidechains have never
affected a campaign number. **The exposure is entirely prospective** — it appears
the moment a sidechain-dependent metric moves onto design files. Two live
candidates: the §3.1d packing scenario, and the new atom-count foot clash, where
a *count* is far more atom-census-sensitive than a minimum distance. The second
feeds the gate (§1.1), so it is the one to pin down first — see §7 item 9.

## 4. Statistical power

Equal-N across stacks is the right allocation for "which stack should we run
next?", because it equalises the error bar per stack. Its limit, at 2,400/stack
(the five-stack case; four stacks gives 3,000 and shifts each row modestly in
your favour):

| baseline hit rate | detectable difference (80% power, α = .05) |
|---|---|
| 1% | 1% vs 2% — right at the edge |
| 0.5% | ≈ 2.5× |
| 0.1% | nothing — would need ≈ 23k/stack |

**If the expected hit rate is nearer 0.1%, the stack comparison this library is
designed to support cannot be made**, and the budget should be reconsidered.
Decide this consciously, before the order.

## 5. Pipeline

```
40k gated set
  │
  ① assemble          join axes; missing-value policy logged, never imputed
  ② normalize         rank → quantile per axis; zero-inflated → binary
  ③ diagnose          ties, correlation, effective dim → sets grid shape
  ④ stack allocation  equal base quota; shortfall redistributed by quality, capped
  ⑤ select 12k        tercile buckets, ipTM-ranked, per stack
  ⑥ compare           overlap first
  ⑦ verify            marginals, occupancy, stack table, per-uuid distribution
```

### 5.1 ① Assemble
Join the three bucketing axes plus `iptm_mean` onto the gated set. **Missing-value policy is logged as a
column, never silently imputed** — missingness correlates with difficulty (AF3
timeouts on larger complexes), so imputation pulls hard designs toward the mean
and biases coverage. Default: a design missing any axis is excluded from grid
placement but stays eligible for the §5.4 quality-redistribution pool.

### 5.2 ② Normalize

**Rank-normalise WITHIN STACK, not globally.** Stack is the outer stratifier
(§5.4) and each stack gets its own grid, so within-stack quantiles are the
natural choice — and they cancel any per-stack systematic offset for free
(sidechain provenance §3.1d, backend-win skew §3.1c, or a generator's
characteristic Rg range). A global rank transform would leave those offsets in
place, pushing whole stacks into uneven quantile bands and making their grid
cells unevenly reachable.

Rank → quantile, per axis. Rank rather than z-score because the metrics are
skewed; under z-scoring "spread" is dominated by whichever axis has the fattest
tail and budget goes to resolving outliers.

Axes with a tie fraction above 60% are **declared binary** rather than pretending
to be continuous — a quantile transform on an axis where most designs tie at one
value yields two effective levels, not four bins, and the grid's dimensionality
would be silently lower than advertised.

### 5.3 ③ Diagnose
Decision rules that **fire mechanically**, so nothing downstream is hand-chosen:

- tie fraction > 60% on an axis → that axis becomes binary (feeds ②)
- |Spearman ρ| > 0.7 between two axes → drop one, keep the more reliable
- effective dimensionality (participation ratio of the rank-correlation
  eigenvalues) < 3 → drop to a 3-axis grid

**Note on the ΔSASA pair (D3).** An earlier draft carried raw `dsasa_base`
alongside `dsasa_total` and needed a report-and-proceed exception, because the
two are nested and would have tripped the ρ > 0.7 rule. **That exception is
withdrawn**: `dsasa_base_fraction` divides size out, so the pair should now come
in well below the threshold and the ordinary auto-drop rule applies to it like
any other. If the measured ρ still exceeds 0.7, that is a genuine finding about
these designs and the rule should be allowed to fire.

The auto-drop rule still applies to any other pair.

### 5.4 ④ Stack allocation — the outer layer

Base quota = 12,000 / n_stacks (3,000 at four stacks; 2,400 at five).
**PROTECTED FLOOR = 1,200 designs per stack** (decided 2026-08-20).

**The floor is absolute, not a fraction of the base quota.** Its justification is
statistical power, which does not depend on how many stacks exist. At a 1% hit
rate, 2,400 designs detects a ≈2× difference in stack hit rate and 1,200 detects
≈2.5× — a modest loss, which is what makes 1,200 a *statistically meaningful*
floor rather than a token one. Below ≈500 a stack is present in the library but
mute: at 1% it yields ~5 expected hits and cannot support any comparison.

At five stacks the floor protects 6,000 of 12,000 — **half the library allocated
by fairness, half by quality.** Feasibility: n_stacks × 1,200 ≤ 12,000 holds at
four and five stacks. Assert it; if it ever fails, the floor must drop before
anything else does.

**Three-phase allocation:**

1. **Floor pass.** Each stack receives `min(1200, available_i)`. **Protected** —
   these slots are never redistributed away, whatever happens downstream.
2. **Equal pass.** Each stack tops up toward its base quota, as far as its
   remaining candidates allow.
3. **Quality pass.** Budget still unallocated — because some stacks were short —
   goes to the stacks with the best remaining `sc_RMSD` candidates, **capped at
   1.5× base quota** per stack (counting phases 1–2), subject to the
   cap-yields-to-count rule in §1.2 item 4.

**Why the floor makes the cap-relaxation safe.** §1.2 item 4 allows the 1.5× cap
to relax when exact-12,000 would otherwise be unreachable. Because floors are
allocated in phase 1, a phase-3 relaxation can only consume *unallocated* budget
— it can never eat into any stack's protected minimum. Without the floor, cap
relaxation and quality redistribution could compound and drift the library far
from equal-N.

**Assert `Σ available_i ≥ 12,000`.** If the gated pool is smaller than the target,
12,000 is unreachable and the run must fail loudly rather than silently
under-deliver.

**Rationale for the cap:** redistribution by quality preserves equal power only
while no stack is short. A 100-slot shortfall redistributes harmlessly; a
1,500-slot shortfall would otherwise convert a balanced design into a
quality-weighted one silently. The cap makes the departure a printed number
rather than an implementation detail.

⚠️ **What this does NOT guarantee.** The architecture guarantees *proportional
access to the budget*; it cannot guarantee *presence in the output*. If the gate
(§1.1) has already reduced a stack to near-zero — foot clash at `abs >= 3.4 both`
left RFD3 with **6 of 463** — no allocation scheme can recover it. The floor
protects a stack's claim on the budget, not its survival through the gate. The
only defence there is the **per-stack survival row** required on every candidate
gate threshold, which is a reporting discipline, not a mechanism.

### 5.5 ⑤ Selection — tercile buckets, ipTM-ranked

**PRIMARY METHOD (supervisor's design, adopted 2026-08-20).** Runs **within each
stack**, filling the quota that §5.4 assigned to that stack. §5.4 decides *how
much* each stack gets; this decides *which designs*.

**Per stack:**

1. For each of the three axes, split the stack's gated designs into **terciles**
   on that axis alone (rank-normalised within stack, §5.2). Three axes x three
   buckets = **9 buckets**. Within a single axis the terciles *partition* the
   stack, so overlap arises only *across* axes.
2. From each bucket take the top **N** by `iptm_mean`, descending.
3. **Union all 9 selections and deduplicate.**
4. Because all nine buckets rank by the same ipTM, high-ipTM designs are selected
   repeatedly and dedup bites hard — so **N is found by search**, not computed.
   Binary-search N until the deduplicated union is at or just under the stack
   quota, then **top up with the highest-ipTM unselected designs** to land on the
   quota exactly. Monotone, deterministic, and consistent with maximising ipTM.
5. **Tiebreak on `design_id`** throughout — with K exact (§1.2), ties must break
   reproducibly.

**Why terciles rather than halves.** Halves guarantee representation only in the
top and bottom half of each axis; top-ipTM within a half can cluster anywhere
inside it, so the middle band can vanish. Terciles protect both extremes **and**
the middle explicitly.

**Where a bucket is short.** If a bucket holds fewer than N candidates it
contributes all of them; the search on N absorbs the shortfall. Log any bucket
that cannot fill.

**Optional comparators, cheap to run and worth the number.** The two earlier
candidate methods can still be run and compared via §5.6:
- *quota grid* — 3 axes x 3 bins = 27 cells per stack, filled by ipTM, with
  largest-remainder apportionment and a deterministic tiebreak on cell index
- *deterministic kNN thinning* — fixed quota per density decile in the 3-axis
  rank space

Their value is the **overlap number**: if all three agree to about 85% or more,
the method choice is immaterial at K/N = 30% and the simplest, most explainable
one wins. That is a real possible outcome and is itself worth reporting.

### 5.5b ⑤′ Selection — Foldseek fold clusters (added 2026-08-26)

An alternative to §5.5 that replaces the three metric buckets with **one
structural partition**: sample every fold cluster, take every structural
singleton, rank within cluster by `iptm_mean`.

**Source.** `rna_binder/analyses/foldseek_map/` (project `~/Documents/binder_map`).
Foldseek all-vs-all on the **designed backbone**, `TM = min(qTM,tTM)`, Leiden at
fixed seed. Two completed runs:

| run | params | clusters | structural singletons |
|---|---|---|---|
| `results` (coarse) | TM >= 0.40, top-k | 245 / 258 / 279 (r05/r10/r20) | 182 |
| `results_fine` | TM >= 0.50, **mutual** kNN | 2,284 ... 2,510 (r05->r80) | 1,763 |

Singleton count is invariant to Leiden resolution: isolates have no edges above
the floor, and resolution cannot split what has no neighbours.

**Coverage is exact.** The map was built on `passing.csv` = **71,866 designs**,
which is the §1.1 base gate with **no foot term**. Verified 2026-08-26: 71,866 of
71,866 gated designs carry a cluster label, zero misses. Any foot variant is a
subset, so the labels transfer without recompute.

**Unit mismatch — clusters are over BACKBONES, not designs.** Foldseek ran on
36,907 unique backbones after hash dedup; those carry the 71,866 designs
(LiMPNN sequence variants collapse onto one backbone; `dedup/members.tsv` holds
the expansion). A per-cluster quota therefore buys *backbones*, and each backbone
can spend 1..N sequence variants against the 12,000. Taking more than one design
per backbone buys sequence diversity at **zero** fold diversity — that is a
legitimate trade but it must be a stated choice, not an accident of the quota.

**Resolution: use `results_fine`, not `results`.** The coarse run is unusable for
equal-per-cluster sampling — its largest cluster holds 20,069 of 71,866 designs
at r05, so an equal quota would hand a 20k cluster the same slots as a 3-design
one. Skew, measured as designs in the ten largest clusters:

| run / resolution | clusters | designs in top-10 clusters |
|---|---|---|
| coarse r10 | 258 | 47,139 (65.6%) |
| coarse r20 | 279 | 32,189 (44.8%) |
| fine r10 | 2,311 | 24,219 (33.7%) |
| fine r80 | 2,510 | 9,984 (13.9%) |

**D5 — the cluster axis must NEST INSIDE the stack allocation, never above it.**
Clusters are strongly stack-dominated (coarse cluster 0 is 3,366/3,374
RFD3+LiMPNN), and stacks reach very unequal numbers of clusters (fine r80: 820
and 774 for the top two, 352 for the bottom). Flat equal-per-cluster is therefore
a vote weighted by how many distinct folds a stack produced, and it **violates
the 1,200 floor**. Measured at K=12,000, fine r10:

| method | BD1 | BG+BIF | BG+LiMPNN | HD | RFD3 | clusters | singletons | backbones | ipTM mean |
|---|---|---|---|---|---|---|---|---|---|
| top-12k ipTM (control) | 518 | 2,010 | 4,527 | **436** | 4,509 | 624 (27.0%) | 325/1,763 | 7,493 | 0.7791 |
| flat cluster-equal | 1,668 | 1,750 | 3,766 | **884** | 3,932 | 2,311 (100%) | 1,763/1,763 | 7,510 | 0.7471 |
| **nested (§5.4 then clusters)** | 2,400 | 2,400 | 2,992 | 1,528 | 2,680 | 2,311 (100%) | 1,763/1,763 | **8,123** | 0.7350 |

Bold = below the protected floor. The nested form dominates the flat form on
every diversity measure *and* fixes the stack table, for 0.013 further ipTM.

**Procedure (nested).** §5.4 assigns each stack its quota. Then, per stack:
1. Take **every** design whose backbone sits in a singleton cluster, ranked by
   `iptm_mean` if the stack quota cannot hold them all.
2. Split the remaining quota equally over that stack's non-singleton clusters;
   fill each by `iptm_mean` descending.
3. Clusters smaller than their quota contribute all they have; the shortfall is
   redistributed by global `iptm_mean` within the stack.
4. Tiebreak on `design_id`.

**What this buys that the metric axes cannot see.** The control covers only
**624 of 2,311 clusters (27%)** and reaches 325 of 1,763 structural singletons —
so a purely quality-ranked library silently omits **73% of the fold space** the
campaign produced. Unique-backbone count barely moves (7,493 control vs 7,510
flat), which shows the gain is *not* backbone count: it is fold-space coverage.
No combination of Rg / dSASA / base-fraction detects this, because those are
continuous scalars and fold identity is not.

**Open sub-decisions.**
- Resolution: r10 (2,311 clusters, ~16 designs/cluster after singletons) or r80
  (2,510, ~11.8). Lower resolution = bigger, more heterogeneous clusters.
- Designs per backbone: cap at 1 (max fold diversity, 12,000 distinct backbones
  needed — feasible, 36,907 available) or leave uncapped (current sim gives
  8,123 backbones for 12,000 designs, i.e. ~1.5 sequences each).
- Whether the fold clusters **replace** the three metric axes, become a **fourth
  bucketing axis**, or act as a **coverage constraint** layered on §5.5.
- The **epitope** axis exists too (48-dim contact fingerprint, cosine, 1,111
  isolates on the fine run) and is orthogonal to fold. Not yet evaluated for
  selection.

### 5.5c ⑤″ RESULT — all methods measured, and the recommendation (2026-08-26)

Nine methods run on one footing over the base-gate pool (71,866, no foot term),
K=12,000, axes rank-normalised within stack, quality = `iptm_mean`.
Script: `$RB/analyses/library_selection/all_methods.py`.

| method | within-stack KS | clusters | singles | backbones | ipTM |
|---|---|---|---|---|---|
| control top-12k | 0.050 | 624/2311 | 325 | 7,493 | 0.7791 |
| tercile (nested) | 0.025 | 1147 | 732 | 8,326 | 0.7609 |
| grid 27 (nested) | 0.013 | 1135 | 723 | 8,375 | 0.7599 |
| kNN decile (nested) | 0.029 | 1145 | 732 | 8,327 | 0.7610 |
| cluster flat | 0.038 | 2311 | 1,763 | 7,510 | 0.7471 |
| cluster (nested) | 0.031 | 2311 | 1,763 | 8,123 | 0.7350 |
| **hyb clusterfloor+grid** | **0.011** | **2311** | **1,763** | **8,898** | **0.7511** |
| hyb clusterfloor+terc | 0.022 | 2311 | 1,763 | 8,855 | 0.7520 |
| hyb singlefloor+grid | 0.011 | 2168 | 1,763 | 8,810 | 0.7527 |

**RECOMMENDED: `clusterfloor+grid`.** Per stack, reserve the highest-ipTM design of
**every fold cluster that stack touches**, then apportion the remaining quota over the
27 joint cells by largest remainder, **crediting what the floor already took in each
cell** so the floor cannot distort the occupancy target.

**It strictly dominates `cluster (nested)`** — better metric uniformity, equal fold
coverage, more backbones, AND higher ipTM. No trade-off; drop the pure cluster method.
Against grid-27 it buys 2x the fold coverage and 523 backbones for 0.009 ipTM.

**The floor is cheap.** Stacks touch 334-754 clusters each, so it costs 2,651 of 12,000
slots (22%): BD1 30% of quota, RFD3 28%, HD 22%, BG+BIF 17%, BG+LiMPNN 14%.
Adding the floor slightly *improved* metric uniformity (0.011 vs grid-27's 0.013) — the
floor picks are spread through metric space and act as a well-distributed seed. 15%
difference; treat as marginal, but it is not a cost.

**Overlap:** hybrid family internally 93.1-97.5%; vs the metric family 87.0-91.3% (above
the §5.6 threshold); vs `cluster (nested)` only **79.4%** (below it) — so the hybrid is a
genuinely different library from pure clustering, and effectively grid-27 plus the
clusters grid-27 was missing. The three metric methods are ONE method (91.7-95.6% mutual):
that question is closed.
Interpretable alternative: `clusterfloor+terc`, 94.4% the same designs, costs 0.011 KS,
explains as "one design per fold cluster, then terciles".
`singlefloor+grid` is a false economy — 97.5% identical to the full floor, gives up 143
clusters, saves nothing.

### 5.5d Axis findings that change the axis list

**A1 — global marginals INVERT the ranking, and that is correct.** Bucketing methods look
*worse* than the control on global pool deciles (tercile KS 0.184 vs control 0.151 on rog)
because terciles are within-stack while the deciles are global, and §5.4 deliberately pulls
RFD3 from 52.8% of the pool to 22.3%. RFD3 is the high-Rg stack, so D9-D10 deplete.
**Stack fairness and global axis uniformity are in direct conflict and stack fairness wins
by construction.** Always report uniformity WITHIN STACK; the global view is misleading.

**A2 — `rog` is half stack.** Between-stack variance share **0.545** — Rg substantially
re-encodes stack identity, confirming D2's warning. Nesting it inside the stack layer is
what makes it usable.

**A3 — the base-fraction axis is INERT.** Between-stack variance share 0.016, correlation
with ipTM +0.034. Every method lands at KS 0.02-0.03 including the control, so bucketing on
it guarantees a coverage that plain ipTM ranking already delivers. It is currently dead
weight. ⚠️ This verdict is on the **contacts** proxy (see A5); the real ΔSASA base fraction
may behave differently.

**A4 — design vs refold Rg is a distinction without a difference: r = 0.997.** D2's
reasoning stands in principle, but after rank-normalisation the two are interchangeable.
Design-backbone Rg is now computed for all 36,907 backbones from the Foldseek CA traces
(`$RB/analyses/library_selection/rog_design.py`, median 15.19 Å vs refolds' 13.42/12.98) and
is what the pipeline uses. No recompute was needed: the traces already existed.

**A5 — ⚠️ THE BASE-FRACTION AXIS IS THE ONE §3.2 REJECTED.** `combined_clean_v2v9` has NO
RNA-side ΔSASA split — only `dsasa_total/polar/hphobic`, all protein-side. The store's
`fbase_*` was sourced from **`frac_base_contacts`**, i.e. the *contacts* base fraction,
ICC .57/.46/.61, which §3.2 rejects in favour of ΔSASA base at .78/.63/.87. To get the
specified axis, extend `contacts_sidecar --sasa` to the gated set (~1.2 CPU-h for 42k x 2).
Until then the third axis is a stand-in, and A3 says it is contributing nothing anyway.

### 5.5e ⚠️ STATUS 2026-08-26 (late) — THE POOL AND THE AXIS LIST BOTH MOVED

Everything in §5.5b–d was computed on the **pre-lrfix** pool of 71,866. Two findings
landed afterwards that invalidate those numbers (though not the method ranking).

**(1) The LRMSD register fix shrinks the gate.** `combined_clean_v2v9_fullmetrics_bf_lrfix.csv`
(6,821,957 rows exact) changed `dockq_self_*` AND `dockq_LRMSD_*` on 100% of rows — both are
gate terms. Base gate **71,866 → 65,750** (gained 3,901, lost 10,017).

| stack | pre | post | delta | new quota |
|---|---|---|---|---|
| BD1+LiMPNN | 3,344 | 2,381 | −963 | 2,381 (0.99×) |
| BG+BIF | 9,881 | 8,394 | −1,487 | 2,400 (1.00×) |
| BG+LiMPNN | 19,162 | 16,460 | −2,702 | 3,063 (1.28×) |
| **HD+LiMPNN** | 1,528 | **564** | **−63%** | **564 (0.23×) — BELOW FLOOR** |
| RFD3+LiMPNN | 37,951 | 37,951 | **0** | **3,592 = 1.5× CAP BINDING** |

RFD3 unchanged **to the row** confirms the fix: it was the only in-register stack.
The mechanism is the opposite of what was predicted — `LRMSD<=3 either` fell 219,484→171,666
while mean LRMSD *improved*, i.e. the buggy `_match` computed RMSD over a **truncated atom
set** and flattered the fit. **The bug inflated 4 of 5 stacks; the fix removes false passes.**

**HD at 564 is below the 1,200 protected floor and at the ~500 mute line (§5.4) — the HD arm
is effectively dead.** ⚠️ UNRESOLVED: BD1 and HD are the **18-nt** stacks and lost the largest
fractions (−63%, −29%) against BG's −14/−15%. Whether the fix handles the length mismatch is
NOT verified. **Resolve before accepting HD's number.** Options if real: accept a mute HD, or
drop to four stacks at 3,000 each (statistically cleaner than pretending 564 supports a
comparison).

Post-fix the pool spans only **1,772 of 2,311 clusters**, and 3,901 designs have no cluster
label (the Foldseek map was built from the pre-fix `passing.csv`). The store is also pre-fix,
so store and CSV now disagree on DockQ/LRMSD everywhere — patch the store additively, or read
the CSV, but never both.

**(2) ⚠️ CROSS-BACKEND AGREEMENT INVALIDATES THE ΔSASA AXES.** AF3 vs ESM, same metric, same
design, n≈100k — a test the spec never applied to any axis:

| metric | global r | within-stack r |
|---|---|---|
| `dsasa_total` | 0.590 | 0.36–0.60 |
| `frac_dsasa_base` | 0.293 | 0.17–0.28 |
| **`frac_base_contacts`** | **0.159** | **0.08–0.17** |

`frac_base_contacts` — the axis every method in §5.5c actually used — is essentially
uncorrelated across backends. It was a near-random third coordinate, which **re-explains the
"inert axis" result of §5.5d A3**: not that binding mode does not vary, but that the column
cannot measure it as a property of the design. Within-stack r is below global r for *every*
metric, so the nesting gets the weaker number.

**`rg_backbone` and the fold clusters are the only axes immune** — both computed on the
designed backbone, with no backend to disagree with. D2's justification is now the only one
that survives contact with the data.

⚠️ Do not conflate this with §3.2's ICC table (.78/.63/.87): that was seed-to-seed **within**
a backend. Cross-**model** is a harsher test. And do not overstate it — low cross-backend r
may partly reflect genuine pose differences, with each metric faithfully describing its own
structure. It still fails as a per-design attribute.

**ACTIONS.**
- **§3.1c is now wrong.** Picking ΔSASA from the single `sc_RMSD`-winning backend takes ONE
  noisy draw. Use the **mean of af3+esm**: Spearman-Brown lifts `dsasa_total` 0.59 → ~0.74 free.
- **Proposed axis list: bucket only on design-derived quantities** — `rg_backbone` + fold
  cluster — and demote ΔSASA to a reported descriptor. A real departure from the agreed spec;
  needs sign-off.
- `frac_dsasa_base` overlap with the post-fix pool is **852/65,750 = 1.3%**; filling it costs
  **43 CPU-h** (52 tasks × ~6 min), NOT the "1.2 CPU-h" in §7 — that was a s→h unit slip.
  Worth it for analysis value, not as bucketing axis 3.
- ⚠️ §5.5d A5 cited r=0.54 between `frac_dsasa_base` and `frac_base_contacts`. Measured
  af3-vs-af3 on n=99,246 it is **r = 0.787**; the 0.54 was almost certainly the bare
  **consensus** columns, which are row-wise minima.

### 5.6 ⑥ Compare
The comparison set is now four methods plus the control: §5.5 terciles, §5.5b
fold clusters (nested), and the two optional comparators (quota grid, kNN
thinning), against the top-12k ipTM baseline.

**Compute pairwise overlap / 12,000 first.** Above about 85% agreement the
method choice is immaterial — ship the primary method (§5.5) and stop. Given K/N = 30% and a
stack-first structure that fixes much of the allocation outright, high overlap is
the expected outcome.

If they diverge, compare on:
- per-axis marginal uniformity (KS distance to uniform on quantiles)
- joint occupancy scored on an **independent grid** — scoring a method on its own
  bucketing hands it the win by construction
- quality retained vs. the **top-12k baseline** — the 12,000 highest
  `iptm_mean` designs, ignoring stacks and axes
  entirely. This is the control both methods must be measured against: it prices
  what the diversity cost, and given K/N = 30% it may well land within noise of
  either method, which is itself a finding worth reporting
- distinct sequence clusters and epitope fingerprints covered

### 5.7 ⑦ Verify
Per-axis marginals, joint occupancy, final stack table, per-uuid distribution.

## 6. Decisions taken, with rationale

**D5 — constraints inside the selection loop, not after.**
Stack floors and quotas are feasibility checks applied *during* fill. Selecting
for coverage and then patching balance by swapping designs in destroys the
coverage just paid for.

**D6 — the objective is to MAXIMISE ipTM; bucketing is a representation floor.**
ipTM ranks every bucket. Three superseded proposals — DockQ, `sc_RMSD`, and
epitope spread — are recorded in §3.1b so none is re-introduced.

**D7 — no per-uuid or per-sequence-cluster caps.**
Dropped by user decision. The risk is real: a single RFD3 backbone can contribute
many designs sharing an *identical* pose — `design12_seq0/seq4/seq9` off one
backbone have the same Cα coordinate md5, differing only in sequence. In a
one-shot library those are correlated bets on one fold.

The countervailing reason is also real, and is why the cap is not reinstated:
multiple sequences on a fixed backbone is the **only** way this library can
reveal whether LiMPNN's sequence choice matters independently of the backbone. A
cap would destroy that question.

Mitigation, at zero cost to the decision: **report the distribution without
constraining it** — designs per uuid and per sequence cluster appear in the §7
diagnostics, so what was bought is known even though it was not restricted.

## 7. Open items — pending the gated CSV

### 7.0 STATUS 2026-08-24 — ALL DATA BLOCKERS CLEARED
`combined_clean_v2v9_fullmetrics.csv` (12 GB, 178 cols, 6,821,957 rows = exact store match)
has rog + dsasa + frac_base at ~100% for ALL campaigns v2-v9, per-stack verified (HD and BD1
went 42.8%/51.7% -> 100%). Filter store augmented via `add_rog_fbase_cols.py`;
`placeholder_cols` is now empty. Inodes 7,138,508 -> 6,302,800.
**Remaining: re-export the gated set with the new columns, then implement the §5.5 selection.**

1. **Per-stack gated counts**, and confirmation of which stacks are in the 40k
   (assumed four or five; BD1+LiMPNN was disabled). Sets the base quota and
   determines whether the §5.4 cap ever fires.
2. **The foot-clash gate threshold** — upstream and unsettled. It determines both
   N and the stack composition of the pool. This spec's *structure* is
   independent of it; its *numbers* are not.
3. **Verify `iptm_af3` / `iptm_esm` column names** in the v2–v9 table and confirm
   whether a mean column already exists or must be derived.
4. **Measured tie fractions and the rank-correlation matrix** — these set the
   grid shape via §5.3 and are not guessable in advance.
5. **Which SASA family** (Biotite or Rosetta) the ΔSASA columns come from.
6. **`iptm_mean` spread inside the gated set.** ipTM is *not* gated (§1.1) so it
   should retain full range — confirm that, and check the **per-stack ipTM
   distributions**, since ipTM-ranking fills every bucket and a large per-stack
   offset would skew which designs each stack contributes. A spread of ≳0.15
   gives ~10 noise units and a meaningful tiebreak; much tighter and the
   ranking is largely noise and within-cell choice is effectively arbitrary — an
   acceptable outcome given the gate has already removed the bad designs, but one
   that must be known rather than discovered.
7. **Verify `dockq_LRMSD_af3` / `_esm` coverage** on the gated 40k — §3.1c
   backend selection needs both present; define the fallback where one is
   missing (proposal: fall back to the backend that has a value; log the count).
   Also record whether the gate applied RMSD as `[AND]` or `[OR]` across
   backends, and how many designs had their backend choice forced.
8. **ICC of `dsasa_base_fraction`** — computable from
   `seed_metrics_merged_v2.csv` with no new compute (§3.1 D3, §7.1).
9. **Launch the §3.1d.1 packing run** once the design list exists — separate
   archived output, not integrated.
10. ~~Where the atom-count foot clash is computed~~ — **CLOSED 2026-08-20.**
    The `footn_*` columns are computed on **refolds** (`_af3` / `_esm`), so the
    sidechain-census bias does not apply. Eleven thresholds (0.5–10 Å), 0% NaN
    across all 6,821,957 store rows, monotone, saturating at 129 heavy atoms
    (the six appended nucleotides). This is the graded profile form, which
    avoids the threshold sensitivity that holds H-bond ICC at 0.31.
11. **BOTH MISSING AXES ARE RECOVERABLE WITHOUT ANY RECOMPUTE** (settled
    2026-08-20 on Euler). Root cause chain, verified end to end:
    - The per-backend fullmetrics OUTPUT chunks (`chunk_<be>_<NN>.csv`, 141 cols)
      contain **`radius_of_gyration_ca`** (col 63) and **`frac_base_contacts`**
      (col 16) in **all three** chunk dirs:
      `fullmetrics_chunks` 1566/1566 · `fullmetrics_chunks_v8` 573/573 ·
      `fullmetrics_chunks_v9` **477/592 (81%)**.
    - `consensus_integrate.py` merges **only columns declared in a family**, and
      `FAMILIES` = interface, hbonds, red, dockq, pistack. **There is no
      `structural` family**, so all **20** structural columns were left behind:
      `radius_of_gyration_ca`, `helix_fraction`, `sse_helix_fraction`,
      `sse_strand_fraction`, `rna_contact_nt_1..16`. Verified: the merged master
      has **0** of them, a chunk has **20**. `frac_base_contacts` survived only
      because it belongs to `interface`.
    - The store builder read the merged master, found no rog, and wrote a
      placeholder — which `app.py` handles gracefully, so it never errored.
    - Timing rules out "the merge predates the metric": chunks are 2026-07-28/29,
      the merged master is 2026-08-05.

    **FIX:** add a `structural` family to `consensus_integrate.py`, re-run against
    the existing chunks, rebuild the store mapping `rog_af3 →
    radius_of_gyration_ca_af3`. **No metric recompute.**

    ⚠️ **Design decision this forces:** the five existing families use worst-case
    consensus (min/max). That is meaningless for Rg and SSE fractions — there is
    no "worse" radius of gyration. What is needed is the **per-backend columns
    merged with no consensus column**, which may require a merge-only path rather
    than just a new `FAMILIES` entry.

    **Residual gap:** v9 is at 477/592 chunks, so ≈19% of v9 designs stay
    uncovered — ≈1,830 of the 9,632 gated v9 rows, or **≈4.4% of the 41,890**.
    Handle under the §5.1 missing-value policy, or finish the remaining 115 chunks.

    **Side benefit:** restores `rna_contact_nt_1..16` (the epitope fingerprint) to
    the master. Rejected for selection (§3.1b) but available as a diagnostic.

### 7.2 Producing `combined_clean_v2v9_fullmetrics.csv` — the unblocking work

Verified on Euler 2026-08-20/21. **Nothing needs recomputing except the v9 tail.**

**Precedent:** there is **no v2v8 base table** — v8 was appended directly to the
v2v7 fullmetrics table (`combined_clean_v2v7_fullmetrics.csv` 6.7G →
`combined_clean_v2v8_fullmetrics.csv` 9.1G). The last real base is
`combined_clean_v2v7.csv` (3.6G). So appending v9 to the v2v8 fullmetrics table
follows the established pattern. `merge_masters.py` / `gapfill_merge.py` exist —
**find how v8 was appended and copy it** rather than writing a new path.

**Step 0 — assemble the v2v9 table.** Append v9's rows from `v9_store_table.csv`
(21 cols: identity + af3/esm paths + iptm; b2-stripped because it was built for
the b2-free store) onto `combined_clean_v2v8_fullmetrics.csv`, which already
carries the five families for v2–v8. v9 rows go NaN for b2 path/iptm — cosmetic,
because `consensus_integrate.py` merges chunks **by `design_id`** and will still
pick up b2 values from `chunk_b2_*.csv` (477 exist for v9).

**Step 1 — finish the v9 tail — DONE + VERIFIED 2026-08-21.** af3/b2/esm each 592/592 files, 1,478,433/1,478,433 rows, deficit 0, metric columns confirmed populated. Jobs 11413058/59/71 (indices 477-591) plus 11416145/52/57 (index 476, which the original run left truncated). Original note kept for the record: Missing chunk indices are a
**contiguous tail 477–591**, identical across all three backends — the array
never finished, these are not scattered failures. 115 × 3 = **345 tasks**.
At the measured 4.8 CPU-s/design/backend: 115 × 2,501 × 3 ≈ 863k units ≈
**1,150 CPU-h** (~35 min wall at 2,000 CPUs). Skipping leaves v9 at 81% and
≈4.4% of the gated 41,890 NaN on the new columns.

**Step 2 — add a `structural` family** to `consensus_integrate.py`.
⚠️ Worst-case min/max is meaningless for Rg and SSE fractions. What is needed is
the per-backend columns merged with **no consensus column** — likely a merge-only
path, not just a `FAMILIES` entry.

**Step 3 — two targeted runs, not a full re-merge:**
- `--chunkdir fullmetrics_chunks_v9 --families interface,hbonds,red,dockq,pistack`
  → fills the five families for **v9 rows only**
- `--chunkdir <each of the three> --families structural`
  → fills the 20 structural columns for **all** rows

**Step 4 — verify** row count against the store's 6,821,957 and column
consistency between the v2–v8 and v9 rows.

**Then rebuild the filter store**, mapping `rog_af3 → radius_of_gyration_ca_af3`.

⚠️ **DISK.** `consensus_integrate.py` backs up the master on **every** run. The
table lands near 11 GB and step 3 is four runs ≈ **44 GB of backups**. Check quota
first or disable the backup after the first run. Use `--dry-run` (supported)
before committing to an 11 GB write.

### 7.1 Accepted risk
The reliability of `dsasa_base_fraction` is unmeasured (§3.1 D3). Its two
components are individually reliable, but a ratio can be worse than either part.
Resolvable from `seed_metrics_merged_v2.csv` at zero compute cost — §7 item 8.
If its ICC lands below ≈0.6 it should be replaced, most likely by raw
`dsasa_base` with `dsasa_total` dropped to avoid the nesting problem.

## 8. Deliverables

Interpretability is a stated requirement, so these are part of the spec, not a
byproduct:

1. **Stack allocation table** — base quota → final quota, with cap events flagged
2. **Per-stack occupancy heatmap**
3. **Redistribution ledger** — slots moved, from which cells, to which
4. **Target vs. achieved marginal**, per axis, per stack
5. **Quality-cost table** — what the tercile floors cost against the top-12k ipTM baseline
6. **Overlap matrix** — primary method vs each optional comparator vs the
   top-12k ipTM baseline
7. **Per-uuid and per-sequence-cluster distribution** (reported, not constrained)
8. **Fold-cluster coverage table** — clusters covered / total and singletons
   hit / total, per method and for the control (§5.5b)
9. **ΔSASA backend-win table** — which refolder was selected by `sc_RMSD`, per
   stack, so any stack-encoding in the ΔSASA axes is visible (§3.1c)

**Every design in the final 12k must be traceable to a one-line reason:
stack S, bucket(s) B that selected it, ipTM rank R within each — or "top-up" for
designs added by the §5.5 step-4 fill.**
