# What the 12k library actually selected for — sequence signature vs the full campaign

Run 2026-10-13. Corpus: `csy4_master_foot.csv` (6,992,515 designs, v2–v9) vs the delivered
panel `csy4_binder_library_12k_v6_k30.csv` (12,000).
Code + results: `$RB/analyses/library_selection/seqsig/` (Euler).
Figures: `seqsig/sequence_signature.html`.

---

## The headline

**91% of the library's amino-acid signature was imposed by the quality gate, not by the
selection algorithm.** The selector — stack quotas, fold/epitope cluster quotas, ipTM
ranking — is very nearly sequence-neutral: it moves composition by a mean of **0.035 log2**
and never by more than **0.168**, against gate shifts of up to **1.44**.

And what the gate selects for is **α-helix**. Not as a stated criterion — nothing in
`sc_RMSD ≤ 3 | DockQ ≥ 0.7 | seq_run ≤ 5 | footnb ≤ 5 | ipTM ≥ 0.6` mentions a residue —
but the structural self-consistency filter behaves as an implicit helix filter, and the
sequence composition records it unmistakably.

If you want one number for the thesis: **GGG is 10× depleted (−3.37 log2) and KLE is 3×
enriched (+1.61).** Turns out, loops out, amphipathic helix in.

---

## 1 · The funnel, reproduced exactly

The gate was rebuilt from the npy store and **lands on 55,706 and all five per-stack counts
to the row** — this is the gate that shipped, not `passing_v1` or `passing_v4`.

| cumulative criterion | remaining | % of total | kept from prev |
|---|---|---|---|
| all designs | 6,992,515 | 100% | — |
| sc_RMSD ≤ 3 Å, either backend | 171,843 | 2.46% | 2.5% |
| DockQ ≥ 0.7, either | 133,863 | 1.91% | 77.9% |
| longest single-AA run ≤ 5 | 93,022 | 1.33% | 69.5% |
| footnb @1.0 Å ≤ 5 atoms, both | 64,185 | 0.92% | 69.0% |
| max ipTM ≥ 0.6 | 55,706 | 0.80% | 86.8% |
| → selection | 12,000 | 0.17% | 21.5% |

Gate survival varies **23× across stacks**, which is why the selector had to rebalance:

| stack | campaign | gate | survival | delivered |
|---|---|---|---|---|
| RFD3+LiMPNN | 1,524,796 | 27,881 | 1.83% | 3,132 |
| BG+LiMPNN | 1,519,958 | 14,752 | 0.97% | 3,266 |
| BG+BIF | 1,480,140 | 9,042 | 0.61% | 2,400 |
| BD1+LiMPNN | 1,444,614 | 3,229 | 0.22% | 2,400 |
| HD+LiMPNN | 1,023,007 | 802 | **0.078%** | 802 |

Note the campaign itself is near-balanced (1.0–1.5M per stack). All the imbalance is created
by the gate and then partly undone by the selector.

---

## 2 · The gate did it, not the selector

`log2(library/campaign)` splits **exactly** into `log2(gate/campaign) + log2(library/gate)`,
so this decomposition is arithmetic, not a model.

| stack | mean \|gate\| | mean \|selector\| | gate share | max \|selector\| |
|---|---|---|---|---|
| BD1+LiMPNN | 0.484 | 0.036 | 93% | 0.145 |
| BG+BIF | 0.408 | 0.028 | 94% | 0.085 |
| BG+LiMPNN | 0.432 | 0.048 | 90% | 0.162 |
| HD+LiMPNN | 0.212 | **0.000** | **100%** | **0.000** |
| RFD3+LiMPNN | 0.224 | 0.062 | 78% | 0.168 |
| **all** | **0.352** | **0.035** | **91%** | |

Two internal checks fall out of this table and both pass:

- **HD's selector column is exactly 0.000 across all 20 residues.** HD received its entire
  gated pool (802 of 802) — no selection happened inside it, and the decomposition correctly
  reports zero. A method that produced any other number here would be wrong.
- **RFD3 has the largest selector term (78% gate share).** It was the most heavily
  down-sampled stack. The one place the selector *does* touch sequence is exactly where it
  did the most sub-setting, which is the expected direction.

**Implication.** The selection algorithm you spent the design effort on is not what shaped
the library's chemistry. If you want different sequences next campaign, change the gate —
tuning the selector will not move this.

---

## 3 · The gate is an implicit α-helix filter

Chou–Fasman helix formers (A, E, L, M, Q, K, R) against breakers (G, P, T, V, S), mean
log2(library/campaign):

| stack | helix formers | breakers | gap |
|---|---|---|---|
| BD1+LiMPNN | +0.542 | −0.807 | **+1.349** |
| BG+BIF | +0.356 | −0.392 | +0.748 |
| BG+LiMPNN | +0.284 | −0.386 | +0.670 |
| HD+LiMPNN | +0.143 | −0.113 | +0.257 |
| RFD3+LiMPNN | +0.164 | +0.132 | **+0.033** |

The motif scan says the same thing far more bluntly (3-mers, 2,685 tested, **1,736 at
FDR < 1%**):

| enriched | log2 | | depleted | log2 |
|---|---|---|---|---|
| MSL | +1.78 | | **GGG** | **−3.37** |
| INS | +1.65 | | GGP | −2.75 |
| KLE | +1.61 | | PPP | −2.55 |
| KQL | +1.57 | | GPG | −2.53 |
| LEK | +1.53 | | GTG | −2.53 |
| ELE | +1.48 | | PGG | −2.52 |
| KKL | +1.40 | | GNG | −2.40 |
| LKK | +1.39 | | GSG | −2.38 |

KLE / LEK / ELE / KKL / LKK are textbook charged amphipathic helix. GGG / PPP / GPG / PGG
are loops, turns and disorder. Supporting physicochemical medians: fraction K+R rises in
4/5 stacks, GRAVY falls (more hydrophilic) in 4/5, Gly and Pro fall.

**Interpretation, flagged as such:** the gate's three structural terms (sc_RMSD, DockQ,
ipTM) are all structure-prediction self-consistency measures, and AF3/ESMFold assign their
highest confidence to idealised helical folds. The composition shift is the observation;
"the predictors prefer helices" is the proposed mechanism, not something this analysis
proves.

---

## 4 · Two places the simple story breaks

I nearly reported a universal aromatic depletion. Within-stack numbers say otherwise:

| stack | aromatic (F+W+Y) campaign → library | fold |
|---|---|---|
| BG+LiMPNN | 0.0373 → 0.0188 | **0.50×** |
| BG+BIF | 0.0417 → 0.0263 | **0.63×** |
| BD1+LiMPNN | 0.0608 → 0.0613 | 1.01× |
| HD+LiMPNN | 0.0327 → 0.0379 | 1.16× |
| RFD3+LiMPNN | 0.0125 → 0.0151 | 1.21× |

**The aromatic loss is BoltzGen-specific, not a property of the funnel.** Worth caring about
because aromatic stacking is a major natural protein–RNA binding mode, and in the two
BoltzGen arms the gate halves it — but it would have been wrong to claim the pipeline
globally selects against it.

**RFD3 inverts.** It shows no helix enrichment (+0.033) and its biggest move is alanine
*down* 0.62× (0.354 → 0.219). The reason is where it started: **RFD3's campaign sequences
are 35.4% alanine.** There is no helix propensity left to select for, so the gate instead
removes the most degenerate, poly-Ala-like sequences. That 35% is itself anomalous and
lines up with the known RFD3 sequence-collapse issue (`rna-binder-rfd3-limpnn-seq-collapse`,
10× duplicate sequences). Read positively: **the gate partially rescues RFD3 from its own
sequence degeneracy.**

---

## 5 · Both controls pass

- **Control A — stack reweighting.** The library's stack mix differs hugely from the
  campaign's (RFD3 0.52×, BD1 3.45×, HD 4.64×), so a pooled comparison could be measuring
  stack choice rather than sequence. Re-running against a campaign reweighted to the
  library's stack shares moves **17 of 19 residues by less than 0.08 log2** (largest: Pro
  0.230, Trp 0.145). The signature is genuine within-stack chemistry, not a composition
  artifact.
- **Control B — backbone redundancy.** The 12,000 designs sit on only **9,759 unique
  backbones**, 1,307 carrying 2–10 LiMPNN variants, so a few over-represented backbones
  could be talking loudly. Collapsing to one random sequence per backbone (200 repeats)
  moves the largest of 100 stack×residue cells by **0.079 log2** — against a signature of up
  to 1.4. Not redundancy.

---

## 6 · Two axes are degenerate — don't analyse them again

- **Length is constant.** All 6,992,515 designs are 80 residues except **620 at 79
  (0.009%)**. No length confound exists, and length-matched null sets would be pure
  ceremony.
- **Cysteine is absent everywhere.** 0 designs contain a cysteine, in all 5 stacks, across
  all 6,992,515. Every C-column log2 is therefore 0/0 = NaN by construction.

⚠️ **This contradicts a stored memory.** `rna-binder-cysteine-by-stack` records RFD3-LiMPNN
at 17.4% and BG+BIF at 21.8% cysteine. That measurement was scoped to the **midrun_v2
aggregation (909k designs, 2026-06-25)**; v2 *is* in this master (507,495 designs) and still
reads 0%. So the v2 contamination is not present in the analysis corpus — it was cleaned or
regenerated upstream (consistent with `rna-binder-bgbif-cysregen`). The 17%/22% figures
should not be quoted against the current master. **Not chased further — flagging, not
resolving.**

---

## 7 · What this does not say

- **Nothing here is about binding.** Every number is composition in a *selected set*. "The
  funnel enriches helix" is not "helical binders work better" — the wet-lab assay is what
  answers that, and this analysis is the thing you compare those results against.
- **The gate/selector split is causal only in bookkeeping terms.** It attributes the shift
  to a stage, not to a mechanism within that stage. Which of the five gate terms carries the
  helix preference is a separate, answerable question (re-run the decomposition term by
  term) and was not in scope.
- **Control B is asymmetric.** The library is backbone-collapsed; the campaign is not, since
  the master carries no backbone/coord-hash column. Given the effect is 0.079 log2, this
  does not change anything, but it is not a symmetric test.
- **HD is underpowered and BG+BIF entered on a disputed reference** — both already on
  record. HD's numbers here rest on 802 designs, and its selector column is zero by
  construction rather than by measurement.

---

## 8 · Reproduce

```
$RB/analyses/library_selection/seqsig/
  extract_seqs.sbatch   # awk pass over the 14.5 GB master -> seqs_all.tsv (49 s, bad_nf=0)
  gate_ids.py           # rebuild the shipped gate from the npy store; asserts 55,706
  featurize.py          # 20 AA counts + length per design; k-mers per (stage x stack)
  analyze.py            # composition, physchem, Control A, k-mer FDR
  control_b.py          # backbone-collapse control
  summary.py            # the headline tables above
  plots.py              # self-contained HTML panel
  results/*.tsv         # every number in this document
```

Everything derives from **20 amino-acid counts + length per design** — net charge, K+R
fraction, GRAVY, aromatic fraction and pI are all functions of those, so nothing else is
stored. Run order is the order listed. Total cost: one `awk` pass plus a few minutes of
single-core Python; no GPU.
