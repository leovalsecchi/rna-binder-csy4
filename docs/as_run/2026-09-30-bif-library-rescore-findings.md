# BG+BIF library rescore — findings

**Date:** 2026-09-30 (written up 2026-10-01)
**Question:** of the 2,400 BG+BIF designs in the delivered 12k library, how many fail
the campaign's own gate when scored against the **true BoltzGen backbone** instead of
against BIF's own re-prediction?
**Artefacts:** `analyses/bif_library_rescore/` — `REPORT.md`, `rescore_2400.csv`,
`PROGRESS.md`. Spec + plan: `docs/2026-09-30-bif-library-rescore-{design,plan}.md`.

---

## 1. The defect

`code/scripts/stages/boltzgen_multif.sh` runs **one BoltzGen job that feeds both inverse
folders**, so BG+BIF and BG+LiMPNN share a `uuid8` and the same `spec_<D>.cif` backbones.
The two stacks then diverge in what is stored as "the design":

| stack | stored design | what it actually is |
|---|---|---|
| BG+LiMPNN | `spec_<D>` + MPNN sequence | the **true backbone** — LigandMPNN moves no atoms |
| BG+BIF | harvested from `bg_output/.../refold_cif/` | a **Boltz prediction** of the designed sequence |

`--steps inverse_folding folding` runs the *folding* step, and `boltzgen_patch/writer.py:92`
writes that file from `pred_dict["coords"]` — best of N samples by `0.8·ipTM + 0.2·pTM`,
with pLDDT in the B-factor column.

So for BG+BIF alone, `dockq_self`, `dockq_LRMSD` (= `sc_RMSD`) and `ca_rmsd_vs_designed`
were measured against **a prediction of the same sequence**, not against the generated
backbone. Both are gate terms. The other four stacks are unaffected.

## 2. Result

Scored 2,397 of 2,400. Accounting closes: 2,397 scored + 1 unresolvable + 2 failed = 2,400.

| | campaign reference | true backbone |
|---|---|---|
| passed the gate | 2,397 *(definitional — the library was selected by this gate)* | **434** |
| sc_RMSD median | 2.744 Å | **8.663 Å** |
| DockQ median | 0.709 | **0.406** |

### **1,963 / 2,397 = 81.9% flip PASS → FAIL.**

Uniform across campaigns (v2 83% · v4 86% · v5 78% · v6 80% · v7 83% · v8 80% · v9 83%);
χ² = 6.27 on 6 df, p ≈ 0.39, i.e. indistinguishable from one common rate — consistent with
a single shared mechanism (the same folding step in every campaign).

The earlier 20-pair probe estimated 70%; it was directionally right but understated it.

## 3. ipTM does not rescue the library

`iptm_mean` = mean(AF3, ESMFold2). BG+BIF over the 2,400: **mean 0.6987**, median 0.7081,
range 0.390–0.867 (AF3 0.7268 / ESM 0.6705). Highest of the five stacks; 12k-wide mean 0.6849.

| | n | mean ipTM |
|---|---|---|
| flipped PASS → FAIL | 1,963 | 0.6938 |
| still passing | 434 | 0.7210 |
| difference | | **+0.027** |

That gap is ~1 seed-noise unit (mean-of-2 noise ≈ 0.025), i.e. essentially no signal.
Survival by ipTM decile runs 8.8% (lowest) to 25.6% (highest), non-monotonically —
**even the top ipTM decile flips at 74%.** Selecting harder on ipTM would not have avoided
this, which is unsurprising: ipTM is a model confidence score and never looked at the
designed backbone at all.

## 4. Why the numbers can be trusted

Every link was measured, not assumed:

- **Regression:** the rebuilt pipeline reproduces the earlier validated 20-pair study at
  max |Δ| = **0.000000** (40/40 values).
- **Baseline:** re-deriving the *published* campaign values from scratch on 266 v6
  comparisons agrees to max |Δ| 0.0005 Å / 0.0001 — so the library CSV is a sound baseline,
  and the extracted refolds are demonstrably the campaign's own.
- **Sequence gate:** 0 mismatches between the recorded sequence and the folded structure
  across all 2,399 (this is the check that would catch the v2 BG+LiMPNN/BG+BIF mislabelling
  precedent).
- **Pairing null test** — the one failure mode no gate covered: is each design paired with
  the *right* backbone? Own pairing gives 0.834 Å median binder-Cα fold RMSD (128/133 < 2 Å);
  shuffled pairing gives 12.828 Å (0/133). Total separation.
- **Independent recompute** of the headline from the raw table, without the reporting script:
  1963/2397, 81.894%, 0 gate disagreements over 2,397 rows.
- Campaign-side values are **byte-identical** to `csy4_binder_library_12k_v6_k30.csv`.

## 5. Caveats

- **Not a stack-wide rate.** These 2,400 were selected *because* they scored well under the
  wrong reference. Conditioning on that maximises the expected regression. The BG+BIF-wide
  figure needs the full 1.48M.
- **fnat is not like-for-like.** The corrected reference is protein-backbone-only; the
  campaign one carries sidechains, and only DockQ's fnat term reads sidechains. Stripping
  sidechains *flatters* the corrected side, so a strict comparison gives a **higher** flip
  rate. Bounded: a uniform +0.05 DockQ rescue recovers 52 flips (2.6%), +0.10 recovers 65
  (3.3%). Worst case 81.9% → ~79.7%.
- **2 designs are unscorable, and that is a finding.** DockQ reports no native interface
  because the true backbone puts the binder *out of contact*: min protein–RNA heavy-atom
  distance 5.27 Å and 5.15 Å against a 5.0 Å cutoff (40 random references: median 2.85 Å,
  max 3.62 Å, none > 5). BoltzGen occasionally emits a design whose binder never touches
  the target.

## 6. Open decisions (deliberately out of scope)

1. **The library.** 1,963 of the 2,400 BG+BIF slots are occupied by designs that fail the
   project's own gate on corrected scores. Whether to re-select is a separate call.
2. **The master / filter store.** Corrected values were **not** backfilled, so every
   downstream tool still serves the inflated BG+BIF numbers.
3. **The pipeline.** `boltzgen_multif.sh` is unchanged; a future BoltzGen campaign would
   reproduce this. The fix is to harvest `spec_<D>.cif` (or store both).
4. **Scope.** Extraction cost ~2 min wall for all seven archives, and the kept working set
   is 0.45 GB — so widening this to the full gated pool is cheap.
