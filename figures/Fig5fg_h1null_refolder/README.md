# Fig 5 f-g — no design metric predicts binding, and the structural ones are not reproducible

Two bare panels + `Fig5fg_assembly` (178 mm, one row). Data `inputs/h1null/`, written
by `export_data.py` here. Source: `~/hits32_browse_flowFC/FLOW_FC_ANALYSIS.md`,
appendices 4-7.

**Replaces `Fig5fgh_h1null_h2hold` (4 panels), which oversold a result that does not
survive. Delete that folder when this one is accepted.** What went wrong there, kept
so it is not repeated:
1. Its panel f was a QQ plot of two scorings — methodology dressed as a result, and
   unreadable without explanation.
2. Its panel g annotated ligand-RMSD at P = 0.038. That was **circular**: the value tested
   was the minimum over the two refolders, selected by that same minimum. Unselected it
   is P = 0.85 on AF3 and 0.26 on ESM. Never test a fidelity metric on the refold you
   chose with it.
3. Its panel i showed aromatic fraction as "the survivor". The effect is entirely a
   lower tail — fewer than two aromatic residues: binders 2/15, non 10/17, P = 0.012;
   above four residues: 2/15 vs 4/17, P = 0.66. The most aromatic designs in the set are
   non-binders. And `arom` correlates with `if_arom` at +0.66, so "sequence yes,
   interface no" is noise splitting one weak signal, not two findings.

## The panels
| | what it shows |
|---|---|
| **f** `h1null` | the ten metrics the field ranks designed binders on, each scored on BOTH refolders, as P(binder > non-binder) with 95% bootstrap intervals. 0 of 20 tests reach P<0.05; all 20 intervals cover 0.5. |
| **g** `refolder` | how reproducibly each interface measurement comes out when the same sequence is refolded by AF3 and by ESMFold2, across all 12,000 designs. Acidic-at-interface — the measurement behind the one "significant" structural result — is the least reproducible of the eight. |

## Decisions worth keeping
- **Both backends in f, never the selected one.** Removes the circularity above, and
  makes g's point without arguing it.
- **P(binder > non) not Δmean**: n = 15/17, incomparable units, bounded and
  zero-inflated features. Same rank statistic as the permutation tests, null fixed at
  0.5, reads directly.
- **Spearman over 11,999 designs in g, not over the 32.** The claim is about the
  measurements, not the binders, so it uses every design; the correlations are precise
  and the 32 contribute nothing.
- Panel g's annotation carries the acid–phosphate concordance because that is the
  consequence a reader needs and no chart of it is clearer than the sentence.

## Caption draft
**No metric used to rank designed binders separates binders from non-binders, and the
structural ones are not reproducible between refolders.** (f) Ten metrics commonly used
to rank designed protein–RNA binders, each computed on the AlphaFold3 and the ESMFold2
refold of the same design, as the probability that one of the 15 flow binders
(FC ≥ 1.2) exceeds one of the 17 non-binders, with 95% bootstrap intervals (5,000
resamples). All twenty intervals include 0.5 and none of the twenty permutation tests
reaches P < 0.05. A scan of all 409 features behaves the same way: 16 at P < 0.05 where
20.5 are expected by chance. (g) Spearman correlation between the AlphaFold3 and
ESMFold2 values of each interface composition measurement, over the 11,999 designs with
both refolds. The two models place about 60% of the same residues at the interface
(median Jaccard 0.65), so composition agrees well for abundant residue classes and
poorly for the sparse ones. Acidic residues at the interface, the measurement
underlying the acid–phosphate contact flag, agree least. That flag marked 0 of 15
binders and 7 of 17 non-binders on AlphaFold3 models (P = 0.008) and 1 of 15 versus 3
of 17 on the refold that best reproduces each design (P = 0.60).
