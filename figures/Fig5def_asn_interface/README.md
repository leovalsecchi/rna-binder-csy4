# Fig 5 f-h — three steps: the metrics are flat, only Asn moves, and it moves on both measures

Three bare panels + `Fig5_assembly` (178 mm: f | g on top, h full width). Data
`inputs/asn/`, written by `export_data.py` and `export_aa20.py`. Analysis:
`~/hits32_browse_flowFC/FLOW_FC_ANALYSIS.md` appendices 8-9e.

**Supersedes `code/Fig5fgh_h1null_h2hold/`** (built on an aromatic result that is Asn in
disguise, and annotating a circular ligand-RMSD P value). That folder can be deleted.
`partition.svg`, `backends.svg`, `stress.svg`, `dist_total.svg` and `dist_interface.svg`
are still emitted by the scripts here and are good panels — they are just not in the
three-panel cut.

| | panel | step |
|---|---|---|
| f | `metrics` | **nothing standard works.** Twelve metrics the field ranks designed binders on, as P(binder > non) with 95% bootstrap intervals. All twelve cover 0.5. |
| g | `aminoacid_volcano` | **of the amino acids, only Asn.** All 19 present in these designs as a volcano: effect on x, significance on y, BH threshold as a rule. `aminoacid.svg` (ranked dot rows) is the same data in the alternative form and is still emitted. |
| h | `asndist` | **Asn is enriched and trends on both measures.** Three densities per measure, each group's designs riding its own curve, with the natural-RBP median as the reference line. |

## Reference brief (2026-11-14)
Searched for how this field draws "do computational metrics predict experimental
success" and "many features tested, which moved".

- **Overath et al., bioRxiv 2025.08.14.670059**, a meta-analysis of 3,766
  experimentally characterised de novo binders with 200+ features per design — the
  closest published analogue to this analysis. **Borrowed**: a separation statistic
  per feature as the primary display, with an explicit uncertainty band (their Fig. 3D/3E
  use median average precision with IQR across targets) and direction-of-association
  marked rather than left implicit. **Not borrowed**: their ranked *bars*. They use
  average precision, which has a meaningful zero; our statistic is P(binder > non) whose
  null is 0.5, and a bar encodes length from zero (STYLE_NOTES §8), so a bar would
  misrepresent it. Dots with intervals, against a 0.5 rule, is the faithful translation.
  They also use violin plots with a p-value per panel for binder-vs-non distributions,
  and no multiple-testing correction — we do correct, and panel g shows the threshold.
- **Volcano-plot convention** (VolcaNoseR, Nat. Commun. 2020 and the general
  proteomics/metabolomics literature): effect size on x, −log10 P on y, hits labelled,
  FDR-adjusted threshold drawn and pre-specified. **Borrowed** wholesale for panel g,
  with the effect axis changed from log fold change to P(binder > non) because these are
  fractions with a natural null and no meaningful ratio.
- **Cao et al., Nature 2022, Extended Data Fig. 13** "Power of computational metrics to
  predict binders" is the obvious precedent for panel f but could not be read — the PDF
  exceeds the fetch limit and the Nature figure pages are behind authentication. Noted
  as unverified rather than described from memory.

## Decisions worth keeping
- **ligand-RMSD and DockQ per backend, never the better of the two.** The refold used
  elsewhere is picked by lowest ligand-RMSD; testing the selected value tests the variable
  that did the selecting. The previous figure did that and manufactured P = 0.038.
- **Panel g is 19 tests and says so on the canvas.** Asn reaches q = 0.046 (sequence) and
  0.064 (interface). **His is nominally significant (P 0.049 / 0.021) and does not survive
  correction (q 0.20)** — it is drawn and labelled, not hidden. Two of nineteen at
  P < 0.05 where one is expected.
- **Cys is absent because the design protocol excludes it**, not because it was dropped.
- **Panel h is percentages, not counts**: naturals are 61-1301 aa with ~41-residue
  interfaces against the designs' 80 aa and ~21. For a design 1.25% of sequence = one
  residue, which the footnote states.
- **Large groups as a density, small groups as their points.** With n = 15 a density
  invents shape that is not in the data.
- **Every row carries "% with none"**, because both distributions are zero-inflated and
  the mean misdescribes them (the library's 4.59% interface mean sits in a trough between
  a spike at zero and a broad shoulder).

## The honest shape, which the caption must carry
The binders are **not** unusually Asn-rich: 5.71% of interface residues against a
600-design library sample at 4.59% (P = 0.38 against the library). **The inactive designs are
the anomaly** — 1.49%, 0.32x the library, P = 0.003, with **71% having no interface
asparagine at all** against 41% of the library and 20% of binders. This predicts failure
from a missing contact, not success from a clever one.

## Caption draft
Everything stripped from the canvas lives here. Nothing below is drawn in the panels.

**No metric used to rank designed binders predicts binding; asparagine does, in sequence
and at the interface.** Throughout, the statistic is the probability that a randomly
chosen flow binder exceeds a randomly chosen inactive design (FC >= 1.2; 15 proposed
binders, 17 inactive designs) on that quantity; **0.5 is no information**.
(**f**) Twelve metrics commonly used to rank designed protein-RNA binders, with 95%
bootstrap intervals from 5,000 resamples. **All twelve intervals cover 0.5.** ligand-RMSD and
DockQ are given per structure predictor, never as the better of the two, since the refold
used elsewhere is selected on ligand-RMSD. At this sample size the design has 77% power at a
true value of 0.75 and 34% at 0.67, so this is evidence against a large effect, not
against any effect.
(**g**) The same statistic for each of the **19 amino acids these designs contain**
(cysteine is excluded by the design protocol), for whole-sequence composition (filled) and
for composition among interface residues (open; protein heavy atom within 4.5 A of RNA).
P values are 20,000 label permutations stratified within design stack, corrected across
the 19 by Benjamini-Hochberg. **Two reach P < 0.05 where one is expected by chance.**
Asparagine is the only one to survive correction (q = 0.046 sequence, 0.064 interface);
histidine is nominally significant (P = 0.049 and 0.021) and does not (q = 0.20).
(**h**) Asparagine as a percentage of the sequence (left) and of the interface residues
(right), for the 12k library, the 32 ordered designs split by flow outcome, and 59 natural
RBP-RNA complexes. Densities for the large samples, individual designs for the small ones;
vertical rule, group mean. Designs are 80 aa, so 1.25% of sequence is one residue;
interfaces average ~21 residues in the designs and ~41 in the naturals, so percentages
rather than counts are comparable across groups.
Both distributions are zero-inflated and the fraction with **none at all** is the clearer
summary: sequence 47% of inactive designs, 21% of the library, 13% of binders, 2% of naturals;
interface **71% of inactive designs, 41% of the library, 20% of binders**, 24% of naturals.
Group means, sequence: 0.96 / 2.13 / 2.50 / 4.59%; interface: 1.49 / 4.59 / 5.71 / 6.21%.
Note that the binders are close to the library and the **inactive designs are the outlying
group** — the pattern predicts failure from a missing contact rather than success.

## Known limits, stated in the appendices not on the canvas
- Asn came out of a 247-409 feature scan; q over that larger scan is ~0.34.
- What the asparagines contact is **not** resolved: of 22 distinct Asn-base contacts only
  7 appear in both the AF3 and the ESM prediction, and 4 survive an ligand-RMSD <= 3 A gate.
  Do not add a base-identity panel (Appendix 9d).
- Power: at 15 vs 17 this design has 77% power at AUC 0.75 and 34% at 0.67, so panel f
  means "no large effect", not "no effect".
