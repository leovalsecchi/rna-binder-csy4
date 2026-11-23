# Fig 5 f-i — H1 is null; of H2, only the sequence part holds

Four bare panels + `Fig5fgh_assembly` for the 178 mm block. Data in `inputs/h1h2/`,
written by `export_data.py` here. Source analysis:
`~/hits32_browse_flowFC/FLOW_FC_ANALYSIS.md`, appendices 4-6.

## The two hypotheses
- **H1** (the field's default): a designed binder works if the predicted complex
  scores well — ipTM, DockQ, ligand-RMSD, interface area, contacts, packstat, BUNS.
- **H2**: it works if the interface carries the chemical groups natural RBPs use,
  in natural proportions, whatever any predictor thinks.

## The panels
| | what it shows |
|---|---|
| **f** `scan` | every feature tested binder-vs-non, as a uniform QQ plot. On the diagonal and inside the envelope = no signal anywhere. Drawn for both scorings. |
| **g** `h1metrics` | the ten metrics the field ranks on, as P(binder > non) with a 95% bootstrap interval. Nine of ten cover 0.5. |
| **h** `h2legs` | each leg of H2 under three stresses. Structural legs collapse; sequence legs cannot move. |
| **i** `arom` | the survivor — aromatic fraction, for the library, the 32, and the 59 natural RBPs. |

## Decisions worth keeping
- **QQ, not a volcano, for panel f.** The claim is negative. A volcano shows the
  smallest few P values out of context and invites belief in them; a QQ against
  uniform shows the whole batch is noise. The envelope (2,000 draws of n uniforms,
  2.5/97.5 per order statistic) is what makes it a test rather than an assertion.
  Both series sit ON or BELOW the diagonal: fewer small P than chance predicts.
- **P(binder > non) in panel g**, not a difference of means. n is 15/17, the
  metrics have incomparable units and several are bounded or zero-inflated. The
  probability of superiority is the same rank statistic the permutation tests use,
  is bounded with a fixed null at 0.5, and reads as "how often does a binder beat
  a non-binder here".
- **No direction arrows in g.** The axis is P(binder > non), so direction is already
  in the mark's position; an extra ↓ on ligand-RMSD made the label wrap and confused
  two different directions.
- **AF3 drawn as an open ring in h.** On the two sequence legs the AF3 and
  best-refold marks coincide exactly — that IS the finding (no structure enters
  them), and a filled mark would hide it.
- **Mathtext subscripts are banned in these labels.** `$-\log_{10}P$` renders the
  subscript at 0.7x, which trips `page_assembly.report`'s stray-size check
  ("f: stray sizes {3.85}"). Plain `-log10 P`. Mathtext also swallows the space
  after `$P$`, giving "Punder the null".
- **`variant`, not `pair`.** Panels g and h are y-label heavy; `pair` rewrites
  font-size without moving labels, so right-anchored text grows out of the
  viewBox and the page clip eats it (see Fig5de's README).

## Caption draft
**H1 is null; of H2, only the sequence component holds.** (f) Permutation P values
for every feature tested between the 15 flow binders and 17 non-binders
(FC ≥ 1.2; 20,000 label permutations), against the uniform expectation. Grey band,
95% envelope of the order statistics under the null. Both scorings lie on or below
the diagonal, and both find fewer features at P < 0.05 than chance predicts.
(g) The ten metrics commonly used to rank designed binders, as the probability that
a binder exceeds a non-binder, with 95% bootstrap intervals (5,000 resamples). Nine
of ten intervals cover 0.5. ligand-RMSD is the exception, in the direction of better
self-consistency among binders; one nominal hit in ten is what chance gives.
(h) Each component of H2 under three nested robustness checks: AF3 used for every
design, the refold that best reproduces each design (DockQ against the designed
complex), and that refold restricted to the 23 designs whose design structure is
unambiguous (BG+BIF excluded, its stored design being a re-prediction). Shaded,
P ≥ 0.05. Every structurally measured leg loses significance; the sequence legs are
unchanged because no structure enters them. Only the pre-specified AF3/cut-1.2 cell
was planned in advance — the others are robustness checks, not independent tests.
(i) Aromatic fraction, the one feature that survives all three checks: the 12k
library, the 32 ordered designs by flow outcome (bars, group means), the 59 natural
RBP-RNA complexes, and Csy4. Natural interfaces sit where the designs do not.
