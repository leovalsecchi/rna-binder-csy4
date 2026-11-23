# bif_quasicontrol

The BoltzIF scoring defect read as a **quasi-control** on the campaign gate. Two bare
panels.

```bash
~/.venvs/binder_map/bin/python3 code/supp/bif_quasicontrol/make_figure.py
```

Emits `output/bif_quasicontrol/panel_{a,b}_*.{pdf,png}`, both 85 × 40 mm. Runs in
~8 s; every number below is re-derived and printed on every run.

## The argument

Two of the five campaign gates, **ligand_RMSD and sc_DockQ**, were effectively not applied to
BoltzGen + BoltzIF: for that route both were scored against BoltzIF's own Boltz
re-prediction rather than against the generated backbone. One fifth of the ordered
panel therefore entered under a materially weaker filter. If that arm yields candidates
at a rate indistinguishable from the properly gated routes, **the gate did not predict
experimental outcome** — shown by internal comparison rather than asserted.

It matters because Hypothesis 3 is framed as performance *above the unfiltered
baseline*, and **no unfiltered baseline exists anywhere else in this work**. The defect
accidentally produced the closest thing to one.

🚨 **Both cautions are in the caption, not only here.**

1. Those designs still passed the other three gates and were still ipTM-ranked within
   cluster. This is a **weakened filter, not no filter** — a **quasi-control**, never a
   control.
2. They were selected for scoring *well under the wrong reference*, so **regression is
   expected**. A lower rate in that arm would not have been evidence against the gate.

## Sources

| file | mtime | what it gives |
|---|---|---|
| `inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv` | 2026-11-17 22:05 | the 12,000-design panel with BoltzGen + BoltzIF's `sc_rmsd_*` / `dockq_*` already replaced by the true-backbone values, plus `bif_rescored` |
| `inputs/candidates32_mapped.csv` | 2026-11-09 18:56 | the 32 selected for ordering; `hts_comparison` splits them POS_B 25 / POS_A 7 |
| `inputs/csy4_binder_library_12k_v6_k30.csv` | 2026-11-12 14:07 | the **stored** scores, read only to assert the panel is the gated set |

No out-of-repo path. `inputs/MANIFEST.md` says to use the CORRECTED library instead of
the original wherever 12k metrics are needed, and this follows that.

## Verified against the files

- **Per-route candidate counts and denominators.** RFD3 + LigandMPNN 4/3,132 (0.128%),
  BoltzGen + LigandMPNN 11/3,266 (0.337%), **BoltzGen + BoltzIF 9/2,400 (0.375%)**,
  BoltzDesign1 + LigandMPNN 7/2,400 (0.292%), HalluDesign + LigandMPNN 1/802 (0.125%);
  all routes 32/12,000 (0.267%). Identical to `code/Fig4d_stack_hit_rate/README.md`,
  which draws the same five counts from `candidates32_mapped.csv` — see **Open** below.
- **The 32 are consistent across two independent lists.** `candidates32_mapped.csv` and
  `~/Downloads/design29_742db4d3/an32/m32.pkl` (2026-10-21 11:43; `hit` / `pos` columns)
  agree on all 32 `design_id`s and on the 7 / 25 POS_A / POS_B split. Only the repo file
  is read by the script; the pickle was used once to cross-check and is not a dependency.
- **χ² = 4.6012, df = 4, p = 0.3307**, smallest expected candidate count **2.14**.
  Hand-rolled Pearson χ² with the closed-form even-df upper tail (**scipy is unusable in
  `~/.venvs/binder_map`** — `_propack` dlopen failure — so every statistic here is
  hand-rolled and the script says so).
- **Permutation p = 0.3249**, 20,000 candidate-label permutations over the 12,000
  designs, seed 42. The asymptotic value is confirmed; the two agree to two decimals.
  ⚠️ **Only two decimals of a permutation p are real here**: Monte-Carlo SE at
  n = 20,000 is 0.0033, and a 200,000-draw run of the same estimator gives
  **0.3244 ± 0.0011**. A figure of 0.3306 quoted elsewhere is the same number one draw
  away, not a different result — but 0.324, not 0.331, is where this estimator sits.
- **1,963/2,397 (81.9%) of the BoltzGen + BoltzIF panel members fail the corrected
  gate**, 434 pass. Identical to `code/supp/bif_effect`, same predicate, same
  denominator: `(sc_rmsd_af3 ≤ 3 OR sc_rmsd_esm ≤ 3) AND (dockq_af3 ≥ 0.7 OR
  dockq_esm ≥ 0.7)`, on stored columns and then on corrected ones. 2,397 is the
  denominator because three members have no true-backbone score (one genuinely
  unmeasurable, two `no_native_interface`); **none of the three is a candidate**, so the
  choice cannot move panel (b)'s numerator.
- **The CORRECTED library reproduces bif_effect's own join exactly.** Rebuilding the
  corrected columns from `csy4_binder_library_12k_v6_k30.csv` +
  `inputs/large/bif_true_per_backend.csv.gz` (2026-11-17 15:18) the way bif_effect does
  gives `max |difference| = 0` on all four columns and the same three NaN rows. That is
  what stops the two figures drifting.
- **Panel (b) strata.** Four properly gated routes 23/9,600 (0.240%); BoltzGen + BoltzIF
  passing the corrected gate 1/434 (0.230%); BoltzGen + BoltzIF failing it
  8/1,963 (0.408%). Collapsed, **fails 8/1,963 (0.408%) vs passes 24/10,034 (0.239%),
  ratio 1.70×, Fisher two-sided p = 0.2261**. Three-group χ² = 1.7505, df = 2, p = 0.4168
  (smallest expected 1.16). Within BoltzGen + BoltzIF alone, Fisher p = 1.0000.
  Fisher is used rather than a Monte-Carlo permutation because for a 2 × 2 with fixed
  margins it *is* the exact permutation test.
- **The panel is the gated set.** All 12,000 pass the stored-score predicate, asserted
  in `load()` off the pre-correction library. That is what makes "weakened filter"
  literally true.

## What the result means — state it this way

The weakly gated arm has the **highest** candidate rate of the five, and the five routes
are **statistically indistinguishable**. 🚨 With 32 candidates over 12,000 designs the
test is badly underpowered, so **"indistinguishable" here is a statement about power as
much as about biology**: equivalence was not demonstrated, a difference merely was not
detected. The smallest expected candidate count is 2.14 (HalluDesign + LigandMPNN),
which is why the exact test is reported beside the asymptotic one.

## Panel (c) was tested and is NOT drawn

Candidate recovery against corrected ligand_RMSD and sc_DockQ across the whole panel does not
hold up, in two separate ways, so no third panel is emitted.

- **Deciles are noise.** Binned by corrected best-of-two ligand_RMSD, candidates per decile
  of ~1,200 run 4, 3, 2, 3, 4, 2, 4, 2, 2, 6; by corrected sc_DockQ, 6, 6, 2, 3, 2, 5, 2, 2,
  3, 1. Nothing monotone survives 32 events spread over ten bins.
- **The one nominal signal is a route-mix artefact.** Whole panel, corrected sc_DockQ:
  candidates sit lower than the rest (median 0.732 vs 0.758), permutation p = 0.0372
  (20,000 draws, seed 42). It is Simpson's: BoltzGen + BoltzIF has both the lowest
  corrected sc_DockQ (median 0.553 against 0.761–0.777 elsewhere) and the highest candidate
  rate. Within the four properly gated routes p = 0.1690; within BoltzGen + BoltzIF alone
  p = 0.1081. Corrected ligand_RMSD is flat everywhere (p = 0.6468 / 0.6069 / 0.2947).
  Uncorrected for two metrics and driven entirely by the confound, 0.0372 must not be
  drawn or quoted as a finding.

Worth keeping in prose rather than as a panel: of the nine BoltzGen + BoltzIF
candidates, **eight fail the corrected gate**, and their median corrected scores are
19.3 Å ligand_RMSD and 0.059 sc_DockQ against 5.8 Å and 0.555 for the rest of the route. n = 9.
(The route's median sc_DockQ quoted above, 0.553, is over all 2,397; 0.555 is the 2,388
non-candidates.)

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the n. 193 words.*

**A design route whose structural gates were effectively not applied yielded candidates
at the highest rate of the five.** BoltzGen + BoltzIF's ligand_RMSD and sc_DockQ were scored
against a re-prediction rather than the generated backbone, so two of the five gates did
not bind for that route. (**a**) Fraction of each route's panel members selected for
ordering, labelled candidates / members; the routes are statistically indistinguishable
(χ² = 4.60, df = 4, p = 0.33; 20,000-permutation p = 0.32), with the weakly gated route
highest. (**b**) The same fraction split by **corrected** gate status: the 1,963
BoltzGen + BoltzIF members failing once rescored against their own backbones yield
0.408%, against 0.239% for those passing (Fisher p = 0.23). n = 12,000 panel members and
32 candidates in (**a**); 11,997 in (**b**), three BoltzGen + BoltzIF members lacking a true-backbone score. A
**quasi-control, not a control**: those designs passed the other three gates and were
ipTM-ranked within cluster, so the filter was weakened rather than absent, and they were
selected for scoring well under the wrong reference, so regression is expected. At 32
candidates the comparison is underpowered: no difference was detected, none excluded.

## Open

- ⚠️ **Panel (a) duplicates `Fig4d_stack_hit_rate`.** The same five counts, from the same
  candidate list. Decide which figure owns them before the set ships — the
  metric_reliability trap. The case for keeping it here is that it is the baseline panel
  (b) is read against and it carries the test; the case against is that Fig 4d is in the
  main text. If Fig 4d keeps them, panel (a) becomes a cross-reference and this figure
  is panel (b) alone. **They differ in two ways that must be reconciled either way:**
  Fig 4d sorts by rate (a deliberate break of the fixed route order, on Leo's request)
  while this uses `palette.STACK_ORDER`, and Fig 4d's README calls the quantity a
  *screen hit rate*.
- `NEED:` **Fig 4d's caption describes the 32 as designs "passing the gate
  (>= 75 reads in POS and matched NEG, splicing >= 0.07, ON/OFF fold change >= 10)".**
  That is a different account from the one this figure is written against — a
  **selection** on ranked splicing ratio and fold change, not a threshold hit rate — and
  the two cannot both stand. Not re-derivable from anything in `inputs/`; Leo has to say
  which is right. Nothing here was changed in `Fig4d_stack_hit_rate`.
- `NEED:` **a citing sentence.** This figure only does its job if the Hypothesis 3
  discussion points at it where it claims performance above an unfiltered baseline.
- The caption's figure number is a placeholder. Supplementary numbering is not settled,
  so no S-number appears in any file or directory name here.

## Naming, settled 2026-11-20

**`sc_DockQ` and `ligand_RMSD (Å)`**, matching the main text (Fig 2c), wherever these two
quantities are an axis label or appear in prose. The bare `DockQ` and the hyphenated
`ligand-RMSD` are not used as display forms anywhere in `code/supp`.

Renamed **`sc_RMSD` -> `ligand_RMSD` on 2026-11-21**: the gate's first term is the
RNA-aligned *ligand* RMSD (`dockq_LRMSD_af3` / `_esm`), a pose metric, not the Ca fold
RMSD it was once written up as. Same measurement throughout, new name.

⚠️ **CSV column names are untouched and stay untouched** — `sc_dockq`, `sc_rmsd`,
`dockq_af3`, `dockq_esm`, `rmsd_af3`, `rmsd_esm`, `dockq_self`, `dockq_LRMSD`,
`dockq_fnat`, `dockq_iRMSD` and the `_true` variants all keep their lower-case spelling,
because they are what is on disk. The sweeps only ever matched capitalised display forms,
never a column name.

## Vocabulary, settled 2026-11-19, swept 2026-11-20

**"re-predictor"** for AF3 / Boltz-2 / ESMFold2, never "refolder", "predictor" or
"repredictor". **"refolding seeds"** for the five seeds. **"design route"**, never
"stack". In all prose, labels and legends.

⚠️ **Identifiers are deliberately exempt. Do not "fix" them — each one breaks a read or
an import:** the `refolder` and `stack` columns in the source CSVs, `palette.REFOLDERS`,
`palette.REFOLDER_ORDER`, `palette.STACK_ORDER`, `palette.SHORT`, the source file names
that contain `_by_stack_refolder`, and the directory name `Fig5fg_h1null_refolder`.
