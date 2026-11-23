# Fig5c_interface_chemistry

The two metrics from `~/hits32_browse_flowFC/HANDOVER_METRICS.md` that separate
binders from non-binders, shown as structures plus the numbers behind them.

| Panel | What |
|---|---|
| a | GLU4 carboxylate on A1's phosphate, OE2-OP2 **1.45 A**, OE1-OP1 2.03 A |
| b | ARG12 salt-bridging C6's phosphate at 3.31 A, no carboxylate within 5 A |
| c | the acid-phosphate flag: 0/15 binders, 7/17 non-binders, 9.7% of the 12k |
| d | the five features where natural Csy4 sits away from the library |
| e | flow fold change for all 32, flagged designs marked |

Reproduced from the source data, not copied from the handover: 32 designs,
15 binders / 17 not at the 1.2 cut, flag 0/15 vs 7/17, library 9.7%.

## Which flow number this is, and a correction

`flow_FC` is **LV2.34.1** `ONOFF_avg` = mean(PosA, PosB) / mean(FullNeg, NegA),
from `inputs/lv2_34_1_flow_on_off.csv`. Verified identical to the `flow_FC`
column in the candidate index to 4 decimal places, on all 32.

**This corrects `Fig5b_top3_structures`.** That panel ranked on the LV2.33
arrayed screen (`LV_FlowCytometry/analysis/out/variants.csv`) and its README
calls the index `flow_FC` "of unknown provenance". That was wrong: it is
traceable, to a different and more complete experiment. LV2.33 scored only 20
of 32 because V1 to V12 were not exported; LV2.34.1 has all 32. The handover
treats LV2.34.1 as ground truth and so should we.

Consequence: the top three by LV2.34.1 are twist 23, 15, 16. Fig5b rendered
15, 16 and 26, so **twist 26 should be twist 23** if that panel is kept.

## The flag is AF3-specific, which the handover does not say

`red_acid_phos_per_ifaceres` exists for all three refolders. The split depends
almost entirely on which one is used:

| Refolder | Binders flagged | Non-binders flagged | 12k library |
|---|---|---|---|
| **AF3** (used here and in the handover) | 0/15 | **7/17** | 9.7% |
| Boltz-2 | 1/15 | 3/17 | 9.0% |
| ESMFold2 | 1/15 | 1/17 | 8.9% |

They do not agree design by design: AF3 and Boltz-2 agree on 23 of 32 and flag
the same design only **once**; AF3 and ESM agree on 25 of 32, again one shared
flag. So the seven designs carrying the signal are flagged **by AF3 alone**.

This applies to the structure in panel a as well. `d503930b` reads 0.185 on
AF3, **0.0 on Boltz-2** and 0.059 on ESM: the 1.45 A carboxylate-phosphate
overlap being drawn is a feature of AF3's model of that complex, not a
consensus of the three predictors.

### It does not survive using the refold closest to the design

Rather than fixing one predictor, take per design whichever refold reproduces
the design best (lowest `dockq_LRMSD`, the ligand_RMSD the gate uses) and read the
flag off that model:

| Rule | Binders flagged | Non-binders flagged | Fisher p |
|---|---|---|---|
| AF3 always (as published) | 0/15 | 7/17 | **0.0077** |
| closest refold by LRMSD | 2/15 | 4/17 | **0.659** |
| closest refold by iRMSD | 2/15 | 4/17 | 0.659 |
| flagged by at least 2 of 3 | 0/15 | 2/17 | - |

**The result does not survive.** AF3 is the closest model for only 10 of the 32
(ESM 14, Boltz-2 8), and the closest model averages 1.95 A LRMSD against AF3's
6.57 A. So the published flag is read off a model that is usually not the best
available description of the design.

This hits panel a directly. For `d503930b`, AF3 sits **7.58 A** from the design
while Boltz-2 sits 2.34 A, and Boltz-2 places no carboxylate on a phosphate at
all. The 1.45 A overlap this figure draws comes from AF3's worst-fitting model
of that complex.

That is not proof the metric is wrong. Applied prospectively on AF3 it is
self-consistent, and a filter does not need a physically correct model to
enrich. But the metric's appeal is its mechanism, like charges repel, and a
mechanism claims something about the molecule rather than about one predictor's
error. On the model that best reproduces the design, the repulsion is not there.

**Treat panel a as an illustration of what the AF3-defined flag detects, not as
evidence that this design fails because of electrostatic repulsion.** Before
building a prospective filter on this, decide which model the filter is a
statement about, and check the flag on the closest refold.

`VERIFY:` the unsuffixed `red_acid_phos_per_ifaceres` gives 2/15 vs 9/17 at
22.1% library prevalence, so it is a looser combination across backends. Which
combination it is was not checked.

## The honest headline is the robustness, not the p-values

Straight from the handover, and it must travel with any use of this figure:

- **The 1.2 cut carries both results.** Cutting at 1.0 (26/6) destroys them:
  Csy4 score p 0.010 -> 0.545, acid-phos p 0.0077 -> 0.59. The score is **not
  monotone in FC**; the 1.0-1.2 band is the chemically worst tier.
- **Dropping that band splits the two leads apart.** Acid-phos keeps its
  direction and runs out of n (0/15 vs 2/6, p = 0.071). The Csy4 score halves
  and goes to p = 0.567.
- **Non-binding has more than one cause.** These metrics catch chemically poor
  interfaces. The hardest failures are chemically ordinary and must fail for
  something no structural feature sees. A prospective filter would enrich, not
  guarantee.
- **Never quote the 15-v-11 numbers**: that split was chosen after seeing which
  one maximises the contrast.
- 7 designs carry the whole acid-phos signal; 1 hit in 408 in the exploratory
  scan (q = 0.79); stack-stratified p = 0.075.

## What does not work, so nobody re-runs it
408-feature scan: 15-16 features at p<0.05 where 20.4 are expected, min q 0.79.
ipTM, DockQ, ligand_RMSD, packstat, BUNS, pi-stacking, epitope, geometry: all flat.
Fold clustering is chance (stack-stratified p = 0.24). The 6.99M campaign adds
nothing. `flow_FC` is uncorrelated with the HTS enrichment the 32 were picked
on (rho = -0.155, p = 0.40).

## Structures
Both from `LV_BinderDesign/LV_32Candidates`, AF3 refold (model 2), which is
where the metric is defined and the only model with sidechains for LiMPNN.

| | Design | FC | Why chosen |
|---|---|---|---|
| a | `BG_BIF_len80_d503930b_design39_seq0` | 0.84 | worst acid-phos value of the 32 |
| b | `BG_LiMPNN_len80_e8cb617b_design26_seq2` | 1.56 | top binder, Arg 0.50 / Lys 0.05 at the interface |

**Not claimed: that ARG12 is bidentate.** It makes one NH1-OP2 contact under
3.6 A. The bidentate geometry in the handover is the general rationale for
preferring Arg, not a measurement on this model.

`VERIFY:` panel d plots Csy4's raw feature values against the library
distribution. The handover's score is the mean of five signed z-scores; the
panel shows the inputs, not the score. The score's own numbers (binders +0.373,
non-binders -0.031, stack-stratified p = 0.010) are not recomputed here.

## Colour
`NEED:` the repulsive pair uses the diverging warm end, which COLORS.md rule 4
reserves for signed data. `code/common/CHIMERAX_RENDERING.md` already records
this as an open NEED for clash colouring. Resolve by adding a synbio entry to
`palette.json` under rule 6, then re-render.

## Rebuild
```
python code/Fig5c_interface_chemistry/render.py output/Fig5c_interface_chemistry --width 3000
python code/Fig5c_interface_chemistry/make_figure.py
```
The render must run first: `make_figure.py` reads the two zoom PNGs. Metrics
come from `~/Downloads/design29_742db4d3/an32/{m32,geom_wide}.pkl`, which are
outside the repo.
