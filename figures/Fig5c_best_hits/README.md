# Figure 5 panel c: best three by flow fold change, plus twist 7's asparagine

Wide four-up strip. Supersedes `Fig5b_top3_structures`, which was built on the
LV2.33 arrayed readout and therefore picked the wrong three.

## What is shown

Ranking is `flow_FC` in `LV_32Candidates/index.tsv`, which is the LV2.34.1
`ONOFF_avg` column, mean(PosA, PosB) / mean(FullNeg, NegA), verified identical
on all 32 candidates.

| panel | twist | FC | stack | design |
| --- | --- | --- | --- | --- |
| 1 | 23 | 1.56 | BG+LiMPNN | e8cb617b_design26_seq2 |
| 2 | 15 | 1.54 | BG+BIF | 25a69495_design18_seq0 |
| 3 | 16 | 1.47 | RFD3+LiMPNN | 73d23e78_design49_seq4 |
| 4 | 07 | 1.41 | BD1+LiMPNN | c5072321_design5_seq1 |

Panel 4 is the asparagine panel requested for the right-hand slot.

## Which model each panel uses

The model with the lowest self-consistency RMSD (`dockq_LRMSD_*`), not AF3 by
default. AF3 is the worst model for three of the four.

| twist | sc_rmsd_af3 | sc_rmsd_esm | used |
| --- | --- | --- | --- |
| 23 | **1.548** | 16.031 | AF3 |
| 15 | 2.412 | **2.079** | ESM |
| 16 | 7.857 | **1.830** | ESM |
| 07 | 3.616 | **2.387** | ESM |

Verified 2026-11-16 against `LV_ThesisFigures/inputs/csy4_binder_library_12k_v6_k30.csv`
(columns `sc_rmsd_af3`, `sc_rmsd_esm`). There is no Boltz2 value: the candidate
PDBs hold exactly three models, and their own REMARK says
`model 1 = design, model 2 = AF3 refold, model 3 = ESMFold2 refold`, all
superposed on the RNA C1' in one common frame. An earlier version of this file
listed a `b2` column; those numbers could not be traced to any file on disk and
have been removed. The choice of model is unaffected.

## Why only two panels carry a dash

ChimeraX `hbonds` on the chosen models. Re-run 2026-11-16; the script and its
full log are kept here as `verify_hbonds.cxc` and `verify_hbonds.log`.

`VERIFY:` the structures carry no hydrogens, so ChimeraX relaxed its
constraints by 0.4 A and 20 degrees and could not evaluate donor-H..acceptor
angles. This is a donor-acceptor distance and geometry screen, not a full
hydrogen bond determination. Word any thesis claim accordingly.


| twist | result |
| --- | --- |
| 23 | 3 pass, but ARG12 NH2-G7 N7 is 1.94 A and LYS23 NZ-A11 OP1 is 1.42 A, overlap distances rather than bonds. Only TYR24 OH-A11 OP2 at 2.59 A is clean. |
| 15 | 0. ARG67 sits 2.43 A from G12 N7 but fails the angle test. |
| 16 | 1. LYS42 NZ-U10 OP1 at 3.11 A. The G7 O6 contact (2.44 A) does not pass. |
| 07 | 0. ASN54's closest contact is the A9 phosphate at 2.78 A; its G7 O6 reach is 3.39 A. |

Only two contacts are both geometrically valid and at a real bond distance, and
only those two are dashed: twist 23 TYR24 OH-A11 OP2 (2.59 A) and twist 16
LYS42 NZ-U10 OP1 (3.11 A). Twist 15 and twist 7 carry no dash.

`VERIFY:` the figure legend must state that a dash marks a contact passing the
ChimeraX `hbonds` geometry test. Without that sentence the two dashless panels
read as a rendering oversight instead of as the result they are.

## Run

    python3 render.py ../../output/Fig5c_best_hits --width 2400
    python3 make_figure.py ../../output/Fig5c_best_hits ../../output/Fig5c_best_hits

`--sweep` on render.py writes 12 angles per panel for re-picking cameras.
Chosen: twist 23 y060, twist 15 y060, twist 16 y090, twist 7 y030.

Rendering follows `code/common/CHIMERAX_RENDERING.md`. All four complexes are
aligned by RNA C1' onto twist 23's hairpin and share one canonical camera, so
the per-panel `turn y` angles are comparable between panels.

## What each panel shows

Moved off the figure so the panel carries only the twist number and the fold
change. Residue and nucleotide labels stay on the structures.

| panel | stack | highlighted | contact |
| --- | --- | --- | --- |
| twist 23 | BG+LiMPNN | TYR24, LYS23, ARG12 on G7 and A11 | TYR24 OH to the A11 phosphate, 2.59 A, dashed. LYS23 and ARG12 sit at 1.42 and 1.94 A, overlap distances, not dashed. |
| twist 15 | BG+BIF | ARG67 in the major groove at G12 | ARG67 NH2 to G12 N7, 2.43 A, fails the angle test, not dashed. |
| twist 16 | RFD3+LiMPNN | LYS42, ARG50 on G7 and U10 | LYS42 NZ to the U10 phosphate, 3.11 A, dashed. The G7 O6 contact at 2.44 A does not pass. |
| twist 7 | BD1+LiMPNN | ASN54, ARG53 on A9 and G7 | ASN54 ND2 to the A9 phosphate, 2.78 A, and to G7 O6 at 3.39 A. Neither passes, not dashed. |

Dashes are ChimeraX `distance` pseudobonds with the measured value printed to
one decimal. The legend must say that a dash marks a contact passing the
ChimeraX `hbonds` geometry test, otherwise the two dashless panels read as an
oversight rather than as a result.

## Handover, 2026-11-16

### Files

| file | what |
| --- | --- |
| `output/Fig5c_best_hits/fig5c_best_hits.{pdf,svg,png}` | the panel, 170 x 48 mm |
| `output/Fig5c_best_hits/t{23,15,16,07}.png` | the four renderings, 2400 px |
| `output/Fig5c_best_hits/render_besthits.cxc` | the ChimeraX script the last run executed |
| `render.py`, `make_figure.py` | regenerate everything |
| `verify_hbonds.cxc`, `verify_hbonds.log` | the H-bond check behind the dashes |

Rebuild from scratch:

    python3 render.py ../../output/Fig5c_best_hits --width 2400
    python3 make_figure.py ../../output/Fig5c_best_hits ../../output/Fig5c_best_hits

### What is verified, and against what

- Ranking, FC values, generator and design_id for all four: `LV_32Candidates/index.tsv`,
  columns `flow_FC` and `flow_FC_rank`. Twist 23, 15, 16 are ranks 1 to 3; twist 7 is
  rank 4, included for the asparagine, not because it is top three.
- `flow_FC` is the LV2.34.1 `ONOFF_avg` column, mean(PosA, PosB) / mean(FullNeg, NegA).
  It is NOT the LV2.33 arrayed readout.
- Model choice: `inputs/csy4_binder_library_12k_v6_k30.csv`, `sc_rmsd_af3` and
  `sc_rmsd_esm`. Each panel uses the lower of the two.
- Contacts and their H-bond status: `verify_hbonds.log`, re-run 2026-11-16, reproduces
  the values printed on the figure exactly.

### Open markers

- `VERIFY:` the legend must say a dash marks a contact passing the ChimeraX `hbonds`
  test. Without it, the two dashless panels read as an oversight rather than a result.
- `VERIFY:` the hbonds test ran on structures with no hydrogens, so constraints were
  relaxed by 0.4 A and 20 degrees and donor-H..acceptor angles were not evaluated.
  It is a donor-acceptor distance and geometry screen. Per
  `LV_writing/style/claim-boundaries.md`, do not upgrade this to a binding claim.
- Twist 23's two other passing contacts (1.42 and 1.94 A) are overlap distances in an
  AF3 model. They are deliberately undashed. If a reviewer asks why the AF3 panel has
  the tightest contacts, that is the answer.

### Superseded

`code/Fig5b_top3_structures/` ranked on the LV2.33 arrayed readout and so picked
twist 26 instead of twist 23. Left on disk, not deleted. Do not use it for panel c.

### Trap: do not resample these PNGs with lanczos

`imshow` draws the renderings at `interpolation="antialiased"`. With `lanczos`
(2026-11-16, fixed 2026-11-17) every silhouette and every ChimeraX label came out
edged in a dotted white stipple. The kernel has negative lobes and matplotlib
resamples RGBA channel by channel without premultiplying, so the ringing drives
alpha below zero on a high-contrast edge and those pixels go fully transparent.
The source PNGs are clean, so the artefact only appears after assembly. Any
transparent-background ChimeraX render is exposed to this.
