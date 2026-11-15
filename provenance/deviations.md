# Deviations and mid-campaign changes

The campaign ran from May to August 2026 and the pipeline changed while it ran. A
design's properties therefore depend on **when** it was generated, not only on which
stack generated it. This file records every change that alters what a stored number
means. It is the companion to `S1.md`: S1 says what the software was, this says what it
did differently before a given date.

Dates are the date the change took effect in production.

## Changes to generation

| Date | Change | What it affects |
|---|---|---|
| 2026-06-11 | **BoltzIF Route B.** The BoltzGen worker was patched (`writer.py`, harvest step, `--steps`) to emit ten sequences per backbone instead of one. | BG+BIF designs generated before this date are 1 sequence per backbone, not 10. Per-backbone sequence diversity is not comparable across the boundary. |
| 2026-06-11 | **BoltzDesign1 internal cutoffs disabled.** The bundled ipTM and pLDDT acceptance thresholds were lowered from 0.5 / 0.7 to 0.0, because at 0.5 the stack rejected every design it produced. | The stack itself was not changed — only its self-filter. BD1 went on to produce 1,444,614 designs. BD1 output is therefore unfiltered by its own confidence criteria, unlike the upstream default. Setting the thresholds to 0.0 also suppresses the `pdb/` output, so backbones are taken from `results_final/` CIFs. |
| 2026-06-25 | **Cysteine excluded.** `omit_AA: C` (LigandMPNN) and `inverse_fold_avoid: C` (BoltzIF) were set. | Designs generated **before** this date contain cysteine. Cysteine-containing designs survive in the earlier campaigns and are identifiable only by run timestamp. |
| 2026-06-25 | **RFdiffusion3 + LigandMPNN chain fix.** `--chains_to_design A` had been passed where chain A is the RNA, so the inverse-folding stage was designing nothing. | RFD3+LiMPNN sequences produced before this date are not designed sequences. |
| 2026-06-27 | **HalluDesign switched from the helix route to de novo.** | Only de-novo HalluDesign designs are in the analysed dataset. The helix-era runs are not part of it. |
| as run | **HalluDesign refinement time steps set to 60** (`configs/halludesign/random_init.yaml`), against an upstream default of 50. | Stated for the record. |

## Changes to metrics

| Date | Change | What it affects |
|---|---|---|
| 2026-08-04 | **DockQ convention changed to RNA-aligned.** Superposition is performed on the RNA chain rather than the previous convention. | `dockq_self` / `dockq_LRMSD` are target-fixed from this date. The older `dockq_std_*` columns are retained but are a different quantity. |
| 2026-08-26 to 27 | **DockQ RNA residue-register bug fixed** in `pipeline/scripts/lrmsd_rna_aligned.py`. Generators number the RNA differently — RFdiffusion3 1–16, BoltzGen 2–17, BoltzDesign1 and HalluDesign 1–18 — while every re-prediction numbers it 1–16, so residues were being paired out of register. | **Any DockQ figure produced before 2026-08-27 is stale and must be regenerated.** After the fix the gated pool fell from 29,992 to 28,315 designs, and the HalluDesign arm from 830 to 430 — below the ~500-design mark at which a stack can still support a hit-rate comparison. Claims about the HalluDesign arm should be hedged accordingly. |
| 2026-08-28 | **Campaign v8 re-collected.** The original collection was stale. | HalluDesign v8 rows went from 12,150 to 182,699. Any v8 analysis predating the re-collection undercounts HalluDesign by an order of magnitude. |
| 2026-07-02 | **ipSAE sidecar introduced.** ipSAE requires the full L×L PAE matrix, which is not retained in the standard output. | ipSAE-family columns exist **only for designs re-predicted after this date**, not for the whole set. |
| late, after analysis | **Boltz-2 dropped from selection.** It was excluded as the least reproducible backend while reporting the highest ipTM. | Wherever the analysis says "both backends" or "either backend", it means **AlphaFold3 and ESMFold2**. The filter store carries only `_af3` and `_esm` columns. Boltz-2 numbers exist but did not enter the gate. |

## Structural caveats that are properties of the data, not of a date

**BG+BIF design coordinates are a re-prediction, not a designed backbone.** The BoltzIF
harvest takes its file from `bg_output/.../refold_cif/`, which is a Boltz prediction
selected as best-of-N by `0.8 · ipTM + 0.2 · pTM`, with pLDDT in the B-factor column.
Self-consistency RMSD and DockQ computed against it are therefore prediction-vs-prediction,
not design-vs-prediction. Rescoring all 2,400 delivered BG+BIF members against their
true backbones (2026-09-30) moved the median self-consistency RMSD from 2.744 Å to
8.663 Å and the median DockQ from 0.709 to 0.406 — **1,963 of 2,397 members, 81.9%,
flip from pass to fail**. The top ipTM decile still flips at 74%. Stored BG+BIF
`sc_RMSD` and `DockQ` values should not be used; use the `*_true` columns.

**Side-chain presence differs by stack.** RFdiffusion3, BG+LiMPNN and BoltzDesign1
design CIFs are backbone-only (N, CA, C, O); BG+BIF and HalluDesign are full atom. The
cause is that the external LigandMPNN stage never passes `--pack_side_chains`, which
defaults to off. Any metric that reads side chains is not comparable across stacks
without the retrospective packing pass.

**A comment in the HalluDesign production config is wrong, and is left as it ran.**
`configs/halludesign/random_init.yaml` carries the comment "deterministic -> the SAME
seed yields the SAME trajectory". It does not. The initial 80-mer is reproducible from
the seed, but the vendored inverse-folding wrapper reseeds the global random state from
OS entropy before sequence sampling, so the trajectory diverges from that point. The
comment has been left untouched because this archive records the configuration as it
ran; it should not be taken at face value.

**No hotspot or epitope conditioning was used anywhere.** `hotspots.json` holds an empty
array for all five stacks and is read by nothing. Published claims of hotspot residues
for this work are incorrect.

**RFdiffusion3 saw a 16-nt target; the other four saw 18 nt.** The RFD3 contig is
`B5-20,/0,80-80`, which drops the two terminal cytosines present in the target PDB.
This is the same asymmetry behind the residue-register bug above.

## Core packing

A bespoke monomer core-packing module was written (`build_core_manifests.py`,
`calibrate_core_packing.py`, `core_integrate.py`, `core_score_refs.py`) but its output
was **never merged** into either metric master, and it did not enter the filter cascade.
Neither master carries a `core_*`, `packing_density` or `core_hydrophobic_fraction`
column. Rosetta `packstat` is present and is a different quantity. The module should be
described as built and not applied.
