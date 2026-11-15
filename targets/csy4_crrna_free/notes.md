# csy4_crrna_free

18-nt crRNA from PDB 4AL5 (Haurwitz et al. 2010), refolded in isolation (no Csy4 protein).

## Target

- **Sequence (5'→3'):** `CACUGCCGUAUAGGCAGC` (18 nt)
- **Source:** trimmed from HalluDesign output27/recycle_8/csy4_complex_recycle_8.pdb
- **Reference structure:** RFdiffusion-AA-derived helix80 + this RNA → csy4_complex.pdb

## Files

- `csy4_complex.pdb` — chain A: 80 CA-only Ala helix (placeholder for `random_init=True`); chain B: 18-nt RNA
- `template_protenix.json` — Protenix list-format template (used by --template_path)
- `template_af3.json` — AlphaFold3 dict-format template (used by --extra_json_path for cross_model)
- `hotspots.json` — RNA residues to fix during MPNN design (B4–B21)

## Note on CA-only chain A

The 80 residues of chain A are CA-only placeholder atoms (no full backbone). With `--random_init`, Protenix replaces them with Gaussian-noise diffusion init on cycle 0, so the input coordinates don't act as a structural prior. Paper-faithful equivalent of "random sequence + Gaussian noise coordinates" (HalluDesign paper line 209-211).

## Working campaigns

- `halludesign__ligandmpnn` (cross_model mode): proven working 2026-05-18, job 66952183
  - Cycle 0: Protenix op, eval_iptm 0.95
  - Cycle 1: AF3-subprocess op, eval_iptm 0.91, first real designed sequence
  - See campaign folder for full CSV
