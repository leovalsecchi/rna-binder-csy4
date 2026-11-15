# RNA Binder Pipeline — File Audit: `csy4_crrna_free` 5-stack campaign

_Audited 2026-06-05. Covers every file required to run the full 5-stack design +
refold + analysis campaign for `csy4_crrna_free`._

**Base path** (all paths below are relative to this unless noted):
`/cluster/work/wicky/lvalsecchi/rna_binder/` → locally `/Volumes/eulermount/lvalsecchi/rna_binder/`

**Note on the mount:** `/cluster/work/wicky` *is* the sshfs mount root `/Volumes/eulermount`
(so a path like `…/wicky/containers/x.sif` is `/Volumes/eulermount/containers/x.sif`).

**The "5 stacks"** = `rfd3__ligandmpnn`, `halludesign__ligandmpnn`, `boltzgen__boltzif`,
`boltzgen__ligandmpnn`, `boltzdesign__ligandmpnn`. The GUI defines **4** generator stacks
(`stacks.py`), but **BoltzGen multif emits 2 stacks** (BoltzIF + LiMPNN) from one job → 5 total.

**Status legend:** ✅ exists & correctly wired · ⚠️ exists, with a caveat · ❌ missing/needs attention

---

## 1. Launch entrypoint

| File | Purpose | Status |
|---|---|---|
| `code/scripts/launch.py` | The orchestrator. Loads a config YAML, applies `--override`, builds the sbatch per generator/refolder, submits via `sbatch`, writes `run_meta.json`. Dispatches on `config_kind: multif`/`refold_af3`/`refold_b2`. `--auto-refold` queues AF3 (×2: msaon+msaoff) + B2 with an `afterok` dependency; `--num-runs N` fans out shared-designs runs. | ✅ (49 KB; many `.bak_*` siblings present — harmless) |

---

## 2. Generators + inverse-fold workers

| File | Purpose | Status |
|---|---|---|
| `code/scripts/stages/rfd3_multif.sh` | RFD3 diffusion (foundry.sif) → gunzip + CIF→PDB → external LigandMPNN per backbone. Called by `build_rfd3_multif_sbatch`. | ✅ |
| `code/scripts/stages/boltzgen_multif.sh` | BoltzGen design once → BIF + LigandMPNN on the **same** backbone (the 2-stack job). Called by `build_boltzgen_multif_sbatch`. | ✅ |
| `code/scripts/stages/boltzdesign_multif.sh` | BD1 hallucination-only → collect survivors → external LigandMPNN. Has a `results_final/` CIF-fallback when `i_ptm_cutoff=0` skips `pdb/`. Called by `build_boltzdesign_multif_sbatch`. | ✅ |
| **HalluDesign — what runs HD** | **No stage worker exists.** HD is built **inline** by `build_halludesign_sbatch` in `launch.py`: `apptainer run` of `python ./HalluDesign_run.py` (cross_model) in `/cluster/work/wicky/lvalsecchi/halludesign/`, then calls `postprocess_halludesign.sh`. | ✅ (`HalluDesign_run.py` + `CHANGES.md` exist) |
| `code/scripts/stages/postprocess_halludesign.sh` | Flattens HD's nested `recycle_N/.../ligand_mpnn_S_K.pdb` into canonical `…_recycle<N>_sample<S>_pack<K>.pdb`. | ✅ |
| `code/scripts/stages/postprocess_boltzdesign.sh` | Postprocess for the **OLD bundled** BD1 path (`build_boltzdesign_sbatch`). | ⚠️ Exists but **not on the GUI path** — the GUI BD1 stack points at the multif config, so this script + `build_boltzdesign_sbatch` are dead for GUI launches. |

The standalone `code/scripts/halludesign_crossmodel_*.sbatch` templates are the pre-`launch.py`
manual approach — superseded by the inline builder; not used by the GUI.

---

## 3. Refold workers

| File | Purpose | Status |
|---|---|---|
| `code/scripts/stages/refold_af3.sh` | Reads `<source_run>/designs/`, extracts protein seq (via boltzgen.sif Python), builds per-design AF3 JSON (RNA chain B = `ACUGCCGUAUAGGCAG`, optional a3m MSA), loops `run_alphafold.py` (`--run_data_pipeline=false`), renames to `…_AF3_seed<S>_sample<N>.cif` + confidences. | ✅ |
| `code/scripts/stages/refold_b2.sh` | Same intake → per-design Boltz-2 YAML (protein `msa=empty`; RNA optional a3m + optional template CIF), `boltz predict`, renames to `…_B2_sample<N>.cif` + confidence JSON. | ✅ |

---

## 4. Analysis

| File | Purpose | Status |
|---|---|---|
| `code/scripts/aggregate_metrics.py` | CLI → `rna_binder.analysis.aggregate{_metrics,_to_designs}`; writes wide refold-metrics CSV. Portable (no hardcoded paths). | ✅ |
| `code/scripts/hbonds_pyrosetta.py` | Backfills `n_h_bonds_interface` with PyRosetta (`pyrosetta.sif`); reads CIF paths from CSV `abs_path` col. | ✅ (pyrosetta.sif exists) |
| `code/src/rna_binder/analysis/parsing.py` | Regex filename → `ParsedName` (gen/invfold/refolder/indices). | ✅ |
| `code/src/rna_binder/analysis/structural.py` | Biotite metrics: contacts, BSA, SSE, Rg, **`ca_rmsd_vs_designed`** (the "rmsd" column), binder seq. `n_h_bonds_interface` is intentionally `None` here (filled by the PyRosetta pass). | ⚠️ Refs `run_structural.sbatch` in a comment — that file is **missing** (comment-only; no functional impact). |
| `code/src/rna_binder/analysis/aggregate.py` (+ `__init__.py`, `refolder_metrics.py`) | Core aggregator. `_iter_refold_samples` only accepts plain `.cif` (gzipped refolds skipped — fine, refolds are plain). | ✅ |
| `analyses/Plots/run_all_plots.sh` | Copies the structural CSV into a per-run dir and runs every `plot_*.py` via boltzgen.sif. Hardcoded Euler paths + container. | ✅ |
| `analyses/Plots/plot_gen_pythons/*.py` (9) | `plot_af3_vs_b2_correlation`, `plot_filter_funnel`, `plot_filter_funnel_stacked`, `plot_hd_rts_comparison`, `plot_iptm_distribution`, `plot_iptm_vs_rmsd`, `plot_msaon_msaoff_correlation`, `plot_rmsd_distribution`, `plot_top_designs_table`. **All 9 are wired into `run_all_plots.sh`** (glob `plot_*.py`); no orphans, none missing. | ✅ (see caveats below) |

**Plot caveats:** `plot_hd_rts_comparison.py` has a **hardcoded uuid8→rts map**
(`328a184e`/`eee4909d`/`fb5198e2`) — dataset-specific, silently drops other rows. All scripts
share a hardcoded msaon-label allowlist `['msaon','prod_bg200_bif','prod_hd200_v2']`; any new
refold label silently falls into `AF3_msaoff`. Funnel/table filters on `n_h_bonds_interface`
pass nothing unless `hbonds_pyrosetta.py` ran first.

**Bonus (not in original list, found via `analyze_campaign.sh`):**
`analyses/Large_Gen_App_v1/full_analysis.py` ✅ exists.

---

## 5. Generator configs

| File | Purpose | Status |
|---|---|---|
| `code/configs/rfd3__multif/default.yaml` | `config_kind: multif`, rfd3 + [ligandmpnn]. Points `inputs`/`cwd` at the **sibling** `rfdiffusion3/examples/rna/`. `n_batches=1` (→8 backbones), LiMPNN `batch_size=10`. | ✅ |
| `code/configs/boltzgen__multif/default.yaml` | boltzgen + [boltzif, ligandmpnn]. spec len=80, target `csy4_rna_only.pdb`, bg `num_designs=1`, BIF 10 seqs, LiMPNN 10. | ✅ |
| `code/configs/boltzdesign__multif/default.yaml` | boltzdesign + [ligandmpnn], hallucination-only (`run_ligandmpnn/alphafold=False`). `design_samples=200`, len 80, target `csy4_rna_only.pdb`. | ✅ |
| `code/configs/halludesign/default.yaml` | cross_model. Container `rna_binder/containers/af3_protenix_esm.sif`; input `csy4_helix.pdb`; templates `template_protenix.json` + `template_af3.json`; RNA `ACUGCCGUAUAGGCAG`. | ✅ (all refs exist) |

---

## 6. Refolder configs

| File | Purpose | Status |
|---|---|---|
| `code/configs/refold__af3/default.yaml` | `alphafold3.sif` + `af3_params`. `rna_msa_path` → `targets/.../msa/csy4_crRNA_16nt_msa.a3m`. (launch.py runs **both** msaon + msaoff regardless of the `use_rna_msa: false` default.) | ✅ |
| `code/configs/refold__b2/default.yaml` | `boltz2.sif` + `params/.boltz`. `diffusion_samples=5`, mmcif, optional a3m/template (both default off). | ✅ |

---

## 7. Target files

| File | Purpose | Status |
|---|---|---|
| `targets/csy4_crrna_free/csy4_rna_only.pdb` | RNA-only target for BoltzGen / BD1 / RFD3-input. | ✅ |
| `targets/csy4_crrna_free/csy4_helix.pdb` | HD `input_file`. | ✅ |
| `targets/csy4_crrna_free/template_af3.json` | HD `extra_json_path` (AF3 dialect; 80-A poly-A + RNA `CACUGCCGUAUAGGCAGC`). | ✅ |
| `targets/csy4_crrna_free/template_protenix.json` | HD `template_path` (Protenix dialect). | ✅ |
| `targets/csy4_crrna_free/csy4_crRNA_16nt_msa.a3m` | RNA MSA for msaon refolds. | ⚠️ **Path correction:** actual location is `targets/csy4_crrna_free/msa/csy4_crRNA_16nt_msa.a3m` (in `msa/`, not the target root). The refold configs already point at the correct `msa/…` path. |
| `rfdiffusion3/examples/rna/csy4_crrna.json` | RFD3 contig input (`B5-20,/0,80-80`, length 96). | ⚠️ **Location correction:** lives in the **sibling repo** `/cluster/work/wicky/lvalsecchi/rfdiffusion3/…`, **not** under `rna_binder/`. The rfd3__multif config references it by absolute path (correct). Its `inputs/csy4_rna_only.pdb` also exists. |

---

## External assets referenced (all verified present)

**Containers** (`/cluster/work/wicky/containers/`): `foundry.sif`, `LigandMPNN.sif`,
`boltzgen.sif`, `boltz_design1.sif`, `alphafold3.sif`, `boltz2.sif`, `pyrosetta.sif` — all ✅.
Plus project-local `rna_binder/containers/af3_protenix_esm.sif` (HD) ✅.

**Params / weights** (`/cluster/work/wicky/params/`): `foundry/ligandmpnn_v_32_010_25.pt`,
`af3_params/`, `.boltz/`, `boltz/boltz1_conf.ckpt`, `boltz/ccd.pkl` — all ✅.

---

## How it all connects to the GUI (`~/Documents/rna_binder_gui/`)

- **`stacks.py`** — declares the 4 stack specs and which `default.yaml` each uses (the override
  paths the panel exposes). Source of truth for stack→config wiring.
- **`app.py`** — sidebar `remote_cwd` defaults to `/cluster/work/wicky/lvalsecchi/rna_binder`;
  builds `(label, config, overrides)` per selected stack and calls `euler.sbatch_via_launcher(...)`.
- **`euler.py`** — `build_launch_cmd` assembles
  `module load … && cd <remote_cwd> && python3 code/scripts/launch.py <config> --label … [--override …] [--auto-refold --refold-af3-override …] [--num-runs N]`
  and runs it over `ssh euler`. Reads (`squeue`, `refold_counts`) go over SSH.
- **`inspect_runs.py`** — reads run state directly off the eulermount
  (`EULERMOUNT_ROOT_DEFAULT = /Volumes/eulermount/lvalsecchi/rna_binder`), with a local disk cache.

**Flow:** GUI → SSH → `launch.py` → per-stage `.sh` workers → SLURM, while the GUI's
"Recent runs" panel reads the campaign tree through the mount.

---

## Bottom line

**The campaign is fully runnable — every container, param, weight, config, stage script, target,
and analysis file referenced by the 5-stack design + refold + analysis path exists and resolves.**
No file needs to be *created* to run it. Items worth noting:

1. **Two original-task paths were slightly off** (both files exist at the real location): the a3m
   is in `…/msa/`, and `csy4_crrna.json` is in the sibling `rfdiffusion3/` repo. Configs already
   use the correct paths.
2. **`postprocess_boltzdesign.sh` + `build_boltzdesign_sbatch`** are dead code on the GUI path
   (GUI uses BD1 multif).
3. **HD is the one generator with no `stages/*.sh` worker** — it's built inline in `launch.py`.
4. **`run_structural.sbatch`** (referenced in a `structural.py` comment) is missing — comment-only,
   no functional impact.
5. Analysis caveats: `n_h_bonds_interface` requires the PyRosetta pass before any funnel/table
   filter on it is meaningful; `plot_hd_rts_comparison.py` is pinned to specific uuid8s;
   `make_cif_view.py` still globs `designs/*.cif.gz` only (multif designs are `.pdb`, so its
   `designs/` view comes out empty).
