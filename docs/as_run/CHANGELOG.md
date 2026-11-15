# rna_binder — Changelog

Project-level history of significant changes and milestones.

---

## 2026-05-18 — Cross_model breakthrough

**Status:** First end-to-end success of HalluDesign cross_model (AF3 + Protenix alternating)
on a Csy4 crRNA target. Job 66952183, 19:19 elapsed.

### Source patches to upstream HalluDesign
See `/cluster/work/wicky/lvalsecchi/halludesign/CHANGES.md` for full diffs.

1. **evaluation.py:935** — fixed `None.item()` syntax bug in `process_confidence_metrics_protenix`
2. **models_utility.py:688** — `cross_model_op_protenix_eval` now uses `source_file` instead of MPNN's chain-A-stripped `packed_path`
3. **models_utility.py:893** — corrected `cif_path` construction in AF3-subprocess branch to match AF3's actual on-disk layout (`file_name/tag/tag_model.cif`, not `tag/tag/tag_model.cif`)

### Sbatch and templates
- New required arg: `--extra_json_path` pointing to AF3-format JSON template (alongside existing `--template_path` for Protenix). AF3-subprocess branch will crash without it.
- AF3 template must include `unpairedMsa: ""` on every chain (including RNA), else AF3 silently skips with "RNA chain X is missing unpaired MSA".
- Working template files: `targets/csy4_crrna_free/template_{af3,protenix}.json`

### Verified output
- Cycle 0 (Protenix op): eval_plddt 73.8, eval_iptm 0.95, op_iptm 0.54, mpnn_seq=AAAA (init)
- Cycle 1 (AF3 subprocess op): eval_plddt 58.8, eval_iptm 0.91, op_iptm 0.63, mpnn_seq=EPVPDVVERKEGIMKKDGKDYIISL…
- All CSV/CIF/PDB files written; no silent errors

---

## 2026-05-17 — Cross_model infrastructure stack

**Status:** Cross_model produces valid Protenix output (no longer crashes on CUDA/CUTLASS).
Cycle-2 chain propagation still broken (fixed next day).

10-layer infrastructure resolved:
- CUDA 13.0.2 bound from `/cluster/software/manual/cuda/13.0.2/x86_64` → `/usr/local/cuda-real`
- DS_SKIP_CUDA_CHECK=1 (DeepSpeed v0.17.5 doesn't recognize CUDA 13)
- CUTLASS 4.5.0 shallow clone at `containers/cutlass`, symlinked at `halludesign/cutlass`
- TMPDIR=`containers/runtime_home/tmp` (nvcc tmpxft needs writable dir)
- Apptainer `--home containers/runtime_home/` for kernel build cache

See memories 16-20 for full env var list.

---

## 2026-05-15 — Repo skeleton finalized; HalluDesign legacy archived

- Layout: `targets/`, `campaigns/`, `containers/`, `code/`, `_index/`, `archive/`
- Pre-skeleton HalluDesign runs (outputs 1-33) archived at `archive/csy4_halludesign_pre_skeleton_2026-05/`
- Targets defined: `csy4_crrna_bound`, `csy4_crrna_free`
- Tier A campaign stacks identified: HalluDesign+LigandMPNN (RNA wet-lab precedent), BoltzGen+BoltzIF, BoltzDesign+LigandMPNN, RFD3+LigandMPNN
