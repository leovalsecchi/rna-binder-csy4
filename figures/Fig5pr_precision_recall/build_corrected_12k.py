#!/usr/bin/env python
"""Build the corrected 12k library metric table.

Why: inputs/csy4_binder_library_12k_v6_k30.csv stores PRE-rescore sc_rmsd_* / dockq_* for the
2,400 BG+BIF designs (mean ligand-RMSD ~6.0 stored vs ~17.5 true; DockQ ~0.62 vs ~0.36). This
script substitutes the true-backbone values and merges the recovered interface metrics.

Inputs (all in inputs/):
  csy4_binder_library_12k_v6_k30.csv        12,000 rows x 22 cols (never modified)
  large/bif_true_per_backend.csv.gz         BG+BIF true-backbone rmsd/dockq per backend
  fig5pr_recovered_metrics_12k.csv.gz       heavy contacts, sc_value, packstat, BUNS, pi-stack,
                                            dockq_self, H-bonds, dsasa_total; _af3/_esm/_b2
Outputs:
  inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv   (NAS)
  ~/rna_binder_plots/_data/csy4_binder_library_12k_v6_k30_CORRECTED.csv  (+ README.md there)

Adds bool `bif_rescored` (True = BG+BIF row with a true value). BG+BIF rows with no rescore
get NaN in the four substituted columns. Recovered `stack` is dropped (library `stack` kept).

Run:  ~/.venvs/binder_map/bin/python code/Fig5pr_precision_recall/build_corrected_12k.py
"""
from pathlib import Path
import shutil
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
IN = REPO / "inputs"
NAME = "csy4_binder_library_12k_v6_k30_CORRECTED.csv"
LOCAL = Path.home() / "rna_binder_plots" / "_data"

src = pd.read_csv(IN / "csy4_binder_library_12k_v6_k30.csv")
true = pd.read_csv(IN / "large" / "bif_true_per_backend.csv.gz")
rec = pd.read_csv(IN / "fig5pr_recovered_metrics_12k.csv.gz")
assert true.design_id.is_unique and rec.design_id.is_unique and src.design_id.is_unique

MAP = {"sc_rmsd_af3": "rmsd_af3_true", "sc_rmsd_esm": "rmsd_esm_true",
       "dockq_af3": "dockq_af3_true", "dockq_esm": "dockq_esm_true"}
df = src.merge(true, on="design_id", how="left", validate="1:1")
bif = df["stack"] == "BG+BIF"
df["bif_rescored"] = bif & df[list(MAP.values())].notna().all(axis=1)
before = {c: src.loc[src["stack"] == "BG+BIF", c].mean() for c in MAP}
for c, t in MAP.items():
    df.loc[bif, c] = df.loc[bif, t]          # NaN where no rescore
df = df.drop(columns=list(MAP.values()))
after = {c: df.loc[bif, c].mean() for c in MAP}

df = df.merge(rec.drop(columns=["stack"]), on="design_id", how="left", validate="1:1")

# column order: identity, then per-metric groups (af3, esm, b2 together)
ident = ["n", "design_id", "sequence", "length", "stack", "generator", "inverse_folder",
         "campaign", "backbone_id", "fold_cluster", "epitope_cluster", "bif_rescored"]
metrics = ["iptm", "sc_rmsd", "dockq", "dsasa", "n_hbonds", "sc_value", "packstat",
           "n_heavy_contacts_4A", "n_heavy_contacts_4p5A", "delta_unsat_hbonds",
           "delta_unsat_hbonds_opth", "buns_vsasa", "buns_vsasa_opth", "pi_stack_score",
           "dockq_self", "n_h_bonds_interface", "dsasa_total"]
order = list(ident)
for m in metrics:
    order += [c for c in (f"{m}_af3", f"{m}_esm", f"{m}_b2", f"{m}_max") if c in df.columns]
order += [c for c in df.columns if c not in order]
df = df[order]

# ---- verification
print("rows", len(df), "| unique ids", df.design_id.is_unique, "| same ids as source",
      set(df.design_id) == set(src.design_id), "| cols", df.shape[1])
for c in ("sc_rmsd_af3", "sc_rmsd_esm", "dockq_af3", "dockq_esm"):
    print(f"BG+BIF mean {c}: {before[c]:.3f} -> {after[c]:.3f}")
print("bif_rescored True:", int(df.bif_rescored.sum()), "of", int(bif.sum()),
      "| NaN (no rescore):", int((bif & ~df.bif_rescored).sum()))
assert 5 < before["sc_rmsd_af3"] < 7 and 16 < after["sc_rmsd_af3"] < 19, "ligand-RMSD did not move as expected"
assert 0.55 < before["dockq_af3"] < 0.7 and 0.3 < after["dockq_af3"] < 0.42, "DockQ did not move as expected"
a = src.set_index("design_id").loc[df.design_id[~bif].values]
b = df.set_index("design_id").loc[a.index]
for c in MAP:
    assert a[c].equals(b[c]), c
print("non-BG+BIF rows identical to source for", list(MAP), ":", int((~bif).sum()), "rows")
print("coverage of merged metrics (non-null / 12000):")
merged = [c for c in rec.columns if c not in ("design_id", "stack")]
cov = df[merged].notna().sum()
print(cov.to_string())

# ---- save
df.to_csv(IN / NAME, index=False)
LOCAL.mkdir(parents=True, exist_ok=True)
shutil.copyfile(IN / NAME, LOCAL / NAME)
print("saved", IN / NAME, "and", LOCAL / NAME)
