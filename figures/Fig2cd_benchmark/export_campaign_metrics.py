#!/usr/bin/env python3
"""Export the per-design campaign metrics as one CSV, figure-ready.

Consensus definitions verified against the per-stack means on the campaign
summary slide (all five stacks reproduce to the printed precision):
    iptm     = min(af3, esm)      "worst of the two refolds"
    sc_dockq = min(af3, esm)
    sc_rmsd  = max(af3, esm)
    rg       = mean(af3, esm)

BG+BIF EXCEPTION: its stored dockq/rmsd were measured against a Boltz
re-prediction, not the true BoltzGen backbone. The slide uses the corrected
values, and so does this file: BG+BIF sc_dockq/sc_rmsd come from
bif_true_bb_rescore/full/bif_true_bb_scores.csv.gz (dockq_self_*_true,
dockq_LRMSD_*_true). Store values give 0.096 / 29.77; corrected give
0.042 / 39.00, which is what the slide shows.
The raw per-backend columns are kept as-is so the correction stays visible.
"""
import numpy as np, pandas as pd, os, json

S = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/filter_store_v2v9_hdtop"
BIF = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_true_bb_rescore/full/bif_true_bb_scores.csv.gz"
OUT = "/cluster/work/wicky/lvalsecchi/campaign_metrics_full.csv.gz"

def col(n): return np.load(os.path.join(S, n + ".npy"), allow_pickle=True)
VOC = json.load(open(os.path.join(S, "meta.json")))["vocabs"]["stack"]

df = pd.DataFrame({
    "design_id": col("design_id").astype(str),
    "stack":     [VOC[i] for i in col("stack")],
    "iptm_af3":  col("iptm_af3"),   "iptm_esm":  col("iptm_esm"),
    "dockq_af3": col("dockq_af3"),  "dockq_esm": col("dockq_esm"),
    "rmsd_af3":  col("rmsd_af3"),   "rmsd_esm":  col("rmsd_esm"),
    "rog_af3":   col("rog_af3"),    "rog_esm":   col("rog_esm"),
})
print(f"store rows {len(df):,}")

df["iptm"]     = np.fmin(df.iptm_af3, df.iptm_esm)
df["sc_dockq"] = np.fmin(df.dockq_af3, df.dockq_esm)
df["sc_rmsd"]  = np.fmax(df.rmsd_af3, df.rmsd_esm)
df["rg"]       = (df.rog_af3 + df.rog_esm) / 2.0

# --- BG+BIF: swap in the true-backbone rescore -----------------------------
t = pd.read_csv(BIF, usecols=["design_id", "dockq_self_af3_true", "dockq_self_esm_true",
                              "dockq_LRMSD_af3_true", "dockq_LRMSD_esm_true"])
t["bif_dockq"] = np.fmin(t.dockq_self_af3_true, t.dockq_self_esm_true)
t["bif_rmsd"]  = np.fmax(t.dockq_LRMSD_af3_true, t.dockq_LRMSD_esm_true)
m = df["stack"].eq("BG+BIF")
before = df.loc[m, ["sc_dockq", "sc_rmsd"]].mean().to_dict()
df = df.merge(t[["design_id", "bif_dockq", "bif_rmsd"]], on="design_id", how="left")
m = df["stack"].eq("BG+BIF")
hit = m & df.bif_dockq.notna()
print(f"BG+BIF rows {int(m.sum()):,}; rescore matched {int(hit.sum()):,}")
df.loc[hit, "sc_dockq"] = df.loc[hit, "bif_dockq"]
df.loc[hit, "sc_rmsd"]  = df.loc[hit, "bif_rmsd"]
df["bif_rescored"] = hit.astype(int)
df = df.drop(columns=["bif_dockq", "bif_rmsd"])
print(f"BG+BIF sc_dockq {before['sc_dockq']:.3f} -> {df.loc[m,'sc_dockq'].mean():.3f}"
      f" | sc_rmsd {before['sc_rmsd']:.2f} -> {df.loc[m,'sc_rmsd'].mean():.2f}")

print("\nper-stack means (compare to slide):")
print(df.groupby("stack")[["iptm", "sc_dockq", "sc_rmsd", "rg"]].mean().round(3).to_string())
print("\nper-stack n:", df.groupby("stack").size().to_dict())

df.to_csv(OUT, index=False, float_format="%.4f")
print("\nwrote", OUT, f"{os.path.getsize(OUT)/1e6:.0f} MB")
