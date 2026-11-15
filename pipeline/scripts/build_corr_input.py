#!/usr/bin/env python3
"""build_corr_input.py — assemble the metric table the correlation matrix is computed on.

The 100k app table still carries the STOCK-convention DockQ (`dockq_*_msa_off`), which was
superseded on 2026-08-04 by the RNA-aligned (target-fixed) values. It also predates the
refolder pair-RMSD and groove ΔSASA work. This joins those in and writes one compact CSV so
the plotting stays local and cheap.

Joins onto app_metrics_100k.csv by design_id:
  dockq_rna_v2v7.csv          RNA-aligned DockQ per backend + worst-case consensus
  pairrmsd_v2v7.csv.gz        refolder pose / fold agreement RMSDs
  app_metrics_100k_groove.csv RNA-side ΔSASA split (base vs backbone) on the AF3 refold
"""
import argparse
import os

import numpy as np
import pandas as pd

BE = ["af3", "b2", "esm"]

BASE = ["design_id", "stack", "consensus_iptm", "mean_iptm", "af3_msaoff_ptm",
        "dsasa_total", "dsasa_hphobic_frac", "n_interface_contacts", "sc_value",
        "packstat", "worst_clash", "n_h_bonds_interface", "delta_unsat_hbonds",
        "pi_stack_count_af3", "tip_weighted_score", "frac_base_contacts",
        "n_contact_res_prot", "n_contact_res_rna", "red_total_per_ifaceres",
        "core_packing_score", "core_underpacked_frac", "apo_rmsd_vs_esm",
        "ca_rmsd_vs_designed", "sse_helix_fraction", "seq_max_kmer",
        "dockq_af3_msa_off", "dockq_b2_msa_off", "dockq_esm_msa_off"]
PAIR = ["design_id"] + [f"rmsd_{k}_{p}" for k in ("pose", "prot", "rna")
                        for p in ("af3_b2", "af3_esm", "b2_esm")]
GROOVE = ["design_id", "dsasa_rna_total_af3", "dsasa_rna_base_af3",
          "dsasa_rna_backbone_af3", "dsasa_rna_base_major_af3"]

ap = argparse.ArgumentParser()
ap.add_argument("--app", required=True)
ap.add_argument("--dockq-rna", required=True)
ap.add_argument("--pairrmsd", required=True)
ap.add_argument("--groove")
ap.add_argument("--out", required=True)
a = ap.parse_args()


def read(path, want, **kw):
    head = pd.read_csv(path, nrows=0, **kw).columns.tolist()
    use = [c for c in want if c in head]
    miss = [c for c in want if c not in head]
    if miss:
        print(f"  {os.path.basename(path)}: MISSING {miss}", flush=True)
    return pd.read_csv(path, usecols=use, low_memory=False, **kw)


df = read(a.app, BASE)
print(f"app: {len(df):,} rows, {len(df.columns)} cols", flush=True)

dq = read(a.dockq_rna, ["design_id", "dockq_rna", "dockq_rna_consensus_n"] +
          [f"dockq_rna_{b}" for b in BE])
df = df.merge(dq, on="design_id", how="left")
print(f"  + dockq_rna: {int(df.dockq_rna.notna().sum()):,} designs have the RNA-aligned "
      f"consensus ({100*df.dockq_rna.notna().mean():.1f}%)", flush=True)

pr = read(a.pairrmsd, PAIR)
df = df.merge(pr, on="design_id", how="left")
# worst-case across the three refolder pairs, matching the project's consensus convention
for k in ("pose", "prot", "rna"):
    cols = [c for c in df.columns if c.startswith(f"rmsd_{k}_")]
    if cols:
        df[f"rmsd_{k}_worst"] = df[cols].max(axis=1)
print(f"  + pairrmsd: {int(df.rmsd_pose_worst.notna().sum()):,} with a pose RMSD", flush=True)

if a.groove and os.path.exists(a.groove):
    gr = read(a.groove, GROOVE)
    df = df.merge(gr, on="design_id", how="left")
    tot = df.get("dsasa_rna_total_af3")
    if tot is not None:
        df["dsasa_base_frac"] = df["dsasa_rna_base_af3"] / tot.where(tot > 0)
        print(f"  + groove: {int(df.dsasa_base_frac.notna().sum()):,} with a base-ΔSASA "
              f"fraction", flush=True)

# stock consensus DockQ, kept so the two conventions can be compared in the same matrix
stock = [f"dockq_{b}_msa_off" for b in BE if f"dockq_{b}_msa_off" in df.columns]
if len(stock) == 3:
    df["dockq_stock_consensus"] = df[stock].min(axis=1).where(df[stock].notna().sum(axis=1).eq(3))

df.to_csv(a.out, index=False)
print(f"\nwrote {a.out}: {len(df):,} rows x {len(df.columns)} cols", flush=True)
