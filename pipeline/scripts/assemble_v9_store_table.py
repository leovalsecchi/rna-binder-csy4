#!/usr/bin/env python3
"""assemble_v9_store_table.py — build the single v9 CSV that build_filter_store.py consumes.

v9's metrics live in separate sidecars because each was produced by its own pipeline. This
joins the ones build_filter_store.py reads directly, under the exact column names it expects,
and adds `campaign` (v9 has no such column — the master merge normally supplies it).

The RNA-aligned rename is the load-bearing part:

    dockq_rna_<be>        -> dockq_self_<be>     (sc_DockQ)
    lrmsd_rna_corr_<be>   -> dockq_LRMSD_<be>    (sc_RMSD)

The master has used the RNA-aligned (target-fixed) convention since 2026-08-04. Feeding the
stock columns in under these names would put v9 ~1.45x off on LRMSD and -0.028 on DockQ
against v2-v8, invisibly — every cross-campaign comparison in the tool would be wrong.

foot / dsasa / pose_disagree / fold_disagree are NOT here: they are joined onto the store as
standalone .npy columns afterwards, the same way they were added for v2-v8.

pandas is fine here (1.5M rows x a few columns). It would NOT be fine on the 9.7 GB master.
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

BE = ("af3", "esm")

ap = argparse.ArgumentParser()
ap.add_argument("--clean", required=True, help="collected_clean_metrics.csv")
ap.add_argument("--dockq-rna", required=True, help="dockq_rna_v9_new.csv")
ap.add_argument("--pis-af3", required=True)
ap.add_argument("--pis-esm", required=True)
ap.add_argument("--hb-af3", required=True)
ap.add_argument("--hb-esm", required=True)
ap.add_argument("--campaign", default="v9")
ap.add_argument("--out", required=True)
a = ap.parse_args()

BASE_COLS = ["design_id", "uuid8", "stack", "gen", "invfold", "length", "sequence",
             "af3_msaoff_iptm", "esm_iptm",
             "af3_msaoff_cif_path", "esm_cif_path", "designed_cif_path"]

df = pd.read_csv(a.clean, usecols=BASE_COLS, low_memory=False)
n0 = len(df)
print(f"clean: {n0:,} rows", flush=True)
df["campaign"] = a.campaign


def join(df, path, usecols, rename, label):
    """Left-join a sidecar, asserting it neither loses nor duplicates rows."""
    s = pd.read_csv(path, usecols=usecols, low_memory=False).drop_duplicates("design_id")
    before = len(df)
    df = df.merge(s.rename(columns=rename), on="design_id", how="left")
    if len(df) != before:
        sys.exit(f"FATAL: {label} merge changed row count {before} -> {len(df)}")
    got = rename.get(usecols[1], usecols[1])
    print(f"  {label}: {int(df[got].notna().sum()):,} / {len(df):,} filled", flush=True)
    return df


# RNA-aligned DockQ + its LRMSD, renamed onto the master's column names
df = join(df, a.dockq_rna,
          ["design_id"] + [f"dockq_rna_{b}" for b in BE] + [f"lrmsd_rna_corr_{b}" for b in BE],
          {**{f"dockq_rna_{b}": f"dockq_self_{b}" for b in BE},
           **{f"lrmsd_rna_corr_{b}": f"dockq_LRMSD_{b}" for b in BE}}, "dockq_rna")

df = join(df, a.pis_af3, ["design_id", "pi_stack_score_af3"], {}, "pistack af3")
df = join(df, a.pis_esm, ["design_id", "pi_stack_score_esm"], {}, "pistack esm")
# hbonds_pyrosetta.py writes the BARE column name; the backend is in the filename only
df = join(df, a.hb_af3, ["design_id", "n_h_bonds_interface"],
          {"n_h_bonds_interface": "n_h_bonds_interface_af3"}, "hbond af3")
df = join(df, a.hb_esm, ["design_id", "n_h_bonds_interface"],
          {"n_h_bonds_interface": "n_h_bonds_interface_esm"}, "hbond esm")

need = (["design_id", "uuid8", "stack", "campaign", "gen", "invfold", "length", "sequence",
         "af3_msaoff_iptm", "esm_iptm",
         "af3_msaoff_cif_path", "esm_cif_path", "designed_cif_path"]
        + [f"dockq_self_{b}" for b in BE] + [f"dockq_LRMSD_{b}" for b in BE]
        + [f"pi_stack_score_{b}" for b in BE] + [f"n_h_bonds_interface_{b}" for b in BE])
missing = [c for c in need if c not in df.columns]
if missing:
    sys.exit(f"FATAL: assembled table lacks {missing}")
if len(df) != n0:
    sys.exit(f"FATAL: row count drifted {n0} -> {len(df)}")

df[need].to_csv(a.out, index=False)
print(f"\nwrote {len(df):,} rows x {len(need)} cols -> {a.out}", flush=True)
for c in need[13:]:
    print(f"  {c:28s} {int(df[c].notna().sum()):>9,} non-null", flush=True)
