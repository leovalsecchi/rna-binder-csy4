#!/usr/bin/env python3
"""Join B2 + ESM interface-suite sidecars onto the master and build the worst-case
3-backend consensus, following hbond/red_consensus_integrate: bare column BECOMES the
consensus (min for higher-is-better, max for lower-is-better), original AF3 preserved as
<col>_af3, per-backend as <col>_b2/_esm. SKIP metrics untouched. Backs up; idempotent.

    python3 interface_consensus_integrate.py <master.csv> <be_dir> [--dry-run]

<be_dir> holds iface_be_b2.csv and iface_be_esm.csv, each with 'design_id' + the
recomputed interface columns for that backend.
"""
import argparse
import os
import shutil
import sys

import pandas as pd

MIN_COLS = ["dsasa_total", "dsasa_polar", "dsasa_hphobic", "sc_value", "packstat",
            "n_heavy_contacts_4A", "n_heavy_contacts_4p5A", "n_contact_res_prot",
            "n_contact_res_rna", "n_contacts_base", "n_contacts_backbone",
            "n_contacts_phosphate", "n_contacts_sugar", "frac_base_contacts"]
MAX_COLS = ["delta_unsat_hbonds", "n_interface_clashes", "clash_overlap_sum", "worst_clash"]
MANAGED = MIN_COLS + MAX_COLS
SENTINEL = "sc_value"   # presence => that backend produced a Rosetta score for the design


def load_sidecar(be_dir, be):
    path = os.path.join(be_dir, f"iface_be_{be}.csv")
    if not os.path.exists(path):
        sys.exit(f"missing sidecar {path}")
    d = pd.read_csv(path, low_memory=False).drop_duplicates("design_id")
    have = [c for c in MANAGED if c in d.columns]
    d = d[["design_id"] + have].rename(columns={c: f"{c}_{be}" for c in have})
    print(f"{be}: {len(d)} rows, {len(have)} managed cols")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("be_dir")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    m = pd.read_csv(args.master, low_memory=False)
    drop = [f"{c}_{s}" for c in MANAGED for s in ("af3", "b2", "esm")
            if f"{c}_{s}" in m.columns]
    if drop:
        m = m.drop(columns=drop)
    m = (m.merge(load_sidecar(args.be_dir, "b2"), on="design_id", how="left")
          .merge(load_sidecar(args.be_dir, "esm"), on="design_id", how="left"))

    for col in MANAGED:
        if col not in m.columns:
            print(f"WARN: master missing managed col {col}; skipping")
            continue
        m[f"{col}_af3"] = m[col]
        trio = [f"{col}_af3", f"{col}_b2", f"{col}_esm"]
        for t in trio:
            if t not in m.columns:
                m[t] = pd.NA
        vals = m[trio].apply(pd.to_numeric, errors="coerce")
        m[col] = vals.min(axis=1) if col in MIN_COLS else vals.max(axis=1)

    sent = [f"{SENTINEL}_{s}" for s in ("af3", "b2", "esm") if f"{SENTINEL}_{s}" in m.columns]
    m["consensus_n_structures"] = (m[sent].apply(pd.to_numeric, errors="coerce")
                                   .notna().sum(axis=1)) if sent else 0

    if args.dry_run:
        for col in ("sc_value", "n_contacts_base", "worst_clash"):
            if f"{col}_af3" in m.columns:
                print(f"  {col}: af3_med={m[f'{col}_af3'].median():.4g} "
                      f"consensus_med={m[col].median():.4g}")
        print("consensus_n_structures:\n"
              f"{m['consensus_n_structures'].value_counts().sort_index()}")
        return

    shutil.copy(args.master, args.master + ".bak_pre_iface_consensus")
    m.to_csv(args.master, index=False)
    print(f"wrote {args.master}: {len(m)} rows x {len(m.columns)} cols")


if __name__ == "__main__":
    main()
