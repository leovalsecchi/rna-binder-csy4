#!/usr/bin/env python3
"""Join the af3/b2/esm groove-dSASA sidecars onto the master and build the worst-case
3-backend consensus, following interface_consensus_integrate.py: the bare column BECOMES
the consensus (min for higher-is-better, max for lower-is-better), per-backend kept as
<col>_af3 / <col>_b2 / <col>_esm.

    python3 groove_consensus_integrate.py <master.csv> <sidecar_dir> [--out F] [--dry-run]

Differences from interface_consensus_integrate.py, on purpose:
  * AF3 is a SIDECAR here, not already-bare in the master (the master has no dsasa_rna_*).
  * Writes to a NEW file (--out, default <master>_groove.csv). The master is not modified.
  * Only the SASA/groove columns are taken from the sidecars. The sidecars also carry
    recomputed CONTACT columns, which are deliberately IGNORED -- the master's contact
    columns are already Rosetta-consensused and must not be clobbered by this pass.

NOTE ON UNITS: dsasa_rna_* is the RNA HALF of the interface, from Biotite. The master's
dsasa_total is Rosetta and counts BOTH sides; the two do not reconcile (measured r=0.80,
biotite/rosetta=1.26 on n=40). Do not mix the two families in one ratio.
"""
import argparse
import os
import sys

import pandas as pd

# higher = more buried / more base reading = better -> worst case is the MIN
MIN_COLS = [
    "dsasa_total_biotite", "dsasa_apolar_biotite", "dsasa_polar_biotite",
    "dsasa_apolar_frac_biotite",
    "dsasa_rna_total", "dsasa_rna_backbone", "dsasa_rna_phosphate", "dsasa_rna_sugar",
    "dsasa_rna_o2prime", "dsasa_rna_base", "dsasa_rna_base_major",
    "dsasa_rna_base_minor", "dsasa_rna_base_wc",
    "frac_dsasa_base", "frac_dsasa_major",
]
# diagnostic: base atoms that got no edge assignment (modified bases). lower = better
MAX_COLS = ["dsasa_rna_base_unassigned"]
MANAGED = MIN_COLS + MAX_COLS
SENTINEL = "dsasa_rna_total"     # presence => that backend produced a result
BACKENDS = ("af3", "b2", "esm")


def load_sidecar(d, be):
    path = os.path.join(d, f"groove_{be}.csv")
    if not os.path.exists(path):
        sys.exit(f"missing sidecar {path}")
    s = pd.read_csv(path, low_memory=False).drop_duplicates("design_id")
    have = [c for c in MANAGED if c in s.columns]
    if not have:
        sys.exit(f"{path}: none of the managed columns present")
    s = s[["design_id"] + have].rename(columns={c: f"{c}_{be}" for c in have})
    n_ok = s[f"{SENTINEL}_{be}"].notna().sum() if f"{SENTINEL}_{be}" in s else 0
    print(f"{be}: {len(s)} rows, {len(have)} managed cols, {n_ok} with {SENTINEL}")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("sidecar_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    out = args.out or args.master.replace(".csv", "_groove.csv")

    m = pd.read_csv(args.master, low_memory=False)
    print(f"master: {len(m)} rows, {len(m.columns)} cols")
    drop = [f"{c}_{s}" for c in MANAGED for s in BACKENDS if f"{c}_{s}" in m.columns]
    drop += [c for c in MANAGED if c in m.columns]
    if drop:
        print(f"dropping {len(drop)} pre-existing groove cols (idempotent re-run)")
        m = m.drop(columns=drop)

    for be in BACKENDS:
        m = m.merge(load_sidecar(args.sidecar_dir, be), on="design_id", how="left")

    for col in MANAGED:
        trio = [f"{col}_{be}" for be in BACKENDS]
        for t in trio:
            if t not in m.columns:
                m[t] = pd.NA
        vals = m[trio].apply(pd.to_numeric, errors="coerce")
        m[col] = vals.min(axis=1) if col in MIN_COLS else vals.max(axis=1)

    sent = [f"{SENTINEL}_{be}" for be in BACKENDS]
    m["groove_n_structures"] = m[sent].apply(pd.to_numeric, errors="coerce").notna().sum(axis=1)

    print("\n--- per-backend medians (sanity: should be close across backends) ---")
    for col in ["dsasa_rna_total", "dsasa_rna_base", "frac_dsasa_base", "frac_dsasa_major"]:
        med = {be: pd.to_numeric(m[f"{col}_{be}"], errors="coerce").median() for be in BACKENDS}
        print(f"  {col:22s} " + "  ".join(f"{be}={v:.4g}" for be, v in med.items())
              + f"   consensus={pd.to_numeric(m[col], errors='coerce').median():.4g}")
    print("\ngroove_n_structures:\n"
          + m["groove_n_structures"].value_counts().sort_index().to_string())

    if args.dry_run:
        print("\n--dry-run: not writing")
        return
    m.to_csv(out, index=False)
    print(f"\nwrote {out}  ({len(m)} rows, {len(m.columns)} cols)")


if __name__ == "__main__":
    main()
