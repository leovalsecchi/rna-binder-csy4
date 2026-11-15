#!/usr/bin/env python3
"""opth_integrate.py — merge the optH columns into the master tables.

Adds, per backend (af3/b2/esm):
    n_h_bonds_interface_opth_<be>   interface H-bonds after core::pack::optimizeH
    delta_unsat_hbonds_opth_<be>    IA buried-unsat polars after optimizeH
    buns_vsasa_<be>                 modern-style BUNS (VSASA burial, probe 1.1)
    buns_vsasa_opth_<be>            same, after optimizeH
    n_polar_h_moved_<be>            polar H that moved >0.1 A

plus bare consensus columns following consensus_integrate.py's convention
(n_h_bonds_* = row-wise MIN = worst case; every unsat flavour = row-wise MAX).

The *_raw columns are NOT merged: they are a verbatim recomputation of the
existing production columns. Instead this script ASSERTS that identity on every
overlapping row and reports any mismatch — a non-zero count means a bad join or
a stale master, and the run should be treated as suspect.

Writes <master>.tmp, verifies the row count, then moves the original aside to
<master>.bak_pre_opth and the tmp into place.

    python3 opth_integrate.py <master.csv> [--chunksize 250000]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

CC = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
BES = ["af3", "b2", "esm"]
SRC = {"hb_opth": "n_h_bonds_interface_opth", "buns_ia_opth": "delta_unsat_hbonds_opth",
       "buns_vsasa_raw": "buns_vsasa", "buns_vsasa_opth": "buns_vsasa_opth",
       "n_polar_h_moved": "n_polar_h_moved"}
CHECK = {"hb_raw": "n_h_bonds_interface", "buns_ia_raw": "delta_unsat_hbonds"}
MINCOLS = ["n_h_bonds_interface_opth"]
MAXCOLS = ["delta_unsat_hbonds_opth", "buns_vsasa", "buns_vsasa_opth"]


def load_lookup():
    use = ["design_id"] + list(SRC) + list(CHECK)
    out = None
    for be in BES:
        d = pd.read_csv(CC + f"opth_full_{be}.csv", usecols=use, low_memory=False)
        d = d.drop_duplicates("design_id")
        ren = {k: f"{v}_{be}" for k, v in SRC.items()}
        ren.update({k: f"__chk_{v}_{be}" for k, v in CHECK.items()})
        d = d.rename(columns=ren)
        print("  %s: %d rows" % (be, len(d)), flush=True)
        out = d if out is None else out.merge(d, on="design_id", how="outer")
    for c in MINCOLS:
        out[c] = out[[f"{c}_{b}" for b in BES]].min(axis=1)
    for c in MAXCOLS:
        out[c] = out[[f"{c}_{b}" for b in BES]].max(axis=1)
    print("  lookup: %d unique design_id, %d cols" % (len(out), out.shape[1]), flush=True)
    return out.set_index("design_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("--chunksize", type=int, default=250_000)
    a = ap.parse_args()
    master = a.master if a.master.startswith("/") else CC + a.master
    if not os.path.exists(master):
        sys.exit("no such master: " + master)

    print("loading optH lookup...", flush=True)
    lk = load_lookup()
    newcols = [f"{v}_{b}" for v in SRC.values() for b in BES] + MINCOLS + MAXCOLS
    chkcols = [c for c in lk.columns if c.startswith("__chk_")]

    tmp = master + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    n_rows = n_hit = 0
    mism = {c: 0 for c in chkcols}
    cmp_n = {c: 0 for c in chkcols}
    first = True
    for ch in pd.read_csv(master, chunksize=a.chunksize, low_memory=False):
        dup = [c for c in newcols if c in ch.columns]
        if dup:
            ch = ch.drop(columns=dup)          # idempotent re-run
        j = ch.join(lk, on="design_id", how="left")
        # identity check against the production columns
        for c in chkcols:
            prod = c.replace("__chk_", "")
            if prod in j.columns:
                m = j[c].notna() & j[prod].notna()
                cmp_n[c] += int(m.sum())
                mism[c] += int((np.abs(j.loc[m, c] - j.loc[m, prod]) > 0.01).sum())
        j = j.drop(columns=chkcols)
        n_rows += len(j)
        n_hit += int(j[newcols[0]].notna().sum())
        j.to_csv(tmp, index=False, mode="w" if first else "a", header=first)
        first = False
        print("  %d rows..." % n_rows, flush=True)

    print("\nidentity check (_raw recomputation vs production columns):", flush=True)
    bad = False
    for c in chkcols:
        print("  %-42s compared %8d  mismatches %d"
              % (c.replace("__chk_", ""), cmp_n[c], mism[c]), flush=True)
        if mism[c]:
            bad = True
    print("\nrows: %d, with optH values: %d (%.2f%%)"
          % (n_rows, n_hit, 100 * n_hit / max(n_rows, 1)), flush=True)

    with open(tmp) as fh:
        tmp_rows = sum(1 for _ in fh) - 1
    if tmp_rows != n_rows:
        sys.exit("FATAL: tmp has %d rows, expected %d — leaving %s in place"
                 % (tmp_rows, n_rows, tmp))
    if bad:
        sys.exit("FATAL: identity check failed — NOT replacing the master. tmp at " + tmp)

    bak = master + ".bak_pre_opth"
    os.replace(master, bak)
    os.replace(tmp, master)
    print("OK: %s updated (+%d cols); original -> %s" % (master, len(newcols), bak), flush=True)


if __name__ == "__main__":
    main()
