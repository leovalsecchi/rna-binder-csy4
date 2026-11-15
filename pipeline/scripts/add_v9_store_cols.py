#!/usr/bin/env python3
"""add_v9_store_cols.py — join the four standalone metrics onto the v9 filter store.

These four are not built by build_filter_store.py; for v2-v8 they were joined onto the store
afterwards as their own .npy columns, and v9 must match column-for-column or concat_stores.py
will NaN-fill one side.

    foot_af3/_esm    <- foot_v9_<be>.csv : foot_tp_min_dist
    dsasa_af3/_esm   <- iface_be_<be>.csv : dsasa_total
    pose_disagree    <- pairrmsd sidecars : rmsd_pose_af3_esm
    fold_disagree    <- pairrmsd sidecars : rmsd_prot_af3_esm

`foot` is the TRANSPLANT estimator (`_tp_`), not the A-form continuation (`_af_`) — the two
rank the same but transplant reads ~1.2 A tighter, so mixing them across campaigns would
silently shift every threshold. v2-v8 stores foot_tp_min_dist__<be>; the v9 sidecars are
per-backend files so the column is unsuffixed.

b2 is never read: the tool has no b2 column.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import numpy as np

csv.field_size_limit(10 ** 7)

ap = argparse.ArgumentParser()
ap.add_argument("--store", required=True)
ap.add_argument("--foot-af3", required=True)
ap.add_argument("--foot-esm", required=True)
ap.add_argument("--iface-af3", required=True)
ap.add_argument("--iface-esm", required=True)
ap.add_argument("--pairrmsd-dir", required=True)
ap.add_argument("--compare-store", default="", help="optional v2-v8 store, for a sanity diff")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

did = np.load(os.path.join(a.store, "design_id.npy"))
n = did.shape[0]
pos = {d: i for i, d in enumerate(did.tolist())}
print(f"store: {n:,} rows", flush=True)
del did

out: dict[str, np.ndarray] = {}
unmatched = 0


def fill(name: str, paths: list[str], src: str) -> None:
    """Read `src` from one or more design_id-keyed CSVs into a new store column."""
    global unmatched
    v = np.full(n, np.nan, dtype=np.float32)
    got = 0
    for p in paths:
        with open(p, newline="") as f:
            r = csv.DictReader(f)
            if src not in (r.fieldnames or []):
                sys.exit(f"FATAL: {p} has no column {src!r}; has {r.fieldnames}")
            for row in r:
                i = pos.get(row["design_id"].encode())
                if i is None:
                    unmatched += 1
                    continue
                s = row.get(src, "")
                if s not in ("", "NA", "nan", "NaN"):
                    try:
                        v[i] = float(s)
                        got += 1
                    except ValueError:
                        pass
    out[name] = v
    ok = np.isfinite(v)
    med = float(np.median(v[ok])) if ok.any() else float("nan")
    print(f"  {name:15s} <- {src:22s} {ok.sum():>9,}/{n:,} ({100*ok.sum()/n:5.2f}%)  "
          f"median {med:.3f}", flush=True)


fill("foot_af3", [a.foot_af3], "foot_tp_min_dist")
fill("foot_esm", [a.foot_esm], "foot_tp_min_dist")
fill("dsasa_af3", [a.iface_af3], "dsasa_total")
fill("dsasa_esm", [a.iface_esm], "dsasa_total")
pr = sorted(glob.glob(os.path.join(a.pairrmsd_dir, "pairrmsd_*.csv")))
print(f"  pairrmsd sidecars: {len(pr)}", flush=True)
fill("pose_disagree", pr, "rmsd_pose_af3_esm")
fill("fold_disagree", pr, "rmsd_prot_af3_esm")

print(f"\nunmatched source rows (not in store): {unmatched:,}")

if a.compare_store:
    print("\nsanity — v9 median vs the v2-v8 store (should be the same ballpark):")
    for k in out:
        p = os.path.join(a.compare_store, k + ".npy")
        if not os.path.exists(p):
            print(f"  {k:15s} (absent from compare store)")
            continue
        o = np.load(p)
        of = o[np.isfinite(o)]
        vf = out[k][np.isfinite(out[k])]
        print(f"  {k:15s} v9 {np.median(vf):9.3f}   v2-v8 {np.median(of):9.3f}")

if a.dry_run:
    print("\ndry run — nothing written")
    sys.exit(0)

for k, v in out.items():
    p = os.path.join(a.store, k + ".npy")
    np.save(p, v)
    print(f"wrote {p} ({os.path.getsize(p)/1e6:.1f} MB)")
