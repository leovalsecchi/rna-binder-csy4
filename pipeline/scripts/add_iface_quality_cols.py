#!/usr/bin/env python3
"""add_iface_quality_cols.py — add the interface-quality family to the v2-v9 filter store.

Five Rosetta interface descriptors, per backend, for the filtered-vs-baseline comparison:

    sc_value            shape complementarity
    packstat            packing quality
    delta_unsat_hbonds  buried unsatisfied H-bonds
    dsasa_polar         polar buried area
    dsasa_hphobic       hydrophobic buried area

They come from two places because v9 was never merged into the master:
  v2-v8  combined_clean_v2v8_fullmetrics.csv, already suffixed `_af3` / `_esm`
  v9     iface_v9/.../iface_be_<be>.csv, unsuffixed (backend is the file)

⚠️ Only the `_af3`/`_esm` columns are read. The master also carries BARE `sc_value`,
`packstat`, `dsasa_polar` etc. which fold Boltz-2 in — reading those would put b2 into a tool
that is supposed to have none. The b2-suffixed columns are equally off-limits. Asserted below.

Rows are matched by design_id, never by position: the store is v2-v8 then v9, but the master
is not guaranteed to enumerate in store order.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

csv.field_size_limit(10 ** 7)

METRICS = ["sc_value", "packstat", "delta_unsat_hbonds", "dsasa_polar", "dsasa_hphobic"]
BE = ("af3", "esm")

ap = argparse.ArgumentParser()
ap.add_argument("--store", required=True)
ap.add_argument("--master", required=True, help="combined_clean_v2v8_fullmetrics.csv")
ap.add_argument("--iface-af3", required=True)
ap.add_argument("--iface-esm", required=True)
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

for m in METRICS:
    for b in BE:
        if "b2" in f"{m}_{b}".lower().split("_"):
            sys.exit(f"REFUSING: {m}_{b} names b2")

did = np.load(os.path.join(a.store, "design_id.npy"))
n = did.shape[0]
pos = {d: i for i, d in enumerate(did.tolist())}
print(f"store: {n:,} rows", flush=True)
del did

out = {f"{m}_{b}": np.full(n, np.nan, dtype=np.float32) for m in METRICS for b in BE}

# ---- v2-v8: one streaming pass over the 9.7 GB master (never pd.read_csv this file)
want = [f"{m}_{b}" for m in METRICS for b in BE]
hit = 0
with open(a.master, newline="") as f:
    r = csv.DictReader(f)
    missing = [c for c in want if c not in (r.fieldnames or [])]
    if missing:
        sys.exit(f"FATAL: master lacks {missing}")
    for row in r:
        i = pos.get(row["design_id"].encode())
        if i is None:
            continue
        hit += 1
        for c in want:
            s = row[c]
            if s not in ("", "NA", "nan", "NaN"):
                try:
                    out[c][i] = float(s)
                except ValueError:
                    pass
        if hit % 1_000_000 == 0:
            print(f"  master: {hit:,} matched", flush=True)
print(f"master: {hit:,} rows matched into the store", flush=True)

# ---- v9: per-backend interface sidecars, columns unsuffixed
for be, path in (("af3", a.iface_af3), ("esm", a.iface_esm)):
    got = 0
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        missing = [m for m in METRICS if m not in (r.fieldnames or [])]
        if missing:
            sys.exit(f"FATAL: {path} lacks {missing}")
        for row in r:
            i = pos.get(row["design_id"].encode())
            if i is None:
                continue
            got += 1
            for m in METRICS:
                s = row[m]
                if s not in ("", "NA", "nan", "NaN"):
                    try:
                        out[f"{m}_{be}"][i] = float(s)
                    except ValueError:
                        pass
    print(f"v9 {be}: {got:,} rows matched", flush=True)

print()
for c in want:
    v = out[c]
    ok = np.isfinite(v)
    med = float(np.median(v[ok])) if ok.any() else float("nan")
    print(f"  {c:26s} {ok.sum():>9,}/{n:,} ({100*ok.sum()/n:5.2f}%)  median {med:.3f}")

if a.dry_run:
    print("\ndry run — nothing written")
    sys.exit(0)
for c in want:
    np.save(os.path.join(a.store, c + ".npy"), out[c])
print(f"\nwrote {len(want)} columns into {a.store}")
