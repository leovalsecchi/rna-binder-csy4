#!/usr/bin/env python3
"""add_rog_fbase_cols.py — add radius-of-gyration and base-contact-fraction to the filter store.

    rog_af3 / rog_esm      <- radius_of_gyration_ca_{af3,esm}   (overwrites the all-NaN placeholder)
    fbase_af3 / fbase_esm  <- frac_base_contacts_{af3,esm}      (new column)

Single source: combined_clean_v2v9_fullmetrics.csv now carries both at ~100% for v2-v9, so
unlike add_iface_quality_cols.py this needs no separate v9 sidecar.

Only _af3/_esm are read. Bare `frac_base_contacts` folds Boltz-2 in and is refused, as is any
b2-suffixed column. Rows are matched by design_id, never by position.
Additive: writes only these four .npy and rewrites meta.json. No other column is touched.
"""
from __future__ import annotations
import argparse, csv, json, os, sys
import numpy as np
csv.field_size_limit(10 ** 7)

SRC = {"rog": "radius_of_gyration_ca", "fbase": "frac_base_contacts"}
BE = ("af3", "esm")

ap = argparse.ArgumentParser()
ap.add_argument("--store", required=True)
ap.add_argument("--master", required=True)
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

for k, src in SRC.items():
    for b in BE:
        if "b2" in f"{src}_{b}".lower().split("_"):
            sys.exit(f"REFUSING: {src}_{b} names b2")

did = np.load(os.path.join(a.store, "design_id.npy"))
n = did.shape[0]
pos = {d: i for i, d in enumerate(did.tolist())}
print(f"store: {n:,} rows", flush=True)
del did

out = {f"{k}_{b}": np.full(n, np.nan, dtype=np.float32) for k in SRC for b in BE}
want = {f"{k}_{b}": f"{src}_{b}" for k, src in SRC.items() for b in BE}

hit = 0
with open(a.master, newline="") as f:
    r = csv.DictReader(f)
    for col in want.values():
        if col not in r.fieldnames: sys.exit(f"MASTER LACKS {col}")
    for row in r:
        i = pos.get(row["design_id"].encode() if isinstance(next(iter(pos)), bytes) else row["design_id"])
        if i is None: continue
        hit += 1
        for dst, srccol in want.items():
            v = row[srccol]
            if v not in ("", "nan", "NaN"):
                try: out[dst][i] = float(v)
                except ValueError: pass
print(f"matched {hit:,} / {n:,} store rows", flush=True)
for k, arr in out.items():
    fin = int(np.isfinite(arr).sum())
    print(f"  {k:12s} {fin:>9,} / {n:,}  ({100*fin/n:.2f}%)")

if a.dry_run:
    print("DRY RUN - nothing written"); raise SystemExit

for k, arr in out.items():
    np.save(os.path.join(a.store, k + ".npy"), arr)
mp = os.path.join(a.store, "meta.json")
m = json.load(open(mp))
m["placeholder_cols"] = [c for c in m.get("placeholder_cols", []) if c not in out]
m.setdefault("float_cols", {}).update(want)
json.dump(m, open(mp, "w"), indent=1)
print("WROTE 4 columns + meta.json updated")
