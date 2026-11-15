#!/usr/bin/env python3
"""concat_stores.py — append one columnar filter store onto another.

Used to fold v9 into the v2-v8 store without rebuilding 6.8M rows from a combined CSV that
does not exist. Row order is A-then-B throughout, so every column stays aligned.

Three things here are not memcpy and are where this would silently corrupt data:

1. **Categorical codes are per-store.** `stack`/`campaign`/`gen`/`invfold` are int8 indices
   into that store's own vocab list. Store B's code 3 is not store A's code 3. Every B code is
   remapped through a merged vocab; a value missing from A's vocab is appended, not dropped.

2. **Path templates are per-store.** `tpl_idx` indexes that store's template list. B's
   templates are appended to A's (deduplicated by exact tuple) and B's indices remapped.
   `int16` is asserted to still fit.

3. **Column sets can differ.** A column present in only one store would misalign the rows of
   the other. Missing float columns are filled with NaN (which the filter fails closed on, the
   correct behaviour for "not computed for that campaign"); a missing non-float column is a
   hard error, because there is no safe fill for an id, a code, or a template index.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

import numpy as np

VOCAB_COLS = ("stack", "campaign", "gen", "invfold")

ap = argparse.ArgumentParser()
ap.add_argument("--a", required=True, help="base store (rows come first)")
ap.add_argument("--b", required=True, help="store to append")
ap.add_argument("--out", required=True)
ap.add_argument("--force", action="store_true")
args = ap.parse_args()

if os.path.exists(args.out):
    if not args.force:
        sys.exit(f"FATAL: {args.out} exists (use --force)")
    shutil.rmtree(args.out)
os.makedirs(args.out)

ma = json.load(open(os.path.join(args.a, "meta.json")))
mb = json.load(open(os.path.join(args.b, "meta.json")))
na, nb = int(ma["n"]), int(mb["n"])
print(f"A: {na:,} rows   B: {nb:,} rows   -> {na+nb:,}")

if list(ma["path_cols"]) != list(mb["path_cols"]):
    sys.exit(f"FATAL: path_cols differ:\n  A={ma['path_cols']}\n  B={mb['path_cols']}")

# ---- templates: append B's, dedup by exact tuple, remap B's indices
ta = [tuple(t) for t in ma["templates"]]
tb = [tuple(t) for t in mb["templates"]]
merged = list(ta)
index = {t: i for i, t in enumerate(merged)}
tmap = np.empty(len(tb), dtype=np.int64)
for i, t in enumerate(tb):
    if t not in index:
        index[t] = len(merged)
        merged.append(t)
    tmap[i] = index[t]
print(f"templates: A={len(ta)} B={len(tb)} -> merged={len(merged)} "
      f"({len(merged)-len(ta)} new from B)")
if len(merged) > np.iinfo(np.int16).max:
    sys.exit("FATAL: merged template count exceeds int16")

# ---- vocabs: union, remap B's codes
vocabs, vmap = {}, {}
for c in VOCAB_COLS:
    va = list(ma["vocabs"].get(c, []))
    vb = list(mb["vocabs"].get(c, []))
    out = list(va)
    pos = {v: i for i, v in enumerate(out)}
    m = np.empty(len(vb), dtype=np.int64)
    for i, v in enumerate(vb):
        if v not in pos:
            pos[v] = len(out)
            out.append(v)
        m[i] = pos[v]
    vocabs[c], vmap[c] = out, m
    added = [v for v in vb if v not in va]
    print(f"vocab {c:9s} A={len(va)} B={len(vb)} -> {len(out)}" +
          (f"  new: {added}" if added else ""))
    if len(out) > 127:
        sys.exit(f"FATAL: vocab {c} exceeds int8 range")

# ---- columns
cols_a = {f[:-4] for f in os.listdir(args.a) if f.endswith(".npy")}
cols_b = {f[:-4] for f in os.listdir(args.b) if f.endswith(".npy")}
only_a, only_b = sorted(cols_a - cols_b), sorted(cols_b - cols_a)
if only_a:
    print(f"columns only in A (B filled NaN): {only_a}")
if only_b:
    print(f"columns only in B (A filled NaN): {only_b}")

for name in sorted(cols_a | cols_b):
    pa = os.path.join(args.a, name + ".npy")
    pb = os.path.join(args.b, name + ".npy")
    if os.path.exists(pa) and os.path.exists(pb):
        va, vb_ = np.load(pa), np.load(pb)
        if name == "tpl_idx":
            vb_ = tmap[vb_.astype(np.int64)].astype(va.dtype)
        elif name in VOCAB_COLS:
            vb_ = vmap[name][vb_.astype(np.int64)].astype(va.dtype)
        elif va.dtype.kind == "S" and va.dtype != vb_.dtype:
            w = max(va.dtype.itemsize, vb_.dtype.itemsize)
            va, vb_ = va.astype(f"S{w}"), vb_.astype(f"S{w}")
        out = np.concatenate([va, vb_])
    else:
        have, n_have, n_missing, first = ((np.load(pa), na, nb, True) if os.path.exists(pa)
                                          else (np.load(pb), nb, na, False))
        if have.dtype.kind != "f":
            sys.exit(f"FATAL: {name} missing from one store and is not float "
                     f"({have.dtype}) — no safe fill for an id/code/template index")
        pad = np.full(n_missing, np.nan, dtype=have.dtype)
        out = np.concatenate([have, pad] if first else [pad, have])
    if len(out) != na + nb:
        sys.exit(f"FATAL: {name} length {len(out)} != {na+nb}")
    np.save(os.path.join(args.out, name + ".npy"), out)

meta = dict(ma)
meta["n"] = na + nb
meta["templates"] = [list(t) for t in merged]
meta["vocabs"] = vocabs
meta["source"] = f"concat({ma.get('source','?')}, {mb.get('source','?')})"
meta["concat_parts"] = [{"store": os.path.abspath(args.a), "n": na},
                        {"store": os.path.abspath(args.b), "n": nb}]
json.dump(meta, open(os.path.join(args.out, "meta.json"), "w"), indent=1)
print(f"\nwrote {args.out}: {na+nb:,} rows, {len(sorted(cols_a | cols_b))} columns")
