#!/usr/bin/env python3
"""backfill_v9_pass_b.py — fill the optH family on the v9 half of a merged master IN PLACE.

opth_integrate.py APPENDS its columns via pandas (and rewrites every float in the table).
This master already carries the columns for v2-v8, so here we only fill blanks, streaming
with the stdlib csv module so untouched fields keep their original text verbatim.

Per backend <be> (af3/b2/esm):
    hb_opth          -> n_h_bonds_interface_opth_<be>
    buns_ia_opth     -> delta_unsat_hbonds_opth_<be>
    buns_vsasa_raw   -> buns_vsasa_<be>
    buns_vsasa_opth  -> buns_vsasa_opth_<be>
    n_polar_h_moved  -> n_polar_h_moved_<be>

Bare consensus columns follow consensus_integrate.py: n_h_bonds_interface_opth = row-wise MIN
(worst case), every unsat flavour = row-wise MAX. Computed over the backends present, matching
the pandas min/max(skipna) that produced the v2-v8 values.

The *_raw columns are not merged; they are a verbatim recomputation of production columns, so
this asserts hb_raw == n_h_bonds_interface_<be> and buns_ia_raw == delta_unsat_hbonds_<be> on
every overlapping row and reports mismatches. Non-zero means a bad join or a stale master.
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(10 ** 7)
BE = ["af3", "b2", "esm"]
# optH source col -> master col stem
SRC = [("hb_opth", "n_h_bonds_interface_opth"), ("buns_ia_opth", "delta_unsat_hbonds_opth"),
       ("buns_vsasa_raw", "buns_vsasa"), ("buns_vsasa_opth", "buns_vsasa_opth"),
       ("n_polar_h_moved", "n_polar_h_moved")]
CHECK = [("hb_raw", "n_h_bonds_interface"), ("buns_ia_raw", "delta_unsat_hbonds")]
MINCOLS = ["n_h_bonds_interface_opth"]
MAXCOLS = ["delta_unsat_hbonds_opth", "buns_vsasa", "buns_vsasa_opth"]
NSRC, NCHK = len(SRC), len(CHECK)
WIDTH = NSRC + NCHK


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--opth", nargs="+", required=True,
                help="optH result CSVs, e.g. opth_v9_af3.csv opth_v9gap_af3.csv (backend "
                     "inferred from the filename suffix)")
a = ap.parse_args()

# ---- lookup: design_id -> [af3 vals..., b2 vals..., esm vals...] --------------------------
print("loading optH results...", flush=True)
lk = {}
for path in a.opth:
    be = next((b for b in BE if os.path.basename(path)[:-4].endswith("_" + b)), None)
    if be is None:
        sys.exit(f"ABORT: cannot infer backend from {path}")
    off = BE.index(be) * WIDTH
    n = skip = 0
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if str(r.get("opth_ok", "")).lower() not in ("true", "1"):
                skip += 1
                continue
            row = lk.get(r["design_id"])
            if row is None:
                row = lk[r["design_id"]] = [None] * (WIDTH * len(BE))
            for k, (s, _) in enumerate(SRC):
                row[off + k] = r.get(s, "")
            for k, (s, _) in enumerate(CHECK):
                row[off + NSRC + k] = r.get(s, "")
            n += 1
    print(f"  {os.path.basename(path)}: {n:,} rows ({be}), {skip:,} not opth_ok", flush=True)
print(f"  lookup: {len(lk):,} unique design_id", flush=True)

with open(a.master, newline="") as f:
    rd = csv.reader(f)
    header = [h.strip() for h in next(rd)]
    ix = {c: i for i, c in enumerate(header)}
    need = ([f"{d}_{b}" for _, d in SRC for b in BE] + MINCOLS + MAXCOLS
            + [f"{d}_{b}" for _, d in CHECK for b in BE])
    missing = [c for c in need if c not in ix]
    if missing:
        sys.exit(f"ABORT: master is missing expected columns {missing}")
    i_did = ix["design_id"]

    n = hit = filled = already = 0
    chk_n = {d: 0 for _, d in CHECK}
    chk_bad = {d: 0 for _, d in CHECK}
    tmp = a.out + ".tmp"
    with open(tmp, "w", newline="") as g:
        w = csv.writer(g, lineterminator="\n")
        w.writerow(header)
        for row in rd:
            n += 1
            src = lk.get(row[i_did])
            if src is not None:
                hit += 1
                touched = False
                for bi, b in enumerate(BE):
                    off = bi * WIDTH
                    for k, (_, dst) in enumerate(SRC):
                        v = src[off + k]
                        if v is None or v == "":
                            continue
                        j = ix[f"{dst}_{b}"]
                        if row[j] == "":
                            row[j] = v
                            touched = True
                    for k, (_, dst) in enumerate(CHECK):
                        v = fnum(src[off + NSRC + k])
                        cur = fnum(row[ix[f"{dst}_{b}"]])
                        if v is not None and cur is not None:
                            chk_n[dst] += 1
                            if abs(v - cur) > 1e-6:
                                chk_bad[dst] += 1
                # bare consensus over the backends present
                for stem, agg in [(c, min) for c in MINCOLS] + [(c, max) for c in MAXCOLS]:
                    j = ix[stem]
                    if row[j] != "":
                        continue
                    vals = [fnum(row[ix[f"{stem}_{b}"]]) for b in BE]
                    vals = [v for v in vals if v is not None]
                    if vals:
                        row[j] = f"{agg(vals):g}"
                        touched = True
                filled += touched
                already += (not touched)
            w.writerow(row)
            if n % 500000 == 0:
                print(f"  {n:,} rows | filled {filled:,}", flush=True)

os.replace(tmp, a.out)
print(f"""
wrote {a.out}
  rows                     {n:,}
  matched optH lookup      {hit:,}
  rows filled              {filled:,}
  rows already complete    {already:,}   (re-run safe)
  raw-identity checks      {dict(chk_n)}
  raw-identity MISMATCHES  {dict(chk_bad)}   (expect all 0)
""", flush=True)
