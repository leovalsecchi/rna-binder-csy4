#!/usr/bin/env python3
"""merge_masters.py — concatenate the v2-v7 and v8 masters into one table.

v8's columns are a strict SUBSET of v2-v7's (verified: 144 of 167, zero v8-only), so the
union schema is just v2-v7's. v8 rows get blanks for the 23 it lacks — the PyRosetta
optimised-hydrogen family (`*_opth`, `buns_vsasa*`, `n_polar_h_moved_*`), which was never
computed for v8, plus the legacy bookkeeping columns.

`campaign` is filled in for v8 rows (v2-v7 already carries it). A `source_table` column
records which master each row came from, so the NaN block is attributable rather than
mysterious.

Streams both files — the result is ~5.3M rows / ~9 GB and cannot be held in memory. Output
is LF-terminated regardless of what the inputs used.
"""
import argparse
import csv
import os

csv.field_size_limit(10 ** 7)

ap = argparse.ArgumentParser()
ap.add_argument("--a", required=True, help="v2-v7 master (defines the schema)")
ap.add_argument("--b", required=True, help="v8 master")
ap.add_argument("--b-campaign", default="v8")
ap.add_argument("--a-label", default="v2v7")
ap.add_argument("--b-label", default="v8")
ap.add_argument("--out", required=True)
a = ap.parse_args()


def header_of(p):
    with open(p, newline="") as f:
        return [h.strip() for h in next(csv.reader(f))]

ha, hb = header_of(a.a), header_of(a.b)
extra_b = [c for c in hb if c not in set(ha)]
if extra_b:
    raise SystemExit(f"ABORT: {a.b} has {len(extra_b)} columns not in the schema: {extra_b[:10]}")
out_header = ha + ["source_table"]
missing_in_b = [c for c in ha if c not in set(hb)]
print(f"schema: {len(ha)} cols (+source_table); v8 lacks {len(missing_in_b)}: "
      f"{', '.join(missing_in_b)}", flush=True)

tmp = a.out + ".tmp"
counts = {}
with open(tmp, "w", newline="") as g:
    w = csv.writer(g, lineterminator="\n")
    w.writerow(out_header)

    for path, label, camp in ((a.a, a.a_label, None), (a.b, a.b_label, a.b_campaign)):
        with open(path, newline="") as f:
            rd = csv.reader(f)
            hdr = [h.strip() for h in next(rd)]
            pos = {c: i for i, c in enumerate(hdr)}
            # build the projection once: for each output column, where to read it from
            take = [pos.get(c, -1) for c in ha]
            ci = ha.index("campaign") if "campaign" in ha else -1
            n = 0
            for row in rd:
                n += 1
                out = [(row[j] if 0 <= j < len(row) else "") for j in take]
                if camp and ci >= 0 and not out[ci]:
                    out[ci] = camp
                w.writerow(out + [label])
                if n % 500000 == 0:
                    print(f"  {label}: {n:,} rows", flush=True)
            counts[label] = n
            print(f"  {label}: {n:,} rows done", flush=True)

os.replace(tmp, a.out)
total = sum(counts.values())
print(f"\nwrote {a.out}: {total:,} rows ({' + '.join(f'{k} {v:,}' for k, v in counts.items())})",
      flush=True)
