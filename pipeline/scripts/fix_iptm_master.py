#!/usr/bin/env python3
"""fix_iptm_master.py — put honest ipTM consensus columns on a master table.

The existing aggregates are not trustworthy and differ per table:
  app_metrics_100k   consensus_iptm = min(AF3,B2), mean_iptm = mean(AF3,B2)   [ESM excluded]
  v2-v7 master       consensus_iptm is MIXED — min-of-3 on v6/v7, best-2-of-3 on v2-v5
  v8 master          no consensus_iptm / mean_iptm at all

so every table gets the same three, recomputed from the per-backend columns and gated on
all three refolds being present:

  iptm_min3    min(af3, b2, esm)      strict consensus — all three must agree
  iptm_best2   median of the three    = min of the BEST two (drops the worst)
  mean_iptm3   mean of the three

`consensus_iptm` is then set to iptm_min3 everywhere. Whatever it held before is preserved
as `consensus_iptm_legacy`, and `mean_iptm` as `mean_iptm_legacy`, so the old numbers stay
recoverable and a merged table can tell "was never computed" (blank) from "was computed
under the old rule".

Streams row by row with the stdlib csv module — a 6.8 GB master would not fit in memory as
a DataFrame. Writes LF line endings (csv.writer defaults to CRLF, which is how the v8 master
ended up mixed).
"""
import argparse
import csv
import os

csv.field_size_limit(10 ** 7)
BE = ["af3_msaoff_iptm", "b2_iptm", "esm_iptm"]


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


ap = argparse.ArgumentParser()
ap.add_argument("--table", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

with open(a.table, newline="") as f:
    rd = csv.reader(f)
    header = [h.strip() for h in next(rd)]
    idx = {c: i for i, c in enumerate(header)}
    missing = [c for c in BE if c not in idx]
    if missing:
        raise SystemExit(f"ABORT: missing per-backend ipTM columns {missing}")

    out_header = list(header)
    for old, new in (("consensus_iptm", "consensus_iptm_legacy"),
                     ("mean_iptm", "mean_iptm_legacy")):
        if old in idx:
            out_header[idx[old]] = new
    added = ["iptm_min3", "iptm_best2", "mean_iptm3", "consensus_iptm"]
    out_header += added
    print(f"{os.path.basename(a.table)}: {len(header)} -> {len(out_header)} cols "
          f"(renamed {[c for c in ('consensus_iptm','mean_iptm') if c in idx]})", flush=True)

    n = n_all3 = 0
    tmp = a.out + ".tmp"
    with open(tmp, "w", newline="") as g:
        w = csv.writer(g, lineterminator="\n")
        w.writerow(out_header)
        for row in rd:
            n += 1
            v = [fnum(row[idx[c]]) if idx[c] < len(row) else None for c in BE]
            ok = [x for x in v if x is not None]
            if len(ok) == 3:
                n_all3 += 1
                ok.sort()
                mn3, best2, mean3 = ok[0], ok[1], sum(ok) / 3.0
                extra = [f"{mn3:.6g}", f"{best2:.6g}", f"{mean3:.6g}", f"{mn3:.6g}"]
            else:
                extra = ["", "", "", ""]
            w.writerow(row + extra)
            if n % 500000 == 0:
                print(f"  {n:,} rows", flush=True)
os.replace(tmp, a.out)
print(f"wrote {a.out}: {n:,} rows, {n_all3:,} with all three ipTMs "
      f"({100 * n_all3 / max(n, 1):.2f}%)", flush=True)
