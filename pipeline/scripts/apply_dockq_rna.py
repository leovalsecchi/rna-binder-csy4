#!/usr/bin/env python3
"""apply_dockq_rna.py — make the RNA-aligned (target-fixed) DockQ the default in a master
metrics table, preserving the standard-convention values.

DockQ assigns receptor = larger chain, so for an ~80-residue binder vs a 16-nt target its
LRMSD superimposes on the BINDER and measures the RNA. For a binder campaign the target
should be held fixed instead. Only the LRMSD term and the composite score change — fnat and
iRMSD are frame-independent.

Per backend <be> and for the consensus:
    dockq_std_LRMSD_<be>  <- old dockq_LRMSD_<be>      (preserved)
    dockq_std_self_<be>   <- old dockq_self_<be>       (preserved)
    dockq_LRMSD_<be>      <- lrmsd_rna_corr_<be>       (new default)
    dockq_self_<be>       <- dockq_rna_<be>            (new default)
Consensus is rebuilt with the project convention: min for higher-is-better (dockq_self),
max for lower-is-better (dockq_LRMSD), NaN unless all three backends scored.

Streams row by row (stdlib csv), so a 6 GB master costs little memory.
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(10 ** 7)
BE = ["af3", "b2", "esm"]


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--sidecar", required=True, help="dockq_rna_*.csv from dockq_rna_recombine")
ap.add_argument("--out", required=True)
a = ap.parse_args()

new = {}
with open(a.sidecar, newline="") as f:
    for r in csv.DictReader(f):
        new[r["design_id"]] = r
print(f"sidecar: {len(new):,} designs", flush=True)

with open(a.master, newline="") as f:
    rd = csv.reader(f)
    header = next(rd)

    add = []
    for b in BE:
        for k in ("LRMSD", "self"):
            if f"dockq_{k}_{b}" in header:
                add.append(f"dockq_std_{k}_{b}")
    for k in ("LRMSD", "self"):
        if f"dockq_{k}" in header:
            add.append(f"dockq_std_{k}")
    out_header = header + add
    idx = {c: i for i, c in enumerate(header)}
    print(f"master: {len(header)} cols -> {len(out_header)} (+{len(add)} preserved)", flush=True)

    tmp = a.out + ".tmp"
    n = n_hit = 0
    with open(tmp, "w", newline="") as g:
        w = csv.writer(g)
        w.writerow(out_header)
        for row in rd:
            n += 1
            did = row[idx["design_id"]]
            src = new.get(did)
            extra = {c: "" for c in add}
            if src:
                n_hit += 1
                per_self, per_lr = [], []
                for b in BE:
                    cs, cl = f"dockq_self_{b}", f"dockq_LRMSD_{b}"
                    if cs in idx:
                        extra[f"dockq_std_self_{b}"] = row[idx[cs]]
                        row[idx[cs]] = src.get(f"dockq_rna_{b}", "")
                        per_self.append(fnum(row[idx[cs]]))
                    if cl in idx:
                        extra[f"dockq_std_LRMSD_{b}"] = row[idx[cl]]
                        row[idx[cl]] = src.get(f"lrmsd_rna_corr_{b}", "")
                        per_lr.append(fnum(row[idx[cl]]))
                # consensus: gated on all three, min for DockQ, max for LRMSD
                if "dockq_self" in idx:
                    extra["dockq_std_self"] = row[idx["dockq_self"]]
                    ok = [v for v in per_self if v is not None]
                    row[idx["dockq_self"]] = f"{min(ok):.6f}" if len(ok) == 3 else ""
                if "dockq_LRMSD" in idx:
                    extra["dockq_std_LRMSD"] = row[idx["dockq_LRMSD"]]
                    ok = [v for v in per_lr if v is not None]
                    row[idx["dockq_LRMSD"]] = f"{max(ok):.4f}" if len(ok) == 3 else ""
            w.writerow(row + [extra[c] for c in add])
            if n % 500000 == 0:
                print(f"  {n:,} rows", flush=True)
os.replace(tmp, a.out)
print(f"wrote {a.out}: {n:,} rows, {n_hit:,} updated ({100*n_hit/max(n,1):.2f}%)", flush=True)
