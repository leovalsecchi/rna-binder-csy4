#!/usr/bin/env python3
"""dockq_rna_recombine.py — rescore DockQ in the binder-campaign frame.

DockQ assigns receptor = larger chain, so for an ~80-residue binder vs a 16-nt target its
LRMSD term superimposes on the BINDER and measures the RNA. Holding the target fixed
instead measures where the binder went — the quantity a binder campaign cares about.

Only LRMSD changes: fnat is contact counting and iRMSD superimposes on BOTH chains'
interface residues, so both are frame-independent and are reused verbatim from the stored
columns. The rescored DockQ is therefore exact, not an approximation:

    DockQ_rna = [ fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD_rna/8.5)^2) ] / 3

LRMSD_rna is anchored to DockQ's own output via the RATIO measured by lrmsd_rna_aligned.py
(both directions from one implementation), so that implementation's ~4% offset vs DockQ's
sequence-alignment residue matching cancels:

    LRMSD_rna = dockq_LRMSD_stored * (lrmsd_rna_mine / lrmsd_prot_mine)

Adds columns; never overwrites dockq_*, so the standard score stays the default.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

BE = ["af3", "b2", "esm"]
GATE = 0.4


def dockq(fnat, irms, lrms):
    return (fnat + 1 / (1 + (irms / 1.5) ** 2) + 1 / (1 + (lrms / 8.5) ** 2)) / 3


ap = argparse.ArgumentParser()
ap.add_argument("--sidecar-dir", required=True)
ap.add_argument("--master", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--label", default="")
a = ap.parse_args()

files = sorted(glob.glob(os.path.join(a.sidecar_dir, "lrmsd_*.csv")))
cols = [f"lrmsd_rna_{b}" for b in BE] + [f"lrmsd_prot_{b}" for b in BE]
sc = pd.concat([pd.read_csv(f, usecols=["design_id"] + cols, low_memory=False)
                for f in files], ignore_index=True).drop_duplicates("design_id")
print(f"{len(files)} sidecars, {len(sc):,} designs", flush=True)

need = ["design_id", "stack"]
for b in BE:
    need += [f"dockq_fnat_{b}", f"dockq_iRMSD_{b}", f"dockq_LRMSD_{b}", f"dockq_self_{b}"]
head = pd.read_csv(a.master, nrows=0).columns.tolist()
use = [c for c in need if c in head]
m = pd.read_csv(a.master, usecols=use, low_memory=False)
print(f"master: {len(m):,} rows", flush=True)

d = m.merge(sc, on="design_id", how="left")
for b in BE:
    ratio = d[f"lrmsd_rna_{b}"] / d[f"lrmsd_prot_{b}"].where(d[f"lrmsd_prot_{b}"] > 0.01)
    d[f"lrmsd_rna_corr_{b}"] = d[f"dockq_LRMSD_{b}"] * ratio
    d[f"dockq_rna_{b}"] = dockq(d[f"dockq_fnat_{b}"], d[f"dockq_iRMSD_{b}"],
                                d[f"lrmsd_rna_corr_{b}"])

trio = [f"dockq_rna_{b}" for b in BE]
n = d[trio].notna().sum(axis=1)
d["dockq_rna_consensus_n"] = n
d["dockq_rna"] = d[trio].min(axis=1).where(n.eq(3))          # worst-case, gated on all 3
old = [f"dockq_self_{b}" for b in BE if f"dockq_self_{b}" in d.columns]
d["dockq_self_consensus"] = d[old].min(axis=1).where(d[old].notna().sum(axis=1).eq(3))

out_cols = (["design_id", "stack"] + [f"lrmsd_rna_corr_{b}" for b in BE] + trio +
            ["dockq_rna_consensus_n", "dockq_rna"])
d[out_cols].to_csv(a.out, index=False)
print(f"wrote {a.out}: {len(d):,} rows", flush=True)

lab = a.label or os.path.basename(a.sidecar_dir)
print(f"\n===== {lab} =====", flush=True)
for b in BE:
    s, r = d[f"dockq_self_{b}"], d[f"dockq_rna_{b}"]
    k = s.notna() & r.notna()
    print(f"  {b:4s} n={int(k.sum()):>9,}  DockQ {s[k].median():.4f} -> {r[k].median():.4f}"
          f"  (median shift {(s[k] - r[k]).median():+.4f})"
          f"  pass≥{GATE}: {100 * (s[k] >= GATE).mean():5.2f}% -> {100 * (r[k] >= GATE).mean():5.2f}%")

k = d.dockq_self_consensus.notna() & d.dockq_rna.notna()
p_old, p_new = d.dockq_self_consensus[k] >= GATE, d.dockq_rna[k] >= GATE
print(f"\n  CONSENSUS n={int(k.sum()):,}")
print(f"    median {d.dockq_self_consensus[k].median():.4f} -> {d.dockq_rna[k].median():.4f}")
print(f"    pass: {int(p_old.sum()):,} ({100 * p_old.mean():.2f}%) -> "
      f"{int(p_new.sum()):,} ({100 * p_new.mean():.2f}%)")
print(f"    lose pass: {int((p_old & ~p_new).sum()):,} "
      f"({100 * (p_old & ~p_new).sum() / max(p_old.sum(), 1):.1f}% of current)   "
      f"gain: {int((~p_old & p_new).sum()):,}")

print("\n  by stack (consensus):")
for st, g in d[k].groupby("stack"):
    po, pn = g.dockq_self_consensus >= GATE, g.dockq_rna >= GATE
    print(f"    {st:14s} n={len(g):>8,}  pass {100 * po.mean():5.2f}% -> {100 * pn.mean():5.2f}%"
          f"   lose {100 * (po & ~pn).sum() / max(po.sum(), 1):5.1f}% of its passes")
