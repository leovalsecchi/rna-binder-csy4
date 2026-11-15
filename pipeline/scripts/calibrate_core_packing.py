#!/usr/bin/env python3
"""Calibrate core_packing_score bounds: sample designs across all 20 chunks, compute the
per-core-residue Cβ coordination distribution at a few radii, print percentiles so CN_LO/CN_HI
(and UNDERPACK_CN) can be set to discriminative values. Read-only; per-task TemporaryDirectory."""
import argparse
import csv
import os
import random
import subprocess
import sys
import tempfile

import numpy as np

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, os.path.join(RB, "code/src"))
from rna_binder.analysis.structural import load_structure, split_chains
import biotite.structure as struc
from scipy.spatial import cKDTree

CORE_CUT = 30.0
RADII = [6.5, 8.0, 10.0]


def per_res_cn(protein, radius):
    protein = protein[protein.element != "H"]
    sasa = struc.sasa(protein, point_number=200)
    starts = struc.get_residue_starts(protein, add_exclusive_stop=True)
    n = len(starts) - 1
    rep = np.full((n, 3), np.nan)
    rsasa = np.zeros(n)
    names, coords = protein.atom_name, protein.coord
    for k in range(n):
        s, e = starts[k], starts[k + 1]
        rsasa[k] = np.nansum(sasa[s:e])
        seg = names[s:e]
        cb = np.where(seg == "CB")[0]
        ca = np.where(seg == "CA")[0]
        if len(cb):
            rep[k] = coords[s + cb[0]]
        elif len(ca):
            rep[k] = coords[s + ca[0]]
    valid = ~np.isnan(rep).any(axis=1)
    core = (rsasa < CORE_CUT) & valid
    if valid.sum() < 2 or core.sum() == 0:
        return []
    tree = cKDTree(rep[valid])
    return [len(tree.query_ball_point(rep[k], radius)) - 1 for k in np.where(core)[0]]


def extract(ref, td):
    if "!" not in ref:
        p = ref if os.path.isabs(ref) else os.path.join(RB, ref)
        return p if os.path.exists(p) else None
    arch, mem = ref.split("!", 1)
    arch = arch if os.path.isabs(arch) else os.path.join(RB, arch)
    out = os.path.join(td, os.path.basename(mem))
    with open(out, "wb") as fh:
        r = subprocess.run(["tar", "--zstd", "-xOf", arch, mem], stdout=fh, stderr=subprocess.DEVNULL)
    return out if (r.returncode == 0 and os.path.getsize(out) > 0) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-chunk", type=int, default=25)
    ap.add_argument("--chunkdir", default="red_chunks")
    a = ap.parse_args()
    tasks = []
    for i in range(20):
        f = os.path.join(RB, a.chunkdir, f"chunk_{i:02d}.tsv")
        if not os.path.exists(f):
            continue
        with open(f) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        random.seed(i)
        random.shuffle(rows)
        for r in rows[:a.per_chunk]:
            tasks.append((r["id"], r["ref"], r.get("stack", "")))
    sys.stderr.write(f"{len(tasks)} designs sampled across chunks\n")
    allcn = {R: [] for R in RADII}
    permean = {R: [] for R in RADII}
    stacks = {}
    ok = 0
    for rid, ref, stack in tasks:
        with tempfile.TemporaryDirectory(prefix="cal_") as td:
            cif = extract(ref, td)
            if not cif:
                continue
            try:
                prot, _ = split_chains(load_structure(cif))
            except Exception:
                continue
            if prot is None:
                continue
            ok += 1
            stacks[stack] = stacks.get(stack, 0) + 1
            for R in RADII:
                cns = per_res_cn(prot, R)
                if cns:
                    allcn[R] += cns
                    permean[R].append(float(np.mean(cns)))
    sys.stderr.write(f"scored {ok} designs; stacks={stacks}\n")
    for R in RADII:
        arr = np.array(allcn[R])
        if not len(arr):
            continue
        pcts = {p: float(np.percentile(arr, p)) for p in (5, 10, 25, 50, 75, 90, 95)}
        pm = np.array(permean[R])
        print(f"\n=== radius {R} Å : {len(arr)} core residues, {len(pm)} designs ===")
        print("  per-core-residue CN percentiles: " + "  ".join(f"p{p}={v:.1f}" for p, v in pcts.items()))
        print(f"  per-design mean CN: median {np.median(pm):.2f}  IQR [{np.percentile(pm,25):.2f}, {np.percentile(pm,75):.2f}]  range [{pm.min():.1f}, {pm.max():.1f}]")


if __name__ == "__main__":
    main()
