#!/usr/bin/env python3
"""lrmsd_convention.py — how much does DockQ's receptor/ligand choice change LRMSD?

DockQ assigns receptor = LARGER chain. For an ~80-residue binder vs a 16-nt RNA that makes
the PROTEIN the receptor, so DockQ's LRMSD superimposes on the binder and measures the RNA.
The natural convention for a binder campaign is the reverse: hold the TARGET fixed and
measure where the binder went.

For every structure, versus the designed backbone, this computes both:

  lrmsd_protein_aligned   superimpose on protein backbone -> RMSD of RNA backbone   (DockQ)
  lrmsd_rna_aligned       superimpose on RNA backbone     -> RMSD of protein backbone

No second fit in either case. Atom sets are DockQ's BACKBONE_ATOMS, matched by
(chain_kind, res_id, atom_name), so the only difference between the two numbers is which
partner defines the frame.

Emitting lrmsd_protein_aligned lets us check against the stored dockq_LRMSD: if they agree,
the counterfactual column is trustworthy.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import biotite.structure as struc

from rna_binder.analysis.structural import load_structure, split_chains

# DockQ constants.py BACKBONE_ATOMS
PROT_BB = {"CA", "C", "N", "O"}
RNA_BB = {"P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'",
          "C1'", "C2'", "C3'", "C4'", "C5'"}
PAT = re.compile(r"^(?P<did>.+?)_(?:AF3|B2|ESM)(?:_seed(?P<seed>\d+))?_sample(?P<samp>\d+)\.cif$")


def _bb(atoms, names):
    if atoms is None:
        return None
    return atoms[np.isin(atoms.atom_name, list(names))]


def _match(a, b):
    pos = {}
    for i, (r, n) in enumerate(zip(a.res_id, a.atom_name)):
        pos.setdefault((int(r), str(n)), i)
    ia, ib = [], []
    for j, (r, n) in enumerate(zip(b.res_id, b.atom_name)):
        i = pos.get((int(r), str(n)))
        if i is not None:
            ia.append(i)
            ib.append(j)
    return np.asarray(ia, int), np.asarray(ib, int)


def _apply(tr, atoms):
    if hasattr(tr, "apply"):
        return tr.apply(atoms)
    return struc.superimpose_apply(atoms, tr)


def _one_direction(fix_ref, fix_mod, meas_ref, meas_mod):
    """Superimpose model onto reference using fix_*, then RMSD of meas_* with no re-fit."""
    if any(x is None or len(x) < 3 for x in (fix_ref, fix_mod, meas_ref, meas_mod)):
        return None
    fa, fb = _match(fix_ref, fix_mod)
    ma, mb = _match(meas_ref, meas_mod)
    if len(fa) < 3 or len(ma) < 3:
        return None
    _, tr = struc.superimpose(fix_ref[fa], fix_mod[fb])
    moved = _apply(tr, meas_mod[mb])
    return float(struc.rmsd(meas_ref[ma], moved))


def score(model_path, ref_path):
    m = load_structure(model_path)
    r = load_structure(ref_path)
    m = m[m.element != "H"]
    r = r[r.element != "H"]
    mp, mr = split_chains(m)
    rp, rr = split_chains(r)
    mp_bb, rp_bb = _bb(mp, PROT_BB), _bb(rp, PROT_BB)
    mr_bb, rr_bb = _bb(mr, RNA_BB), _bb(rr, RNA_BB)
    return {
        # DockQ's convention: receptor = larger chain = the protein here
        "lrmsd_protein_aligned": _one_direction(rp_bb, mp_bb, rr_bb, mr_bb),
        # binder-campaign convention: hold the RNA target fixed
        "lrmsd_rna_aligned": _one_direction(rr_bb, mr_bb, rp_bb, mp_bb),
        "n_prot_bb": 0 if rp_bb is None else len(rp_bb),
        "n_rna_bb": 0 if rr_bb is None else len(rr_bb),
    }


_BACKEND = None


def _score_item(it):
    did, seed, model_path, ref_path = it
    try:
        return {"design_id": did, "backend": _BACKEND, "seed": seed,
                **score(model_path, ref_path)}
    except Exception as e:
        return {"design_id": did, "backend": _BACKEND, "seed": seed,
                "error": f"{type(e).__name__}: {e}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="dir of extracted <...>.cif structures")
    ap.add_argument("--refs", required=True, help="design_id,designed_cif_path (loose)")
    ap.add_argument("--backend", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    refs = {r["design_id"]: r["designed_cif_path"] for r in csv.DictReader(open(a.refs))}
    items = []
    for fn in sorted(os.listdir(a.dir)):
        m = PAT.match(fn)
        if not m:
            continue
        did = m.group("did")
        if did in refs:
            items.append((did, int(m.group("seed") or m.group("samp")),
                          os.path.join(a.dir, fn), refs[did]))
    if a.limit:
        items = items[:a.limit]
    print(f"{len(items)} structures", flush=True)

    global _BACKEND
    _BACKEND = a.backend
    with Pool(a.jobs) as pool:
        rows = pool.map(_score_item, items, chunksize=8)

    cols = ["design_id", "backend", "seed", "lrmsd_protein_aligned",
            "lrmsd_rna_aligned", "n_prot_bb", "n_rna_bb", "error"]
    with open(a.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if r.get("lrmsd_rna_aligned") is not None)
    print(f"wrote {a.output}: {len(rows)} rows, {ok} scored", flush=True)


if __name__ == "__main__":
    main()
