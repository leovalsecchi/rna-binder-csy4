#!/usr/bin/env python3
"""seed_lrmsd_rna.py — LRMSD in BOTH frames for per-SEED refolds vs the designed backbone.

The seed self-consistency run scored each seed against the design with stock DockQ, which
superimposes on the BINDER (receptor = larger chain) and measures the RNA. This emits the
binder-campaign direction as well, so the per-seed DockQ can be rescored target-fixed:

  lrmsd_rna    superimpose on RNA backbone     -> RMSD of protein backbone   (target fixed)
  lrmsd_prot   superimpose on protein backbone -> RMSD of RNA backbone       (= DockQ's)

The ratio of the two anchors the correction to DockQ's own LRMSD, so this implementation's
residue-matching offset cancels (see dockq_rna_recombine.py).

Flat one-row-in/one-row-out over (design_id, backend, seed) — unlike lrmsd_rna_aligned.py,
which keys by design_id and carries three backend columns per row.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import biotite.structure as struc

from rna_binder.analysis.structural import load_structure, split_chains

csv.field_size_limit(10 ** 7)

PROT_BB = ("CA", "C", "N", "O")
RNA_BB = ("P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'",
          "C1'", "C2'", "C3'", "C4'", "C5'")
KEYS = ["design_id", "backend", "seed"]
COLUMNS = ["lrmsd_rna", "lrmsd_prot", "lrmsd_error"]


def _bb(atoms, names):
    return None if atoms is None else atoms[np.isin(atoms.atom_name, list(names))]


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


def _direction(fix_ref, fix_mod, meas_ref, meas_mod):
    """Superimpose model on reference via fix_*, then RMSD of meas_* with NO second fit."""
    if any(x is None or len(x) < 3 for x in (fix_ref, fix_mod, meas_ref, meas_mod)):
        return None
    fa, fb = _match(fix_ref, fix_mod)
    ma, mb = _match(meas_ref, meas_mod)
    if len(fa) < 3 or len(ma) < 3:
        return None
    _, tr = struc.superimpose(fix_ref[fa], fix_mod[fb])
    return float(struc.rmsd(meas_ref[ma], _apply(tr, meas_mod[mb])))


def _parts(path):
    a = load_structure(path)
    a = a[a.element != "H"]
    p, r = split_chains(a)
    return _bb(p, PROT_BB), _bb(r, RNA_BB)


def _score_design(item):
    """All seeds of one design share the designed reference — load it once."""
    ref, rows = item
    out = []
    try:
        rp, rr = _parts(ref)
    except Exception as e:                                  # ref unreadable -> whole group
        return [(r, {"lrmsd_error": f"ref {type(e).__name__}: {e}"}) for r in rows]
    for r in rows:
        rec = {}
        try:
            mp, mr = _parts(r["model_path"])
            rec["lrmsd_rna"] = _direction(rr, mr, rp, mp)     # target fixed, binder measured
            rec["lrmsd_prot"] = _direction(rp, mp, rr, mr)    # DockQ's convention
        except Exception as e:
            rec["lrmsd_error"] = f"{type(e).__name__}: {e}"
        out.append((r, rec))
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()

    with open(a.csv, newline="") as f:
        rows = list(csv.DictReader(f))

    groups = {}
    for r in rows:
        groups.setdefault(r["designed_cif_path"], []).append(r)
    items = sorted(groups.items(), key=lambda kv: kv[0])
    print(f"{len(rows)} structures, {len(items)} designs, jobs={a.jobs}",
          file=sys.stderr, flush=True)

    tmp = a.output + ".tmp"
    n_err = n = 0
    with open(tmp, "w", newline="") as f, Pool(a.jobs) as pool:
        w = csv.DictWriter(f, fieldnames=KEYS + COLUMNS, extrasaction="ignore")
        w.writeheader()
        for res in pool.imap_unordered(_score_design, items):
            for r, rec in res:
                n += 1
                n_err += int("lrmsd_error" in rec)
                w.writerow({k: r[k] for k in KEYS} | rec)
    os.replace(tmp, a.output)
    print(f"wrote {a.output}: {n} rows, {n_err} errors", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
