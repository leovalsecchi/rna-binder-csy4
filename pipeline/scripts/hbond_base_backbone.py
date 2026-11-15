#!/usr/bin/env python3
"""hbond_base_backbone.py — interface H-bonds split by which part of the RNA they reach.

The stored H-bond columns (n_h_bonds_interface*, and the optH variants) are per-design
TOTALS, and the aa x nt map resolves the nucleotide but not the ATOM. Neither can answer
"is the binder hydrogen-bonding to the BASE or to the BACKBONE" — i.e. is it reading
sequence, or just gripping the phosphate cage.

Same engine as hbond_map_pyrosetta.py / hbonds_pyrosetta.py (Rosetta builds hydrogens and
does orientation-dependent detection, so the counts are comparable to the production
columns). For every protein<->RNA H-bond it classifies the RNA-side heavy atom:

    phosphate   P, OP1, OP2, O5', O3'
    sugar       C1'-C5', O4', O2'          (O2' is the 2'-OH — sugar, not base)
    base        everything else (ring N/C, exocyclic O/N)

Both directions are counted: protein donor -> RNA acceptor, and RNA donor -> protein
acceptor (the RNA heavy atom is then the donor base atom, not the hydrogen).

Emits per design: hb_base, hb_sugar, hb_phosphate, hb_backbone (= sugar+phosphate),
hb_total, and hb_base_frac. Run inside pyrosetta.sif.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive_read import resolve

csv.field_size_limit(10 ** 7)

PHOSPHATE = {"P", "OP1", "OP2", "OP3", "O5'", "O3'", "O1P", "O2P", "O3P"}
SUGAR = {"C1'", "C2'", "C3'", "C4'", "C5'", "O4'", "O2'"}
COLUMNS = ["hb_total", "hb_base", "hb_sugar", "hb_phosphate", "hb_backbone",
           "hb_base_frac", "hb_error"]
_SF = None


def _worker_init():
    global _SF
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false")
    _SF = pyrosetta.get_fa_scorefxn()


def _cls(name):
    n = name.strip()
    if n in PHOSPHATE:
        return "phosphate"
    if n in SUGAR:
        return "sugar"
    return "base"


def _one(path):
    import pyrosetta
    from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
    try:
        pose = pyrosetta.pose_from_file(resolve(path))
        _SF(pose)
        pose.update_residue_neighbors()
        hb = HBondSet()
        fill_hbond_set(pose, False, hb)
        c = defaultdict(int)
        for i in range(1, hb.nhbonds() + 1):
            b = hb.hbond(i)
            dr, ar = b.don_res(), b.acc_res()
            rd, ra = pose.residue(dr), pose.residue(ar)
            if rd.is_protein() and (ra.is_RNA() or ra.is_DNA()):
                c[_cls(ra.atom_name(b.acc_atm()))] += 1
            elif ra.is_protein() and (rd.is_RNA() or rd.is_DNA()):
                # RNA is the donor: the H's parent heavy atom is the RNA partner
                c[_cls(rd.atom_name(rd.atom_base(b.don_hatm())))] += 1
        base, sug, pho = c["base"], c["sugar"], c["phosphate"]
        tot = base + sug + pho
        return path, {"hb_total": tot, "hb_base": base, "hb_sugar": sug,
                      "hb_phosphate": pho, "hb_backbone": sug + pho,
                      "hb_base_frac": round(base / tot, 4) if tot else ""}
    except Exception as e:
        return path, {"hb_error": f"{type(e).__name__}: {e}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--path-col", default="af3_msaoff_cif_path")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()

    with open(a.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    refs = sorted({(r.get(a.path_col) or "").strip() for r in rows} - {""})
    print(f"{len(rows)} rows, {len(refs)} structures, jobs={a.jobs}", file=sys.stderr, flush=True)

    got = {}
    with Pool(a.jobs, initializer=_worker_init) as pool:
        for ref, m in pool.imap_unordered(_one, refs, chunksize=4):
            got[ref] = m
    n_err = sum(1 for m in got.values() if m.get("hb_error"))

    tmp = a.output + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["design_id"] + COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            ref = (r.get(a.path_col) or "").strip()
            if ref in got:
                w.writerow({"design_id": r["design_id"], **got[ref]})
    os.replace(tmp, a.output)
    print(f"wrote {a.output}: {len(got)} structures, {n_err} errors", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
