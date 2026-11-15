#!/usr/bin/env python3
"""opth_metrics.py — paired H-optimisation experiment.

For every structure in a chunk CSV, computes the interface H-bond count and the
buried-unsatisfied-polars count TWICE on the same pose:

  *_raw   production settings: hydrogens as Rosetta builds them from the
          heavy-atom CIF, never optimised (packing:no_optH is true by default)
  *_opth  after core::pack::optimizeH(), i.e. polar-H torsions
          optimised and Asn/Gln/His flipped (-flip_HNQ)

and additionally the FIELD-STANDARD BUNS via BuriedUnsatHbondFilter
(ddG-style, VSASA burial, probe 1.1 A) on both H states, so the legacy
InterfaceAnalyzer number can be compared to the metric the binder-design
field actually reports.

Runs INSIDE pyrosetta.sif. Stdlib-only IO.

  python3 opth_metrics.py --csv chunk.csv --path-col af3_msaoff_cif_path \
      --out chunk_out.csv [--jobs 8] [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive_read import resolve

csv.field_size_limit(10 ** 7)

OUT_COLS = [
    "hb_raw", "hb_opth",
    "buns_ia_raw", "buns_ia_opth",
    "buns_vsasa_raw", "buns_vsasa_opth",
    "dsasa_total_raw", "n_polar_h_moved", "opth_ok",
]

_SF = None


def _worker_init():
    global _SF
    import pyrosetta
    # production flags + -flip_HNQ so optimize_H can flip Asn/Gln/His.
    # -flip_HNQ only affects the optH routine; the *_raw pass is computed before
    # optimize_H is ever called, so it stays identical to production.
    pyrosetta.init("-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false "
                   "-flip_HNQ")
    _SF = pyrosetta.get_fa_scorefxn()


def _iface_string(pose):
    info = pose.pdb_info()
    cp, cr = set(), set()
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        ch = info.chain(i)
        if r.is_protein():
            cp.add(ch)
        elif r.is_RNA() or r.is_DNA():
            cr.add(ch)
    if not cp or not cr:
        return None
    return "".join(sorted(cp)) + "_" + "".join(sorted(cr))


def _count_iface_hbonds(pose):
    from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
    _SF(pose)
    pose.update_residue_neighbors()
    hb = HBondSet()
    fill_hbond_set(pose, False, hb)
    n = 0
    for k in range(1, hb.nhbonds() + 1):
        b = hb.hbond(k)
        rd, ra = pose.residue(b.don_res()), pose.residue(b.acc_res())
        dn = rd.is_RNA() or rd.is_DNA()
        an = ra.is_RNA() or ra.is_DNA()
        if (rd.is_protein() and an) or (ra.is_protein() and dn):
            n += 1
    return n


def _ia(pose, iface):
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    iam = InterfaceAnalyzerMover(iface)
    iam.set_compute_packstat(False)
    iam.set_compute_interface_sc(False)
    iam.set_pack_separated(False)
    iam.set_scorefunction(_SF)
    iam.apply(pose)
    d = iam.get_all_data()
    return int(d.delta_unsat_hbonds), float(list(d.dSASA)[0])


def _buns_vsasa(pose):
    """Modern-style BuriedUnsatHbonds: VSASA burial, probe radius 1.1 A, all
    heavy-atom unsats. This is the burial definition the binder-design field
    uses, in contrast to InterfaceAnalyzer's legacy absolute-SASA < 0.3 A^2.

    Two deliberate deviations from the published `buns_heavy_ball_1.1D`:
      * DAlphaBall is off - the binary is not in the container, and
        set_dalphaball_sasa() takes no argument in this build (calling it would
        ENABLE it).
      * ddG-style is off - set_use_ddG_style(True) returns 0.0 for every
        protein-RNA pose in this PyRosetta build (verified on several), so this
        is a whole-pose ABSOLUTE count, not a bound-minus-unbound delta. Use it
        for the optH delta and for cross-design comparison, NOT as a substitute
        for delta_unsat_hbonds."""
    from pyrosetta.rosetta.protocols.simple_filters import BuriedUnsatHbondFilter
    f = BuriedUnsatHbondFilter()
    f.set_use_vsasa(True)
    f.set_probe_radius(1.1)
    f.set_jump_num(1)
    f.set_report_all_heavy_atom_unsats(True)
    f.set_ignore_surface_res(False)
    return float(f.report_sm(pose))


def _one(row_path):
    idx, path = row_path
    import pyrosetta
    from pyrosetta.rosetta.core.pack import optimizeH
    out = {c: "" for c in OUT_COLS}
    try:
        pose = pyrosetta.pose_from_file(resolve(path))
        iface = _iface_string(pose)
        if iface is None:
            return idx, out
        # ---- raw (production) ----
        out["hb_raw"] = _count_iface_hbonds(pose)
        b, sa = _ia(pose, iface)
        out["buns_ia_raw"] = b
        out["dsasa_total_raw"] = round(sa, 2)
        try:
            out["buns_vsasa_raw"] = _buns_vsasa(pose)
        except Exception:
            pass
        # ---- optimised hydrogens ----
        p2 = pyrosetta.Pose(pose)
        xyz_before = [(i, j, p2.residue(i).xyz(j))
                      for i in range(1, p2.total_residue() + 1)
                      for j in range(1, p2.residue(i).natoms() + 1)
                      if p2.residue(i).atom_is_polar_hydrogen(j)]
        optimizeH(p2, _SF)
        moved = sum(1 for i, j, v in xyz_before
                    if (p2.residue(i).xyz(j) - v).norm() > 0.1)
        out["n_polar_h_moved"] = moved
        out["hb_opth"] = _count_iface_hbonds(p2)
        b2, _ = _ia(p2, iface)
        out["buns_ia_opth"] = b2
        try:
            out["buns_vsasa_opth"] = _buns_vsasa(p2)
        except Exception:
            pass
        out["opth_ok"] = 1
        return idx, out
    except Exception as e:
        out["opth_ok"] = 0
        return idx, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--path-col", default="af3_msaoff_cif_path")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--progress", action="store_true")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.csv)))
    if a.limit:
        rows = rows[:a.limit]
    jobs = [(i, r[a.path_col]) for i, r in enumerate(rows)]
    print(f"{len(jobs)} structures, {a.jobs} workers", flush=True)
    if not rows:
        # a trailing chunk can be header-only when N does not divide the manifest;
        # emit the header so the merge's header-consistency check still passes
        with open(a.csv) as fh, open(a.out, "w") as out:
            out.write(fh.readline().rstrip("\n") + "," + ",".join(OUT_COLS) + "\n")
        print(f"wrote {a.out}: 0 rows (empty chunk)", flush=True)
        return

    done = 0
    with Pool(a.jobs, initializer=_worker_init) as pool:
        for idx, out in pool.imap_unordered(_one, jobs, chunksize=1):
            rows[idx].update(out)
            done += 1
            if a.progress and done % 25 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)

    fields = list(rows[0].keys())
    for c in OUT_COLS:
        if c not in fields:
            fields.append(c)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if str(r.get("opth_ok")) == "1")
    print(f"wrote {a.out}: {len(rows)} rows, {ok} ok", flush=True)


if __name__ == "__main__":
    main()
