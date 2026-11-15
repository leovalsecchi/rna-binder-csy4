#!/usr/bin/env python3
"""interface_aa_nt_map.py — protein↔RNA interface at (amino-acid type × RNA nt position)
resolution, for BOTH H-bonds and heavy-atom contacts, in ONE PyRosetta pass.

H-bonds via Rosetta HBondSet (same engine as hbond_map_pyrosetta.py); contacts = heavy-
atom pairs <=4.5 A (same 4.5 A definition as analysis/structural.py rna_contact). nt_pos
= rank of the RNA residue along the chain (1..N). Emits a LONG CSV:
    stack,refolder,msa_variant,abs_path,aa,nt_pos,n_hbonds,n_contacts
(one row per (structure, aa, nt_pos) with >=1 H-bond OR >=1 contact). Runs inside pyrosetta.sif.

Usage (inside pyrosetta.sif):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        pyrosetta.sif python3 interface_aa_nt_map.py --csv <in.csv> --output <map.csv> \\
        [--path-col abs_path] [--limit N] [--jobs 8] --progress
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from archive_read import resolve
except ImportError:  # allow running from plots_tools; archive_read lives in code/scripts
    sys.path.insert(0, "/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts")
    from archive_read import resolve

_SF = None
CONTACT_A = 4.5


def _worker_init():
    global _SF
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false")
    _SF = pyrosetta.get_fa_scorefxn()


def _map_one(path):
    """Return (path, ({(aa,nt): n_hbonds}, {(aa,nt): n_contacts})) or (path, None)."""
    import numpy as np
    import pyrosetta
    from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
    try:
        pose = pyrosetta.pose_from_file(resolve(path))
        _SF(pose)
        pose.update_residue_neighbors()
        nuc = sorted(i for i in range(1, pose.total_residue() + 1)
                     if pose.residue(i).is_RNA() or pose.residue(i).is_DNA())
        ntpos = {ri: k + 1 for k, ri in enumerate(nuc)}

        # --- H-bonds (aa, nt) ---
        hb = HBondSet()
        fill_hbond_set(pose, False, hb)
        hbonds: dict = defaultdict(int)
        for i in range(1, hb.nhbonds() + 1):
            b = hb.hbond(i)
            dr, ar = b.don_res(), b.acc_res()
            rd, ra = pose.residue(dr), pose.residue(ar)
            if rd.is_protein() and (ra.is_RNA() or ra.is_DNA()):
                hbonds[(rd.name1(), ntpos[ar])] += 1
            elif ra.is_protein() and (rd.is_RNA() or rd.is_DNA()):
                hbonds[(ra.name1(), ntpos[dr])] += 1

        # --- heavy-atom contacts <=4.5 A (aa, nt) ---
        pc, paa = [], []
        for i in range(1, pose.total_residue() + 1):
            r = pose.residue(i)
            if not r.is_protein():
                continue
            aa1 = r.name1()
            for a in range(1, r.nheavyatoms() + 1):
                x = r.xyz(a); pc.append((x.x, x.y, x.z)); paa.append(aa1)
        rc, rnt = [], []
        for j in nuc:
            r = pose.residue(j); ntp = ntpos[j]
            for a in range(1, r.nheavyatoms() + 1):
                x = r.xyz(a); rc.append((x.x, x.y, x.z)); rnt.append(ntp)
        contacts: dict = defaultdict(int)
        if pc and rc:
            P = np.asarray(pc); R = np.asarray(rc)
            d2 = ((P[:, None, :] - R[None, :, :]) ** 2).sum(-1)
            ii, jj = np.where(d2 <= CONTACT_A * CONTACT_A)
            paa = np.asarray(paa); rnt = np.asarray(rnt)
            for a, b in zip(paa[ii], rnt[jj]):
                contacts[(str(a), int(b))] += 1
        return path, (dict(hbonds), dict(contacts))
    except Exception:  # noqa: BLE001
        return path, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--path-col", default="abs_path")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    pc = args.path_col
    meta = {}
    for r in rows:
        p = r.get(pc, "")
        if p and p not in meta:
            meta[p] = (r.get("stack", ""), r.get("refolder", ""), r.get("msa_variant", ""))
    paths = sorted(meta)
    if args.limit:
        paths = paths[:args.limit]
    print(f"{len(rows)} rows, {len(paths)} unique structures, jobs={args.jobs}", file=sys.stderr)

    results = {}
    n_err = 0
    with Pool(args.jobs, initializer=_worker_init) as pool:
        for i, (path, m) in enumerate(pool.imap_unordered(_map_one, paths, chunksize=4)):
            results[path] = m; n_err += int(m is None)
            if args.progress and i % 200 == 0:
                print(f"  {i}/{len(paths)}", file=sys.stderr)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stack", "refolder", "msa_variant", "abs_path", "aa", "nt_pos", "n_hbonds", "n_contacts"])
        nrows = 0
        for path in paths:
            m = results.get(path)
            if not m:
                continue
            hbonds, contacts = m
            stack, refolder, msa = meta[path]
            for key in set(hbonds) | set(contacts):
                aa, nt = key
                w.writerow([stack, refolder, msa, path, aa, nt, hbonds.get(key, 0), contacts.get(key, 0)])
                nrows += 1
    print(f"wrote {nrows} (structure,aa,nt) rows -> {args.output}  (errors: {n_err}/{len(paths)})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
