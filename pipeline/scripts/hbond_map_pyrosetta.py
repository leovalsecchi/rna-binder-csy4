#!/usr/bin/env python3
"""Protein-RNA H-bond contact MAP (amino-acid type x nucleotide position) via PyRosetta.

Runs INSIDE pyrosetta.sif (CIFs are heavy-atom only; Rosetta builds hydrogens and
does orientation-dependent H-bond detection — same engine as hbonds_pyrosetta.py).
For every interface H-bond (one partner protein, the other a nucleic acid) it records
(aa 1-letter type, nt position 1..N) where nt position = rank of the RNA residue along
the chain (5'->3'). Emits a LONG CSV: stack,refolder,msa_variant,abs_path,aa,nt_pos,n_hbonds
(one row per (structure, aa, nt_pos) with >=1 H-bond). Aggregate -> heatmap downstream.

Usage (inside pyrosetta.sif):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        pyrosetta.sif python3 code/scripts/hbond_map_pyrosetta.py \\
        --csv <structural.csv> --output <hbond_map.csv> [--limit N] [--jobs 8] --progress
"""
from __future__ import annotations
import argparse, csv, os, sys
from collections import defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive_read import resolve   # archive-aware (.tar.zst!member) path resolution

_SF = None


def _worker_init():
    global _SF
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false")
    _SF = pyrosetta.get_fa_scorefxn()


def _map_one(path):
    """Return (path, {(aa1, nt_pos): n_hbonds}) or (path, None) on failure."""
    import pyrosetta
    from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
    try:
        pose = pyrosetta.pose_from_file(resolve(path))
        _SF(pose)
        pose.update_residue_neighbors()
        # nucleic residue index -> nt position (rank along chain, 1-based)
        nuc = sorted(i for i in range(1, pose.total_residue() + 1)
                     if pose.residue(i).is_RNA() or pose.residue(i).is_DNA())
        ntpos = {ri: k + 1 for k, ri in enumerate(nuc)}
        hb = HBondSet()
        fill_hbond_set(pose, False, hb)
        counts: dict = defaultdict(int)
        for i in range(1, hb.nhbonds() + 1):
            b = hb.hbond(i)
            dr, ar = b.don_res(), b.acc_res()
            rd, ra = pose.residue(dr), pose.residue(ar)
            if rd.is_protein() and (ra.is_RNA() or ra.is_DNA()):
                counts[(rd.name1(), ntpos[ar])] += 1
            elif ra.is_protein() and (rd.is_RNA() or rd.is_DNA()):
                counts[(ra.name1(), ntpos[dr])] += 1
        return path, dict(counts)
    except Exception as e:
        return path, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="structural CSV (provides paths + stack/refolder/msa per path)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--path-col", default="abs_path")
    ap.add_argument("--limit", type=int, default=0, help="process only first N unique structures (sampling)")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    pc = args.path_col
    # path -> (stack, refolder, msa_variant) from the first row that has it
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
    runner = Pool(args.jobs, initializer=_worker_init).imap_unordered if args.jobs > 1 else None
    if runner:
        with Pool(args.jobs, initializer=_worker_init) as pool:
            for i, (path, m) in enumerate(pool.imap_unordered(_map_one, paths, chunksize=4)):
                results[path] = m; n_err += int(m is None)
                if args.progress and i % 200 == 0:
                    print(f"  {i}/{len(paths)}", file=sys.stderr)
    else:
        _worker_init()
        for i, path in enumerate(paths):
            _, m = _map_one(path); results[path] = m; n_err += int(m is None)
            if args.progress and i % 200 == 0:
                print(f"  {i}/{len(paths)}", file=sys.stderr)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stack", "refolder", "msa_variant", "abs_path", "aa", "nt_pos", "n_hbonds"])
        nrows = 0
        for path in paths:
            m = results.get(path)
            if not m:
                continue
            stack, refolder, msa = meta[path]
            for (aa, pos), n in sorted(m.items()):
                w.writerow([stack, refolder, msa, path, aa, pos, n]); nrows += 1
    print(f"Wrote {nrows} (structure,aa,nt) rows for {len(paths)-n_err}/{len(paths)} structures "
          f"({n_err} failed) -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
