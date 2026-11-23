#!/usr/bin/env python3
"""
nt_contacts.py  (runs on Euler)
===============================

Per-nucleotide protein->RNA heavy-atom contacts for the SAME design structures
the density renders are built from, full-atom.

    heavy-atom PAIRS within 4.5 A, counted per nucleotide of the 16-mer core

4.5 A and pair-counting to match the epitope precedent in the thesis figures
(Fig 2d's funnel and Fig 4g's hit footprint, both `geom.py`), not the 4.0 A of
the campaign's `n_heavy_contacts_4A`.

## Why this script exists at all

The design files three of the five stacks emit are BACKBONE ONLY (RFD3+LiMPNN,
BG+LiMPNN, BD1+LiMPNN ~ 4.0 atoms/residue; HD+LiMPNN and BG+BIF are full-atom).
Counting "heavy-atom contacts" straight off them is therefore not a measurement
of the design, it is a measurement of what each generator happened to write to
disk - the union fingerprint behind Fig 2d/4e shows exactly that, 55 contacts
per gated RFD3 design against 277 for HalluDesign, a 5x split that falls on the
backbone-only / full-atom line and not on anything biological.

The fix already exists: `analyses/library_selection/packed_designs/` holds
LigandMPNN side-chain retro-packs of the three backbone-only stacks, sequence
frozen, backbone preserved (CA exact, N/C/O transplanted back). This reads those
for the stacks that need them and the design file itself for the ones that do
not, so every structure here is full-atom and the four groups are comparable.

Side chains are a model's guess either way; what the packed route preserves and
a refold does not is the design's own POSE, which is what the render shows.

    python3 nt_contacts.py --designs <dir> --packed <dir> --out nt_contacts.csv

Writes one row per design: design_id, group, source, n_res, atoms_per_res,
nt01..nt16.  `source` is `packed` or `design`; `atoms_per_res` is the check that
the structure really is full-atom, and anything under 4.5 is refused rather than
silently counted.
"""
import argparse
import csv
import importlib.util
import io
import os
import sys
import tarfile

import numpy as np
import zstandard as zstd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "_extract_binders", os.path.join(HERE, "extract_binders.py"))
_eb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_eb)

CUTOFF = 4.5
MIN_ATOMS_PER_RES = 4.5      # backbone alone is 4.0; anything at/below is refused
GROUPS = ("RFD3", "BoltzGen", "BoltzDesign", "HalluDesign")


def heavy(atoms):
    return atoms[np.array([str(e).upper() != "H" for e in atoms.element])]


def packed_index(packed_dir, wanted):
    """{design_id: cif bytes} for the wanted ids, in one pass over the archives."""
    out = {}
    archives = sorted(f for f in os.listdir(packed_dir) if f.endswith(".tar.zst"))
    dctx = zstd.ZstdDecompressor()
    for i, name in enumerate(archives, 1):
        with open(os.path.join(packed_dir, name), "rb") as fh:
            with dctx.stream_reader(fh) as r:
                blob = io.BytesIO(r.read())
        with tarfile.open(fileobj=blob, mode="r:") as tf:
            for m in tf.getmembers():
                did = os.path.basename(m.name)
                if not did.endswith(".cif"):
                    continue
                did = did[:-4]
                if did in wanted and did not in out:
                    out[did] = tf.extractfile(m).read()
        if i % 50 == 0 or i == len(archives):
            print("  scanned %d/%d archives, %d/%d found"
                  % (i, len(archives), len(out), len(wanted)), flush=True)
    return out


def per_nt(path_or_bytes, load_structure, split_chains, tmp):
    """(16 contact counts, n_res, atoms_per_res) for one full-atom complex."""
    if isinstance(path_or_bytes, bytes):
        with open(tmp, "wb") as fh:
            fh.write(path_or_bytes)
        path = tmp
    else:
        path = path_or_bytes
    atoms = load_structure(path)
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        raise ValueError("no protein/RNA pair")
    ids = _eb.core_residues(rna)
    if ids is None:
        raise ValueError("core 16-mer not found")
    prot = heavy(protein)
    n_res = len(set(int(r) for r in protein.res_id))
    apr = len(prot) / max(n_res, 1)
    if apr <= MIN_ATOMS_PER_RES:
        raise ValueError("backbone-only structure (%.2f atoms/residue)" % apr)
    tree = cKDTree(prot.coord)
    counts = np.zeros(16, dtype=int)
    for k, rid in enumerate(ids):
        sel = heavy(rna[rna.res_id == rid])
        if len(sel) == 0:
            continue
        counts[k] = sum(len(p) for p in tree.query_ball_point(sel.coord, CUTOFF))
    return counts, n_res, apr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", required=True, help="dir with <group>/<design_id>.cif")
    ap.add_argument("--packed", required=True, help="packed_designs dir (tar.zst)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--src", default=_eb.DEFAULT_SRC)
    ap.add_argument("--tmp", default="/tmp/nt_contacts_work.cif")
    args = ap.parse_args()

    load_structure, split_chains = _eb.load_helpers(args.src)

    members = {}
    for g in GROUPS:
        d = os.path.join(args.designs, g)
        members[g] = sorted(f[:-4] for f in os.listdir(d) if f.endswith(".cif"))
        print("%-12s %d designs" % (g, len(members[g])))

    # Which ones need a packed structure: try the design file first, and fall
    # back only when it turns out to be backbone-only. No stack is hardcoded.
    rows, need, kept = [], set(), {}
    for g in GROUPS:
        for did in members[g]:
            p = os.path.join(args.designs, g, did + ".cif")
            try:
                counts, n_res, apr = per_nt(p, load_structure, split_chains, args.tmp)
            except ValueError:
                need.add(did)
                kept[did] = g
                continue
            rows.append([did, g, "design", n_res, "%.2f" % apr] + list(counts))
    print("full-atom straight from the design file: %d" % len(rows))
    print("needing the side-chain retro-pack: %d" % len(need))

    if need:
        blobs = packed_index(args.packed, need)
        missing = sorted(need - set(blobs))
        print("packed found %d, missing %d" % (len(blobs), len(missing)))
        for did, blob in blobs.items():
            try:
                counts, n_res, apr = per_nt(blob, load_structure, split_chains, args.tmp)
            except ValueError as e:
                print("  SKIP %s: %s" % (did, e))
                continue
            rows.append([did, kept[did], "packed", n_res, "%.2f" % apr] + list(counts))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design_id", "group", "source", "n_res", "atoms_per_res"]
                   + ["nt%02d" % i for i in range(1, 17)])
        w.writerows(rows)
    print("wrote %s (%d designs)" % (args.out, len(rows)))
    for g in GROUPS:
        sel = [r for r in rows if r[1] == g]
        if sel:
            v = np.array([r[5:] for r in sel], dtype=float)
            print("%-12s n=%4d  mean total %6.1f  per-nt %s"
                  % (g, len(sel), v.sum(1).mean(),
                     " ".join("%.0f" % x for x in v.mean(0))))


if __name__ == "__main__":
    main()
