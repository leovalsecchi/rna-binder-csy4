#!/usr/bin/env python3
"""Extract the full complex for in-browser display: protein backbone + RNA target.

Why more than CA:
  * 3Dmol cannot draw a real cartoon from alpha carbons alone - it needs the
    carbonyl to orient the ribbon - so a CA-only page can only ever show a
    featureless tube, with no helices or sheets.
  * A binder shown without its RNA target is close to useless on an epitope
    map, which is precisely a map of where the protein touches the RNA.

What is stored, per structure:
  protein  N, CA, C, O          80 res x 4 = 320 atoms
  RNA      P, C4', C1'          16 nt  x 3 =  48 atoms
  sse      one byte per residue (H / E / C), from biotite's P-SEA annotation

Encoding matches stage7b: 0.1 A grid, int16, delta-encoded along the atom
order so neighbouring atoms cost a few bits instead of a full coordinate.

Usage:  python scripts/stage7c_complex_traces.py --layer design
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import biotite.structure as struc

from binder_map.chains import load_structure, split_chains
from binder_map.rnaframe import canonical_positions

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
FM = RB / "analyses" / "foldseek_map"

PROT_ATOMS = ("N", "CA", "C", "O")
RNA_ATOMS = ("P", "C4'", "C1'")
N_RES, N_NT = 80, 16
N_PROT = N_RES * len(PROT_ATOMS)      # 320
N_RNA = N_NT * len(RNA_ATOMS)         # 48
N_ATOMS = N_PROT + N_RNA              # 368
SCALE = 10.0                          # 0.1 A grid


def ordered_atoms(chain, names, per_unit, n_units, positions=None):
    """(n_units * len(names), 3) coords, NaN where an atom is absent."""
    out = np.full((n_units * len(names), 3), np.nan, dtype=np.float64)
    if positions is None:                       # protein: order of appearance
        keys = list(zip(chain.chain_id, chain.res_id, chain.ins_code))
        seen, ordinal = {}, np.empty(len(keys), dtype=int)
        for i, k in enumerate(keys):
            if k not in seen:
                seen[k] = len(seen)
            ordinal[i] = seen[k]
    else:
        ordinal = positions - 1                 # RNA: canonical index 1..16
    names_arr = np.array([str(a).strip() for a in chain.atom_name])
    for j, want in enumerate(names):
        sel = names_arr == want
        for idx, coord in zip(ordinal[sel], chain.coord[sel]):
            if 0 <= idx < n_units:
                out[idx * len(names) + j] = coord
    return out


def sse_string(protein) -> str:
    """H/E/C per residue via P-SEA, which works from CA positions alone."""
    try:
        codes = struc.annotate_sse(protein)
    except Exception:
        return "C" * N_RES
    table = {"a": "H", "b": "E", "c": "C"}
    out = "".join(table.get(str(c), "C") for c in codes)
    return (out + "C" * N_RES)[:N_RES]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    ap.add_argument("--build", type=Path, default=FM / "build")
    ap.add_argument("--results", type=Path, default=FM / "results")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()

    out = args.results / args.layer
    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    want = {n: i for i, n in enumerate(nodes["struct_name"])}
    n = len(nodes)
    print(f"{args.layer}: {n:,} complexes to extract")

    # Only the archives holding a canonical structure...
    manifest = manifest[manifest["design_id"].isin(want)]
    # ...and only this shard's slice of them. Whole archives per shard, so no
    # archive is opened twice. Single-threaded this walks ~8,900 archives and
    # takes hours; stage 1 does the same traversal 64-way in ten minutes.
    archives = sorted(manifest["archive_path"].unique())
    mine = set(archives[args.shard::args.n_shards])
    manifest = manifest[manifest["archive_path"].isin(mine)]
    print(f"  shard {args.shard}/{args.n_shards}: "
          f"{manifest['archive_path'].nunique():,} archives, "
          f"{len(manifest):,} members")

    coords = np.zeros((n, N_ATOMS, 3), dtype=np.int16)
    sse = np.full((n, N_RES), ord("C"), dtype=np.uint8)
    found = np.zeros(n, dtype=bool)

    with tempfile.TemporaryDirectory() as tmpdir:
        scratch = Path(tmpdir)
        for archive_path, group in manifest.groupby("archive_path", sort=True):
            archive = RB / archive_path
            if not archive.exists():
                continue
            members = group["member_path"].tolist()
            work = scratch / "x"
            if work.exists():
                for f in work.rglob("*"):
                    if f.is_file():
                        f.unlink()
            work.mkdir(parents=True, exist_ok=True)
            subprocess.run(["tar", "--zstd", "-x", "-C", str(work), "-f",
                            str(archive)] + members,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for _, row in group.iterrows():
                idx = want.get(row["design_id"])
                if idx is None or found[idx]:
                    continue
                path = work / row["member_path"]
                if not path.exists():
                    continue
                try:
                    protein, rna = split_chains(load_structure(path))
                    p = ordered_atoms(protein, PROT_ATOMS, len(PROT_ATOMS), N_RES)
                    r = ordered_atoms(rna, RNA_ATOMS, len(RNA_ATOMS), N_NT,
                                      canonical_positions(rna))
                    allc = np.vstack([p, r])
                    if np.isnan(allc[:N_PROT]).all(axis=1).any():
                        continue                      # incomplete protein backbone
                    # Missing atoms take the previous atom's position, so a gap
                    # is a zero-length step rather than a wild spike.
                    bad = np.isnan(allc).any(axis=1)
                    if bad.any():
                        good_idx = np.where(~bad)[0]
                        if len(good_idx) == 0:
                            continue
                        for b in np.where(bad)[0]:
                            allc[b] = allc[good_idx[np.argmin(np.abs(good_idx - b))]]
                    centre = allc[:N_PROT].mean(axis=0)
                    q = np.clip(np.round((allc - centre) * SCALE), -32767, 32767)
                    q[1:] -= q[:-1].copy()
                    coords[idx] = q.astype(np.int16)
                    sse[idx] = np.frombuffer(sse_string(protein).encode(),
                                             dtype=np.uint8)
                    found[idx] = True
                except Exception:
                    continue

    print(f"  extracted {int(found.sum()):,} / {n:,}")
    if (miss := int((~found).sum())):
        print(f"  WARNING: {miss:,} complexes could not be extracted")

    dest = out / "complex"
    dest.mkdir(parents=True, exist_ok=True)
    if args.n_shards == 1:
        cpacked = gzip.compress(coords.tobytes(), 6)
        spacked = gzip.compress(sse.tobytes(), 6)
        (dest / "complex_int16.bin.gz").write_bytes(cpacked)
        (dest / "sse.bin.gz").write_bytes(spacked)
        print(f"  coords {len(cpacked)/1e6:.1f} MB gz (raw {coords.nbytes/1e6:.1f} MB)")
    else:
        # Per-shard partials; only this shard's rows are populated, so the
        # merge is an OR over `found` rather than a concatenation.
        part = dest / "parts"
        part.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(part / f"shard_{args.shard:02d}.npz",
                            coords=coords, sse=sse, found=found)
        print(f"  wrote partial -> {part / f'shard_{args.shard:02d}.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
