#!/usr/bin/env python3
"""Export N singleton structures per stack as PDB files.

A singleton here is a backbone with no reciprocal structural neighbour above
TM 0.50 among the 36,907 passing designs - structurally unique within this
filtered set. Note that says nothing about quality: every one of these passed
DockQ >= 0.7 and RMSD <= 3 on at least one backend.

Chains are normalised to protein = A, RNA = B, since the generators disagree
(RFD3 design files put the protein in chain B).
"""
from __future__ import annotations

import argparse
import io as _io
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import biotite.structure.io.pdb as pdbio
from binder_map.chains import load_structure, split_chains

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", type=Path, default=Path("local_results"))
    ap.add_argument("--fine", type=Path, default=Path("local_results_fine"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--out", type=Path, default=Path("dist/singletons"))
    args = ap.parse_args()

    nodes = pd.read_csv(args.coarse / args.layer / "dedup" / "nodes.tsv", sep="\t")
    members = pd.read_csv(args.coarse / args.layer / "dedup" / "members.tsv", sep="\t")
    nodes = nodes.merge(members, on="struct_name", how="left")
    meta = pd.read_parquet(args.build / "meta" / "designs.parquet").set_index("design_id")
    stack = meta.reindex(nodes["struct_name"])["stack"].to_numpy()
    f = pd.read_parquet(args.fine / args.layer / "clusters" / "fold" /
                        "clusters.parquet").sort_values("node_id")[args.resolution].to_numpy()
    sizes = pd.Series(f).value_counts()
    single = np.array([sizes[x] == 1 for x in f])

    rng = np.random.default_rng(42)
    picks = []
    for st in ["RFD3+LiMPNN", "BG+LiMPNN", "BG+BIF", "BD1+LiMPNN", "HD+LiMPNN"]:
        cand = np.flatnonzero((stack == st) & single)
        if len(cand) == 0:
            print(f"  {st}: no singletons")
            continue
        take = rng.choice(cand, size=min(args.n, len(cand)), replace=False)
        # Spread rather than cluster: sort so the sample is reproducible.
        for k, node in enumerate(sorted(int(x) for x in take)):
            picks.append({"stack": st, "k": k, "node": node,
                          "design_id": nodes["struct_name"].iloc[node]})
        print(f"  {st:12s} {len(cand):>4,} singletons available, taking {len(take)}")

    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    ref_of = (manifest.assign(r=manifest["archive_path"] + "!" + manifest["member_path"])
              .set_index("design_id")["r"])
    listing = "\n".join(ref_of[p["design_id"]].replace("!", "\t") for p in picks)

    remote = f"""
set -u
cd {RB} || exit 1
WORK=$(mktemp -d); n=0
while IFS=$'\\t' read -r arch member; do
  [ -n "$arch" ] || continue
  tar --zstd -xOf "$arch" "$member" > "$WORK/$(printf '%03d' $n).cif" 2>/dev/null
  n=$((n+1))
done
cd "$WORK" && tar -cf - .
rm -rf "$WORK"
"""
    print(f"\nfetching {len(picks)} structures ...")
    proc = subprocess.run(["ssh", "-o", "ConnectTimeout=60", "euler", remote],
                          input=listing.encode(), capture_output=True)
    if proc.returncode != 0:
        print("FATAL: fetch failed")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    written = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with tarfile.open(fileobj=_io.BytesIO(proc.stdout)) as tf:
            tf.extractall(tmp)
        for n, p in enumerate(picks):
            src = tmp / f"{n:03d}.cif"
            if not src.exists() or src.stat().st_size < 200:
                continue
            protein, rna = split_chains(load_structure(src))
            protein, rna = protein.copy(), rna.copy()
            protein.chain_id[:] = "A"
            rna.chain_id[:] = "B"
            tag = p["stack"].replace("+", "-")
            dest = args.out / f"{tag}_{p['k']+1}_{p['design_id'][:44]}.pdb"
            pf = pdbio.PDBFile()
            pf.set_structure(protein + rna)
            pf.write(str(dest))
            p["pdb"] = dest.name
            written.append(p)

    (args.out / "singletons.json").write_text(json.dumps(written, indent=1))
    # A ready-to-run ChimeraX script, all superposed on the RNA target.
    lines = ["set bgColor white", "graphics silhouettes true width 1.2",
             "lighting soft"]
    COL = {"RFD3+LiMPNN": "#4a7ebb", "BG+LiMPNN": "#e07b39", "BG+BIF": "#4f9d69",
           "BD1+LiMPNN": "#b5514f", "HD+LiMPNN": "#8264a5"}
    for k, p in enumerate(written, start=1):
        lines.append(f"open {(args.out / p['pdb']).resolve()}")
        if k > 1:
            lines.append(f"matchmaker #{k}/B to #1/B")
        lines += [f"hide #{k} atoms", f"show #{k} cartoon",
                  f"color #{k}/A {COL.get(p['stack'], '#888888')}"]
        if k == 1:
            lines += [f"nucleotides #{k}/B tube/slab", f"color #{k}/B #2a9d8f"]
        else:
            lines.append(f"hide #{k}/B cartoon")
        lines.append(f"rename #{k} {p['stack'].replace('+','-')}_{k}")
    lines += ["view", "zoom 0.85"]
    (args.out / "open_all.cxc").write_text("\n".join(lines) + "\n")

    print(f"  wrote {len(written)} PDBs -> {args.out}")
    print(f"  ChimeraX: open {args.out / 'open_all.cxc'}  (all RNA-aligned)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
