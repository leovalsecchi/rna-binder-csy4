#!/usr/bin/env python3
"""Open 5 designs from a coarse bucket and 5 from a fine bucket, side by side.

The coarse five are deliberately drawn from FIVE DIFFERENT fine sub-clusters,
because that is exactly what the coarse bucket lumps together - sampling them
at random could accidentally pick five near-identical structures and hide the
very heterogeneity being examined.

All ten are superposed on the RNA target, so differences are in the binder.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import biotite.structure.io.pdb as pdbio
from binder_map.chains import load_structure, split_chains

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CX = "/Applications/ChimeraX-1.11.1.app/Contents/bin/ChimeraX"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", type=Path, default=Path("local_results"))
    ap.add_argument("--fine", type=Path, default=Path("local_results_fine"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--out", type=Path, default=Path("dist/bucket_compare"))
    args = ap.parse_args()

    nodes = pd.read_csv(args.coarse / args.layer / "dedup" / "nodes.tsv", sep="\t")
    c = pd.read_parquet(args.coarse / args.layer / "clusters" / "fold" /
                        "clusters.parquet").sort_values("node_id")[args.resolution].to_numpy()
    f = pd.read_parquet(args.fine / args.layer / "clusters" / "fold" /
                        "clusters.parquet").sort_values("node_id")[args.resolution].to_numpy()

    big = pd.Series(c).value_counts().index[0]
    inside = pd.Series(f[c == big]).value_counts()
    fine_id = int(inside.index[0])

    # Coarse set: one member from each of the n largest fine sub-clusters, so
    # the sample spans what the coarse bucket merged.
    coarse_pick = []
    for sub in inside.index[:args.n]:
        cand = np.flatnonzero((c == big) & (f == sub))
        coarse_pick.append(int(cand[len(cand) // 2]))
    # Fine set: evenly spaced through the single fine cluster.
    cand = np.flatnonzero(f == fine_id)
    fine_pick = [int(cand[i]) for i in
                 np.linspace(0, len(cand) - 1, args.n).astype(int)]

    print(f"coarse cluster {big} ({int((c==big).sum()):,} backbones) "
          f"-> {len(inside)} fine clusters")
    print(f"  COARSE set: 5 members from 5 different fine sub-clusters "
          f"{list(inside.index[:args.n])}")
    print(f"  FINE set:   5 members of fine cluster {fine_id} "
          f"({int(inside.iloc[0]):,} backbones)")

    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    ref_of = (manifest.assign(r=manifest["archive_path"] + "!" + manifest["member_path"])
              .set_index("design_id")["r"])

    args.out.mkdir(parents=True, exist_ok=True)
    jobs = ([("COARSE", i) for i in coarse_pick] + [("FINE", i) for i in fine_pick])
    listing = "\n".join(
        ref_of[nodes["struct_name"].iloc[i]].replace("!", "\t") for _, i in jobs)

    remote = f"""
set -u
cd {RB} || exit 1
WORK=$(mktemp -d); n=0
while IFS=$'\\t' read -r arch member; do
  [ -n "$arch" ] || continue
  tar --zstd -xOf "$arch" "$member" > "$WORK/$(printf '%02d' $n).cif" 2>/dev/null
  n=$((n+1))
done
cd "$WORK" && tar -cf - .
rm -rf "$WORK"
"""
    proc = subprocess.run(["ssh", "-o", "ConnectTimeout=60", "euler", remote],
                          input=listing.encode(), capture_output=True)
    if proc.returncode != 0:
        print("FATAL: fetch failed"); return 1

    import io as _io, tarfile
    paths = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with tarfile.open(fileobj=_io.BytesIO(proc.stdout)) as tf:
            tf.extractall(tmp)
        for n, (group, node) in enumerate(jobs):
            src = tmp / f"{n:02d}.cif"
            if not src.exists():
                continue
            protein, rna = split_chains(load_structure(src))
            protein, rna = protein.copy(), rna.copy()
            protein.chain_id[:] = "A"
            rna.chain_id[:] = "B"
            name = nodes["struct_name"].iloc[node]
            dest = args.out / f"{group}_{n:02d}.pdb"
            pf = pdbio.PDBFile(); pf.set_structure(protein + rna); pf.write(str(dest))
            paths.append((group, dest, name))

    # Warm colours for the coarse set, cool for the fine set: at a glance the
    # question is whether the cool group is visibly tighter than the warm one.
    WARM = ["#b5514f", "#e07b39", "#d1495b", "#c1663a", "#a33b3b"]
    COOL = ["#4a7ebb", "#2a9d8f", "#3f6f9c", "#2f8f83", "#5a93cc"]
    lines = ["set bgColor white", "graphics silhouettes true width 1.2",
             "lighting soft"]
    for k, (group, path, name) in enumerate(paths, start=1):
        lines.append(f"open {path.resolve()}")
        if k > 1:
            lines.append(f"matchmaker #{k}/B to #1/B")
        lines += [f"hide #{k} atoms", f"show #{k} cartoon"]
        col = (WARM if group == "COARSE" else COOL)[(k - 1) % 5]
        lines.append(f"color #{k}/A {col}")
        # only the reference RNA is shown; ten superposed copies is noise
        lines.append(f"color #{k}/B #999999")
        if k > 1:
            lines.append(f"hide #{k}/B cartoon")
        else:
            lines.append(f"nucleotides #{k}/B tube/slab")
            lines.append(f"color #{k}/B #2a9d8f")
        lines.append(f"rename #{k} {group}_{name[:34]}")
    lines += ["view", "zoom 0.85",
              f"save {(args.out / 'compare.png').resolve()} width 1500 "
              "height 1150 supersample 3"]
    script = args.out / "compare.cxc"
    script.write_text("\n".join(lines) + "\n")

    print(f"\nopening {len(paths)} structures in ChimeraX ...")
    subprocess.Popen([CX, "--cmd", f"open {script.resolve()}"])
    print(f"  warm = coarse bucket {big} (5 different fine sub-clusters)")
    print(f"  cool = fine bucket {fine_id} (one tight cluster)")
    print(f"  script: {script}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
