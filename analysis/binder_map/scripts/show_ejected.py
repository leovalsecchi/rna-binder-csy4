#!/usr/bin/env python3
"""Open designs the stricter clustering KEPT together vs ones it EJECTED.

Coarse cluster 0 held 3,374 backbones. The refined run (TMalign scoring,
TM >= 0.50, mutual k-NN) kept 1,752 of them together and split the other 1,622
into 147 separate clusters, 112 of them singletons.

Everything is superposed on the RNA target, so any visible difference is in
the binder and not in how the pair happens to sit in space.
"""
from __future__ import annotations

import argparse
import io as _io
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
CX = "/Applications/ChimeraX-1.11.1.app/Contents/bin/ChimeraX"
KEPT = ["#2a6fb5", "#3f87c9", "#5a9bd6", "#2f6f9c", "#1d5c94"]      # blues
EJECT = ["#c0392b", "#e07b39", "#d1495b", "#a33b3b", "#e0a020"]     # warm


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", type=Path, default=Path("local_results"))
    ap.add_argument("--fine", type=Path, default=Path("local_results_fine"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-open", action="store_true")
    ap.add_argument("--coarse-cluster", type=int, default=None,
                    help="which coarse cluster to dissect (default: largest)")
    args = ap.parse_args()

    if args.out is None:
        tag = ("largest" if args.coarse_cluster is None
               else f"c{args.coarse_cluster}")
        args.out = Path(f"dist/ejected_{tag}")
    nodes = pd.read_csv(args.coarse / args.layer / "dedup" / "nodes.tsv", sep="\t")
    c = pd.read_parquet(args.coarse / args.layer / "clusters" / "fold" /
                        "clusters.parquet").sort_values("node_id")[args.resolution].to_numpy()
    f = pd.read_parquet(args.fine / args.layer / "clusters" / "fold" /
                        "clusters.parquet").sort_values("node_id")[args.resolution].to_numpy()

    big = (args.coarse_cluster if args.coarse_cluster is not None
           else int(pd.Series(c).value_counts().index[0]))
    members = np.flatnonzero(c == big)
    core = int(pd.Series(f[members]).value_counts().index[0])
    fine_sizes = pd.Series(f).value_counts()

    kept = members[f[members] == core]
    ejected = members[f[members] != core]
    singles = np.array([i for i in ejected if fine_sizes[f[i]] == 1], dtype=int)

    kept_pick = [int(kept[i]) for i in
                 np.linspace(0, len(kept) - 1, args.n).astype(int)]
    # Prefer singletons among the ejected - those are the structures the strict
    # criterion found NO reciprocal partner for, so they show the split most
    # clearly; top up from the largest ejected pieces if there are too few.
    ej_pick = [int(x) for x in singles[np.linspace(0, len(singles) - 1,
                                                   min(3, len(singles))).astype(int)]]
    others = pd.Series(f[ejected]).value_counts().index
    for sub in others:
        if len(ej_pick) >= args.n:
            break
        cand = ejected[f[ejected] == sub]
        pick = int(cand[len(cand) // 2])
        if pick not in ej_pick:
            ej_pick.append(pick)

    print(f"coarse cluster {big}: {len(members):,} backbones")
    print(f"  KEPT together (fine cluster {core}): {len(kept):,}")
    print(f"  EJECTED: {len(ejected):,} into {len(set(f[ejected]))} clusters, "
          f"{len(singles):,} singletons")
    jobs = ([("KEPT", i) for i in kept_pick] + [("EJECTED", i) for i in ej_pick])
    for g, i in jobs:
        tag = f"fine={f[i]}" + (" SINGLETON" if fine_sizes[f[i]] == 1 else
                                f" (size {fine_sizes[f[i]]})")
        print(f"    {g:8s} {nodes['struct_name'].iloc[i][:46]:46s} {tag}")

    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    ref_of = (manifest.assign(r=manifest["archive_path"] + "!" + manifest["member_path"])
              .set_index("design_id")["r"])
    listing = "\n".join(ref_of[nodes["struct_name"].iloc[i]].replace("!", "\t")
                        for _, i in jobs)

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
        print("FATAL: fetch failed")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    made = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with tarfile.open(fileobj=_io.BytesIO(proc.stdout)) as tf:
            tf.extractall(tmp)
        for n, (group, node) in enumerate(jobs):
            src = tmp / f"{n:02d}.cif"
            if not src.exists() or src.stat().st_size < 200:
                continue
            protein, rna = split_chains(load_structure(src))
            protein, rna = protein.copy(), rna.copy()
            protein.chain_id[:] = "A"
            rna.chain_id[:] = "B"
            dest = args.out / f"{n:02d}_{group}.pdb"
            pf = pdbio.PDBFile()
            pf.set_structure(protein + rna)
            pf.write(str(dest))
            made.append((group, dest, nodes["struct_name"].iloc[node], f[node],
                         int(fine_sizes[f[node]])))

    lines = ["set bgColor white", "graphics silhouettes true width 1.2",
             "lighting soft"]
    ki = ei = 0
    for k, (group, path, name, fid, fsize) in enumerate(made, start=1):
        lines.append(f"open {path.resolve()}")
        # EVERY model superposed on the reference RNA, so the target is the
        # fixed frame and all visible difference is in the binder.
        if k > 1:
            lines.append(f"matchmaker #{k}/B to #1/B")
        lines += [f"hide #{k} atoms", f"show #{k} cartoon"]
        if group == "KEPT":
            col = KEPT[ki % len(KEPT)]; ki += 1
        else:
            col = EJECT[ei % len(EJECT)]; ei += 1
        lines.append(f"color #{k}/A {col}")
        if k == 1:
            lines += [f"nucleotides #{k}/B tube/slab", f"color #{k}/B #2a9d8f"]
        else:
            lines.append(f"hide #{k}/B cartoon")     # 10 copies of one RNA is noise
        suffix = "SINGLETON" if fsize == 1 else f"n{fsize}"
        lines.append(f"rename #{k} {group}_f{fid}_{suffix}")
    lines += ["view", "zoom 0.8",
              f"save {(args.out / 'ejected_compare.png').resolve()} "
              "width 1600 height 1200 supersample 3"]
    script = args.out / "show.cxc"
    script.write_text("\n".join(lines) + "\n")

    print(f"\n  blue = kept in fine cluster {core}")
    print(f"  warm = ejected from the coarse cluster by the strict criterion")
    print(f"  all superposed on the RNA target (shown once, teal)")
    if not args.no_open:
        subprocess.Popen([CX, "--cmd", f"open {script.resolve()}"])
        print("  opening in ChimeraX ...")
    print(f"  script: {script}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
