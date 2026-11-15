#!/usr/bin/env python3
"""Fetch the real CIFs for the deck and normalise them for rendering.

The compressed in-page blob carries three atoms per nucleotide, which is enough
for a tube but not for a nucleic cartoon - ChimeraX draws the RNA as loose dots
from it. A deck is only a few hundred structures, so fetch the actual files and
render those instead.

Normalising matters too: chain IDs are not consistent across generators (RFD3
design files put the protein in chain B), so every structure is rewritten with
protein = A and RNA = B. The ChimeraX script can then treat them uniformly and
superpose on chain B without special-casing.

Usage:
    python scripts/fetch_deck_structures.py --layer design --axis fold \
        --resolution cluster_r10 --min-size 4 --per-cluster 4
"""
from __future__ import annotations

import argparse
import json
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("local_results"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--axis", default="fold")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--min-size", type=int, default=4)
    ap.add_argument("--per-cluster", type=int, default=4)
    ap.add_argument("--out", type=Path, default=Path("dist/deck_structures"))
    args = ap.parse_args()

    base = args.results / args.layer
    nodes = pd.read_csv(base / "dedup" / "nodes.tsv", sep="\t")
    members = pd.read_csv(base / "dedup" / "members.tsv", sep="\t")
    nodes = nodes.merge(members, on="struct_name", how="left")
    nodes["n_designs"] = nodes["n_designs"].fillna(1).astype(int)
    meta = pd.read_parquet(args.build / "meta" / "designs.parquet").set_index("design_id")
    stack = meta.reindex(nodes["struct_name"])["stack"].to_numpy()
    labels = pd.read_parquet(
        base / "clusters" / args.axis / "clusters.parquet"
    ).sort_values("node_id")[args.resolution].to_numpy()
    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    ref_of = (manifest.assign(ref=manifest["archive_path"] + "!" + manifest["member_path"])
              .set_index("design_id")["ref"])

    sizes = pd.Series(labels).value_counts()
    clusters = [c for c, k in sizes.items() if k >= args.min_size]
    clusters.sort(key=lambda c: -sizes[c])

    plan, chosen = [], []
    for cid in clusters:
        idx = np.flatnonzero(labels == cid)
        pick = idx[np.linspace(0, len(idx) - 1,
                               min(args.per_cluster, len(idx))).astype(int)]
        for slot, node in enumerate(pick):
            name = nodes["struct_name"].iloc[node]
            plan.append({
                "cluster": int(cid), "slot": slot, "node": int(node),
                "design_id": name, "stack": str(stack[node]),
                "ref": str(ref_of.get(name, "")),
            })
        chosen.append({
            "cluster": int(cid),
            "backbones": int(sizes[cid]),
            "designs": int(nodes["n_designs"].to_numpy()[idx].sum()),
            "pct": 100.0 * sizes[cid] / len(nodes),
            "stacks": pd.Series(stack[idx]).value_counts().to_dict(),
        })
    print(f"{len(clusters)} clusters, {len(plan)} structures to fetch")

    args.out.mkdir(parents=True, exist_ok=True)
    by_archive: dict[str, list[dict]] = {}
    for item in plan:
        if not item["ref"]:
            continue
        arch, member = item["ref"].split("!", 1)
        item["member"] = member
        by_archive.setdefault(arch, []).append(item)
    print(f"  spread over {len(by_archive)} archives")

    # One ssh for the whole deck: a remote helper extracts every wanted
    # member into a scratch dir and streams it back as a single tar. Per-file
    # round-trips would be ~172 ssh handshakes for no benefit.
    want = [(i["ref"].split("!", 1)[0], i["ref"].split("!", 1)[1], i)
            for i in plan if i["ref"]]
    listing = "\n".join(f"{a}\t{m}" for a, m, _ in want)

    remote = f"""
set -u
RB={RB}
WORK=$(mktemp -d)
cd "$RB" || exit 1
n=0
while IFS=$'\\t' read -r arch member; do
  [ -n "$arch" ] || continue
  tar --zstd -xOf "$arch" "$member" > "$WORK/$(printf '%05d' $n).cif" 2>/dev/null
  n=$((n+1))
done
cd "$WORK" && tar -cf - . 
rm -rf "$WORK"
"""
    print(f"  requesting {len(want)} structures in one call...")
    proc = subprocess.run(["ssh", "-o", "ConnectTimeout=60", "euler", remote],
                          input=listing.encode(), capture_output=True)
    if proc.returncode != 0:
        print("FATAL: remote extraction failed")
        print(proc.stderr.decode()[:500])
        return 1

    ok = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        import tarfile
        with tarfile.open(fileobj=__import__("io").BytesIO(proc.stdout)) as tf:
            tf.extractall(tmp)
        for n, (_, _, item) in enumerate(want):
            src = tmp / f"{n:05d}.cif"
            if not src.exists() or src.stat().st_size < 200:
                continue
            try:
                protein, rna = split_chains(load_structure(src))
            except Exception:
                continue
            protein, rna = protein.copy(), rna.copy()
            # Normalise chain ids: RFD3 design files put the protein in B, and
            # the ChimeraX script superposes on chain B as the RNA.
            protein.chain_id[:] = "A"
            rna.chain_id[:] = "B"
            dest = args.out / f"c{item['cluster']:04d}_s{item['slot']}.pdb"
            f = pdbio.PDBFile()
            f.set_structure(protein + rna)
            f.write(str(dest))
            item["pdb"] = dest.name
            ok += 1

    meta_out = {
        "layer": args.layer, "axis": args.axis, "resolution": args.resolution,
        "n_backbones": int(len(nodes)),
        "n_designs": int(nodes["n_designs"].sum()),
        "n_clusters_total": int(len(sizes)),
        "clusters": chosen,
        "panels": [p for p in plan if "pdb" in p],
    }
    (args.out / "deck_meta.json").write_text(json.dumps(meta_out, indent=1))
    print(f"  fetched {ok}/{len(plan)} structures -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
