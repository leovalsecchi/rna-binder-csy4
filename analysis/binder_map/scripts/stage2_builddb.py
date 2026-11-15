#!/usr/bin/env python3
"""Stage 2: elect one canonical structure per coordinate hash.

Stage 1 dedups within a shard; sequence variants of one backbone share an
archive so most redundancy dies there. This stage catches hashes that recur
across shards, which naming would never reveal.

Usage:  python scripts/stage2_builddb.py --layer esm
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    files = sorted((out / "hashes").glob("shard_*.tsv"))
    if not files:
        print("FATAL: no hash tables; did stage 1 run?")
        return 1

    frames = [pd.read_csv(f, sep="\t") for f in files]
    hashes = pd.concat(frames, ignore_index=True)
    print(f"{args.layer}: {len(hashes):,d} structures, "
          f"{hashes['coord_hash'].nunique():,d} unique backbones")

    tally = hashes.groupby("protein_chain").size().to_dict()
    print(f"  protein chain tally: {tally}")

    # Elect the lexicographically first design_id per hash: deterministic, and
    # independent of shard order.
    hashes = hashes.sort_values(["coord_hash", "design_id"], kind="stable")
    canonical = hashes.groupby("coord_hash", sort=True)["design_id"].first()
    canonical = canonical.rename("struct_name").reset_index()

    members = (hashes.groupby("coord_hash", sort=True)["design_id"]
               .apply(lambda s: ",".join(sorted(s)))
               .rename("design_ids").reset_index())
    members = canonical.merge(members, on="coord_hash")

    nodes = canonical.copy()
    nodes["n_designs"] = members["design_ids"].str.count(",") + 1
    nodes.insert(0, "node_id", range(len(nodes)))

    dedup = out / "dedup"
    dedup.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(dedup / "canonical.tsv", sep="\t", index=False)
    members[["struct_name", "design_ids"]].to_csv(
        dedup / "members.tsv", sep="\t", index=False)
    nodes[["node_id", "struct_name", "n_designs"]].to_csv(
        dedup / "nodes.tsv", sep="\t", index=False)

    # Foldseek createsubdb takes a plain list of the names to keep.
    (dedup / "keep_names.txt").write_text(
        "\n".join(f"{n}.pdb" for n in nodes["struct_name"]) + "\n")

    print(f"  wrote {len(nodes):,d} nodes -> {dedup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
