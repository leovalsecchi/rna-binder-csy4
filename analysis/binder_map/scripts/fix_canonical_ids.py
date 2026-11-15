#!/usr/bin/env python3
"""Rewrite canonical.parquet with design_id restored, without re-running stage 4.

stage4_graph wrote this table with `index=False` while design_id was the index,
so the key vanished. Everything else in the file is correct and expensive to
recompute, so patch the column back rather than redoing the graph.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--layer", required=True)
    ap.add_argument("--results", type=Path, default=FM / "results")
    args = ap.parse_args()

    out = args.results / args.layer
    path = out / "fp" / "canonical.parquet"
    df = pd.read_parquet(path)
    if "design_id" in df.columns:
        print(f"{args.layer}: already has design_id, nothing to do")
        return 0

    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    if len(nodes) != len(df):
        print(f"{args.layer}: FATAL length mismatch {len(nodes)} vs {len(df)}")
        return 1
    # Row order is node_id order in both, which stage 4 guarantees.
    if not (df["node_id"].to_numpy() == nodes["node_id"].to_numpy()).all():
        print(f"{args.layer}: FATAL node_id order differs")
        return 1

    df.insert(1, "design_id", nodes["struct_name"].to_numpy())
    df.to_parquet(path, index=False)
    print(f"{args.layer}: restored design_id on {len(df):,} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
