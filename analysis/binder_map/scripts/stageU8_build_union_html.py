#!/usr/bin/env python3
"""Stage U8: bake the union fold map into one self-contained page.

Reads the four outputs of the union pipeline, joins them onto one row per
backbone, and hands them to `binder_map.unionhtml` to pack and bake.

The join is the delicate part. `embed/fold/embed.parquet`, `clusters/fold/
clusters.parquet` and `dedup/nodes.tsv` all key on `node_id`, but
`meta/layers.parquet` keys on `coord_hash`, and `dedup/canonical.tsv` is the
only bridge between them. `unionhtml.join_layers` does it with LEFT joins and
refuses to lose or multiply a row: an inner join here would answer a broken key
by quietly dropping nodes, and the page would render perfectly while wrong.

The per-cluster medoid backbones come from stage 7 (`--layer design --results
results_union`), which bundles every medoid across every resolution. Only the
operative resolution's medoids are embedded; see `--no-medoids` to build the
points-only page.

Usage:
    python scripts/stage7_medoids.py --layer design --results results_union
    python scripts/stageU8_build_union_html.py \
        --results results_union --out dist/binder_map_union.html
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.unionhtml import RESOLUTION, bake, join_layers, union_payload

ROOT = Path(__file__).resolve().parent.parent

#: Measured on 2026-09-28 and asserted, never assumed. The whole point of the
#: union tree is that these three populations are in one picture; a silent
#: change to any of them means the picture is of something else.
EXPECTED_NODES = 526_467
EXPECTED_LAYERS = {"in_sample500k": 499_999, "in_passing60k": 31_763, "in_12k": 9_759}
EXPECTED_CONSERVED = {
    "sample_passing": 5_295, "sample_library": 1_553,
    "passing_library": 9_759, "all_three": 1_553,
}


def load(results: Path) -> pd.DataFrame:
    """One row per backbone: node_id, x, y, cluster, struct_name, layer flags."""
    base = results / "design"
    embed = pd.read_parquet(base / "embed" / "fold" / "embed.parquet")
    clusters = pd.read_parquet(base / "clusters" / "fold" / "clusters.parquet")
    nodes = pd.read_csv(base / "dedup" / "nodes.tsv", sep="\t")
    canonical = pd.read_csv(base / "dedup" / "canonical.tsv", sep="\t")
    layers = pd.read_parquet(base / "meta" / "layers.parquet")

    frame = join_layers(nodes, canonical, layers)
    for name, other in (("embed", embed), ("clusters", clusters[["node_id", RESOLUTION]])):
        before = len(frame)
        frame = frame.merge(other, on="node_id", how="left", validate="one_to_one")
        if len(frame) != before:
            raise ValueError(f"{name} join changed the row count: {before} -> {len(frame)}")
    blank = frame[["x", "y", RESOLUTION]].isna().any(axis=1)
    if bool(blank.any()):
        raise ValueError(f"{int(blank.sum())} node(s) have no embedding or cluster")
    return frame.sort_values("node_id", kind="stable").reset_index(drop=True)


def load_medoids(results: Path) -> tuple[pd.DataFrame, dict]:
    """Stage 7's medoid table and CA-trace bundle, or (None, None) if unbuilt.

    Missing files are not fatal: the page degrades to points-only and says so,
    which is a better failure than refusing to build a map that is otherwise
    complete. A HALF-present pair is fatal, because that is a broken stage 7
    run rather than one that never happened.
    """
    base = results / "design"
    table = base / "clusters" / "fold" / "medoids.parquet"
    bundle = base / "medoids" / "medoids.json.gz"
    if not table.exists() and not bundle.exists():
        return None, None
    for path in (table, bundle):
        if not path.exists():
            raise SystemExit(f"FATAL: {path} is missing but its partner is not; "
                             "rerun scripts/stage7_medoids.py")
    return pd.read_parquet(table), json.loads(gzip.open(bundle).read())


def check(frame: pd.DataFrame) -> dict:
    """Assert the measured counts and return the conserved breakdown."""
    if len(frame) != EXPECTED_NODES:
        raise SystemExit(f"FATAL: {len(frame):,d} nodes, expected {EXPECTED_NODES:,d}")
    for column, expected in EXPECTED_LAYERS.items():
        got = int(frame[column].sum())
        if got != expected:
            raise SystemExit(f"FATAL: {column} = {got:,d}, expected {expected:,d}")

    s = frame["in_sample500k"].to_numpy(bool)
    p = frame["in_passing60k"].to_numpy(bool)
    k = frame["in_12k"].to_numpy(bool)
    if not (s | p | k).all():
        raise SystemExit(f"FATAL: {int((~(s | p | k)).sum()):,d} node(s) carry no layer flag")
    conserved = {"sample_passing": int((s & p).sum()),
                 "sample_library": int((s & k).sum()),
                 "passing_library": int((p & k).sum()),
                 "all_three": int((s & p & k).sum())}
    if conserved != EXPECTED_CONSERVED:
        raise SystemExit(f"FATAL: conserved counts {conserved}, expected {EXPECTED_CONSERVED}")
    conserved["any"] = int((s.astype(int) + p + k > 1).sum())
    return conserved


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "results_union")
    parser.add_argument("--out", type=Path, default=ROOT / "dist" / "binder_map_union.html")
    parser.add_argument("--template", type=Path, default=ROOT / "web_union" / "template.html")
    parser.add_argument("--fflate", type=Path,
                        default=ROOT / "web_union" / "vendor" / "fflate.min.js")
    parser.add_argument("--threedmol", type=Path,
                        default=ROOT / "web" / "vendor" / "3Dmol-min.js")
    parser.add_argument("--no-medoids", action="store_true",
                        help="build the points-only page: no CA traces, no viewer")
    parser.add_argument("--max-mb", type=float, default=25.0,
                        help="size target; over it, campaign-ONLY nodes are "
                             "decimated and the page says by how much")
    args = parser.parse_args()

    frame = load(args.results)
    conserved = check(frame)
    medoids, traces = (None, None) if args.no_medoids else load_medoids(args.results)
    print(f"{len(frame)} nodes")
    print(f"  campaign {int(frame['in_sample500k'].sum()):,d} · "
          f"passing {int(frame['in_passing60k'].sum()):,d} · "
          f"library {int(frame['in_12k'].sum()):,d}")
    print(f"  conserved: sample∩passing {conserved['sample_passing']:,d} · "
          f"sample∩library {conserved['sample_library']:,d} · "
          f"passing∩library {conserved['passing_library']:,d} · "
          f"all three {conserved['all_three']:,d} · any {conserved['any']:,d}")
    if traces is None:
        print("  no medoid bundle: building the points-only page")
        threedmol = None
    else:
        print(f"  medoids: {len(traces):,d} CA traces in the stage 7 bundle, "
              f"{int((medoids['resolution'] == RESOLUTION).sum()):,d} at {RESOLUTION}")
        threedmol = args.threedmol

    decimation, payload, size = 1, None, 0
    while True:
        payload = union_payload(frame, conserved, decimation, medoids, traces)
        payload["payload_mb"] = None            # placeholder, filled after baking
        size = bake(payload, args.template, args.fflate, args.out, threedmol)
        if size <= args.max_mb * 1e6 or decimation >= 32:
            break
        # Over target. Thin the campaign-only layer -- never the passing pool,
        # the library, or anything conserved -- and say so on the page.
        decimation = max(2, int(np.ceil(size / (args.max_mb * 1e6))) * decimation)
        print(f"  {size / 1e6:.1f} MB over the {args.max_mb:.0f} MB target -> "
              f"decimating the campaign-only layer {decimation}x")

    payload["payload_mb"] = round(size / 1e6, 1)
    size = bake(payload, args.template, args.fflate, args.out, threedmol)

    if decimation > 1:
        print(f"  campaign-only layer decimated {decimation}x "
              f"({payload['n']:,d} of {len(frame):,d} nodes drawn)")
    if size > args.max_mb * 1e6:
        print(f"WARNING: {size / 1e6:.1f} MB is over the {args.max_mb:.0f} MB target")
    print(f"wrote {args.out} ({size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
