#!/usr/bin/env python3
"""Clustering overview: how many fold clusters, how big, and what is in them.

This is the direct answer to "cluster my binders by backbone similarity and
show me by stack" - the cluster table and composition, not just a scatter.

Usage:
    python scripts/cluster_overview.py --layer design --resolution cluster_r10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = {
    "RFD3+LiMPNN": "#4a7ebb", "BG+LiMPNN": "#e07b39", "BG+BIF": "#4f9d69",
    "BD1+LiMPNN": "#b5514f", "HD+LiMPNN": "#8264a5",
}
STACKS = list(PALETTE)


def load(results: Path, build: Path, layer: str, axis: str, resolution: str):
    base = results / layer
    nodes = pd.read_csv(base / "dedup" / "nodes.tsv", sep="\t")
    members = pd.read_csv(base / "dedup" / "members.tsv", sep="\t")
    meta = pd.read_parquet(build / "meta" / "designs.parquet").set_index("design_id")
    clusters = pd.read_parquet(
        base / "clusters" / axis / "clusters.parquet").sort_values("node_id")
    xy = pd.read_parquet(base / "embed" / axis / "embed.parquet").sort_values("node_id")

    df = nodes.merge(members, on="struct_name", how="left")
    df["cluster"] = clusters[resolution].to_numpy()
    df["x"] = xy["x"].to_numpy()
    df["y"] = xy["y"].to_numpy()
    df["stack"] = meta.reindex(df["struct_name"])["stack"].to_numpy()
    # n_designs: how many original designs collapse onto this backbone
    df["n_designs"] = df["n_designs"].fillna(1).astype(int)
    return df


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """One row per cluster: size, dominant stack, and per-stack counts."""
    rows = []
    for cid, g in df.groupby("cluster", sort=True):
        counts = g["stack"].value_counts()
        top = counts.index[0]
        rows.append({
            "cluster": int(cid),
            "backbones": len(g),
            "designs": int(g["n_designs"].sum()),
            "n_stacks": int((counts > 0).sum()),
            "dominant_stack": top,
            "dominant_pct": round(100 * counts.iloc[0] / len(g), 1),
            **{s: int(counts.get(s, 0)) for s in STACKS},
        })
    out = pd.DataFrame(rows).sort_values("backbones", ascending=False)
    return out.reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("local_results"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--axis", default="fold")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--out", type=Path, default=Path("dist/overview"))
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    df = load(args.results, args.build, args.layer, args.axis, args.resolution)
    summary = summarise(df)
    tag = f"{args.layer}_{args.axis}_{args.resolution}"
    summary.to_csv(args.out / f"clusters_{tag}.csv", index=False)

    sizes = summary["backbones"].to_numpy()
    singletons = int((sizes == 1).sum())
    print(f"\n=== {args.layer} / {args.axis} / {args.resolution} ===")
    print(f"  {len(df):,} backbones ({df['n_designs'].sum():,} designs) "
          f"in {len(summary):,} clusters")
    print(f"  largest {sizes.max():,} | median {int(np.median(sizes))} | "
          f"singletons {singletons:,}")
    covered = sizes[sizes >= 10].sum()
    print(f"  clusters with >=10 backbones: {(sizes >= 10).sum():,} "
          f"covering {covered:,} backbones ({100*covered/len(df):.1f}%)")
    pure = summary[summary["dominant_pct"] >= 90]
    print(f"  single-stack clusters (>=90% one stack): {len(pure):,} "
          f"of {len(summary):,}")

    print(f"\n  top {args.top} clusters:")
    show = summary.head(args.top)[
        ["cluster", "backbones", "designs", "dominant_stack", "dominant_pct"] + STACKS]
    print(show.to_string(index=False))

    # ---- figure: map by cluster, size distribution, stack composition
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.28, wspace=0.22)

    ax0 = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(0)
    order = summary["cluster"].to_numpy()
    rank = {c: i for i, c in enumerate(order)}
    cmap = plt.get_cmap("tab20")
    colours = np.array([cmap(rank[c] % 20) if rank[c] < 60 else (0.85, 0.85, 0.85, 1.0)
                        for c in df["cluster"]])
    ax0.scatter(df["x"], df["y"], s=1.2, c=colours, linewidths=0)
    ax0.set_title(f"{args.layer} fold map — coloured by cluster "
                  f"({len(summary):,} clusters, top 60 coloured)", fontsize=10)
    ax0.set_xticks([]); ax0.set_yticks([])

    ax1 = fig.add_subplot(gs[0, 1])
    for s, c in PALETTE.items():
        m = df["stack"].to_numpy() == s
        ax1.scatter(df["x"].to_numpy()[m], df["y"].to_numpy()[m], s=1.2, c=c,
                    label=s, linewidths=0, alpha=0.7)
    ax1.legend(fontsize=8, frameon=False, markerscale=6, loc="best")
    ax1.set_title(f"{args.layer} fold map — coloured by stack", fontsize=10)
    ax1.set_xticks([]); ax1.set_yticks([])

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.hist(sizes, bins=np.logspace(0, np.log10(max(sizes.max(), 2)), 40),
             color="#4a7ebb", edgecolor="white", linewidth=0.4)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("cluster size (backbones)", fontsize=9)
    ax2.set_ylabel("clusters", fontsize=9)
    ax2.set_title("cluster size distribution", fontsize=10)

    ax3 = fig.add_subplot(gs[1, 1])
    top = summary.head(20)
    bottom = np.zeros(len(top))
    idx = np.arange(len(top))
    for s in STACKS:
        vals = top[s].to_numpy()
        ax3.bar(idx, vals, bottom=bottom, color=PALETTE[s], label=s, width=0.85)
        bottom += vals
    ax3.set_xticks(idx)
    ax3.set_xticklabels(top["cluster"].astype(str), fontsize=6, rotation=90)
    ax3.set_xlabel("cluster id (20 largest)", fontsize=9)
    ax3.set_ylabel("backbones", fontsize=9)
    ax3.set_title("stack composition of the largest clusters", fontsize=10)

    fig.suptitle(
        f"Fold clustering overview — {args.layer} layer, {len(df):,} backbones, "
        f"resolution {args.resolution}", fontsize=13)
    dest = args.out / f"overview_{tag}.png"
    fig.savefig(dest, dpi=150, facecolor="white", bbox_inches="tight")
    print(f"\n  wrote {dest}")
    print(f"  wrote {args.out / f'clusters_{tag}.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
