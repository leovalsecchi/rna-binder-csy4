#!/usr/bin/env python3
"""Static preview figures for the binder map. Runs locally after results are pulled.

The HTML explorer is the real deliverable; these are for putting in a deck or a
message without asking anyone to open a 30 MB page.

Usage:
    python scripts/make_preview.py --results local_results --build build \
        --out dist/preview
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

LAYERS = ("design", "af3", "esm")
AXES = (("fold", "Fold — Foldseek backbone similarity"),
        ("epitope", "Epitope — RNA contact fingerprint"))

# Same palette as the explorer, so figures and app agree.
PALETTE = {
    "RFD3+LiMPNN": "#4a7ebb", "BG+LiMPNN": "#e07b39", "BG+BIF": "#4f9d69",
    "BD1+LiMPNN": "#b5514f", "HD+LiMPNN": "#8264a5",
}


def load_layer(results: Path, build: Path, layer: str):
    base = results / layer
    if not (base / "dedup" / "nodes.tsv").exists():
        return None
    nodes = pd.read_csv(base / "dedup" / "nodes.tsv", sep="\t")
    meta = pd.read_parquet(build / "meta" / "designs.parquet").set_index("design_id")
    stack = meta.reindex(nodes["struct_name"])["stack"].to_numpy()

    out = {"nodes": nodes, "stack": stack, "axes": {}}
    for axis, _ in AXES:
        embed = base / "embed" / axis / "embed.parquet"
        clust = base / "clusters" / axis / "clusters.parquet"
        if embed.exists() and clust.exists():
            out["axes"][axis] = {
                "xy": pd.read_parquet(embed).sort_values("node_id"),
                "cl": pd.read_parquet(clust).sort_values("node_id"),
            }
    return out


def scatter(ax, xy, stack, title, point_size):
    for name, colour in PALETTE.items():
        m = stack == name
        if not m.any():
            continue
        ax.scatter(xy["x"].to_numpy()[m], xy["y"].to_numpy()[m],
                   s=point_size, c=colour, label=name, alpha=0.7, linewidths=0)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor("#dddddd")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("local_results"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--out", type=Path, default=Path("dist/preview"))
    ap.add_argument("--resolution", default="cluster_r10")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    loaded = {l: d for l in LAYERS if (d := load_layer(args.results, args.build, l))}
    if not loaded:
        print("no layers found under", args.results)
        return 1

    # One row per layer, one column per axis.
    n_rows = len(loaded)
    fig, axes = plt.subplots(n_rows, 2, figsize=(13, 5.6 * n_rows), squeeze=False)
    for r, (layer, d) in enumerate(loaded.items()):
        n = len(d["nodes"])
        size = 1.2 if n > 20000 else (4 if n > 5000 else 14)
        for c, (axis, label) in enumerate(AXES):
            ax = axes[r][c]
            if axis not in d["axes"]:
                ax.set_visible(False)
                continue
            scatter(ax, d["axes"][axis]["xy"], d["stack"],
                    f"{layer} — {label}  (n={n:,})", size)
        if r == 0:
            axes[0][0].legend(loc="best", fontsize=8, frameon=False, markerscale=4)
    fig.suptitle("Binder map — fold vs epitope, coloured by generation stack",
                 fontsize=13)
    fig.tight_layout()
    dest = args.out / "maps_by_stack.png"
    fig.savefig(dest, dpi=150, facecolor="white")
    print("wrote", dest)

    # Fold vs epitope agreement, per layer.
    rows = []
    for layer, d in loaded.items():
        if len(d["axes"]) < 2:
            continue
        a = d["axes"]["fold"]["cl"][args.resolution].to_numpy()
        b = d["axes"]["epitope"]["cl"][args.resolution].to_numpy()
        rows.append({
            "layer": layer,
            "n": len(a),
            "fold_clusters": int(len(np.unique(a))),
            "epitope_clusters": int(len(np.unique(b))),
            "AMI": adjusted_mutual_info_score(a, b),
            "ARI": adjusted_rand_score(a, b),
            "AMI_fold_stack": adjusted_mutual_info_score(a, d["stack"]),
            "AMI_epitope_stack": adjusted_mutual_info_score(b, d["stack"]),
        })
    if rows:
        table = pd.DataFrame(rows)
        table.to_csv(args.out / "axis_agreement.csv", index=False)
        print("\n=== fold vs epitope, at", args.resolution, "===")
        print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print("\nA low AMI means backbone fold does not predict where a binder")
        print("lands on the RNA. That is a result, not a failure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
