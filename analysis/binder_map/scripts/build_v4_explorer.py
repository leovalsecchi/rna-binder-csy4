#!/usr/bin/env python3
"""Bake the binder_map explorer against the passing_v4 fold map (`results_v4`).

Differences from build_library_explorer.py:
  * no staging/symlinking -- results_v4 has its own dedup / traces / complex /
    medoids, because v4 is a DIFFERENT design set from the run that built `results`,
  * no `library` / `lib_*` columns -- those 12k selections were computed on
    passing_v1 and would be wrong here,
  * the coarse ~20-cluster partition goes in as its OWN categorical `fold20`,
    not as a sixth cluster_r* resolution: it comes from a different graph
    (TM floor 0.30, non-mutual), so presenting it as another resolution of the
    same clustering would misrepresent it.

    python3 scripts/build_v4_explorer.py
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from binder_map.htmlbuild import bake, build_payload, encode_categorical, encode_floats

STAGE = Path.home() / "binder_map_v4"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=STAGE / "results")
    ap.add_argument("--build", type=Path, default=STAGE / "build")
    ap.add_argument("--out", type=Path, default=STAGE / "binder_map_v4.html")
    args = ap.parse_args()

    payload = build_payload(args.results, args.build)
    entry = payload["layers"]["design"]
    nodes = pd.read_csv(args.results / "design" / "dedup" / "nodes.tsv", sep="\t")

    cpath = args.results / "design" / "clusters" / "fold_coarse" / "clusters.parquet"
    if cpath.exists():
        coarse = pd.read_parquet(cpath).set_index("node_id")
        lab = coarse["cluster_coarse"].reindex(nodes["node_id"]).to_numpy()
        iso = coarse["is_isolate"].reindex(nodes["node_id"]).to_numpy()
        # Rank clusters by size so c00 is the biggest -- a stable, readable order.
        order = pd.Series(lab).value_counts().index.tolist()
        rank = {c: i for i, c in enumerate(order)}
        names = ["isolate" if b else f"c{rank[c]:02d}" for c, b in zip(lab, iso)]
        lookup, codes = encode_categorical(names)
        entry.setdefault("categorical", {})["fold20"] = {"lookup": lookup, "codes": codes}
        size = pd.Series(lab).value_counts()
        entry.setdefault("numeric", {})["fold20_size"] = encode_floats(
            size.reindex(lab).to_numpy().astype(float))
        n_iso = int(np.asarray(iso).sum())
        print(f"  fold20: {len(set(names)) - (1 if n_iso else 0)} clusters "
              f"+ {n_iso} isolates over {len(nodes):,} backbones")

    res = ", ".join(sorted(entry["axes"]["fold"]["clusters"]))
    print(f"  design: {entry['n']:,} nodes   axes: {', '.join(entry['axes'])}   "
          f"resolutions: {res}")
    print(f"  colour-by: {', '.join(sorted(entry['categorical']))}")

    size = bake(payload, ROOT / "web" / "template.html", ROOT / "web" / "app.js",
                ROOT / "web" / "style.css", args.out)
    print(f"wrote {args.out}  ({size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
