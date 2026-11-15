#!/usr/bin/env python3
"""Bake the binder_map explorer against the FINE Foldseek run, with the 12k library
selections layered on.

The stock explorer (stage8_build_html.py) bakes the coarse run and knows nothing about
the library. This one:

  * stages a results tree that reuses the coarse run's dedup / traces / complexes
    (those depend only on the design SET, which is identical between runs) but swaps in
    the fine run's fold clusters (r05..r80) and its UMAP,
  * injects a `library` categorical per node -- base only / cluster only / both /
    in pool / not in pool -- so it shows up in the colour dropdown,
  * injects numeric `lib_base_n`, `lib_cluster_n`, `pool_n` so you can also colour by
    how many of a backbone's designs made each library.

Everything else -- the 3Dmol backbone+RNA viewer, the CA traces, cluster tabs, and the
"open in ChimeraX" bridge calls -- is the existing app.js untouched.

    python3 scripts/build_library_explorer.py --embed <fine embed.parquet>
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from binder_map.htmlbuild import bake, build_payload, encode_categorical, encode_floats

EXPORTS = Path("/Users/leo/Documents/binder_filter/exports")
LIBS = {"base": EXPORTS / "library_12k_base_iptm.csv",
        "cluster": EXPORTS / "library_12k_cluster_rr_r80.csv"}


def stage(coarse: Path, fine_clusters: Path, fine_embed: Path, dest: Path,
          assets: Path | None = None) -> Path:
    """Fold axis only: mixing a coarse epitope axis with fine fold clusters would offer
    resolutions on one axis that do not exist on the other."""
    d = dest / "design"
    (d / "clusters" / "fold").mkdir(parents=True, exist_ok=True)
    (d / "embed" / "fold").mkdir(parents=True, exist_ok=True)
    # complex/ carries the backbone+RNA cartoon (the good viewer); traces/ is the CA
    # fallback; medoids/ backs the cluster tab. Symlink rather than copy: these are
    # tens of MB and live in iCloud, where a copy stalls on dataless files.
    for sub in ("dedup", "traces", "fp", "complex", "medoids"):
        src = (assets / "design" / sub) if (assets and (assets / "design" / sub).exists())               else (coarse / "design" / sub).resolve()
        if src.exists():
            (d / sub).symlink_to(src, target_is_directory=True)
    src_graph = (coarse / "design" / "graph" / "fold").resolve()
    if src_graph.exists():
        (d / "graph").mkdir(exist_ok=True)
        (d / "graph" / "fold").symlink_to(src_graph, target_is_directory=True)
    shutil.copy(fine_clusters, d / "clusters" / "fold" / "clusters.parquet")
    shutil.copy(fine_embed, d / "embed" / "fold" / "embed.parquet")
    return dest


def library_columns(results: Path, build: Path) -> tuple[list[str], np.ndarray, dict]:
    """Per-NODE library membership. A node is a backbone; the libraries pick designs, so
    a backbone counts as 'in' when any of its designs was selected."""
    nodes = pd.read_csv(results / "design" / "dedup" / "nodes.tsv", sep="\t")
    members = pd.read_csv(results / "design" / "dedup" / "members.tsv", sep="\t")
    d2s: dict[str, str] = {}
    for sn, ids in zip(members["struct_name"], members["design_ids"]):
        for i in str(ids).split(","):
            d2s[i] = sn
    idx = pd.Series(np.arange(len(nodes)), index=nodes["struct_name"])

    counts = {}
    for tag, path in LIBS.items():
        c = np.zeros(len(nodes), dtype=np.float32)
        if path.exists():
            sel = pd.read_csv(path, usecols=["design_id"])["design_id"]
            hit = idx.reindex(pd.Index([d2s.get(x) for x in sel])).dropna()
            for n, k in hit.value_counts().items():
                c[int(n)] = k
        else:
            print(f"  WARNING: {path.name} missing -> {tag} column is all zero")
        counts[tag] = c

    pool = np.zeros(len(nodes), dtype=np.float32)
    exp = EXPORTS / "passing_v1.csv"
    if exp.exists():
        ids = pd.read_csv(exp, skiprows=2, usecols=["design_id"])["design_id"]
        hit = idx.reindex(pd.Index([d2s.get(x) for x in ids])).dropna()
        for n, k in hit.value_counts().items():
            pool[int(n)] = k

    a, b = counts["base"] > 0, counts["cluster"] > 0
    lab = np.where(a & b, "both",
          np.where(b, "cluster only",
          np.where(a, "base only",
          np.where(pool > 0, "in pool, unpicked", "not in pool"))))
    return lab.tolist(), pool, counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", type=Path, default=ROOT / "local_results")
    ap.add_argument("--fine-clusters", type=Path,
                    default=ROOT / "local_results_fine/design/clusters/fold/clusters.parquet")
    ap.add_argument("--embed", type=Path, required=True,
                    help="fine run design/embed/fold/embed.parquet")
    ap.add_argument("--assets", type=Path, default=None,
                    help="dir with design/{complex,traces,medoids} pulled off Euler; "
                         "iCloud evicts these locally and a read just hangs")
    ap.add_argument("--build", type=Path, default=ROOT / "build")
    ap.add_argument("--out", type=Path, default=ROOT / "dist" / "library_explorer.html")
    args = ap.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="libexp_"))
    try:
        results = stage(args.coarse, args.fine_clusters, args.embed, tmp / "results",
                        args.assets)
        payload = build_payload(results, args.build)
        entry = payload["layers"]["design"]

        lab, pool, counts = library_columns(results, args.build)
        lookup, codes = encode_categorical(lab)
        entry.setdefault("categorical", {})["library"] = {"lookup": lookup, "codes": codes}
        entry.setdefault("numeric", {})["lib_base_n"] = encode_floats(counts["base"])
        entry["numeric"]["lib_cluster_n"] = encode_floats(counts["cluster"])
        entry["numeric"]["pool_n"] = encode_floats(pool)

        n = len(lab)
        for k in lookup:
            print(f"  {k:22s} {lab.count(k):>7,} backbones "
                  f"({100 * lab.count(k) / n:5.1f}%)")
        axes = ", ".join(entry["axes"])
        res = ", ".join(sorted(entry["axes"]["fold"]["clusters"]))
        print(f"  design: {entry['n']:,} nodes  axes: {axes}  resolutions: {res}")

        size = bake(payload, ROOT / "web" / "template.html", ROOT / "web" / "app.js",
                    ROOT / "web" / "style.css", args.out)
        print(f"wrote {args.out} ({size / 1e6:.1f} MB)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
