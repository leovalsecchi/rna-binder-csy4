#!/usr/bin/env python3
"""Stage 9: compare the three fold maps.

Three clusterings built under identical conditions, differing only in sample
size: the 500k campaign sample, the 31,763-node filtered pool, and the
9,759-node ordered library.

Nothing here recomputes a gate. `in_passing` and `in_12k` were set in stage 0 by
a join against exports that already exist and are already audited.

Cluster IDS ARE NOT COMPARABLE ACROSS MAPS -- three independent Leiden runs on
three node sets. The crosswalk reports ARI/AMI on designs present in more than
one map as a RESULT, not as a correspondence to be relied on.

Usage:
    python scripts/stage9_compare.py --out compare
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")

#: The Leiden resolution the 12k selection consumed (v6hier_q.py:61).
RESOLUTION = "cluster_r80"

#: name -> (results dir, build dir or None for the pre-existing reference map)
MAPS = {
    "sample500k": (FM / "results_full" / "design", FM / "build_full"),
    "passing60k": (FM / "results_sw_k30_f50_mut" / "design", None),
    "library12k": (FM / "results_lib12k" / "design", FM / "build_lib12k"),
}


def expand_members(members: pd.DataFrame) -> pd.DataFrame:
    """One row per design: (struct_name, design_id).

    A graph node is a coordinate hash standing for one or more design_ids, so
    per-design labels must be expanded before they can be joined.
    """
    out = members.assign(design_id=members["design_ids"].str.split(",")).explode("design_id")
    out["design_id"] = out["design_id"].str.strip()
    return out[["struct_name", "design_id"]].reset_index(drop=True)


def design_to_cluster(results: Path, resolution: str = RESOLUTION) -> pd.Series:
    """design_id -> fold cluster id, for one completed map."""
    clusters = pd.read_parquet(results / "clusters" / "fold" / "clusters.parquet")
    nodes = pd.read_csv(results / "dedup" / "nodes.tsv", sep="\t")
    members = expand_members(pd.read_csv(results / "dedup" / "members.tsv", sep="\t"))
    joined = (members.merge(nodes[["struct_name", "node_id"]], on="struct_name")
                     .merge(clusters[["node_id", resolution]], on="node_id"))
    return joined.set_index("design_id")[resolution]


#: Collapse the dedup key's known over-splitting: BG+BIF seq variants are the same
#: fold rigid-body displaced (0.94 A superposed RMSD), and HD sample variants
#: differ by 0.022 A of packer jitter. Both measured 2026-09-14.
_SEQ = re.compile(r"_seq\d+$")
_SAMPLE_PACK = re.compile(r"_sample\d+_pack\d+$")


def fold_key(design_id: str) -> str:
    """The distinct FOLD a node belongs to, as opposed to its coord-hash node."""
    if "_HD_LiMPNN_" in design_id:
        return _SAMPLE_PACK.sub("", design_id)
    return _SEQ.sub("", design_id)


def summarise(labels: pd.DataFrame, name: str) -> dict:
    sizes = labels.groupby("cluster").size()
    return {
        "map": name,
        "nodes": int(len(labels)),
        "designs": int(labels["n_designs"].sum()),
        "clusters": int(labels["cluster"].nunique()),
        "singletons": int((sizes == 1).sum()),
        "largest_cluster": int(sizes.max()),
        "median_cluster": float(sizes.median()),
        "p90_cluster": float(sizes.quantile(0.90)),
        "n_passing": int(labels["in_passing"].sum()),
        "n_12k": int(labels["in_12k"].sum()),
    }


def cluster_table(labels: pd.DataFrame) -> pd.DataFrame:
    """One row per fold cluster, with node count AND distinct-fold count.

    Both size columns are mandatory and must be computed over the same rows.
    The dedup key over-splits BG+BIF ~9.8x and HD ~10x (rigid-body offset and
    0.022 A packer jitter respectively), so `size` alone is an inflated
    denominator for any cluster-size or fold-coverage claim; `folds` is what a
    fold-level reading needs. Aggregating `struct_name` here instead of
    `fold_key` would make the two columns identical and silently hide exactly
    the inflation this table exists to expose.
    """
    return (
        labels.groupby("cluster")
        .agg(
            size=("struct_name", "size"),
            folds=("fold_key", "nunique"),
            n_designs=("n_designs", "sum"),
            n_passing=("in_passing", "sum"),
            n_12k=("in_12k", "sum"),
            top_stack=("stack", lambda x: x.value_counts().idxmax()),
        )
        .reset_index()
        .sort_values("size", ascending=False)
    )


def crosswalk_scores(old: pd.Series, new: pd.Series) -> dict:
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

    shared = old.index.intersection(new.index)
    if len(shared) == 0:
        return {"n_shared": 0, "ari": float("nan"), "ami": float("nan")}
    return {
        "n_shared": int(len(shared)),
        "ari": float(adjusted_rand_score(old.loc[shared], new.loc[shared])),
        "ami": float(adjusted_mutual_info_score(old.loc[shared], new.loc[shared])),
    }


def load_labels(results: Path, build: Path | None) -> pd.DataFrame:
    clusters = pd.read_parquet(results / "clusters" / "fold" / "clusters.parquet")
    nodes = pd.read_csv(results / "dedup" / "nodes.tsv", sep="\t")
    # nodes.tsv already carries n_designs and per_node recomputes it from the
    # stage-0 metadata below; merging both emits n_designs_x/_y and the fillna
    # further down dies with KeyError. The aggregated one is authoritative.
    nodes = nodes.drop(columns=["n_designs"], errors="ignore")
    members = expand_members(pd.read_csv(results / "dedup" / "members.tsv", sep="\t"))

    if build is not None:
        meta = pd.read_parquet(build / "meta" / "nodes.parquet")[
            ["design_id", "stack", "in_passing", "in_12k", "n_designs"]]
    else:
        # The reference map predates stage 0, so rebuild its labels the same way:
        # every design in it is by definition in `passing`, and 12k membership
        # comes from the library export.
        lib = set(pd.read_csv(
            FM.parent / "library_selection" / "csy4_binder_library_12k_v6_k30.csv",
            usecols=["design_id"])["design_id"])
        meta = pd.DataFrame({"design_id": members["design_id"].unique()})
        meta["stack"] = "unknown"
        meta["in_passing"] = True
        meta["in_12k"] = meta["design_id"].isin(lib)
        meta["n_designs"] = 1

    joined = members.merge(meta, on="design_id", how="left")
    per_node = joined.groupby("struct_name").agg(
        stack=("stack", lambda x: x.dropna().iloc[0] if x.notna().any() else "unknown"),
        in_passing=("in_passing", "any"),
        in_12k=("in_12k", "any"),
        n_designs=("n_designs", "sum"),
    ).reset_index()

    labels = (nodes.merge(per_node, on="struct_name", how="left")
                   .merge(clusters[["node_id", RESOLUTION]], on="node_id", how="left")
                   .rename(columns={RESOLUTION: "cluster"}))
    labels[["in_passing", "in_12k"]] = labels[["in_passing", "in_12k"]].fillna(False)
    labels["n_designs"] = labels["n_designs"].fillna(0)
    labels["fold_key"] = [fold_key(d) for d in labels["struct_name"]]
    return labels.dropna(subset=["cluster"]).astype({"cluster": int})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=FM / "compare")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    summaries, per_design = [], {}
    for name, (results, build) in MAPS.items():
        if not (results / "clusters" / "fold" / "clusters.parquet").exists():
            print(f"SKIP {name}: not built yet")
            continue
        labels = load_labels(results, build)
        summaries.append(summarise(labels, name))

        table = cluster_table(labels)
        table.to_parquet(args.out / f"cluster_table_{name}.parquet", index=False)
        per_design[name] = design_to_cluster(results)

    summary = pd.DataFrame(summaries)
    summary.to_parquet(args.out / "map_summary.parquet", index=False)
    print(summary.to_string(index=False))

    rows = []
    names = list(per_design)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            scores = crosswalk_scores(per_design[a], per_design[b])
            scores.update({"map_a": a, "map_b": b})
            rows.append(scores)
            print(f"{a} vs {b}: n={scores['n_shared']:,d} "
                  f"ARI={scores['ari']:.3f} AMI={scores['ami']:.3f}")
    pd.DataFrame(rows).to_parquet(args.out / "crosswalk.parquet", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
