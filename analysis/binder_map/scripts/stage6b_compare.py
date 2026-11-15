#!/usr/bin/env python3
"""Stage 6b: compare the fold and epitope partitions, and epitopes across layers.

This is the payoff of having two axes. A low AMI is a real result — backbone
fold does not predict where a binder lands — not a failure.

Usage:  python scripts/stage6b_compare.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.clustering import RESOLUTIONS
from binder_map.compare import (
    contingency, fingerprint_agreement, partition_scores, stack_profiles,
)
from binder_map import RNA_LEN

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
LAYERS = ("design", "af3", "esm")
PER_NT = [f"per_nt_{i}" for i in range(1, RNA_LEN + 1)]


def load_axis_clusters(results: Path, layer: str, axis: str):
    return pd.read_parquet(results / layer / "clusters" / axis / "clusters.parquet")


def design_fingerprints(results: Path, layer: str) -> pd.DataFrame:
    """Fingerprints keyed by design_id, fanned out from canonical structures.

    The design layer collapses sequence variants onto one backbone, so one
    canonical row stands for many design_ids.
    """
    canonical = pd.read_parquet(results / layer / "fp" / "canonical.parquet")
    members = pd.read_csv(results / layer / "dedup" / "members.tsv", sep="\t")
    canonical = canonical.merge(
        members, left_on="design_id", right_on="struct_name", how="left")
    canonical["design_ids"] = canonical["design_ids"].fillna(canonical["design_id"])
    expanded = canonical.assign(
        design_id=canonical["design_ids"].str.split(",")).explode("design_id")
    return expanded.drop(columns=["design_ids", "struct_name"]).set_index("design_id")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=FM / "results")
    parser.add_argument("--build", type=Path, default=FM / "build")
    args = parser.parse_args()

    out = args.results / "compare"
    out.mkdir(parents=True, exist_ok=True)

    # ---- 1. fold x epitope, per layer and resolution
    scores = []
    for layer in LAYERS:
        try:
            fold = load_axis_clusters(args.results, layer, "fold")
            epi = load_axis_clusters(args.results, layer, "epitope")
        except FileNotFoundError:
            print(f"  {layer}: clusters missing, skipping")
            continue
        merged = fold.merge(epi, on="node_id", suffixes=("_fold", "_epi"))
        for resolution in RESOLUTIONS:
            a = merged[f"{resolution}_fold"].to_numpy()
            b = merged[f"{resolution}_epi"].to_numpy()
            s = partition_scores(a, b)
            s.update({"layer": layer, "resolution": resolution})
            scores.append(s)
            print(f"  {layer}/{resolution}: AMI={s['ami']:.3f} ARI={s['ari']:.3f} "
                  f"n={s['n']:,d}")

            mat, la, lb = contingency(a, b)
            pd.DataFrame(mat, index=la, columns=lb).to_parquet(
                out / f"contingency_{layer}_{resolution}.parquet")
    if scores:
        pd.DataFrame(scores).to_parquet(out / "partition_scores.parquet", index=False)

    # ---- 2/3. cross-layer epitope agreement (and design-vs-refold)
    fps = {}
    for layer in LAYERS:
        try:
            fps[layer] = design_fingerprints(args.results, layer)
        except FileNotFoundError:
            print(f"  {layer}: fingerprints missing, skipping agreement")

    from binder_map.fingerprint import FINGERPRINT_COLUMNS
    agreement = None
    pairs = [("af3", "esm"), ("design", "af3"), ("design", "esm")]
    for a, b in pairs:
        if a not in fps or b not in fps:
            continue
        common = fps[a].index.intersection(fps[b].index)
        va = fps[a].loc[common, FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
        vb = fps[b].loc[common, FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
        col = pd.DataFrame({f"epi_agree_{a}_{b}": fingerprint_agreement(va, vb)},
                           index=common)
        agreement = col if agreement is None else agreement.join(col, how="outer")
        print(f"  epitope agreement {a} vs {b}: n={len(common):,d}, "
              f"median={np.median(col.iloc[:, 0]):.3f}")
    if agreement is not None:
        agreement.reset_index(names="design_id").to_parquet(
            out / "epitope_agreement.parquet", index=False)

    # ---- 4. per-stack nucleotide profiles
    meta = pd.read_parquet(args.build / "meta" / "designs.parquet")
    profiles = []
    for layer, fp in fps.items():
        joined = fp.join(meta.set_index("design_id")["stack"], how="inner")
        frame = stack_profiles(joined[PER_NT].to_numpy(dtype=np.float64),
                               joined["stack"].to_numpy())
        frame.insert(0, "layer", layer)
        profiles.append(frame.reset_index())
    if profiles:
        pd.concat(profiles, ignore_index=True).to_parquet(
            out / "stack_profiles.parquet", index=False)

    print(f"comparisons written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
