#!/usr/bin/env python3
"""clean.py — filter the raw collection and draw a per-stack iteration sample.

Step 2 of the analyses/ pipeline (collect -> clean -> compute_metrics).

Takes collect.py's raw CSV and applies, in order:
    1. drop designs whose binder sequence contains cysteine (C)
    2. drop designs missing any of the AF3-msaoff / B2 / ESM refolds
    3. drop exact sequence duplicates, keeping the highest AF3-msaoff iPTM

Order note: the spec listed dedup before the all-three requirement, but that can
silently discard a *complete* design (AF3+B2+ESM) in favour of a sequence-duplicate
that scores higher on AF3 yet lacks B2/ESM — which then gets dropped, losing both.
Requiring all three refolds *before* dedup means dedup only ever chooses among
fully-refolded designs. (Confirmed with the user.)

Writes:
    collected_clean.csv   all passing, deduplicated designs
    sample10k.csv          a random, seed-fixed sample of N per stack (default 600;
                          5 stacks -> 3000) drawn from the passing designs, for fast
                          iteration in compute_metrics.py.

Rows with no sequence at all can't be cys-checked or deduped, so they're dropped in
stage 1 and reported separately.

Usage
-----
    python3 clean.py --input collected_raw.csv --outdir analyses/<name> \\
        [--per-stack 600] [--seed 0]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def _funnel(label: str, before: int, after: int):
    print(f"  {label:<34} {before:>7} -> {after:>7}  ({before - after:+d})")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", "-i", type=Path, required=True,
                    help="raw CSV from collect.py")
    ap.add_argument("--outdir", "-d", type=Path, default=None,
                    help="output dir (default: input's dir). Writes "
                         "collected_clean.csv + sample10k.csv here.")
    ap.add_argument("--per-stack", type=int, default=2000,
                    help="designs sampled per stack for sample10k.csv (default 600)")
    ap.add_argument("--seed", type=int, default=0,
                    help="RNG seed for the per-stack sample (default 0)")
    args = ap.parse_args()

    outdir = args.outdir or args.input.parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, low_memory=False)
    n0 = len(df)
    print(f"Loaded {n0} designs from {args.input}")
    print("Funnel:")

    # --- 1. cysteine (drops null-sequence rows too: they can't be cys-checked) ---
    has_seq = df["sequence"].notna()
    n_no_seq = int((~has_seq).sum())
    no_cys = has_seq & ~df["sequence"].str.contains("C", case=True, na=True)
    before = len(df)
    df = df[no_cys].copy()
    _funnel("1. remove cysteine / no-seq", before, len(df))
    if n_no_seq:
        print(f"     ({n_no_seq} of those had no sequence and were dropped)")

    # --- 2. require AF3-msaoff + B2 + ESM ---
    for c in ("has_af3_msaoff", "has_b2", "has_esm"):
        df[c] = df[c].astype(str).str.lower().isin(("true", "1"))
    before = len(df)
    df = df[df["has_af3_msaoff"] & df["has_b2"] & df["has_esm"]].copy()
    _funnel("2. require AF3+B2+ESM refolds", before, len(df))

    # --- 3. dedup exact sequences, keep highest AF3-msaoff iPTM ---
    before = len(df)
    df = (df.sort_values("af3_msaoff_iptm", ascending=False, na_position="last")
            .drop_duplicates(subset="sequence", keep="first"))
    _funnel("3. dedup sequences (keep best AF3)", before, len(df))

    df = df.sort_values(["stack", "uuid8", "design_id"]).reset_index(drop=True)
    clean_path = outdir / "collected_clean.csv"
    df.to_csv(clean_path, index=False)
    print(f"\nWrote {len(df)} clean designs -> {clean_path}")
    print(f"  per stack: {df.groupby('stack').size().to_dict()}")

    # --- per-stack random sample (seed-fixed) ---
    # Plain stack loop (no groupby.apply / include_groups, which changed across
    # pandas 2.x->3.x) so the 'stack' column is always retained.
    parts = [g.sample(n=min(args.per_stack, len(g)), random_state=args.seed)
             for _, g in df.groupby("stack", sort=True)]
    sample = pd.concat(parts) if parts else df.iloc[0:0]
    sample = sample.sort_values(["stack", "design_id"]).reset_index(drop=True)
    sample_path = outdir / "sample10k.csv"
    sample.to_csv(sample_path, index=False)
    print(f"Wrote {len(sample)} sampled designs "
          f"(<= {args.per_stack}/stack, seed={args.seed}) -> {sample_path}")
    print(f"  per stack: {sample.groupby('stack').size().to_dict()}")


if __name__ == "__main__":
    main()
