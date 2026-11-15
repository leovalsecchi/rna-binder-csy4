#!/usr/bin/env python3
"""Plot per-nucleotide protein->RNA contact patterns from the aggregated structural CSV.

Reads a structural-aggregated CSV that carries rna_contact_nt_1..N columns (added by
analysis/structural.py: per-nucleotide heavy-atom contact counts, 4.5 A). Produces:

  (1) rna_contact_heatmap.png   - per-stack heatmap: designs (rows, sorted by AF3 msaoff
                                  iPTM) x nucleotide position (cols), color = contacts.
  (2) rna_contact_mean_bar.png  - mean contacts per nucleotide position, grouped by stack.

NOTE on sorting: the spec asks to sort by AF3 msaoff iPTM. If the CSV has no AF3 iPTM
column (e.g. a B2-only aggregation while AF3 is un-canonicalized), it falls back to
B2 interface iPTM (b2_protein_rna_iptm) and prints a warning. Override with --sort-col.

Usage:
    plot_rna_contact_heatmap.py STRUCTURAL_CSV OUTDIR [--sort-col COL] [--refolder b2|af3]
"""
import argparse
import os
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", f"/tmp/mplconfig_{os.getpid()}")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def pick_sort_col(df, requested):
    """Return (column, fell_back_bool). Prefer AF3 msaoff iPTM, else B2 iPTM."""
    if requested:
        if requested in df.columns:
            return requested, False
        sys.exit(f"--sort-col {requested!r} not in CSV")
    if "af3_protein_rna_iptm" in df.columns and df["af3_protein_rna_iptm"].notna().any():
        return "af3_protein_rna_iptm", False        # AF3 interface iPTM (msaoff handled by --refolder)
    if "b2_protein_rna_iptm" in df.columns and df["b2_protein_rna_iptm"].notna().any():
        return "b2_protein_rna_iptm", True
    return None, True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("outdir")
    ap.add_argument("--sort-col", default=None, help="column to sort designs by (desc)")
    ap.add_argument("--refolder", default=None, help="filter to one refolder (e.g. b2, af3); default: all rows")
    ap.add_argument("--msa", default=None, help="filter to one msa_variant (e.g. msaoff); default: all")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv, low_memory=False)
    ccols = sorted([c for c in df.columns if re.fullmatch(r"rna_contact_nt_\d+", c)],
                   key=lambda c: int(c.rsplit("_", 1)[1]))
    N = len(ccols)
    if N == 0:
        sys.exit("no rna_contact_nt_* columns in CSV — run the structural aggregation first")

    # optional refolder / msa filtering (contacts are per refold CIF)
    if args.refolder and "refolder" in df.columns:
        df = df[df["refolder"].astype(str).str.lower() == args.refolder.lower()]
    if args.msa and "msa_variant" in df.columns:
        df = df[df["msa_variant"].astype(str).str.lower() == args.msa.lower()]
    df = df[df[ccols].notna().any(axis=1)].copy()
    if df.empty:
        sys.exit("no rows with contact data after filtering")
    df[ccols] = df[ccols].fillna(0.0)

    sortcol, fell_back = pick_sort_col(df, args.sort_col)
    if fell_back and sortcol:
        print(f"WARNING: AF3 msaoff iPTM not available in CSV; sorting designs by {sortcol} instead.",
              file=sys.stderr)
    stacks = [s for s in pd.unique(df["stack"].dropna())]
    positions = np.arange(1, N + 1)

    # (1) per-stack heatmap: designs (sorted) x position
    fig, axes = plt.subplots(1, len(stacks), figsize=(3.2 * len(stacks) + 1.5, 6.5), squeeze=False)
    im = None
    vmax = float(np.nanpercentile(df[ccols].to_numpy(), 99)) or 1.0  # robust color scale
    for ax, s in zip(axes[0], stacks):
        sub = df[df["stack"] == s]
        if sortcol:
            sub = sub.sort_values(sortcol, ascending=False)
        mat = sub[ccols].to_numpy(dtype=float)
        im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_title(f"{s}\n(n={len(sub)})", fontsize=9)
        ax.set_xlabel("nt position", fontsize=8)
        ax.set_xticks(np.arange(N)); ax.set_xticklabels(positions, fontsize=6)
        ax.set_yticks([])
    axes[0][0].set_ylabel(f"designs (sorted by {sortcol} desc)" if sortcol else "designs")
    if im is not None:
        fig.colorbar(im, ax=axes[0].tolist(), label="heavy-atom contacts", fraction=0.025)
    fig.suptitle("Per-nucleotide protein-RNA contacts by stack (4.5 A heavy-atom)", fontsize=11)
    p1 = os.path.join(args.outdir, "rna_contact_heatmap.png")
    fig.savefig(p1, dpi=150, bbox_inches="tight")

    # (2) mean contacts per position, grouped by stack
    fig2, ax2 = plt.subplots(figsize=(max(8, N * 0.6), 5))
    w = 0.8 / max(len(stacks), 1)
    for i, s in enumerate(stacks):
        sub = df[df["stack"] == s]
        means = sub[ccols].mean().to_numpy()
        ax2.bar(positions + i * w - 0.4 + w / 2, means, width=w, label=f"{s} (n={len(sub)})")
    ax2.set_xlabel("nucleotide position"); ax2.set_ylabel("mean contacts / design")
    ax2.set_xticks(positions)
    ax2.set_title("Mean per-nucleotide protein-RNA contacts, by stack")
    ax2.legend(fontsize=8)
    p2 = os.path.join(args.outdir, "rna_contact_mean_bar.png")
    fig2.savefig(p2, dpi=150, bbox_inches="tight")

    print(f"wrote {p1}\nwrote {p2}\nstacks={stacks} N_nt={N} sort={sortcol} rows={len(df)}")


if __name__ == "__main__":
    main()
