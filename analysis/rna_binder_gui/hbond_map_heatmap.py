#!/usr/bin/env python3
"""Heatmap of protein->RNA H-bonds: amino-acid type (y) x nucleotide position (x).

Reads the long CSV from hbond_map_pyrosetta.py (cols: stack,refolder,msa_variant,
abs_path,aa,nt_pos,n_hbonds) and draws, per stack, a heatmap with:
    x = nucleotide position (1..N, optionally labeled with the RNA sequence)
    y = amino-acid type (20)
    cell = mean H-bonds of that aa->nt per design (sum over designs / n designs in stack)
A design here = one unique structure (abs_path). Filter refolder/msa if desired.

Usage:
    hbond_map_heatmap.py HBOND_MAP_CSV OUTDIR [--refolder af3] [--msa msaoff]
                         [--rna-seq ACUGCCGUAUAGGCAG] [--metric mean|total]
"""
import argparse, os, sys
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/mplconfig_{os.getpid()}")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

AA_ORDER = list("RHKDESTNQCGPAVILMFYW")   # charged/polar -> hydrophobic, fixed 20-row axis


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv"); ap.add_argument("outdir")
    ap.add_argument("--refolder", default=None); ap.add_argument("--msa", default=None)
    ap.add_argument("--rna-seq", default="ACUGCCGUAUAGGCAG")
    ap.add_argument("--metric", choices=["mean", "total"], default="mean")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)
    if args.refolder: df = df[df["refolder"].astype(str).str.lower() == args.refolder.lower()]
    if args.msa: df = df[df["msa_variant"].astype(str).str.lower() == args.msa.lower()]
    if df.empty: sys.exit("no rows after filtering")
    N = int(df["nt_pos"].max())
    positions = list(range(1, N + 1))
    seq = args.rna_seq if len(args.rna_seq) == N else None
    stacks = [s for s in pd.unique(df["stack"].dropna())]

    fig, axes = plt.subplots(1, len(stacks), figsize=(3.4 * len(stacks) + 1.5, 6.5),
                             squeeze=False, sharey=True)
    # shared color scale across stacks
    mats = {}
    for s in stacks:
        sub = df[df["stack"] == s]
        ndes = sub["abs_path"].nunique()
        piv = (sub.groupby(["aa", "nt_pos"])["n_hbonds"].sum()
               .unstack(fill_value=0).reindex(index=AA_ORDER, columns=positions, fill_value=0))
        m = piv.to_numpy(dtype=float)
        if args.metric == "mean" and ndes:
            m = m / ndes
        mats[s] = (m, ndes)
    vmax = max((mat.max() for mat, _ in mats.values()), default=1.0) or 1.0
    im = None
    for ax, s in zip(axes[0], stacks):
        m, ndes = mats[s]
        im = ax.imshow(m, aspect="auto", cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_title(f"{s}\n(n={ndes} designs)", fontsize=9)
        ax.set_xticks(range(N))
        ax.set_xticklabels([f"{p}\n{seq[p-1]}" if seq else str(p) for p in positions], fontsize=6)
        ax.set_xlabel("nucleotide position", fontsize=8)
    axes[0][0].set_yticks(range(len(AA_ORDER))); axes[0][0].set_yticklabels(AA_ORDER, fontsize=7)
    axes[0][0].set_ylabel("amino acid")
    fig.colorbar(im, ax=axes[0].tolist(),
                 label=("mean H-bonds / design" if args.metric == "mean" else "total H-bonds"),
                 fraction=0.025)
    fig.suptitle("Protein->RNA H-bond map (PyRosetta): amino acid x nucleotide, by stack", fontsize=11)
    p1 = os.path.join(args.outdir, "hbond_map_aa_x_nt.png")
    fig.savefig(p1, dpi=150, bbox_inches="tight")
    print(f"wrote {p1}\nstacks={stacks} N_nt={N} metric={args.metric}")


if __name__ == "__main__":
    main()
