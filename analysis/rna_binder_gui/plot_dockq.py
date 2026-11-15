#!/usr/bin/env python3
"""
plot_dockq.py
DockQ (binder-RNA pose quality) plots, AF3 msaoff, best sample per design.

Two figures:
  1. dockq_lrmsd_vs_carmsd.png  — pose (DockQ ligand-RMSD) vs fold (Cα-RMSD).
     Exposes designs that fold as intended (low Cα-RMSD) but place the RNA wrong
     (high LRMSD) — the trap that Cα-RMSD alone misses.
  2. dockq_quality_by_stack.png — CAPRI quality breakdown per stack from the
     DockQ score (incorrect <0.23, acceptable, medium, high >0.80).

Usage:
    python3 plot_dockq.py --csv <path_to_csv> --outdir <output_dir>
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STACK_ORDER = ['RFD3+LiMPNN', 'BG+LiMPNN', 'BG+BIF', 'BD1+LiMPNN', 'HD+LiMPNN']
STACK_COLORS = {
    'RFD3+LiMPNN': '#9C27B0',
    'BG+LiMPNN':   '#2196F3',
    'BG+BIF':      '#64B5F6',
    'BD1+LiMPNN':  '#FF9800',
    'HD+LiMPNN':   '#F44336',
}
# CAPRI DockQ-score quality bins.
CAPRI_BINS = [('incorrect', 0.00, 0.23, '#BDBDBD'),
              ('acceptable', 0.23, 0.49, '#FFB74D'),
              ('medium', 0.49, 0.80, '#4FC3F7'),
              ('high', 0.80, 1.01, '#388E3C')]


def normalize_refold(row):
    if row['refolder'] == 'B2':
        return 'B2'
    if row['refold_label'] in ['msaon', 'prod_bg200_bif', 'prod_hd200_v2']:
        return 'AF3_msaon'
    return 'AF3_msaoff'


def best_af3_msaoff(df):
    """AF3 msaoff, best sample per design by af3_iptm."""
    df = df.copy()
    df['refold_type'] = df.apply(normalize_refold, axis=1)
    af3 = df[(df['refolder'] == 'AF3') & (df['refold_type'] == 'AF3_msaoff')].copy()
    af3['_design_key'] = af3['file'].str.replace(
        r'_AF3_seed\d+_sample\d+.*', '', regex=True)
    return af3.loc[af3.groupby(['stack', '_design_key'])['af3_iptm'].idxmax()]


def plot_lrmsd_vs_carmsd(best, outdir, rmsd_max=15.0, lrmsd_max=30.0):
    fig, ax = plt.subplots(figsize=(8, 6))
    for stack in STACK_ORDER:
        sub = best[best['stack'] == stack][['ca_rmsd_vs_designed', 'dockq_lrmsd']].dropna()
        if sub.empty:
            continue
        ax.scatter(sub['ca_rmsd_vs_designed'].clip(upper=rmsd_max),
                   sub['dockq_lrmsd'].clip(upper=lrmsd_max),
                   color=STACK_COLORS[stack], alpha=0.4, s=12,
                   label=f'{stack} (n={len(sub)})', zorder=3)
    ax.axvline(2.0, color='darkred', linestyle='--', linewidth=1, alpha=0.6)   # fold filter
    ax.axhline(10.0, color='navy', linestyle='--', linewidth=1, alpha=0.6)     # rough pose cutoff
    ax.fill_between([0, 2.0], 0, 10.0, alpha=0.06, color='green', zorder=1)     # folds + docks
    ax.text(0.30, 0.06, 'folds right\n+ docks right', fontsize=8, ha='center',
            color='darkgreen', alpha=0.8, transform=ax.transAxes)
    ax.text(0.30, 0.62, 'folds right\nbut wrong pose', fontsize=8, ha='center',
            color='crimson', alpha=0.8, transform=ax.transAxes)
    ax.set_xlabel('Cα RMSD vs designed (Å)  — fold self-consistency', fontsize=11)
    ax.set_ylabel('DockQ ligand-RMSD (Å)  — RNA pose vs designed', fontsize=11)
    ax.set_title('Pose quality vs fold quality (AF3 msaoff, best sample/design)',
                 fontsize=12, fontweight='bold')
    ax.set_xlim(0, rmsd_max * 1.05)
    ax.set_ylim(0, lrmsd_max * 1.05)
    ax.legend(fontsize=8, markerscale=2, loc='upper right')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = outdir / 'dockq_lrmsd_vs_carmsd.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def plot_quality_by_stack(best, outdir):
    stacks = [s for s in STACK_ORDER if (best['stack'] == s).any()]
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(stacks))
    bottom = np.zeros(len(stacks))
    for name, lo, hi, color in CAPRI_BINS:
        counts = []
        for s in stacks:
            d = best[best['stack'] == s]['dockq'].dropna()
            counts.append(int(((d >= lo) & (d < hi)).sum()))
        ax.bar(x, counts, bottom=bottom, color=color, label=name, zorder=3)
        bottom += np.array(counts)
    ax.set_xticks(x)
    ax.set_xticklabels(stacks, rotation=20, ha='right', fontsize=9)
    ax.set_ylabel('designs (best AF3 msaoff sample)', fontsize=11)
    ax.set_title('DockQ CAPRI quality by stack', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, title='DockQ quality')
    ax.grid(alpha=0.3, axis='y')
    fig.tight_layout()
    out = outdir / 'dockq_quality_by_stack.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    if 'dockq' not in df.columns:
        print("WARNING: no 'dockq' column in CSV — skipping DockQ plots")
        return
    best = best_af3_msaoff(df)
    n_dockq = int(best['dockq'].notna().sum())
    print(f"{len(best)} designs (AF3 msaoff best), {n_dockq} with DockQ")
    plot_lrmsd_vs_carmsd(best, args.outdir)
    plot_quality_by_stack(best, args.outdir)


if __name__ == '__main__':
    main()
