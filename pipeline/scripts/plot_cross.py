#!/usr/bin/env python3
"""plot_cross.py — interactive cross-campaign plots from the combined sample10k.

Self-contained (no dependency on the v2-era analysis.py); consumes the new
per-design sample10k_metrics schema (af3_msaoff_iptm / b2_iptm / esm_iptm /
ca_rmsd_vs_designed / helix_fraction / rna_contact_nt_1..16 / stack ...).

Produces under OUTDIR:
    af3_vs_b2_iptm.html        scatter, colour=stack, y=x consensus diagonal
    consensus_hist.html        histogram of min(AF3,B2) per stack
    iptm_by_refolder.html      AF3 / B2 / ESM iPTM histograms (subplots)
    violins_by_stack.html      iPTM + Ca-RMSD + helix per stack
    rna_contact_mean_bar.png   mean protein->RNA contacts per nt, grouped by stack
    rna_contact_heatmap.png    per-stack design x nt contact heatmap (sorted by consensus)
    index.html

Usage:  plot_cross.py COMBINED_SAMPLE3K_CSV OUTDIR
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", f"/tmp/mplconfig_{os.getpid()}")
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AF3, B2, ESM = "af3_msaoff_iptm", "b2_iptm", "esm_iptm"
STACK_COLORS = {
    "RFD3+LiMPNN": "#9C27B0", "HD+LiMPNN": "#F44336",
    "BG+LiMPNN": "#2196F3", "BG+BIF": "#64B5F6", "BD1+LiMPNN": "#FF9800",
}
_FALLBACK = ["#26A69A", "#EC407A", "#7E57C2", "#8D6E63", "#789262", "#FF7043"]


def stacks_present(df: pd.DataFrame):
    present = list(df["stack"].dropna().unique())
    known = [s for s in STACK_COLORS if s in present]
    extra = sorted(set(present) - set(STACK_COLORS))
    out = [(s, STACK_COLORS[s]) for s in known]
    out += [(s, _FALLBACK[i % len(_FALLBACK)]) for i, s in enumerate(extra)]
    return out


def scatter_af3_b2(df, out: Path):
    sub = df.dropna(subset=[AF3, B2])
    if sub.empty:
        print("  (skip af3_vs_b2: no data)"); return
    fig = go.Figure()
    for stk, col in stacks_present(sub):
        g = sub[sub["stack"] == stk]
        fig.add_trace(go.Scattergl(
            x=g[B2], y=g[AF3], mode="markers", name=f"{stk} (n={len(g)})",
            marker=dict(size=6, color=col, opacity=0.55),
            customdata=g["design_id"],
            hovertemplate="<b>%{customdata}</b><br>B2=%{x:.3f}<br>AF3=%{y:.3f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", name="AF3=B2",
                             line=dict(dash="dash", color="gray"), showlegend=True))
    fig.update_layout(title="AF3-msaoff vs B2 interface iPTM (consensus = closeness to top-right)",
                      xaxis_title="B2 iPTM", yaxis_title="AF3-msaoff iPTM",
                      template="plotly_white", height=680, legend_title="stack")
    fig.update_xaxes(range=[0, 1]); fig.update_yaxes(range=[0, 1])
    fig.write_html(out / "af3_vs_b2_iptm.html", include_plotlyjs="cdn")
    print(f"  wrote af3_vs_b2_iptm.html ({len(sub)} pts)")


def consensus_hist(df, out: Path):
    if "consensus_iptm" not in df.columns:
        df = df.assign(consensus_iptm=df[[AF3, B2]].min(axis=1))
    sub = df.dropna(subset=["consensus_iptm"])
    if sub.empty:
        print("  (skip consensus_hist)"); return
    fig = go.Figure()
    for stk, col in stacks_present(sub):
        g = sub[sub["stack"] == stk]
        fig.add_trace(go.Histogram(x=g["consensus_iptm"], name=f"{stk} (n={len(g)})",
                                   marker_color=col, opacity=0.6, nbinsx=40))
    fig.update_layout(title="Consensus iPTM  min(AF3, B2)  by stack",
                      xaxis_title="min(AF3, B2) iPTM", yaxis_title="designs",
                      barmode="overlay", template="plotly_white", height=560)
    fig.write_html(out / "consensus_hist.html", include_plotlyjs="cdn")
    print(f"  wrote consensus_hist.html ({len(sub)} designs)")


def iptm_by_refolder(df, out: Path):
    cols = [(AF3, "AF3-msaoff"), (B2, "B2"), (ESM, "ESM")]
    cols = [(c, lab) for c, lab in cols if c in df.columns and df[c].notna().any()]
    if not cols:
        print("  (skip iptm_by_refolder)"); return
    fig = make_subplots(rows=1, cols=len(cols), subplot_titles=[l for _, l in cols])
    seen = set()
    for ci, (c, lab) in enumerate(cols, 1):
        sub = df.dropna(subset=[c])
        for stk, col in stacks_present(sub):
            g = sub[sub["stack"] == stk]
            fig.add_trace(go.Histogram(x=g[c], name=stk, legendgroup=stk,
                                       showlegend=stk not in seen, marker_color=col,
                                       opacity=0.6, nbinsx=35), row=1, col=ci)
            seen.add(stk)
        fig.update_xaxes(title_text=f"{lab} iPTM", row=1, col=ci)
    fig.update_layout(title="Interface iPTM by refolder", barmode="overlay",
                      template="plotly_white", height=520, legend_title="stack")
    fig.write_html(out / "iptm_by_refolder.html", include_plotlyjs="cdn")
    print("  wrote iptm_by_refolder.html")


def violins(df, out: Path):
    metrics = [(AF3, "AF3 iPTM", None), (B2, "B2 iPTM", None), (ESM, "ESM iPTM", None),
               ("ca_rmsd_vs_designed", "Ca-RMSD vs designed (A)", 15),
               ("helix_fraction", "helix fraction", None)]
    metrics = [(c, l, cl) for c, l, cl in metrics if c in df.columns]
    if not metrics:
        print("  (skip violins)"); return
    fig = make_subplots(rows=1, cols=len(metrics), subplot_titles=[l for _, l, _ in metrics])
    seen = set()
    for ci, (c, lab, clip) in enumerate(metrics, 1):
        sub = df.dropna(subset=[c])
        for stk, col in stacks_present(sub):
            g = sub[sub["stack"] == stk]
            if g.empty:
                continue
            y = g[c].clip(upper=clip) if clip else g[c]
            fig.add_trace(go.Violin(y=y, name=stk, legendgroup=stk, showlegend=stk not in seen,
                                    scalegroup=c, line_color=col, fillcolor=col, opacity=0.6,
                                    box_visible=True, meanline_visible=True, points=False),
                          row=1, col=ci)
            seen.add(stk)
        fig.update_yaxes(title_text=lab, row=1, col=ci)
        fig.update_xaxes(showticklabels=False, row=1, col=ci)
    fig.update_layout(title="Distributions by stack", template="plotly_white",
                      height=600, legend_title="stack", violinmode="group")
    fig.write_html(out / "violins_by_stack.html", include_plotlyjs="cdn")
    print(f"  wrote violins_by_stack.html ({len(seen)} stacks)")


def contact_cols(df):
    cs = [c for c in df.columns if c.startswith("rna_contact_nt_")]
    return sorted(cs, key=lambda c: int(c.rsplit("_", 1)[1]))


def rna_contact_plots(df, out: Path):
    cc = contact_cols(df)
    if not cc:
        print("  (skip rna_contact: no rna_contact_nt_* cols)"); return
    nts = [int(c.rsplit("_", 1)[1]) for c in cc]
    order = [s for s, _ in stacks_present(df)]
    # (1) mean contacts per nt, grouped by stack
    fig, ax = plt.subplots(figsize=(11, 5))
    w = 0.8 / max(len(order), 1)
    for i, stk in enumerate(order):
        g = df[df["stack"] == stk]
        means = [g[c].mean() for c in cc]
        ax.bar(np.arange(len(nts)) + i * w, means, width=w,
               label=f"{stk} (n={len(g)})", color=STACK_COLORS.get(stk, _FALLBACK[i % len(_FALLBACK)]))
    ax.set_xticks(np.arange(len(nts)) + 0.4 - w / 2)
    ax.set_xticklabels(nts)
    ax.set_xlabel("RNA nucleotide position"); ax.set_ylabel("mean heavy-atom contacts (4.5 A)")
    ax.set_title("Mean protein->RNA contacts per nucleotide, by stack")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "rna_contact_mean_bar.png", dpi=130); plt.close(fig)
    print("  wrote rna_contact_mean_bar.png")
    # (2) per-stack heatmap, designs sorted by consensus (or AF3)
    sortcol = "consensus_iptm" if "consensus_iptm" in df.columns else AF3
    fig, axes = plt.subplots(1, len(order), figsize=(3.1 * len(order), 6), squeeze=False)
    for ax, stk in zip(axes[0], order):
        g = df[df["stack"] == stk].sort_values(sortcol, ascending=False)
        M = g[cc].to_numpy(dtype=float)
        im = ax.imshow(M, aspect="auto", cmap="viridis", interpolation="nearest")
        ax.set_title(stk, fontsize=9); ax.set_xlabel("nt")
        ax.set_xticks(range(len(nts))); ax.set_xticklabels(nts, fontsize=6)
        ax.set_yticks([])
    axes[0][0].set_ylabel(f"designs (sorted by {sortcol} desc)")
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.6, label="contacts")
    fig.suptitle("Per-nucleotide contact map by stack")
    fig.savefig(out / "rna_contact_heatmap.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    print("  wrote rna_contact_heatmap.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("outdir")
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.csv, low_memory=False)
    print(f"{len(df):,} sampled designs; stacks: {sorted(df['stack'].dropna().unique())}")
    scatter_af3_b2(df, out)
    consensus_hist(df, out)
    iptm_by_refolder(df, out)
    violins(df, out)
    rna_contact_plots(df, out)
    links = "\n".join(f'<li><a href="{p.name}">{p.name}</a></li>'
                      for p in sorted(out.glob("*.html")) if p.name != "index.html")
    imgs = "\n".join(f'<h3>{p.name}</h3><img src="{p.name}" style="max-width:100%">'
                     for p in sorted(out.glob("*.png")))
    (out / "index.html").write_text(
        f"<h2>Cross-campaign plots</h2><p>HD+RFD3 (v3) + BG+BD (v2), "
        f"ranked by AF3&amp;B2 consensus.</p><ul>{links}</ul>{imgs}")
    print(f"index -> {out/'index.html'}")


if __name__ == "__main__":
    main()
