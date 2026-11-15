#!/usr/bin/env python3
"""plot_21 — ipTM distribution per refolder (three violins: AF3 / Boltz-2 / ESMFold2).

One violin per refolder over all designs scored by all three refolders; each design contributes
one point per refolder = that refolder's ipTM. Interactive points are click-to-open (the design's
3 refolds, RNA-aligned, in ChimeraX). Complements plot_02 (per-stack consensus min) and plot_07
(per-refolder means): here you see the full shape of each refolder's ipTM distribution + the dots.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_23_violin_iptm_refolder"
LABEL, TMIN = "ipTM", "iptm_min"
CLICK_JS = H.CLICK_JS   # violin points are click-to-open (design's 3 refolds in ChimeraX)

_IDCIF = ["design_id", "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]


def generate(df, cfg):
    t = D.refolder_triplet("iptm")          # per-design af3/b2/esm + cif paths (all-three only)
    if t.empty:
        return H.pending_pair(f"{LABEL} per refolder",
                              f"{LABEL} not available for all three refolders yet")
    refs = [r for r in cfg.REFOLDER_ORDER if r in t.columns]     # af3, b2, esm
    labels = [cfg.REFOLDER_LABELS[r] for r in refs]
    thr = cfg.THRESHOLDS[TMIN]

    # ---- static (matplotlib) ----
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    parts = ax.violinplot([t[r].dropna().to_numpy() for r in refs], showextrema=False, widths=0.85)
    for i, r in enumerate(refs):
        body = parts["bodies"][i]
        body.set_facecolor(cfg.REFOLDER_COLORS[r]); body.set_edgecolor(cfg.MARK_EDGE)
        body.set_alpha(0.85); body.set_linewidth(1.0)
        vals = t[r].dropna().to_numpy()
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.vlines(i + 1, q1, q3, color=cfg.INK_PRIMARY, lw=5, zorder=3)
        ax.scatter(i + 1, med, color="white", edgecolor=cfg.INK_PRIMARY, s=24, zorder=4)
    ax.set_xticks(range(1, len(refs) + 1))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(LABEL)
    ax.set_title(f"{LABEL} distribution per refolder  (n = {len(t):,} designs)")
    ax.grid(axis="x", visible=False)
    H.add_mpl_threshold(ax, thr, "y", f"{LABEL} {thr}")

    # ---- interactive (plotly): one violin per refolder; points click-to-open ----
    pfig = go.Figure()
    for r in refs:
        g = t[t[r].notna()].copy()
        cd = g[_IDCIF + ["af3", "b2", "esm"]].copy()
        cd[_IDCIF] = cd[_IDCIF].fillna("")               # keep af3/b2/esm numeric for hover
        pfig.add_trace(go.Violin(
            y=g[r], name=cfg.REFOLDER_LABELS[r], line_color=cfg.MARK_EDGE,
            fillcolor=cfg.REFOLDER_COLORS[r], opacity=0.85, box_visible=True, meanline_visible=True,
            points="all", jitter=0.5, pointpos=0,
            marker=dict(size=3, opacity=0.35, color=cfg.REFOLDER_COLORS[r]),
            customdata=cd,
            hovertemplate=(f"{cfg.REFOLDER_LABELS[r]}<br>{LABEL}=%{{y:.3f}}"
                           "<br>AF3=%{customdata[4]:.3f} · B2=%{customdata[5]:.3f}"
                           " · ESM=%{customdata[6]:.3f}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds</i><extra></extra>")))
    pfig.update_layout(
        title=f"{LABEL} per refolder — AF3 / Boltz-2 / ESMFold2 · click a point → 3 refolds",
        yaxis_title=LABEL, xaxis_title="refolder", showlegend=False)
    H.add_plotly_threshold(pfig, thr, "y", f"{LABEL} {thr}")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
