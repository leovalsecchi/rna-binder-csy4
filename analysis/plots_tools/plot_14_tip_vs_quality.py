#!/usr/bin/env python3
"""plot_14 — tip score vs ipTM×DockQ quality (2-D precursor to the 3-D cube).

X = tip-weighted contact score; Y = consensus ipTM × consensus DockQ (the two
structural-agreement metrics combined into one quality axis — high only if both are high).
One point per design, colored by stack. Click a point → its 3 refolds open in ChimeraX.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_14_tip_vs_quality"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.quality_axes(stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("Tip vs ipTM×DockQ", "needs ipTM + DockQ + tip score")
    d = d.copy()
    d["quality"] = d["iptm"] * d["dockq"]
    stacks = H.ordered_stacks(d)
    imin, dmin, tmin = cfg.THRESHOLDS["iptm_min"], cfg.THRESHOLDS["dockq_min"], cfg.THRESHOLDS["tip_min"]
    qmin = imin * dmin                       # ipTM×DockQ ≥ this is NECESSARY for both to pass
    n_pass = int(((d["iptm"] >= imin) & (d["dockq"] >= dmin) & (d["tip"] >= tmin)).sum())
    xmax, ymax = float(d["tip"].max()), float(d["quality"].max())

    # ---- static ----
    fig, ax = plt.subplots(figsize=(8.6, 7))
    # light "pass all 3" box: tip ≥ tip_min & ipTM×DockQ ≥ iptm_min·dockq_min (approx on the product)
    ax.add_patch(Rectangle((tmin, qmin), xmax - tmin, ymax - qmin, facecolor="#2ca02c",
                           alpha=0.10, edgecolor="#2ca02c", linewidth=0.8, zorder=0))
    for s in stacks:
        g = d[d["stack"] == s]
        ax.scatter(g["tip"], g["quality"], s=12, alpha=0.5, color=cfg.color(s),
                   edgecolor="none", label=f"{H.display(s)} (n={len(g)})", zorder=3)
    ax.set_xlabel("tip-weighted contact score")
    ax.set_ylabel("ipTM × DockQ  (consensus)")
    ax.set_title(f"Tip score vs ipTM×DockQ quality\n(green box ≈ pass-all-3 zone; exact pass count = {n_pass})",
                 fontsize=11)
    ax.legend(fontsize=8, framealpha=0.7, loc="upper right")

    # ---- interactive ----
    pfig = go.Figure()
    pfig.add_shape(type="rect", x0=tmin, x1=xmax, y0=qmin, y1=ymax, layer="below",
                   fillcolor="rgba(44,160,44,0.10)", line=dict(color="rgba(44,160,44,0.5)", width=1))
    pfig.add_annotation(x=tmin, y=ymax, text=f"≈ pass all 3 (n={n_pass})", showarrow=False,
                        xanchor="left", yanchor="top", font=dict(color="#2ca02c", size=11))
    for s in stacks:
        g = d[d["stack"] == s]
        pfig.add_trace(go.Scattergl(
            x=g["tip"], y=g["quality"], mode="markers", name=H.display(s),
            marker=dict(color=cfg.color(s), size=6, opacity=0.6,
                        line=dict(color=cfg.MARK_EDGE, width=0.4)),
            customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path",
                          "esm_cif_path"]].fillna(""),
            hovertemplate=(f"{H.display(s)}<br>tip=%{{x:.1f}}  ipTM×DockQ=%{{y:.3f}}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds (RNA-aligned)"
                           " in ChimeraX</i><extra></extra>")))
    pfig.update_layout(title="Tip score vs ipTM×DockQ quality  (click a point → 3 refolds)",
                       xaxis_title="tip-weighted contact score",
                       yaxis_title="ipTM × DockQ (consensus)")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
