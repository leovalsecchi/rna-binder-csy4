#!/usr/bin/env python3
"""plot_04 — consensus ipTM vs consensus DockQ scatter, colored by stack.

ONE point per design at its CONSENSUS values: ipTM = min(AF3, Boltz-2, ESMFold2),
DockQ = min across the available refolders. Same definition as the strict funnel
(plot_10) and quality cube (plot_15), so the top-right "pass" quadrant holds the
same designs (n_pass matches). Only designs refolded by all conditions appear.
Points top-right of both thresholds (ipTM ≥ iptm_min, DockQ ≥ dockq_min) are
full-opacity, the rest faded. Pending until DockQ exists.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_04_scatter_iptm_dockq"
CLICK_JS = H.CLICK_JS   # points are click-to-open: a click opens the design's 3 refolds


def generate(df, cfg):
    if not (D.has_metric(df, "iptm") and D.has_metric(df, "dockq")):
        return H.pending_pair("Consensus ipTM vs DockQ", "needs both ipTM and DockQ (DockQ pending)")
    imin, dmin = cfg.THRESHOLDS["iptm_min"], cfg.THRESHOLDS["dockq_min"]
    # consensus per design: iptm = min(af3,b2,esm), dockq = min across refolders (same as cube)
    d = D.quality_axes(stacks=list(df["stack"].unique()), require=("iptm", "dockq"))
    if d.empty:
        return H.pending_pair("Consensus ipTM vs DockQ", "needs both ipTM and DockQ")
    d = d.copy()
    d["pass"] = (d["iptm"] >= imin) & (d["dockq"] >= dmin)
    for c in ("af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"):
        if c not in d.columns:
            d[c] = ""
    stacks = H.ordered_stacks(d)
    n_pass = int(d["pass"].sum())

    # ---- static ----
    fig, ax = plt.subplots(figsize=(8.2, 7))
    ax.axhspan(dmin, 1.0, xmin=0, xmax=1, color="#f4f4f0", zorder=0)
    for s in stacks:
        for passed, alpha, sz in [(False, 0.18, 12), (True, 0.95, 26)]:
            sub = d[(d["stack"] == s) & (d["pass"] == passed)]
            ax.scatter(sub["iptm"], sub["dockq"], s=sz, alpha=alpha, color=cfg.color(s),
                       edgecolor=cfg.MARK_EDGE if passed else "none",
                       linewidth=0.5, zorder=3 if passed else 2)
    H.add_mpl_threshold(ax, imin, "x", f"ipTM {imin}")
    H.add_mpl_threshold(ax, dmin, "y", f"DockQ {dmin}")
    ax.set_xlabel("consensus ipTM  (min AF3 / Boltz-2 / ESMFold2)")
    ax.set_ylabel("consensus DockQ  (min across refolders)")
    ax.set_title(f"Consensus ipTM vs DockQ  (top-right = pass both, n={n_pass})")
    handles = [plt.Line2D([], [], marker="o", ls="", color=cfg.color(s), label=H.display(s))
               for s in stacks]
    ax.legend(handles=handles, fontsize=8, framealpha=0.7, loc="lower left")

    # ---- interactive ----
    # ONE point per design at its consensus (min) ipTM/DockQ — matches plot_10 / plot_15.
    pfig = go.Figure()
    for s in stacks:
        sub = d[d["stack"] == s]
        pfig.add_trace(go.Scattergl(
            x=sub["iptm"], y=sub["dockq"], mode="markers", name=H.display(s),
            marker=dict(color=cfg.color(s), size=[9 if p else 5 for p in sub["pass"]],
                        opacity=0.8, line=dict(color=cfg.MARK_EDGE, width=0.5)),
            customdata=sub[["design_id", "af3_msaoff_cif_path", "b2_cif_path",
                            "esm_cif_path"]].fillna(""),
            hovertemplate=(f"{H.display(s)}<br>consensus ipTM=%{{x:.3f}}  "
                           "DockQ=%{{y:.3f}}<br>%{customdata[0]}"
                           "<br><i>click → open 3 refolds (RNA-aligned) in ChimeraX</i><extra></extra>")))
    H.add_plotly_threshold(pfig, imin, "x", f"ipTM {imin}")
    H.add_plotly_threshold(pfig, dmin, "y", f"DockQ {dmin}")
    pfig.update_layout(
        title=f"Consensus ipTM vs DockQ  (top-right = pass both, n={n_pass} · click a point → 3 refolds)",
        xaxis_title="consensus ipTM", yaxis_title="consensus DockQ")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
