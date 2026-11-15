#!/usr/bin/env python3
"""plot_15 — quality cube: consensus ipTM × consensus DockQ × tip score, as a 3-D scatter.

The interactive HTML is rotatable; a click on any point opens that design's 3 refolds in
ChimeraX (RNA-aligned). The PNG is a fixed-view 3-D projection. One point per design, by stack.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_15_quality_cube"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.quality_axes(stacks=list(df["stack"].unique()), require=("iptm", "dockq"))
    if d.empty:
        return H.pending_pair("Quality cube (ipTM × DockQ × tip)", "needs ipTM + DockQ")
    stacks = H.ordered_stacks(d)
    # pass box: ipTM ≥ iptm_min AND DockQ ≥ dockq_min. tip is shown as an axis but does NOT
    # gate passing — keeps this equal-strict with the strict funnel (plot_10), which filters
    # on ipTM + DockQ only. Box spans the full tip range.
    imin, dmin = cfg.THRESHOLDS["iptm_min"], cfg.THRESHOLDS["dockq_min"]
    _has_tip = bool(d["tip"].notna().any())
    tlo = float(d["tip"].min()) if _has_tip else 0.0
    x1, y1, z1 = 1.0, 1.0, (float(d["tip"].max()) if _has_tip else 1.0)
    passmask = (d["iptm"] >= imin) & (d["dockq"] >= dmin)
    n_pass = int(passmask.sum())
    box_name = f"pass both ({n_pass}): ipTM≥{imin:g}·DockQ≥{dmin:g}"

    # ---- static (fixed-view 3-D) ----
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    for s in stacks:
        g = d[d["stack"] == s]
        ax.scatter(g["iptm"], g["dockq"], g["tip"], s=8, alpha=0.5,
                   color=cfg.color(s), edgecolor="none", label=H.display(s))
    # shaded pass box (6 faces, light green)
    c = [(imin, dmin, tlo), (imin, dmin, z1), (imin, y1, tlo), (imin, y1, z1),
         (x1, dmin, tlo), (x1, dmin, z1), (x1, y1, tlo), (x1, y1, z1)]
    faces = [[c[0], c[1], c[3], c[2]], [c[4], c[5], c[7], c[6]], [c[0], c[1], c[5], c[4]],
             [c[2], c[3], c[7], c[6]], [c[0], c[2], c[6], c[4]], [c[1], c[3], c[7], c[5]]]
    ax.add_collection3d(Poly3DCollection(faces, facecolor="#2ca02c", alpha=0.10,
                                         edgecolor="#2ca02c", linewidths=0.5))
    ax.set_xlabel("consensus ipTM"); ax.set_ylabel("consensus DockQ"); ax.set_zlabel("tip score")
    ax.set_title(f"Quality cube — ipTM × DockQ × tip\n(green box = pass ipTM+DockQ, n={n_pass})", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.view_init(elev=22, azim=-60)

    # ---- interactive (rotatable, clickable) ----
    pfig = go.Figure()
    for s in stacks:
        g = d[d["stack"] == s]
        pfig.add_trace(go.Scatter3d(
            x=g["iptm"], y=g["dockq"], z=g["tip"], mode="markers", name=H.display(s),
            marker=dict(color=cfg.color(s),
                        size=[7 if p else 3 for p in passmask.loc[g.index]],  # only passers enlarged
                        opacity=0.7, line=dict(color=cfg.MARK_EDGE, width=0.3)),
            customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path",
                          "esm_cif_path"]].fillna(""),
            hovertemplate=(f"{H.display(s)}<br>ipTM=%{{x:.3f}}  DockQ=%{{y:.3f}}  tip=%{{z:.1f}}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds</i><extra></extra>")))
    # "pass all three" box as a WIREFRAME (a filled mesh would swallow point clicks in 3-D)
    ex, ey, ez = H.box_edges(imin, x1, dmin, y1, tlo, z1)
    pfig.add_trace(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines",
                                line=dict(color="#1a7d1a", width=7),
                                name=box_name, hoverinfo="skip", showlegend=True))
    pfig.update_layout(
        title="Quality cube — ipTM × DockQ × tip  (green box = pass ipTM+DockQ · drag to rotate · click → 3 refolds)",
        scene=dict(xaxis_title="consensus ipTM", yaxis_title="consensus DockQ",
                   zaxis_title="tip score"),
        margin=dict(l=0, r=0, t=50, b=0), height=760)
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
