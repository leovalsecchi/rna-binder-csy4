#!/usr/bin/env python3
"""plot_02 — CONSENSUS ipTM distribution per stack (one violin per stack, 5 total).

Each design contributes ONE point = min(AF3, Boltz-2, ESMFold2) ipTM — the conservative
"all three refolders must agree" read that matches the consensus ranking, NOT a per-refolder
pool or average. Only designs with all three refolds are shown. Static matplotlib +
interactive plotly (points click-to-open the design's 3 refolds).
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_02_violin_iptm"
KIND, LABEL, TMIN = "iptm", "ipTM", "iptm_min"
CLICK_JS = H.CLICK_JS   # violin points are click-to-open (design's 3 refolds in ChimeraX)

_IDCIF = ["design_id", "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]


def generate(df, cfg):
    t = D.refolder_triplet(KIND)            # per-design af3/b2/esm + cif paths (all-three only)
    if t.empty:
        return H.pending_pair(f"consensus {LABEL} by stack",
                              f"{LABEL} not available for all three refolders yet")
    t = t.copy()
    t["value"] = t[["af3", "b2", "esm"]].min(axis=1)     # consensus = the worst refolder
    stacks = [s for s in H.ordered_stacks(t) if (t["stack"] == s).any()]
    data = {s: t[t["stack"] == s]["value"].to_numpy() for s in stacks}
    thr = cfg.THRESHOLDS[TMIN]

    # ---- static (matplotlib) ----
    fig, ax = plt.subplots(figsize=(9, 5.4))
    parts = ax.violinplot([data[s] for s in stacks], showextrema=False, widths=0.85)
    for i, s in enumerate(stacks):
        body = parts["bodies"][i]
        body.set_facecolor(cfg.color(s)); body.set_edgecolor(cfg.MARK_EDGE)
        body.set_alpha(0.85); body.set_linewidth(1.0)
        q1, med, q3 = np.percentile(data[s], [25, 50, 75])
        ax.vlines(i + 1, q1, q3, color=cfg.INK_PRIMARY, lw=5, zorder=3)
        ax.scatter(i + 1, med, color="white", edgecolor=cfg.INK_PRIMARY, s=24, zorder=4)
    ax.set_xticks(range(1, len(stacks) + 1))
    ax.set_xticklabels([H.display(s) for s in stacks], rotation=18, ha="right", fontsize=9)
    ax.set_ylabel(f"consensus {LABEL}  =  min(AF3, Boltz-2, ESMFold2)")
    ax.set_title(f"Consensus {LABEL} by stack  (min across the 3 refolders)")
    ax.grid(axis="x", visible=False)
    H.add_mpl_threshold(ax, thr, "y", f"min {thr}")

    # ---- interactive (plotly): one point per design = its min; hover shows all 3 refolders ----
    pfig = go.Figure()
    for s in stacks:
        g = t[t["stack"] == s].copy()
        cd = g[_IDCIF + ["af3", "b2", "esm"]].copy()
        cd[_IDCIF] = cd[_IDCIF].fillna("")                # keep af3/b2/esm numeric for hover
        pfig.add_trace(go.Violin(
            y=g["value"], name=H.display(s), line_color=cfg.MARK_EDGE, fillcolor=cfg.color(s),
            opacity=0.85, box_visible=True, meanline_visible=True,
            points="all", jitter=0.5, pointpos=0, marker=dict(size=3, opacity=0.35),
            customdata=cd,
            hovertemplate=(f"{H.display(s)}<br>consensus {LABEL}=%{{y:.3f}}"
                           "<br>AF3=%{customdata[4]:.3f} · B2=%{customdata[5]:.3f}"
                           " · ESM=%{customdata[6]:.3f}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds</i><extra></extra>")))
    pfig.update_layout(
        title=f"Consensus {LABEL} by stack — min(AF3, Boltz-2, ESMFold2) · click a point → 3 refolds",
        yaxis_title=f"consensus {LABEL}", showlegend=False)
    H.add_plotly_threshold(pfig, thr, "y", f"min {thr}")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
