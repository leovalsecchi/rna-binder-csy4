#!/usr/bin/env python3
"""plot_08 — secondary-structure content per stack.

Mean helix_frac and sheet_frac per stack, computed twice: over ALL designs (solid)
and over only threshold-passing designs (hatched). Two panels (helix | sheet),
grouped bars per stack; texture is the secondary channel for all-vs-passing.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_08_secondary_structure"
PANELS = [("helix_frac", "Helix fraction"), ("sheet_frac", "Sheet fraction")]


def generate(df, cfg):
    if not (D.has_metric(df, "helix_frac") or D.has_metric(df, "sheet_frac")):
        return H.pending_pair("Secondary-structure content", "helix/sheet fractions not available")
    uniq = df.drop_duplicates("design_id")
    keep = D.passing_designs(df, "avg")
    stacks = H.ordered_stacks(uniq)
    npass = len(keep)

    def means(col):
        allm = {s: uniq.loc[uniq["stack"] == s, col].mean() for s in stacks}
        passm = {s: uniq[(uniq["stack"] == s) & uniq["design_id"].isin(keep)][col].mean()
                 for s in stacks}
        return allm, passm

    # ---- static ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    x = np.arange(len(stacks))
    for ax, (col, lab) in zip(axes, PANELS):
        allm, passm = means(col)
        ax.bar(x - 0.2, [allm[s] for s in stacks], width=0.38, color=[cfg.color(s) for s in stacks],
               edgecolor=cfg.MARK_EDGE, linewidth=0.4, label="all designs")
        ax.bar(x + 0.2, [passm[s] for s in stacks], width=0.38, color=[cfg.color(s) for s in stacks],
               edgecolor=cfg.MARK_EDGE, linewidth=0.6, hatch="///", label="threshold-passing")
        ax.set_xticks(x); ax.set_xticklabels([H.display(s) for s in stacks], rotation=20,
                                             ha="right", fontsize=8)
        ax.set_title(lab); ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("fraction")
    h = [plt.Rectangle((0, 0), 1, 1, fc="#bbb", ec=cfg.MARK_EDGE),
         plt.Rectangle((0, 0), 1, 1, fc="#bbb", ec=cfg.MARK_EDGE, hatch="///")]
    axes[1].legend(h, ["all designs", f"threshold-passing (n={npass})"], fontsize=8, framealpha=0.7)

    # ---- interactive ----
    pfig = make_subplots(rows=1, cols=2, subplot_titles=[lab for _, lab in PANELS])
    for j, (col, lab) in enumerate(PANELS, start=1):
        allm, passm = means(col)
        disp = [H.display(s) for s in stacks]
        pfig.add_trace(go.Bar(x=disp, y=[allm[s] for s in stacks], name="all",
                              marker_color=[cfg.color(s) for s in stacks],
                              marker_line_color=cfg.MARK_EDGE, marker_line_width=0.4,
                              showlegend=False), row=1, col=j)
        pfig.add_trace(go.Bar(x=disp, y=[passm[s] for s in stacks], name="passing",
                              marker_color=[cfg.color(s) for s in stacks],
                              marker_pattern_shape="/", marker_line_color=cfg.MARK_EDGE,
                              marker_line_width=0.4, showlegend=False), row=1, col=j)
    pfig.update_layout(title=f"Secondary-structure content  (solid = all, hatched = passing, n={npass})",
                       barmode="group")
    pfig.update_yaxes(title_text="fraction", row=1, col=1)
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
