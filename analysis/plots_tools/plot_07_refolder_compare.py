#!/usr/bin/env python3
"""plot_07 — refolder ipTM comparison per stack.

Mean ipTM from each MSA-off refolder (AF3, Boltz-2, ESMFold2) per stack, grouped
bars colored by refolder; the highest-scoring refolder per stack is marked (★ +
bold edge) so 'which refolder scores highest' is obvious at a glance.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_07_refolder_compare"


def generate(df, cfg):
    d = D.refolder_iptm(stacks=list(df["stack"].unique()))
    if d.empty or not d["iptm"].notna().any():
        return H.pending_pair("Refolder ipTM comparison", "no refolder ipTM available")
    piv = d.groupby(["stack", "refolder"])["iptm"].mean().unstack("refolder")
    stacks = [s for s in cfg.STACK_ORDER if s in piv.index]
    refs = [r for r in cfg.REFOLDER_ORDER if r in piv.columns and piv[r].notna().any()]
    best = {s: max(refs, key=lambda r: (piv.loc[s, r] if not np.isnan(piv.loc[s, r]) else -1))
            for s in stacks}
    title = "Mean ipTM per refolder, by stack   (★ = highest-scoring refolder)"

    # ---- static ----
    x = np.arange(len(stacks))
    w = 0.8 / len(refs)
    fig, ax = plt.subplots(figsize=(11, 5.6))
    for j, r in enumerate(refs):
        vals = [piv.loc[s, r] for s in stacks]
        bars = ax.bar(x + (j - (len(refs) - 1) / 2) * w, vals, width=w,
                      color=cfg.REFOLDER_COLORS[r], edgecolor=cfg.MARK_EDGE, linewidth=0.5,
                      label=cfg.REFOLDER_LABELS[r])
        for s, b, v in zip(stacks, bars, vals):
            if np.isnan(v):
                continue
            ax.text(b.get_x() + b.get_width() / 2, v + 0.006, f"{v:.2f}", ha="center",
                    va="bottom", fontsize=7.5, color=cfg.INK_SECONDARY)
            if best[s] == r:
                b.set_edgecolor(cfg.INK_PRIMARY); b.set_linewidth(2)
                ax.text(b.get_x() + b.get_width() / 2, v + 0.045, "★", ha="center",
                        va="bottom", fontsize=13, color=cfg.INK_PRIMARY)
    ax.set_xticks(x); ax.set_xticklabels([H.display(s) for s in stacks], rotation=18,
                                         ha="right", fontsize=9)
    ax.set_ylabel("mean ipTM"); ax.set_title(title); ax.grid(axis="x", visible=False)
    ax.legend(title="Refolder", fontsize=9, framealpha=0.9, loc="upper left",
              bbox_to_anchor=(1.005, 1.0))
    ax.set_ylim(0, min(1.0, float(np.nanmax(piv.values)) + 0.14))
    fig.subplots_adjust(right=0.82)

    # ---- interactive ----
    pfig = go.Figure()
    for r in refs:
        vals = [piv.loc[s, r] for s in stacks]
        text = ["★" if best[s] == r else "" for s in stacks]
        pfig.add_trace(go.Bar(
            x=[H.display(s) for s in stacks], y=vals, name=cfg.REFOLDER_LABELS[r],
            marker_color=cfg.REFOLDER_COLORS[r], marker_line_color=cfg.MARK_EDGE,
            marker_line_width=0.5, text=text, textposition="outside", textfont_size=15,
            hovertemplate=f"{cfg.REFOLDER_LABELS[r]}<br>%{{x}}: mean ipTM %{{y:.3f}}<extra></extra>"))
    pfig.update_layout(title=title, yaxis_title="mean ipTM", barmode="group",
                       legend_title="Refolder")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
