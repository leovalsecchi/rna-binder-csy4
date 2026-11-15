#!/usr/bin/env python3
"""plot_05 — RNA nucleotide interface frequency.

x = RNA nt position (16-nt numbering), y = % of designs where that nt appears in
interface_contacts, grouped bars per stack. Denominator = ALL designs per stack
(design-level metric, deduped to one row per design).
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_05_nt_interface_freq"
CONTACT_COL, LABEL = "interface_contacts", "in interface"


def _freq_by_stack(df):
    """Return {stack: np.array(len=RNA_NT_COUNT) of % designs with contact at nt}."""
    nts = np.arange(1, cfg.RNA_NT_COUNT + 1)
    out = {}
    uniq = df.drop_duplicates("design_id")  # design-level -> one row per design
    for s in H.ordered_stacks(uniq):
        sub = uniq[uniq["stack"] == s]
        n = len(sub)
        if not n:
            continue
        counts = np.zeros(cfg.RNA_NT_COUNT)
        for contacts in sub[CONTACT_COL]:
            for nt in contacts:
                if 1 <= nt <= cfg.RNA_NT_COUNT:
                    counts[nt - 1] += 1
        out[s] = 100.0 * counts / n
    return nts, out


def generate(df, cfg):
    if not (df[CONTACT_COL].map(len).sum()):
        return H.pending_pair("Nucleotide interface frequency", f"{CONTACT_COL} not available")
    nts, freq = _freq_by_stack(df)
    stacks = list(freq)
    title = "RNA nucleotide interface frequency  (% of all designs, per stack)"

    # ---- static ----
    fig, ax = plt.subplots(figsize=(12, 5.4))
    w = 0.8 / len(stacks)
    for i, s in enumerate(stacks):
        ax.bar(nts + (i - (len(stacks) - 1) / 2) * w, freq[s], width=w, color=cfg.color(s),
               edgecolor=cfg.MARK_EDGE, linewidth=0.4, label=H.display(s))
    ax.set_xticks(nts)
    ax.set_xlabel("RNA nucleotide position"); ax.set_ylabel("% of designs")
    ax.set_title(title); ax.grid(axis="x", visible=False)
    ax.legend(fontsize=8, framealpha=0.7, ncol=2)

    # ---- interactive ----
    pfig = go.Figure()
    for s in stacks:
        pfig.add_trace(go.Bar(x=nts, y=freq[s], name=H.display(s), marker_color=cfg.color(s),
                              marker_line_color=cfg.MARK_EDGE, marker_line_width=0.4,
                              hovertemplate=f"{H.display(s)}<br>nt %{{x}}: %{{y:.1f}}%<extra></extra>"))
    pfig.update_layout(title=title, xaxis_title="RNA nucleotide position", yaxis_title="% of designs",
                       barmode="group", xaxis=dict(tickmode="linear"))
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
