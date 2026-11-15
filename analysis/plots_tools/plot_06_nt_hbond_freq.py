#!/usr/bin/env python3
"""plot_06 — RNA nucleotide H-bond frequency.

Same as plot_05 but y = % of designs where that nt appears in hbond_contacts
(RNA nts specifically involved in an H-bond). Pending until a per-nt H-bond
sidecar exists (same hbond definition as tip_weighted_contact_score.py, type=='hbond').
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_06_nt_hbond_freq"
CONTACT_COL = "hbond_contacts"


def _freq_by_stack(df):
    uniq = df.drop_duplicates("design_id")
    nts = np.arange(1, cfg.RNA_NT_COUNT + 1)
    out = {}
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
        return H.pending_pair("Nucleotide H-bond frequency",
                              "per-nt H-bond contacts not computed yet (needs hbond_nts sidecar)")
    nts, freq = _freq_by_stack(df)
    stacks = list(freq)
    title = "RNA nucleotide H-bond frequency  (% of all designs, per stack)"

    fig, ax = plt.subplots(figsize=(12, 5.4))
    w = 0.8 / len(stacks)
    for i, s in enumerate(stacks):
        ax.bar(nts + (i - (len(stacks) - 1) / 2) * w, freq[s], width=w, color=cfg.color(s),
               edgecolor=cfg.MARK_EDGE, linewidth=0.4, label=H.display(s))
    seq = cfg.RNA_SEQ if len(cfg.RNA_SEQ) == cfg.RNA_NT_COUNT else None
    ax.set_xticks(nts)
    ax.set_xticklabels([f"{q}\n{seq[q - 1]}" if seq else str(q) for q in nts])
    ax.set_xlabel("RNA nucleotide position (5'→3')"); ax.set_ylabel("% of designs")
    ax.set_title(title); ax.grid(axis="x", visible=False)
    ax.legend(fontsize=8, framealpha=0.7, ncol=2)

    letters = [seq[q - 1] if seq else "" for q in nts]
    pfig = go.Figure()
    for s in stacks:
        pfig.add_trace(go.Bar(x=nts, y=freq[s], name=H.display(s), marker_color=cfg.color(s),
                              marker_line_color=cfg.MARK_EDGE, marker_line_width=0.4,
                              customdata=letters,
                              hovertemplate=f"{H.display(s)}<br>nt %{{x}} (%{{customdata}}): %{{y:.1f}}%<extra></extra>"))
    xt = [f"{q}<br>{seq[q - 1]}" if seq else str(q) for q in nts]
    pfig.update_layout(title=title, xaxis_title="RNA nucleotide position (5'→3')", yaxis_title="% of designs",
                       barmode="group", xaxis=dict(tickmode="array", tickvals=list(nts), ticktext=xt))
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
