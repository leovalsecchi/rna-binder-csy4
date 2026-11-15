#!/usr/bin/env python3
"""plot_01 — per-stack overview table.

Number of designs generated per stack, and mean ipTM / pTM / DockQ / Cα-RMSD / PAE
(averaged across the in-scope MSA-off refolds). Metrics not yet computed show '—'.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_01_overview"
COLS = [("iptm", "ipTM"), ("ptm", "pTM"), ("dockq", "DockQ"),
        ("protein_rmsd", "Cα-RMSD"), ("pae", "PAE")]


def _table(df):
    stacks = H.ordered_stacks(df)
    header = ["Stack", "N designs"] + [lab for _, lab in COLS]
    rows = []
    for s in stacks:
        sub = df[df["stack"] == s]
        n = sub["design_id"].nunique()
        cells = [H.display(s), f"{n:,}"]
        for col, _ in COLS:
            m = sub[col].mean()
            cells.append("—" if (m is None or np.isnan(m)) else f"{m:.3f}")
        rows.append(cells)
    return header, rows, stacks


def generate(df, cfg):
    header, rows, stacks = _table(df)

    # ---- static ----
    fig, ax = plt.subplots(figsize=(11, 0.55 * len(rows) + 1.4))
    ax.set_axis_off()
    colw = [0.30, 0.11] + [0.118] * len(COLS)
    tbl = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center", colWidths=colw)
    tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d9d8d2")
        if r == 0:
            cell.set_facecolor("#f0efec"); cell.set_text_props(weight="bold", color=cfg.INK_PRIMARY)
        elif c == 0:  # stack name cell — tint with the stack color for identity
            cell.set_facecolor(cfg.color(stacks[r - 1]) + "22")
            cell.set_text_props(weight="bold", ha="left")
            cell.PAD = 0.03
    ax.set_title("Campaign overview — designs and mean metrics per stack", pad=16)

    # ---- interactive ----
    cols_t = list(zip(*rows)) if rows else [[] for _ in header]
    # per-column fill: stack-name column tinted by stack color, rest white
    fill = [[H.rgba(cfg.color(s), 0.13) for s in stacks]]
    for _ in range(len(header) - 1):
        fill.append(["#ffffff"] * len(rows))
    pfig = go.Figure(go.Table(
        header=dict(values=[f"<b>{h}</b>" for h in header], fill_color="#f0efec",
                    font=dict(color=cfg.INK_PRIMARY, size=13), align="center"),
        cells=dict(values=cols_t, fill_color=fill[:len(header)], align="center",
                   font=dict(size=12), height=28)))
    pfig.update_layout(title="Campaign overview — designs and mean metrics per stack")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
