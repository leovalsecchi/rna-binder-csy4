#!/usr/bin/env python3
"""plot_20 — pairwise refolder DockQ agreement.

Three joint scatterplots (AF3 vs Boltz-2, AF3 vs ESMFold2, Boltz-2 vs ESMFold2), each with
marginal DockQ distributions on both axes, an identity diagonal, threshold lines at dockq_min,
and Pearson r. One point per design (all three refolders present). Click a point → its 3
refolds open (RNA-aligned) in ChimeraX. Shows how much the refolders agree on DockQ.
Pending until an ESM DockQ column exists in the sidecar (needs all three refolders).
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_20_refolder_dockq_pairs"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.refolder_triplet("dockq", stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("Pairwise refolder DockQ", "needs AF3 + Boltz-2 + ESMFold2 DockQ")
    return H.refolder_pairgrid(d, metric_label="DockQ", value_min=cfg.THRESHOLDS["dockq_min"])


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
