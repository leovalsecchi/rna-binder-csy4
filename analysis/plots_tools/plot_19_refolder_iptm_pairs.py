#!/usr/bin/env python3
"""plot_19 — pairwise refolder ipTM agreement.

Three joint scatterplots (AF3 vs Boltz-2, AF3 vs ESMFold2, Boltz-2 vs ESMFold2), each with
marginal ipTM distributions on both axes, an identity diagonal, threshold lines at iptm_min,
and Pearson r. One point per design (all three refolders present). Click a point → its 3
refolds open (RNA-aligned) in ChimeraX. Shows how much the refolders agree on ipTM.
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_19_refolder_iptm_pairs"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.refolder_triplet("iptm", stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("Pairwise refolder ipTM", "needs AF3 + Boltz-2 + ESMFold2 ipTM")
    return H.refolder_pairgrid(d, metric_label="ipTM", value_min=cfg.THRESHOLDS["iptm_min"])


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
