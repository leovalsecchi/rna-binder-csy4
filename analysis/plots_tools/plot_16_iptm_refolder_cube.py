#!/usr/bin/env python3
"""plot_16 — ipTM refolder-agreement cube: AF3 ipTM × Boltz-2 ipTM × ESMFold2 ipTM.

One point per design (each axis is that refolder's ipTM). The green box in the upper corner
is where ALL THREE refolders clear iptm_min — i.e. full cross-refolder agreement. Rotatable;
click a point → its 3 refolds open (RNA-aligned) in ChimeraX.
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_16_iptm_refolder_cube"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.refolder_triplet("iptm", stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("ipTM refolder cube", "needs AF3 + B2 + ESM ipTM")
    return H.refolder_cube(d, value_min=cfg.THRESHOLDS["iptm_min"], metric_label="ipTM")


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
