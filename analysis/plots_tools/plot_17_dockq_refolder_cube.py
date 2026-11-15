#!/usr/bin/env python3
"""plot_17 — DockQ refolder-agreement cube: AF3 DockQ × Boltz-2 DockQ × ESMFold2 DockQ.

One point per design (each axis is that refolder's DockQ). The green box in the upper corner
is where ALL THREE refolders clear dockq_min. Rotatable; click a point → its 3 refolds open
(RNA-aligned) in ChimeraX. Needs the ESM DockQ sidecar column (else pending).
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_17_dockq_refolder_cube"
CLICK_JS = H.CLICK_JS


def generate(df, cfg):
    d = D.refolder_triplet("dockq", stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("DockQ refolder cube",
                              "needs AF3 + B2 + ESM DockQ (run ESM DockQ sidecar)")
    return H.refolder_cube(d, value_min=cfg.THRESHOLDS["dockq_min"], metric_label="DockQ")


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
