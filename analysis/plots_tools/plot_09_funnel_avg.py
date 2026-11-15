#!/usr/bin/env python3
"""plot_09 — filtering funnel, AVERAGE across refolders.

Cumulative filters (iptm → dockq → disorder → clashes), % of designs still
passing after each stage, per stack. 'Passing' = the design's MEAN value across the
refold conditions clears the threshold. Filters whose metric isn't computed yet
are skipped and named in the subtitle (they slot in automatically once available).
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_09_funnel_avg"


def generate(df, cfg):
    labels, counts, totals, filters = D.funnel(df, mode="avg")
    if len(labels) <= 1:
        return H.pending_pair("Filter funnel — average across refolders",
                              "no threshold metrics available yet")
    pending = [m for (m, _thrs, _d) in cfg.FUNNEL_STAGES if not D.has_metric(df, m)]
    sub = "staggered cumulative · passing = design's MEAN across the two refold conditions clears each threshold"
    if pending:
        sub += "   ·   pending filters: " + ", ".join(pending)
    return H.funnel_figure(labels, counts, totals,
                           "Filter funnel — average across refolders (share per stack)", sub)


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
