#!/usr/bin/env python3
"""plot_10 — filtering funnel, STRICT (every refolder).

Same filter order + layout as plot_09 (avg), but only designs refolded by ALL conditions
(AF3 + Boltz-2 + ESMFold2, MSA-off) are eligible, and 'passing' requires the design to
independently clear each threshold in EVERY one of them, not just on average. Filters on
ipTM + DockQ only — the same criterion as the quality cube (plot_15), so the two report
the same pass count.
"""
from __future__ import annotations
import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_10_funnel_strict"


def generate(df, cfg):
    labels, counts, totals, filters = D.funnel(df, mode="strict")
    if len(labels) <= 1:
        return H.pending_pair("Filter funnel — strict (every refolder)",
                              "no threshold metrics available yet")
    pending = [m for (m, _thrs, _d) in cfg.FUNNEL_STAGES if not D.has_metric(df, m)]
    sub = "staggered cumulative · eligible = refolded by all conditions · passing = clears each threshold in every one"
    if pending:
        sub += "   ·   pending filters: " + ", ".join(pending)
    return H.funnel_figure(labels, counts, totals,
                           "Filter funnel — strict (every refolder, share per stack)", sub)


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
