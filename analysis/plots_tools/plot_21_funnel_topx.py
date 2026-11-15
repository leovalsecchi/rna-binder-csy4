#!/usr/bin/env python3
"""plot_21 — ELITE selection funnel: consensus DockQ hard-filter + triple top-X ipTM.

Per-stack stacked funnel (same style as 09/10). Stages, cumulative top→bottom:
  1. Total                    — designs refolded by all three (AF3, Boltz-2, ESMFold2)
  2. [optional] Consensus ipTM ≥ iptm_min   — min(AF3,B2,ESM); shown only if iptm_min > 0
  3. Consensus DockQ ≥ dockq_min            — min(AF3,B2,ESM); the HARD filter
  4. ∩ each refolder's top-X ipTM           — in AF3's top-X AND Boltz-2's top-X AND ESM's
                                              top-X (each ranked by its OWN ipTM, global)

All knobs in config.TOPX_FUNNEL (top_x, dockq_min, iptm_min). ipTM/DockQ thresholds also
live in config.THRESHOLDS. Runs on whatever MASTER_CSV / DOCKQ_SIDECAR config points at
(the 10k sample by default).
"""
from __future__ import annotations
import numpy as np

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_21_funnel_topx"


def generate(df, cfg):
    k = cfg.TOPX_FUNNEL
    top_x = int(k["top_x"]); dqmin = float(k["dockq_min"]); ipmin = float(k.get("iptm_min") or 0)

    ti = D.refolder_triplet("iptm")
    td = D.refolder_triplet("dockq")
    title = f"Elite funnel — consensus DockQ + triple top-{top_x} ipTM (per stack)"
    if ti.empty or td.empty:
        return H.pending_pair(title, "needs ipTM AND DockQ for all three refolders (af3/b2/esm)")

    m = (ti.rename(columns={"af3": "ip_af3", "b2": "ip_b2", "esm": "ip_esm"})
           [["design_id", "stack", "ip_af3", "ip_b2", "ip_esm"]]
         .merge(td.rename(columns={"af3": "dq_af3", "b2": "dq_b2", "esm": "dq_esm"})
                  [["design_id", "dq_af3", "dq_b2", "dq_esm"]], on="design_id", how="inner"))
    if m.empty:
        return H.pending_pair(title, "no designs carry all three ipTM AND all three DockQ")

    m["dq_min"] = m[["dq_af3", "dq_b2", "dq_esm"]].min(axis=1)
    m["ip_min"] = m[["ip_af3", "ip_b2", "ip_esm"]].min(axis=1)
    # global top-X per refolder by its OWN ipTM; a design must be in ALL THREE lists
    tri = (set(m.nlargest(top_x, "ip_af3")["design_id"])
           & set(m.nlargest(top_x, "ip_b2")["design_id"])
           & set(m.nlargest(top_x, "ip_esm")["design_id"]))
    m["in_triple"] = m["design_id"].isin(tri)

    stages = [("Total (3 refolds)", None)]
    if ipmin > 0:
        stages.append((f"Consensus ipTM ≥ {ipmin:g}", m["ip_min"] >= ipmin))
    stages.append((f"Consensus DockQ ≥ {dqmin:g}", m["dq_min"] >= dqmin))
    stages.append((f"∩ each refolder's top-{top_x} ipTM", m["in_triple"]))

    stacks = [s for s in H.ordered_stacks(m) if (m["stack"] == s).any()]
    labels, totals = [], []
    counts = {s: [] for s in stacks}
    cum = np.ones(len(m), bool)
    for lab, mask in stages:
        if mask is not None:
            cum = cum & mask.to_numpy()
        labels.append(lab)
        totals.append(int(cum.sum()))
        vc = m[cum]["stack"].value_counts()
        for s in stacks:
            counts[s].append(int(vc.get(s, 0)))

    sub = (f"{len(m):,} designs with all 3 refolds  ·  DockQ filter = min(AF3,B2,ESM) ≥ {dqmin:g}"
           f"  ·  top-{top_x} = global per-refolder ipTM, design in ALL THREE"
           + (f"  ·  ipTM pre-filter min ≥ {ipmin:g}" if ipmin > 0 else ""))
    return H.funnel_figure(labels, counts, totals, title, sub)


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_png}\n      {_html}")
