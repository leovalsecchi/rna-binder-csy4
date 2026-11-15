#!/usr/bin/env python3
"""plot_01_2 — best-N designs PER stack.

Same consensus ipTM = min(AF3, B2) ranking as plot_01_1, but the top N of each stack
(N = cfg.TOP_N_PER_STACK), grouped in fixed stack order — so every generator's best
binders are visible even if one stack dominates the global top-N. Feeds
best_binders/best<N>_per_stack/.
"""
from __future__ import annotations
import re
import numpy as np

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_01_2_best_per_stack"
DOCKQ_COL = "dockq_af3_msa_off"


def _fmt(v, p=3):
    return "—" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{p}f}"


def _short(did):
    """Drop the redundant 'csy4free_<stack>_len80_' prefix -> 'uuid_designN_seqN'."""
    d = str(did)
    m = re.sub(r"^csy4free_.*?_len\d+_", "", d)
    return m if m != d else d.replace("csy4free_", "")


def generate(df, cfg):
    r = D.ranked_wide(stacks=list(df["stack"].unique()))
    if r.empty:
        return H.pending_pair("Best designs per stack",
                              "no designs with BOTH an AF3-msaoff and a Boltz-2 refold")
    r = r[r["rank_in_stack"] <= cfg.TOP_N_PER_STACK].copy()
    order = {s: i for i, s in enumerate(cfg.STACK_ORDER)}
    r["_so"] = r["stack"].map(order)
    r = r.sort_values(["_so", "rank_in_stack"], kind="stable")
    header = ["Stack", "In-stack #", "Design", "consensus", "AF3", "B2", "ESM", "DockQ", "global #"]
    rows, row_stacks, row_refs = [], [], []
    for _, x in r.iterrows():
        rows.append([
            H.short_stack(x["stack"]), str(int(x["rank_in_stack"])), _short(x["design_id"]),
            _fmt(x.get("consensus_iptm")), _fmt(x.get("af3_msaoff_iptm")),
            _fmt(x.get("b2_iptm")), _fmt(x.get("esm_iptm")),
            _fmt(x.get(DOCKQ_COL)) if DOCKQ_COL in r.columns else "—",
            str(int(x["rank"])),
        ])
        row_stacks.append(x["stack"])
        row_refs.append([x["design_id"], x.get("af3_msaoff_cif_path", ""),
                         x.get("b2_cif_path", ""), x.get("esm_cif_path", "")])
    title = f"Best {cfg.TOP_N_PER_STACK} designs per stack — by consensus ipTM = min(AF3, B2, ESM)"
    sub = "click a row → its 3 refolds open in ChimeraX  ·  extracted to best_binders/best%d_per_stack/" % cfg.TOP_N_PER_STACK
    fig_mpl, fig_ply = H.designs_table(header, rows, row_stacks, title, sub)
    return fig_mpl, fig_ply, H.clickable_table_html(header, rows, row_refs, title, sub)


if __name__ == "__main__":
    _df = D.load()
    res = generate(_df, cfg)
    _png, _html = H.save_pair(res[0], res[1], NAME, html_override=res[2] if len(res) > 2 else None)
    print(f"wrote {_png}\n      {_html}")
