#!/usr/bin/env python3
"""plot_01_1 — top-N designs across ALL stacks.

Global ranking table by consensus ipTM = min(AF3, B2) (a design ranks high only if BOTH
refolders agree; matches combine_and_rank.py). One row per design with its per-refolder
ipTMs and DockQ. N = cfg.TOP_N_ACROSS. Feeds best_binders/top<N>_across_all/.
"""
from __future__ import annotations
import re
import numpy as np

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_01_1_top_designs"
DOCKQ_COL = "dockq_af3_msa_off"


def _fmt(v, p=3):
    return "—" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{p}f}"


def _short(did):
    """Drop the redundant 'csy4free_<stack>_len80_' prefix (stack shown separately) ->
    keeps 'uuid_designN_seqN'."""
    d = str(did)
    m = re.sub(r"^csy4free_.*?_len\d+_", "", d)
    return m if m != d else d.replace("csy4free_", "")


def generate(df, cfg):
    r = D.ranked_wide(stacks=list(df["stack"].unique()))
    if r.empty:
        return H.pending_pair("Top designs across all stacks",
                              "no designs with BOTH an AF3-msaoff and a Boltz-2 refold")
    r = r.head(cfg.TOP_N_ACROSS)
    header = ["Rank", "Stack", "Design", "consensus", "AF3", "B2", "ESM", "DockQ"]
    rows, row_stacks, row_refs = [], [], []
    for _, x in r.iterrows():
        rows.append([
            str(int(x["rank"])), H.short_stack(x["stack"]), _short(x["design_id"]),
            _fmt(x.get("consensus_iptm")), _fmt(x.get("af3_msaoff_iptm")),
            _fmt(x.get("b2_iptm")), _fmt(x.get("esm_iptm")),
            _fmt(x.get(DOCKQ_COL)) if DOCKQ_COL in r.columns else "—",
        ])
        row_stacks.append(x["stack"])
        row_refs.append([x["design_id"], x.get("af3_msaoff_cif_path", ""),
                         x.get("b2_cif_path", ""), x.get("esm_cif_path", "")])
    title = f"Top {len(r)} designs across all stacks — by consensus ipTM = min(AF3, B2, ESM)"
    sub = "click a row → its 3 refolds open in ChimeraX  ·  extracted to best_binders/top%d_across_all/" % cfg.TOP_N_ACROSS
    fig_mpl, fig_ply = H.designs_table(header, rows, row_stacks, title, sub)
    return fig_mpl, fig_ply, H.clickable_table_html(header, rows, row_refs, title, sub)


if __name__ == "__main__":
    _df = D.load()
    res = generate(_df, cfg)
    _png, _html = H.save_pair(res[0], res[1], NAME, html_override=res[2] if len(res) > 2 else None)
    print(f"wrote {_png}\n      {_html}")
