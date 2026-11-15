#!/usr/bin/env python3
"""plot_12 — protein amino-acid x RNA nucleotide interface-contact heatmap, per stack.

y = amino-acid type, x = RNA nt position (16-nt numbering), cell = mean heavy-atom
contacts (<=4.5 A) per design at that (aa, nt), from the AF3 refold. Pending until the
aa×nt map (interface_aa_nt_map.py) exists.
"""
from __future__ import annotations
from pathlib import Path

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_12_interface_heatmap"


def generate(df, cfg):
    mp = Path(cfg.INTERFACE_AA_NT_MAP)
    if not mp.exists():
        return H.pending_pair("Interface-contact map (protein aa x RNA nt)",
                              "aa x nt map not computed yet (interface_aa_nt_map.py)")
    n_designs = df.drop_duplicates("design_id").groupby("stack").size().to_dict()
    return H.aa_nt_heatmap(mp, "n_contacts", n_designs,
                           "Protein aa x RNA nt — interface-contact map",
                           "AF3 refold · cell = mean heavy-atom contacts (<=4.5 A) per design")


if __name__ == "__main__":
    _df = D.load()
    _p, _h = H.save_pair(*generate(_df, cfg), NAME)
    print(f"wrote {_p}\n      {_h}")
