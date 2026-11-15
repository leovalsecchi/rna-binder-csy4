#!/usr/bin/env python3
"""plot_18 — scaffold fit: does the binder dock the csy4 hairpin WITHIN the full 121-nt
construct? Two criteria on one plot:

    X = scaffold_interface_frac  (fraction of the binder's RNA interface ON the hairpin;
                                  →1 = interface on-target, <1 = spilling onto the rest)
    Y = scaffold_clashes         (binder heavy atoms overlapping the REST of the construct;
                                  0 = fits cleanly, high = sterically blocked by scaffold)

The green box (bottom-right) = passes BOTH: interface on the hairpin AND ~no clashes. One point
per design, colored by stack; click → its 3 refolds open (RNA-aligned) in ChimeraX. The scaffold
was folded once (Boltz-2 winner); this metric is pure geometry (no per-design folding).
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import plotly.graph_objects as go

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_18_scaffold_fit"
CLICK_JS = H.CLICK_JS


def _fit_pct(d):
    """per-stack '(n=…, fit …%)' label, fit = passes both criteria."""
    cmax, imin = cfg.THRESHOLDS["scaffold_clash_max"], cfg.THRESHOLDS["scaffold_interface_min"]
    out = {}
    for s, g in d.groupby("stack"):
        ok = ((g["scaffold_clashes"] <= cmax) & (g["scaffold_interface_frac"] >= imin)).sum()
        out[s] = (len(g), 100 * ok / len(g) if len(g) else 0.0)
    return out


def generate(df, cfg):
    d = D.scaffold_fit(stacks=list(df["stack"].unique()))
    if d.empty:
        return H.pending_pair("Scaffold fit (interface on hairpin × clashes)",
                              "needs the scaffold_fit sidecar (run scaffold_fit.py)")
    stacks = H.ordered_stacks(d)
    cmax, imin = cfg.THRESHOLDS["scaffold_clash_max"], cfg.THRESHOLDS["scaffold_interface_min"]
    n_pass = int(((d["scaffold_clashes"] <= cmax) & (d["scaffold_interface_frac"] >= imin)).sum())
    lab = _fit_pct(d)
    # clashes span 0..hundreds; clip the display so the pass sliver is legible (outliers pinned to cap)
    ycap = float(max(cmax + 5, np.percentile(d["scaffold_clashes"], 97)))
    d = d.copy()
    d["cy"] = d["scaffold_clashes"].clip(upper=ycap)

    # ---- static ----
    fig, ax = plt.subplots(figsize=(8.8, 7))
    ax.add_patch(Rectangle((imin, -0.5), 1.02 - imin, cmax + 0.5, facecolor="#2ca02c",
                           alpha=0.12, edgecolor="#2ca02c", linewidth=0.8, zorder=0))
    for s in stacks:
        g = d[d["stack"] == s]
        n, pct = lab[s]
        ax.scatter(g["scaffold_interface_frac"], g["cy"], s=12, alpha=0.5, color=cfg.color(s),
                   edgecolor="none", label=f"{H.display(s)} (n={n}, fit {pct:.0f}%)", zorder=3)
    ax.axhline(cmax, color="#2ca02c", lw=0.7, ls="--", zorder=1)
    ax.set_xlabel("interface fraction on the hairpin  →  on-target")
    ax.set_ylabel(f"clashes with the rest of the construct  (clipped at {ycap:.0f})")
    ax.set_ylim(-0.5, ycap + 1)
    ax.set_title("Scaffold fit — binder vs the full 121-nt construct\n"
                 f"(green box = interface ≥{imin:g} & clashes ≤{cmax}: {n_pass}/{len(d)} designs)",
                 fontsize=11)
    ax.legend(fontsize=7.5, framealpha=0.7, loc="upper left")

    # ---- interactive ----
    pfig = go.Figure()
    pfig.add_shape(type="rect", x0=imin, x1=1.02, y0=-0.5, y1=cmax, layer="below",
                   fillcolor="rgba(44,160,44,0.12)", line=dict(color="rgba(44,160,44,0.5)", width=1))
    pfig.add_annotation(x=1.02, y=-0.5, text=f"fits: interface≥{imin:g} & clashes≤{cmax} (n={n_pass})",
                        showarrow=False, xanchor="right", yanchor="bottom",
                        font=dict(color="#2ca02c", size=11))
    for s in stacks:
        g = d[d["stack"] == s]
        n, pct = lab[s]
        pfig.add_trace(go.Scattergl(
            x=g["scaffold_interface_frac"], y=g["cy"], mode="markers",
            name=f"{H.display(s)} · fit {pct:.0f}%",
            marker=dict(color=cfg.color(s), size=6, opacity=0.6,
                        line=dict(color=cfg.MARK_EDGE, width=0.4)),
            customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path",
                          "scaffold_clashes"]].fillna(""),
            hovertemplate=(f"{H.display(s)}<br>interface on hairpin=%{{x:.2f}}"
                           "<br>clashes=%{customdata[4]}<br>%{customdata[0]}"
                           "<br><i>click → open 3 refolds (RNA-aligned)</i><extra></extra>")))
    pfig.update_layout(
        title=f"Scaffold fit — binder vs full 121-nt construct  ({n_pass}/{len(d)} fit · click → 3 refolds)",
        xaxis_title="interface fraction on the hairpin  →  on-target",
        yaxis_title=f"clashes with the rest of the construct (clipped at {ycap:.0f})")
    return fig, pfig


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
