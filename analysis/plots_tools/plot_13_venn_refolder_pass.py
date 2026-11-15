#!/usr/bin/env python3
"""plot_13 — Venn: designs clearing a threshold, per refolder.

One Venn panel per metric (ipTM, DockQ). Each circle is a refolder (AF3,
Boltz-2, ESMFold2) and holds the designs whose refold clears that metric's
`config` threshold; the overlaps are the designs that pass in several refolders
at once — i.e. refolder *agreement*. The centre of a 3-circle panel is the
designs that pass in all three.

Universe per panel = designs scored by *every* refolder shown, so set membership
is well-defined (a design missing an ESM refold is not silently counted as an ESM
failure — it is excluded and reported as `n`).

DockQ has no ESM sidecar column yet (only af3/b2), so its panel is a 2-circle
Venn. If an ESM DockQ column ever appears it auto-grows to 3 with no code change.
Thresholds/colors come only from config.py.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_13_venn_refolder_pass"
CLICK_JS = H.VENN_CLICK_JS   # click a region → open a RANDOM design from it in ChimeraX
SAMPLE_N = 30                # designs embedded per region (JS random-picks one per click)

# metrics to render, in order. `source` = where the per-refolder value comes from:
#   "refolder_iptm" -> D.refolder_iptm() (af3/b2/esm, MASTER_CSV)
#   "long"          -> the passed long df, column m["key"], keyed by refold_condition
METRICS = [
    {"key": "iptm",  "label": "ipTM",  "thr_key": "iptm_min",  "source": "refolder_iptm"},
    {"key": "dockq", "label": "DockQ", "thr_key": "dockq_min", "source": "dockq_sidecar"},
]

# DockQ sidecar columns per refolder (dockq_refolds.py). ESM column lands once its
# DockQ is computed → the DockQ panel auto-grows from 2 to 3 circles.
DOCKQ_SIDECAR_COLS = {"af3": "dockq_af3_msa_off", "b2": "dockq_b2_msa_off", "esm": "dockq_esm_msa_off"}

# geometry per set-count: circle centres (in refolder order) + label anchors.
_GEO = {
    3: {
        "R": 0.72,
        "centres": [(0.0, 0.52), (-0.52, -0.36), (0.52, -0.36)],  # top, bl, br
        "name_xy": [(0.0, 1.40), (-1.00, -1.14), (1.00, -1.14)],
        "name_va": ["bottom", "top", "top"],
        "regions": {                       # keyed by tuple of circle indices
            (0,): (0.0, 0.86), (1,): (-0.86, -0.56), (2,): (0.86, -0.56),
            (0, 1): (-0.48, 0.12), (0, 2): (0.48, 0.12), (1, 2): (0.0, -0.60),
            (0, 1, 2): (0.0, -0.04),
        },
        "xlim": (-1.55, 1.55), "ylim": (-1.65, 1.7),
    },
    2: {
        "R": 0.70,
        "centres": [(-0.37, 0.0), (0.37, 0.0)],
        "name_xy": [(-0.58, 0.92), (0.58, 0.92)],
        "name_va": ["bottom", "bottom"],
        "regions": {(0,): (-0.68, 0.0), (1,): (0.68, 0.0), (0, 1): (0.0, 0.0)},
        "xlim": (-1.35, 1.35), "ylim": (-1.1, 1.25),
    },
}


def _metric_wide(m, df):
    """dict refolder_key -> Series(design_id -> value) for refolders that have data,
    in cfg.REFOLDER_ORDER order."""
    if m["source"] == "dockq_sidecar":
        import pandas as pd
        p = Path(cfg.DOCKQ_SIDECAR)
        if not p.exists():
            return {}
        sc = pd.read_csv(p, low_memory=False)
        valid = set(df["design_id"])            # honor any --stacks / --master filtering
        out = {}
        for r in cfg.REFOLDER_ORDER:
            col = DOCKQ_SIDECAR_COLS.get(r)
            if col and col in sc.columns:
                s = sc[["design_id", col]].dropna()
                s = s[s["design_id"].isin(valid)].drop_duplicates("design_id").set_index("design_id")[col]
                if len(s):
                    out[r] = s
        return out
    # per-refolder iptm from MASTER_CSV (af3/b2/esm), same source as plot 7
    d = D.refolder_iptm(stacks=list(df["stack"].unique()))
    col = "iptm"
    out = {}
    for r in cfg.REFOLDER_ORDER:
        s = d[(d["refolder"] == r) & d[col].notna()]
        if len(s):
            out[r] = s.drop_duplicates("design_id").set_index("design_id")[col]
    return out


def _panel(m, df):
    """Compute the Venn payload for one metric, or None if <2 refolders have data."""
    wide = _metric_wide(m, df)
    refs = [r for r in cfg.REFOLDER_ORDER if r in wide]
    if len(refs) < 2:
        return None
    thr = cfg.THRESHOLDS.get(m["thr_key"])
    if thr is None:
        return None
    common = set.intersection(*(set(wide[r].index) for r in refs))
    passsets = {r: set(wide[r].index[wide[r] >= thr]) & common for r in refs}
    return {"m": m, "refs": refs, "thr": thr, "N": len(common), "passsets": passsets}


def _regions(refs, passsets):
    """{tuple-of-indices: count} for the 3 (2-set) or 7 (3-set) Venn regions, + union size."""
    S = [passsets[r] for r in refs]
    union = len(set().union(*S)) if S else 0
    if len(refs) == 3:
        A, B, C = S
        out = {(0,): len(A - B - C), (1,): len(B - A - C), (2,): len(C - A - B),
               (0, 1): len((A & B) - C), (0, 2): len((A & C) - B), (1, 2): len((B & C) - A),
               (0, 1, 2): len(A & B & C)}
    else:
        A, B = S
        out = {(0,): len(A - B), (1,): len(B - A), (0, 1): len(A & B)}
    return out, union


def _region_sets(refs, passsets):
    """{tuple-of-indices: set of design_ids} for each Venn region (parallels _regions)."""
    S = [passsets[r] for r in refs]
    if len(refs) == 3:
        A, B, C = S
        return {(0,): A - B - C, (1,): B - A - C, (2,): C - A - B,
                (0, 1): (A & B) - C, (0, 2): (A & C) - B, (1, 2): (B & C) - A,
                (0, 1, 2): A & B & C}
    A, B = S
    return {(0,): A - B, (1,): B - A, (0, 1): A & B}


def _region_label(idx, refs):
    names = [cfg.REFOLDER_LABELS[refs[i]] for i in idx]
    return (" ∩ ".join(names)) if len(idx) > 1 else (names[0] + " only")


def _region_sample(ids, refs_lut):
    """Deterministic even spread of <=SAMPLE_N designs (as [id,af3,b2,esm]) from a region."""
    ids = sorted(d for d in ids if d in refs_lut)
    if len(ids) > SAMPLE_N:
        step = len(ids) / SAMPLE_N
        ids = [ids[int(j * step)] for j in range(SAMPLE_N)]
    return [refs_lut[d] for d in ids]


def _draw_mpl(ax, panel):
    refs, k = panel["refs"], len(panel["refs"])
    g = _GEO[k]
    regions, union = _regions(refs, panel["passsets"])
    for (cx, cy), r in zip(g["centres"], refs):
        ax.add_patch(Circle((cx, cy), g["R"], facecolor=cfg.REFOLDER_COLORS[r], alpha=0.42,
                            edgecolor=cfg.MARK_EDGE, lw=1.3, zorder=2))
    for idx, (rx, ry) in g["regions"].items():
        ax.text(rx, ry, f"{regions.get(idx, 0):,}", ha="center", va="center",
                fontsize=15 if len(idx) == k else 13, fontweight="bold",
                color=cfg.INK_PRIMARY, zorder=5)
    for (nx, ny), va, r in zip(g["name_xy"], g["name_va"], refs):
        ax.text(nx, ny, cfg.REFOLDER_LABELS[r], ha="center", va=va, fontsize=10,
                fontweight="bold", color=cfg.REFOLDER_COLORS[r], zorder=5)
    ax.set_xlim(*g["xlim"]); ax.set_ylim(*g["ylim"])
    ax.set_aspect("equal"); ax.axis("off")
    scored_by = "both refolders" if k == 2 else f"all {k} refolders"
    ax.set_title(f"{panel['m']['label']} ≥ {panel['thr']:g}\n"
                 f"n = {panel['N']:,} designs scored by {scored_by}",
                 fontsize=12, fontweight="bold")
    ax.text(0.5, 0.005, f"{union:,} pass in ≥1 refolder", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9, color=cfg.INK_MUTED)


def _draw_plotly(panels, refs_lut):
    n = len(panels)
    titles = [f"{p['m']['label']} ≥ {p['thr']:g}  ·  n = {p['N']:,}" for p in panels]
    pfig = make_subplots(rows=1, cols=n, subplot_titles=titles, horizontal_spacing=0.06)
    for i, panel in enumerate(panels, 1):
        refs, k = panel["refs"], len(panel["refs"])
        g = _GEO[k]
        regions, union = _regions(refs, panel["passsets"])
        for (cx, cy), r in zip(g["centres"], refs):
            pfig.add_shape(type="circle", x0=cx - g["R"], y0=cy - g["R"],
                           x1=cx + g["R"], y1=cy + g["R"], layer="below",
                           fillcolor=H.rgba(cfg.REFOLDER_COLORS[r], 0.38),
                           line=dict(color=cfg.MARK_EDGE, width=1.3), row=1, col=i)
        # region counts as CLICKABLE points (near-invisible marker under the number); a click
        # opens a random design from that region. customdata = [label, JSON sample of refs].
        rsets = _region_sets(refs, panel["passsets"])
        rxs, rys, txts, cds, hts = [], [], [], [], []
        for idx, (rx, ry) in g["regions"].items():
            cnt = regions.get(idx, 0)
            rxs.append(rx); rys.append(ry); txts.append(f"<b>{cnt:,}</b>")
            cds.append([_region_label(idx, refs),
                        json.dumps(_region_sample(rsets.get(idx, set()), refs_lut))])
            hts.append(f"{_region_label(idx, refs)}: {cnt:,}<br>click → open a random design")
        pfig.add_trace(go.Scatter(
            x=rxs, y=rys, mode="markers+text", text=txts, textposition="middle center",
            textfont=dict(size=15, color=cfg.INK_PRIMARY),
            marker=dict(size=36, color="rgba(0,0,0,0.001)"),
            customdata=cds, hovertext=hts, hoverinfo="text", showlegend=False), row=1, col=i)
        for (nx, ny), r in zip(g["name_xy"], refs):
            pfig.add_annotation(x=nx, y=ny, text=f"<b>{cfg.REFOLDER_LABELS[r]}</b>",
                                showarrow=False, row=1, col=i,
                                font=dict(size=12, color=cfg.REFOLDER_COLORS[r]))
        pfig.add_annotation(x=sum(g["xlim"]) / 2, y=g["ylim"][0], yanchor="bottom",
                            text=f"{union:,} pass in ≥1 refolder", showarrow=False,
                            row=1, col=i, font=dict(size=10, color=cfg.INK_MUTED))
        xref = "x" if i == 1 else f"x{i}"
        pfig.update_xaxes(visible=False, range=list(g["xlim"]), row=1, col=i)
        pfig.update_yaxes(visible=False, range=list(g["ylim"]), scaleanchor=xref,
                          scaleratio=1, row=1, col=i)
    pfig.update_layout(showlegend=False,
                       title="Designs clearing each metric, by refolder  "
                             "(overlap = refolder agreement)")
    return pfig


def generate(df, cfg):
    panels = [p for p in (_panel(m, df) for m in METRICS) if p]
    if not panels:
        return H.pending_pair("Refolder pass Venn",
                              "no metric has ≥2 refolders scored yet")
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(5.7 * n, 6.0))
    axes = np.atleast_1d(axes)
    for ax, panel in zip(axes, panels):
        _draw_mpl(ax, panel)
    fig.suptitle("Designs clearing each metric, by refolder  (overlap = refolder agreement)",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.86, bottom=0.05, wspace=0.08)
    cp = D.cif_paths(stacks=list(df["stack"].unique()))
    refs_lut = {r["design_id"]: [r["design_id"], r.get("af3_msaoff_cif_path", ""),
                                 r.get("b2_cif_path", ""), r.get("esm_cif_path", "")]
                for _, r in cp.iterrows()}
    return fig, _draw_plotly(panels, refs_lut)


if __name__ == "__main__":
    _df = D.load()
    _png, _html = H.save_pair(*generate(_df, cfg), NAME, click_js=CLICK_JS)
    print(f"wrote {_png}\n      {_html}")
