#!/usr/bin/env python3
"""plot_22 — consensus ipTM × DockQ scatter with a clickable table of the ELITE designs below.

Axes = one point/design at consensus values (ipTM = min(AF3,B2,ESM), DockQ = min across refolders).
"Passing" (big dots + table) uses the SAME elite mechanism as plot_21 (config.TOPX_FUNNEL): a design
passes iff it is in EACH refolder's global top-`top_x` by that refolder's OWN ipTM (triple
intersection) AND its consensus DockQ ≥ dockq_min. BOTH the scatter points and the table rows are
click-to-open (design's 3 refolds → ChimeraX). Static PNG = scatter over a rendered table; interactive
HTML = plotly scatter + clickable HTML table.
"""
from __future__ import annotations
import json
import re
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import plotly.graph_objects as go
import plotly.io as pio

import config as cfg
import data_loader as D
import plot_helpers as H

NAME = "plot_22_scatter_table"
CLICK_JS = H.CLICK_JS
MAX_TABLE = 40          # cap table rows (sorted by consensus ipTM) when many pass


def _short(did):
    d = str(did); m = re.sub(r"^csy4free_.*?_len\d+_", "", d)
    return m if m != d else d.replace("csy4free_", "")


def generate(df, cfg):
    top_x = int(cfg.TOPX_FUNNEL["top_x"]); dmin = float(cfg.TOPX_FUNNEL["dockq_min"])
    title = f"Consensus ipTM × DockQ — scatter + elite table (triple top-{top_x} ipTM)"

    ti = D.refolder_triplet("iptm")
    td = D.refolder_triplet("dockq")
    if ti.empty or td.empty:
        return H.pending_pair(title, "needs ipTM AND DockQ for all three refolders")
    d = (ti.rename(columns={"af3": "ip_af3", "b2": "ip_b2", "esm": "ip_esm"})
           [["design_id", "stack", "ip_af3", "ip_b2", "ip_esm",
             "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]]
         .merge(td.rename(columns={"af3": "dq_af3", "b2": "dq_b2", "esm": "dq_esm"})
                  [["design_id", "dq_af3", "dq_b2", "dq_esm"]], on="design_id", how="inner"))
    if d.empty:
        return H.pending_pair(title, "no designs carry all three ipTM AND all three DockQ")
    d["iptm"] = d[["ip_af3", "ip_b2", "ip_esm"]].min(axis=1)     # x-axis (consensus)
    d["dockq"] = d[["dq_af3", "dq_b2", "dq_esm"]].min(axis=1)    # y-axis (consensus)
    # elite ipTM criterion = in EACH refolder's global top-X by its OWN ipTM
    tri = (set(d.nlargest(top_x, "ip_af3")["design_id"])
           & set(d.nlargest(top_x, "ip_b2")["design_id"])
           & set(d.nlargest(top_x, "ip_esm")["design_id"]))
    d["in_triple"] = d["design_id"].isin(tri)
    d["pass"] = d["in_triple"] & (d["dockq"] >= dmin)

    stacks = H.ordered_stacks(d)
    passing = d[d["pass"]].sort_values("iptm", ascending=False)
    n_pass = len(passing)
    sub = (f"axes = consensus min(AF3,B2,ESM)  ·  big dots + table = ELITE: in all-3 top-{top_x} "
           f"ipTM AND consensus DockQ ≥ {dmin:g}  ·  {len(tri)} in triple top-{top_x} → {n_pass} also pass DockQ")

    # ---------- table rows (elite designs) ----------
    header = ["#", "Stack", "Design", "cons. ipTM", "cons. DockQ"]
    rows, row_refs = [], []
    for i, (_, x) in enumerate(passing.head(MAX_TABLE).iterrows(), 1):
        rows.append([str(i), H.short_stack(x["stack"]), _short(x["design_id"]),
                     f"{x['iptm']:.3f}", f"{x['dockq']:.3f}"])
        row_refs.append([x["design_id"], x.get("af3_msaoff_cif_path", ""),
                         x.get("b2_cif_path", ""), x.get("esm_cif_path", "")])
    ntab = len(rows)

    # ---------- static (matplotlib): scatter over table ----------
    fig = plt.figure(figsize=(9.5, 7.2 + 0.30 * max(ntab, 1)))
    gs = GridSpec(2, 1, height_ratios=[7.0, max(1.3, 0.32 * max(ntab, 1))], hspace=0.30)
    ax = fig.add_subplot(gs[0])
    ax.axhspan(dmin, 1.0, color="#f4f4f0", zorder=0)     # DockQ hard-filter band
    for s in stacks:
        for passed, alpha, sz in [(False, 0.16, 11), (True, 0.95, 46)]:
            g = d[(d["stack"] == s) & (d["pass"] == passed)]
            ax.scatter(g["iptm"], g["dockq"], s=sz, alpha=alpha, color=cfg.color(s),
                       edgecolor=cfg.MARK_EDGE if passed else "none", linewidth=0.6,
                       zorder=3 if passed else 2)
    H.add_mpl_threshold(ax, dmin, "y", f"DockQ {dmin:g}")
    ax.set_xlabel("consensus ipTM  (min AF3 / Boltz-2 / ESMFold2)")
    ax.set_ylabel("consensus DockQ")
    ax.set_title(f"Consensus ipTM × DockQ   (big dots = triple top-{top_x} ipTM & DockQ ≥ {dmin:g}, n={n_pass})")
    handles = [plt.Line2D([], [], marker="o", ls="", color=cfg.color(s), label=H.display(s))
               for s in stacks]
    ax.legend(handles=handles, fontsize=8, framealpha=0.7, loc="lower left")

    axt = fig.add_subplot(gs[1]); axt.set_axis_off()
    if rows:
        tbl = axt.table(cellText=rows, colLabels=header, loc="center", cellLoc="center",
                        colWidths=[0.07, 0.20, 0.48, 0.125, 0.125])
        tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.35)
        for (rr, cc), cell in tbl.get_celld().items():
            cell.set_edgecolor("#d9d8d2")
            if rr == 0:
                cell.set_facecolor("#f0efec"); cell.set_text_props(weight="bold")
            elif cc == 1 and (rr - 1) < ntab:
                cell.set_facecolor(cfg.color(passing.iloc[rr - 1]["stack"]) + "22")
                cell.set_text_props(weight="bold")
            if cc == 2:
                cell.set_text_props(ha="left")
        axt.set_title(f"Elite designs (triple top-{top_x} ipTM ∩ DockQ ≥ {dmin:g})"
                      + (f" — top {MAX_TABLE} of {n_pass}" if n_pass > MAX_TABLE else ""),
                      fontsize=9, pad=6)
    else:
        axt.text(0.5, 0.5, f"no designs in all-3 top-{top_x} ipTM with DockQ ≥ {dmin:g}",
                 ha="center", va="center", fontsize=11, color=cfg.INK_MUTED)

    # ---------- interactive scatter (plotly) ----------
    pfig = go.Figure()
    for s in stacks:
        g = d[d["stack"] == s]
        pfig.add_trace(go.Scattergl(
            x=g["iptm"], y=g["dockq"], mode="markers", name=H.display(s),
            marker=dict(color=cfg.color(s), size=[13 if p else 5 for p in g["pass"]],
                        opacity=0.82, line=dict(color=cfg.MARK_EDGE, width=0.5)),
            customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path",
                          "esm_cif_path"]].fillna(""),
            hovertemplate=(f"{H.display(s)}<br>consensus ipTM=%{{x:.3f}}  DockQ=%{{y:.3f}}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds</i><extra></extra>")))
    H.add_plotly_threshold(pfig, dmin, "y", f"DockQ {dmin:g}")
    pfig.update_layout(xaxis_title="consensus ipTM", yaxis_title="consensus DockQ",
                       height=560, margin=dict(t=12, b=44, l=64, r=20))

    # ---------- combined HTML: plotly scatter + clickable HTML table ----------
    scatter_div = pio.to_html(pfig, include_plotlyjs="cdn", full_html=False, div_id="scat")
    th = "".join(f"<th>{c}</th>" for c in header)
    trs = []
    for r, refs in zip(rows, row_refs):
        tds = "".join(f"<td>{c}</td>" for c in r)
        trs.append(f"<tr onclick='window.openRefolds({json.dumps(refs)})' "
                   f"title='click → open 3 refolds in ChimeraX'>{tds}</tr>")
    table_html = (f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(trs)}</tbody></table>"
                  if rows else f"<p><em>No designs in all-3 top-{top_x} ipTM with DockQ ≥ {dmin:g}.</em></p>")
    tcap = (f"Elite designs — {n_pass} (triple top-{top_x} ipTM ∩ DockQ ≥ {dmin:g})"
            + (f", showing top {MAX_TABLE}" if n_pass > MAX_TABLE else "") + " · click a row → 3 refolds")
    combined = f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<style>body{{font:14px/1.55 system-ui,-apple-system,Segoe UI,sans-serif;margin:22px;color:#0b0b0b}}
h2{{font-size:19px;margin:0 0 3px}} .sub{{color:#6b6a66;font-size:12.5px;margin:0 0 8px}}
h3{{font-size:15px;margin:20px 0 8px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{padding:6px 10px;border-bottom:1px solid #e1e0d9;text-align:center}}
td:nth-child(3){{text-align:left;font-family:ui-monospace,Menlo,monospace;font-size:12px}}
th{{background:#f0efec;font-weight:700}}
tr[onclick]{{cursor:pointer}} tr[onclick]:hover{{background:#eaf3ff}}
.hint{{color:#898781;font-size:12px;margin-top:10px}}</style></head><body>
<h2>{title}</h2><div class="sub">{sub}</div>
{scatter_div}
<h3>{tcap}</h3>
{table_html}
<p class="hint">Click a scatter point <b>or</b> a table row → its 3 refolds open (RNA-aligned) in
ChimeraX (needs <code>openrefolds</code> running on :45001).</p>
<script>{H.CLICK_JS}</script>
</body></html>"""
    return fig, pfig, combined


if __name__ == "__main__":
    _df = D.load()
    res = generate(_df, cfg)
    _png, _html = H.save_pair(res[0], res[1], NAME, html_override=res[2] if len(res) > 2 else None)
    print(f"wrote {_png}\n      {_html}")
