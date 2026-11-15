"""plot_helpers.py — shared house style + save/degrade utilities so all 10 plots
read as one system. Static = matplotlib, interactive = plotly, both driven by
config.py (never hardcode color/threshold here either).
"""
from __future__ import annotations
import json
import warnings
from html import escape as _escape
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

import config as cfg

# ------------------------------------------------------------------ matplotlib style
plt.rcParams.update({
    "font.family": "DejaVu Sans",          # system-ui isn't available to mpl on Euler
    "font.size": 11,
    "figure.facecolor": cfg.SURFACE,
    "axes.facecolor": cfg.SURFACE,
    "axes.edgecolor": "#c3c2b7",
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": cfg.GRID,
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": cfg.INK_SECONDARY,
    "ytick.color": cfg.INK_SECONDARY,
    "text.color": cfg.INK_PRIMARY,
    "axes.labelcolor": cfg.INK_PRIMARY,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "figure.dpi": 110,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

# ------------------------------------------------------------------ plotly template
pio.templates["binder"] = go.layout.Template(layout=dict(
    font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=13, color=cfg.INK_PRIMARY),
    paper_bgcolor=cfg.SURFACE,
    plot_bgcolor=cfg.SURFACE,
    colorway=list(cfg.STACK_COLORS.values()),
    xaxis=dict(gridcolor=cfg.GRID, zerolinecolor=cfg.GRID, linecolor="#c3c2b7", ticks="outside",
               tickcolor="#c3c2b7", title_font_color=cfg.INK_PRIMARY),
    yaxis=dict(gridcolor=cfg.GRID, zerolinecolor=cfg.GRID, linecolor="#c3c2b7", ticks="outside",
               tickcolor="#c3c2b7", title_font_color=cfg.INK_PRIMARY),
    legend=dict(bgcolor="rgba(255,255,255,0.7)", bordercolor=cfg.GRID, borderwidth=1),
    title=dict(font=dict(size=17, color=cfg.INK_PRIMARY)),
    margin=dict(l=70, r=30, t=70, b=60),
))
pio.templates.default = "binder"


# ------------------------------------------------------------------ helpers
def display(stack: str) -> str:
    return cfg.STACK_DISPLAY.get(stack, stack)


_SHORT_STACK = {v: k for k, v in cfg.STACK_ALIASES.items()}


def short_stack(stack: str) -> str:
    """Compact stack label for dense tables: 'RFD3+LiMPNN' vs 'RFD3 + LigandMPNN'."""
    return _SHORT_STACK.get(stack, display(stack))


def ordered_stacks(df) -> list[str]:
    present = set(df["stack"].unique())
    return [s for s in cfg.STACK_ORDER if s in present]


def rgba(hexcol: str, alpha: float = 1.0) -> str:
    """'#RRGGBB' -> 'rgba(r,g,b,a)' (plotly rejects 8-digit hex; matplotlib accepts it)."""
    h = hexcol.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def threshold_note(metric_min_key: str) -> str:
    v = cfg.THRESHOLDS.get(metric_min_key)
    return "" if v is None else f"{metric_min_key} = {v}"


def add_mpl_threshold(ax, value, axis="y", label=None):
    if value is None:
        return
    fn = ax.axhline if axis == "y" else ax.axvline
    fn(value, color=cfg.THRESHOLD_LINE, ls="--", lw=1.4, zorder=1)
    if label:
        if axis == "y":
            ax.text(0.995, value, f" {label}", color=cfg.THRESHOLD_LINE, va="bottom",
                    ha="right", transform=ax.get_yaxis_transform(), fontsize=9)
        else:
            ax.text(value, 0.995, f" {label}", color=cfg.THRESHOLD_LINE, va="top",
                    ha="left", rotation=90, transform=ax.get_xaxis_transform(), fontsize=9)


def add_plotly_threshold(fig, value, axis="y", label=None):
    if value is None:
        return
    kw = dict(line=dict(color=cfg.THRESHOLD_LINE, dash="dash", width=1.4),
              annotation_text=label, annotation_position="top left" if axis == "y" else "top right")
    (fig.add_hline if axis == "y" else fig.add_vline)(value, **kw)


# ------------------------------------------------------------------ degradation
def pending_pair(title: str, message: str):
    """A clean 'metric pending' placeholder (static + interactive) for missing data."""
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.set_axis_off()
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=15, fontweight="bold",
            color=cfg.INK_PRIMARY)
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11, color=cfg.INK_SECONDARY,
            wrap=True)
    ax.text(0.5, 0.28, "pending — auto-populates once the metric is computed",
            ha="center", va="center", fontsize=9.5, color=cfg.INK_MUTED, style="italic")
    pfig = go.Figure()
    pfig.update_layout(title=title,
                       annotations=[dict(text=message + "<br><i>pending — auto-populates once "
                                         "the metric is computed</i>", showarrow=False,
                                         x=0.5, y=0.5, font=dict(size=14, color=cfg.INK_SECONDARY))],
                       xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig, pfig


# ------------------------------------------------------------------ save
# ------------------------------------------------------------------ click-to-open
# Clicking a plot element opens a design's 3 refolds (RNA-aligned) in ChimeraX via the local
# open_refolds bridge (:45001). Shared window.openRefolds([id,af3,b2,esm]) fires the request;
# each plot family adds its own handler (scatter/violin points, venn regions, table rows).
_BRIDGE = "http://127.0.0.1:45001/open"
_BRIDGE_HINT = "open_refolds bridge not reachable on :45001 - start it with:  openrefolds"

OPEN_REFOLDS_JS = r"""
window.openRefolds = function(refs){
  if(!refs || !refs.length) return;
  var q = 'id=' + encodeURIComponent(refs[0]||'');
  if(refs[1]) q += '&af3=' + encodeURIComponent(refs[1]);
  if(refs[2]) q += '&b2='  + encodeURIComponent(refs[2]);
  if(refs[3]) q += '&esm=' + encodeURIComponent(refs[3]);
  fetch('__BRIDGE__?' + q, {mode:'no-cors'})
    .catch(function(){ alert('__HINT__'); });
};
""".replace("__BRIDGE__", _BRIDGE).replace("__HINT__", _BRIDGE_HINT)

# scatter / violin: click a point -> open its design (customdata = [id, af3, b2, esm])
CLICK_JS = OPEN_REFOLDS_JS + r"""
(function(){
  function attach(){
    var gd = document.querySelector('.plotly-graph-div');
    if(!gd || !gd.on){ setTimeout(attach, 300); return; }
    gd.on('plotly_click', function(data){
      // pick the first clicked point that carries a design (skips the pass-box mesh etc.)
      for(var i=0;i<data.points.length;i++){
        var cd = data.points[i].customdata;
        if(cd && cd.length){ window.openRefolds(cd); return; }
      }
    });
  }
  attach();
})();
"""

# venn: click a region -> open a RANDOM design from that region. The region's clickable point
# carries customdata = [label, JSON list of [id,af3,b2,esm]]; each click re-samples.
VENN_CLICK_JS = OPEN_REFOLDS_JS + r"""
(function(){
  function attach(){
    var gd = document.querySelector('.plotly-graph-div');
    if(!gd || !gd.on){ setTimeout(attach, 300); return; }
    gd.on('plotly_click', function(data){
      var cd = data.points[0].customdata; if(!cd || !cd[1]) return;
      var lst; try { lst = JSON.parse(cd[1]); } catch(e){ return; }
      if(!lst.length){ alert('no designs in region ' + (cd[0]||'')); return; }
      window.openRefolds(lst[Math.floor(Math.random()*lst.length)]);
    });
  }
  attach();
})();
"""


def clickable_table_html(header, rows, row_refs, title, subtitle=None):
    """Standalone HTML table whose rows are click-to-open. row_refs[i] = [id,af3,b2,esm]
    for row i (a falsy entry -> that row is not clickable)."""
    th = "".join(f"<th>{_escape(str(c))}</th>" for c in header)
    body = []
    for r, refs in zip(rows, row_refs):
        tds = "".join(f"<td>{_escape(str(c))}</td>" for c in r)
        if refs and refs[0]:
            body.append(f"<tr onclick='window.openRefolds({json.dumps(refs)})' "
                        f'title="click - open 3 refolds in ChimeraX">{tds}</tr>')
        else:
            body.append(f"<tr>{tds}</tr>")
    sub = f'<div class="sub">{_escape(subtitle)}</div>' if subtitle else ""
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{_escape(title)}</title>
<style>body{{font:14px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:26px;color:#0b0b0b}}
h2{{font-size:18px;margin:0 0 4px}} .sub{{color:#6b6a66;font-size:12px;margin-bottom:12px}}
table{{border-collapse:collapse;width:100%}} th,td{{padding:6px 10px;border-bottom:1px solid #e1e0d9;text-align:center}}
th{{background:#f0efec;font-weight:700}} tr[onclick]{{cursor:pointer}} tr[onclick]:hover{{background:#eaf3ff}}
.hint{{color:#898781;font-size:12px;margin-top:12px}}</style></head><body>
<h2>{_escape(title)}</h2>{sub}
<table><thead><tr>{th}</tr></thead><tbody>
{chr(10).join(body)}
</tbody></table>
<p class="hint">Click a row - its 3 refolds open (RNA-aligned) in ChimeraX (needs <code>openrefolds</code> running).</p>
<script>{OPEN_REFOLDS_JS}</script>
</body></html>"""


def save_pair(mpl_fig, plotly_fig, name: str, png_dir=None, html_dir=None,
              click_js=None, html_override=None) -> tuple[Path, Path]:
    """Write one PNG (matplotlib) and one HTML, named after the script. `click_js` injects a
    click handler before </body> (scatter/violin/venn). `html_override` writes that HTML string
    verbatim instead of the plotly figure (used by the clickable tables)."""
    png_dir = Path(png_dir or cfg.PNG_DIR)
    html_dir = Path(html_dir or cfg.HTML_DIR)
    png_dir.mkdir(parents=True, exist_ok=True)
    html_dir.mkdir(parents=True, exist_ok=True)
    png = png_dir / f"{name}.png"
    html = html_dir / f"{name}.html"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mpl_fig.savefig(png)
        plt.close(mpl_fig)
    if html_override is not None:
        html.write_text(html_override)
    else:
        plotly_fig.write_html(html, include_plotlyjs="cdn", full_html=True)
        if click_js:
            txt = html.read_text()
            if "</body>" in txt:
                txt = txt.replace("</body>", f"<script>\n{click_js}\n</script>\n</body>", 1)
                html.write_text(txt)
    return png, html


def funnel_figure(labels, counts, totals, title, subtitle=None):
    """Horizontal stacked funnel (reference style): stages top→bottom, x = Number of
    Designs, each bar segmented by stack; total count centered + 'count (pct%)' to the
    right; legend lists each stack's Final count. plots 09 (avg) / 10 (strict) share it."""
    stacks = list(counts)
    n = len(labels)
    start = totals[0] if totals[0] else 1
    y = list(range(n))[::-1]            # first stage on top

    fig, ax = plt.subplots(figsize=(12, 0.52 * n + 2.6))
    for i in range(n):
        left = 0
        for s in stacks:
            v = counts[s][i]
            if v:
                ax.barh(y[i], v, left=left, height=0.72, color=cfg.color(s),
                        edgecolor=cfg.SURFACE, linewidth=0.7)
                left += v
        tot = totals[i]
        if tot:
            ax.text(tot / 2, y[i], f"{tot}", ha="center", va="center", color="white",
                    fontweight="bold", fontsize=9)
        ax.text(tot + start * 0.012, y[i], f"{tot}  ({100.0 * tot / start:.1f}%)",
                ha="left", va="center", fontsize=9, color=cfg.INK_PRIMARY)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Number of Designs"); ax.set_xlim(0, start * 1.16)
    ax.grid(axis="y", visible=False); ax.spines["left"].set_visible(False)

    finals = {s: counts[s][-1] for s in stacks}
    handles = [plt.Rectangle((0, 0), 1, 1, fc=cfg.color(s), ec=cfg.MARK_EDGE) for s in stacks]
    labs = [f"{display(s)}    {finals[s]}" for s in stacks]
    ax.legend(handles, labs, title="Stack        Final", loc="upper left",
              bbox_to_anchor=(1.01, 1.0), fontsize=9, title_fontsize=10, framealpha=0.9)
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.97)
    if subtitle:
        fig.text(0.5, 0.895, subtitle, ha="center", fontsize=8.5, color=cfg.INK_MUTED)
    fig.subplots_adjust(left=0.19, right=0.8, top=0.83, bottom=0.12)

    # plotly: horizontal stacked, stages reversed so the first is on top
    ry = labels[::-1]
    pfig = go.Figure()
    for s in stacks:
        pfig.add_trace(go.Bar(y=ry, x=counts[s][::-1], name=f"{display(s)} · {finals[s]}",
                              orientation="h", marker_color=cfg.color(s),
                              marker_line_color=cfg.SURFACE, marker_line_width=0.7,
                              hovertemplate=f"{display(s)}<br>%{{y}}: %{{x}} designs<extra></extra>"))
    rtot = totals[::-1]
    ann = [dict(x=rtot[i], y=ry[i], text=f"  {rtot[i]} ({100.0 * rtot[i] / start:.1f}%)",
                showarrow=False, xanchor="left", font=dict(size=10, color=cfg.INK_PRIMARY))
           for i in range(n)]
    pfig.update_layout(title=title + (f"<br><sub>{subtitle}</sub>" if subtitle else ""),
                       xaxis_title="Number of Designs", xaxis_range=[0, start * 1.18],
                       barmode="stack", legend_title="Stack · Final", annotations=ann)
    return fig, pfig


def aa_nt_heatmap(map_csv, value_col, n_designs, title, subtitle=None):
    """Per-stack (amino-acid x RNA-nt) heatmap; cell = mean count/design of value_col.
    Static = matplotlib panels, interactive = plotly subplots (shared color scale)."""
    import pandas as pd
    from plotly.subplots import make_subplots
    m = pd.read_csv(map_csv, low_memory=False)
    m["stack"] = m["stack"].map(lambda x: cfg.STACK_ALIASES.get(x, x))
    m = m[m["stack"].isin(cfg.STACK_ORDER)]
    stacks = [s for s in cfg.STACK_ORDER if s in set(m["stack"])]
    nts = list(range(1, cfg.RNA_NT_COUNT + 1))
    mats, vmax = {}, 0.0
    for s in stacks:
        sub = m[m["stack"] == s]
        piv = (sub.groupby(["aa", "nt_pos"])[value_col].sum().unstack(fill_value=0)
               .reindex(index=cfg.AA_ORDER, columns=nts, fill_value=0))
        mat = piv.to_numpy(float) / max(1, n_designs.get(s, 1))
        mats[s] = mat; vmax = max(vmax, float(mat.max()))
    vmax = vmax or 1.0
    seq = cfg.RNA_SEQ if len(cfg.RNA_SEQ) == cfg.RNA_NT_COUNT else None
    xlab = [f"{q}\n{seq[q-1]}" if seq else str(q) for q in nts]

    fig, axes = plt.subplots(1, len(stacks), figsize=(2.9 * len(stacks) + 1.3, 6.3), sharey=True)
    if len(stacks) == 1:
        axes = [axes]
    im = None
    for ax, s in zip(axes, stacks):
        im = ax.imshow(mats[s], aspect="auto", cmap=cfg.HEATMAP_CMAP, vmin=0, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(display(s), fontsize=8.5)
        ax.set_xticks(range(len(nts))); ax.set_xticklabels(xlab, fontsize=6)
        ax.set_xlabel("RNA nt", fontsize=8); ax.grid(False)
    axes[0].set_yticks(range(len(cfg.AA_ORDER))); axes[0].set_yticklabels(cfg.AA_ORDER, fontsize=7)
    axes[0].set_ylabel("amino acid")
    fig.colorbar(im, ax=axes, fraction=0.015, pad=0.01, label="mean / design")
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)
    if subtitle:
        fig.text(0.5, 0.945, subtitle, ha="center", fontsize=8.5, color=cfg.INK_MUTED)

    pfig = make_subplots(rows=1, cols=len(stacks), subplot_titles=[display(s) for s in stacks],
                         shared_yaxes=True, horizontal_spacing=0.015)
    for c, s in enumerate(stacks, start=1):
        pfig.add_trace(go.Heatmap(z=mats[s], x=[str(n) for n in nts], y=cfg.AA_ORDER,
                                  coloraxis="coloraxis",
                                  hovertemplate=f"{display(s)}<br>aa %{{y}} nt %{{x}}<br>mean %{{z:.3f}}<extra></extra>"),
                       row=1, col=c)
    pfig.update_layout(title=title + (f"<br><sub>{subtitle}</sub>" if subtitle else ""),
                       coloraxis=dict(cmin=0, cmax=vmax, colorscale="Magma",
                                      colorbar=dict(title="mean/design")))
    pfig.update_yaxes(autorange="reversed")
    return fig, pfig


def designs_table(header, rows, row_stacks, title, subtitle=None):
    """Ranking table (static + interactive), one row per design; the 'Stack' cell is tinted
    by stack color. row_stacks[i] = canonical stack of row i. Used by plots 01_1 / 01_2."""
    n = len(rows)
    sc = header.index("Stack") if "Stack" in header else 1

    def _w(h):  # relative column widths — text columns wide, numbers/ranks narrow
        if h == "Design":
            return 2.6
        if h == "Stack":
            return 1.7
        if h in ("Rank", "In-stack #", "global #"):
            return 0.62
        return 0.9
    weights = [_w(h) for h in header]
    colw = [w / sum(weights) for w in weights]

    fig, ax = plt.subplots(figsize=(13, 0.5 * max(n, 1) + 1.5))
    ax.set_axis_off()
    tbl = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center", colWidths=colw)
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if c in (sc, header.index("Design") if "Design" in header else -1):
            cell.set_text_props(ha="left"); cell.PAD = 0.02
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d9d8d2")
        if r == 0:
            cell.set_facecolor("#f0efec"); cell.set_text_props(weight="bold", color=cfg.INK_PRIMARY)
        elif c == sc and (r - 1) < len(row_stacks):
            cell.set_facecolor(cfg.color(row_stacks[r - 1]) + "22")
            cell.set_text_props(weight="bold")
    ax.set_title(title, pad=14)
    if subtitle:
        fig.text(0.5, 0.02, subtitle, ha="center", fontsize=8.5, color=cfg.INK_MUTED)

    cols_t = list(zip(*rows)) if rows else [[] for _ in header]
    fill = []
    for h in header:
        if h == "Stack":
            fill.append([rgba(cfg.color(s), 0.13) for s in row_stacks])
        else:
            fill.append(["#ffffff"] * n)
    pfig = go.Figure(go.Table(
        header=dict(values=[f"<b>{h}</b>" for h in header], fill_color="#f0efec",
                    font=dict(color=cfg.INK_PRIMARY, size=12), align="center"),
        cells=dict(values=cols_t, fill_color=fill, align="center", font=dict(size=11), height=26)))
    pfig.update_layout(title=title + (f"<br><sub>{subtitle}</sub>" if subtitle else ""))
    return fig, pfig


def box_edges(x0, x1, y0, y1, z0, z1):
    """The 12 edges of the cuboid [x0,x1]×[y0,y1]×[z0,z1] as (x,y,z) lists with None separators —
    for a clickable-through WIREFRAME box (go.Scatter3d mode='lines'); a filled Mesh3d would
    occlude the 3-D click picker and swallow point clicks."""
    C = {"000": (x0, y0, z0), "100": (x1, y0, z0), "110": (x1, y1, z0), "010": (x0, y1, z0),
         "001": (x0, y0, z1), "101": (x1, y0, z1), "111": (x1, y1, z1), "011": (x0, y1, z1)}
    edges = [("000", "100"), ("100", "110"), ("110", "010"), ("010", "000"),
             ("001", "101"), ("101", "111"), ("111", "011"), ("011", "001"),
             ("000", "001"), ("100", "101"), ("110", "111"), ("010", "011")]
    x, y, z = [], [], []
    for a, b in edges:
        pa, pb = C[a], C[b]
        x += [pa[0], pb[0], None]; y += [pa[1], pb[1], None]; z += [pa[2], pb[2], None]
    return x, y, z


def refolder_cube(d, *, value_min, metric_label, hi=1.0):
    """3-D scatter of a metric across the 3 refolders (x=AF3, y=Boltz-2, z=ESMFold2) + a
    lightly-colored 'all three pass' box at [value_min, hi]^3. `d` has columns af3/b2/esm +
    stack + design_id + *_cif_path. Returns (mpl_fig, plotly_fig); points are click-to-open."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    lo = value_min
    stacks = ordered_stacks(d)
    passmask = (d["af3"] >= lo) & (d["b2"] >= lo) & (d["esm"] >= lo)
    n_pass = int(passmask.sum())
    xl, yl, zl = f"AF3 {metric_label}", f"Boltz-2 {metric_label}", f"ESMFold2 {metric_label}"

    # ---- static ----
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    for s in stacks:
        g = d[d["stack"] == s]
        ax.scatter(g["af3"], g["b2"], g["esm"], s=8, alpha=0.5, color=cfg.color(s),
                   edgecolor="none", label=display(s))
    c = [(lo, lo, lo), (lo, lo, hi), (lo, hi, lo), (lo, hi, hi),
         (hi, lo, lo), (hi, lo, hi), (hi, hi, lo), (hi, hi, hi)]
    faces = [[c[0], c[1], c[3], c[2]], [c[4], c[5], c[7], c[6]], [c[0], c[1], c[5], c[4]],
             [c[2], c[3], c[7], c[6]], [c[0], c[2], c[6], c[4]], [c[1], c[3], c[7], c[5]]]
    ax.add_collection3d(Poly3DCollection(faces, facecolor="#2ca02c", alpha=0.10,
                                         edgecolor="#2ca02c", linewidths=0.5))
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_zlabel(zl)
    ax.set_title(f"{metric_label} across refolders (AF3 × Boltz-2 × ESMFold2)\n"
                 f"(green box = all three ≥ {lo:g}, n={n_pass})", fontsize=11)
    ax.legend(fontsize=8, loc="upper left"); ax.view_init(elev=22, azim=-60)

    # ---- interactive ----
    pfig = go.Figure()
    for s in stacks:
        g = d[d["stack"] == s]
        sizes = [7 if p else 3 for p in passmask.loc[g.index]]   # enlarge ONLY the passers (in the box)
        pfig.add_trace(go.Scatter3d(
            x=g["af3"], y=g["b2"], z=g["esm"], mode="markers", name=display(s),
            marker=dict(color=cfg.color(s), size=sizes, opacity=0.7,
                        line=dict(color=cfg.MARK_EDGE, width=0.3)),
            customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path",
                          "esm_cif_path"]].fillna(""),
            hovertemplate=(f"{display(s)}<br>AF3=%{{x:.3f}}  B2=%{{y:.3f}}  ESM=%{{z:.3f}}"
                           "<br>%{customdata[0]}<br><i>click → open 3 refolds</i><extra></extra>")))
    # BOLD wireframe pass box — a filled surface would occlude clicks on the points inside it,
    # so it stays a wireframe (thick + solid so it's clearly visible).
    ex, ey, ez = box_edges(lo, hi, lo, hi, lo, hi)
    pfig.add_trace(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines",
                                line=dict(color="#1a7d1a", width=7),
                                name=f"all three ≥ {lo:g} ({n_pass})",
                                hoverinfo="skip", showlegend=True))
    pfig.update_layout(
        title=f"{metric_label} across refolders — AF3 × Boltz-2 × ESMFold2  "
              f"(green box = all three ≥ {lo:g} · drag to rotate · click → 3 refolds)",
        scene=dict(xaxis_title=xl, yaxis_title=yl, zaxis_title=zl),
        margin=dict(l=0, r=0, t=50, b=0), height=760)
    return fig, pfig

# ------------------------------------------------------------------ refolder pair grid
def _pairgrid_titles(pairs, corrs):
    """subplot_titles for the plotly grid: pair name + Pearson r / Spearman ρ above each
    top-marginal cell (row 1). Later rows get no titles."""
    top = []
    for i, (_xk, _yk, xd, yd) in enumerate(pairs):
        pr, sr = corrs[i]
        top.append(f"{xd} vs {yd}  (r={pr:.2f}, ρ={sr:.2f})")
        top.append("")
    return top


def _corr_pair(a, xk, yk):
    """(pearson, spearman, n) for columns xk,yk of frame `a`; NaNs dropped. Spearman = Pearson
    on average ranks (tie-aware, no scipy). Stats are None if fewer than 3 finite points."""
    b = a[[xk, yk]].dropna()
    if len(b) < 3:
        return None, None, len(b)
    pear = float(np.corrcoef(b[xk].to_numpy(float), b[yk].to_numpy(float))[0, 1])
    spear = float(np.corrcoef(b[xk].rank().to_numpy(), b[yk].rank().to_numpy())[0, 1])
    return pear, spear, len(b)


def refolder_pairgrid(d, *, metric_label, value_min, lo=0.0, hi=1.0):
    """3 pairwise refolder scatters (AF3 vs Boltz-2, AF3 vs ESMFold2, Boltz-2 vs ESMFold2) of
    `metric_label`, each a JOINT plot: main scatter (by stack) + marginal-distribution histograms
    on both axes (PER GENERATOR, stacked; in the interactive figure they are legend-linked, so
    toggling a generator drops it from the scatter AND the marginals), an identity diagonal (grey)
    and threshold lines (red dashed at `value_min`), plus Pearson r + Spearman ρ per pair. A table
    at the bottom breaks r / ρ / n down for every generator × refolder pair. `d` has af3/b2/esm +
    stack + design_id + *_cif_path. Interactive points are click-to-open. Returns (mpl, plotly)."""
    from plotly.subplots import make_subplots
    pairs = [("af3", "b2", "AF3", "Boltz-2"),
             ("af3", "esm", "AF3", "ESMFold2"),
             ("b2", "esm", "Boltz-2", "ESMFold2")]
    stacks = ordered_stacks(d)
    pad = (hi - lo) * 0.02
    rng = [lo - pad, hi + pad]

    overall = [_corr_pair(d, xk, yk) for (xk, yk, _x, _y) in pairs]
    corrs = [(p, s) for (p, s, _n) in overall]

    # ---- correlation table data: rows = generators + "All", cols = pair×(r, ρ) + n ----
    def _fmt(v):
        return f"{v:.2f}" if isinstance(v, float) and v == v else "–"
    header = ["Generator"]
    for (_x, _y, xd, yd) in pairs:
        header += [f"{xd}–{yd} r", f"{xd}–{yd} ρ"]
    header += ["n"]
    table_rows, row_stacks = [], []
    for s in stacks:
        sub = d[d["stack"] == s]
        row = [short_stack(s)]
        for (xk, yk, _x, _y) in pairs:
            pr, sr, _n = _corr_pair(sub, xk, yk)
            row += [_fmt(pr), _fmt(sr)]
        row.append(f"{len(sub):,}")
        table_rows.append(row); row_stacks.append(s)
    row = ["All generators"]
    for (pr, sr, _n) in overall:
        row += [_fmt(pr), _fmt(sr)]
    row.append(f"{len(d):,}")
    table_rows.append(row); row_stacks.append(None)

    # ---------------- static (matplotlib) ----------------
    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.7], hspace=0.34)
    topgs = gs[0].subgridspec(1, 3, wspace=0.3)
    for i, (xk, yk, xd, yd) in enumerate(pairs):
        inner = topgs[i].subgridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                                     wspace=0.06, hspace=0.06)
        ax_top = fig.add_subplot(inner[0, 0])
        ax_main = fig.add_subplot(inner[1, 0], sharex=ax_top)
        ax_right = fig.add_subplot(inner[1, 1], sharey=ax_main)
        for s in stacks:
            g = d[d["stack"] == s]
            ax_main.scatter(g[xk], g[yk], s=9, alpha=0.45, color=cfg.color(s), edgecolor="none")
        ax_main.plot(rng, rng, color=cfg.INK_MUTED, ls="-", lw=1.0, zorder=1)   # identity
        add_mpl_threshold(ax_main, value_min, "x")
        add_mpl_threshold(ax_main, value_min, "y")
        ax_main.set_xlim(rng); ax_main.set_ylim(rng)
        ax_main.set_xlabel(f"{xd} {metric_label}"); ax_main.set_ylabel(f"{yd} {metric_label}")
        pr, sr = corrs[i]
        ax_main.text(0.04, 0.96, f"Pearson r = {pr:.2f}\nSpearman ρ = {sr:.2f}\nn = {len(d):,}",
                     transform=ax_main.transAxes, va="top", ha="left", fontsize=9,
                     color=cfg.INK_SECONDARY,
                     bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=cfg.GRID, alpha=0.85))
        _xs = [d[d["stack"] == s][xk].dropna() for s in stacks]
        _ys = [d[d["stack"] == s][yk].dropna() for s in stacks]
        _cols = [cfg.color(s) for s in stacks]
        ax_top.hist(_xs, bins=40, range=(lo, hi), stacked=True, color=_cols, edgecolor="white", linewidth=0.15)
        ax_right.hist(_ys, bins=40, range=(lo, hi), orientation="horizontal", stacked=True,
                      color=_cols, edgecolor="white", linewidth=0.15)
        ax_top.set_title(f"{xd} vs {yd}", fontsize=11, fontweight="bold")
        for a in (ax_top, ax_right):
            a.grid(False)
        ax_top.tick_params(labelbottom=False, left=False, labelleft=False)
        ax_right.tick_params(labelleft=False, bottom=False, labelbottom=False)
        for sp in ("left", "bottom"):
            ax_top.spines[sp].set_visible(False); ax_right.spines[sp].set_visible(False)
    # bottom correlation table
    ax_tbl = fig.add_subplot(gs[1]); ax_tbl.axis("off")
    ax_tbl.set_title("Pearson r / Spearman ρ per generator × refolder pair", fontsize=10, pad=8)
    tbl = ax_tbl.table(cellText=table_rows, colLabels=header, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d9d8d2")
        if r == 0:
            cell.set_facecolor("#f0efec"); cell.set_text_props(weight="bold", color=cfg.INK_PRIMARY)
        elif c == 0 and (r - 1) < len(row_stacks) and row_stacks[r - 1]:
            cell.set_facecolor(cfg.color(row_stacks[r - 1]) + "33"); cell.set_text_props(weight="bold")
    handles = [plt.Line2D([], [], marker="s", ls="", color=cfg.color(s), label=display(s)) for s in stacks]
    fig.legend(handles=handles, loc="lower center", ncol=len(stacks), fontsize=9, framealpha=0.85,
               bbox_to_anchor=(0.5, 0.008))
    fig.suptitle(f"{metric_label} — pairwise refolder agreement   "
                 f"(grey = identity, red dashed = {metric_label} {value_min:g})",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.subplots_adjust(top=0.89, bottom=0.075, left=0.05, right=0.98)

    # ---------------- interactive (plotly) ----------------
    pfig = make_subplots(rows=3, cols=6,
                         specs=[[{}, {}, {}, {}, {}, {}],
                                [{}, {}, {}, {}, {}, {}],
                                [{"type": "table", "colspan": 6}, None, None, None, None, None]],
                         column_widths=[0.26, 0.066, 0.26, 0.066, 0.26, 0.066],
                         row_heights=[0.16, 0.52, 0.32], horizontal_spacing=0.03, vertical_spacing=0.06,
                         subplot_titles=_pairgrid_titles(pairs, corrs))
    _bin = (hi - lo) / 40.0
    for i, (xk, yk, xd, yd) in enumerate(pairs):
        mcol, rcol = 2 * i + 1, 2 * i + 2
        first = (i == 0)
        for s in stacks:
            g = d[d["stack"] == s]
            pfig.add_trace(go.Scattergl(
                x=g[xk], y=g[yk], mode="markers", name=display(s), legendgroup=s, showlegend=first,
                marker=dict(color=cfg.color(s), size=5, opacity=0.6, line=dict(color=cfg.MARK_EDGE, width=0.3)),
                customdata=g[["design_id", "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]].fillna(""),
                hovertemplate=(f"{display(s)}<br>{xd}=%{{x:.3f}}  {yd}=%{{y:.3f}}<br>%{{customdata[0]}}"
                               "<br><i>click → open 3 refolds</i><extra></extra>")), row=2, col=mcol)
        pfig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=cfg.INK_MUTED, width=1),
                                  showlegend=False, hoverinfo="skip"), row=2, col=mcol)
        pfig.add_trace(go.Scatter(x=[value_min, value_min], y=rng, mode="lines",
                                  line=dict(color=cfg.THRESHOLD_LINE, width=1.2, dash="dash"),
                                  showlegend=False, hoverinfo="skip"), row=2, col=mcol)
        pfig.add_trace(go.Scatter(x=rng, y=[value_min, value_min], mode="lines",
                                  line=dict(color=cfg.THRESHOLD_LINE, width=1.2, dash="dash"),
                                  showlegend=False, hoverinfo="skip"), row=2, col=mcol)
        for s in stacks:
            gm = d[d["stack"] == s]
            pfig.add_trace(go.Histogram(x=gm[xk], xbins=dict(start=lo, end=hi, size=_bin),
                                        marker_color=cfg.color(s), legendgroup=s, showlegend=False,
                                        hoverinfo="skip"), row=1, col=mcol)
            pfig.add_trace(go.Histogram(y=gm[yk], ybins=dict(start=lo, end=hi, size=_bin),
                                        marker_color=cfg.color(s), legendgroup=s, showlegend=False,
                                        hoverinfo="skip"), row=2, col=rcol)
        pfig.update_xaxes(range=rng, title_text=f"{xd} {metric_label}", row=2, col=mcol)
        pfig.update_yaxes(range=rng, title_text=f"{yd} {metric_label}", row=2, col=mcol)
        pfig.update_xaxes(range=rng, showticklabels=False, row=1, col=mcol)
        pfig.update_yaxes(showticklabels=False, row=1, col=mcol)
        pfig.update_yaxes(range=rng, showticklabels=False, row=2, col=rcol)
        pfig.update_xaxes(showticklabels=False, row=2, col=rcol)
        pfig.update_xaxes(visible=False, row=1, col=rcol)
        pfig.update_yaxes(visible=False, row=1, col=rcol)
    # bottom table
    _genfill = [rgba(cfg.color(s), 0.18) if s else "#f4f4f0" for s in row_stacks]
    _fills = [_genfill] + [["#ffffff"] * len(table_rows) for _ in range(len(header) - 1)]
    pfig.add_trace(go.Table(
        columnwidth=[2.4] + [1.0] * (len(header) - 1),
        header=dict(values=[f"<b>{h}</b>" for h in header], fill_color="#f0efec",
                    font=dict(color=cfg.INK_PRIMARY, size=11), align="center"),
        cells=dict(values=list(map(list, zip(*table_rows))), fill_color=_fills,
                   align="center", font=dict(size=10), height=22)), row=3, col=1)
    pfig.update_layout(
        title=f"{metric_label} — pairwise refolder agreement   "
              f"(identity line · red dashed = {value_min:g} · click a point → 3 refolds)",
        barmode="stack", height=980,
        legend=dict(orientation="h", y=-0.02, x=0.5, xanchor="center", yanchor="top",
                    groupclick="togglegroup"),
        margin=dict(l=60, r=20, t=90, b=70))
    return fig, pfig
