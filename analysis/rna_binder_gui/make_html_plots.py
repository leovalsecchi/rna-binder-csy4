#!/usr/bin/env python3
"""Interactive (Plotly) HTML versions of the analysis plots.

Per-design SCATTERS are click-to-open: clicking a point fires `open <4 variants>
+ matchmaker` into a RUNNING ChimeraX via its REST server (the same one `opend`
starts on port 45000), using a no-CORS fetch — so the browser doesn't need to
launch anything itself. Histograms are interactive (hover/zoom) but not clickable
(no per-design point).

Prereq for click-to-open: ChimeraX running with REST on :45000 — start it once
with any `opend ...` (or `opendesign.py`). The structures are read from the mount
by ChimeraX locally, so the eulermount must be up.

Usage:
    make_html_plots.py STRUCTURAL_CSV OUTDIR [--mount BASE]
(STRUCTURAL_CSV = an aggregate_metrics --structural CSV, ideally with dockq cols.)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A
import opendesign as OD

REST_PORT = OD.REST_PORT
STACK_COLORS = {
    "RFD3+LiMPNN": "#9C27B0", "BG+LiMPNN": "#2196F3", "BG+BIF": "#64B5F6",
    "BD1+LiMPNN": "#FF9800", "HD+LiMPNN": "#F44336",
}
# Colours for stacks not in STACK_COLORS, so every plot defaults to EVERY stack
# present in the CSV — not just the five known ones.
_FALLBACK_COLORS = ["#26A69A", "#EC407A", "#7E57C2", "#8D6E63", "#789262",
                    "#FF7043", "#5C6BC0", "#9E9D24"]


def stacks_present(sub: pd.DataFrame) -> list[tuple[str, str]]:
    """[(stack, colour), ...] for every stack in `sub`: the known stacks in
    canonical order first, then any extras (alphabetical) with fallback colours."""
    present = set(sub["stack"].dropna().unique())
    out = [(s, STACK_COLORS[s]) for s in STACK_COLORS if s in present]
    extra = sorted(present - set(STACK_COLORS))
    out += [(s, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
            for i, s in enumerate(extra)]
    return out

# JS click handler: read the point's embedded variants, build the ChimeraX
# command (mirrors opendesign's .cxc), fire it at the REST server (no-CORS).
_CLICK_JS = """
var gd = document.getElementById('{plot_id}');
gd.on('plotly_click', function(data){{
  var cd = data.points[0].customdata;
  if(!cd || !cd[0]) return;
  var v; try {{ v = JSON.parse(cd[0]); }} catch(e) {{ return; }}
  if(!v.length) {{ alert('No designed-source path for this point.'); return; }}
  var c = ['close session'], n = v.length;
  for(var i=0;i<n;i++){{ c.push('open ' + JSON.stringify(v[i][2])); c.push('rename #'+(i+1)+' '+v[i][0]); }}
  if(n>1) c.push('matchmaker #2-'+n+' to #1');
  c.push('hide atoms'); c.push('show cartoon'); c.push('color #1-'+n+' & protein gray');
  for(var i=0;i<n;i++){{ c.push('color #'+(i+1)+' & nucleic '+v[i][1]); }}
  c.push('view');
  var url = 'http://127.0.0.1:%d/run?command=' + encodeURIComponent(c.join('; '));
  fetch(url, {{mode:'no-cors'}}).catch(function(){{
    alert('ChimeraX not reachable on :%d — start it once with `opend ...`');
  }});
}});
""" % (REST_PORT, REST_PORT)


def _variants_json(designed_source, mount: str) -> str:
    """JSON [[label,colour,abspath],...] for a design's 4 variants (no fs check)."""
    if not designed_source or pd.isna(designed_source):
        return "[]"
    p = str(designed_source)
    if not p.startswith("/"):
        p = os.path.join(mount, p)
    run = Path(p).parent.parent
    stem = OD._strip_ext(Path(p).name)
    paths = {
        "design": p,
        "af3_msaoff": str(run / "refold_af3/msaoff/designs" / f"{stem}_AF3_top.cif"),
        "af3_msaon": str(run / "refold_af3/msaon/designs" / f"{stem}_AF3_top.cif"),
        "b2": str(run / "refold_b2/msaoff/designs" / f"{stem}_B2_sample0.cif"),
        "esm": str(run / "refold_esm/msaoff/designs" / f"{stem}_ESM_sample0.cif"),
    }
    return json.dumps([[label, colour, paths[label]] for label, colour in OD._VARIANTS])


def clickable_scatter(designs: pd.DataFrame, xcol, ycol, *, title, xlab, ylab,
                      outpath, mount, xclip=None, yclip=None):
    missing = [c for c in (xcol, ycol) if c not in designs.columns]
    if missing:
        print(f"  (skip {outpath.name}: missing column(s) {missing})")
        return
    sub = designs.dropna(subset=[xcol, ycol])
    if sub.empty:
        print(f"  (skip {outpath.name}: no {xcol}/{ycol} data)")
        return
    fig = go.Figure()
    for stack, color in stacks_present(sub):
        g = sub[sub["stack"] == stack]
        x = g[xcol].clip(upper=xclip) if xclip else g[xcol]
        y = g[ycol].clip(upper=yclip) if yclip else g[ycol]
        cd = [[_variants_json(src, mount), key]
              for src, key in zip(g.get("designed_source", [None] * len(g)), g["design_key"])]
        fig.add_trace(go.Scattergl(
            x=x, y=y, mode="markers", name=f"{stack} (n={len(g)})",
            marker=dict(size=7, color=color, opacity=0.55),
            customdata=cd,
            hovertemplate=(f"<b>%{{customdata[1]}}</b><br>{xlab}=%{{x:.2f}}<br>"
                           f"{ylab}=%{{y:.2f}}<br><i>click → open in ChimeraX</i><extra></extra>"),
        ))
    fig.update_layout(title=title + "  (click a point → open its 4 variants in ChimeraX)",
                      xaxis_title=xlab, yaxis_title=ylab, template="plotly_white",
                      legend_title="stack", height=650)
    fig.write_html(str(outpath), include_plotlyjs="cdn", post_script=_CLICK_JS)
    print(f"  wrote {outpath}  ({len(sub)} points)")


def histogram(designs: pd.DataFrame, col, *, title, xlab, outpath, nbins=40):
    if col not in designs.columns:
        print(f"  (skip {outpath.name}: missing column {col})")
        return
    sub = designs.dropna(subset=[col])
    if sub.empty:
        print(f"  (skip {outpath.name}: no {col})")
        return
    fig = go.Figure()
    for stack, color in stacks_present(sub):
        g = sub[sub["stack"] == stack]
        fig.add_trace(go.Histogram(x=g[col], name=f"{stack} (n={len(g)})",
                                   marker_color=color, opacity=0.6, nbinsx=nbins))
    fig.update_layout(title=title, xaxis_title=xlab, yaxis_title="designs",
                      barmode="overlay", template="plotly_white", height=550)
    fig.write_html(str(outpath), include_plotlyjs="cdn")
    print(f"  wrote {outpath}  ({len(sub)} designs)")


def violins_by_stack(designs: pd.DataFrame, metrics, *, title, outpath):
    """Violin distributions per stack — defaults to EVERY stack present in the
    CSV, one subplot per metric.

    `metrics` = [(col, ylabel, clip_or_None), ...]. Metric columns absent from
    the CSV are dropped (that panel is skipped, no crash); a stack with no data
    for a metric is simply omitted from that panel. One shared, click-to-toggle
    legend entry per stack spans all panels.
    """
    metrics = [(c, lab, clip) for c, lab, clip in metrics if c in designs.columns]
    if not metrics:
        print(f"  (skip {outpath.name}: none of the metric columns present)")
        return
    fig = make_subplots(rows=1, cols=len(metrics),
                        subplot_titles=[lab for _, lab, _ in metrics])
    seen: set[str] = set()
    for ci, (col, ylab, clip) in enumerate(metrics, start=1):
        sub = designs.dropna(subset=[col])
        for stack, color in stacks_present(sub):
            g = sub[sub["stack"] == stack]
            if g.empty:
                continue
            y = g[col].clip(upper=clip) if clip else g[col]
            fig.add_trace(go.Violin(
                y=y, name=f"{stack} (n={len(g)})", legendgroup=stack,
                showlegend=stack not in seen, scalegroup=f"{col}",
                line_color=color, fillcolor=color, opacity=0.6,
                box_visible=True, meanline_visible=True, points=False,
                hovertemplate=f"<b>{stack}</b><br>{ylab}=%{{y:.2f}}<extra></extra>",
            ), row=1, col=ci)
            seen.add(stack)
        fig.update_yaxes(title_text=ylab,
                         range=[0, clip * 1.05] if clip else None, row=1, col=ci)
        fig.update_xaxes(showticklabels=False, row=1, col=ci)
    fig.update_layout(title=title, template="plotly_white", height=600,
                      legend_title="stack", violinmode="group")
    fig.write_html(str(outpath), include_plotlyjs="cdn")
    print(f"  wrote {outpath}  ({len(seen)} stacks)")


def pick_refolder(designs: pd.DataFrame, choice: str) -> str:
    """Refolder prefix ('af3'/'b2') for the violin. 'auto' picks whichever
    refolder's Cα-RMSD covers more stacks (ties -> af3), so the default plot
    shows every available stack instead of silently dropping un-refolded ones."""
    if choice in ("af3", "b2", "esm"):
        return choice
    def coverage(prefix: str) -> int:
        col = f"{prefix}_rmsd"
        return designs.dropna(subset=[col])["stack"].nunique() if col in designs else 0
    cov = {p: coverage(p) for p in ("af3", "b2", "esm")}
    pick = max(("af3", "b2", "esm"), key=lambda p: cov[p])   # ties -> af3, then b2, then esm
    print(f"  violin refolder=auto -> {pick.upper()} "
          f"(AF3 covers {cov['af3']} stacks, B2 {cov['b2']}, ESM {cov['esm']})")
    return pick


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("outdir")
    ap.add_argument("--mount", default=OD.MOUNT_BASE, help="mount root for relative designed_source")
    ap.add_argument("--violin-refolder", choices=["auto", "af3", "b2", "esm"], default="auto",
                    help="Refolder basis for the RMSD-by-stack violin. 'auto' (default) "
                         "uses whichever refolder covers the most stacks, so every "
                         "available stack appears.")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    designs = A.to_per_design(A.load_metrics(args.csv))
    designs = A.add_composite(A.add_composite(A.add_composite(designs, "af3"), "b2"), "esm")
    print(f"{len(designs)} designs from {args.csv}")

    # Clickable per-design scatters
    clickable_scatter(designs, "af3_rmsd", "af3_iptm",
                      title="AF3 iPTM vs Cα-RMSD", xlab="af3_rmsd", ylab="af3_iptm",
                      outpath=out / "iptm_vs_rmsd.html", mount=args.mount, xclip=15)
    clickable_scatter(designs, "af3_rmsd", "af3_lrmsd",
                      title="DockQ pose (ligand-RMSD) vs fold (Cα-RMSD)",
                      xlab="af3_rmsd (fold)", ylab="af3_lrmsd (pose)",
                      outpath=out / "dockq_lrmsd_vs_carmsd.html", mount=args.mount,
                      xclip=15, yclip=30)
    clickable_scatter(designs, "b2_iptm", "af3_iptm",
                      title="AF3 vs B2 iPTM", xlab="b2_iptm", ylab="af3_iptm",
                      outpath=out / "af3_vs_b2_iptm.html", mount=args.mount)
    # Interactive (non-clickable) distributions
    histogram(designs, "af3_iptm", title="AF3 iPTM distribution", xlab="af3_iptm",
              outpath=out / "iptm_distribution.html")
    histogram(designs, "af3_rmsd", title="AF3 Cα-RMSD distribution", xlab="af3_rmsd (Å)",
              outpath=out / "rmsd_distribution.html")
    histogram(designs, "af3_dockq", title="DockQ score distribution", xlab="af3_dockq",
              outpath=out / "dockq_distribution.html")

    # RMSD distributions by stack (all available stacks): fold (Cα) + pose (ligand).
    # Default refolder = whichever covers the most stacks, so every stack shows.
    rp = pick_refolder(designs, args.violin_refolder)
    violins_by_stack(
        designs,
        [(f"{rp}_rmsd", "Cα-RMSD vs designed (Å)", 15),
         (f"{rp}_lrmsd", "DockQ ligand-RMSD (Å)", 30)],
        title=f"RMSD by stack — fold (Cα) vs pose (ligand), {rp.upper()} best sample/design",
        outpath=out / "rmsd_violins_by_stack.html")

    # index
    links = "\n".join(f'<li><a href="{p.name}">{p.name}</a></li>'
                      for p in sorted(out.glob("*.html")) if p.name != "index.html")
    (out / "index.html").write_text(
        f"<h2>Interactive plots</h2><p>Click a scatter point to open its 4 variants "
        f"in ChimeraX (needs ChimeraX running with REST on :{REST_PORT} — run any "
        f"<code>opend ...</code> once).</p><ul>{links}</ul>")
    print(f"index: {out/'index.html'}")


if __name__ == "__main__":
    main()
