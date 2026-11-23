"""Fig 3, assembled on an A4 page from each panel's newest file.

    python code/Fig3_assembly/make_figure.py [outdir]     # default working/Fig3_assembly

Layout (Leo, 2026-11-13; supersedes the deck's a / b / d):
    a            filter funnel, criterion renders over their columns
    b            secondary-structure cluster map, campaign | gated
    c | d        3d's two plots side by side: c the bars (backbones and fold
                 clusters each rule keeps), d the columns (stack mix per set)
    e            GPU cost against pass rate (Fig3e), full width at the bottom

V1 "as_sketched"  that arrangement with every panel at the size it was drawn,
                  text untouched; shrunk as a whole if it overflows the page.
V1 "adjusted"     V1's arrangement and panel sizes (Leo's pick for the book),
                  with V2's text and strokes: every label at the house size,
                  every stroke at least LINE_MIN. If re-set text makes the page
                  overflow, every panel comes down by the same factor.
V2 "adjusted"     a and b fill the block, c and d share a row at one height,
                  e is redrawn by its own script to fill what is left - every
                  label at the house size, every stroke at least LINE_MIN.

Panels
  a  output/Fig3a_filter_funnel/Fig3a_funnel_3to1_box.svg with the five criterion
     renders (selection/Fig3a_box_*.png) set into the empty top of their columns,
     as on the deck. Composited here, into _a_composite/.
  b  output/Fig3b_sse_cluster_map/sse_cluster_map_compare.svg
  c, d  output/Fig3d_selection_difference/selection_difference.svg, cut in two
     at the blank band between its plots (found on every run) - _cd_split/.
  e  output/Fig3e_gpu_efficiency/Fig3e_gpu_efficiency.svg; V2 redraws it with
     its script's size variant (W_MM H_MM FONT_PT OUTDIR) into _e_variant/.
"""
import base64
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common.page_assembly import (GAP_H, GAP_V, MM, PAGE_H, TOP, W, T_HOUSE,  # noqa: E402
                                  Panel, fit, floor, ink, page, pair, report,
                                  variant, vb, window, write)
import numpy as np                                            # noqa: E402
from PIL import Image                                         # noqa: E402

# 2026-11-17 (Leo): "Add the venn as 3a". The funnel panel is replaced by the
# four-set Venn, which already carries its own criterion renders and the
# foot-clash reduction in a right-hand column, so build_a's render-compositing
# is no longer needed for panel a. The funnel source and build_a are kept below
# so the old panel can be restored by pointing A_SRC back at build_a().
VENN = REPO / "output/Fig3b_filter_venn/Fig3b_filter_venn.svg"
PICKS = REPO / "output/Fig3cd_pools/panel_picks.svg"
MATRIX = REPO / "output/Fig3g_fold_coverage/panel.svg"
BACKEND = REPO / "output/Fig3a_backend_agreement/Fig3a_backend_agreement.svg"
FUNNEL = REPO / "output/Fig3a_filter_funnel/Fig3a_funnel_3to1_box.svg"
RENDERS = [REPO / f"selection/Fig3a_box_{k}.png"            # left to right, as the
           for k in ("scrmsd", "dockq", "aarun", "iptm", "footclash")]   # columns
# 2026-11-17 (Leo): the SSE cluster map is replaced by Alpha/Beta/Loop ternary
# density plots, two per backbone design method (before and after filtering).
# The old source is kept so the cluster map can be restored by swapping B back.
B = REPO / "output/Fig3c_sse_ternary/panel_c.svg"
B_WAS = REPO / "output/Fig3b_sse_cluster_map/sse_cluster_map_compare.svg"
D = REPO / "output/Fig3d_selection_difference/selection_difference.svg"
E = REPO / "output/Fig3e_gpu_efficiency/Fig3e_gpu_efficiency.svg"
E_SCRIPT = REPO / "code/Fig3e_gpu_efficiency/make_figure.py"
SPLIT_PAD = 4.0                # pt of air kept round each half of 3d
BOX_FILL = "#f0f6fa"           # the funnel's pale criterion columns
PAD = 3.0                      # pt, render to its column's top / to its label
RENDER_PX = 700                # embedded width: ~600 dpi at the printed size


def build_a(workdir):
    """The funnel with each criterion's render in the empty top of its column."""
    t = FUNNEL.read_text()
    boxes = []
    for m in re.finditer(r'<path d="([^"]+)"[^>]*fill: ?%s' % BOX_FILL, t, re.S):
        n = [float(v) for v in re.findall(r"-?\d+\.?\d*", m.group(1))]
        boxes.append((min(n[0::2]), min(n[1::2]), max(n[0::2]), max(n[1::2])))
    boxes.sort()
    assert len(boxes) == len(RENDERS) == 5, f"{len(boxes)} columns, {len(RENDERS)} renders"
    lab = [(float(m.group(1)), float(m.group(2))) for m in re.finditer(
        r'font-size: ([\d.]+)px[^>]*?y="([\d.]+)"[^>]*>[^<]*[≤≥]', t)]
    assert len(lab) == 5, f"{len(lab)} criterion labels"
    size, base = lab[0]
    top_of_label = base - 0.716 * size
    imgs = []
    with tempfile.TemporaryDirectory() as d:
        for (x0, y0, x1, _), png in zip(boxes, RENDERS):
            small = Path(d) / png.name
            subprocess.run(["magick", png, "-trim", "+repage",          # the renders carry
                            "-bordercolor", "none", "-border", "1%",          # wide clear margins
                            "-resize", f"{RENDER_PX}x", small], check=True)
            w_px, h_px = (int(v) for v in subprocess.run(
                ["magick", "identify", "-format", "%w %h", small],
                check=True, capture_output=True, text=True).stdout.split())
            sw, sh = x1 - x0 - 2 * PAD, top_of_label - y0 - 2 * PAD
            k = min(sw / w_px, sh / h_px)
            w, h = w_px * k, h_px * k
            x, y = (x0 + x1 - w) / 2, y0 + PAD + (sh - h) / 2
            b64 = base64.b64encode(small.read_bytes()).decode()
            imgs.append(f'<image x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" '
                        f'xlink:href="data:image/png;base64,{b64}"/>')
    workdir.mkdir(exist_ok=True)
    out = workdir / "Fig3a_funnel_with_renders.svg"
    i = t.rindex("</svg>")
    out.write_text(t[:i] + "".join(imgs) + t[i:], encoding="utf-8")
    return out


def split_3d(workdir):
    """3d's two plots as two files: cut at the widest blank band between them."""
    v = vb(D.read_text())
    px = 2400
    with tempfile.TemporaryDirectory() as d:
        png = Path(d) / "d.png"
        subprocess.run(["rsvg-convert", "-w", str(px), "-b", "white", D, "-o", png], check=True)
        dark = np.asarray(Image.open(png).convert("L")) < 245
    k = v[2] / px
    xs, rows = np.where(dark.any(0))[0], dark.any(1)
    ys = np.where(rows)[0]
    gaps, start = [], None
    for y in range(ys.min(), ys.max() + 1):
        if not rows[y] and start is None:
            start = y
        elif rows[y] and start is not None:
            gaps.append((y - start, start, y)); start = None
    _, g0, g1 = max(gaps)
    assert (g1 - g0) * k > 10, f"no clear band between 3d's plots ({(g1 - g0) * k:.1f} pt)"
    cut, x0 = (g0 + g1) / 2 * k, max(xs.min() * k - SPLIT_PAD, v[0])
    workdir.mkdir(exist_ok=True)
    c = window(D, (x0, v[1], v[0] + v[2] - x0, cut - v[1]), workdir / "Fig3c_bars.svg")
    d = window(D, (x0, cut, v[0] + v[2] - x0, v[1] + v[3] - cut), workdir / "Fig3d_columns.svg")
    return c, d


def v1(a_src, c_src, d_src):
    a = Panel("a", a_src).at(1); a.at(W / a.iw)
    b = Panel("b", B).at(1); b.at(W / b.iw)
    c, d = Panel("c", c_src).at(1.0), Panel("d", d_src).at(1.0)   # as drawn
    e = Panel("e", E).at(1.0)
    return [[a], [b], [c, d], [e]]


def v2(a_src, c_src, d_src, T, workdir):
    """2026-11-17 (Leo, from his original list): the AF3-vs-ESMFold2 backend
    agreement scatters become panel a and everything else shifts one letter —
    "New panel A. Rename all the old ones one letter later."

    a backend agreement / b Venn / c ternaries / d coverage matrix /
    e pies + f GPU.
    """
    a = fit(lambda T_, s: Panel("a", BACKEND, T_, s), T, W)
    b = fit(lambda T_, s: Panel("b", a_src, T_, s), T, W)
    c = fit(lambda T_, s: Panel("c", B, T_, s), T, W)
    d = fit(lambda T_, s: Panel("d", MATRIX, T_, s), T, W)
    e = Panel("e", d_src, T, 1.0)                     # the pies, on the left
    h_f = PAGE_H - 2 * TOP - (a.h + b.h + c.h + d.h) - 4 * GAP_V - 1.5 * MM
    F_W, F_H_MAX = 82 * MM, 45 * MM
    if h_f > F_H_MAX:
        print(f"f capped at {F_H_MAX / MM:.0f} mm "
              f"({(h_f - F_H_MAX) / MM:.0f} mm returned to the foot)")
        h_f = F_H_MAX
    assert h_f >= 30 * MM, f"only {h_f / MM:.1f} mm left for f (floor 30)"
    f = variant("f", E_SCRIPT, E.name, F_W, h_f, T, workdir, sys.executable)
    e.at(min(1.0, (W - GAP_H - f.w) / e.iw))          # pies left, GPU right
    rows = [[a], [b], [c], [d], [e, f]]
    for q, lo in ((a, 150 * MM), (d, 150 * MM), (e, 60 * MM)):
        assert q.w >= lo, (
            f"{q.key} is only {q.w / MM:.1f} mm wide - it is being sized as "
            "part of a pair it no longer belongs to")
    for p in sum(rows, []):              # last, once every scale is final
        p.retouch(lambda u, s=p.s: floor(u, s))
    return rows


def v1_adjusted(a_src, c_src, d_src, T):
    """V1's layout at the house standard: each panel at its V1 printed size
    (page fit included), text re-set to T there, strokes floored."""
    probe = v1(a_src, c_src, d_src)
    page(probe, "probe")                 # applies V1's page fit to each scale
    s1 = {p.key: p.s for p in sum(probe, [])}
    src = {"a": a_src, "b": B, "c": c_src, "d": d_src, "e": E}
    keys = (("a",), ("b",), ("c", "d"), ("e",))
    room = (PAGE_H - 2 * TOP) * 0.995    # a hair of slack for the stroke floor
    k = 1.0
    for _ in range(4):
        rows = [[Panel(q, src[q], T, s1[q] * k).at(s1[q] * k) for q in r] for r in keys]
        height = sum(max(p.h for p in r) for r in rows) + GAP_V * (len(rows) - 1)
        if height <= room:
            break
        k *= room / height
    for p in sum(rows, []):
        p.retouch(lambda u, s=p.s: floor(u, s))
    print(f"V1 adjusted: panels at {k:.3f} of their V1 size")
    return rows


def main(outdir):
    outdir.mkdir(exist_ok=True)
    a_src = VENN          # was build_a(outdir / "_a_composite"), the funnel
    assert a_src.exists(), (
        f"{a_src} is missing: run code/Fig3b_filter_venn/make_figure.py first")
    # 2026-11-17 (Leo): c and d are replaced. c = the 15 largest k30 fold
    # clusters, column = gate-passing backbones with an inset for the
    # ipTM-ranked selection; d = the three pools as stack-mix pies with arrows
    # from the gate pool. Both come from code/Fig3cd_pools/ as separate SVGs, so
    # split_3d() (which cut the old single selection_difference.svg in two) is
    # no longer called. It and D are kept so the old pair can be restored.
    c_src = REPO / "output/Fig3cd_pools/panel_c.svg"
    d_src = REPO / "output/Fig3cd_pools/panel_d.svg"
    for q in (c_src, d_src, MATRIX, BACKEND):
        assert q.exists(), f"{q} is missing: run code/Fig3cd_pools/make_figure.py first"
    rows1 = v1(a_src, c_src, d_src)
    svg1, fit1, bot1 = page(rows1, "Fig3 V1 as sketched")
    rows2 = v2(a_src, c_src, d_src, T_HOUSE, outdir / "_e_variant")
    svg2, fit2, bot2 = page(rows2, "Fig3 V2 adjusted")
    rows3 = v1_adjusted(a_src, c_src, d_src, T_HOUSE)
    svg3, fit3, bot3 = page(rows3, "Fig3 V1 adjusted")
    report("V1 as sketched", rows1, fit1, bot1)
    report(f"V1 adjusted, T = {T_HOUSE:.2f} pt", rows3, fit3, bot3, T_HOUSE)
    assert fit3 == 1.0, "V1 adjusted must fit without shrinking"
    report(f"V2 adjusted, T = {T_HOUSE:.2f} pt", rows2, fit2, bot2, T_HOUSE)
    assert fit2 == 1.0, "V2 must fit without shrinking, or its type is no longer T"
    write(outdir, "Fig3_V1_as_sketched", svg1)
    write(outdir, "Fig3_V2_adjusted", svg2)
    write(outdir, "Fig3_V1_adjusted", svg3)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "working" / "Fig3_assembly")
