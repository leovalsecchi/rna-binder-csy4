"""Fig 2, assembled on an A4 page from each panel's newest file.

    python code/Fig2_assembly/make_figure.py [outdir]     # default working/Fig2_assembly

Same method as code/Fig1_assembly (copied rather than imported while that one
is still being finished; fold the two into code/common once it settles).

V1 "as_built"   no sketch exists for Fig 2, so every panel sits at the size it
                was drawn at (b, drawn 234 mm wide, comes down to the block),
                stacked a / b / c / d, text untouched. Shrunk as a whole if
                it overflows the page.

V2 "adjusted"   every label at ONE size, T = 5.5 pt, the same as Fig 1, so the
                two figures match. a and b fill the block;
                c and d share the last row at one height. Text is re-set to T
                in each file's own units; b also goes from DejaVu Sans to Arial.
                Every stroke under LINE_MIN is raised to it.

Panels
  a  output/Fig2a_pipeline/Fig2a_pipeline.svg
  b  map composite with the 8 paired callouts. The canonical file is PDF/PNG
     only, with glyphs as paths, so its own script is re-run with an OUTDIR
     to get an SVG with live text (output/Fig2b_map_composite/symmetric/zoom/
     make_paired.py OUTDIR) - same drawing, byte-identical PNG/PDF
  c  output/Fig2cd_benchmark/metrics.svg          (violins)
  d  output/Fig2cd_benchmark/gpu_per_design.svg   (bars; _log is the dot version)
  e  output/Fig2d_epitope_funnel/epitope_funnel_columns.svg
     per-nucleotide contact share, campaign / gate / library
  f  output/Fig2d_epitope_funnel/u10_rise.svg
     square side panel: U10 against the other 15 nt, same three stages

REMOVED 2026-11-13 (Leo): e and f are off the page for now - they do not fit
there. Their files and PANELS entries are kept; put them back by restoring the
[e, f] row in v1() and v2().

NOTE on e's letter: it is the FIFTH panel and is lettered e by position, but
its folder is `Fig2d_epitope_funnel` because Leo asked for it as "figure 2
panel d". Panel d is currently the benchmark's GPU-cost bars. Either that
becomes e and this becomes d, or this stays e and the folder gets renamed -
not decided, so nothing else has been renumbered.

Page: A4, 16 mm sides = the 178 mm block in config/style.yaml, 15 mm top.
Panel letters: 8 pt bold, outside the panel, cap-top level with its ink.
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette, style                            # noqa: E402
from common.retext import allowed, retext                    # noqa: E402

MM = 72 / 25.4
PAGE_W, PAGE_H = 210 * MM, 297 * MM
SIDE, TOP = 16 * MM, 15 * MM
W = PAGE_W - 2 * SIDE                    # 504.567 pt, the 178 mm block
GAP_V, GAP_H = 5 * MM, 6 * MM
LETTER_PT, LETTER_DX = 8.0, 4.5 * MM
INK = palette.INK

T_FIG1 = float(style.CFG["font_size_pt"])   # pt, the house size - one across the
                       # thesis (5.17 pt version in working/Fig2_assembly/_5p17/)
LINE_MIN = 0.5         # pt: the strictest common journal floor (BMC, Elsevier)

B_SCRIPT = REPO / "output/Fig2b_map_composite/symmetric/zoom/make_paired.py"
D_SCRIPT = REPO / "code/Fig2cd_benchmark/make_figure.py"
PYTHON = sys.executable

# file, riders ({size in the file: multiple of T}); every other label -> T
PANELS = {
    "a": ("output/Fig2a_pipeline/Fig2a_pipeline.svg", {}),
    "b": (None, {}),                                       # set by _b_svg()
    "c": ("output/Fig2cd_benchmark/metrics.svg", {}),
    # "d" is now the binder-density row (2026-11-18): the target with each
    # backbone generator's design-backbone density over it, one panel per
    # generator. code/Fig2d_binder_density/.
    "d": ("output/Fig2d_binder_density/binder_density_row.svg", {}),
    # the old "d" (gpu_per_design) came off the page 2026-11-17; the script
    # still emits it and Fig2cd_benchmark keeps it, so it can be restored by
    # naming this entry in a row.
    "d_gpu": ("output/Fig2cd_benchmark/gpu_per_design.svg", {}),
    "e": ("output/Fig2d_epitope_funnel/epitope_funnel_columns.svg", {}),
    "f": ("output/Fig2d_epitope_funnel/u10_rise.svg", {}),
}
SIZE = r'font-size(?::\s*|=")([0-9.]+)'
STROKE = r'(stroke-width(?::\s*|="))([0-9.]+)'
ARIAL = "Arial, Helvetica, sans-serif"


def _vb(t):
    return [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]


def _retext(key, t, T, s):
    """Every label -> T (common/retext.py); b also goes to the house face."""
    t = retext(t, T, s, PANELS[key][1])
    if key == "b":                       # matplotlib default -> the house face
        t = re.sub(r"font-family: ?'DejaVu Sans'[^;\"]*", f"font-family: {ARIAL}", t)
    return t


def _floor(t, s):
    lo = LINE_MIN / s
    return re.sub(STROKE, lambda m: m.group(0) if not 0 < float(m.group(2)) < lo
                  else f"{m.group(1)}{lo:.4f}", t)


def _ink(t):
    """The drawn extent, in the file's units (render, trim, convert back)."""
    vb = _vb(t)
    px = 2400
    with tempfile.TemporaryDirectory() as d:
        src, png = Path(d) / "p.svg", Path(d) / "p.png"
        src.write_text(t)
        subprocess.run(["rsvg-convert", "-w", str(px), "-b", "white", src, "-o", png],
                       check=True)
        w, h, x, y = (int(v) for v in subprocess.run(
            ["magick", png, "-fuzz", "1%", "-trim", "-format", "%w %h %X %Y", "info:"],
            check=True, capture_output=True, text=True).stdout.split())
    u = vb[2] / px
    return (vb[0] + (x - 1) * u, vb[1] + (y - 1) * u, (w + 2) * u, (h + 2) * u)


def _body(key, t):
    """Content safe beside other panels: ids prefixed, CSS scoped, metadata gone."""
    t = re.sub(r"<metadata.*?</metadata>|<title>.*?</title>", "", t, flags=re.S)
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    ids = set(re.findall(r'\sid="([^"]+)"', body))
    body = re.sub(r'(\sid=")([^"]+)"', rf'\g<1>{key}-\2"', body)
    body = re.sub(r'(url\(#|href="#)([^)"]+)',
                  lambda m: m.group(1) + (f"{key}-" if m.group(2) in ids else "")
                  + m.group(2), body)

    def scope(m):
        css = re.sub(r'([^{}]+)\{', lambda r: ",".join(
            f"#P{key} {sel.strip()}" for sel in r.group(1).split(",")) + "{", m.group(2))
        return m.group(1) + css + m.group(3)
    return re.sub(r"(<style[^>]*>)(.*?)(</style>)", scope, body, flags=re.S)


class Panel:
    src_of = {}

    def __init__(self, key, T=None, s=None):
        self.key = key
        self.src = Path(self.src_of[key]).read_text()
        self.text = self.src if T is None else _retext(key, self.src, T, s)
        self.ix, self.iy, self.iw, self.ih = _ink(self.text)

    def at(self, s):
        self.s, self.w, self.h = s, self.iw * s, self.ih * s
        return self

    def svg(self, x, y):
        """Placed and clipped to its own ink (a matplotlib figure patch is the
        size of its canvas and would paint over the neighbours)."""
        c = f"clip-{self.key}"
        return (f'<clipPath id="{c}"><rect x="{x:.3f}" y="{y:.3f}" width="{self.w:.3f}" '
                f'height="{self.h:.3f}"/></clipPath><g clip-path="url(#{c})">'
                f'<g id="P{self.key}" transform="translate({x:.3f} {y:.3f}) '
                f'scale({self.s:.6f}) translate({-self.ix:.3f} {-self.iy:.3f})">'
                f"{_body(self.key, self.text)}</g></g>")

    def sizes(self):
        return sorted({round(float(v) * self.s, 2) for v in re.findall(SIZE, self.text)})

    def rule(self):
        return min(float(m.group(2)) for m in re.finditer(STROKE, self.text)
                   if float(m.group(2)) > 0) * self.s


def _fit(key, T, width):
    """Scaled to a width, text re-set to T; the ink moves with the text, so iterate."""
    s = width / Panel(key).iw
    for _ in range(3):
        s = width / Panel(key, T, s).iw
    return Panel(key, T, s).at(s)


def _pair(k1, k2, T=None):
    """Two panels at one height, filling the block."""
    def solve(p1, p2):
        h = (W - GAP_H) / (p1.iw / p1.ih + p2.iw / p2.ih)
        return h / p1.ih, h / p2.ih
    p1, p2 = Panel(k1), Panel(k2)
    s1, s2 = solve(p1, p2)
    if T is not None:
        for _ in range(3):
            p1, p2 = Panel(k1, T, s1), Panel(k2, T, s2)
            s1, s2 = solve(p1, p2)
        p1, p2 = Panel(k1, T, s1), Panel(k2, T, s2)
    return p1.at(s1), p2.at(s2)


def _d_variant(w, h, T, workdir):
    """2d redrawn by its own script into its slot, text at T. Scaled into the
    slot it shrank to ~0.75x, so text re-set to T grew ~1.3x against the bars:
    the longest stack name was clipped and the value labels ran off the axis
    (2026-11-13). bbox_inches='tight' adds the labels to the requested canvas,
    so iterate the request until the ink lands on w x h. Placed at scale 1."""
    workdir.mkdir(exist_ok=True)
    Panel.src_of["d"] = workdir / "gpu_per_design.svg"
    rw, rh = w / MM, h / MM
    for _ in range(4):
        subprocess.run([PYTHON, D_SCRIPT, f"{rw:.3f}", f"{rh:.3f}", f"{T:.4f}", workdir],
                       check=True, capture_output=True)
        p = Panel("d").at(1.0)
        rw, rh = rw * w / p.w, rh * h / p.h
    return p


def _letter(ch, x, y):
    return (f'<text x="{x - LETTER_DX:.2f}" y="{y + LETTER_PT * 0.716:.2f}" '
            f'font-family="{ARIAL}" font-size="{LETTER_PT}" '
            f'font-weight="700" fill="{INK}">{ch}</text>')


def _page(rows, name):
    height = sum(max(p.h for p in r) for r in rows) + GAP_V * (len(rows) - 1)
    fit = min(1.0, (PAGE_H - 2 * TOP) / height)
    x0 = (PAGE_W - W * fit) / 2
    out, y = [], TOP
    for r in rows:
        x = x0
        for p in r:
            p.at(p.s * fit)
            # 2026-11-18 (Leo): capitals, as the manuscript cites them
            # (Fig. 1A, 2B, 3D ...) and as the journals set them. The
            # panel KEYS stay lowercase - they name files and dict entries.
            out += [p.svg(x, y), _letter(p.key.upper(), x, y)]
            x += p.w + GAP_H * fit
        y += max(p.h for p in r) + GAP_V * fit
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'xmlns:xlink="http://www.w3.org/1999/xlink" width="210mm" height="297mm" '
           f'viewBox="0 0 {PAGE_W:.3f} {PAGE_H:.3f}"><title>{name}</title>'
           f'<rect width="100%" height="100%" fill="#ffffff"/>' + "".join(out) + "</svg>")
    return svg, fit, y - GAP_V * fit


def v1():
    a = Panel("a").at(1); a.at(W / a.iw)
    b = Panel("b").at(1); b.at(W / b.iw)
    c = Panel("c").at(1.0)               # matplotlib: 1 unit = 1 pt as drawn
    d = Panel("d").at(1); d.at(W / d.iw)
    return [[a], [b], [c], [d]]          # e, f removed 2026-11-13


def v2(T, workdir):
    a = _fit("a", T, W)
    b = _fit("b", T, W)
    # 2026-11-17 (Leo): "you can remove panel d now. its data is in the panel c
    # now". d was GPU-hours per million designs; c now carries Designs/GPU-h as
    # its second sub-panel, which is the same comparison the other way up. c
    # takes the full block, and the diffusion/hallucination split that d showed
    # in its hatched legend moved into c as brackets on the stack labels.
    c = _fit("c", T, W)
    # d (2026-11-18): the binder-density row, full block. Its renders are NOT
    # grown to fill the width - the render is taller than it is wide, so a
    # full-width set of four overflows the page and _page would shrink
    # everything below T. They sit at a fixed width (IMG_W in
    # code/Fig2d_binder_density/make_figure.py) with the slack as gaps.
    d = _fit("d", T, W)
    rows = [[a], [b], [c], [d]]          # e, f removed 2026-11-13
    for p in sum(rows, []):              # last, once every scale is final
        p.text = _floor(p.text, p.s)
        p.ix, p.iy, p.iw, p.ih = _ink(p.text)
        p.at(p.s)
    return rows


def _report(name, rows, fit, bottom, T=None):
    print(f"\n{name}: scale-to-page {fit:.3f}, block {W * fit / MM:.1f} mm, "
          f"figure ends {bottom / MM:.1f} mm down a 297 mm page")
    for r in rows:
        for p in r:
            print(f"  {p.key}  {p.w / MM:6.1f} x {p.h / MM:5.1f} mm   type "
                  + ", ".join(f"{v:g}" for v in p.sizes())
                  + f" pt   thinnest rule {p.rule():.2f} pt")
            if T is not None:
                ok = allowed(T, PANELS[p.key][1])
                assert all(any(abs(v - w) < 0.02 for w in ok) for v in p.sizes()), f"{p.key} off T"
                assert p.rule() >= LINE_MIN - 0.005, f"{p.key} rule {p.rule():.3f}"


def _b_svg(workdir):
    """2b's own script, re-run into workdir for an SVG with live text."""
    workdir.mkdir(exist_ok=True)
    svg = workdir / "map_composite_paired.svg"
    subprocess.run([PYTHON, B_SCRIPT, workdir], check=True, capture_output=True)
    return svg


def main(outdir):
    outdir.mkdir(exist_ok=True)
    Panel.src_of = {k: REPO / f for k, (f, _) in PANELS.items() if f}
    Panel.src_of["b"] = _b_svg(outdir / "_b_variant")
    rows1 = v1()
    svg1, fit1, bot1 = _page(rows1, "Fig2 V1 as built")
    rows2 = v2(T_FIG1, outdir / "_d_variant")
    svg2, fit2, bot2 = _page(rows2, "Fig2 V2 adjusted")
    assert fit2 == 1.0, "V2 must fit without shrinking, or its type is no longer T"
    _report("V1 as built", rows1, fit1, bot1)
    _report(f"V2 adjusted, T = {T_FIG1:.2f} pt", rows2, fit2, bot2, T_FIG1)
    for stem, svg in (("Fig2_V1_as_built", svg1), ("Fig2_V2_adjusted", svg2)):
        p = outdir / f"{stem}.svg"
        p.write_text(svg, encoding="utf-8")
        subprocess.run(["rsvg-convert", "-f", "pdf", p, "-o", p.with_suffix(".pdf")], check=True)
        subprocess.run(["rsvg-convert", "--dpi-x", "200", "--dpi-y", "200", p,
                        "-o", p.with_suffix(".png")], check=True)
        print("wrote", outdir / stem, "{svg,pdf,png}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "working" / "Fig2_assembly")
