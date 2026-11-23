"""Panels onto an A4 page, for the thesis figure assemblers.

The generic half of code/Fig1_assembly and code/Fig2_assembly (which still
carry their own copies): read a panel's SVG, crop it to its drawn ink, make it
safe beside other panels (ids prefixed, CSS scoped), re-set its text to the
house size, raise thin strokes to a floor, and stack rows on an A4 page with
the panel letters outside.

Page: A4, 16 mm sides = the 178 mm block in config/style.yaml, 15 mm top.
"""
import re
import subprocess
import tempfile
from pathlib import Path

from common import palette, style
from common.retext import allowed, retext

MM = 72 / 25.4
PAGE_W, PAGE_H = 210 * MM, 297 * MM
SIDE, TOP = 16 * MM, 15 * MM
W = PAGE_W - 2 * SIDE                    # 504.567 pt, the 178 mm block
GAP_V, GAP_H = 5 * MM, 6 * MM
LETTER_PT, LETTER_DX = float(style.CFG["title_size_pt"]), 4.5 * MM
T_HOUSE = float(style.CFG["font_size_pt"])
LINE_MIN = 0.5                           # pt: the strictest common journal floor
ARIAL = "Arial, Helvetica, sans-serif"
SIZE = r'font-size(?::\s*|=")([0-9.]+)'
STROKE = r'(stroke-width(?::\s*|="))([0-9.]+)'


def vb(t):
    return [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]


def ink(t):
    """The drawn extent, in the file's units (render, trim, convert back)."""
    v = vb(t)
    px = 2400
    with tempfile.TemporaryDirectory() as d:
        src, png = Path(d) / "p.svg", Path(d) / "p.png"
        src.write_text(t)
        subprocess.run(["rsvg-convert", "-w", str(px), "-b", "white", src, "-o", png], check=True)
        w, h, x, y = (int(q) for q in subprocess.run(
            ["magick", png, "-fuzz", "1%", "-trim", "-format", "%w %h %X %Y", "info:"],
            check=True, capture_output=True, text=True).stdout.split())
    u = v[2] / px
    return (v[0] + (x - 1) * u, v[1] + (y - 1) * u, (w + 2) * u, (h + 2) * u)


def floor(t, s, skip=()):
    """Strokes that would print under LINE_MIN at scale s raised to it; CSS
    classes named in skip are left alone (e.g. ladder rungs that would fuse)."""
    lo = LINE_MIN / s
    fix = lambda u: re.sub(STROKE, lambda m: m.group(0) if not 0 < float(m.group(2)) < lo
                           else f"{m.group(1)}{lo:.4f}", u)
    if not skip:
        return fix(t)
    parts = re.split(r"(\.(?:%s)\{[^}]*\})" % "|".join(map(re.escape, skip)), t)
    return "".join(p if i % 2 else fix(p) for i, p in enumerate(parts))


def body(key, t):
    """Content safe beside other panels: ids prefixed, CSS scoped, metadata gone."""
    t = re.sub(r"<metadata.*?</metadata>|<title>.*?</title>", "", t, flags=re.S)
    b = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    ids = set(re.findall(r'\sid="([^"]+)"', b))
    b = re.sub(r'(\sid=")([^"]+)"', rf'\g<1>{key}-\2"', b)
    b = re.sub(r'(url\(#|href="#)([^)"]+)',
               lambda m: m.group(1) + (f"{key}-" if m.group(2) in ids else "") + m.group(2), b)

    def scope(m):
        css = re.sub(r'([^{}]+)\{', lambda r: ",".join(
            f"#P{key} {sel.strip()}" for sel in r.group(1).split(",")) + "{", m.group(2))
        return m.group(1) + css + m.group(3)
    return re.sub(r"(<style[^>]*>)(.*?)(</style>)", scope, b, flags=re.S)


class Panel:
    """One panel file. With T and s, its text is re-set to T (riders aside).
    canvas_x: place by the canvas's width, not the ink's - for a panel whose
    empty side is meant (a column held for something still to come)."""

    def __init__(self, key, src, T=None, s=None, riders=None, canvas_x=False, prep=None):
        self.key, self.riders = key, riders or {}
        self.src = Path(src).read_text()
        t = self.src if T is None else retext(self.src, T, s, self.riders)
        self.canvas_x = canvas_x
        self.retouch(lambda u: prep(u) if prep else u, t)

    def retouch(self, fn, t=None):
        """Change the text (e.g. floor its strokes) and re-measure its extent."""
        self.text = fn(self.text if t is None else t)
        self.ix, self.iy, self.iw, self.ih = ink(self.text)
        if self.canvas_x:
            v = vb(self.text)
            self.ix, self.iw = v[0], v[2]
        if hasattr(self, "s"):
            self.at(self.s)
        return self

    def at(self, s):
        self.s, self.w, self.h = s, self.iw * s, self.ih * s
        return self

    def svg(self, x, y):
        c = f"clip-{self.key}"
        return (f'<clipPath id="{c}"><rect x="{x:.3f}" y="{y:.3f}" width="{self.w:.3f}" '
                f'height="{self.h:.3f}"/></clipPath><g clip-path="url(#{c})">'
                f'<g id="P{self.key}" transform="translate({x:.3f} {y:.3f}) '
                f'scale({self.s:.6f}) translate({-self.ix:.3f} {-self.iy:.3f})">'
                f"{body(self.key, self.text)}</g></g>")

    def sizes(self):
        return sorted({round(float(v) * self.s, 2) for v in re.findall(SIZE, self.text)})

    def rule(self, skip=()):
        """The thinnest stroke as printed, or None for a panel that has none -
        an embedded raster, say. min() over nothing used to raise."""
        t = re.sub(r"\.(?:%s)\{[^}]*\}" % "|".join(map(re.escape, skip)), "", self.text) if skip else self.text
        w = [float(m.group(2)) for m in re.finditer(STROKE, t) if float(m.group(2)) > 0]
        return min(w) * self.s if w else None


def fit(make, T, width=None, height=None):
    """make(T, s) -> Panel; scaled to width (or to height) with text at T.
    Re-set text moves the ink, so iterate the scale."""
    goal, dim = (width, lambda q: q.iw) if height is None else (height, lambda q: q.ih)
    s = goal / dim(make(None, None))
    for _ in range(3):
        s = goal / dim(make(T, s))
    return make(T, s).at(s)


def pair(make1, make2, T=None, width=None):
    """Two panels side by side at one height, filling width (the block by
    default). make(T, s) -> Panel. Text re-set to T moves the ink: iterate."""
    width = width or W

    def solve(p1, p2):
        h = (width - GAP_H) / (p1.iw / p1.ih + p2.iw / p2.ih)
        return h / p1.ih, h / p2.ih
    p1, p2 = make1(None, None), make2(None, None)
    s1, s2 = solve(p1, p2)
    if T is not None:
        for _ in range(3):
            p1, p2 = make1(T, s1), make2(T, s2)
            s1, s2 = solve(p1, p2)
        p1, p2 = make1(T, s1), make2(T, s2)
    return p1.at(s1), p2.at(s2)


def variant(key, script, svg_name, w, h, T, workdir, python):
    """A matplotlib panel redrawn by its own script to w x h pt with its text
    at T (script CLI: W_MM H_MM FONT_PT OUTDIR). It saves with
    bbox_inches='tight', so labels add to the canvas: iterate the request until
    the ink lands on w x h. Placed at scale 1, so its type is exactly T."""
    workdir.mkdir(exist_ok=True)
    rw, rh = w / MM, h / MM
    for _ in range(4):
        subprocess.run([python, script, f"{rw:.3f}", f"{rh:.3f}", f"{T:.4f}", workdir],
                       check=True, capture_output=True, cwd=Path(script).parents[1])
        p = Panel(key, workdir / svg_name).at(1.0)
        rw, rh = rw * w / p.w, rh * h / p.h
    return p


def window(src, box, out):
    """A copy of an SVG showing only box = (x, y, w, h) of it, in its units -
    for taking one sub-plot out of a file that holds several."""
    t = Path(src).read_text()
    x, y, w, h = box
    t = re.sub(r'width="[^"]*" height="[^"]*" viewBox="[^"]*"',
               f'width="{w:.3f}pt" height="{h:.3f}pt" viewBox="{x:.3f} {y:.3f} {w:.3f} {h:.3f}"',
               t, count=1)
    Path(out).write_text(t, encoding="utf-8")
    return Path(out)


def letter(ch, x, y):
    return (f'<text x="{x - LETTER_DX:.2f}" y="{y + LETTER_PT * 0.716:.2f}" '
            f'font-family="{ARIAL}" font-size="{LETTER_PT}" font-weight="700" '
            f'fill="{palette.INK}">{ch}</text>')


def page(rows, name, letters=None):
    """Rows stacked from the top margin, shrunk as a whole only if they overflow."""
    height = sum(max(p.h for p in r) for r in rows) + GAP_V * (len(rows) - 1)
    fit_ = min(1.0, (PAGE_H - 2 * TOP) / height)
    x0 = (PAGE_W - W * fit_) / 2
    out, y = [], TOP
    for r in rows:
        x = x0
        for p in r:
            p.at(p.s * fit_)
            # 2026-11-18: p.dx widens the gap before a panel, in the page's own
            # units, for a row whose panels do not add up to the block - Fig 5's
            # a | b, where b is sized by HEIGHT to match a and so leaves slack.
            # Without it the row packs left and the block looks short against
            # the full-width rows below.
            x += getattr(p, "dx", 0.0) * fit_
            # 2026-11-18 (Leo): capitals, as the manuscript cites them
            # (Fig. 1A, 2B, 3D ...) and as the journals set them. The
            # panel KEYS stay lowercase - they name files and dict entries.
            out += [p.svg(x, y),
                    letter((letters or {}).get(p.key, p.key.upper()), x, y)]
            x += p.w + GAP_H * fit_
        y += max(p.h for p in r) + GAP_V * fit_
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
           f'width="210mm" height="297mm" viewBox="0 0 {PAGE_W:.3f} {PAGE_H:.3f}">'
           f'<title>{name}</title><rect width="100%" height="100%" fill="#ffffff"/>'
           + "".join(out) + "</svg>")
    return svg, fit_, y - GAP_V * fit_


def report(name, rows, fit_, bottom, T=None, skip=None):
    """Print sizes; with T, assert every label is T (riders aside) and no rule
    is under LINE_MIN."""
    print(f"\n{name}: scale-to-page {fit_:.3f}, block {W * fit_ / MM:.1f} mm, "
          f"figure ends {bottom / MM:.1f} mm down a 297 mm page")
    for r in rows:
        for p in r:
            sk = (skip or {}).get(p.key, ())
            r = p.rule(sk)
            print(f"  {p.key}  {p.w / MM:6.1f} x {p.h / MM:5.1f} mm   type "
                  + (", ".join(f"{v:g}" for v in p.sizes()) or "none (raster)")
                  + ("   thinnest rule none" if r is None else f" pt   thinnest rule {r:.2f} pt"))
            if T is not None:
                ok = allowed(T, p.riders)
                stray = {v for v in p.sizes() if not any(abs(v - w) <= 0.02 for w in ok)}
                assert not stray, f"{p.key}: stray sizes {stray}"
                assert r is None or r >= LINE_MIN - 0.005, f"{p.key}: rule {r:.3f} pt"


def write(outdir, stem, svg):
    p = outdir / f"{stem}.svg"
    p.write_text(svg, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-f", "pdf", p, "-o", p.with_suffix(".pdf")], check=True)
    subprocess.run(["rsvg-convert", "--dpi-x", "200", "--dpi-y", "200", p,
                    "-o", p.with_suffix(".png")], check=True)
    print("wrote", outdir / stem, "{svg,pdf,png}")
