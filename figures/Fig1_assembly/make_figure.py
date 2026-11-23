"""Fig 1, assembled on an A4 page from each panel's newest file.

    python code/Fig1_assembly/make_figure.py [outdir]     # default working/Fig1_assembly

Two versions of the same six panels; nothing is redrawn.

V1 "as_sketched"  Leo's arrangement: a / b c / d / e / f. Each panel is scaled
                  as a whole to fill its slot and its own text is left alone,
                  so every panel prints at a different type size. If the stack
                  is taller than the page it is shrunk as a whole.

V2 "adjusted"     every caption in every panel at ONE size, T. T is what Fig 1d
                  prints at when it fills the text block, so T is set by d, not
                  chosen. V2 uses d's narrowed copy ('Full Panel V3 compact.svg',
                  code/Fig1d_integrated_construct/compact_full_panel.py), which
                  is what lifts T over 5 pt.
                  b and d are drawings built around their labels (block labels
                  sit inside blocks), so they keep the scale that puts their
                  label tier at T. a, c, e and f are scaled to their slot and
                  their text is re-set to T in their own units. e moves up
                  beside f so the page fits - and since e's labels sit at
                  fixed offsets that collide when its text grows, e is
                  regenerated at its slot size with T as its font instead
                  (code/Fig1e_assay_bridge/make_panel.py W H PT OUTDIR).

Page: A4, 16 mm sides = the 178 mm block in config/style.yaml, 15 mm top.
Panel letters: 8 pt bold (title_size_pt), outside the panel, cap-top level
with the panel's ink (STYLE_NOTES section 1).
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
GAP_V = 5 * MM                           # between rows
GAP_H = 6 * MM                           # between panels in a row: room for a letter
LETTER_PT, LETTER_DX = 8.0, 4.5 * MM
INK = palette.INK

# file, riders. In V2 every label is set to T except the riders, named by
# their size in the file: {size: multiple of T}.
PANELS = {
    "a": ("output/Fig1a_screen_pipeline/Fig1a_screen_pipeline.svg", {}),
    "b": ("output/Fig1b_sms_assay/Fluorescence V1.svg",
          {12.0: 12.0 / 9.997}),                        # '*' stop codon rides with 5'/3'
    "c": ("output/Fig1c_sms_titration/sms_titration.svg", {}),
    # 2026-11-17 (Leo): scaffold optimisation moves UP, so the concept is
    # introduced before the integrated construct. The letters are the page
    # order, so the sources were re-keyed rather than the rows re-sorted:
    #   d = the scaffold workflow   (was the right-hand half of old f)
    #   e = the integrated construct (was d)
    #   f = both readouts, now carrying the scaffold bars too (was e)
    "d": ("output/Fig1f_scaffold_splicing/scaffold_workflow.svg", {}),
    "e": ("selection/Fig1d_integrated_construct.svg",
          {34.995: 34.995 / 27.496}),                   # 'Splicing ratio' heading
    "f": ("selection/Fig1e_assay_bridge.svg", {}),
}
LABEL = {"b": 8.0, "e": 27.496}          # the tier a drawing was built around

# V1 only: e and f as fractions of the block, read off Leo's sketch
E_FRAC, F_FRAC = 0.36, 0.76


# V2: e is regenerated and gets the width its labels need; f is scaled into the
# rest, so its text grows against its own drawing. Its tightest gap - the
# 'Circular scaffold' label against the wet-lab arrow - sets how far that can go.
# f (the readouts) is regenerated full width at the foot of the page. 38 mm is
# what the page budget leaves once a, b/c and the size-locked construct are in;
# six bars across 178 mm are comfortable at that height.
F_H_MM = 38.0
# d (the workflow) takes the height still unspent, at its own ink aspect, rather
# than a fixed size: it is the only panel on the page free to give.
D_H_MIN = 26.0         # below this its three-line wet-lab label crowds the cell glyph
LINE_MIN = 0.5         # pt, V2 only: the strictest common journal floor (BMC,
                       # Elsevier). Science's own figure page could not be read.
                       # Thinner strokes are raised to it; thicker ones are kept.
# Exempt from the floor: d's mechanism-inset base-pair rungs (light grey, 16 of
# them) sit 0.73 pt apart, so at 0.5 pt the ladder fuses into a solid bar.
FLOOR_SKIP = {"e": ("cls-38",)}
T_TARGET = float(style.CFG["font_size_pt"])   # pt, every label in V2: the house
                       # size, 5.5 since 2026-11-13. The 5.17 pt version this
                       # replaced is kept in working/Fig1_assembly/_5p17/
D_V2 = REPO / "working" / "Fig1d_integrated_construct" / "Full Panel V4 compact.svg"
E_SCRIPT = REPO / "code" / "Fig1e_assay_bridge" / "make_panel.py"
C_SCRIPT = REPO / "code" / "Fig1c_sms_titration" / "make_panel.py"
PYTHON = sys.executable
SIZE = r'font-size(?::\s*|=")([0-9.]+)'


def _vb(t):
    return [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]


def _retext(key, t, T, s):
    """Every label -> T, riders aside (common/retext.py)."""
    return retext(t, T, s, PANELS[key][1])


STROKE = r'(stroke-width(?::\s*|="))([0-9.]+)'


def _split(t, skip):
    """t cut into (text, exempt?) pieces: the CSS rules named in skip are exempt."""
    if not skip:
        return [(t, False)]
    parts = re.split(r"(\.(?:%s)\{[^}]*\})" % "|".join(map(re.escape, skip)), t)
    return [(p, i % 2 == 1) for i, p in enumerate(parts)]


def _floor(t, s, skip=()):
    """Every stroke that would print under LINE_MIN at scale s, raised to it."""
    lo = LINE_MIN / s
    fix = lambda u: re.sub(STROKE, lambda m: m.group(0) if not 0 < float(m.group(2)) < lo
                           else f"{m.group(1)}{lo:.4f}", u)
    return "".join(p if exempt else fix(p) for p, exempt in _split(t, skip))


def _thinnest(t, skip=()):
    return min(float(m.group(2)) for p, exempt in _split(t, skip) if not exempt
               for m in re.finditer(STROKE, p) if float(m.group(2)) > 0)


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
    u = vb[2] / px                       # units per pixel
    return (vb[0] + (x - 1) * u, vb[1] + (y - 1) * u, (w + 2) * u, (h + 2) * u)


def _body(key, t):
    """The file's content, safe to sit beside five others: ids prefixed, CSS
    scoped to this panel's group, metadata (with its rdf: prefix) dropped."""
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
    def __init__(self, key, T=None, s=None, src=None):
        self.key = key
        self.src = Path(src or REPO / PANELS[key][0]).read_text()
        self.text = self.src if T is None else _retext(key, self.src, T, s)
        self.ix, self.iy, self.iw, self.ih = _ink(self.text)

    def at(self, s):
        self.s, self.w, self.h = s, self.iw * s, self.ih * s
        return self

    def svg(self, x, y):
        """Placed at (x, y) and clipped to its own ink. A matplotlib file
        carries a white figure patch the size of its CANVAS, which reaches past
        its ink and paints over whatever sits next to it."""
        c = f"clip-{self.key}"
        return (f'<clipPath id="{c}"><rect x="{x:.3f}" y="{y:.3f}" width="{self.w:.3f}" '
                f'height="{self.h:.3f}"/></clipPath><g clip-path="url(#{c})">'
                f'<g id="P{self.key}" transform="translate({x:.3f} {y:.3f}) '
                f'scale({self.s:.6f}) translate({-self.ix:.3f} {-self.iy:.3f})">'
                f"{_body(self.key, self.text)}</g></g>")

    def sizes(self, text=None):
        """Every font size in the placed panel, in printed pt."""
        return {round(float(v) * self.s, 2) for v in re.findall(SIZE, text or self.text)}


def _slot(key, T, width_of):
    """A panel scaled to a slot whose width depends on its own ink, with its
    text re-set to T. Re-setting the text moves the ink, so iterate."""
    s = width_of(Panel(key))
    for _ in range(3):
        s = width_of(Panel(key, T, s))
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


def _f_variant(w, h, T, workdir):
    """The readout panel redrawn by its own script at T, ink landing on w x h pt.
    The script saves with bbox_inches='tight', so labels add to the requested
    canvas; iterate the request until the ink fits. Placed at scale 1, so its
    type is exactly T."""
    workdir.mkdir(exist_ok=True)
    rw, rh = w / MM, h / MM
    for _ in range(4):
        subprocess.run([PYTHON, E_SCRIPT, f"{rw:.3f}", f"{rh:.3f}", f"{T:.4f}", workdir],
                       check=True, capture_output=True)
        p = Panel("f", src=workdir / "Fig1e_assay_bridge.svg").at(1.0)
        rw, rh = rw * w / p.w, rh * h / p.h
    return p


def _c_variant(w, T, workdir):
    """Fig 1c redrawn by its own script with every label at one size, then
    scaled to width w. Its geometry is fixed in mm, so the label size the
    script needs is T / scale; matplotlib then places its y-label clear of
    the ticks, which re-setting text in place cannot do. Iterate: the labels
    hang outside the grid, so the ink - and the scale - move with them."""
    workdir.mkdir(exist_ok=True)
    s = w / Panel("c").iw
    for _ in range(4):
        subprocess.run([PYTHON, C_SCRIPT, f"{T / s:.4f}", workdir],
                       check=True, capture_output=True)
        p = Panel("c", src=workdir / "sms_titration.svg")
        s = w / p.iw
    return p.at(s)


def _letter(ch, x, y):
    return (f'<text x="{x - LETTER_DX:.2f}" y="{y + LETTER_PT * 0.716:.2f}" '
            f'font-family="Arial, Helvetica, sans-serif" font-size="{LETTER_PT}" '
            f'font-weight="700" fill="{INK}">{ch}</text>')


def _page(rows, name):
    """Stack rows top-down from the top margin; shrink the lot if it overflows."""
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
    a = Panel("a").at(0); a.at(W / a.iw)
    b, c = _pair("b", "c")
    d = Panel("d").at(0); d.at(W / d.iw)
    e = Panel("e").at(0); e.at(W / e.iw)
    f = Panel("f").at(0); f.at(F_FRAC * W / f.iw)
    d = Panel("d").at(0); d.at(E_FRAC * W / d.iw)
    return [[a], [b, c], [d], [e], [f]]


def v2(workdir):
    T = T_TARGET
    e = Panel("e", src=D_V2).at(T / LABEL["e"])     # construct: scale pinned by its type tier
    assert e.w <= W + 0.01, f"e is {e.w / MM:.1f} mm at {T} pt: narrow it further"
    b = Panel("b", T, T / LABEL["b"]).at(T / LABEL["b"])
    a = _slot("a", T, lambda p: W / p.iw)
    c = _c_variant(W - GAP_H - b.w, T, workdir.parent / "_c_variant")
    f = _f_variant(W, F_H_MM * MM, T, workdir)
    # The construct's width is locked by its text tier, so it is the one panel
    # that cannot give. Everything else is fixed too, which leaves exactly one
    # free quantity: the workflow's height. Spend what is left, no more.
    spent = a.h + max(b.h, c.h) + e.h + f.h + 4 * GAP_V
    d_h = (PAGE_H - 2 * TOP) - spent
    print(f"d gets {d_h / MM:.1f} mm of height "
          f"({(PAGE_H - 2 * TOP - spent) / MM:.1f} mm unspent before it)")
    assert d_h >= D_H_MIN * MM, (
        f"only {d_h / MM:.1f} mm left for the workflow (floor {D_H_MIN} mm). "
        "The page is full: shrink another panel or move one to supplementary.")
    d0 = Panel("d")
    s_d = d_h / d0.ih
    for _ in range(2):                   # re-set text moves the ink, so iterate
        s_d = d_h / Panel("d", T, s_d).ih
    # Width can bind before height does, and when it does the panel is PLACED
    # smaller than the scale its text was re-set for, so its type prints under T
    # and _report's stray-size assert fires. Settle the final scale first, then
    # re-set the text for that scale (2026-11-17: a wide, short workflow made
    # width the binding constraint for the first time).
    s_d = min(s_d, W / d0.iw)
    d = Panel("d", T, s_d).at(s_d)
    rows = [[a], [b, c], [d], [e], [f]]
    for k, skip in FLOOR_SKIP.items():   # the exemption must still name the rungs
        p = {q.key: q for q in sum(rows, [])}[k]
        for c in skip:
            assert re.search(rf"\.{c}\{{[^}}]*stroke:#b3b8bc", p.text), f"{k}: {c} is not the rungs"
    for p in sum(rows, []):              # last, once every scale is final
        p.text = _floor(p.text, p.s, FLOOR_SKIP.get(p.key, ()))
        p.ix, p.iy, p.iw, p.ih = _ink(p.text)
        p.at(p.s)
    return rows, T


def _report(name, rows, fit, bottom, T=None):
    print(f"\n{name}: scale-to-page {fit:.3f}, block {W * fit / MM:.1f} mm, "
          f"figure ends {bottom / MM:.1f} mm down a 297 mm page")
    for r in rows:
        for p in r:
            rule = _thinnest(p.text, FLOOR_SKIP.get(p.key, ()) if T else ())
            print(f"  {p.key}  {p.w / MM:6.1f} x {p.h / MM:5.1f} mm   type "
                  + ", ".join(f"{v:g}" for v in sorted(p.sizes()))
                  + f" pt   thinnest rule {rule * p.s:.2f} pt")
            if T is None:
                continue
            assert rule * p.s >= LINE_MIN - 0.005, f"{p.key}: rule {rule * p.s:.3f} pt"
            ok = allowed(T, PANELS[p.key][1])
            stray = {v for v in p.sizes() if not any(abs(v - w) <= 0.02 for w in ok)}
            assert not stray, f"{p.key}: stray sizes {stray}"


def main(outdir):
    outdir.mkdir(exist_ok=True)
    rows1 = v1()
    svg1, fit1, bot1 = _page(rows1, "Fig1 V1 as sketched")
    rows2, T = v2(outdir / "_f_variant")   # c's lands in _c_variant
    svg2, fit2, bot2 = _page(rows2, "Fig1 V2 adjusted")
    assert fit2 == 1.0, "V2 must fit without shrinking, or its type is no longer T"
    _report("V1 as sketched", rows1, fit1, bot1)
    _report(f"V2 adjusted, T = {T:.2f} pt", rows2, fit2, bot2, T)
    for stem, svg in (("Fig1_V1_as_sketched", svg1), ("Fig1_V2_adjusted", svg2)):
        p = outdir / f"{stem}.svg"
        p.write_text(svg, encoding="utf-8")
        subprocess.run(["rsvg-convert", "-f", "pdf", p, "-o", p.with_suffix(".pdf")], check=True)
        subprocess.run(["rsvg-convert", "--dpi-x", "200", "--dpi-y", "200", p,
                        "-o", p.with_suffix(".png")], check=True)
        print("wrote", outdir / stem, "{svg,pdf,png}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "working" / "Fig1_assembly")
