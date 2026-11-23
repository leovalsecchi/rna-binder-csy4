"""All five thesis figures in one SVG, added to as each figure is settled.

    python code/Fig_book/make_book.py               # rebuild from what is in selection/
    python code/Fig_book/make_book.py --promote 3   # copy Fig 3's assembled page in, then rebuild
    python code/Fig_book/make_book.py --promote 3:Fig3_V1_as_sketched   # a named version instead

The book is a function of selection/: page N is selection/FigN.svg when it
exists, otherwise an empty A4 page carrying its title. Promoting is the
deliberate step - a figure joins the book when its assembled page
(working/FigN_assembly/FigN_V2_adjusted.*) is copied into selection/ as FigN.* -
so re-running an assembler never changes the book by itself.

Writes to selection/ (README: things arrive there freely; final/ stays manual):
  Thesis_figures.svg  the five A4 pages side by side, each fully embedded with
                      its ids prefixed so the pages cannot repaint each other
  Thesis_figures.pdf  one A4 page per figure, in order
  Thesis_figures.png  preview of the SVG
Sits in code/Fig_book/ so run_all.sh, which takes Fig*/ in order, builds it last.
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette, style                            # noqa: E402

SEL = REPO / "selection"
N_FIGS = 5
MM = 72 / 25.4
PW, PH = 210 * MM, 297 * MM            # A4, in pt
GAP, MARGIN, BAND = 12 * MM, 12 * MM, 9 * MM   # between pages, round the sheet, title band
FONT = f'{style.CFG["font_family"]}, Helvetica, sans-serif'
DESK, EDGE, INK, MUT = (palette.NEUTRAL["300"], palette.NEUTRAL["400"],
                        palette.INK, palette.NEUTRAL["600"])


def titles():
    """Figure titles from FIGURES.md, '## Fig N — title'."""
    t = (REPO / "FIGURES.md").read_text()
    return {int(n): s.strip() for n, s in re.findall(r"^## Fig (\d+) — (.+)$", t, re.M)}


def promote(spec):
    """'N' promotes FigN_V2_adjusted; 'N:STEM' promotes working/FigN_assembly/STEM."""
    n, _, stem = spec.partition(":")
    stem = stem or f"Fig{n}_V2_adjusted"
    src = REPO / "working" / f"Fig{n}_assembly"
    for ext in ("svg", "pdf", "png"):
        f = src / f"{stem}.{ext}"
        if not f.exists():
            raise SystemExit(f"no assembled page to promote: {f}")
        shutil.copy2(f, SEL / f"Fig{n}.{ext}")
    print(f"promoted Fig {n}: {src.name}/{stem}.* -> selection/Fig{n}.*")


def _vb(t):
    return [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]


def _embed(n, t):
    """A page's content with every id prefixed - in id=, url(#), href="#" and
    the #id selectors of its style blocks - so five pages can share one file."""
    t = re.sub(r"<metadata.*?</metadata>|<title>.*?</title>", "", t, flags=re.S)
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    ids = set(re.findall(r'\sid="([^"]+)"', body))
    p = f"f{n}-"
    body = re.sub(r'(\sid=")([^"]+)"', lambda m: f'{m.group(1)}{p}{m.group(2)}"', body)
    body = re.sub(r'(url\(#|href="#)([^)"]+)',
                  lambda m: m.group(1) + (p if m.group(2) in ids else "") + m.group(2), body)
    body = re.sub(r"(<style[^>]*>)(.*?)(</style>)", lambda m: m.group(1) + re.sub(
        r"#([A-Za-z_][\w-]*)", lambda s: "#" + (p if s.group(1) in ids else "") + s.group(1),
        m.group(2)) + m.group(3), body, flags=re.S)
    return body


def placeholder(n, title):
    """An empty A4 page holding figure n's place, as its own SVG."""
    cx = PW / 2
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="210mm" height="297mm" '
            f'viewBox="0 0 {PW:.3f} {PH:.3f}"><rect width="100%" height="100%" fill="#ffffff"/>'
            f'<text x="{cx:.1f}" y="{PH / 2 - 8:.1f}" font-family="{FONT}" font-size="14" '
            f'font-weight="700" fill="{MUT}" text-anchor="middle">Figure {n}</text>'
            f'<text x="{cx:.1f}" y="{PH / 2 + 10:.1f}" font-family="{FONT}" font-size="9" '
            f'fill="{MUT}" text-anchor="middle">{title}</text>'
            f'<text x="{cx:.1f}" y="{PH / 2 + 26:.1f}" font-family="{FONT}" font-size="7" '
            f'fill="{EDGE}" text-anchor="middle">not yet in selection/ - '
            f'make_book.py --promote {n}</text></svg>')


def pages():
    """[(n, title, page svg text, promoted?)] for figures 1..N_FIGS."""
    tt = titles()
    out = []
    for n in range(1, N_FIGS + 1):
        f = SEL / f"Fig{n}.svg"
        title = tt.get(n, "")
        out.append((n, title, f.read_text() if f.exists() else placeholder(n, title), f.exists()))
    return out


def book(pp):
    W = 2 * MARGIN + N_FIGS * PW + (N_FIGS - 1) * GAP
    H = 2 * MARGIN + BAND + PH
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
         f'width="{W / MM:.1f}mm" height="{H / MM:.1f}mm" viewBox="0 0 {W:.3f} {H:.3f}">'
         f'<title>Thesis figures</title><rect width="100%" height="100%" fill="{DESK}"/>']
    for i, (n, title, t, done) in enumerate(pp):
        x, y = MARGIN + i * (PW + GAP), MARGIN + BAND
        vb = _vb(t)
        k = PW / vb[2]
        head = f"Figure {n}" + (f" — {title}" if title else "")
        o.append(f'<text x="{x:.2f}" y="{MARGIN + BAND * 0.62:.2f}" font-family="{FONT}" '
                 f'font-size="11" font-weight="{700 if done else 400}" '
                 f'fill="{INK if done else MUT}">{head}</text>')
        o.append(f'<clipPath id="page{n}"><rect x="{x:.3f}" y="{y:.3f}" width="{PW:.3f}" '
                 f'height="{PH:.3f}"/></clipPath>'
                 f'<g id="F{n}" clip-path="url(#page{n})"><rect x="{x:.3f}" y="{y:.3f}" '
                 f'width="{PW:.3f}" height="{PH:.3f}" fill="#ffffff"/>'
                 f'<g transform="translate({x:.3f} {y:.3f}) scale({k:.6f}) '
                 f'translate({-vb[0]:.3f} {-vb[1]:.3f})">{_embed(n, t)}</g></g>'
                 f'<rect x="{x:.3f}" y="{y:.3f}" width="{PW:.3f}" height="{PH:.3f}" '
                 f'fill="none" stroke="{EDGE}" stroke-width="0.6"/>')
    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--promote", action="append", default=[], metavar="N[:STEM]",
                    help="copy working/FigN_assembly/FigN_V2_adjusted.* (or STEM.*) to selection/FigN.*")
    for n in ap.parse_args().promote:
        promote(n)
    pp = pages()
    svg = SEL / "Thesis_figures.svg"
    svg.write_text(book(pp), encoding="utf-8")
    with tempfile.TemporaryDirectory() as d:     # PDF: one A4 page per figure
        files = []
        for n, title, t, done in pp:
            f = Path(d) / f"p{n}.svg"
            f.write_text(t, encoding="utf-8")
            files.append(str(f))
        subprocess.run(["rsvg-convert", "-f", "pdf", *files, "-o", svg.with_suffix(".pdf")],
                       check=True)
    subprocess.run(["rsvg-convert", "-w", "3000", svg, "-o", svg.with_suffix(".png")], check=True)
    have = [n for n, *_, done in pp if done]
    print(f"wrote selection/Thesis_figures.{{svg,pdf,png}}: figures {have} in, "
          f"{[n for n in range(1, N_FIGS + 1) if n not in have]} still placeholders")


if __name__ == "__main__":
    main()
