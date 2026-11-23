"""Fig 4a - the library cassette before integration.

This IS the promoterless donor plasmid from the cassette-integration figure
(`output/cartoon_library/integration/Promoterless-donor-plasmid.svg`), not a
copy of it and not a redrawing. The file is read, and exactly two things are
done to it:

1. **The plasmid outline is left alone**, byte for byte. Matteo's backbone is a
   rounded rectangle 0.2 to 150.329 wide and 11.78 to 43.78 tall.
2. **Everything else is translated straight up by 32**, which is the distance
   between that rectangle's bottom edge and its top edge. The cassette, which
   sat on the bottom edge, now sits on the top edge with the loop hanging
   below. Nothing is scaled, mirrored or re-spaced: every block keeps its own
   geometry, its width, its colour and its distance from its neighbours, and
   the brackets and label ticks travel with the blocks they belong to.

Then the two additions: the oligo pool dropping into the binder, and the
target scaffold blown up into the four constructs that were screened.

The only deletion is Matteo's "Promoterless donor plasmid" heading, which is a
title and belongs in the caption (STYLE_NOTES section 1).

    python code/Fig4a_library_cassette/make_figure.py
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                    # noqa: E402
from common.paths import ROOT                                 # noqa: E402

DONOR = ROOT / "output/cartoon_library/integration/Promoterless-donor-plasmid.svg"
GLYPHS = ROOT / "output/Fig1f_scaffold_splicing/elements"
OUT = ROOT / "output/Fig4a_library_cassette"

INK, OUTLINE = palette.INK, palette.SYNBIO["outline"]
MUT = palette.NEUTRAL["500"]

# The four screened targets: mod1 (Y array, "A") and mod4 (circular, "B"),
# sites live or dead. Fig1f already emits these, so they are reused as files.
# One flat row rather than a 2x2 block: it keeps the call-out short, which is
# what stops it driving the panel's height. Order matches Fig 4b's 2x2 read
# row-major, so the two panels list the four in the same sequence.
CELLS = ["f30_variation_positive", "circular_scaffold_positive",
         "f30_variation_negative", "circular_scaffold_negative"]
BOX = (81.0, -34.0, 168.0, -6.0)              # the call-out, x0 y0 x1 y1
TITLE_H = 8.0                                 # the title's band inside the box
BOX_R = 2.0                                   # its corner radius
# Arial advance widths, enough for the module captions - used to find where
# those captions actually sit, so a leader can be checked against them.
_ADV = {"S": .667, "e": .556, "n": .556, "s": .5, "o": .556, "r": .333,
        " ": .278, "m": .833, "d": .556, "u": .556, "l": .222, "T": .611,
        "a": .556, "g": .556, "t": .278, "R": .722, "c": .5, "i": .222}
CAP_UP, CAP_DN = 3.9, 0.0                     # cap height and descender, 5.5px
# the oligo pool: three strands, its label above them, its arrow below
WAVE_TOP, WAVE_DY, WAVE_N = -24.5, 4.2, 3
LABEL_GAP = 5.4                               # stack to label baseline
ARROW_GAP = 2.6                               # stack to where the arrow starts


def read_donor():
    """His style block, his backbone, and everything else as one string."""
    t = DONOR.read_text()
    style = re.search(r"<style>(.*?)</style>", t, re.S).group(1)
    body = t[t.index("</defs>") + len("</defs>"):t.rindex("</svg>")]
    body = re.sub(r"<title>.*?</title>", "", body, flags=re.S)
    back = re.search(r'<path id="Donor_plasmid_backbone"[^>]*/>', body).group(0)
    body = body.replace(back, "")
    # his heading is a title; it belongs to the caption
    body = re.sub(r'<text id="Donor_heading".*?</text>', "", body, flags=re.S)
    return style, back, body


LOOP_WIDEN = 4.0        # each side wall, so the att captions clear it inside
LOOP_LIFT = 14.0        # the bottom wall, to stop the loop eating the panel


def reshape(back, pad=LOOP_WIDEN, lift=LOOP_LIFT):
    """Widen the plasmid and raise its floor, without distorting its corners.

    Two problems, one fix. The att captions move inside the loop when the
    cassette goes to the top edge, and they are wider than the gap between
    their site and his corner, so they graze the side walls. And his rectangle
    is 32 deep because it used to sit above a cassette; below one it is mostly
    empty and takes a third of the panel for nothing.

    Scaling would stretch the corner arcs into ellipses. Moving every point on
    a given wall by the same amount leaves each corner exactly as he drew it:
    left wall left, right wall right, bottom wall up. The top wall is the
    strand and never moves.
    """
    n = [float(v) for v in re.findall(r"-?[\d.]+", back)]
    midx = (min(n[0::2]) + max(n[0::2])) / 2
    midy = (min(n[1::2]) + max(n[1::2])) / 2
    i = [0]

    def shift(m):
        i[0] += 1
        v = float(m.group(0))
        if i[0] % 2:                                   # an x
            return f"{v - pad:g}" if v < midx else f"{v + pad:g}"
        return f"{v - lift:g}" if v > midy else m.group(0)

    head, _, tail = back.partition(' d="')
    d, _, rest = tail.partition('"')
    moved = re.sub(r"-?[\d.]+", shift, d)
    return head + ' d="' + moved + '"' + rest


def captions(body, dy):
    """Where the module captions actually sit, as boxes."""
    out = []
    for m in re.finditer(r'<text[^>]*matrix\(1 0 0 1 ([\d.-]+) ([\d.-]+)\)'
                         r'[^>]*>\s*<tspan[^>]*>([^<]*)', body):
        cx, cy, txt = float(m.group(1)), float(m.group(2)) + dy, m.group(3)
        if ("module" not in txt
                and txt not in ("Recorder", "Sensor", "Target array")):
            continue
        w = sum(_ADV.get(c, .556) for c in txt) * 5.5
        out.append((cx - w / 2, cx + w / 2, cy - CAP_UP, cy + CAP_DN))
    return out


def check_leader(a, b, boxes, n=60):
    """A leader must miss every caption over its WHOLE height.

    Sampling only the baseline is how the left leader came to run through
    "Sensor module" twice: it cleared at y = 1 and cut the text's top edge.
    """
    for i in range(n + 1):
        t = i / n
        x, y = a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t
        for x0, x1, y0, y1 in boxes:
            if x0 <= x <= x1 and y0 <= y <= y1:
                sys.exit(f"leader {a}->{b} crosses a caption at "
                         f"({x:.1f}, {y:.1f}) - move BOX")


def box_of(body, pid):
    """The x span and the top/bottom y of one of his elements."""
    d = re.search(rf'id="{re.escape(pid)}"[^>]*\sd="([^"]+)"', body).group(1)
    n = [float(v) for v in re.findall(r"-?[\d.]+", d)]
    return min(n[0::2]), max(n[0::2]), min(n[1::2]), max(n[1::2])


def glyph(name, x0, y0, x1, y1):
    """One Fig1f scaffold glyph, fitted to the box by its INK.

    Fitting these by their viewBox is what kept the targets looking small: a
    matplotlib figure is mostly margin, and the Y array's drawing is 62 x 21
    inside an 88 x 59 box - three quarters of what was being scaled was empty
    padding. Measure the paths instead.
    """
    t = (GLYPHS / f"{name}.svg").read_text()
    # matplotlib's last </g> closes figure_1, so slicing to it leaves axes_1 whole
    body = t[t.index('<g id="axes_1">'):t.rindex("</g>")]
    body = re.sub(r'\sclip-path="url\(#[^)]*\)"', "", body)
    body = re.sub(r'\sid="[^"]*"', "", body)
    xs, ys = [], []
    for d in re.findall(r'\sd="([^"]+)"', body):
        n = [float(v) for v in re.findall(r"-?[\d.]+", d)]
        xs += n[0::2]
        ys += n[1::2]
    pad = 0.7                                  # half the heaviest stroke
    ix0, iy0 = min(xs) - pad, min(ys) - pad
    ix1, iy1 = max(xs) + pad, max(ys) + pad
    k = min((x1 - x0) / (ix1 - ix0), (y1 - y0) / (iy1 - iy0))
    tx = (x0 + x1) / 2 - k * (ix0 + ix1) / 2
    ty = (y0 + y1) / 2 - k * (iy0 + iy1) / 2
    return (f'<g transform="translate({tx:.3f} {ty:.3f}) scale({k:.5f})">'
            f"{body}</g>")


def main():
    style, back, body = read_donor()
    _, _, rect_top, rect_bot = box_of(back, "Donor_plasmid_backbone")
    dy = rect_top - rect_bot                   # -32: bottom edge up to the top

    bx0, by0, bx1, by1 = BOX
    b0, b1, btop, _ = box_of(body, "Binder")
    s0, s1, _, _ = box_of(body, "RNA_scaffold")
    btop += dy                                 # the binder's top edge, moved

    o = [reshape(back), f'<g transform="translate(0 {dy:g})">{body}</g>']

    # --- the oligo pool, dropping straight into the middle of the binder.
    # The sensor-module bracket runs right over the binder, so the arrow has
    # to cross it. It crosses the way Matteo's own exchange curves cross: a
    # white casing under the ink, so the arrow reads as passing in front.
    px = (b0 + b1) / 2
    py = WAVE_TOP + (WAVE_N - 1) * WAVE_DY + ARROW_GAP
    tipy = btop - 0.5
    for i in range(WAVE_N):
        y = WAVE_TOP + i * WAVE_DY
        o.append(
            f'<path class="wave" stroke="'
            f'{palette.THEME["600" if i % 2 else "400"]}" '
            f'd="M{px - 13:g} {y:g} C{px - 8:g} {y - 2.2:g} {px - 5:g} '
            f'{y + 2.2:g} {px:g} {y:g} C{px + 5:g} {y - 2.2:g} {px + 8:g} '
            f'{y + 2.2:g} {px + 13:g} {y:g}"/>')
    # above the stack, not beside it: set to its right the label sat between
    # the first and second wave and crowded them. Shares the call-out title's
    # baseline, so the two headings read as one row.
    o.append(f'<text class="lbl" x="{px:g}" y="{WAVE_TOP - LABEL_GAP:g}">'
             # 2026-11-18 (Leo, from his supervisor): name the pool size
             # here. Leo's wording, chosen over "12,000-design": shorter, and
             # "binder" says what the 12k ARE, which the long form did not.
             f"12,000 binder oligo pool</text>")
    o.append(f'<path class="casing" d="M{px:g} {py:g} L{px:g} {tipy - 2:g}"/>'
             f'<path class="arrow" d="M{px:g} {py:g} L{px:g} {tipy - 2.6:g}"/>'
             f'<path fill="{INK}" d="M{px:g} {tipy:g} '
             f'L{px - 1.7:g} {tipy - 2.8:g} L{px + 1.7:g} {tipy - 2.8:g} Z"/>')

    # --- the target call-out. Its leaders fan OUTWARD from the scaffold, so
    # they too pass outside the target-module bracket instead of over it.
    # The leaders land on the TARGET MODULE BRACKET's ends, not on the block
    # itself. Running them down to the block means crossing the bracket row,
    # and with the call-out this wide the left one cuts straight through the
    # sensor-module bracket. Stopping at the bracket also says the right
    # thing: what the call-out expands is the module that bracket names.
    tb0, tb1, tbt, _ = box_of(body, "Target_module_bracket")
    tby = tbt + dy
    # Start on the straight part of the bottom edge, not at the rect's corner
    # point: the corners are rounded, so a leader from the corner itself
    # springs out of the curve and reads as misaligned.
    lead = [((bx0 + BOX_R, by1), (tb0, tby)), ((bx1 - BOX_R, by1), (tb1, tby))]
    for a, b in lead:
        check_leader(a, b, captions(body, dy))
    o.append("".join(f'<path class="lead" d="M{a[0]:g} {a[1]:g} '
                     f'L{b[0]:g} {b[1]:g}"/>' for a, b in lead))
    o.append(f'<rect class="callout" x="{bx0:g}" y="{by0:g}" '
             f'width="{bx1 - bx0:g}" height="{by1 - by0:g}" rx="2"/>')
    # the title sits inside the box, on its own band, rather than floating
    # above it - that gap was the panel's tallest piece of nothing
    tb = by0 + TITLE_H
    # 2026-11-18 (Leo): "don't call it target scaffold but RNA presenting
    # module" - the box shows the four things the RNA is presented ON, and
    # "scaffold" collided with the scaffold vocabulary used elsewhere. Left
    # aligned inside the box, so the longer string just runs further right
    # into the band rather than needing the box resized.
    # 2026-11-23 (Leo): the thing itself is now the TARGET ARRAY, which the
    # bracket under the cassette already names, so the band says what the four
    # cells are instead of repeating it - the caption's own words.
    o.append(f'<text class="ttl" x="{bx0 + 3.5:g}" y="{tb - 2.4:g}">'
             f'Four array architectures</text>')
    cw = (bx1 - bx0) / len(CELLS)
    o.append(f'<path class="grid" d="M{bx0 + 1.5:g} {tb:g} '
             f'L{bx1 - 1.5:g} {tb:g}"/>')
    for i, name in enumerate(CELLS):
        cx0 = bx0 + i * cw
        if i:
            o.append(f'<path class="grid" d="M{cx0:g} {tb + 1.5:g} '
                     f'L{cx0:g} {by1 - 1.5:g}"/>')
        o.append(glyph(name, cx0 + 2.5, tb + 1.6, cx0 + cw - 2.5, by1 - 1.8))

    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="-10 -38 188 72">'
        f"<style>{style}"
        f".lbl{{font-family:Arial;font-size:5.5px;fill:{INK};text-anchor:middle}}"
        f".ttl{{font-family:Arial;font-size:5.5px;fill:{INK}}}"
        f".wave{{fill:none;stroke-width:1.1;stroke-linecap:round}}"
        f".arrow{{fill:none;stroke:{INK};stroke-width:0.5}}"
        f".casing{{fill:none;stroke:{palette.SURFACE};stroke-width:1.9}}"
        f".lead{{fill:none;stroke:{palette.NEUTRAL['400']};stroke-width:0.3}}"
        f".callout{{fill:none;stroke:{OUTLINE};stroke-width:0.4}}"
        f".grid{{fill:none;stroke:{palette.NEUTRAL['300']};stroke-width:0.25}}"
        "</style>" + "".join(o) + "</svg>")

    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "library_cassette.svg"
    p.write_text(svg, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-w", "2000", "-b", "white", str(p),
                    "-o", str(p.with_suffix(".png"))], check=False)
    subprocess.run(["inkscape", str(p), "--export-type=pdf",
                    f"--export-filename={p.with_suffix('.pdf')}"], check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    stray = sorted({h.lower() for h in re.findall(r"#[0-9a-fA-F]{6}", svg)} -
                   {c.lower() for c in (*palette.THEME.values(),
                                        *palette.NEUTRAL.values(),
                                        *palette.SYNBIO.values(),
                                        *palette.DIVERGING, palette.INK,
                                        palette.GRID, palette.SURFACE)})
    print(f"cassette moved up {-dy:g} onto the plasmid's top edge; "
          f"backbone untouched" + (f"  LEFTOVER {stray}" if stray else
                                   "  palette ok"))


if __name__ == "__main__":
    main()
