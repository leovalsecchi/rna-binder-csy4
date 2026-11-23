"""Fig 1d, right-hand panel: the mechanism, and what it does to the recorder.

Three pieces, none of them redrawn here:

  * `working/Fig1d_integrated_construct/Fig1d_integrated_construct.svg` - the
    construct, used TWICE: once as drawn, once with the recorder spliced;
  * `output/cartoon_library/svg/Hairpin Splicing Mechanism V2 compact.svg` -
    Leo's hand-built mechanism with the apo state and the released intron
    dropped, in a pale container on the left. The full three-state version is
    kept beside it as `Hairpin Splicing Mechanism V2.svg` and is not used
    here;
  * `hairpin.build(spliced=True)` and `(intron=True)` - the cut hairpin and
    the piece it loses.

A swung arrow runs from the first row to the second and passes BEHIND the
container, so the mechanism sits inside the step it causes rather than beside
it. Draw order carries that: arrow, then container, then mechanism, then rows.

**The two rows stay in register.** Only the hairpin differs; the ATG, binder,
SMS and target keep the same x in both, so the reader compares them by looking
straight down. The spliced hairpin is much the narrower, and the slack is left
as a longer 5' lead rather than being closed up.

    python code/Fig1d_integrated_construct/make_panel.py [outdir]
"""
import math
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
REPO = _HERE.parents[2]
sys.path.insert(0, str(REPO / "code"))
sys.path.insert(0, str(REPO / "code" / "cartoon_library"))

import hairpin as AT                                         # noqa: E402
from common import palette                                   # noqa: E402

WORK = REPO / "working" / "Fig1d_integrated_construct"
LIB = REPO / "output" / "cartoon_library" / "svg"
CONSTRUCT = WORK / "Fig1d_integrated_construct.svg"
MECHANISM = LIB / "Hairpin Splicing Mechanism V2 compact.svg"

# --- the construct's own frame, read off the file it comes from ------------
STEM_X = 18.476        # the hairpin's stem, and so the recorder's centre
STRAND_Y = 30.28
HP_SCALE = 0.42        # construct units per hairpin-module unit
EXON_W, RUNG_W = 1.714, 1.19     # local widths that land on Matteo's 0.72/0.5

# --- layout ----------------------------------------------------------------
BOX_W = 100.0          # the container holding the mechanism
BOX_PAD = 7.0
GAP = 34.0             # container to constructs: wide enough that the swung
                       # arrow has somewhere to be seen
ROW_DY = 80.0          # first row to second. The compact mechanism drops the
                       # apo state and the released intron, so the container
                       # is shorter and the two rows can close up.
MARGIN = 4.0
ARC_PULL = 92.0        # how far left the swung arrow reaches, i.e. how much
                       # of it the container hides

# The expressed product, under the TOP row only. PROD_X is the lane between
# the Binder and SMS captions, which is the one place the arrow can drop out
# of the row without crossing a caption.
MONO_VB, MONO_W, MONO_H = (32.3, 41.2), 26.1, 50.5
PROD_X = 81.3
PROD_TOP, PROD_H = 48.0, 31.0

NEUT = palette.NEUTRAL
MUT = NEUT["500"]


def _head(x, y, ang, size, col):
    """A filled wedge at (x, y) pointing along ang."""
    pts = [(0.0, 0.0), (-size, -size * 0.42), (-size, size * 0.42)]
    c, s = math.cos(ang), math.sin(ang)
    d = " ".join(f"{'M' if i == 0 else 'L'}{x + a * c - b * s:.2f},"
                 f"{y + a * s + b * c:.2f}" for i, (a, b) in enumerate(pts))
    return f'<path d="{d} Z" fill="{col}"/>'


def _inner(svg, tag="svg"):
    return svg[svg.index(">", svg.index(f"<{tag}")) + 1:svg.rindex(f"</{tag}>")]


def _construct():
    """The construct's <style> and the body of its one group."""
    t = CONSTRUCT.read_text()
    style = re.search(r"<style>.*?</style>", t, re.S).group(0)
    body = _inner(t[t.index('<g id="Gene_circuit"'):] + "", "g")
    return style, body


def _spliced_hairpin():
    """The cut hairpin, at the construct's scale and line weights.

    hairpin.py draws it for a standalone figure, where one unit is 7 across;
    here one unit is 0.42, so the strokes are reset to the widths that land on
    the construct's own 0.72 and 0.5 rather than being scaled down with it.
    """
    s = AT.build(spliced=True)
    body = _inner(s[s.index("<g"):], "g")
    body = re.sub(r'(stroke="#0678a6"[^/]*?stroke-width=")[\d.]+',
                  lambda m: m.group(1) + f"{EXON_W}", body)
    body = re.sub(r'(stroke="#b3b8bc"[^/]*?stroke-width=")[\d.]+',
                  lambda m: m.group(1) + f"{RUNG_W}", body)
    return (f'<g data-name="AT hairpin | spliced" '
            f'transform="translate({STEM_X} {STRAND_Y}) scale({HP_SCALE}) '
            f'translate(0 {-AT.SP_BOT})">{body}</g>')


def _intron(cx, cy, scale):
    s = AT.build(intron=True)
    body = _inner(s[s.index("<g"):], "g")
    body = re.sub(r'(stroke-width=")[\d.]+', lambda m: m.group(1) + "1.5", body)
    return (f'<g data-name="Excised intron" transform="translate({cx} {cy}) '
            f'scale({scale}) translate({-AT.DRIFT_LEN / 2:.2f} 0)">{body}</g>')


def _splice(body):
    """Swap the uncut hairpin for the cut one, and let the bracket follow."""
    out = re.sub(r'<g id="AT_hairpin".*?</g>', _spliced_hairpin(), body, flags=re.S)
    if out == body:
        raise SystemExit("the hairpin group moved - update the pattern")
    half = (AT.SP_R + EXON_W / 2) * HP_SCALE
    a, b = STEM_X - half - 0.45, STEM_X + half + 0.45
    out = re.sub(r'(<path d=")M1\.96 13\.58[^"]*(")',
                 lambda m: (f'{m.group(1)}M{a:.2f} 13.58 L{a:.2f} 11.00 '
                            f'L{b:.2f} 11.00 L{b:.2f} 13.58{m.group(2)}'), out)
    # the piece the recorder has just lost, drifting in the slack it left
    return out + _intron(32.0, 24.6, 0.30)


def _expression():
    """A production arrow from the sensor CDS down to what it makes.

    SBOL Visual's way of saying a coding sequence is expressed: an arrow from
    the CDS to its product. The product here is the binder-SMS chain, drawn
    from the same library asset the mechanism uses, so the reader meets the
    same object twice rather than two drawings of it. It carries no caption -
    its three colours are the Binder and SMS blocks directly above it.

    The arrow drops at PROD_X, which is the gap between the Binder and SMS
    captions, so it reaches the band below without crossing either.
    """
    svg = _inner((LIB / "sms_crown_monomer.svg").read_text())
    style = re.search(r"<style>(.*?)</style>", svg, re.S).group(1)
    body = re.sub(r"<defs>.*?</defs>", "", svg, flags=re.S)
    style = f"<style>{style.replace('.cls-', '.prot-cls-')}</style>"
    body = body.replace('class="cls-', 'class="prot-cls-')
    k = PROD_H / MONO_H
    g = (f'<g data-name="Binder-SMS, expressed" '
         f'transform="translate({PROD_X - MONO_W * k / 2:.2f},{PROD_TOP:.2f}) '
         f'scale({k:.4f}) translate({-MONO_VB[0]},{-MONO_VB[1]})">{body}</g>')
    arrow = (f'<path data-name="Production" d="M{PROD_X} 34.4 '
             f'L{PROD_X} {PROD_TOP - 3.4:.2f}" fill="none" stroke="{MUT}" '
             f'stroke-width="0.5" stroke-linecap="round"/>'
             + _head(PROD_X, PROD_TOP - 1.2, math.pi / 2, 2.0, MUT))
    return style, arrow + g


def _mechanism(x, y, h):
    t = MECHANISM.read_text()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]
    k = h / vb[3]
    style = re.search(r"<style>(.*?)</style>", t, re.S).group(1)
    body = re.sub(r"<defs>.*?</defs>", "", _inner(t), flags=re.S)
    style = f"<style>{style.replace('.cls-', '.mech-cls-')}</style>"
    body = body.replace('class="cls-', 'class="mech-cls-')
    w = vb[2] * k
    g = (f'<g transform="translate({x:.2f},{y:.2f}) scale({k:.5f}) '
         f'translate({-vb[0]},{-vb[1]})">{body}</g>')
    return style, g, w


def build():
    c_style, c_body = _construct()
    rows_x = BOX_W + GAP
    box_h = ROW_DY + 44.0 - 2.30
    m_style, mech, m_w = _mechanism(0, 0, box_h - 2 * BOX_PAD)
    mech = mech.replace("translate(0.00,0.00)",
                        f"translate({(BOX_W - m_w) / 2:.2f},{BOX_PAD:.2f})", 1)

    # The swung arrow, emitted first so the container covers its middle. It
    # starts and ends clear of the 5' captions, and its head follows the
    # curve's own tangent rather than a guessed angle.
    x0 = rows_x - 9.0
    y0, y1 = STRAND_Y - 1.0, STRAND_Y + ROW_DY + 1.0
    c1, c2 = (ARC_PULL, y0 + 26.0), (ARC_PULL, y1 - 26.0)
    ang = math.atan2(y1 - c2[1], x0 - c2[0])
    arc = (f'<path d="M{x0:.2f},{y0:.2f} C{c1[0]:.2f},{c1[1]:.2f} '
           f'{c2[0]:.2f},{c2[1]:.2f} {x0:.2f},{y1:.2f}" fill="none" '
           f'stroke="{MUT}" stroke-width="1.6" stroke-linecap="round"/>'
           + _head(x0, y1, ang, 5.2, MUT))
    box = (f'<rect x="0" y="0" width="{BOX_W}" height="{box_h:.2f}" rx="4" '
           f'fill="{NEUT["100"]}" stroke="{NEUT["300"]}" stroke-width="0.4"/>')

    p_style, product = _expression()
    rows = (f'<g transform="translate({rows_x},0)">{c_body}{product}</g>'
            f'<g transform="translate({rows_x},{ROW_DY})">{_splice(c_body)}</g>')

    w = rows_x + 174.0 + MARGIN
    y_min = -MARGIN
    h = box_h + 2 * MARGIN
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w * 6:.0f}" '
            f'height="{h * 6:.0f}" viewBox="{-MARGIN} {y_min} '
            f'{w + MARGIN:.2f} {h:.2f}">'
            + c_style + m_style + p_style + arc + box + mech + rows
            + "</svg>")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "Fig1d_panel_mechanism.svg"
    p.write_text(build(), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else WORK)
