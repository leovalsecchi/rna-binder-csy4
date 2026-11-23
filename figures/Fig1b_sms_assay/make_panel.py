"""Fig1 panel B: the SMS reporter assay, from bound target to mScarlet.

Three zones: the sensor on the left, the three transfected constructs on the
right, and the two reporter states below.

Nothing here draws its own construct bars or its own molecules. The bars come
from `common.constructs` (so they match panel d), the hairpins come from
`cartoon_library.hairpin`, the SMS module is Matteo's drawing through
`common.cartoon`, and the scaffold is
`cartoon_library.f30`, imported rather than copied so that panel picks up any
change to the F30 model without being touched.

No title, no standfirst, no panel letter (STYLE_NOTES.md items 1 and 8).

Writes output/Fig1b_sms_assay/Fig1b_sms_assay.svg.
"""
import re
import sys
from pathlib import Path
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parents[1] / "cartoon_library"))

import math
import yaml
from common import paths, palette, cartoon
from common import constructs as K
import f30
import hairpin as HP

out = paths.fig_output("Fig1b_sms_assay")

_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
DPI = int(_cfg["dpi"])
WIDTH_MM = float(_cfg["width_double_mm"])

S = palette.SYNBIO
LINE = S["outline"]
INK = palette.INK
MUT = palette.NEUTRAL["600"]
RULE = palette.NEUTRAL["300"]
BG = palette.SURFACE

W, H = 1600, 600
SENSOR_CX, SENSOR_CY = 330, 205      # the scaffold sits here
SENSOR_K = 1.3 / f30.K     # f30 is drawn at the AT hairpin's size now,
                           # so divide out the factor that bought it
CONSTRUCT_X, CONSTRUCT_Y = 780, 78   # the ruled block
STATE_Y = 530                        # the two reporter mRNAs


def T(x, y, s, sz=13, col=INK, w=400, anch="middle", ls=0):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{sz}" '
            f'fill="{col}" font-weight="{w}" text-anchor="{anch}" '
            f'letter-spacing="{ls}">{s}</text>')


def stroke(d, col=LINE, w=2.0, dash=None):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{w}" '
            f'stroke-linecap="round" stroke-linejoin="round"{extra}/>')


def arrow(x1, y1, x2, y2, col=palette.NEUTRAL["400"], w=1.9, head=9):
    a = math.atan2(y2 - y1, x2 - x1)
    bx, by = x2 - head * math.cos(a), y2 - head * math.sin(a)
    px, py = -math.sin(a) * head * 0.6, math.cos(a) * head * 0.6
    return (stroke(f"M{x1:.1f},{y1:.1f} L{bx:.1f},{by:.1f}", col, w)
            + f'<path d="M{x2:.1f},{y2:.1f} L{bx+px:.1f},{by+py:.1f} '
              f'L{bx-px:.1f},{by-py:.1f} Z" fill="{col}"/>')


_SVG = re.compile(r'<svg[^>]*width="([\d.]+)"[^>]*height="([\d.]+)"[^>]*'
                  r'viewBox="([^"]+)"[^>]*>(.*)</svg>\s*$', re.S)


def embed(svg_text, cx, cy, scale=1.0):
    """Place another module's finished SVG, centred on (cx, cy).

    Returns the group and a function mapping that module's own coordinates into
    this panel's, so things can be anchored to it rather than to guessed pixels.
    """
    m = _SVG.match(svg_text)
    if not m:
        raise ValueError("cannot read that SVG's header")
    w, h = float(m.group(1)), float(m.group(2))
    vx, vy = (float(v) for v in m.group(3).split()[:2])
    dx = cx - (vx + w / 2) * scale
    dy = cy - (vy + h / 2) * scale
    to_panel = lambda pt: (dx + pt[0] * scale, dy + pt[1] * scale)
    box = (vx, vy, w, h)                 # the source's own viewBox
    return (f'<g transform="translate({dx:.1f},{dy:.1f}) scale({scale})">'
            + m.group(4) + '</g>'), to_panel, box


# ----------------------------------------------------------------- sensor ---
def _arm_tips():
    """Where the Csy4 hairpins end, in the scaffold module's own coordinates.

    Read off f30's constants rather than hardcoded, so the receptors stay on the
    arms if that model changes. It is better to fail here than to drift quietly.
    """
    try:
        reach = (f30.JUNC_GAP + f30.JUNC_LEN + f30.BULGE_LEN + f30.ARM_LEN)
        return [f30._along(f30.J, d, reach) for d in (f30._DL, f30._DR)]
    except AttributeError as e:          # f30 was renamed under us
        raise RuntimeError(
            f"cannot find the arm geometry in cartoon_library.f30 ({e}); "
            f"panel b anchors its receptors to it") from e


def receptor(tip, out_dir):
    """Binder clamping the hairpin, with the SMS module hanging below it.

    Returns the svg and where the SMS module landed, so the caption can sit
    under it rather than at a guessed height.
    """
    # 28, not the 40 the wide helix bundle needed: the crown is a compact
    # cup, so at 40 it floated off the hairpin it is meant to be gripping
    bx, by = tip[0] + out_dir[0] * 28, tip[1] + out_dir[1] * 28
    sx, sy = bx + out_dir[0] * 16, by + 94
    return "".join([
        stroke(f"M{bx:.0f},{by+20:.0f} L{sx:.0f},{sy-30:.0f}", LINE, 2.2),
        # the crown (his VHH cup), not the helix bundle: a binder should not
        # read as "some protein". 3.6 keeps its width close to the bundle's,
        # which the linker and the caption below are spaced for.
        cartoon.fragment("crown_monomer.svg", "Nanobody v7.9", scale=3.6,
                         at=(bx, by), anchor="center"),
        cartoon.fragment("Selected_SMS_v6_5_1x.svg", scale=3.4,
                         at=(sx, sy), anchor="center"),
    ]), (sx, sy)


def sensor():
    svg, to_panel, box = embed(f30.build(), SENSOR_CX, SENSOR_CY, SENSOR_K)
    o = [svg]
    vx, vy, vw, vh = box

    sms_at = []
    for tip, d in zip([to_panel(t) for t in _arm_tips()],
                      ((-0.57, -0.82), (0.57, -0.82))):
        part, c = receptor(tip, d)
        o.append(part)
        sms_at.append(c)
    for (sx, sy) in sms_at:
        o.append(T(sx, sy + 58, "SMS receptor", 12, MUT))

    # mNeonGreen is translated from the same transcript (construct 1). No
    # leader to the strand: the construct map already states where it sits, and
    # a line to a point this drawing does not expose would only drift.
    gx, gy = to_panel((vx - 18, vy + vh - 30))
    o.append(cartoon.fragment("Selected_GFP_v1_1x.svg", scale=3.0,
                              at=(gx, gy), anchor="center"))
    o.append(T(gx, gy + 46, "mNeonGreen", 11.5, MUT))

    bottom = max(to_panel((vx, vy + vh))[1], max(c[1] for c in sms_at) + 58)
    o.append(T(SENSOR_CX, bottom + 34, "Target array", 13.5, INK, 600))
    return "".join(o), sms_at[1]


# ------------------------------------------------------------- constructs ---
ROWS = [("reporter_rna", "target RNA"),
        ("sensor", "SMS receptor"),
        ("readout", "reporter")]


def constructs(x0=CONSTRUCT_X, y0=CONSTRUCT_Y, width=760):
    o = [stroke(f"M{x0},{y0} L{x0+width},{y0}", RULE, 1.0),
         T(x0, y0 - 12, "CONSTRUCTS", 10.5, MUT, 600, "start", 1.4)]
    for i, (name, note) in enumerate(ROWS):
        y = y0 + 54 + i * 66
        o.append(f'<circle cx="{x0+13}" cy="{y}" r="11" fill="none" '
                 f'stroke="{RULE}" stroke-width="1.4"/>')
        o.append(T(x0 + 13, y + 4.4, str(i + 1), 11.5, MUT, 600))
        o.append(K.bar(x0 + 46, y, K.BARS[name]))
        o.append(T(x0 + 46 + K.width(K.BARS[name]) + 30, y + 4, note, 11,
                   MUT, 400, "start"))
    o.append(stroke(f"M{x0},{y0+228} L{x0+width},{y0+228}", RULE, 1.0))
    return "".join(o)


# ----------------------------------------------------------------- states ---
HP_H = 74
HP_X, STOP_X, SCAR_X = 190, 272, 330


def hairpin(spliced, cx, baseline):
    """The AT hairpin, from its own builder.

    Not `cartoon.fragment("v19.18.svg", ...)`, which this used to call: that is
    a different, compact hairpin in Matteo's own colour assignment, where the
    exon is blue and the intron grey - the opposite way round from the drawing
    this panel is meant to carry.
    """
    svg = HP.build(spliced=spliced)
    w, h = (float(v) for v in
            re.search(r'width="([\d.]+)" height="([\d.]+)"', svg).groups())
    k = HP_H / h
    body = re.sub(r"^<svg[^>]*>|</svg>$", "", svg)
    return (f'<g transform="translate({cx - w * k / 2:.2f},'
            f'{baseline - HP_H:.2f}) scale({k:.4f})">{body}</g>')


def state(x0, spliced):
    """One reporter mRNA. The hairpin IS the AT module, so it is not boxed too."""
    y = STATE_Y
    o = [T(x0, y - 104, "Spliced" if spliced else "Unspliced", 13.5, INK,
           600, "start"),
         stroke(f"M{x0},{y} L{x0+600},{y}", LINE, 2.0),
         K.bar(x0 + 40, y, ["aug"]),
         hairpin(spliced, x0 + HP_X, y),
         K.bar(x0 + SCAR_X, y, ["mscarlet" if spliced else "mscarlet_off"])]
    if not spliced:                      # the stop sits ahead of the reporter
        o.append(K.bar(x0 + STOP_X, y, ["stop"]))
    note = ("cleaved: mScarlet now in frame" if spliced
            else "stop codon in frame, no reporter")
    o.append(T(x0, y + 58, note, 11, MUT, 400, "start"))
    return "".join(o)


def build():
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" '
         f'height="{WIDTH_MM*H/W:.2f}mm" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>',
         ]
    sensor_svg, right_sms = sensor()
    o += [sensor_svg, constructs(), state(80, False), state(890, True)]
    o.append(arrow(726, STATE_Y, 866, STATE_Y))
    o.append(T(796, STATE_Y - 14, "SMS RNase", 11, MUT))
    # the activated receptor is what cleaves; the curve starts at its RNase
    rx, ry = right_sms
    o.append(stroke(f"M{rx+26:.0f},{ry+22:.0f} C{rx+140:.0f},{ry+90:.0f} "
                    f"{680:.0f},{STATE_Y-80:.0f} {726:.0f},{STATE_Y-16:.0f}",
                    RULE, 1.6, dash="5 6"))
    o.append('</svg>')
    return "\n".join(o)


if __name__ == "__main__":
    svg_path = out / "Fig1b_sms_assay.svg"
    svg_path.write_text(build(), encoding="utf-8")
    pdf, png = out / "Fig1b_sms_assay.pdf", out / "Fig1b_sms_assay.png"
    try:
        import cairosvg
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf))
        cairosvg.svg2png(url=str(svg_path), write_to=str(png),
                         output_width=int(W * DPI / 96))
    except ImportError:
        for stale in (pdf, png):
            stale.unlink(missing_ok=True)
        print("SVG only (no cairosvg).")
    print("wrote", svg_path)
