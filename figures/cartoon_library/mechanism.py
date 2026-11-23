"""The sensing mechanism: target binding -> dimerization -> cleavage.

Nothing is drawn from scratch; the library already holds every piece:

  * `f30.build()`             the target - F30 junction, two Csy4 hairpins
  * `sms_crown_dimer.svg`     two chains, dimerized through their kinase lobes
  * `sms_crown_monomer.svg`   one chain, for the unbound state
  * `hairpin.build()`         the AT hairpin, uncut, and the released intron

Arrangement follows Matteo's concept sketch - target on top, chains hanging
from its arms, substrate below - with two deliberate departures:

  * the AT hairpin sits **on the recorder mRNA**, a 5'-3' line, so the cut
    visibly yields the spliced recorder the readout rows show; his floats free;
  * the cut is marked by **arrowheads at the two sites**, the RNA-processing
    convention, rather than by scissors, which do not survive 7 pt.

The F30 is flipped, so its arms splay downward and the bottom stem - the rest
of the transcript - leaves at the top.

**Why the dimer is one asset and not two monomers.** Two straight chains tilted
far enough inward for their kinase lobes to touch send their RNase lobes past
each other, because the RNase sits below the kinase and a straight chain cannot
bend. A real IRE1 dimer bends there. `sms_crown_dimer` already holds that pose,
so the panel is built in ITS units and the F30 is scaled until the span of its
two arms equals the span of the dimer's two crowns. Change the F30's arm angle
and the fit still holds.

    python code/cartoon_library/mechanism.py [outdir]
"""
import math
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parent))

import f30                                                   # noqa: E402
import hairpin as AT                                         # noqa: E402
from common import palette                                   # noqa: E402

REPO = _HERE.parents[2]
LIB = REPO / "output/cartoon_library/svg"

INK = palette.INK
MUT = palette.NEUTRAL["500"]
RNA = palette.NEUTRAL["700"]          # the recorder mRNA backbone
FONT = "Arial, Helvetica, sans-serif"

# --- measured once off the crown assets, in their own viewport units --------
# Both files were cropped by crown_receptor.frame(), so a viewBox that does not
# start at the origin still has to be shifted out before anything is placed.
DIMER_VB = (17.3, 27.7)
DIMER_W, DIMER_H = 43.2, 50.5
CROWN_L, CROWN_R = 11.29, 31.93       # the two crown mouths
CROWN_Y = 4.94                        # and their lip
MONO_VB = (32.3, 41.2)
MONO_W, MONO_H = 26.1, 50.5
MONO_CROWN = (11.30, 5.00)
MONO_RNASE = (16.09, 41.75)   # the lobe the chain hangs its weight on

GRIP_SINK = 1.1        # how far the hairpin loop settles into the crown
REC_SCALE = 2.15       # the recorder hairpin, against the target's arm span:
                       # it is the same kind of object as the rest, so it is
                       # drawn at the same size rather than as an inset
STEP = 13.0            # between one state and the next
ARROW = 9.0
CLEAR = 4.2            # empty space above a step arrow
DROP = 11.0            # and below it. Always the larger of the two: what a
                       # step arrow points at is an RNA stem that rises into
                       # the gap, so an even split still reads as crowding.
LEAN = 12.0            # how far an approaching chain tilts its crown inward
APPROACH = (24.0, 21.0)   # a waiting chain's crown, from the hairpin tip it
                          # will dock on: out by x, UP by y, so its SMS hangs
                          # low and outboard and the wedge between opens up
ARC_FROM = (7.5, 10.0)    # the motion arc leaves the crown's inner shoulder
ARC_TO = (-2.5, 2.2)      # and stops short of the tip, so its head reads
ARC_BOW = 4.5             # bowed outward and down, under the chord


def _inner(svg):
    return svg[svg.index(">", svg.index("<svg")) + 1:svg.rindex("</svg>")]


def _g(body, *xf):
    return f'<g transform="{" ".join(x for x in xf if x)}">{body}</g>'


def _asset(name, tag):
    """An asset's <style>, namespaced, and its body.

    Both crown files number their classes from cls-1, so dropping the two in
    one document unchanged would let the second file's rules repaint the
    first. Each gets its own prefix instead, exactly as the integration
    assembler does.
    """
    svg = _inner((LIB / f"{name}.svg").read_text())
    style = re.search(r"<style>(.*?)</style>", svg, re.S)
    body = re.sub(r"<defs>.*?</defs>", "", svg, flags=re.S)
    css = (style.group(1).replace(".cls-", f".{tag}-cls-") if style else "")
    return f"<style>{css}</style>", body.replace('class="cls-',
                                                 f'class="{tag}-cls-')


def _arrowhead(x, y, ang, size, col=INK):
    """A filled wedge at (x, y) pointing along ang (radians, 0 = +x)."""
    pts = [(0.0, 0.0), (-size, -size * 0.42), (-size, size * 0.42)]
    c, s = math.cos(ang), math.sin(ang)
    d = " ".join(f"{'M' if i == 0 else 'L'}{x + a * c - b * s:.2f},"
                 f"{y + a * s + b * c:.2f}" for i, (a, b) in enumerate(pts))
    return f'<path d="{d} Z" fill="{col}"/>'


def _curve_arrow(p0, p1, bulge, size=1.9, col=MUT):
    """An arc from p0 to p1, bowed `bulge` off the chord, head at p1.

    One quadratic, so the head's direction is just the tangent at the end -
    the control point - and the arc cannot kink.
    """
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = math.hypot(dx, dy) or 1.0
    c = ((p0[0] + p1[0]) / 2 - dy / n * bulge,
         (p0[1] + p1[1]) / 2 + dx / n * bulge)
    return (f'<path d="M{p0[0]:.2f},{p0[1]:.2f} Q{c[0]:.2f},{c[1]:.2f} '
            f'{p1[0]:.2f},{p1[1]:.2f}" fill="none" stroke="{col}" '
            f'stroke-width="0.6" stroke-linecap="round"/>'
            + _arrowhead(p1[0], p1[1],
                         math.atan2(p1[1] - c[1], p1[0] - c[0]), size, col))


def _step_arrow(x, y, length=ARROW, col=MUT):
    head = 2.2
    return (f'<path d="M{x:.2f},{y:.2f} L{x:.2f},{y + length - head:.2f}" '
            f'fill="none" stroke="{col}" stroke-width="0.6" '
            f'stroke-linecap="round"/>'
            + _arrowhead(x, y + length, math.pi / 2, head, col))


# --- the target -------------------------------------------------------------
def _target(span):
    """The F30, flipped and scaled so its two arm tips are `span` apart.

    Returns the placed body, the two tips and the bbox, all in panel units.
    """
    f30._pts.clear()
    body = _inner(f30.build())
    reach = f30.JUNC_GAP + f30.JUNC_LEN + f30.BULGE_LEN + f30.ARM_LEN
    tl = f30._along(f30.J, f30._DL, reach)
    tr = f30._along(f30.J, f30._DR, reach)
    k = span / (tr[0] - tl[0])
    xs, ys = [p[0] for p in f30._pts], [p[1] for p in f30._pts]
    flip = [(x * k, -y * k) for x, y in zip(xs, ys)]
    tlf, trf = (tl[0] * k, -tl[1] * k), (tr[0] * k, -tr[1] * k)
    box = (min(p[0] for p in flip), min(p[1] for p in flip),
           max(p[0] for p in flip), max(p[1] for p in flip))
    return (body, k, tlf, trf, box)


def _place_target(span, tip_to):
    """The F30 placed so its left arm tip lands on `tip_to`."""
    body, k, tlf, trf, box = _target(span)
    dx, dy = tip_to[0] - tlf[0], tip_to[1] - tlf[1]
    placed = _g(body, f"translate({dx:.3f},{dy:.3f})", f"scale({k:.5f})",
                "scale(1,-1)")
    box = (box[0] + dx, box[1] + dy, box[2] + dx, box[3] + dy)
    return placed, box


# --- the substrate ----------------------------------------------------------
def _recorder(cx, cy, k, cuts=True):
    """The AT hairpin on the recorder mRNA, with the two cut sites marked.

    Both the half-width and the cut sites come from hairpin.py: the exon meets
    the intron at (+-(SIDE_CX + SIDE_R), ARM_Y), so this follows that module
    rather than restating its numbers.
    """
    body = _inner(AT.build())
    body = body[body.index(">", body.index("<g")) + 1:body.rindex("</g>")]
    s = k * AT.SCALE
    far, foot = AT.SIDE_CX + AT.SIDE_R, AT.STEM_BOT
    half = far * s
    lead = half + 11.0
    o = [f'<path d="M{cx - lead:.2f},{cy:.2f} L{cx + lead:.2f},{cy:.2f}" '
         f'fill="none" stroke="{RNA}" stroke-width="0.45" '
         f'stroke-linecap="round"/>',
         _g(body, f"translate({cx:.2f},{cy:.2f})", f"scale({s:.5f})",
            f"translate(0,{-foot:.2f})")]
    if cuts:
        y = cy - (foot - AT.ARM_Y) * s
        for sgn in (-1, 1):
            o.append(_arrowhead(cx + sgn * (half + 1.6), y,
                                0.0 if sgn < 0 else math.pi, 2.0))
    return "".join(o), half


def _intron(cx, cy, k):
    body = _inner(AT.build(intron=True))
    body = body[body.index(">", body.index("<g")) + 1:body.rindex("</g>")]
    return _g(body, f"translate({cx:.2f},{cy:.2f})",
              f"scale({k * AT.SCALE:.5f})")


# --- the panel --------------------------------------------------------------
def build(states=3):
    d_style, d_body = _asset("sms_crown_dimer", "d")
    m_style, m_body = _asset("sms_crown_monomer", "m")
    cx = (CROWN_L + CROWN_R) / 2
    span = CROWN_R - CROWN_L

    # the bound state: the dimer where it was drawn, the target on its crowns
    dimer = _g(d_body, f"translate({-DIMER_VB[0]},{-DIMER_VB[1]})")
    tgt, tbox = _place_target(span, (CROWN_L, CROWN_Y + GRIP_SINK))
    bound = tgt + dimer
    b_box = (min(tbox[0], 0.0), tbox[1], max(tbox[2], DIMER_W), DIMER_H)

    o, y = [], 0.0

    def band(body, box):
        """Drop a state in, its top edge at the running y."""
        nonlocal y
        o.append(_g(body, f"translate(0,{y - box[1]:.2f})"))
        y += box[3] - box[1]

    if states >= 3:
        # apo: the target apart from two free chains, so that dimerization
        # reads as something the target causes rather than something assumed
        #
        # The target goes in the MIDDLE and a chain closes on it from either
        # side, the right one a mirror of the left. Three things in a row read
        # as a list; two things converging on a third read as an encounter,
        # which is what the next state resolves. Each chain also leans its
        # crown inward, by the same angle the bound state hangs them at.
        a_tgt, a_box = _place_target(span, (CROWN_L, CROWN_Y + GRIP_SINK))
        tcx = (a_box[0] + a_box[2]) / 2
        # Each chain waits out to the side with its crown high - level with the
        # target - so the SMS hangs low and outboard and leaves an open wedge
        # between the chain and the target. The motion arrow sweeps up through
        # that wedge, from beside the SMS to the hairpin tip the crown will
        # take. Starting it at the crown instead gave an arc too short to read.
        tip = (CROWN_L, CROWN_Y + GRIP_SINK)
        crown = (tip[0] - APPROACH[0], tip[1] - APPROACH[1])
        chain = _g(m_body,
                   f"translate({crown[0]:.2f},{crown[1]:.2f})",
                   f"rotate({LEAN:.1f})",
                   f"translate({-MONO_CROWN[0]:.2f},{-MONO_CROWN[1]:.2f})",
                   f"translate({-MONO_VB[0]},{-MONO_VB[1]})")
        # The arc runs from the CROWN to the hairpin head, because that pair is
        # the event: this binder takes that hairpin. Running it from the SMS
        # drew the right sweep between the wrong two things.
        chain += _curve_arrow((crown[0] + ARC_FROM[0], crown[1] + ARC_FROM[1]),
                              (tip[0] + ARC_TO[0], tip[1] + ARC_TO[1]),
                              ARC_BOW, 2.3)
        apo = a_tgt + chain + _g(chain, f"translate({2 * tcx:.2f},0)",
                                 "scale(-1,1)")
        span_x = max(a_box[2] - tcx,
                     tcx - (crown[0] - MONO_W * 0.62))
        apo_box = (tcx - span_x, min(a_box[1], crown[1] - 3.0),
                   tcx + span_x, max(a_box[3], crown[1] + MONO_H * 0.94))
        # centre the row on the column the other states sit in
        apo = _g(apo, f"translate({cx - tcx:.2f},0)")
        apo_box = (cx - span_x, apo_box[1], cx + span_x, apo_box[3])
        band(apo, apo_box)
        o.append(_step_arrow(cx, y + CLEAR))
        y += ARROW + CLEAR + DROP

    band(bound, b_box)
    o.append(_step_arrow(cx, y + CLEAR))

    # The arrow keeps CLEAR above it and the larger DROP below, because the
    # hairpin's own stem and junction loop rise into that gap and the arrow
    # points straight down the same centre line. Its height above the strand is the stem foot less the
    # top of the small junction loop, both from hairpin.py.
    k = span / (2 * (AT.SIDE_CX + AT.SIDE_R) * AT.SCALE) * REC_SCALE
    s = k * AT.SCALE
    rise = (AT.STEM_BOT - (AT.MID_CY - AT.MID_R)) * s
    cy = y + CLEAR + ARROW + DROP + rise
    rec, half = _recorder(cx, cy, k)
    o.append(rec)

    # the excised piece leaves the cut, so an arrow carries it out of the site
    ix = cx + half + 7.0
    o.append(f'<path d="M{cx + half + 2.6:.2f},{cy - rise * 0.52:.2f} '
             f'L{ix - 1.2:.2f},{cy - rise * 0.52:.2f}" fill="none" '
             f'stroke="{MUT}" stroke-width="0.6" stroke-linecap="round"/>')
    o.append(_arrowhead(ix, cy - rise * 0.52, 0.0, 2.2, MUT))
    o.append(_intron(ix + 2.0, cy - rise * 0.72, k))
    y = cy + 3.0
    half += 20.0

    xs = [b_box[0] - 2.0, b_box[2] + 2.0, cx - half - 12.0, cx + half + 26.0]
    if states >= 3:
        xs += [apo_box[0] - 2.0, apo_box[2] + 2.0]
    x0, x1, m = min(xs), max(xs), 1.5
    w, h = x1 - x0 + 2 * m, y + 2 * m
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w * 8:.0f}" '
            f'height="{h * 8:.0f}" viewBox="{x0 - m:.2f} {-m:.2f} '
            f'{w:.2f} {h:.2f}">' + d_style + m_style
            + "".join(o) + "</svg>")


def flatten(svg):
    """Fold the <style> rules onto the elements as plain attributes.

    The crown assets carry their colours in CSS classes. Illustrator opens
    that, but every shape sharing a class then moves as one, which is no use
    to someone who wants to nudge a single lobe. Writing the declarations out
    per element makes each shape independent. CSS beats a presentation
    attribute in SVG, so a rule overwrites an attribute already on the tag.
    """
    rules = {}
    for block in re.findall(r"<style>(.*?)</style>", svg, re.S):
        for sel, body in re.findall(r"\.([\w-]+)\s*\{([^}]*)\}", block):
            rules[sel] = [tuple(x.strip() for x in p.split(":", 1))
                          for p in body.split(";") if ":" in p]
    out = re.sub(r"<style>\s*</style>|<style>.*?</style>", "", svg, flags=re.S)

    def fix(m):
        raw = m.group(0)
        close = "/>" if raw.endswith("/>") else ">"
        tag = raw[:-len(close)].rstrip()
        tag = re.sub(r'\sclass="[^"]*"', "", tag)
        for name in m.group(1).split():
            for prop, val in rules.get(name, []):
                tag = re.sub(rf'\s{re.escape(prop)}="[^"]*"', "", tag)
                tag += f' {prop}="{val}"'
        return tag + close

    return re.sub(r'<[^>]*\sclass="([^"]*)"[^>]*>', fix, out)


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for name, n in (("mechanism_3state", 3), ("mechanism_2state", 2)):
        (outdir / f"{name}.svg").write_text(flatten(build(n)),
                                             encoding="utf-8")
        print(f"wrote {name}.svg")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else LIB)
