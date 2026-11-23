"""Fig 1d, PCR panel: the primer pair, and the two products it gives.

Four rows. Each recorder state - uncut, then spliced - with the primer pair on
it, and directly under each the amplicon that pair produces. The two amplicons
are left-aligned on the forward primer, so their right ends carry the result:
the spliced template is short by exactly the piece the recorder lost.

The target module is not amplified, so it is not drawn. The construct stops
after the SMS and the strand runs on as slack, which is where the reverse
primer sits.

**Primer convention.** Published figures draw a primer as a short arrow on the
template: its position is where it anneals, its direction is which way
polymerase extends - forward right, reverse left. The arrows are neutral,
because a primer annotates the construct rather than being one of its modules
(COLORS.md rule 3), and they sit BELOW the strand: above it the uncut
hairpin's own loops reach over the 5' lead, and a forward primer drawn there
crosses the recorder it is meant to sit upstream of. Sources in the README.

For a cassette exon this is the standard design: ONE pair in the flanking
sequence amplifies both states in one reaction, the retained-piece product
running longer. That is why both rows carry the pair at identical positions.

    python code/Fig1d_integrated_construct/make_pcr.py [outdir]
"""
import math
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))

import make_panel as MP                                      # noqa: E402

WORK = MP.WORK
NEUT = MP.NEUT

# --- the construct, trimmed to what the reaction actually sees -------------
SHIFT = 3.0            # the sensor module slid left, closing the dead lead
                       # between the recorder and the ATG
SMS_END = 120.6 - SHIFT
SLACK = 18.0           # strand past it, so the reverse primer has somewhere
END = SMS_END + SLACK
STRAND_Y = 30.28

# --- primers ---------------------------------------------------------------
PRIMER_Y = 35.3        # below the strand, in the lane the label ticks used
LABEL_Y = 41.6
FWD = (6.0, 17.5)      # anneals upstream of the recorder, extends right
REV = (END - 3.0, END - 14.5)   # in the slack, extends left
HEAD = 2.4

# --- amplicons -------------------------------------------------------------
AMP_Y = 50.5           # the product, directly under the template it came from
AMP_H = 3.4
# The bar's segments. The recorder's share is set here rather than taken from
# the hairpin drawing, because a hairpin is drawn wide and sits on a few bases
# of strand - its picture is not its sequence. The four add to 36, which is
# where the binder starts on the template, so everything downstream lines up
# with the row above.
FLANK_W, EXON_W, INTRON_W, LEAD_W = 3.0, 5.0, 17.0, 3.0
EXON = "#0678a6"        # theme 600, the hairpin's exons
INTRON = MP.palette.SYNBIO["intron"]   # the piece the recorder loses
# widths straight off the construct: binder, linker, kinase, RNase, slack
MODULES = [(22.28, "#015d82"), (4.95, NEUT["300"]), (31.47, "#3d8f85"),
           (19.90, "#ffc233"), (SLACK - 3.0, NEUT["300"])]
ROW_GAP = 60.0
# The rule is centred in the gap: from the upper amplicon's lower edge to the
# spliced hairpin's top, which is the highest thing in the lower half.
SPLICED_TOP = 12.56

# --- the ratio, gathered to the right ---------------------------------------
BR_X = 151.0           # the bracket's spine, clear of the 3' captions
BR_ARM = 3.5
BR_GAP = 9.0           # spine to the fraction
BR_TOP, BR_BOT = 13.0, 115.0
MINI_K = 0.20          # the fraction's bars: WIDTH only (see _ratio)
MINI_H = 2.8
PLUS_W = 7.5           # room for the plus, which sits between two bars
BASE = 6.4             # bar underside to the caption baseline below it
TITLE_GAP = 3.5        # box top to the title's baseline above it
TITLE_PT = 7.0         # the title carries the panel, so it is set
                       # larger than the captions inside the fraction
RATIO_PAD = 4.5        # the card's inset around the whole block

# the annotated construct's furniture, which this panel does not carry
STRIP = ("Recorder | bracket and label",
         "Sensor module | bracket and label",
         "Target module | bracket and label",
         "ATG site", "Start codon",
         "Label tick | Binder", "Label | Binder",
         "Label tick | SMS", "Label | SMS",
         "RNA scaffold")


def _bare(body):
    """Drop the captions, brackets and the target, and pull the strand in."""
    for name in STRIP:
        # the tail has to stay inside the element: an unanchored .*?</text>
        # lets a self-closing <path> run on and eat the next <text>
        pat = (rf'<g[^>]*data-name="{re.escape(name)}".*?</g>|'
               rf'<(?:path|text)[^>]*data-name="{re.escape(name)}"'
               rf'[^>]*?(?:/>|>.*?</text>)')
        body, n = re.subn(pat, "", body, count=1, flags=re.S)
        if not n:
            raise SystemExit(f"nothing stripped for {name!r}")
    body, n = re.subn(r'(data-name="RNA baseline" d="M41\.4 30\.28 L)[\d.]+',
                      lambda m: m.group(1) + f"{END:g}", body)
    assert n, "the baseline moved - update the pattern"
    body, n = re.subn(r'(data-name="RNA 3 prime" x=")[\d.]+',
                      lambda m: m.group(1) + f"{END + 2.0:g}", body)
    assert n, "the 3' caption moved - update the pattern"
    body, n = re.subn(r'(<g id="Sensor_module"[^>]*)>',
                      rf'\1 transform="translate({-SHIFT},0)">', body)
    assert n, "the sensor module group moved - update the pattern"
    return body


def _primer(x0, x1, label, y=PRIMER_Y, col=NEUT["700"]):
    d = 1.0 if x1 > x0 else -1.0
    return (f'<path data-name="Primer {label}" d="M{x0:.2f},{y:.2f} '
            f'L{x1 - d * HEAD * 0.8:.2f},{y:.2f}" fill="none" stroke="{col}" '
            f'stroke-width="0.7" stroke-linecap="round"/>'
            + MP._head(x1, y, 0.0 if d > 0 else math.pi, HEAD, col)
            + f'<text data-name="Label | {label}" text-anchor="middle" '
              f'x="{(x0 + x1) / 2:.2f}" y="{LABEL_Y}" class="lbl">{label}</text>')


def _outline(x0, x1, y, r=1.7, h=None):
    h = AMP_H if h is None else h
    r = min(r, h / 2)
    return (f'M{x0 + r:.3f} {y:.3f} L{x1 - r:.3f} {y:.3f} '
            f'A{r} {r} 0 0 1 {x1 - r:.3f} {y + h:.3f} '
            f'L{x0 + r:.3f} {y + h:.3f} '
            f'A{r} {r} 0 0 1 {x0 + r:.3f} {y:.3f} Z')


def _segments(spliced):
    """The product's colour map, 5' to 3'.

    The recorder is the point. Uncut it reads exon - intron - exon, the red
    flanked by blue on both sides; splicing joins the two exons, so the
    spliced product carries ONE continuous blue and ends short by the red's
    width. Everything downstream keeps the template's widths, so the protein
    part of the bar lines up with the row above it.
    """
    segs = ([(FLANK_W, NEUT["300"]), (2 * EXON_W, EXON)] if spliced else
            [(FLANK_W, NEUT["300"]), (EXON_W, EXON),
             (INTRON_W, INTRON), (EXON_W, EXON)])
    return segs + [(LEAD_W, NEUT["300"])] + MODULES


def amp_width(spliced):
    return sum(w for w, _ in _segments(spliced))


def _amplicon(x0, y, spliced, uid, k=1.0, h=AMP_H, r=1.7):
    """One product as a segmented bar, at full size or shrunk for the fraction.

    `k` scales width only. The height is given separately, because scaling a
    bar uniformly down to fraction size turns it into a hairline and the
    segments stop reading.
    """
    x = x0
    rects = []
    for w, col in _segments(spliced):
        rects.append(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w * k:.3f}" '
                     f'height="{h:.3f}" fill="{col}"/>')
        x += w * k
    d = _outline(x0, x, y, r, h)
    return (f'<clipPath id="amp{uid}"><path d="{d}"/></clipPath>'
            f'<g data-name="Amplicon" clip-path="url(#amp{uid})">'
            + "".join(rects) + "</g>"
            + f'<path d="{d}" fill="none" stroke="{NEUT["700"]}" '
              f'stroke-width="0.4"/>')


def _ratio(uid):
    """The bracket gathering both states, and the ratio written out in them.

    The fraction uses the amplicon glyphs rather than the words "spliced" and
    "unspliced", because the reader has just met those two bars full size an
    inch to the left: the short one and the long one with red in it. Writing
    the formula in the same objects means the definition needs no key.

    MINI_K scales width only. Shrinking the bars uniformly would reduce them
    to hairlines and the segments - which is the whole identification - would
    stop reading.
    """
    o = []
    mid = (BR_TOP + BR_BOT) / 2
    # square bracket: arms toward the states it gathers, nub toward the result
    o.append(f'<path data-name="Ratio bracket" '
             f'd="M{BR_X - BR_ARM:.2f},{BR_TOP:.2f} L{BR_X:.2f},{BR_TOP:.2f} '
             f'L{BR_X:.2f},{BR_BOT:.2f} L{BR_X - BR_ARM:.2f},{BR_BOT:.2f} '
             f'M{BR_X:.2f},{mid:.2f} L{BR_X + BR_ARM:.2f},{mid:.2f}" '
             f'fill="none" stroke="{MP.INK if hasattr(MP, "INK") else "#1d2428"}" '
             f'stroke-width="0.4" stroke-linecap="round" '
             f'stroke-linejoin="round"/>')

    w_spl = amp_width(True) * MINI_K
    w_unc = amp_width(False) * MINI_K
    den_w = w_spl + PLUS_W + w_unc
    x0 = BR_X + BR_ARM + BR_GAP
    rule_y = mid

    # Every term is captioned, so each one is a bar with its name under it and
    # the stack above the rule has to make room for a caption too. BASE is the
    # drop from a bar's underside to the caption's baseline: it is the
    # ascender, not the baseline, that approaches the bar from below, which is
    # why it is larger than the TITLE gap, where a baseline sits under text.
    num_lab_y = rule_y - 2.0
    num_y = num_lab_y - BASE - MINI_H
    o.append(_amplicon(x0 + (den_w - w_spl) / 2, num_y,
                       True, uid, MINI_K, MINI_H, 0.9))
    o.append(f'<text data-name="Species | spliced (numerator)" '
             f'text-anchor="middle" x="{x0 + den_w / 2:.2f}" '
             f'y="{num_lab_y:.2f}" class="lbl">spliced</text>')
    o.append(f'<path data-name="Fraction rule" d="M{x0:.2f},{rule_y:.2f} '
             f'L{x0 + den_w:.2f},{rule_y:.2f}" fill="none" stroke="#1d2428" '
             f'stroke-width="0.5" stroke-linecap="round"/>')
    den_y = rule_y + 3.0
    o.append(_amplicon(x0, den_y, True, uid + 1, MINI_K, MINI_H, 0.9))
    o.append(_amplicon(x0 + w_spl + PLUS_W, den_y, False, uid + 2,
                       MINI_K, MINI_H, 0.9))
    # the plus, drawn rather than set, so it keeps its weight at this size
    px, py, a = x0 + w_spl + PLUS_W / 2, den_y + MINI_H / 2, 1.3
    o.append(f'<path data-name="Plus" d="M{px - a:.2f},{py:.2f} '
             f'L{px + a:.2f},{py:.2f} M{px:.2f},{py - a:.2f} '
             f'L{px:.2f},{py + a:.2f}" fill="none" stroke="#1d2428" '
             f'stroke-width="0.5" stroke-linecap="round"/>')
    lab_y = den_y + MINI_H + BASE
    for cx, txt in ((x0 + w_spl / 2, "spliced"),
                    (x0 + w_spl + PLUS_W + w_unc / 2, "unspliced")):
        o.append(f'<text data-name="Species | {txt}" text-anchor="middle" '
                 f'x="{cx:.2f}" y="{lab_y:.2f}" class="lbl">{txt}</text>')
    # The box holds the fraction only; its title stands above it, naming the
    # box the way a panel letter names a panel. It takes the panel's own
    # ground colour, so the outline alone does the separating - a lighter or
    # darker fill would read as a second surface sitting on the first.
    bx0, bx1 = x0 - RATIO_PAD, x0 + den_w + RATIO_PAD
    by0, by1 = num_y - RATIO_PAD, lab_y + 2.0
    card = (f'<rect data-name="Ratio card" x="{bx0:.2f}" y="{by0:.2f}" '
            f'width="{bx1 - bx0:.2f}" height="{by1 - by0:.2f}" rx="2.5" '
            f'fill="{NEUT["100"]}" stroke="{NEUT["300"]}" stroke-width="0.4"/>')
    o.append(f'<text data-name="Ratio title" text-anchor="middle" '
             f'x="{x0 + den_w / 2:.2f}" y="{by0 - TITLE_GAP:.2f}" class="lbl" '
             f'style="font-size:{TITLE_PT}px">Splicing ratio</text>')
    return card + "".join(o), bx1


def build():
    style, raw = MP._construct()
    uncut = _bare(raw)
    spliced = re.sub(r'<g data-name="Excised intron".*?</g>', "",
                     _bare(MP._splice(MP._construct()[1])), flags=re.S)
    pair = _primer(*FWD, "Fwd") + _primer(*REV, "Rev")

    rows = (f'<g data-name="Uncut">{uncut}{pair}{_amplicon(FWD[0], AMP_Y, False, 1)}</g>'
            f'<g data-name="Spliced" transform="translate(0,{ROW_GAP})">'
            f'{spliced}{pair}{_amplicon(FWD[0], AMP_Y, True, 2)}</g>')
    m, top, bot = 4.0, 11.0, AMP_Y + AMP_H + 2.0
    w = END + 12.0
    h = ROW_GAP + (bot - top) + 2 * m
    pad = 3.0
    x0, y0 = -m + 0.6, top - m + 0.6
    bw, bh = w + 2 * m - 1.2, h - 1.2
    # a container behind everything, and a hairline rule between the two
    # states, so the pair reads as one comparison rather than two drawings
    box = (f'<rect data-name="Panel ground" x="{x0:.2f}" y="{y0:.2f}" '
           f'width="{bw:.2f}" height="{bh:.2f}" rx="3" '
           f'fill="{NEUT["100"]}" stroke="{NEUT["300"]}" stroke-width="0.35"/>')
    dy = ((AMP_Y + AMP_H) + (ROW_GAP + SPLICED_TOP)) / 2
    rule = (f'<path data-name="State divider" d="M{x0 + pad:.2f},{dy:.2f} '
            f'L{x0 + bw - pad:.2f},{dy:.2f}" fill="none" '
            f'stroke="{NEUT["400"]}" stroke-width="0.35" '
            f'stroke-dasharray="2 1.6" stroke-linecap="butt"/>')
    ratio, ratio_x1 = _ratio(10)
    w = ratio_x1 + 5.0 - (-m)
    bw = w + 2 * m - 1.2
    box = (f'<rect data-name="Panel ground" x="{x0:.2f}" y="{y0:.2f}" '
           f'width="{bw:.2f}" height="{bh:.2f}" rx="3" '
           f'fill="{NEUT["100"]}" stroke="{NEUT["300"]}" stroke-width="0.35"/>')
    rule = (f'<path data-name="State divider" d="M{x0 + pad:.2f},{dy:.2f} '
            f'L{BR_X - BR_ARM - 2.0:.2f},{dy:.2f}" fill="none" '
            f'stroke="{NEUT["400"]}" stroke-width="0.35" '
            f'stroke-dasharray="2 1.6" stroke-linecap="butt"/>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{(w + 2 * m) * 6:.0f}" height="{h * 6:.0f}" '
            f'viewBox="{-m} {top - m} {w + 2 * m:.2f} {h:.2f}">'
            + style + box + rule + rows + ratio + "</svg>")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "Fig1d_panel_pcr.svg"
    p.write_text(build(), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else WORK)
