"""Fig 4a - the screen, end to end, in one horizontal panel.

Three drawings this repo already emits, placed side by side. NOTHING is
redrawn here; this script only reads, rescales and positions:

  * output/Fig4a_library_cassette/library_cassette.svg - the donor cassette
    with the oligo pool going into the binder and the target-scaffold
    call-out;
  * output/Fig4b_cell_populations/cell_populations.svg - the four populations
    the library lands in, which this panel puts in a box of their own;
  * output/readout_workflow/readout_workflow.svg - RNA out of the cells
    through to filtered reads.

Plus the screen's last step, drawn here: the hit gate, splicing ratio against
fold change. That one is NOT a cartoon - its points are Fig 4c's own rows,
read from the same CSV under the same gate, so the cloud's shape and where
the 32 hits sit are the screen's. Only the out-of-gate cloud is thinned, and
deterministically.

**Every label in the panel prints at the same size.** Each source is scaled by
`5.5 / its own label size` rather than by eye, so a piece can be re-emitted at
a different internal scale and still land right. The cassette and the
populations are already drawn at 5.5; the workflow is drawn at 7 pt in its own
85 mm frame, which is 13.65 of its units, so it comes in at about 0.40. That
one number is also why the workflow's circles end up about the size of a cell
cluster's zoom - the match follows from the type, it was not set by eye.

**The plot frame sits in the L's own empty corner**, and where that corner
starts is read off the workflow script itself (its column labels' right edge,
its row labels' bottom edge) rather than typed in here. Move a step or
re-word a label there and the frame still clears them.

**The readout workflow is switched OFF** (Leo, 2026-11-13: "remove this whole
part for now"). With WITH_WORKFLOW = False the panel is the left column
(cassette over the populations box) and the gate plot beside it, PLOT_ASPECT
wide for its height, the same height as the column. Set it True and the
original L-shaped layout comes back unchanged - nothing of it was deleted.

Source order in run_all.sh: Fig4b runs AFTER this, so a changed Fig 4b needs
its own script run first, or a second pass.

    python code/Fig4a_screen_overview/make_figure.py [--check]
"""
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette, style                                      # noqa: E402
from common.paths import ROOT, OUTPUT, fig_output               # noqa: E402

OUT = fig_output("Fig4a_screen_overview")

CASSETTE = OUTPUT / "Fig4a_library_cassette" / "library_cassette.svg"
POPS = OUTPUT / "Fig4b_cell_populations" / "cell_populations.svg"
READOUT = OUTPUT / "readout_workflow" / "readout_workflow.svg"

NEUT = palette.NEUTRAL
FRAME = palette.SYNBIO["outline"]       # the cassette's own call-out stroke
FONT = "Arial"

LABEL = 5.5             # the panel's one label size, in panel units
# 2026-11-15 (Leo, font audit): the axis title was written with a literal
# U+2082 SUBSCRIPT TWO. Arial has no subscript digits, so rsvg-convert silently
# fell back to HELVETICA for that one glyph - the only non-Arial face in the
# whole book besides the mathtext italics. Drawn as a real tspan instead, so
# the 2 stays Arial and the subscript is typeset rather than borrowed from a
# codepoint. The dy is undone on the next run or the rest of the line rides low.
def sub2(inner):
    """log, a typeset subscript 2, then (inner)."""
    return ('log<tspan dy="{:.2f}" font-size="{:.2f}px">2</tspan>'
            '<tspan dy="{:.2f}">({})</tspan>').format(
                LABEL * 0.20, LABEL * 0.72, -LABEL * 0.20, inner)


SUB2_LABEL = sub2("Fold Change")
RO_LABEL = 13.65        # what the workflow draws its 7 pt labels at

# ONE gap between blocks, everywhere: plasmid down to the populations box,
# and that box across to the workflow. Leo spaced those by hand at 16.2 and
# 20.0 and asked for them even, so they are both the larger of the two - the
# smaller would not hold the integration arrow and its label.
GAP = 20.0
BOX_PAD = 5.0           # inside the populations box
CORNER_PAD = 7.0        # the L's corner to the plot frame. NOT the block gap:
                        # this is clearance from a label to a frame tucked
                        # inside the L, and at 20 it would eat 13 units of
                        # the plot to line up with a space nothing abuts.
MARGIN = 3.0            # a hair of air, so strokes are not clipped
LINE_H = 7.6            # between stacked label lines

# Plasmid into the cells. INT_X is where the arrow drops; the label sits to
# its right because the gap is 20 units and two stacked lines leave an arrow
# no room to be an arrow in.
INT_X = 40.0
INTEGRATION = ["Bxb1 landing pad", "transfection and integration"]
CONN = palette.NEUTRAL["400"]           # the workflow's own connector grey
SW_CONN = 0.67                          # its stroke, in panel units
HEAD = 2.0
WITH_WORKFLOW = False                   # see the docstring; True restores the L
# The gate plot left panel a on 2026-11-14 (Leo: split 4a into "the whole cell
# and integration part" and a screen-output panel with two graphs). With
# WITH_GATE False this script emits TWO files: Fig4a_screen_overview.svg is
# the left column alone, and Fig4b_screen_gates.svg is the real gate beside
# the negative-against-negative null. Set it True and panel a is what it was.
WITH_GATE = False
PLOT_ASPECT = 1.2                       # gate plot width per height, workflow off


# --- reading a source ------------------------------------------------------
def _parts(path):
    """A source's viewBox, its <style> text and its body."""
    t = path.read_text()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]
    m = re.search(r"<style>(.*?)</style>", t, re.S)
    style = m.group(1) if m else ""
    body = re.sub(r"<style>.*?</style>", "", t, flags=re.S)
    body = body[body.index(">", body.index("<svg")) + 1:body.rindex("</svg>")]
    return vb, style, body


def _namespace(style, body, tag):
    """Prefix every class this source defines, so two sources can both own a
    `.lbl` without one silently restyling the other."""
    names = set(re.findall(r"\.([A-Za-z][\w-]*)\s*\{", style))
    for n in sorted(names, key=len, reverse=True):
        style = re.sub(rf"\.{re.escape(n)}(?=\s*[,{{])", f".{tag}-{n}", style)
        body = body.replace(f'class="{n}"', f'class="{tag}-{n}"')
        body = body.replace(f'class="{n} ', f'class="{tag}-{n} ')
    return style, body


def conn_arrow(x1, y1, x2, y2):
    """The workflow's own connector, at panel scale. The two links this panel
    adds are process steps like the workflow's, so they are drawn as its
    arrows are and not invented."""
    import math
    a = math.atan2(y2 - y1, x2 - x1)
    c, s = math.cos(a), math.sin(a)
    bx, by = x2 - HEAD * c, y2 - HEAD * s
    wg = HEAD * 0.46
    return (f'<path d="M{x1:.2f},{y1:.2f} L{bx:.2f},{by:.2f}" fill="none" '
            f'stroke="{CONN}" stroke-width="{SW_CONN}" stroke-linecap="round"/>'
            f'<path d="M{x2:.2f},{y2:.2f} L{bx - wg * s:.2f},{by + wg * c:.2f} '
            f'L{bx + wg * s:.2f},{by - wg * c:.2f} Z" fill="{CONN}"/>')


def ink_x0(path):
    """The leftmost drawn x in a source's own units.

    The cassette's viewBox is padded well left of its plasmid, so laying the
    left column out on viewBoxes leaves the plasmid inset from the
    populations box below it by about a fifth of the box's width. Text is not
    measured: in this source nothing textual reaches past the plasmid's own
    left edge, and `_check()` re-tests that the two left edges still agree.
    """
    _, _, body = _parts(path)
    xs = []
    for d in re.findall(r'\sd="([^"]+)"', body):
        xs += [float(v) for v in re.findall(r"-?[\d.]+", d)][0::2]
    xs += [float(v) for v in re.findall(r'<rect[^>]*\sx="(-?[\d.]+)"', body)]
    return min(xs)


def piece(path, tag, k, x, y, drop_bg=False, crop=None):
    """One source, scaled by k and placed with its top left ink at (x, y).

    `crop` is (x0, y0, w, h) in the source's own units, used in place of its
    viewBox. The workflow is cropped: it is a panel in its own right and
    carries a panel's margin, which inside this figure is just a hole between
    it and the populations box.
    """
    vb, style, body = _parts(path)
    if drop_bg:        # the workflow paints its own canvas white
        body = re.sub(r'<rect width="[\d.]+" height="[\d.]+" fill="#\w+"/>',
                      "", body, count=1)
    style, body = _namespace(style, body, tag)
    box = crop or vb
    g = (f'<g data-name="{tag}" transform="translate({x:.2f},{y:.2f}) '
         f'scale({k:.5f}) translate({-box[0]:.2f},{-box[1]:.2f})">{body}</g>')
    return style, g, box[2] * k, box[3] * k


# --- the workflow's own geometry, so the plot frame can clear it -----------
def _readout_module():
    spec = importlib.util.spec_from_file_location(
        "ro_workflow", ROOT / "code" / "readout_workflow" / "make_figure.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


RO = _readout_module()
K_RO = LABEL / RO_LABEL

# Where the L stops being an L: right of the column's longest label, below the
# row's deepest one. Both come from the workflow script, not from this one.
COL_LAB_R = (RO.COL_X + RO.EDGE + RO.LAB_DX
             + max(RO.label_w(lines) for i, (_, lines) in enumerate(RO.STEPS)
                   if RO.beside(i)))
ROW_LAB_B = max(RO.COL_TOP + RO.EDGE + RO.LAB_DY + len(lines) * RO.LINE_H
                for i, (_, lines) in enumerate(RO.STEPS) if not RO.beside(i))

# The workflow's ink, in its own units - what is left once its panel margin
# goes. Its row labels are centred and none of them reaches past the last
# circle, so the right edge is a container edge; the deepest thing is the
# bottom circle, not a label, since the column wears its labels at its side.
_RO_POS = RO.layout()
RO_INK = (RO.COL_X - RO.EDGE,
          RO.COL_TOP - RO.EDGE,
          max(cx for cx, _ in _RO_POS) + RO.EDGE - (RO.COL_X - RO.EDGE),
          max(cy for _, cy in _RO_POS) + RO.EDGE - (RO.COL_TOP - RO.EDGE))


# --- the gate, as the panel's last step -----------------------------------
# A reduction of Fig 4c, not a cartoon of one: the points are that figure's
# own rows, read from the same CSV under the same gate, so the cloud's shape
# and where the hits sit are the screen's and not a drawing. Only the
# out-of-gate cloud is thinned, deterministically, and the thinning is said
# out loud rather than left to look like the whole screen.
VOLC = ROOT / "inputs" / "screen_points.csv"
R_GATE, S_GATE, FC_GATE = 75, 0.07, 10.0        # Fig 4c's gate, unchanged
SPL_CAP, Y_TOP = 0.25, 8.0                      # Fig 4c's display limits
BG_N = 2600                                     # out-of-gate points kept
# Fig 4c's own encoding, kept so the two read as one plot: the two screen
# comparisons on the theme's dark/light pair, everything out of gate in ONE
# neutral grey rather than as faded series colours.
ORDER = ["POS_A_vs_NEG_A", "POS_B_vs_NEG_B"]
HIT_COL = {"POS_A_vs_NEG_A": palette.THEME["700"],
           "POS_B_vs_NEG_B": palette.THEME["400"]}
# 2026-11-23: drawn labels follow the figure convention (PosA/NegA, as the
# supplementary flow panels use); the raw column values stay the dict keys.
HIT_LAB = {"POS_A_vs_NEG_A": "PosA / NegA", "POS_B_vs_NEG_B": "PosB / NegB"}
XTICKS = (0, 0.05, 0.10, 0.15, 0.20, 0.25)
YTICKS = (0, 2, 4, 6, 8)
PAD = (20.0, 3.0, 12.0, 14.0)                   # left, right, top, bottom
                                                # top carries the legend row
LEG_GAP = 7.0                                   # swatch to the next entry


def _points():
    """Fig 4c's points: (x, y, comparison, in-gate), under its own reads cut."""
    import csv
    import math
    out = []
    with open(VOLC) as f:
        for r in csv.DictReader(f):
            fc = r["fold_change_raw"]
            if not fc:
                continue
            if float(r["reads"]) < R_GATE or float(r["neg_reads"]) < R_GATE:
                continue
            fc, spl = float(fc), float(r["splicing_ratio"])
            lfc = max(math.log2(fc), 0.0) if fc > 0 else 0.0
            out.append((min(spl, SPL_CAP), min(lfc, Y_TOP), r["comparison"],
                        spl >= S_GATE and fc >= FC_GATE))
    return out


# The null: the same gate run on one negative against the other. NO LONGER
# DRAWN (2026-11-16) - signal_plot above replaced it - but kept because the
# reason it was limited is the reason it was replaced.
# WHY IT IS BUILT HERE AND NOT READ. screen_points.csv carries the POSITIVE
# side's splicing ratio and the POS/NEG fold change, not the negatives' own
# ratios, so each negative's ratio is recovered as splicing_ratio /
# fold_change_raw. That division needs a non-zero numerator, so the null can
# only be built for designs that spliced in BOTH positive conditions - 2,557
# of the 9,094 with reads on both negatives. It is therefore a null over the
# screen's most active designs, which is the subset where a leaky gate would
# show, not a null over the whole library. The per-condition table that would
# give all of them (matched_comparisons.csv, or the report it was extracted
# from) is not on the NAS or on Euler.
NEG_ORDER = ["A_over_B", "B_over_A"]
NEG_LAB = {"A_over_B": "NegA / NegB", "B_over_A": "NegB / NegA"}


def _neg_points():
    """(x, y, direction, in-gate) for NEG_A against NEG_B and the reverse."""
    import csv
    import math
    by = {}
    with open(VOLC) as f:
        for r in csv.DictReader(f):
            fc = r["fold_change_raw"]
            if not fc or float(fc) <= 0 or float(r["neg_reads"]) < R_GATE:
                continue
            # the negative's own splicing ratio, recovered
            by.setdefault(r["transcript"].split("::")[-1], {})[r["comparison"]] = \
                float(r["splicing_ratio"]) / float(fc)
    out = []
    for d, c in by.items():
        na = c.get("POS_A_vs_NEG_A")
        nb = c.get("POS_B_vs_NEG_B")
        if na is None or nb is None:
            continue
        for key, num, den in (("A_over_B", na, nb), ("B_over_A", nb, na)):
            fc = num / den
            out.append((min(num, SPL_CAP), min(max(math.log2(fc), 0.0), Y_TOP),
                        key, num >= S_GATE and fc >= FC_GATE))
    return out


def gate_plot(x, y, w, h, pts=None, order=None, col=None, lab=None,
              xlab="Splicing ratio (fraction spliced)", ylab=None,
              fclab="FC 10"):
    """Fig 4c, in the panel: splicing ratio against fold change, with gates.

    Same data, same gate, same encoding and the same legend as Fig 4c, so a
    reader meets one plot twice rather than two. What is NOT the same is the
    type: Fig 4c sets its own 9-15 pt for a figure of its own, and dropping
    that file in here whole would put three further type sizes into a panel
    whose whole premise is one. So it is redrawn in panel units at the
    panel's 5.5, and landscape, since the corner it fills is about 2:1 while
    Fig 4c is square.
    """
    import math
    import random
    pl, pr, pt, pb = PAD
    fx, fy = x + pl, y + pt
    fw, fh = w - pl - pr, h - pt - pb

    def px(v):
        return fx + v / SPL_CAP * fw

    def py(v):
        return fy + fh - v / Y_TOP * fh

    def T(tx, ty, s, col=palette.INK, anch="middle", rot=None):
        r = f' transform="rotate(-90 {tx:.2f} {ty:.2f})"' if rot else ""
        return (f'<text x="{tx:.2f}" y="{ty:.2f}" font-family="{FONT}" '
                f'font-size="{LABEL}px" fill="{col}" text-anchor="{anch}"'
                f'{r}>{s}</text>')

    pts = _points() if pts is None else pts
    order = ORDER if order is None else order
    col = HIT_COL if col is None else col
    lab = HIT_LAB if lab is None else lab
    hits = [p for p in pts if p[3]]
    out_of = [p for p in pts if not p[3]]
    # the SAME number of background marks in both plots, so the two read as
    # one plot twice and a difference in density is a difference in the data
    rest = (random.Random(0).sample(out_of, BG_N) if len(out_of) > BG_N
            else out_of)

    o = [f'<g data-name="Hit gate">']
    o += [f'<circle cx="{px(a):.2f}" cy="{py(b):.2f}" r="0.45" '
          f'fill="{NEUT["400"]}" fill-opacity="0.4"/>' for a, b, _, _ in rest]
    gx, gy = px(S_GATE), py(math.log2(FC_GATE))
    o.append(f'<path d="M{gx:.2f},{fy:.2f} L{gx:.2f},{fy + fh:.2f} '
             f'M{fx:.2f},{gy:.2f} L{fx + fw:.2f},{gy:.2f}" fill="none" '
             f'stroke="{palette.INK}" stroke-width="0.35" '
             f'stroke-dasharray="1.4 1.4"/>')
    for comp in order:          # A under B, as Fig 4c draws them
        o += [f'<circle cx="{px(a):.2f}" cy="{py(b):.2f}" r="0.9" '
              f'fill="{col[comp]}"/>' for a, b, c, _ in hits if c == comp]

    # frame and ticks, over the marks as in Fig 4c
    o.append(f'<rect x="{fx:.2f}" y="{fy:.2f}" width="{fw:.2f}" '
             f'height="{fh:.2f}" fill="none" stroke="{palette.INK}" '
             f'stroke-width="0.5"/>')
    tk = []
    for v in XTICKS:
        tk.append(f'M{px(v):.2f},{fy + fh:.2f} l0,1.6')
        o.append(T(px(v), fy + fh + LABEL + 1.0, f"{v:g}"))
    for v in YTICKS:
        tk.append(f'M{fx:.2f},{py(v):.2f} l-1.6,0')
        o.append(T(fx - 2.8, py(v) + LABEL * 0.35, f"{v:g}", anch="end"))
    o.append(f'<path d="{" ".join(tk)}" fill="none" stroke="{palette.INK}" '
             f'stroke-width="0.4"/>')

    # The gates are named on the lines themselves. Fig 4c leaves them to its
    # caption; in a panel there is no caption to leave them to, and both sit
    # in the empty bottom-right quadrant rather than over any mark.
    o += [T(gx + 2.0, fy + fh - 2.2, "0.07", col=NEUT["600"], anch="start"),
          T(fx + fw - 2.0, gy + LABEL + 1.2, fclab, col=NEUT["600"],
            anch="end")]

    # axis titles, Fig 4c's wording
    o += [T(fx + fw / 2, y + h, xlab),
          T(x + LABEL, fy + fh / 2, ylab or SUB2_LABEL, rot=True)]

    # Legend above the frame, one row, with each comparison's hit count.
    # Fig 4c's third entry, "outside gate", is dropped: it would push the row
    # past the frame, and one neutral grey against two coloured series needs
    # no swatch to be read as everything else.
    # It starts at the box's own left edge, not the frame's: the row is 12
    # units wider than the frame, and the alternative was to lengthen the
    # workflow's arrows until the plot grew to fit a legend, which is the
    # layout driving the wrong thing.
    lx = x + 2.0
    for comp in order:
        n = sum(1 for _, _, c, _ in hits if c == comp)
        o.append(f'<circle cx="{lx + 1.3:.2f}" cy="{y + 5.0:.2f}" r="1.3" '
                 f'fill="{col[comp]}"/>')
        o.append(T(lx + 4.4, y + 5.0 + LABEL * 0.35,
                   f"{lab[comp]} ({n})", anch="start"))
        lx += LEG_GAP + style.text_w(f"{lab[comp]} ({n})", LABEL)
    assert lx <= x + w + 0.5, (
        f"the legend is {lx - x - w:.0f} units wider than the plot box")
    o.append("</g>")
    return "".join(o)


def build():
    styles, body = [], []

    # The left column has two containers - the plasmid and the populations
    # box - and they share a left edge. The cassette is placed by its ink,
    # not by its viewBox, which is padded well left of the plasmid.
    vb_cas = _parts(CASSETTE)[0]
    cas_x = -(ink_x0(CASSETTE) - vb_cas[0])
    box_w = _parts(POPS)[0][2] + 2 * BOX_PAD
    left_w = max(cas_x + vb_cas[2], box_w)

    s, g_cas, _, h_cas = piece(CASSETTE, "cas", 1.0, cas_x, 0.0)
    styles.append(s)

    s, g_pop, w_pop, h_pop = piece(POPS, "pop", 1.0, BOX_PAD,
                                   h_cas + GAP + BOX_PAD)
    styles.append(s)
    box_h = h_pop + 2 * BOX_PAD
    pop_box = (f'<rect data-name="Cell populations" x="0" '
               f'y="{h_cas + GAP:.2f}" width="{box_w:.2f}" '
               f'height="{box_h:.2f}" rx="2" fill="none" stroke="{FRAME}" '
               f'stroke-width="0.4"/>')

    left_h = h_cas + GAP + box_h

    # Plasmid down into the cells. The label sits beside the arrow rather
    # than under it: the gap is 26 units and two lines of type would leave
    # the arrow no length to be an arrow in.
    gap_mid = h_cas + GAP / 2
    y0 = gap_mid - (len(INTEGRATION) - 1) * LINE_H / 2 + LABEL * 0.35
    step = [conn_arrow(INT_X, h_cas + 2.0, INT_X, h_cas + GAP - 1.5)]
    step += [f'<text x="{INT_X + 6.0:.2f}" y="{y0 + j * LINE_H:.2f}" '
             f'font-family="{FONT}" font-size="{LABEL}px" '
             f'fill="{palette.INK}">{ln}</text>'
             for j, ln in enumerate(INTEGRATION)]
    int_lab_r = INT_X + 6.0 + max(style.text_w(l, LABEL) for l in INTEGRATION)
    cas_ink = cas_x + ink_x0(CASSETTE) - vb_cas[0]

    if not WITH_WORKFLOW and not WITH_GATE:
        # Panel a is the left column alone: the cassette, the populations and
        # the arrow between them. The gate moved to its own panel.
        W, H = left_w, left_h
        geo = dict(plot_w=None, plot_h=None, box_w=box_w, h_cas=h_cas,
                   box_y=h_cas + GAP, int_lab_r=int_lab_r, left_w=left_w,
                   cas_x=cas_ink, box_x=0.0, ro_x=None, ro_y=None,
                   plot_x=None, plot_y=None)
        body += [g_cas, pop_box, g_pop] + step
    elif not WITH_WORKFLOW:
        # The gate plot beside the left column, as tall as it, nothing between.
        H = left_h
        plot_x, plot_y, plot_h = left_w + GAP, 0.0, H
        plot_w = PLOT_ASPECT * H
        plot = gate_plot(plot_x, plot_y, plot_w, plot_h)
        W = plot_x + plot_w
        geo = dict(plot_w=plot_w, plot_h=plot_h, box_w=box_w, h_cas=h_cas,
                   box_y=h_cas + GAP, int_lab_r=int_lab_r, left_w=left_w,
                   cas_x=cas_ink, box_x=0.0, ro_x=None, ro_y=None,
                   plot_x=plot_x, plot_y=plot_y)
        body += [g_cas, pop_box, g_pop, plot] + step
    else:
        ro_x = left_w + GAP
        H = max(left_h, RO_INK[3] * K_RO)
        ro_y = (H - RO_INK[3] * K_RO) / 2       # centred against the left column
        s, g_ro, w_ro, h_ro = piece(READOUT, "ro", K_RO, ro_x, ro_y,
                                    drop_bg=True, crop=RO_INK)
        styles.append(s)

        # The plot fills the L's corner exactly: left edge clear of the column's
        # labels, top clear of the row's, right edge flush with the workflow's.
        plot_x = ro_x + (COL_LAB_R - RO_INK[0]) * K_RO + CORNER_PAD
        plot_y = ro_y + (ROW_LAB_B - RO_INK[1]) * K_RO + CORNER_PAD
        plot_h = H - plot_y
        plot_w = ro_x + RO_INK[2] * K_RO - plot_x
        plot = gate_plot(plot_x, plot_y, plot_w, plot_h)

        # The cell box into the first workflow step - the workflow's OWN
        # connector, level with the step it feeds.
        rna_cx, rna_cy = _RO_POS[0]
        rna = (ro_x + (rna_cx - RO_INK[0]) * K_RO,
               ro_y + (rna_cy - RO_INK[1]) * K_RO)
        step.append(conn_arrow(box_w + 3.0, rna[1],
                               rna[0] - RO.EDGE * K_RO - 3.0, rna[1]))
        W = plot_x + plot_w
        geo = dict(plot_w=plot_w, plot_h=plot_h, box_w=box_w, h_cas=h_cas,
                   box_y=h_cas + GAP, int_lab_r=int_lab_r, left_w=left_w,
                   cas_x=cas_ink, box_x=0.0, ro_x=ro_x, ro_y=ro_y,
                   plot_x=plot_x, plot_y=plot_y)
        body += [g_cas, pop_box, g_pop, g_ro, plot] + step

    return geo, (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{(W + 2 * MARGIN) * 4:.0f}" '
            f'height="{(H + 2 * MARGIN) * 4:.0f}" '
            f'viewBox="{-MARGIN} {-MARGIN} {W + 2 * MARGIN:.2f} '
            f'{H + 2 * MARGIN:.2f}">'
            f'<style>{"".join(styles)}</style>'
            + "".join(body) + "</svg>"), W, H


# Panel b: the screen's output, the hit gate once per positive comparison.
# Two frames of the same size, STACKED - POS_A above, POS_B below - so they
# share an x axis by eye and the panel comes out about square. Each carries
# its own direct label (comparison and hit count), which is what tells them
# apart; no titles (STYLE_NOTES section 1).
GATES_GAP = 14.0                # between the two plot boxes
GATES_H = 100.0                 # each box, in panel units
GATES_W = 210.0                 # wide enough for the legend row, which is
                                # checked against the box, not guessed


def _split(comp):
    """The gate's points for ONE comparison; hits counted per frame."""
    return [p for p in _points() if p[2] == comp]


def build_gates():
    # 2026-11-17 (Leo): both frames are the hit gate; POS_A only above, POS_B
    # only below. Same limits, thresholds and thinning (gate_plot's own).
    body = [gate_plot(0.0, i * (GATES_H + GATES_GAP), GATES_W, GATES_H,
                      pts=_split(c), order=[c],
                      xlab="Splicing ratio (fraction spliced)")
            for i, c in enumerate(ORDER)]
    W, H = GATES_W, 2 * GATES_H + GATES_GAP
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{(W + 2 * MARGIN) * 4:.0f}" '
            f'height="{(H + 2 * MARGIN) * 4:.0f}" '
            f'viewBox="{-MARGIN} {-MARGIN} {W + 2 * MARGIN:.2f} '
            f'{H + 2 * MARGIN:.2f}">' + "".join(body) + "</svg>"), W, H


def _check():
    import xml.etree.ElementTree as ET
    for p in (CASSETTE, POPS) + ((READOUT,) if WITH_WORKFLOW else ()):
        assert p.exists(), f"missing source: {p} - run its own script first"
    geo, svg, W, H = build()
    ET.fromstring(svg)                       # it has to parse

    # the left column's two containers share a left edge
    assert abs(geo["cas_x"] - geo["box_x"]) < 0.01, \
        f"plasmid and populations box are {geo['cas_x']:.2f} apart"
    if WITH_WORKFLOW:
        # the workflow starts clear of the left column, ink to ink
        assert geo["ro_x"] >= geo["left_w"], "the workflow overlaps the left column"
        # the plot frame sits in the L's corner, clear of the labels on both sides
        assert geo["plot_x"] > geo["ro_x"] + (COL_LAB_R - RO_INK[0]) * K_RO, \
            "the plot frame runs into the workflow's column labels"
        assert geo["plot_y"] > geo["ro_y"] + (ROW_LAB_B - RO_INK[1]) * K_RO, \
            "the plot frame runs into the workflow's row labels"
    elif WITH_GATE:
        assert geo["plot_x"] >= geo["left_w"] + GAP - 0.01, "the plot overlaps the left column"
    if geo["plot_x"] is not None:
        assert geo["plot_y"] < H - 20, "the plot frame has almost no height"
    # the plot has to have room for its own axis labels and both gate lines
    assert geo["int_lab_r"] < geo["box_w"], \
        "the integration label runs past the populations box"
    # the two block gaps have to stay the SAME gap, which they only do while
    # the populations box is the widest thing in the left column
    assert not WITH_WORKFLOW or abs(geo["ro_x"] - geo["box_w"] - GAP) < 0.01, \
        (f"the box-to-workflow gap is {geo['ro_x'] - geo['box_w']:.1f}, not "
         f"{GAP} - something in the left column is now wider than the box")
    assert geo["box_y"] - geo["h_cas"] == GAP, "the cassette-to-box gap moved"
    assert GAP - 2 * LINE_H > 3.0, \
        "the integration label no longer clears the gap it sits in"
    # the gate plots, wherever they live now
    gsvg, GW, GH = (svg, geo["plot_w"], geo["plot_h"]) if WITH_GATE \
        else (build_gates()[0], GATES_W, GATES_H)
    assert GW - PAD[0] - PAD[1] > 60, "the plot is too narrow to read"
    assert GH - PAD[2] - PAD[3] > 40, "the plot is too short to read"
    # each frame has to yield its share of Fig 4c's 7 + 25. A frame draws its
    # hits in its comparison's fill and the legend swatch is the -1.
    if not WITH_GATE:
        ET.fromstring(gsvg)
        frames = [f for f in gsvg.split('<g data-name="Hit gate">')[1:]]
        got = {c: len(re.findall(rf'<circle[^>]*fill="{HIT_COL[c]}"', f)) - 1
               for c, f in zip(ORDER, frames)}
        print(f"  hits per frame: {got}")
        assert got == {"POS_A_vs_NEG_A": 7, "POS_B_vs_NEG_B": 25}, \
            f"the frames yield {got}, not Fig 4c's 7 + 25 = 32 ordered hits"

    raw = set(re.findall(r"font-size:([\d.]+)px", svg))
    raw |= set(re.findall(r'font-size="([\d.]+)', svg))
    # Every label in the panel resolves to 5.5. A source writes its own
    # number and the transform above it decides what that prints at, so the
    # three the panel can contain are named here; a fourth means a source has
    # changed its type and the scale that brings it in has to be rechecked.
    PRINTS_AT = {"5.5": LABEL,           # cassette and populations, as drawn
                 str(RO_LABEL): LABEL,   # the workflow, under K_RO
                 "18.0664": 4.5}         # the populations' "or", inside their
                                         # own 0.249 group - deliberately the
                                         # one smaller mark in the panel
    assert raw <= set(PRINTS_AT), f"an unknown label size: {raw - set(PRINTS_AT)}"
    got = sorted({PRINTS_AT[s] for s in raw})
    assert LABEL in got, f"no {LABEL} label: {got}"

    print(f"checks pass: {W:.0f} x {H:.0f} units, aspect {W / H:.2f}, "
          f"workflow {'at %.3f' % K_RO if WITH_WORKFLOW else 'OFF'}, labels {sorted(got)}")


def _emit(svg, name):
    p = OUT / name
    p.write_text(svg, encoding="utf-8")
    for ext in ("pdf", "png"):
        subprocess.run(["rsvg-convert", "-f", ext, "-o",
                        str(p.with_suffix("." + ext)), str(p)], check=False)
    print("wrote", p)


def main():
    _, svg, _, _ = build()
    _emit(svg, "Fig4a_screen_overview.svg")
    if not WITH_GATE:
        _emit(build_gates()[0], "Fig4b_screen_gates.svg")


if __name__ == "__main__":
    _check()
    if "--check" not in sys.argv:
        main()
