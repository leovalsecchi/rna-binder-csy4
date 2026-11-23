"""Target-presentation optimisation WORKFLOW: one wide strip, six steps.

The order and the named things are Methods 2.12.2-2.12.3, not a guess:

  1 module design      a sketch of the intended RNA geometry, converted to a
                       dot-bracket structure that preserves the Csy4 hairpin
  2 sequence design    PyFurnace 1.2.2, optimising the sequence against that
                       geometry with AlphaFold3 and ViennaRNA
  3 re-prediction      the fold is checked; a design that misses its geometry
                       goes back to step 2. Methods frames 2 and 3 as ONE
                       iterative loop ("iterating until each fragment satisfied
                       the synthesis constraints"), which is why the revise arrow
                       returns to the design step rather than to the sketch.
  4 cloning            Golden Gate with PaqCI into the ccdB entry vector
                       (pTWIST-Kan), transformed into E. coli Top10
  5 landing-pad bb     a second Golden Gate, BbsI, joining the insert to the
                       ccdB-weak cassette part pART177 into the pTWIST-Amp
                       landing-pad-compatible backbone, in ccdB-tolerant DB3.1
  6 integration        Bxb1 attB-GT into the landing pad of HEK293-H C8 Spartan,
                       integrase supplied in trans (pMLS596), Lipofectamine 3000
  7 readout            RNeasy Plus -> SuperScript IV -> KAPA HiFi (~850 bp) ->
                       ONT; splicing ratio = removed / (intact + removed)

2026-11-17 (Leo): "space out figure d more, more horizontal less vertical" and
"we also had more steps". The panel was 100 x 62 mm with three nodes; it is now
178 x 44 and carries all seven. It sits alone on its row, where the old width
left ~112 mm of the block empty.

## Two things the Methods changed, beyond the layout

**All four geometries are redrawn from Leo's own schematics**, at
`~/rna_binder_plots/thesis_figures/RNA_Scaff_Mod/schematics/` (`schematics.py`,
which asserts 8 hairpin rings per design). Coordinates here are that script's,
scaled into the panel; only the palette and the line weights are this figure's.

  mod1   linear backbone, four three-way junctions, TALL stems (24 u) and the
         hairpin pair rising 22 u above each junction
  mod2   a balanced binary tree: trunk, two branches, four junctions, eight
         hairpins - the branched design
  mod3   the same linear array as mod1 with SHORT stems (15 u, rise 12 u);
         that stem height is the only difference between them
  mod4   a closed ring, eight hairpins radiating outward, 5' trunk at the foot

Spacing along the array is EVEN. An earlier pass drew it uneven, reading
`F30_mod1_variablespac` and the Methods phrase "sharing no spacer sequence" as
variable spacing; the schematics settle it - the spacer SEQUENCES differ, the
spacing does not, which is also what "near-uniformly spaced" says.

**All four are drawn, though only three were ever built.** Of ART243-ART254
(mod1 = 243-245, mod2 = 246-248, mod3 = 249-251, mod4 = 252-254), mod3 never
existed: ART249 was never delivered and pART250/251 failed to assemble. It is
still drawn, because the panel's first step is what was DESIGNED, and four
architectures were designed; what survived assembly belongs in the caption.
Bare panel, no title; see README.md.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch
from common import style, paths, palette, scaffold_glyphs

style.apply()
plt.rcParams["hatch.linewidth"] = 0.6

# Construct variants, not design stacks -> theme steps + neutrals (rules 1, 3).
INK, MUTED, LINEG, TEAL = (palette.INK, palette.NEUTRAL["600"],
                           palette.NEUTRAL["400"], palette.THEME["700"])
# the one warm mark in the panel: the diverging ramp's warm pole (rule 4 reserves
# the RAMP for signed data; a single pole as an annotation colour is not a ramp)
LOOP = palette.DIVERGING[-1]
# Display names. The CSV, the keys and the README provenance all keep mod1 /
# mod4; only the printed labels change, so nothing downstream has to be re-mapped.
DISPLAY = {"mod1": "F30 variation", "mod4": "Circular scaffold"}
disp = lambda k: DISPLAY.get(k, k)
W, H = 178.0, 38.0                       # mm, full block; was 100 x 62 with three nodes
fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
pa = fig.add_axes([0, 0, 1, 1]); pa.set_xlim(0, W); pa.set_ylim(H, 0); pa.axis("off"); pa.set_zorder(-1)

def seg(p, q, col, lw=1.3):
    pa.plot([p[0], q[0]], [p[1], q[1]], color=col, lw=lw, solid_capstyle="round", zorder=2)

def dot(x, y, r=0.8):                           # Csy4 hairpin: one filled dot
    pa.add_patch(Circle((x, y), r, color=LOOP, lw=0, zorder=3))

def arrow(p, q, col=MUTED, lw=0.8):
    pa.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="-|>", lw=lw, color=col,
                mutation_scale=5.5, shrinkA=0, shrinkB=0))

def twisted_hairpin(cx, ybase, sc, col, cross_at=0.72, loop_h=10.0, ybot=43.0, ytop=22.0,
                    amp=7.2, taper=0.32, op=0.45, lw=1.1):
    """Same RNA stem-loop glyph as the pipeline figure (common/glyphs.py), drawn in mm."""
    pw = math.log(0.5) / math.log(cross_at)
    def rail(sg, u):
        a = amp * (1.0 - taper * u ** 3)
        return (cx + sc * sg * a * math.cos(math.pi * u ** pw), ybase - sc * (ybot - ytop) * u)
    def run(sg, u0, u1, n=22):
        pts = [rail(sg, u0 + (u1 - u0) * i / n) for i in range(n + 1)]
        pa.plot([p[0] for p in pts], [p[1] for p in pts], color=col, lw=lw,
                solid_capstyle="round", solid_joinstyle="round", zorder=2)
    run(1, 0.0, cross_at + 0.005); run(-1, 0.0, 1.0); run(1, cross_at - 0.005, 1.0)
    for u in (0.10, 0.30, 0.50, 0.66):                   # base pairs in the stem
        (x1, y1), (x2, y2) = rail(1, u), rail(-1, u)
        pa.plot([x1, x2], [y1, y2], color=col, lw=0.6, alpha=op, solid_capstyle="round", zorder=2)
    (lx, ly), (rx, ry) = rail(1, 1.0), rail(-1, 1.0)
    h = loop_h * sc                                # cubic bezier closing the loop
    t = [i / 24 for i in range(25)]
    pa.plot([(1-u)**3*lx + 3*(1-u)**2*u*lx + 3*(1-u)*u**2*rx + u**3*rx for u in t],
            [(1-u)**3*ly + 3*(1-u)**2*u*(ly-h) + 3*(1-u)*u**2*(ry-h) + u**3*ry for u in t],
            color=col, lw=lw, solid_capstyle="round", zorder=2)
    return (lx + rx) / 2, ly                       # foot of the stem


def mini_ring(cx, cy, col, bad=False):
    """The circular scaffold at node size. `bad` is the mis-predicted variant:
    the ring opens, the arms splay unevenly and one hairpin is left dangling."""
    r = 1.8
    if bad:                                   # open arc instead of a closed ring
        th = [math.radians(a) for a in range(35, 330, 5)]
        pa.plot([cx + r * math.cos(t) for t in th], [cy - r * math.sin(t) for t in th],
                color=col, lw=1.1, solid_capstyle="round", zorder=2)
    else:
        pa.add_patch(Circle((cx, cy), r, fill=False, ec=col, lw=1.1, zorder=2))
    angs = ([90, 140, 190, 240, 290, 340, 30] if bad          # uneven, one missing
            else [90 + 360 * k / 7 for k in range(7)])
    for j, a in enumerate(angs):
        rad = math.radians(a)
        ca, sa = math.cos(rad), -math.sin(rad)
        arm = 1.3 + (0.8 if (bad and j in (1, 4)) else 0.0)    # two arms over-extended
        pa.plot([cx + r * ca, cx + (r + arm) * ca], [cy + r * sa, cy + (r + arm) * sa],
                color=col, lw=1.1, solid_capstyle="round", zorder=2)
        pa.add_patch(Circle((cx + (r + arm + 0.7) * ca, cy + (r + arm + 0.7) * sa),
                            0.45, color=LOOP, lw=0, zorder=3))


def badge(cx, cy, ok):
    col = TEAL if ok else MUTED
    pa.add_patch(Circle((cx, cy), 1.7, color=col, lw=0, zorder=3))
    if ok:
        pa.plot([cx - 0.85, cx - 0.2, cx + 0.95], [cy, cy + 0.68, cy - 0.8],
                color="white", lw=0.8, solid_capstyle="round",
                solid_joinstyle="round", zorder=4)
    else:
        for sg in (1, -1):
            pa.plot([cx - 0.72, cx + 0.72], [cy - 0.72 * sg, cy + 0.72 * sg],
                    color="white", lw=0.8, solid_capstyle="round", zorder=4)


def blob(cx, cy, r, col, a1=0.07, a2=0.05, p1=0.0, p2=1.1, k1=3, k2=5, z=2, lw=1.2):
    t = [2 * math.pi * i / 72 for i in range(73)]
    rr = [r * (1 + a1 * math.sin(k1 * u + p1) + a2 * math.sin(k2 * u + p2)) for u in t]
    pa.fill([cx + q * math.cos(u) for q, u in zip(rr, t)],
            [cy + q * math.sin(u) for q, u in zip(rr, t)],
            facecolor="white", edgecolor=col, lw=lw, joinstyle="round", zorder=z)


def glyph_cell(cx, cy, col, k=1.0):
    blob(cx - 4.9 * k, cy - 4.3 * k, 2.0 * k, col, p1=0.9, p2=2.4, z=2, lw=1.0)   # two smaller cells behind
    pa.add_patch(Circle((cx - 4.9 * k, cy - 4.3 * k), 0.62 * k, color=col, lw=0, zorder=2.2))
    blob(cx - 3.8 * k, cy + 3.7 * k, 2.2 * k, col, p1=2.1, p2=0.4, z=2, lw=1.0)
    pa.add_patch(Circle((cx - 3.8 * k, cy + 3.7 * k), 0.78 * k, color=col, lw=0, zorder=2.2))
    blob(cx + 0.9 * k, cy - 0.3 * k, 4.3 * k, col, p1=0.3, p2=1.7, z=3)           # the cell in front
    pa.add_patch(Circle((cx + 0.9 * k, cy - 0.3 * k), 1.75 * k, color=col, lw=0, zorder=3.2))

# Construct variants, not design stacks -> theme steps + neutrals (rules 1, 3).
PLUM, GREY = palette.THEME["400"], palette.NEUTRAL["600"]

# ---------------------------------------------------------------- strip geometry
X0, X1 = 3.0, 175.0
N = 6            # was 7; the two cloning steps are one glyph (2026-11-19)
# Step 1 carries four architectures side by side, so it gets its own width and
# the other six share what is left; equal slots would squeeze them to 6 mm each.
GEO_W = 36.0
SLOT = (X1 - X0 - GEO_W) / (N - 1)
CX = [X0 + GEO_W / 2] + [X0 + GEO_W + SLOT / 2 + i * SLOT for i in range(N - 1)]
ICON_CY = 14.0            # glyph centre (step 1 stacks two rows)
LAB_Y = 25.5              # first label line, drawn downward
LINE_H = 3.3
LAB_PT = 6.0
REVISE_Y = 3.4            # the loop runs OVER the strip (Leo, 2026-11-17)


def chevron(x, y, col=LINEG, s=1.9, lw=0.9):
    """Step-to-step mark. Fig1a's idiom: a chevron, not an arrow per step - the
    strip reads as one flow and only the revise loop needs a head."""
    pa.plot([x - s * 0.55, x + s * 0.35, x - s * 0.55],
            [y - s, y, y + s], color=col, lw=lw,
            solid_capstyle="round", solid_joinstyle="round", zorder=2)


def label(cx, lines):
    for j, ln in enumerate(lines):
        pa.text(cx, LAB_Y + j * LINE_H, ln, ha="center", va="top",
                fontsize=LAB_PT, color=INK)


# Panel background (Leo, 2026-11-17). NEUTRAL 100 is the house panel tint - the
# same fill Fig 1a's first band and Fig 1d's inner boxes carry. palette.SURFACE
# is white, so it is not the one to reach for here.
BG_FILL = palette.NEUTRAL["100"]
# The scaffold strand tone, shared with Fig 1c and fig04 panel a
# (`common/scaffold_glyphs.py`), so a target scaffold is the same colour wherever
# it appears (Leo, 2026-11-17). All four variants take it - they are one kind of
# object, told apart by shape and by their labels, not by hue.
# CONTRAST: 2.24:1 against this panel's background and 2.44:1 even on white, so
# it is under the 3:1 COLORS.md asks for. The strands are drawn a little heavier
# here to hold up; if it still reads weak in print, either drop the background
# behind step 1 or move the strands to THEME 700 (6.67:1).
SCAF = palette.SYNBIO["rna"]
BG_X0, BG_X1 = X0 - 2.0, X1 + 2.0
BG_Y0, BG_Y1 = 0.6, LAB_Y + 2 * LINE_H + 2.2
pa.add_patch(FancyBboxPatch((BG_X0, BG_Y0), BG_X1 - BG_X0, BG_Y1 - BG_Y0,
                            boxstyle="round,pad=0,rounding_size=2.5",
                            fc=BG_FILL, ec="none", zorder=-2))


# ---------------------------------------------------------------- the seven glyphs
# The four architectures are drawn by `common/scaffold_glyphs.py` - the SAME
# module that makes the target-scaffold glyphs in Fig 4a, via
# `Fig1f_scaffold_splicing/make_glyphs.py`, whose SVGs Fig 4a imports from
# `output/Fig1f_scaffold_splicing/elements/`. Drawing them here from the shared
# glyphs rather than from a local copy is the point (Leo, 2026-11-17): a target
# scaffold is the same object in both figures, so it is built by one function.
#
# mod1 and mod3 are the same `array` with a long and a short stem, and they
# share ARRAY_K so that arm length - the only thing separating them - survives.
# mod2 is `tree`, added to the module for this. mod4 is `ring`.
# Numeric order, reading across then down: mod1 mod2 / mod3 mod4.
GEOMS = ("mod1", "mod2", "mod3", "mod4")
GEOM_LAB = {k: k for k in GEOMS}

ARRAY_K, TREE_K, RING_K = 0.62, 0.43, 0.46
# One red dot size for all four: the Csy4 hairpin tip is the SAME feature in
# every architecture, so it must not shrink because mod2 is drawn smaller.
DOT_R = 0.8 * ARRAY_K


def mod2_tree():
    """mod2's REAL branching, read from its dot-bracket target.

    Returned as (helix_bp, [children]); a leaf is a Csy4 hairpin. Parsed rather
    than transcribed so the glyph cannot drift from the design: mod2 is a
    CASCADE, not a balanced tree - trunk 8 bp, then 7/6, 7/6, 13/13, with 16 bp
    hairpin arms - and that asymmetry is the shape."""
    import re
    txt = (paths.INPUTS / "csy4_dotbracket" / "F30_mod2_pos.dbn").read_text()
    ss = [l for l in txt.split("\n") if l and set(l) <= set("().")][0]
    unit = ss[: len(ss) // 3]                       # one repeat = the 8-site form
    st, pt = [], [-1] * len(unit)
    for i, c in enumerate(unit):
        if c == "(":
            st.append(i)
        elif c == ")":
            j = st.pop(); pt[i] = j; pt[j] = i

    def kids(i, j):
        out, m = [], i + 1
        while m < j:
            if pt[m] > m:
                out.append((m, pt[m])); m = pt[m] + 1
            else:
                m += 1
        return out

    def node(i, j):
        n, a, b = 0, i, j
        while a < b and pt[a] == b:                 # count consecutive bp
            n += 1; a += 1; b -= 1
            while a < b and unit[a] == ".":
                a += 1
            while a < b and unit[b] == ".":
                b -= 1
        ch = [c for c in kids(a - 1, b + 1) if c != (i, j)] if a <= b else []
        return (n, [node(u, v) for u, v in ch])

    top = kids(-1, len(unit))
    assert len(top) == 1, f"mod2 should have one top-level helix, got {len(top)}"
    return node(*top[0])


MOD2 = mod2_tree()
# Every glyph's STRAND lands on one baseline per row. The shared glyphs each
# measure their strand from `y1` by a different offset (4.5k for array and
# tree, 3.5k for ring), so y1 is solved per kind rather than guessed - which
# is what put mod4's strand below the step labels on the first pass.
STRAND_DY = 2.8
OVERHANG = 1.1      # 5' and 3' strand either side of the scaffold's foot
STEM_LONG, STEM_SHORT = 3.6, 1.8
# The four architectures are four different heights, so a label hung at a fixed
# offset from the CELL sat 3.2 mm clear of mod1 and 0.7 mm clear of mod2. Each
# label is placed off the ink actually drawn instead, all at this one gap.
LAB_GAP = 0.9
GLYPH_LW = 0.9


def glyph_geoms(cx, cy):
    """Step 1: the four designed architectures, 2 x 2, from the shared glyphs."""
    CW, CH = GEO_W / 2, 11.5
    for i, kind in enumerate(GEOMS):
        row, col = divmod(i, 2)
        gx = cx - GEO_W / 2 + CW / 2 + col * CW
        gy = cy - CH / 2 + row * CH
        # The 5' and 3' overhangs are SHORT (Leo, 2026-11-17). The strand is
        # sized from each scaffold's own FOOT plus OVERHANG, not set to one
        # width for all four: the array stands on four junctions and needs to
        # span them, while the tree and the ring meet the transcript at a single
        # point and need almost nothing either side.
        if kind in ("mod1", "mod3"):
            half = 2.5 * ARRAY_K * 3 + OVERHANG    # spans the four junctions
        else:
            half = 2.6 + OVERHANG                  # one foot: a short strand
        x0, x1 = gx - half, gx + half
        base = gy + STRAND_DY                      # where the transcript sits
        n_l, n_p = len(pa.lines), len(pa.patches)
        if kind in ("mod1", "mod3"):
            scaffold_glyphs.array(
                pa, gx, base + 4.5 * ARRAY_K, SCAF, LOOP, x0, x1, n=4,
                stem=STEM_LONG if kind == "mod1" else STEM_SHORT,
                lw=GLYPH_LW, k=ARRAY_K, dot_r=DOT_R)
        elif kind == "mod2":
            scaffold_glyphs.tree(pa, gx, base + 4.5 * TREE_K, SCAF, LOOP,
                                 x0, x1, MOD2, lw=GLYPH_LW, k=TREE_K,
                                 dot_r=DOT_R)
        else:
            scaffold_glyphs.ring(pa, gx, base + 3.5 * RING_K, SCAF, LOOP,
                                 x0, x1, n=8, r=3.2, lw=GLYPH_LW, k=RING_K,
                                 dot_r=DOT_R)
        top = min([min(l.get_ydata()) for l in list(pa.lines)[n_l:]]
                  + [c.center[1] - c.radius for c in list(pa.patches)[n_p:]])
        pa.text(gx, top - LAB_GAP, GEOM_LAB[kind], ha="center", va="bottom",
                fontsize=4.6, color=MUTED)


def glyph_seq(cx, cy):
    """Step 2: sequence design - bases laid under a held structure."""
    seg((cx - 6.2, cy - 4.4), (cx + 6.2, cy - 4.4), MUTED, 0.9)
    for k, (x0, x1) in enumerate(((-6.0, 1.2), (-3.0, 5.8), (-6.2, 3.0))):
        seg((cx + x0, cy - 1.2 + k * 2.6), (cx + x1, cy - 1.2 + k * 2.6), TEAL, 1.0)


# The FAIL badge is on top and the pass below (Leo, 2026-11-17), so the revise
# loop leaves from the badge it belongs to and runs over the strip rather than
# crossing it.
BADGE_DX, BADGE_FAIL_DY, BADGE_PASS_DY = 4.6, -2.6, 3.0


def glyph_predict(cx, cy):
    """Step 3: the fold is re-predicted, accepted or sent back.

    The structure shown is the SAME scaffold glyph step 1 draws, not a bespoke
    ring (Leo, 2026-11-17): it is the same molecule being checked."""
    k = 0.50
    scaffold_glyphs.ring(pa, cx - 2.6, cy + 4.0 + 3.5 * k, SCAF, LOOP,
                         cx - 8.2, cx + 3.0, n=8, r=3.2, lw=GLYPH_LW, k=k)
    badge(cx + BADGE_DX, cy + BADGE_FAIL_DY, False)
    badge(cx + BADGE_DX, cy + BADGE_PASS_DY, True)


# Panel e draws its donor plasmid as a ROUNDED RECTANGLE with the cassette on
# the bottom edge, and panel b's three component plasmids are the same shape.
# These two match it (Leo, 2026-11-17) instead of being plain circles, so a
# plasmid looks like a plasmid everywhere in the figure.
PL_W, PL_H, PL_R = 13.0, 9.2, 2.8


# The insert sits ON the bottom edge and attB hangs below it, so a box centred
# on cy leaves the glyph's INK centre ~0.5 mm low against every other step.
PL_LIFT = 0.5


def _plasmid_box(cx, cy, col=TEAL, lw=0.9):
    cy -= PL_LIFT
    pa.add_patch(FancyBboxPatch((cx - PL_W / 2, cy - PL_H / 2), PL_W, PL_H,
                                boxstyle=f"round,pad=0,rounding_size={PL_R}",
                                fc="none", ec=col, lw=lw, zorder=2))
    return cy + PL_H / 2                      # the bottom edge: y grows downward


def glyph_clone(cx, cy):
    """Step 4: cloning, as ONE glyph.

    2026-11-19 (Leo): "one cloning glyph is enough (make them into one)". This
    was two steps - the insert closed into the ccdB entry vector (PaqCI), then
    the same insert moved into the landing-pad-compatible backbone (BbsI), each
    with its own plasmid ring. Two near-identical rings side by side read as a
    repeat rather than as two assemblies, and the strip has one slot per idea.

    The merged glyph carries BOTH facts the two carried: the insert's cut sites
    (what cloning does) and the flanking attB squares (what the landing-pad
    backbone adds, and what step 6 then integrates). The two Golden Gates stay
    in the Methods and in the caption; only the picture merges.
    """
    y = _plasmid_box(cx, cy)
    seg((cx - 2.6, y), (cx + 2.6, y), LOOP, 1.8)           # the insert
    for sx in (-2.6, 2.6):                                 # its cut sites
        seg((cx + sx, y - 1.0), (cx + sx, y + 1.0), LINEG, 0.7)
    for sx in (-4.6, 4.6):                                 # attB, flanking it
        pa.add_patch(FancyBboxPatch((cx + sx - 1.15, y - 0.9), 2.3, 1.8,
                                    boxstyle="round,pad=0,rounding_size=0.4",
                                    fc="white", ec=GREY, lw=0.8, zorder=4))


def glyph_reads(cx, cy):
    """Step 7: long reads spanning the construct, which is what the ratio counts."""
    # TEAL, not INK: the construct is a molecule like every other one here,
    # and black was the only pure-black mark on the strip.
    seg((cx - 7.0, cy - 4.6), (cx + 7.0, cy - 4.6), TEAL, 1.0)
    for k, (x0, x1) in enumerate(((-6.4, 1.6), (-2.6, 6.6), (-5.2, 3.2))):
        seg((cx + x0, cy - 1.0 + k * 2.7), (cx + x1, cy - 1.0 + k * 2.7), MUTED, 1.0)


# ---------------------------------------------------------------- draw
STEPS = [
    (glyph_geoms,   ["Array design"]),
    (glyph_seq,     ["Sequence design", "(PyFurnace)"]),
    (glyph_predict, ["Structure", "re-prediction"]),
    (glyph_clone,   ["Cloning into", "landing-pad backbone"]),
    # 0.86: the cluster is the widest-spread glyph here and read as the
    # heaviest step on a strip where every step is one step.
    (lambda cx, cy: glyph_cell(cx, cy, TEAL, 0.86), ["Integration", "into cells"]),
    (glyph_reads,   ["Splicing readout", "(ONT)"]),
]
assert len(STEPS) == N, "N must match STEPS"

# Half-width of each step's drawn content, for placing the chevrons. Step 1 is
# the wide one; the rest are single glyphs.
HALF = [GEO_W / 2] + [8.6] * (N - 1)

for i, (draw, lines) in enumerate(STEPS):
    draw(CX[i], ICON_CY)
    label(CX[i], lines)
    if i < N - 1:
        chevron((CX[i] + HALF[i] + CX[i + 1] - HALF[i + 1]) / 2, ICON_CY)

# The revise loop runs ABOVE the strip and leaves the FAIL badge: a fold that
# misses its geometry goes back to the DESIGN step, not to the sketch - Methods
# frames design and re-prediction as one iteration.
BX = CX[2] + BADGE_DX                       # the fail badge's own x
seg((BX, ICON_CY + BADGE_FAIL_DY - 1.7), (BX, REVISE_Y), TEAL, 0.8)
seg((BX, REVISE_Y), (CX[1], REVISE_Y), TEAL, 0.8)
arrow((CX[1], REVISE_Y), (CX[1], ICON_CY - 6.6), TEAL)
pa.text((CX[1] + BX) / 2, REVISE_Y - 0.8, "revise", ha="center", va="bottom",
        fontsize=5.0, color=MUTED)

out = paths.fig_output("Fig1f_scaffold_splicing")
fig.savefig(out / "scaffold_workflow.svg")
style.save(fig, out, "scaffold_workflow")
print("wrote", out / "scaffold_workflow.svg", f"{W}x{H} mm, {N} steps")
