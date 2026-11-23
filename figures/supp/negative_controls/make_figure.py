"""Negative controls: what is replaced, what that does to the fold at two
scales, which scaffolds were never obtained, and the colony screen that
produced the clone actually used.

SIX bare panels, from the tables `make_inputs.py` writes. Rearranged and split
2026-11-19 on Leo's note, "rearrange the pannels (maybe also split b into two
panels) / also would be good to see two failed scaffolds and how much shorter
they were", and corrected 2026-11-20 on his clarification of which failure he
meant: "with the failed scaffolds, I meant the negatives that failed and had
missing repeats because of the spike in fail."

Order, and why:

  a  the swap at one hairpin position        sequence scale
  b  the two elements, ViennaRNA MFE fold    element scale
  c  the whole insert, folding energy        construct scale
  d  the negatives that lost whole repeats   construct scale, as obtained
  e  the colony screen                       construction, clone level
  f  the two scaffolds never obtained        construction, design level

Every panel emits its AUTHORED 85 mm canvas - `style.save(..., tight=False)` -
and that is a layout requirement, not a detail. `make_supp_book.py` packs panels
greedily into rows no wider than the 182 mm text block, so at the 90-92 mm that
`tight` padded these panels out to, 85 + 6 + 85 = 176 fitted but 90 + 6 + 90 =
186 did not, and all six landed ONE PER ROW, 299 mm of content on a 283 mm page
body. Leo, 2026-11-20: "Align the panels. Should be two per horizontal level
always." They now emit 85.0 mm each and the book pairs them. `check_emitted()`
re-derives the packing with make_supp_book's OWN `rows_for` at the end of every
build, so the figure cannot quietly go back to one per row.

The cost of `tight=False` is that ink outside the canvas is CLIPPED rather than
padded around, so four places that used to overflow were brought inside: (a)'s
left labels (2.1 mm off the left edge, now placed off their measured width),
(c)'s axis caption (0.2 mm, canvas 25 -> 26 mm), (f)'s caption row (2.5 mm,
42 -> 45 mm) and (e)'s right spine (1.8 mm, right=0.99 -> 0.96).

a-b-c is one ascending scale, which is the distinction this figure exists to
keep straight. The claim is made at (b) - *the replacement cannot form THE
hairpin* - and bounded immediately at (c), where the negative insert still
folds 101 pairs over 319 nt, so no reader carries "the negative is
unstructured" past the panel that disproves it. The old panel (b) carried both
scales at once, which is exactly the conflation the split removes.

(d) then states what was actually obtained - pART257 at SEVEN F30 units against
pART158's eight, and, since 2026-11-20, pLV16 at FOUR against pLV13's eight -
and (e) is the screen that explains the first of those: 0/56 (0%) colonies
carried eight units. The pLV16 pair is there on Leo's note "for d add the
negative with less copies in the 1 by 1 resulting in pLV15"; it is pLV16, not
pLV15, confirmed with him, and the README records what pLV15 is. (f) is a DIFFERENT failure mode, two scaffolds whose
assembly was never obtained at all, and it is last so it cannot be read as part
of the spike-in story. It was panel (d) until 2026-11-20; it is a real measured
result (mod3's repeat unit is the longest of the four, so pART250 and pART251
are 74 to 134 nt LONGER than the siblings that assembled at the same copy
number) and is kept as its own panel rather than dropped to the caption.

The caution this figure has to carry: the replacement removes secondary
structure as well as sequence at every hairpin position, so a binder that loses
signal on the negative has not been shown to be reading sequence rather than
shape.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name

INK = palette.INK
MUTED = palette.NEUTRAL["600"]
LINEG = palette.NEUTRAL["400"]
MAIN = palette.MAIN
LOOP = palette.DIVERGING[-1]            # the Csy4 site, COLORS rule 4
NEUT = palette.NEUTRAL["300"]           # the negative: a reference, rule 3
BGFILL = palette.NEUTRAL["100"]
# Palette v0.9 (2026-11-20), which SUPERSEDES COLORS rule 9: the four designed
# scaffold architectures are ONE colour, and since v0.9 that colour is theme
# step 600 `#0678a6` rather than `palette.MAIN` `#015d82` (v0.8) - lighter, on
# Leo's instruction. The indigo ramp was withdrawn at v0.8 for being too dark to
# separate (adjacent dE 9.8-10.9).
#
# Nothing in this file names the colour: it is read from
# `palette.SCAFFOLDS[arch]`, so a palette revision lands without a code change.
# The consequence of one colour for four architectures is unchanged and is the
# reason this figure reads as it does: COLOUR NO LONGER NAMES AN ARCHITECTURE -
# text does. Every mark in (a), (c), (d) and (f) carries its architecture in its
# own label.
#
# F30 is NOT one of the four and keeps `palette.F30`, the Fig 1 construct-map
# blue, so the F30 array still reads as a different kind of object from the mod
# rings. Csy4 sites stay in the diverging warm end (rule 4).
MOD4 = palette.SCAFFOLDS["mod4"]
F30 = palette.F30
FS = 5.5

# The flank either side of one hairpin position, from the construct files. Both
# constructs carry it unchanged, which is what makes the swap a swap.
LEFT = "GGTTAGCTG"
RIGHT = "CAGCACTGCG"


def load():
    d = pd.read_csv(paths.INPUTS / "supp_negative_controls_design.csv")
    f = pd.read_csv(paths.INPUTS / "supp_negative_controls_fold.csv").set_index("insert")
    c = pd.read_csv(paths.INPUTS / "supp_negative_controls_colonies.csv",
                    dtype={"plate": str})
    xy = pd.read_csv(paths.INPUTS / "supp_negative_controls_foldxy.csv")
    sc = pd.read_csv(paths.INPUTS / "supp_negative_controls_scaffolds.csv")
    return d.set_index("negative"), f, c, xy, sc


def family(design):
    """{construct: fill}, from the table rather than from a hardcoded dict.

    Needed from 2026-11-20, when panel (d) gained a third pair: the constructs a
    panel draws are whatever the table holds, so the fill a construct gets has
    to be looked up by its architecture and its role, not keyed by name in this
    script. A positive takes its family colour (F30 is not one of the four
    designed architectures and keeps the Fig 1 construct-map blue); a negative
    is a reference, not a series, so it takes NEUTRAL 300 - COLORS rule 3.
    """
    out = {}
    for neg, r in design.iterrows():
        out[r.positive] = F30 if r.arch == "F30" else MOD4
        out[neg] = NEUT
    return out


def arch_of(design):
    """{construct: architecture}, from the table rather than from this script.

    Needed from palette v0.8 (2026-11-20): one colour now serves all four
    scaffold architectures, so a panel that draws them has to NAME them, and
    every panel here that draws a construct names its architecture beside it.
    """
    return {p: r.arch for _, r in design.iterrows() for p in (r.positive, r.name)}


def w_mm(s, size=FS):
    """Ink width of a label in mm, MEASURED - see common/style.text_w."""
    return style.text_w(s, size) / 72.0 * 25.4


# --- reading a fold --------------------------------------------------------
# scipy is unusable in this venv (`_propack` dlopen), so the pairing and the
# geometry below are hand-rolled; none of it needs more than numpy.
def pairs_of(db):
    """0-based (i, j) base pairs from a dot-bracket string."""
    st, out = [], []
    for k, c in enumerate(db):
        if c == "(":
            st.append(k)
        elif c == ")":
            out.append((st.pop(), k))
    return sorted(out)


# How many hydrogen bonds hold a pair, which is the whole difference between
# these two stems. A G.U wobble is two bonds but NOT a Watson-Crick pair, so it
# is kept as its own class rather than lumped in with A-U.
def pair_class(a, b):
    s = frozenset((a, b))
    return {frozenset("GC"): "GC", frozenset("AU"): "AU",
            frozenset("GU"): "GU"}.get(s, "other")


# Rung weight IS the pair chemistry - STYLE_NOTES section 4, class by texture
# and not by colour, and COLORS rule 8, base-pairing ticks stay neutral.
RUNG = {"GC": dict(lw=1.2, ls="-"),      # three bonds
        "AU": dict(lw=0.5, ls="-"),      # two
        "GU": dict(lw=0.5, ls=(0, (0.8, 1.0)))}   # two, but a wobble


def oriented(p, db):
    """Rotate a naview layout stem-horizontal, hairpin loop to the RIGHT.

    naview hands back whatever orientation its own layout fell into, which for
    these two differs by ~40 degrees. Two folds compared side by side have to
    be in the same pose or the reader is comparing poses, so both are rotated
    onto the axis that runs from the outermost pair to the loop.
    """
    pr = pairs_of(db)
    i, j = pr[0]                                    # outermost pair: stem foot
    foot = (p[i] + p[j]) / 2
    loop = p[pr[-1][0] + 1:pr[-1][1]].mean(0)       # the bases inside the stem
    a = np.arctan2(*(loop - foot)[::-1])            # angle of foot -> loop
    c, s = np.cos(-a), np.sin(-a)
    q = (p - foot) @ np.array([[c, s], [-s, c]])
    return q * [1, -1] if q[0, 1] < 0 else q        # 5' tail always on top


def draw_fold(ax, sub, x0, yaxis, scale, col):
    """One MFE fold, drawn from its ViennaRNA layout. Returns (width, half-height).

    `x0` is the fold's left edge and `yaxis` the height of its STEM AXIS, which
    `oriented` has put at y = 0, so both folds sit on one line however much
    their unpaired tails fan out around it.

    `scale` is mm per naview unit and is passed in rather than fitted per fold,
    so the 28 nt repeat and the 19 nt replacement are drawn at ONE scale and
    the reader sees they are not the same length.
    """
    seq = "".join(sub.base)
    db = "".join(sub.dot_bracket)
    q = oriented(sub[["x", "y"]].to_numpy(float), db) * scale
    q = q - [q[:, 0].min(), 0] + [x0, yaxis]
    ax.plot(q[:, 0], q[:, 1], color=col, lw=0.5, zorder=2,
            solid_capstyle="round")
    for i, j in pairs_of(db):
        v = q[j] - q[i]
        # INK, not MUTED: the negative's BACKBONE is MUTED, so a muted rung
        # and a muted backbone segment were the same mark, and the loop arc
        # read as if it carried a fifth pair.
        ax.plot(*zip(q[i] + v * 0.25, q[j] - v * 0.25), color=INK, zorder=3,
                **RUNG[pair_class(seq[i], seq[j])])
    for (x, y), b in zip(q, seq):
        ax.text(x, y, b, ha="center", va="center", fontsize=FS,
                color=col, zorder=4,
                bbox=dict(fc="white", ec="none", pad=0.3))
    # Step OUTWARD along the chain from base 1, not a fixed offset: on the
    # replacement a fixed one landed the marker on top of base 1 itself.
    d = q[0] - q[1]
    ax.text(*(q[0] + d / np.hypot(*d) * 1.7), "5'", ha="center", va="center",
            fontsize=FS - 1, color=MUTED, zorder=4)
    return np.ptp(q[:, 0]), np.ptp(q[:, 1]) / 2


def _selfcheck(design, fold, xy, sc):
    """The drawn fold must BE the fold the table reports, and the stems must be
    what the sequences allow. Cheap, and it is the one place a silent regression
    in the dot-bracket parse or the naview orientation would otherwise show up
    only as a wrong picture that still looks like a hairpin.
    """
    want = {"Csy4 repeat alone": dict(GC=4, AU=1, GU=0),
            "replacement alone": dict(GC=1, AU=2, GU=1)}
    for name, exp in want.items():
        sub = xy[xy.element == name].sort_values("i")
        seq, db = "".join(sub.base), "".join(sub.dot_bracket)
        pr = pairs_of(db)
        assert len(sub) == fold.loc[name, "nt"], f"{name}: layout/length mismatch"
        assert len(pr) == fold.loc[name, "predicted_pairs"], \
            f"{name}: drawing {len(pr)} pairs, table says " \
            f"{fold.loc[name, 'predicted_pairs']}"
        got = {k: sum(pair_class(seq[i], seq[j]) == k for i, j in pr)
               for k in ("GC", "AU", "GU")}
        assert got == exp, f"{name}: stem chemistry {got}, expected {exp}"
        q = oriented(sub[["x", "y"]].to_numpy(float), db)
        i, j = pr[0]
        assert abs(q[i, 1] + q[j, 1]) < 1e-6 and q[pr[-1][0], 0] > q[i, 0], \
            f"{name}: stem is not horizontal with the loop to the right"
    # The claim the panel exists to protect, checked rather than trusted: the
    # replacement's best WATSON-CRICK stem is 3 bp and all A-U, which is not
    # the same statement as the MFE fold's 4 pairs (that one uses a G.U
    # wobble). Both numbers are right under their own rule; the README says so.
    wc = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
    s = "".join(xy[xy.element == "replacement alone"].sort_values("i").base)
    best = (0, "")
    for i in range(len(s)):
        for j in range(i + 4, len(s)):       # >= 3 nt between the strands
            n = 0
            while i + n < j - n and (s[i + n], s[j - n]) in wc:
                n += 1
            best = max(best, (n, s[i:i + n]))
    assert best == (3, "AUU"), f"best Watson-Crick stem is {best}, expected (3, 'AUU')"
    # The negative scaffold is NOT unstructured, which is the other half of the
    # claim and the reason (c) exists as its own panel.
    assert fold.loc["pART256 insert", "predicted_pairs"] > 100, \
        "pART256's insert no longer folds; panel (c)'s point has gone"
    _check_missing_repeat(design)
    _check_failed(sc)


def _check_missing_repeat(design):
    """Panel (d)'s arithmetic, re-closed from the table it draws.

    `make_inputs.design()` asserts this against the construct files; this
    re-asserts it against the CSV, so a hand-edited or stale table cannot put a
    number on the canvas that its own bars contradict. Exactly ONE of the two
    pairs lost a whole repeat, and each shortfall is one lost unit plus 9 nt at
    every hairpin position the negative kept.
    """
    lost = {n: int(r.positive_units - r.negative_units)
            for n, r in design.iterrows() if r.positive_units > r.negative_units}
    assert lost == {"pART257": 1, "pLV16": 4}, \
        f"the pairs short a repeat are {lost}, not pART257 by 1 and pLV16 by 4"
    for n, r in design.iterrows():
        assert r.delta_nt == (r.positive_units - r.negative_units) \
            * r.positive_unit_nt + r.nt_per_position * r.loops_in_negative, \
            f"{n}: the shortfall does not close"
        assert r.delta_nt == r.positive_array_nt - r.negative_array_nt, \
            f"{n}: the length difference is not all inside the hairpin array"
    # The F30 units are drawn at the POSITIVE's period and the negatives' own
    # unit is 2 x 9 nt shorter; both decompositions of each shortfall close on
    # those two numbers, and the README explains why mixing them gives 224.
    for n in ("pART257", "pLV16"):
        r = design.loc[n]
        assert r.positive_unit_nt - r.negative_unit_nt == 2 * r.nt_per_position, \
            f"{n}: the two F30 periods no longer differ by 2 x 9 nt"
    # pLV16 is an AS-BUILT sequence and the panel prints "ordered at 8 units and
    # 16 loops; sequenced at 4 and 8" beside it. That has to be the SAME event
    # as the four unfilled blocks drawn above it, so the as-built shortfall is
    # re-closed here on the negative's own 98 nt unit.
    r = design.loc["pLV16"]
    assert r.negative_design_nt - r.negative_nt == \
        (r.negative_design_units - r.negative_units) * r.negative_unit_nt, \
        "pLV16's as-built shortfall does not close on its own 98 nt unit"
    assert r.negative_design_units == r.positive_units and \
        r.negative_design_loops == r.csy4_in_positive, \
        "pLV16's design is no longer valency-matched to pLV13"
    assert design.loc["pART256", "positive_unit_nt"] == 0, \
        "the mod4 ring has acquired a repeat period; it cannot lose a unit"
    # pART257 and pART256 ARE their own as-built record, so nothing but pLV16
    # may print an "ordered at" line.
    asbuilt = [n for n, r in design.iterrows()
               if r.negative_design_units > r.negative_units]
    assert asbuilt == ["pLV16"], f"as-built pairs are {asbuilt}, not just pLV16"


def _check_failed(sc):
    """The two scaffolds that failed are LONGER than every sibling that did not.

    Leo's note expected "how much shorter they were". The files say the other
    way round, and panel (f) prints the word "longer" beside a measured number,
    so a table edit that flipped the direction would otherwise leave the panel
    asserting the opposite of its own bars.
    """
    assert len(sc) == 12 and sc.csy4_sites.eq(sc.copies * 8).all()
    assert sorted(sc[sc.outcome == "assembly not obtained"].part) == \
        ["pART250", "pART251"], "the two failed assemblies are not the mod3 pair"
    for copies in (2, 3):
        g = sc[sc.copies == copies]
        lost, kept = g[g.outcome != "assembled"], g[g.outcome == "assembled"]
        assert len(lost) == 1 and len(kept) == 3
        assert int(lost.array_nt.iloc[0]) > int(kept.array_nt.max()), \
            f"{copies}x: the failed construct is no longer the longest"
    # The repeat unit, per architecture, from the 1x -> 2x step. mod3's is the
    # LONGEST of the four, which is the measured form of Leo's question.
    per = {a: int(g[g.copies == 2].array_nt.iloc[0] - g[g.copies == 1].array_nt.iloc[0])
           for a, g in sc.groupby("arch")}
    assert per == {"mod1": 392, "mod2": 412, "mod3": 438, "mod4": 405}, per
    assert max(per, key=per.get) == "mod3", per
    # pART254's third unit is a mod3 unit, not a mod4 one (its map name says so
    # and its 2x -> 3x step is mod3's 438, not mod4's 405). Asserted so the
    # caption cannot quietly start calling it three copies of one module.
    m4 = sc[sc.arch == "mod4"].set_index("copies").array_nt
    assert int(m4[3] - m4[2]) == per["mod3"], "pART254's third unit is not mod3's"


def stem_words(sub):
    """"4 bp stem: 1 G-C, 2 A-U, 1 G.U wobble", counted off the fold itself."""
    seq, db = "".join(sub.base), "".join(sub.dot_bracket)
    pr = pairs_of(db)
    n = {k: sum(pair_class(seq[i], seq[j]) == k for i, j in pr)
         for k in ("GC", "AU", "GU")}
    bits = [f"{n['GC']} G–C" if n["GC"] else "",
            f"{n['AU']} A–U" if n["AU"] else "",
            f"{n['GU']} G·U wobble" if n["GU"] else ""]
    return f"{len(pr)} bp stem: " + ", ".join(b for b in bits if b)


def emit(fig, out, stem):
    """Save the AUTHORED canvas, having checked nothing hangs off it.

    Every panel here goes out with `tight=False` so that all six emit the same
    85 mm box and `make_supp_book` pairs them two per row. The cost is that ink
    outside the canvas is **clipped instead of padded around**, and clipped ink
    is silent - a label simply loses its last word and the build still passes.

    So the tight bounding box is compared with the canvas before saving, and a
    panel that overflows stops the build and says by how much and on which side.
    This is the check `tight=False` needs: without it the layout fix trades a
    visible page overflow for an invisible missing label.

    The 0.05 mm tolerance is for rounding, not for slack - 0.03 mm is what
    panel (e)'s right spine sat outside by at `right=0.96`, which is a hairline
    and still a clip.
    """
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    MM = 25.4
    over = {"left": -bb.x0 * MM, "right": (bb.x1 - W) * MM,
            "bottom": -bb.y0 * MM, "top": (bb.y1 - H) * MM}
    bad = {k: round(v, 2) for k, v in over.items() if v > 0.05}
    assert not bad, (f"{stem}: ink hangs off the {W * MM:.0f} x {H * MM:.0f} mm "
                     f"canvas by {bad} mm and tight=False would CLIP it - widen "
                     f"the canvas or move the label, do not switch tight back on")
    style.save(fig, out, stem, tight=False)


def mm_axes(w, h):
    fig = plt.figure(figsize=(style.mm(w), style.mm(h)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return fig, ax


# --- (a) one hairpin position, positive against negative -------------------
def panel_a(out, design):
    """One hairpin position, positive against negative.

    Leo's "dont add to much text. thats for the caption to do" is a principle,
    not a note about panel (d), so (a) was put through the same test on
    2026-11-20: a label earns the canvas only if the DRAWING CANNOT SAY IT.

      gone    `pART252 6,194 nt, pART256 6,122 nt, 72 nt apart` - pure
              arithmetic, and the caption already carries it.
      cut     `8 positions replaced, 0 Csy4 sites left` -> `the same at all 8
              positions`. The "0 Csy4 hairpins left" half is arithmetic and goes to
              the caption. The scope does NOT: this panel draws ONE position,
              so a reader cannot count eight off the marks, and without the line
              the swap reads as happening at a single site.
      kept    `shared flank, mod4 scaffold` - a row name, and under palette v0.9
              the colour can no longer say it.
      kept    the GC line, shortened to name its window FIRST. It is the one
              place the figure says WHICH window it means - the 28 nt replaced,
              not the 16 nt hairpin core - and both figures are right under
              their own rule (README). A caption can say it; on the panel it
              stops the misreading at the moment it would happen.
    """
    # 27 mm, not the 33 it took with three lines of text under the tracks.
    # The rows did not spread: the tracks and their scale are where they were
    # and the canvas lost the line that went to the caption. Lowest ink 26.1 mm.
    W, H = 85.0, 27.0
    fig, ax = mm_axes(W, H)
    r = design.loc["pART256"]
    nt_w = 0.34                           # mm per nucleotide, one scale for both
    # x0 off the MEASURED label widths (common/style.text_w), not a round
    # number: at the old 16.0 the right-anchored "pART256, negative" ran 2.1 mm
    # off the left edge of the canvas, which `tight` absorbed by padding the
    # panel out to 92 mm - wide enough that make_supp_book could only fit ONE
    # panel per row. Every panel in this figure now emits its authored 85 mm
    # (tight=False), so ink that leaves the canvas is CLIPPED, not padded.
    labels = [f"{r.positive}, target", f"{r.name}, negative"]
    x0 = max(w_mm(t) for t in labels) + 1.5 + 0.8

    def track(y, mid, mid_col, mid_lab, label):
        x = x0
        # The flank is the mod4 scaffold itself, so it takes the mod4 colour and
        # the Csy4 site the warm end - blue scaffold, red sites, the Fig 1
        # construct-map idiom (COLORS rule 9). It also makes the family visible
        # here, not only on (c).
        for seg, col, txt in ((LEFT, MOD4, None), (mid, mid_col, mid_lab),
                              (RIGHT, MOD4, None)):
            w = len(seg) * nt_w * 3.4
            ax.add_patch(Rectangle((x, y), w, 3.4, facecolor=col,
                                   edgecolor=INK, lw=0.3, zorder=2))
            if txt:
                ax.text(x + w / 2, y + 1.7, txt, ha="center", va="center",
                        fontsize=FS, color="white", zorder=3)
            x += w
        ax.text(x0 - 1.5, y + 1.7, label, ha="right", va="center", fontsize=FS,
                color=INK)
        return x

    track(7.0, r.csy4_motif, LOOP, f"Csy4 hairpin, {int(r.csy4_nt)} nt", labels[0])
    track(14.0, r.replacement, NEUT, "", labels[1])
    ax.text(x0 + len(LEFT) * nt_w * 3.4 + len(r.replacement) * nt_w * 3.4 / 2,
            15.7, f"replacement, {int(r.replacement_nt)} nt", ha="center",
            va="center", fontsize=FS, color=INK, zorder=3)
    ax.text(x0 + len(LEFT) * nt_w * 1.7, 5.6,
            f"shared flank, {design.loc['pART256', 'arch']} array",
            ha="center", va="bottom", fontsize=FS, color=MUTED)
    # Name the WINDOW. "26% against 56%" is also right, and is the replacement
    # against the 16 nt hairpin core; this panel is about the swap, so it uses
    # the 28 nt that actually comes out, and says so rather than leaving a
    # reader to find two different GC figures for the same element.
    # The window named FIRST, which is the whole reason this line is on the
    # canvas rather than in the caption alone.
    ax.text(x0, 24.5,
            f"{int(r.csy4_nt)} nt replaced: {r.csy4_gc:.0f}% GC, "
            f"against {r.replacement_gc:.0f}%",
            ha="left", va="top", fontsize=FS, color=MUTED)
    # Scope, not arithmetic: the panel draws ONE position and the swap is the
    # same at every one of them, which no mark on this canvas can say.
    ax.text(x0, 21.0, f"the same at all {int(r.loops_in_negative)} positions",
            ha="left", va="top", fontsize=FS, color=MUTED)
    emit(fig, out, "panel_a_replacement")
    plt.close(fig)


# --- (b) the two elements, at one scale -------------------------------------
def panel_b(out, fold, xy):
    """The MFE fold DRAWN, from ViennaRNA's own naview layout, not asserted as a
    dot-bracket string. The strings this replaced were `(((((.....)))))`
    against `((((....))))`, which at a glance say the two elements fold alike -
    the opposite of the truth and of the MFE numbers beside them.

    Both folds take ONE mm-per-unit scale, so the 19 nt element is drawn smaller
    than the 28 nt it replaces. Fitting each to its own box would have hidden
    the 9 nt per site that identifies the crossed pair.
    """
    W, H = 85.0, 30.0
    fig, ax = mm_axes(W, H)
    SCALE, YAXIS = 0.235, 17.0
    # Draw first, label second. naview fans the unpaired tails into an arc
    # whose height differs between the two folds, so labels have to clear the
    # TALLER one - offsetting each from its own extent put the two captions at
    # two different heights and read as a misalignment rather than a fit.
    drawn, x0 = [], 2.0
    for name, col in (("Csy4 repeat alone", LOOP), ("replacement alone", MUTED)):
        sub = xy[xy.element == name].sort_values("i")
        w, half = draw_fold(ax, sub, x0, YAXIS, SCALE, col)
        drawn.append((name, sub, x0, half))
        x0 += w + 7.0
    hmax = max(h for *_, h in drawn)
    for name, sub, x0, _ in drawn:
        r = fold.loc[name]
        ax.text(x0, YAXIS - hmax - 3.6, name.replace(" alone", ""), ha="left",
                va="baseline", fontsize=FS, color=INK)
        ax.text(x0, YAXIS - hmax - 1.4,
                f"{int(r.nt)} nt, {r.mfe_kcal_mol:+.1f} kcal/mol", ha="left",
                va="baseline", fontsize=FS, color=MUTED)
        ax.text(x0, YAXIS + hmax + 3.0, stem_words(sub), ha="left",
                va="baseline", fontsize=FS, color=MUTED)
    emit(fig, out, "panel_b_element_fold")
    plt.close(fig)


# --- (c) the whole insert, where the negative is still folded ---------------
def panel_c(out, fold, design):
    """The construct scale. Its own panel since 2026-11-19: it used to share one
    with the element folds, and a reader who met "-0.9 kcal/mol" and a bar chart
    in the same frame could carry "unstructured" from the first to the second.
    The negative insert folds 101 pairs over 319 nt. Splitting them is the point.
    """
    W, H = 85.0, 26.0
    fig, ax = mm_axes(W, H)
    ins = ["pART252 insert", "pART256 insert", "pART158 insert", "pART257 insert"]
    y0, bh = 4.0, 3.2
    # 32, not 40: the longest bar's right label is 21.0 mm wide (measured, not
    # estimated - style.text_w), so a 40 mm bar ran it to x = 92 on an 85 mm
    # canvas and `tight` padded the panel out to 97. Ratios are unchanged,
    # every bar being scaled by the same vmax.
    xa, xw = 30.0, 32.0
    arch = arch_of(design)            # v0.9: colour no longer says which it is
    FAMILY = family(design)
    vmax = max(abs(fold.loc[i, "mfe_per_nt"]) for i in ins)
    for i, name in enumerate(ins):
        r = fold.loc[name]
        part = name.split()[0]
        y = y0 + i * (bh + 1.4)
        w = abs(r.mfe_per_nt) / vmax * xw
        ax.add_patch(Rectangle((xa, y), w, bh, facecolor=FAMILY[part],
                               edgecolor=INK, lw=0.3, zorder=2))
        ax.text(xa - 1.2, y + bh / 2, f"{part}  {arch[part]}", ha="right",
                va="center", fontsize=FS, color=INK)
        ax.text(xa + w + 1.0, y + bh / 2,
                f"{r.mfe_per_nt:.2f}   {int(r.predicted_pairs)} pairs / "
                f"{int(r.nt)} nt", ha="left", va="center", fontsize=FS,
                color=INK)
    ax.plot([xa, xa], [y0 - 0.6, y0 + 4 * (bh + 1.4) - 0.8], color=INK, lw=0.5)
    ax.text(xa, y0 + 4 * (bh + 1.4) + 1.2, "|MFE| per nucleotide (kcal/mol)",
            ha="left", va="top", fontsize=FS, color=INK)
    emit(fig, out, "panel_c_insert_fold")
    plt.close(fig)


# --- (d) the negatives that came back whole repeats short -------------------
def panel_d(out, design):
    """Each negative against its own positive, and the F30 repeats two of them
    came back short.

    Leo, 2026-11-20: "with the failed scaffolds, I meant the negatives that
    failed and had missing repeats because of the spike in fail", and then "for
    d add the negative with less copies in the 1 by 1 resulting in pLV15". That
    last one is pLV16, confirmed with him - pLV15 exists as a map but was never
    transfected, and the negative the arrayed one-by-one screen used is pLV16.
    The README says what pLV15 is and why it is not a bar here.

    THREE pairs now, at ONE mm-per-nt scale, in ASCENDING shortfall: the mod4
    ring lost nothing (it has no repeat unit to lose), pART257 came back one F30
    unit of eight short, and pLV16 four of eight. The pair that lost nothing is
    what makes a missing repeat visible at all, which is why it is still drawn.

    What a bar is: the hairpin array, first position to the end of the last,
    measured off the construct files. The whole length difference between a
    negative and its positive lies inside that span - asserted in
    `make_inputs.design()` - so one drawn quantity carries each comparison.

    WHAT IS ON THE PANEL, AND WHAT IS IN THE CAPTION. Leo, 2026-11-20: "dont
    add to much text. thats for the caption to do." He is right, and the test he
    set is a good one: a label earns its place only if the DRAWING CANNOT SAY IT.
    So the panel keeps the row names, the value labels (`888 nt  8 units, 16
    Csy4 sites`) and a two-or-three-word marker on each unfilled group. The four
    sentences that used to sit under the bars - each pair's shortfall
    arithmetic, and pLV16's ordered-against-sequenced line - are arithmetic and
    provenance, which a reader does not need in order to read the bars. They are
    in the caption draft in the README, and nothing was dropped on the way.

    What the caption carries, per `make_inputs.design()`:
      pART252 -> pART256   72 nt =                            8 x 9 nt
      pART158 -> pART257  242 nt = 1 x 116 nt F30 unit   +   14 x 9 nt
      pLV13   -> pLV16    536 nt = 4 x 116 nt F30 units  +    8 x 9 nt
    9 nt is the per-position step: the 19 nt replacement for a 28 nt repeat. The
    unfilled blocks are the POSITIVE's 116 nt unit, which is why they are drawn
    wider than the negatives' own 98 nt blocks. The other decomposition of each
    shortfall is equally true - 242 = one 98 nt negative unit + 16 x 9, and
    536 = four of them + 16 x 9 - and the README says why mixing them gives 224.

    The one thing that STAYS on the canvas is that the two F30 shortfalls are
    NOT the same event, because that distinction is what Leo asked to have
    visualised ("Exactly, this was known. Thats why i want it visualized") and a
    caption cannot do it - the reader has to be able to tell the two unfilled
    groups apart while looking at them. Cut to the shortest phrasing that still
    separates them:

      pART257 IS the seven-unit clone the spike-in screen adopted, so its own
        map is its as-built record and the eighth unit was `never placed`.
        Panel (e) is that screen: 0/56 (0%) colonies came back at eight units.
      pLV16 was ordered at eight units and sixteen loops and the sequenced
        plasmid carries four and eight, so its four are `designed, not built` -
        a 392 nt deletion, 4 x its own 98 nt unit.

    `negative_design_units` in the table is what separates the two, and the
    marker each one gets is read off that column rather than decided here.

    No cause is claimed for either, on the panel or in the caption.
    """
    # 52 mm, not the 68 it took while the four prose lines were on the canvas.
    # The rows CLOSED UP rather than spreading: bar height, scale and row
    # spacing are untouched and the canvas simply lost the band the text
    # occupied. Measured, not guessed - the lowest ink sits at 50.8 mm.
    W, H = 85.0, 52.0
    fig, ax = mm_axes(W, H)
    bh, gap, ggap, line = 3.2, 1.2, 2.4, 3.2
    FAMILY = family(design)
    # Ascending shortfall: none, one unit, four. The pair that lost nothing
    # establishes the pattern the other two depart from.
    pairs = [design.loc["pART256"], design.loc["pART257"], design.loc["pLV16"]]
    # Every mark carries its ARCHITECTURE as text, not as a hue. Under palette
    # v0.9 all four scaffold architectures are one colour (COLORS.md rule 9,
    # superseded), so the mod4 ring and the F30 array cannot be told apart by
    # fill, and the name goes in the left label where it groups the pair as
    # well. `arch` comes from the table, not from this script.
    left, detail = {}, {}
    for r in pairs:
        for role in ("positive", "negative"):
            part = r.positive if role == "positive" else r.name
            left[part] = f"{part}  {r.arch}"
            u, n = int(r[f"{role}_units"]), int(r[f"{role}_array_nt"])
            what = ("Csy4 hairpins" if role == "positive" else "loops")
            k = int(r.csy4_in_positive if role == "positive" else r.loops_in_negative)
            detail[part] = f"{n:,} nt   " + (f"{u} units, " if u else "") + f"{k} {what}"
    # Fit to MEASURED label widths, not guessed ones (common/style.text_w): the
    # bars have to clear the widest right-hand label on any row, which since the
    # third pair arrived is pLV13's, not pART158's.
    vmax = float(design.positive_array_nt.max())
    right = detail
    xa = max(w_mm(t) for t in left.values()) + 1.8
    LGAP = 1.8
    xw = W - xa - LGAP - max(w_mm(t) for t in right.values()) - 0.5
    s = xw / vmax                                     # mm per nt, one scale
    # The marker on each unfilled group - two or three words, the shortest
    # phrasing that still separates the two events. Never placed, or missing
    # from a plasmid that was designed with them? Read off the table:
    # `negative_design_units` is the negative's own ORDERED unit count, which
    # for pART257 equals its as-built count and for pLV16 does not.
    #
    # The unit count and the 116 nt are NOT repeated here - the blocks are drawn
    # one per unit and the row's value label already gives the nt, so a marker
    # that said "4 F30 units, 116 nt each" was saying three times over what the
    # drawing says once.
    block = {}
    for r in pairs:
        if r.positive_units > r.negative_units:
            block[r.name] = ("designed, not built"
                             if int(r.negative_design_units) > int(r.negative_units)
                             else "never placed")
    y = 4.6
    for r in pairs:
        # The pair's own name, from the table. Two of the three pairs are F30
        # arrays and all four designed architectures are one colour, so without
        # this the panel's six bars could be read as one series instead of three
        # pairs - and pLV13/pLV16 are the arrayed screen's constructs, not more
        # library pARTs.
        ax.text(0.6, y, r.pair, ha="left", va="top", fontsize=FS, color=MUTED)
        y += line
        top = y
        short = int(r.positive_units - r.negative_units)
        for role, col in (("positive", FAMILY[r.positive]), ("negative", NEUT)):
            part = r.positive if role == "positive" else r.name
            n = int(r[f"{role}_array_nt"])
            ax.add_patch(Rectangle((xa, y), n * s, bh, facecolor=col,
                                   edgecolor=INK, lw=0.3, zorder=2))
            # Unit boundaries where there ARE units: hairlines inside the bar,
            # so eight blocks against seven - or against four - is readable as a
            # count and not only as a length. The mod4 ring has no repeat unit
            # (its eight positions sit at irregular gaps), which is why it
            # cannot lose one. The last block of each F30 bar is SHORT because
            # the span stops at the last hairpin position, inside the final
            # unit - not a drawing error.
            per = int(r[f"{role}_unit_nt"])
            for k in range(1, int(r[f"{role}_units"] or 1)):
                ax.plot([xa + k * per * s] * 2, [y, y + bh], color=INK, lw=0.25,
                        zorder=3)
            ax.text(xa - 1.2, y + bh / 2, left[part], ha="right", va="center",
                    fontsize=FS, color=INK)
            # Both rows of a pair label in ONE column, placed off the POSITIVE's
            # length: it makes the pair read as a pair, and it keeps the
            # reference line below out of the negative's own label.
            x = xa + int(r.positive_array_nt) * s
            if role == "negative" and short:
                # The units the array came back without, at the positive's own
                # period, where they would have sat. ONE BLOCK PER UNIT, so four
                # short reads as a count and not as one long gap. No fill, no
                # hatch - absence of fill IS absence, and it survives a
                # greyscale print. Each group's own label names it, so no key.
                g = int(r.positive_unit_nt) * s
                x0 = xa + n * s
                for k in range(short):
                    ax.add_patch(Rectangle((x0 + k * g, y), g, bh,
                                           facecolor="none",
                                           edgecolor=FAMILY[r.positive], lw=0.7,
                                           zorder=2))
                # The leader sits at the FIRST missing block's left edge, not
                # the group's centre: the label is wider than the group, so a
                # centred leader came down in the middle of the string and read
                # as a stray mark. At the left edge it meets the label's own
                # start in both anchorings below.
                ax.plot([x0] * 2, [y + bh + 0.3, y + bh + 1.4],
                        color=LINEG, lw=0.3, zorder=1)
                # Anchored left at the blocks, or right at the canvas edge where
                # the label would otherwise run off it - MEASURED, not guessed
                # (common/style.text_w). Panels emit their authored 85 mm, so
                # ink that leaves the canvas is clipped.
                t = block[r.name]
                tx, ha = (x0, "left") if x0 + w_mm(t) <= W - 0.6 else (W - 0.6, "right")
                ax.text(tx, y + bh + 1.6, t, ha=ha, va="top", fontsize=FS,
                        color=MUTED)
            ax.text(x + LGAP, y + bh / 2, right[part], ha="left", va="center",
                    fontsize=FS, color=INK)
            bottom = y + bh
            y += bh + gap
        # Where the POSITIVE's array ends, carried down across its negative, so
        # each shortfall is a bounded gap and not an unanchored one. Without it
        # the unfilled units float: a reader can see they are 116 nt each but
        # not that they close 464 of pLV16's 536 and leave 72 - the 8 x 9 nt.
        ax.plot([xa + int(r.positive_array_nt) * s] * 2, [top + bh + 0.4, bottom],
                color=LINEG, lw=0.3, zorder=1)
        # One baseline per pair, not one down the whole panel: a single line ran
        # through the first pair's own summary text.
        ax.plot([xa, xa], [top - 0.6, bottom + 0.2], color=INK, lw=0.5)
        if short:
            y += 3.4                                  # the marker's own line
        y += ggap
    ax.text(xa, 3.4, "Hairpin-array length (nt)", ha="left", va="bottom",
            fontsize=FS, color=INK)
    emit(fig, out, "panel_d_missing_repeat")
    plt.close(fig)


# --- (f) the two scaffolds that were never obtained -------------------------
def deltas(sc, copies):
    """(smallest, largest) nt by which the failed construct exceeds the three
    siblings that assembled at the same copy number."""
    g = sc[sc.copies == copies]
    lost = int(g[g.outcome != "assembled"].array_nt.iloc[0])
    d = sorted(lost - int(v) for v in g[g.outcome == "assembled"].array_nt)
    return d[0], d[-1]


def panel_f(out, sc):
    """Designed Csy4-array length, 2x and 3x, with the two mod3 failures.

    Leo asked for "two failed scaffolds and how much shorter they were", and
    on 2026-11-20 said he meant the spike-in shortfall, which is now panel (d).
    This panel keeps the OTHER failure, which is a real measured result. The
    files answer in the other direction and the panel draws the measured one:
    mod3 carries the LONGEST repeat unit of the four architectures (438 nt
    against mod1's 392, mod2's 412, mod4's 405), so pART250 and pART251 are
    74 to 134 nt LONGER than the constructs that did assemble at the same copy
    number. Nothing here says that length is why they failed - n is one
    architecture, and the panel reports a size, not a cause.

    1x is left out: mod3's one-copy member pART249 had its fragment never
    delivered, so that copy number has no assembly outcome to compare with.
    Length is the DESIGNED length, from the SnapGene map, which exists for all
    twelve - including the three constructs that were never obtained.
    """
    W, H = 85.0, 45.0
    fig, ax = mm_axes(W, H)
    bh, gap, ggap = 3.0, 1.3, 3.2
    rows = [(c, a) for c in (2, 3) for a in palette.SCAFFOLD_ORDER]
    by = sc.set_index(["copies", "arch"])
    left = {k: f"{by.loc[k, 'part']}  {k[1]} {k[0]}×" for k in rows}
    # The value and the delta are two labels, not one string: the length is the
    # number Leo asked for and stays INK on every row, and only the comparison
    # beside it is an annotation. The delta names its own baseline, since the
    # three bars it is measured against are not the only bars on the panel.
    value = {k: f"{int(by.loc[k, 'array_nt']):,} nt" for k in rows}
    delta = {(c, a): "+%d to +%d nt vs mod1/2/4" % deltas(sc, c) for c, a in rows
             if by.loc[(c, a), "outcome"] != "assembled"}
    GAP = 1.8
    # Fit the bars to the MEASURED label widths rather than guessing: the left
    # column has to clear its widest pART label and the longest bar its widest
    # value, or `tight` quietly pads the panel past the 85 mm column.
    wide = {k: w_mm(value[k]) + (GAP + w_mm(delta[k]) if k in delta else 0.0)
            for k in rows}
    xa = max(w_mm(t) for t in left.values()) + 1.8
    xw = W - xa - 1.2 - max(wide.values()) - 0.8
    vmax = sc.array_nt.max()
    y = 5.0
    for i, (c, a) in enumerate(rows):
        if i == 4:
            y += ggap
        r = by.loc[(c, a)]
        lost = r.outcome != "assembled"
        w = r.array_nt / vmax * xw
        # Absence of fill IS "never obtained": texture, not colour (STYLE_NOTES
        # section 4), so it survives a greyscale print, and no hatch is used.
        # Under palette v0.8 all four architectures are one colour, so fill is
        # the ONLY thing separating these marks - which is all this panel needed
        # from the encoding anyway, and why the colour change costs it nothing.
        # Every row names its architecture in its own label.
        ax.add_patch(Rectangle((xa, y), w, bh,
                               facecolor="none" if lost else palette.SCAFFOLDS[a],
                               edgecolor=palette.SCAFFOLDS[a] if lost else INK,
                               lw=0.7 if lost else 0.3, zorder=2))
        ax.text(xa - 1.2, y + bh / 2, left[(c, a)], ha="right", va="center",
                fontsize=FS, color=INK)
        ax.text(xa + w + 1.2, y + bh / 2, value[(c, a)], ha="left",
                va="center", fontsize=FS, color=INK)
        if lost:
            ax.text(xa + w + 1.2 + w_mm(value[(c, a)]) + GAP, y + bh / 2,
                    delta[(c, a)], ha="left", va="center", fontsize=FS,
                    color=MUTED)
        y += bh + gap
    ax.plot([xa, xa], [4.4, y - gap + 0.2], color=INK, lw=0.5)
    ax.text(xa, y - gap + 1.6, "Designed Csy4-array length (nt)", ha="left",
            va="top", fontsize=FS, color=INK)
    # One frameless line, outside the bars, rather than a legend box: it is the
    # only mark on the panel that is not named by its own label.
    ax.text(W - 0.5, y - gap + 1.6, "unfilled: assembly not obtained",
            ha="right", va="top", fontsize=FS, color=MUTED)
    emit(fig, out, "panel_f_failed_scaffolds")
    plt.close(fig)


# --- (e) the colony screen -------------------------------------------------
def panel_e(out, col):
    plates = [("49761", "LV1.24, 24 colonies"), ("49780", "LV1.28, 32 colonies")]
    units = list(range(1, 10))
    fig, axes = plt.subplots(2, 1, figsize=(style.mm(85), style.mm(46)),
                             sharex=True)
    for ax, (plate, lab) in zip(axes, plates):
        d = col[col.plate == plate]
        clean = d[d.call != "ambiguous"]
        counts = [int((clean.units == u).sum()) for u in units]
        # MAIN, not palette.F30: the mark here is a COUNT OF COLONIES, a single
        # series (COLORS rule 1), not the F30 construct itself. Family colour is
        # spent where a mark IS a construct - panels (a), (c) and (d).
        ax.bar(units, counts, 0.66, color=MAIN, edgecolor=INK, linewidth=0.3)
        amb = int((d.call == "ambiguous").sum())
        for u, c in zip(units, counts):
            if c:
                ax.annotate(str(c), (u, c), textcoords="offset points",
                            xytext=(0, 1.6), ha="center", fontsize=FS)
        ax.set_ylim(0, max(counts) * 1.45)
        ax.set_yticks([0, 5, 10])
        ax.set_ylabel("Colonies")
        ax.tick_params(axis="both", length=2, width=0.5)
        ax.annotate(f"{lab}; {amb} read ambiguously", (0.5, 1.0),
                    xycoords="axes fraction", xytext=(0, 1),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=FS, color=MUTED)
        if plate == "49780":
            # "4 clones at seven units" was the bar's own value label said
            # twice - the bar is already annotated with its count - so only the
            # PROVENANCE is left, and it now names the construct rather than
            # saying "one of them". That link is why this panel is in the
            # figure at all: without it (e) is a histogram and nothing on the
            # page says which clone became the negative that panel (d) draws.
            n7 = counts[units.index(7)]
            ax.annotate("one became pART257", (7, n7),
                        textcoords="offset points", xytext=(6, 6),
                        ha="left", va="bottom", fontsize=FS, color=MUTED,
                        arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED))
    axes[-1].set_xlabel("F30 repeat units in the clone")
    axes[-1].set_xticks(units)
    # The whole point of the screen: nobody landed on eight.
    for ax in axes:
        ax.axvspan(7.5, 8.5, color=BGFILL, zorder=0)
    # -26, not -16: since the axis label became "F30 repeat units in the
    # clone" (2026-11-23) the two sit on one line with 2.5 pt between them and
    # read as one run-on sentence. A second line keeps the tick association.
    axes[-1].annotate("eight units, the target", (8, 0),
                      xycoords=("data", "axes fraction"),
                      textcoords="offset points", xytext=(0, -26), ha="center",
                      va="top", fontsize=FS, color=MUTED, annotation_clip=False)
    # right=0.955, not 0.99: at 0.99 the axes' own right spine sat 1.8 mm beyond
    # the 85 mm canvas, which tight absorbed as padding and tight=False clips.
    # 0.96 left 0.03 mm still outside - a hairline, and still a clip, which is
    # why `emit()` checks rather than trusting the arithmetic.
    fig.subplots_adjust(left=0.13, right=0.955, top=0.92, bottom=0.26, hspace=0.34)
    emit(fig, out, "panel_e_colony_screen")
    plt.close(fig)


def report(design, fold, col, sc):
    print("--- what the negative puts where the Csy4 repeat was ---")
    print(design.reset_index()[["pair", "positive", "negative", "positive_nt",
                                "negative_nt", "delta_nt", "csy4_in_positive",
                                "csy4_in_negative", "loops_in_negative"]]
          .to_string(index=False))
    print("\n--- ViennaRNA MFE fold ---")
    print(fold[["nt", "predicted_pairs", "mfe_kcal_mol", "mfe_per_nt"]].to_string())
    print("\n--- panel (d): the repeat the spike-in did not place ---")
    for n, r in design.iterrows():
        u = int(r.positive_units - r.negative_units)
        print(f"  {r.positive} ({r.arch}) {int(r.positive_array_nt):,} nt array, "
              f"{int(r.positive_units)} units, {int(r.csy4_in_positive)} Csy4 hairpins"
              f"  ->  {n} {int(r.negative_array_nt):,} nt, "
              f"{int(r.negative_units)} units, {int(r.loops_in_negative)} loops")
        print(f"     {int(r.delta_nt)} nt shorter than {r.positive} = "
              + (f"{u} x {int(r.positive_unit_nt)} nt unit + " if u else "")
              + f"{int(r.loops_in_negative)} x {int(r.nt_per_position)} nt"
              + ("" if u else ", no unit lost"))
        if r.negative_design_units > r.negative_units:
            print(f"     {n} is an AS-BUILT sequence: ordered at "
                  f"{int(r.negative_design_nt):,} nt, "
                  f"{int(r.negative_design_units)} units, "
                  f"{int(r.negative_design_loops)} loops, from "
                  f"{r.source_negative_design}; the "
                  f"{int(r.negative_design_nt) - int(r.negative_nt)} nt it came "
                  f"back short is "
                  f"{int(r.negative_design_units) - int(r.negative_units)} x "
                  f"{int(r.negative_unit_nt)} nt")
    print("\n--- the twelve designed scaffolds, and which were obtained ---")
    print(sc[["part", "arch", "copies", "plasmid_nt", "csy4_sites", "array_nt",
              "outcome"]].to_string(index=False))
    per = {a: int(g[g.copies == 2].array_nt.iloc[0] - g[g.copies == 1].array_nt.iloc[0])
           for a, g in sc.groupby("arch")}
    print("  repeat unit per architecture (nt):", per,
          "-> longest:", max(per, key=per.get))
    n = sc.outcome.value_counts().to_dict()
    tot = len(sc)
    for k, v in n.items():
        print(f"  {k}: {v}/{tot} ({100 * v / tot:.0f}%)")
    for copies in (2, 3):
        g = sc[sc.copies == copies]
        lo, hi = deltas(sc, copies)
        f = g[g.outcome != "assembled"].iloc[0]
        print(f"  {copies}x: {f.part} ({f.arch}) {int(f.array_nt)} nt is "
              f"+{lo} to +{hi} nt LONGER than the three that assembled "
              f"({', '.join(str(int(v)) for v in sorted(g[g.outcome == 'assembled'].array_nt))})")
    print("\n--- colony screen ---")
    for plate in ("49761", "49780"):
        d = col[col.plate == plate]
        clean = d[d.call != "ambiguous"]
        hist = clean.units.value_counts().sort_index().to_dict()
        print(f"  plate {plate} ({d.journal_entry.iloc[0]}): n = {len(d)}, "
              f"units {hist}, ambiguous {int((d.call == 'ambiguous').sum())}")
        eight = clean[clean.units == 8]
        print(f"     colonies at exactly eight units: {len(eight)}")
        seven = clean[clean.units == 7]
        if len(seven):
            print(f"     seven-unit colonies: {sorted(seven.colony)}")
        if len(clean):
            top = clean.loc[clean.units.idxmax()]
            print(f"     largest clean clone: colony {int(top.colony)}, "
                  f"{int(top.units)} units ({int(top.loops)} loops)")
    n8 = int((col[col.call != 'ambiguous'].units == 8).sum())
    print(f"  eight-unit clones across all {len(col)} colonies: "
          f"{n8}/{len(col)} ({100 * n8 / len(col):.0f}%)")


def check_emitted(out):
    """The six panels have to come out as THREE ROWS OF TWO in the supp book.

    Leo, 2026-11-20: "Align the panels. Should be two per horizontal level
    always." The book places panels by their EMITTED size, so this is not
    checked by eye or by what the canvases were authored at - it is re-derived
    with `make_supp_book.rows_for` and `make_supp_book.svg_size`, the same two
    functions the book itself uses, reading the SVGs this run just wrote.

    It is the only check here that can fail on a change that leaves every number
    in the figure correct: a longer tick label, one more line of annotation, or
    a `tight=True` slipping back in pads a panel past 88 mm, and the book then
    silently drops to one panel per row and overflows the page.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import make_supp_book as book
    panels = sorted(q for q in out.glob("*.svg") if not q.name.startswith("._"))
    sizes = [book.svg_size(q) for q in panels]
    avail = book.A4_W - 2 * book.MARGIN
    rows = book.rows_for(sizes, avail)
    y = book.MARGIN + book.TITLE_H
    print("\n--- page layout, from make_supp_book's own rows_for ---")
    for row in rows:
        h = max(sizes[i][1] for i in row)
        y += h + book.GUTTER
        print("  " + " | ".join(f"{panels[i].stem[6:]} {sizes[i][0]:.1f} x "
                                f"{sizes[i][1]:.1f} mm" for i in row))
    body = book.A4_H - book.MARGIN
    print(f"  {len(panels)} panels in {len(rows)} rows; content runs "
          f"{y:.0f} mm of the {body:.0f} mm page body")
    # ONE emitted width, and two of them fit a row. Not "== 85.0": matplotlib
    # rounds the SVG canvas to whole pixels at the save dpi, so an 85 mm canvas
    # comes out at 84.8 mm (85/25.4 x 100 = 334.6 px -> 334 px -> 240.5 pt).
    # What the layout needs is that the six agree with each other and that a
    # pair clears the text block, which is what is asserted.
    w = {round(x, 1) for x, _ in sizes}
    assert len(w) == 1, f"panels emit different widths: {sorted(w)}"
    pair = 2 * w.pop() + book.GUTTER
    print(f"  one emitted width; a pair plus the {book.GUTTER:.0f} mm gutter is "
          f"{pair:.1f} mm of the {avail:.0f} mm text block")
    assert pair <= avail, f"two panels plus a gutter is {pair:.1f} mm, over {avail}"
    assert all(len(r) == 2 for r in rows), \
        f"not two per row: {[[panels[i].stem for i in r] for r in rows]}"
    assert y <= body, f"content runs {y:.0f} mm, past the {body:.0f} mm page body"


def main():
    style.apply()
    out = paths.fig_output(FIG)
    design, fold, col, xy, sc = load()
    assert len(col) == 56
    _selfcheck(design, fold, xy, sc)
    panel_a(out, design)
    panel_b(out, fold, xy)
    panel_c(out, fold, design)
    panel_d(out, design)
    panel_e(out, col)
    panel_f(out, sc)
    report(design, fold, col, sc)
    check_emitted(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
