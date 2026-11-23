"""RNA scaffold design: the F30 starting point, four architectures, their
dot-bracket targets, and the design space with its fold checks.

Four bare panels. Every number is re-derived by `make_inputs.py` from the
SnapGene maps, the dot-bracket targets and a ViennaRNA re-fold on Euler; this
script only draws what those tables say.

All five scaffolds - the F30 unit and the four designed architectures - are
drawn the SAME way, as folded secondary structures laid out by ViennaRNA's
naview, the layout RNAplot itself uses: stems as ladders, loops opened out.
Panel (a) is the F30 unit as the eight-unit array of pLV009 actually folds it;
panel (c) is one repeat of each dot-bracket design target. Both come from
`inputs/supp_scaffold_design_layouts.csv`, so nothing here invents a fold.
(Panel (c) used to be an arc diagram and panel (a) a hand-drawn junction;
replaced 2026-11-19 at Leo's request - "all 5 scaffolds ... use the dotbrace
representation to get the nice folds".)

Colour is palette v0.9. `palette.SCAFFOLDS` is ONE colour for all four
architectures - theme 600 `#0678a6` since 2026-11-20 (Leo: "adjust the rna
presentation modules blue to a lighter one"), `MAIN` theme 700 under v0.8, and
the four-step indigo ramp of v0.7 before that having been withdrawn - so
colour no longer separates them anywhere in this figure - the drawn fold and the
printed NAME do, on (b), (c) and (d) alike. The F30 array keeps `palette.F30`,
the blue the Fig 1 construct maps give it, so panel (a) still reads as F30 and
not as a designed architecture. Csy4 sites are the warm pole of the diverging
ramp in both (COLORS.md rule 4) and base-pair rungs stay neutral, structure
inside a cartoon never being coloured (rule 8). Panel (b)'s glyphs still come
from `common/scaffold_glyphs.py`, the module Fig 1f and Fig 4a share.

Panel (d)'s twelve cells are the twelve ORDERED constructs, folded at their real
copy number: the D-series of 22-JUL-2026 (D1=mod1 ... D4=mod4), where a
construct at k copies is the first k of its three DISTINCT fragments assembled
head to tail. Nothing is truncated and no unit is repeated - an earlier version
of this panel cut a 3x design into thirds, which both used the superseded
first-generation targets and reconstructed a repeat that does not exist.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from common import paths, style, palette, scaffold_glyphs

FIG = Path(__file__).resolve().parent.name

INK = palette.INK
MUTED = palette.NEUTRAL["600"]
LINEG = palette.NEUTRAL["400"]
LOOP = palette.DIVERGING[-1]          # the Csy4 hairpin tip, as in Fig 1f
PAIR = palette.NEUTRAL["300"]         # base-pair rungs stay neutral (COLORS r8)
TEAL = palette.THEME["700"]
BGFILL = palette.NEUTRAL["100"]
FS = 5.5

# Display names. RENAMED 2026-11-20, when the panels moved to the ORDERED
# designs: every arm is 5 bp there (the low-structure redesign's budget of at
# most 30 bp per 120 nt window - `csy4_redesigns_lowstructure.txt`, "Tandem 4
# modules, stems 5/5/5"), so "long arms" / "short arms" described the
# SUPERSEDED generation and, worse, inverted on the new one: mod3 now draws
# LARGER than mod1. What separates them now is the linker between modules
# (`41 57 ...` against `38 75 ...`) and the exterior loop it opens, which is
# what the design files themselves call it - `D1_tandem_lean` against
# `D3_splayed_tandem`. The tables, the files and the journal keep mod1..mod4, so
# the token stays beside the name rather than replacing it.
ARCH = ["mod1", "mod2", "mod3", "mod4"]
NAME = {"mod1": "tandem, compact",
        "mod2": "branched",
        "mod3": "tandem, splayed",
        "mod4": "ring"}


def _mm(pt):
    """Points to mm, for laying a panel out from MEASURED type widths."""
    return pt * 25.4 / 72


def load():
    space = pd.read_csv(paths.INPUTS / "supp_scaffold_design_space.csv")
    targ = pd.read_csv(paths.INPUTS / "supp_scaffold_design_targets.csv")
    folds = pd.read_csv(paths.INPUTS / "supp_scaffold_design_folds.csv")
    lay = pd.read_csv(paths.INPUTS / "supp_scaffold_design_layouts.csv")
    cells = pd.read_csv(paths.INPUTS / "supp_scaffold_design_ordered.csv")
    layo = pd.read_csv(paths.INPUTS / "supp_scaffold_design_layouts_ordered.csv")
    targ["arch"] = targ["target"].str.extract(r"(mod\d)")
    assert len(cells) == 12, f"panel (d) needs 12 cells, got {len(cells)}"
    return (space, targ.set_index("arch"), folds.set_index("run"),
            dict(tuple(lay.groupby("scaffold"))), cells.set_index(["arch", "copies"]),
            dict(tuple(layo.groupby("scaffold"))))


def mm_axes(w, h):
    fig = plt.figure(figsize=(style.mm(w), style.mm(h)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)            # y grows downward, as the shared glyphs assume
    ax.axis("off")
    return fig, ax


# --- the shared secondary-structure drawing --------------------------------
# One representation for all five scaffolds. Three strokes and no more: the
# backbone, the base-pair rungs that make a stem read as a stem, and the Csy4
# repeats drawn over the backbone in the warm colour. No per-nucleotide dots -
# at 400-odd nt in 40 mm they are a texture, not information.
SCALE = 0.0395          # mm per naview unit, set by the widest target (mod2)


def _runs(flag):
    """[(start, stop)] of each contiguous run of 1s."""
    out, i = [], 0
    while i < len(flag):
        if flag[i]:
            j = i
            while j < len(flag) and flag[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def draw_fold(ax, g, cx, cy, strand, s=SCALE, lw=0.5):
    """Draw one fold centred on (cx, cy) at `s` mm per layout unit.

    Returns the placed coordinates, so a caller can hang a leader off a
    nucleotide it cares about.
    """
    gx, gy = g["x"].to_numpy(), g["y"].to_numpy()
    X = cx + (gx - (gx.min() + gx.max()) / 2) * s
    Y = cy - (gy - (gy.min() + gy.max()) / 2) * s      # the axes are y-down
    pair = g["pair"].to_numpy()
    for i, j in enumerate(pair):
        if j > i:
            ax.plot([X[i], X[j]], [Y[i], Y[j]], color=PAIR, lw=lw * 0.55,
                    solid_capstyle="butt", zorder=1)
    ax.plot(X, Y, color=strand, lw=lw, solid_capstyle="round", zorder=2)
    for a, b in _runs(g["csy4"].to_numpy()):
        ax.plot(X[a:b], Y[a:b], color=LOOP, lw=lw * 1.8, solid_capstyle="round",
                zorder=3)
    return X, Y


def fold_size(g, s=SCALE):
    return ((g["x"].max() - g["x"].min()) * s, (g["y"].max() - g["y"].min()) * s)


# --- (a) the F30 three-way junction we started from ------------------------
def panel_a(out, lay):
    """One F30 unit as the pLV009 array folds it: a 3WJ with two Csy4 sites.

    Two sites per unit is what makes pART158's eight tandem units sixteen
    sites, and it is the only place in the figure where that ratio is drawn.
    Drawn at twice panel (c)'s scale - it is 138 nt against their 400-odd, and
    at one scale it would be a smudge.
    """
    W, H = 66.0, 50.0
    fig, ax = mm_axes(W, H)
    g = lay["f30"]
    # palette.F30, not a SCAFFOLDS step: the F30 array is the thing the designs
    # START from, drawn as the Fig 1 construct maps draw it (COLORS.md rule 9).
    X, Y = draw_fold(ax, g, 26.0, 24.0, palette.F30, s=SCALE * 2, lw=0.75)

    # 5' and 3' are the two ends of the SAME closing helix, a base pair apart,
    # so a fixed left/right offset lands them on top of each other. Push each
    # one along the line joining them (away from the other) and a little
    # further along its own last step (out of the helix).
    n = len(X) - 1
    for i, k, lab in ((0, 1, "5′"), (n, n - 1, "3′")):
        vx, vy = X[i] - X[n - i], Y[i] - Y[n - i]
        m = math.hypot(vx, vy) or 1.0
        ox, oy = X[i] - X[k], Y[i] - Y[k]
        mo = math.hypot(ox, oy) or 1.0
        ax.text(X[i] + 2.0 * vx / m + 1.6 * ox / mo,
                Y[i] + 2.0 * vy / m + 1.6 * oy / mo,
                lab, ha="center", va="center", fontsize=FS - 0.8, color=MUTED)

    # one leader, onto the midpoint of the right-hand Csy4 repeat
    runs = _runs(g["csy4"].to_numpy())
    r = max(runs, key=lambda ab: X[ab[0]:ab[1]].mean())
    m = (r[0] + r[1]) // 2
    ax.annotate("Csy4 hairpin", xy=(X[m], Y[m]), xytext=(X[m] + 6.5, Y[m] + 4.0),
                ha="left", va="center", fontsize=FS, color=MUTED,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED,
                                shrinkA=0, shrinkB=2))
    ax.text(26.0, H - 2.0, f"{len(X)} nt, 2 sites", ha="center", va="bottom",
            fontsize=FS, color=MUTED)
    style.save(fig, out, "panel_a_f30_unit")
    plt.close(fig)


# --- (b) the four architectures --------------------------------------------
def panel_b(out, cells, layo):
    """The four architectures as glyphs, on the ORDERED generation.

    Brought over from the superseded first generation on 2026-11-20, with (c).
    Two things changed that (c) did not need. The per-repeat length now comes
    from the ordered table, not `supp_scaffold_design_targets.csv`; and mod2's
    branching is parsed from the ORDERED target, whose tree really is different
    - first generation `(8, [(7, ...), (6, ...)])` against the ordered
    `(4, [(5, ...), (4, ...)])`, so the glyph drew the wrong number of stacked
    base pairs per branch.

    ⚠️ mod1 and mod3 now draw the SAME glyph. `stem` is what used to separate
    them, and in the ordered designs their arms are both 5 bp, so asserting a
    difference here would be drawing something that is not true. What separates
    them now is how far apart the modules sit, which `scaffold_glyphs.array`
    has no parameter for - adding one is a change to a module Fig 1f and Fig 4a
    share, so it is flagged in the README rather than made here.
    """
    W, H = 178.0, 34.0
    fig, ax = mm_axes(W, H)
    slot = W / 4
    k = 0.70
    for i, a in enumerate(ARCH):
        gx = slot * (i + 0.5)
        base = 20.0
        n0, p0 = len(ax.lines), len(ax.patches)
        col = palette.SCAFFOLDS[a]
        if a in ("mod1", "mod3"):
            half = 2.5 * k * 3 + 1.4
            scaffold_glyphs.array(ax, gx, base + 4.5 * k, col, LOOP,
                                  gx - half, gx + half, n=4,
                                  stem=2.7,          # 5 bp arms in BOTH, now
                                  lw=1.0, k=k, dot_r=0.8 * k)
        elif a == "mod2":
            scaffold_glyphs.tree(ax, gx, base + 4.5 * k, col, LOOP,
                                 gx - 11, gx + 11, mod2_tree(layo), lw=1.0,
                                 k=k, dot_r=0.8 * k)
        else:
            scaffold_glyphs.ring(ax, gx, base + 3.5 * k, col, LOOP,
                                 gx - 9, gx + 9, n=8, r=3.2, lw=1.0, k=k,
                                 dot_r=0.8 * k)
        top = min([min(l.get_ydata()) for l in list(ax.lines)[n0:]]
                  + [c.center[1] - c.radius for c in list(ax.patches)[p0:]])
        ax.text(gx, top - 1.4, a, ha="center", va="bottom", fontsize=FS,
                color=INK)
        ax.text(gx, base + 6.0, NAME[a], ha="center", va="top", fontsize=FS,
                color=MUTED)
        ax.text(gx, base + 9.3,
                f"{int(cells.loc[(a, 1), 'nt'])} nt, 8 sites",
                ha="center", va="top", fontsize=FS, color=MUTED)
    style.save(fig, out, "panel_b_architectures")
    plt.close(fig)


def mod2_tree(layo):
    """mod2's real branching, from the ORDERED 1x design target.

    Read back out of the layout table's `pair` column rather than from a
    separate file, so the glyph and the fold drawn in (c) are the same
    structure by construction.
    """
    g = layo["mod2_target"].sort_values("i")
    unit = "".join("(" if q > i else ")" if 0 <= q < i else "."
                   for i, q in zip(g["i"], g["pair"]))
    st, pt = [], [-1] * len(unit)
    for i, c in enumerate(unit):
        if c == "(":
            st.append(i)
        elif c == ")":
            j = st.pop()
            pt[i] = j
            pt[j] = i

    def kids(i, j):
        out, m = [], i + 1
        while m < j:
            if pt[m] > m:
                out.append((m, pt[m]))
                m = pt[m] + 1
            else:
                m += 1
        return out

    def node(i, j):
        n, a, b = 0, i, j
        while a < b and pt[a] == b:
            n += 1
            a += 1
            b -= 1
            while a < b and unit[a] == ".":
                a += 1
            while a < b and unit[b] == ".":
                b -= 1
        ch = [c for c in kids(a - 1, b + 1) if c != (i, j)] if a <= b else []
        return (n, [node(u, v) for u, v in ch])

    top = kids(-1, len(unit))
    return node(*top[0])


# --- (c) the ordered fragments: design target against RNAfold prediction ---
def panel_c(out, cells, layo):
    """Each ordered 1x fragment twice: what it was designed to be, and what
    ViennaRNA predicts it folds into.

    Replaced the first-generation panel on 2026-11-20 (Leo: *"Can you recreate
    the c panel (dot brace based) with the real ordered designs?"*), then gained
    the second row (*"we also have the actual predictions of the sequences by
    viennafold. Can you add them to figure."*). The retired first-generation
    drawing is in `archive/`.

    **Both rows come from one `RNA.fold` call in `make_inputs.py`**, the same
    call whose distance panel (d) reports, and that distance is asserted
    against `supp_scaffold_design_ordered.csv` when the layout is built - so the
    two panels cannot disagree about what a sequence folds into.

    Layout: architectures in columns, target above and prediction below, all
    EIGHT folds at one `SCALE` - unchanged from the rest of the figure, which
    is why the panel is 88 mm tall rather than redrawn smaller. One scale
    matters more here than anywhere else in the figure, because the comparison
    being asked for is structural: a target and a prediction drawn at different
    sizes cannot be read against each other.

    🚨 **Area is not magnitude, and the caption has to say so.** naview places
    an exterior loop from its pairing context, so a handful of pairs opening or
    closing near a multiloop reorganises the WHOLE drawing. mod3's two rows
    look dramatically different where mod1's look near-identical, and the
    measured distances are 11 / 240 (4.6%) against 2 / 240 - a reader given the
    drawings alone takes mod3 as catastrophically wrong. The panel shows WHERE
    a design departs from its target; the `bp off` figure under each column is
    what says HOW MUCH. Do not drop that sentence from the caption, and do not
    "fix" this by scaling or re-laying-out the folds: the layout is the thing
    that reorganises, not the molecule. mod3's eleven pairs are in the linker
    between modules - all eight Csy4 hairpins form in the prediction, asserted
    in `main()`, as they do in all twelve constructs of panel (d). The reading
    the panel supports is "the linker reorganises, the sites do not".

    mod2 and mod4 are at distance 0, so their two rows are the SAME drawing.
    They are still drawn twice: a 2 x 4 grid with two cells left empty reads as
    missing data rather than as agreement, and the coincidence is the check
    that the exact cells really are exact. The per-column label carries the
    distance so a reader knows what to look for.
    """
    W, H = 178.0, 88.0
    fig, ax = mm_axes(W, H)
    ROWS = (("target", "design target", 20.0),
            ("mfe", "RNAfold prediction", 58.0))
    gut = _mm(max(style.text_w(lab, FS) for _, lab, _ in ROWS)) + 3.0
    col = (W - gut) / len(ARCH)
    assert max(fold_size(layo[f"{a}_{k}"])[0] for a in ARCH
               for k, _, _ in ROWS) < col - 2.0, "folds overlap"
    LABEL_TOP = H - 9.6
    for a in ARCH:
        bottom = ROWS[-1][2] + fold_size(layo[f"{a}_{ROWS[-1][0]}"])[1] / 2
        assert LABEL_TOP - bottom > 1.2, (
            f"{a}: fold reaches {bottom:.2f} mm, label baseline {LABEL_TOP:.2f} mm")
    for k, lab, cy in ROWS:
        # direct row label at the left, not a legend box; a reader must never
        # have to guess which drawing is the design and which the prediction
        ax.text(gut - 3.0, cy, lab, ha="right", va="center", fontsize=FS,
                color=MUTED)
        for i, a in enumerate(ARCH):
            draw_fold(ax, layo[f"{a}_{k}"], gut + col * (i + 0.5), cy,
                      palette.SCAFFOLDS[a])
    for i, a in enumerate(ARCH):
        cx = gut + col * (i + 0.5)
        r = cells.loc[(a, 1)]
        ax.text(cx, LABEL_TOP, a, ha="center", va="baseline", fontsize=FS,
                color=INK)
        ax.text(cx, H - 6.1, NAME[a], ha="center", va="baseline", fontsize=FS,
                color=MUTED)
        ax.text(cx, H - 2.6,
                f"{int(r.nt)} nt, 8 sites \u00b7 {int(r.bp_distance)} bp off",
                ha="center", va="baseline", fontsize=FS, color=MUTED)
    style.save(fig, out, "panel_c_ordered_fragments")
    plt.close(fig)


# --- (d) the design space, and at what scale each fold was checked ----------
def panel_d(out, space, cells):
    """Architecture against copy number, every cell carrying both facts.

    The fold check used to sit in a second tier below the grid, with its own
    copy of the four row labels and columns that did not line up with the ones
    above. It is now the second line of each cell, so the scale a fold was
    checked at IS the column it sits in, and the panel cannot misregister.

    All twelve cells carry a distance, and every one is a REAL ORDERED
    construct folded at its real copy number (2026-11-20). The cells used to
    carry five (a cell matched a ROAD run by name), then twelve truncations of a
    3x design - wrong twice over, because those were the superseded 20-JUL
    targets and because a 3x construct is three DISTINCT fragments, not one unit
    tripled. `make_inputs.py` now assembles the D-series GenBank files and folds
    f1, f1+f2 and f1+f2+f3 against the matching prefix of the design target.

    ⚠️ **The panel has to explain itself** (Leo, 2026-11-20: "explain what panel
    d does"). It is met before its caption, and a bare grid of `2 / 240` is
    unreadable. Three things were added, all as headings and keys rather than
    as prose on the canvas (STYLE_NOTES section 0 keeps notes in the caption,
    but an axis heading is not a note):

    * the two grid axes are NAMED - `architecture` over the row labels,
      `fragments assembled` over the columns - so the reader is told what the
      12 cells are indexed by;
    * both cell lines are keyed as numerator / denominator beside the first
      row, `Csy4 hairpins formed / designed` and `bp off target / 2 x target
      bp`, which says what the slash means instead of leaving it to be guessed;
    * the TOP line now draws `hairpins_ok / n_csy4` (8/8, 16/16, 24/24) where
      it used to draw the site count alone. That turns the figure's actual
      result into something on the panel: **every Csy4 hairpin forms in all
      twelve**, so where mod1 and mod3 depart from target it is the scaffold
      between the sites and never a site. `main()` asserts the two columns are
      equal, so the key cannot outlive the fact.

    `cw` dropped 14.0 -> 13.0 to pay for the wider right-hand key: the gutter
    needs 25.1 mm for `Csy4 hairpins formed / designed` at 5 pt (MEASURED, and
    asserted below), and at cw 14.0 only 24.2 mm was left.
    """
    W, H = 85.0, 42.0
    fig, ax = mm_axes(W, H)
    cols = [1, 2, 3]
    cw, ch = 13.0, 7.4
    # The row label now carries the architecture's NAME as well as its modN
    # token, so the gutter is set from the MEASURED widest name - 0.52 em per
    # character is out by a tenth to a third per string (memory: measure, do
    # not estimate text).
    gx0 = _mm(max(style.text_w(NAME[a], FS) for a in ARCH)) + 2.9
    gy0 = 10.0

    # pART254's architecture reads "mod4+mod3": it is the mod4 row's 3-copy
    # construct with a mod3 unit in its middle position (see README).
    built = {(r.architecture.split("+")[0], r.copies): r
             for _, r in space.iterrows()}
    for j, c in enumerate(cols):
        ax.text(gx0 + cw * (j + 0.5), gy0 - 1.4, f"{c}×", ha="center",
                va="bottom", fontsize=FS, color=INK)
    # The two grid axes, named. Without them the reader is given twelve boxes
    # of numbers and left to infer what indexes them.
    ax.text(gx0 + cw * 1.5, gy0 - 4.6, "fragments assembled, 8 Csy4 hairpins each",
            ha="center", va="bottom", fontsize=FS, color=MUTED)
    ax.text(gx0 - 1.5, gy0 - 1.4, "architecture", ha="right", va="bottom",
            fontsize=FS, color=MUTED)

    for i, a in enumerate(ARCH):
        y = gy0 + ch * i
        ax.text(gx0 - 1.5, y + ch / 2 - 1.2, a, ha="right", va="center",
                fontsize=FS, color=INK)
        ax.text(gx0 - 1.5, y + ch / 2 + 1.2, NAME[a], ha="right", va="center",
                fontsize=FS, color=MUTED)
        for j, c in enumerate(cols):
            x = gx0 + cw * j
            r = built[(a, c)]
            cloned = r.cloned == "yes"
            # A dashed cell is a map that was designed but never cloned; the
            # caption names it. A diagonal strike used to do this and ran
            # straight through the numeral.
            ax.add_patch(Rectangle(
                (x + 0.5, y + 0.6), cw - 1.0, ch - 1.2, facecolor=BGFILL,
                edgecolor=LINEG, lw=0.4, zorder=1,
                linestyle="solid" if cloned else (0, (1.6, 1.2))))
            f = cells.loc[(a, c)]
            # Hairpins FORMED over sites designed, not the site count alone:
            # 8/8, 16/16, 24/24 everywhere, which is the panel's result.
            ax.text(x + cw / 2, y + 3.2,
                    f"{int(f.hairpins_ok)} / {int(f.n_csy4)}", ha="center",
                    va="center", fontsize=FS, color=INK if cloned else MUTED)
            ax.text(x + cw / 2, y + 5.7,
                    f"{int(f.bp_distance)} / {2 * int(f.target_bp)}",
                    ha="center", va="center", fontsize=FS, color=TEAL)

    # The key for both cell lines, on the first row, rather than a legend box.
    # Each one says what its slash means, in the same numerator / denominator
    # form the cell is printed in - so "8 / 8" and "2 / 240" are read off the
    # panel instead of off the caption. FS - 0.5: MEASURED at 25.1 mm for the
    # longer of the two against a 27.2 mm gutter (asserted).
    xr = gx0 + cw * len(cols) + 1.8
    keys = [(3.2, "Csy4 hairpins formed / designed"),
            (5.7, "bp off target / 2 × target bp")]
    for dy, lab in keys:
        assert xr + _mm(style.text_w(lab, FS - 0.5)) < W, f"clipped: {lab}"
        ax.text(xr, gy0 + dy, lab, ha="left", va="center", fontsize=FS - 0.5,
                color=MUTED)
    style.save(fig, out, "panel_d_design_space")
    plt.close(fig)


def report(space, targ, folds, lay, cells, layo):
    print("\n--- drawn folds (naview layout) ---")
    for n in ("f30", *ARCH):
        g = lay[n]
        w, h = fold_size(g)
        print(f"{n:5s} {len(g):4d} nt  {int(g['pair'].gt(g['i']).sum()):3d} bp  "
              f"{len(_runs(g['csy4'].to_numpy()))} Csy4 sites  "
              f"{w:.1f} x {h:.1f} mm at SCALE")
    print("\n--- ORDERED 1x fragments, as panel (c) draws them ---")
    for a in sorted(layo):
        g = layo[a]
        w, h = fold_size(g)
        print(f"{a:5s} {len(g):4d} nt  {int(g['pair'].gt(g['i']).sum()):3d} bp  "
              f"{len(_runs(g['csy4'].to_numpy()))} Csy4 sites  "
              f"{w:.1f} x {h:.1f} mm at SCALE")
    print("\n--- design space (SnapGene maps) ---")
    print(space[["construct", "architecture", "copies", "insert_nt", "n_csy4",
                 "cloned"]].to_string(index=False))
    print("\n--- dot-bracket targets ---")
    print(targ[["target", "total_nt", "unit_nt", "target_bp", "n_csy4",
                "top_level_helices"]].to_string())
    print("\n--- full-length fold check ---")
    for a in ARCH:
        r = folds.loc[f"F30_{a}_pos"]
        print(f"F30_{a}_pos  {int(r.target_nt)} nt  {int(r.target_bp)} target bp  "
              f"bp distance {int(r.bp_distance)}  {r.pct_pairs}% of pairs  "
              f"{int(r.n_csy4)} Csy4")
    print("\n--- 1x fragment runs present on Euler ---")
    frag = folds[folds["target_nt"] < 600]
    print(frag[["best_nt", "target_bp", "bp_distance", "pct_pairs", "n_csy4",
                "n_neg_loop", "mfe"]].to_string())
    print("\n--- panel (d): the twelve ORDERED constructs ---")
    print(cells[["nt", "target_bp", "bp_distance", "pct_pairs", "n_csy4",
                 "hairpins_ok", "kind"]].to_string())
    bad = cells[cells["hairpins_ok"] != cells["n_csy4"]]
    print(f"\nCsy4 hairpins formed in the MFE: "
          f"{'ALL, in all twelve' if bad.empty else bad.index.tolist()}")
    print("engine:", sorted(set(cells["engine"])))


def main():
    style.apply()
    out = paths.fig_output(FIG)
    space, targ, folds, lay, cells, layo = load()
    assert len(space) == 12, "12 designed constructs expected"
    assert set(lay) == {"f30", *ARCH}, f"layouts has {sorted(lay)}"
    assert set(layo) == {f"{a}_{k}" for a in ARCH for k in ("target", "mfe")}, \
        f"ordered layouts has {sorted(layo)}"
    # A palette edit must not silently re-key the architectures. Under palette
    # v0.9 the four still share ONE colour - only its step moved, MAIN -> theme
    # 600 - so the printed name is the ONLY thing telling them apart; assert
    # every architecture has one.
    assert palette.SCAFFOLD_ORDER == ARCH, palette.SCAFFOLD_ORDER
    assert set(NAME) == set(ARCH)
    # Panel (d) draws "hairpins formed / designed" and the whole point of that
    # line is that the two are equal in all twelve - it is what lets a reader
    # see that mod1's and mod3's drift is in the scaffold BETWEEN the sites and
    # never at a site. Assert it, so a regenerated table cannot leave the panel
    # claiming something the data stopped saying.
    bad = cells[cells["hairpins_ok"] != cells["n_csy4"]]
    assert bad.empty, f"a Csy4 hairpin failed to fold:\n{bad}"
    # The Csy4 runs are what the warm colour draws, so count them here: a
    # broken _runs(), or a layout table rebuilt against the wrong sequence,
    # would otherwise show up only as a fold with the colour in odd places.
    for n, want in (("f30", 2), *((a, 8) for a in ARCH)):
        got = len(_runs(lay[n]["csy4"].to_numpy()))
        assert got == want, f"{n}: {got} Csy4 runs, expected {want}"
    for a in ARCH:                       # the ordered fragments carry 8 too
        for k in ("target", "mfe"):
            got = len(_runs(layo[f"{a}_{k}"]["csy4"].to_numpy()))
            assert got == 8, f"{a}_{k}: {got} Csy4 runs, expected 8"
    panel_a(out, lay)
    panel_b(out, cells, layo)
    panel_c(out, cells, layo)
    panel_d(out, space, cells)
    report(space, targ, folds, lay, cells, layo)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
