"""Scaffold construction and selection: what was built, and what spliced.

Three bare panels, from `inputs/supp_scaffold_selection_ont.csv`
(`make_inputs.py`, which reads the ONT counter's own output, not the chart
pasted into the journal).

Revised 2026-11-19 (Leo, second round): "a is good but rename sequence checked
to something like signal strength measured / b is good but just put the ratios
in 0.xx format / c same for this also why is there a range for the ratios?
should only be one fixed number".

* **Colour is the architecture**, from the `scaffolds` block (palette v0.9,
  COLORS.md rule 9): ONE colour for all four, because they are four variants of
  the same kind of object. F30 is NOT one of them - pART158 keeps
  `palette.F30`, the blue the main figures' construct maps give the F30 array,
  which is what takes it out of the mods' hue family. Before today this panel
  gave F30 theme 900, i.e. the same family, which is what Leo asked to undo.
  The `stacks` block is reserved for design-route identity (rule 2) and an RNA
  scaffold is not a design route, so `palette.stack_color` is never called here.
* **Every mark is labelled.** The four architectures are the SAME colour, so
  colour groups but never identifies: each bar carries its pART number, each
  block its architecture name IN FULL - mod1 tandem compact, mod2 branched,
  mod3 tandem splayed, mod4 ring - and panel (a) keys all four over its top bar.
* **Ratios are printed as 0.xx**, with the raw counts under them, in (b) and
  (c). The axis carried the derived value and the bar carried only its inputs;
  now the reader gets the number without reading it off the axis.
* **The fold-change brackets print ONE number.** They used to read
  `1.10× (0.89-1.36)`, and the parenthesis - a 95% confidence interval - was
  read as a range of measured values. There is one measurement per construct
  per run, so there is no spread across replicates to show; the interval is an
  interval ESTIMATE of that single number. It therefore moves to the caption,
  labelled as a confidence interval, and the panel prints the point estimate.
  Nothing is lost: the binomial standard errors are still drawn on every bar
  and `report()` prints every interval on every run.
* **The sample inversion at LV1.17 is no longer drawn.** The bow-tie strip under
  panel (b) is gone; each bar now carries its corrected pART number and the
  correction is asserted in `_check_inversion` instead of illustrated. The
  caption carries the rest.
* **Counting error.** At ~300 evaluable reads the binomial standard error on a
  ratio near 0.4 is about 3 percentage points, so panel (c) draws it.

A splicing ratio is a reporter proxy for receptor dimerisation on a multivalent
scaffold, and so for avidity. The binder here is Csy4-H29A, which is
catalytically dead, so the readout cannot be cleavage of the scaffold: binding
the hairpin array recruits the MgIRE1 halves together and it is their
reconstitution that splices. It is not an affinity and there is no Kd anywhere
in this figure.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name

INK = palette.INK
MUTED = palette.NEUTRAL["600"]
NEG = palette.NEUTRAL["300"]              # negative controls: a reference, rule 3
FS = 5.5

# One colour per scaffold architecture, from the `scaffolds` block (COLORS.md
# rule 9, palette v0.9, 2026-11-20 - theme 600 `#0678a6`, a lighter step than
# the `MAIN` of v0.8). NOT the `stacks` block, which rule 2 reserves for
# design-route identity - these are RNA scaffolds, not design routes. ONE
# colour for all four, deliberately: they are four variants of the same kind of
# object, and four unrelated hues would say they are four unrelated things.
#
# F30 (pART158) is NOT one of them. It keeps `palette.F30`, the blue the main
# figures' construct maps give the F30 array, which is dE 69 from the four
# architectures - the point of the change. Until today this panel gave it
# theme 900, i.e. the mods' own hue family, which is what Leo asked to undo.
#
# The four architectures are the SAME colour, so COLORS.md rule 9 and
# STYLE_NOTES.md section 9 both require a LABEL ON EVERY MARK: every bar here
# carries its pART number, every block its architecture name in full, and
# panel (a) keys all four over its top bar.
ARCH_ORDER = ["mod1", "mod2", "mod3", "mod4"]
ARCH = dict(palette.SCAFFOLDS, F30=palette.F30)
assert [a for a in ARCH_ORDER] == palette.SCAFFOLD_ORDER, "scaffold order moved"

# What each architecture IS, spelled out wherever it is labelled (Leo,
# 2026-11-19), RENAMED 2026-11-20. mod1 and mod3 are the same tandem array
# differing in how far apart their modules sit, which is exactly the
# distinction a bare "mod1 / mod3" hides. They used to be named "F30, long
# arms" / "F30, short arms"; in the designs that were actually ORDERED every
# arm is 5 bp, so that difference no longer exists and the old names inverted
# (mod3 draws larger than mod1). See COLORS.md rule 9 and
# supp/scaffold_design/README.md.
ARCH_NAME = {"mod1": "tandem, compact", "mod2": "branched",
             "mod3": "tandem, splayed", "mod4": "ring",
             "F30": "F30 array"}

# The first run's ten samples, in the order the panel draws them: by
# architecture, then by copy number. pART158 is the pre-existing tenth.
RUN1_ORDER = ["pART243", "pART244", "pART245",
              "pART246", "pART247", "pART248",
              "pART252", "pART253", "pART254", "pART158"]
RUN1_GROUP = {"pART243": "mod1", "pART244": "mod1", "pART245": "mod1",
              "pART246": "mod2", "pART247": "mod2", "pART248": "mod2",
              "pART252": "mod4", "pART253": "mod4", "pART254": "mod4",
              "pART158": "F30"}
RUN1_UNITS = {"pART243": 1, "pART244": 2, "pART245": 3,
              "pART246": 1, "pART247": 2, "pART248": 3,
              "pART252": 1, "pART253": 2, "pART254": 3,
              "pART158": 8}
# Csy4 hairpins presented by one unit of each scaffold: a designed mod unit
# carries 8 sites, an F30 repeat 2. The panel DERIVES its hairpin count from
# these and the copy number rather than carrying a typed list of twelve
# numbers, and `_check_sites()` asserts the result against the counts taken
# off the plasmid records themselves.
SITES_PER_UNIT = {"mod1": 8, "mod2": 8, "mod3": 8, "mod4": 8, "F30": 2}
GAP = 0.6          # extra x spacing between architecture blocks


def copy_label(c):
    """What the first tick line says: copies for a mod, units for the array."""
    return (f"{RUN1_UNITS[c]} units" if RUN1_GROUP[c] == "F30"
            else f"{RUN1_UNITS[c]}×")


def csy4_sites(c):
    """Csy4 hairpins the construct presents (Leo, 2026-11-20: "add to b how
    many csy4 hairpins each column includes"). Derived, never typed."""
    return RUN1_UNITS[c] * SITES_PER_UNIT[RUN1_GROUP[c]]

# Panel (a). Each stage is a count over the SAME 12 designed combinations -
# four architectures at one, two and three copies - and every one of them is
# traceable to a file, not to a summary. The per-architecture split is read off
# the 12 SnapGene map names in LV_Cloning/Projects/F30_mod/Assembled_Final:
#   designed            12 maps, three per architecture (243-245 mod1,
#                           246-248 mod2, 249-251 mod3, 252-254 mod4)
#   fragment delivered  11; pART249 is the map marked `_notcloned`, so mod3
#                           loses its one-copy member
#   assembled           9;  the 9 submitted to LV1.16 Sanger sequencing - a
#                           plasmid has to exist to be sequenced - and the whole
#                           mod3 series (pART250, pART251) is absent
#   splicing measured   9;  the 9 designed constructs in the LV1.15 ONT run,
#                           which is panel (b). Renamed from "Sequence checked"
#                           2026-11-19 (Leo: "rename sequence checked to
#                           something like signal strength measured"). CHECKED
#                           FIRST, and the two records are NOT the same thing:
#                           LV1.16_SangerSeq holds a .dna per assembled
#                           construct and the LV1.15 ONT run holds read counts.
#                           They cover the IDENTICAL nine (243-248, 252-254,
#                           plus pART158 in both), so drawing one bar loses no
#                           construct - but the stage now counts the readout
#                           this figure is built on, and the caption says the
#                           same nine were Sanger-verified.
#   carried forward     1;  pART252 (mod4), from the LV1.27 repeat readout
STAGES = [("Designed", dict(mod1=3, mod2=3, mod3=3, mod4=3)),
          ("Fragment delivered", dict(mod1=3, mod2=3, mod3=2, mod4=3)),
          ("Assembled", dict(mod1=3, mod2=3, mod3=0, mod4=3)),
          ("Splicing measured", dict(mod1=3, mod2=3, mod3=0, mod4=3)),
          ("Carried forward", dict(mod1=0, mod2=0, mod3=0, mod4=1))]
N_DESIGNED = 12


def load():
    d = pd.read_csv(paths.INPUTS / "supp_scaffold_selection_ont.csv")
    return d[d.run == "LV1.15"].set_index("construct"), d[d.run == "LV1.27.1"]


# --- statistics, all hand-rolled: scipy is broken in this venv (_propack) ----
def se(k, n):
    """Binomial standard error sqrt(p(1-p)/n)."""
    p = k / n
    return math.sqrt(p * (1 - p) / n)


def rule_of_three(n):
    """One-sided 95% upper bound on a rate when the count is zero."""
    return 3.0 / n


def fc_ci(k1, n1, k2, n2, z=1.96):
    """Fold change p1/p2 with a two-sided 95% CI.

    Delta method on log(p1/p2): var = (1-p1)/(n1 p1) + (1-p2)/(n2 p2), which is
    the usual risk-ratio interval. z = 1.96 is the normal quantile, written out
    because there is no `scipy.stats.norm` to ask. Both counts must be > 0; a
    zero denominator makes the ratio undefined and is handled with
    `rule_of_three` instead, never printed as an infinity.
    """
    assert k1 > 0 and k2 > 0, "fold change is undefined with a zero count"
    p1, p2 = k1 / n1, k2 / n2
    fc = p1 / p2
    s = math.sqrt((1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2))
    return fc, fc * math.exp(-z * s), fc * math.exp(z * s)


def _selfcheck():
    """Smallest check that fails if the ratio statistics break."""
    fc, lo, hi = fc_ci(100, 300, 100, 300)      # identical rates
    assert abs(fc - 1.0) < 1e-12 and lo < 1 < hi, (fc, lo, hi)
    fc, lo, hi = fc_ci(200, 300, 100, 300)      # a real doubling
    assert abs(fc - 2.0) < 1e-12 and lo > 1, (fc, lo, hi)
    assert abs(se(150, 300) - 0.5 / math.sqrt(300)) < 1e-12
    assert rule_of_three(300) == 0.01


def _check_sites():
    """The Csy4 hairpin count on every (b) column, against the plasmid records.

    The panel derives it (`csy4_sites`: units x sites per unit). These two
    tables COUNTED it, by searching the 28 nt Csy4 extended repeat in the
    SnapGene maps and the GenBank files - `supp_scaffold_design_space.csv` for
    the nine designed constructs, `supp_negative_controls_design.csv` for
    pART158, which is not a mod and so is in neither mod table. If a construct
    is ever relabelled or the per-unit site count changes, the number drawn on
    the panel and the number counted off the plasmid stop agreeing here rather
    than in the thesis.
    """
    space = pd.read_csv(paths.INPUTS / "supp_scaffold_design_space.csv")
    counted = dict(zip(space.construct, space.n_csy4))
    neg = pd.read_csv(paths.INPUTS / "supp_negative_controls_design.csv")
    counted.update(zip(neg.positive, neg.csy4_in_positive))
    for c in RUN1_ORDER:
        assert csy4_sites(c) == int(counted[c]), (
            f"{c}: panel draws {csy4_sites(c)} Csy4 hairpins, the plasmid "
            f"record counts {counted[c]}")


def _check_inversion(run1):
    """The LV1.17 correction, asserted rather than drawn.

    Panel (b) used to carry a bow-tie strip showing sample k -> construct
    11 - k. Leo asked for the lines out and the right pART number on each bar
    instead, so the only thing left guarding the labels is this: in the
    corrected assignment the ten samples run the design order BACKWARDS, so
    `RUN1_ORDER[i]` must be sample 10 - i. If the input CSV ever reverted to
    the uncorrected file, every bar would be mislabelled silently and this is
    what catches it.
    """
    for i, c in enumerate(RUN1_ORDER):
        k = int(run1.loc[c, "sample"].split(".")[-1])
        assert k == 10 - i, f"{c} is sample {k}, expected {10 - i}"


def _positions(groups):
    """x for each bar, with GAP inserted between architecture blocks."""
    xs, x = [], 0.0
    for i, g in enumerate(groups):
        if i and g != groups[i - 1]:
            x += GAP
        xs.append(x)
        x += 1.0
    return np.array(xs)


# --- (a) what survived from design to the screen ---------------------------
def panel_a(out):
    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(38)))
    y = np.arange(len(STAGES))[::-1].astype(float)
    for yy, (_, split) in zip(y, STAGES):
        left = 0.0
        for a in ARCH_ORDER:
            w = split[a]
            if w:
                # SURFACE, not INK: the four architectures are ONE colour, so
                # a near-black hairline between two
                # adjacent segments is invisible and the split - the whole
                # point of segmenting (a) - disappears. A white divider is
                # what makes it read. Panels (b) and (c) keep the INK edge:
                # their bars stand alone on white, and the pale F30 fill
                # needs it (STYLE_NOTES section 9).
                ax.barh(yy, w, 0.62, left=left, color=ARCH[a],
                        edgecolor=palette.SURFACE, linewidth=0.5)
                # EVERY segment is named, not just the top row's. As soon as
                # one architecture goes to zero the ones after it slide left,
                # so the key above the top bar stops applying: in the two
                # 9/12 rows the mod4 segment sits under the "tandem, splayed"
                # key. With all four the same colour a reader cannot recover
                # the mapping at all, so a key that stops applying to three of
                # five rows is worse than repeating the label. 4.0 pt is what
                # the narrowest segment takes - MEASURED with style.text_w,
                # "mod4" is 3.37 mm at 4.0 pt against the 4.43 mm one-unit
                # segment in "Carried forward" (3.79 mm at 4.5). The full
                # names stay on the key above the top row, where each sits
                # directly over its own `mod` id, which is what maps them.
                ax.annotate(a, (left + w / 2, yy), ha="center", va="center",
                            fontsize=FS - 1.5, color=palette.SURFACE)
            left += w
        ax.annotate(f"{left:.0f}/{N_DESIGNED}", (left, yy),
                    textcoords="offset points", xytext=(3, 0), va="center",
                    ha="left", fontsize=FS)
    # The architecture key: the full names, direct over the top row's
    # segments (no legend box), where each sits over its own `mod` id and so
    # names it. Every segment in every row carries that id, so identity never
    # rests on colour - which it could not, the four being one colour.
    # 5.0 pt because that is what fits - MEASURED with style.text_w, the
    # widest name is 12.1 mm at 5.0 pt against a 13.3 mm segment, and 13.4 mm
    # at the house 5.5.
    for i, a in enumerate(ARCH_ORDER):
        ax.annotate(ARCH_NAME[a], (3 * i + 1.5, y[0] + 0.33),
                    textcoords="offset points", xytext=(0, 1), ha="center",
                    va="bottom", fontsize=FS - 0.5, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([s for s, _ in STAGES])
    ax.set_xlabel("Designed combinations")
    ax.set_xlim(0, N_DESIGNED * 1.2)
    ax.set_xticks([0, 4, 8, 12])
    ax.set_ylim(-0.6, len(STAGES) - 0.15)
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.annotate("plus pART158, which already existed",
                (1, 0), textcoords="offset points", xytext=(16, 0),
                va="center", ha="left", fontsize=FS, color=MUTED)
    fig.subplots_adjust(left=0.24, right=0.99, top=0.93, bottom=0.20)
    style.save(fig, out, "panel_a_assembly_outcome")
    plt.close(fig)


# --- (b) first readout, constructs labelled with the CORRECTED pART ---------
def panel_b(run1, out):
    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(44)))
    d = run1.reindex(RUN1_ORDER)
    groups = [RUN1_GROUP[c] for c in RUN1_ORDER]
    x = _positions(groups)
    # A light band behind each architecture's bars (Leo, 2026-11-20: "you could
    # do a light background behind three bars of the same mod in b to bunch
    # them together more visually"). NEUTRAL 100, a tint and not a hue: the
    # band groups, it does not encode, and COLORS.md rule 3 keeps anything that
    # is a reference rather than a series neutral - a fifth colour here would
    # read as a fifth data series. zorder 0 with the bars at 2 puts it under
    # the bars, and matplotlib draws the spines at 2.5, so neither the axis nor
    # the bars' INK edges are tinted by it. The band stops 0.5 either side of
    # the end bars, leaving the full GAP of white between blocks.
    for g in dict.fromkeys(groups):
        xs = x[[i for i, gg in enumerate(groups) if gg == g]]
        ax.axvspan(xs[0] - 0.5, xs[-1] + 0.5, facecolor=palette.NEUTRAL["100"],
                   edgecolor="none", zorder=0)
    ax.bar(x, d.splicing_ratio, 0.84, color=[ARCH[g] for g in groups],
           edgecolor=INK, linewidth=0.3, zorder=2)
    # Leo 2026-11-19: "just put the ratios in 0.xx format". The derived value
    # on the upper line at the house size, the raw counts it came from under it
    # (STYLE_NOTES section 5 wants both). The counts are set at 4.0 pt because
    # that is what FITS: measured with style.text_w, "408/1009" is 5.80 mm at
    # 4.0 pt against a 6.03 mm bar pitch, and 6.53 mm at 4.5 - which is why the
    # old single-line labels of two neighbouring bars nearly touched.
    for xx, (_, r) in zip(x, d.iterrows()):
        ax.annotate(f"{r.removed}/{r.evaluable}", (xx, r.splicing_ratio),
                    textcoords="offset points", xytext=(0, 1.8), ha="center",
                    fontsize=FS - 1.5, color=MUTED)
        ax.annotate(f"{r.splicing_ratio:.2f}", (xx, r.splicing_ratio),
                    textcoords="offset points", xytext=(0, 7.0), ha="center",
                    fontsize=FS)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{copy_label(c)}\n{c.replace('pART', '')}"
                        for c in RUN1_ORDER])
    # A third tick row: how many Csy4 hairpins the column's construct presents.
    # Drawn as annotations rather than a third line in the tick label so the
    # row can carry its own name in the gutter - a bare 8 / 16 / 24 under a bar
    # would be read as anything. The number is DERIVED (copy number x sites per
    # unit) and asserted against the plasmid records; see `_check_sites`.
    for xx, c in zip(x, RUN1_ORDER):
        ax.annotate(f"{csy4_sites(c)}", (xx, 0),
                    xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(0, -19), ha="center",
                    va="top", fontsize=FS, color=MUTED, annotation_clip=False)
    # The row's name, in the left gutter, on TWO lines. MEASURED: the gutter
    # between the left spine and the rotated y label is 7.7 mm, and
    # "Csy4 hairpins" on one line is 11.6 mm at the house 5.5 pt - setting it
    # on one line pushed the tight bounding box out to 262 pt, so this panel
    # would have been scaled to the 85 mm column 9% harder than its siblings
    # (245-248 pt) and its type would have come out visibly smaller than
    # theirs. Broken in two the widest line is "hairpins", 6.6 mm, and the
    # panel stays at its sibling width.
    for dy, lab in ((-19, "Csy4"), (-24.5, "hairpins")):
        ax.annotate(lab, (0, 0), xycoords="axes fraction",
                    textcoords="offset points", xytext=(-3, dy), ha="right",
                    va="top", fontsize=FS, color=MUTED, annotation_clip=False)
    ax.set_ylabel("Splicing ratio")
    ax.set_ylim(0, 0.62)          # headroom for the two-line bar labels
    ax.set_yticks([0, 0.2, 0.4])
    ax.spines["left"].set_bounds(0, 0.4)     # as in (c): the axis stops at the
    ax.set_xlim(x[0] - 0.8, x[-1] + 0.8)     # last tick, not in the label space
    ax.tick_params(axis="both", length=2, width=0.5)

    # architecture names under their own block, direct, no legend, in full
    for g in dict.fromkeys(groups):
        xs = x[[i for i, gg in enumerate(groups) if gg == g]]
        ax.annotate(ARCH_NAME[g], ((xs[0] + xs[-1]) / 2, 0),
                    xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(0, -31), ha="center",
                    va="top", fontsize=FS, color=MUTED, annotation_clip=False)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.95, bottom=0.38)
    style.save(fig, out, "panel_b_first_readout")
    plt.close(fig)


def _bracket(ax, x0, x1, y, text, sub=None, drop=0.022, inset=0.10):
    """A comparison bracket: one value, and optionally the basis for it.

    `sub` used to carry the fold change's 95% CI as `(0.89-1.36)`, and Leo read
    that parenthesis as a RANGE OF VALUES rather than an interval estimate -
    "why is there a range for the ratios? should only be one fixed number".
    He is right that it misreads, so the bracket now prints the point estimate
    alone and the confidence intervals are quoted in the caption, labelled as
    confidence intervals, where a word of prose can say what they are. The
    uncertainty is not dropped: the binomial standard errors are still drawn on
    every bar, which is the same information in the form that cannot be
    mistaken for a range of measurements, and `report()` prints every CI on
    every run.

    `sub` survives for the control bracket, whose "vs 95% control bound" is a
    BASIS and not a range. `inset` pulls both ends in so two brackets sharing a
    bar do not run into one line.
    """
    x0, x1 = x0 + inset, x1 - inset
    ax.plot([x0, x0, x1, x1], [y - drop, y, y, y - drop], color=MUTED,
            lw=0.5, solid_capstyle="butt", clip_on=False)
    mid = (x0 + x1) / 2
    if sub:
        ax.annotate(sub, (mid, y), textcoords="offset points", xytext=(0, 1.5),
                    ha="center", va="bottom", fontsize=FS - 1.0, color=MUTED)
    ax.annotate(text, (mid, y), textcoords="offset points",
                xytext=(0, 8.0 if sub else 1.5), ha="center", va="bottom",
                fontsize=FS)


# --- (c) the repeat, which is what the selection rests on ------------------
def panel_c(run2, out):
    order = ["pART158", "pART252", "pART243", "pART215", "pART242"]
    arch = ["F30", "mod4", "mod1", None, None]
    d = run2.set_index("construct").reindex(order)
    x = _positions(["t", "t", "t", "n", "n"])
    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(52)))
    ax.bar(x, d.splicing_ratio, 0.72,
           color=[ARCH[a] if a else NEG for a in arch],
           edgecolor=INK, linewidth=0.3)
    err = [se(r.removed, r.evaluable) if r.removed else 0.0
           for _, r in d.iterrows()]
    ax.errorbar(x, d.splicing_ratio, yerr=err, fmt="none", ecolor=INK,
                elinewidth=0.5, capsize=1.2, capthick=0.5)
    # As in (b): the ratio in 0.xx on the upper line, its raw counts under it.
    # There is room here - 12.1 mm of bar pitch against a 5.4 mm "117/300" - so
    # the counts keep 4.5 pt rather than (b)'s 4.0. A zero bar gets its
    # one-sided bound on the ratio line instead of a "0.00" the bar already
    # says, in the same 0.xxx format rather than the old percent, so every
    # number in the column reads on one scale.
    for xx, (c, r) in zip(x, d.iterrows()):
        top = r.splicing_ratio + (se(r.removed, r.evaluable) if r.removed else 0)
        ax.annotate(f"{r.removed}/{r.evaluable}", (xx, top),
                    textcoords="offset points", xytext=(0, 2.0), ha="center",
                    fontsize=FS - 1.0, color=MUTED)
        ax.annotate(f"{r.splicing_ratio:.2f}" if r.removed
                    else f"< {rule_of_three(r.evaluable):.3f}", (xx, top),
                    textcoords="offset points", xytext=(0, 7.4), ha="center",
                    fontsize=FS, color=INK if r.removed else MUTED)

    # --- the fold changes between conditions, which is what (c) is for -----
    # Between target scaffolds: the point estimate, one number per bracket.
    # The 95% CI is still computed - it is asserted on below and printed by
    # report() - but it is QUOTED IN THE CAPTION rather than drawn, because
    # "(0.89-1.36)" beside a ratio reads as a range of measured values and was
    # read that way. Both intervals contain 1, so the panel does not rank the
    # target scaffolds; the overlapping error bars are what says so on the
    # canvas, and the caption says it in words.
    for i, j in ((0, 1), (1, 2)):
        a, b = d.iloc[i], d.iloc[j]
        fc, lo, hi = fc_ci(a.removed, a.evaluable, b.removed, b.evaluable)
        assert lo < 1 < hi, (
            f"{d.index[i]}/{d.index[j]} CI {lo:.2f}-{hi:.2f} no longer "
            "contains 1 - the caption claims it does, so fix the caption")
        _bracket(ax, x[i], x[j], 0.55, f"{fc:.2f}×")
    # Against the negative controls the denominator is 0/329 and 0/267, so the
    # ratio is undefined and an infinity would be a lie. The honest statement
    # is a LOWER bound: divide by the weaker of the two one-sided 95% bounds
    # (the larger number, 3/267) and take the weakest target, so the figure
    # printed holds for every target against both controls.
    floor = max(rule_of_three(r.evaluable) for _, r in d.iterrows() if not r.removed)
    lower = min(r.splicing_ratio / floor for _, r in d.iterrows() if r.removed)
    _bracket(ax, x[0], x[4], 0.70, f"> {math.floor(lower):.0f}×",
             "vs 95% control bound")

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("pART", "") for c in order])
    ax.set_ylabel("Splicing ratio")
    ax.set_ylim(0, 0.80)
    ax.set_yticks([0, 0.2, 0.4])
    ax.spines["left"].set_bounds(0, 0.4)     # the axis stops where the data do
    ax.set_xlim(x[0] - 0.8, x[-1] + 0.8)
    ax.tick_params(axis="both", length=2, width=0.5)
    for xs, lab in (((0, 2), "target arrays"), ((3, 4), "negative controls")):
        ax.annotate(lab, ((x[xs[0]] + x[xs[1]]) / 2, 0),
                    xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(0, -13), ha="center",
                    va="top", fontsize=FS, color=MUTED, annotation_clip=False)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.97, bottom=0.21)
    style.save(fig, out, "panel_c_repeat_readout")
    plt.close(fig)


def report(run1, run2):
    print("--- LV1.15, corrected construct assignment ---")
    for c in RUN1_ORDER:
        r = run1.loc[c]
        print(f"  {c}  {RUN1_GROUP[c]:5s} {copy_label(c):8s} "
              f"{csy4_sites(c):2d} Csy4  "
              f"{r.removed:4d}/{r.evaluable:5d} = {r.splicing_ratio:.3f} "
              f"({r.splicing_ratio:.1%})  (sample {r['sample']})")
    lo, hi = run1.splicing_ratio.min(), run1.splicing_ratio.max()
    print(f"  span {lo:.3f} to {hi:.3f} "
          f"({run1.splicing_ratio.idxmin()} to {run1.splicing_ratio.idxmax()})")

    print("--- LV1.27.1 repeat ---")
    d = run2.set_index("construct")
    for c in ["pART158", "pART252", "pART243", "pART215", "pART242"]:
        r = d.loc[c]
        pm = se(r.removed, r.evaluable) if r.removed else rule_of_three(r.evaluable)
        kind = "+/-" if r.removed else "< (rule of three)"
        print(f"  {c}  {r.removed:4d}/{r.evaluable:4d} = {r.splicing_ratio:.3f} "
              f"({r.splicing_ratio:.1%})  {kind} {pm:.3f}   {r.note}")

    print("--- fold changes drawn on panel (c) ---")
    for a, b in (("pART158", "pART252"), ("pART252", "pART243"),
                 ("pART158", "pART243")):
        ra, rb = d.loc[a], d.loc[b]
        fc, lo, hi = fc_ci(ra.removed, ra.evaluable, rb.removed, rb.evaluable)
        flag = "" if lo < 1 < hi else "   <- 1 outside the CI"
        print(f"  {a} / {b}  {fc:.2f}x  95% CI {lo:.2f}-{hi:.2f}{flag}")
    floor = max(rule_of_three(r.evaluable) for _, r in d.iterrows() if not r.removed)
    print(f"  control 95% bounds: " + ", ".join(
        f"{c} < {rule_of_three(r.evaluable):.4f}"
        for c, r in d.iterrows() if not r.removed))
    for c, r in d.iterrows():
        if r.removed:
            print(f"  {c} / weaker control bound ({floor:.4f}) "
                  f"> {r.splicing_ratio / floor:.1f}x")

    print("--- between the two runs (NOT drawn; see README) ---")
    for c in ["pART158", "pART252", "pART243"]:
        print(f"  {c}  run1 {run1.loc[c].splicing_ratio:.3f} -> "
              f"run2 {d.loc[c].splicing_ratio:.3f}  "
              f"({d.loc[c].splicing_ratio / run1.loc[c].splicing_ratio:.2f}x)")


def main():
    _selfcheck()
    style.apply()
    out = paths.fig_output(FIG)
    run1, run2 = load()
    assert len(run1) == 10 and len(run2) == 5
    assert sum(sum(s.values()) for _, s in STAGES[:1]) == N_DESIGNED
    _check_inversion(run1)
    _check_sites()
    panel_a(out)
    panel_b(run1, out)
    panel_c(run2, out)
    report(run1, run2)
    print("wrote", out)


if __name__ == "__main__":
    main()
