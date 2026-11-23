"""Fig 1a drafts, round 2: one horizontal row kept, glyphs and the dry/wet
separation reworked.

    python code/Fig1a_screen_pipeline/make_drafts2.py [outdir]
    default outdir: working/Fig1a_drafts/

Round 1 moved the pipeline onto two rows; one row was preferred, so the row
stays. What changes here is (a) the glyphs, which now draw the campaign's
molecules in the campaign's colours (see rich_glyphs.py), and (b) how the
in-silico half is told apart from the experimental half, which in the shipped
panel is only two words on one hairline.

Three separations, everything else identical, so the device is what is being
compared:

    gap       the two halves pushed apart by EXTRA units of air. Nothing drawn.
    ground    a pale NEUTRAL-100 band behind the in-silico pair, as the other
              panels use for their containers. Additive: the rules stay.
    divider   a vertical hairline on the boundary, running the panel's height.
    bands     BOTH halves on a ground: neutral behind in silico, theme tint
              behind experimental, so the split is a change of colour rather
              than an edge. Three of these, see BANDS.
    bands_soft   the same, both tints a step paler.
    bands_named  the same as bands, with the phase names set into the left of
              each band instead of into a rule, since the band already marks
              the phase and the rule then says it twice.

A variant whose name ends in `_old` is drawn with the SHIPPED one-colour line
glyphs from code/common/glyphs.py rather than the richer set in
rich_glyphs.py, so the separation can be judged without changing the icons.

COLORS.md sanctions both fills: rule 1 reserves theme 100-300 for "fills and
backgrounds only", rule 3 the neutrals for reference bands. The bands carry no
stroke and no label of their own - a stroked pale box is the slide-deck tell
STYLE_NOTES warns about, a tint behind a phase is not.

All three keep the house rules: no box round a step, chevrons not arrows,
sizes in points, no panel letter. The shipped panel and Fig 2a are untouched.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml                                                   # noqa: E402
import rich_glyphs                                            # noqa: E402
from common import paths, palette, glyphs                     # noqa: E402

REPO = paths.ROOT
_cfg = yaml.safe_load((REPO / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
WIDTH_MM = float(_cfg["width_double_mm"])
LABEL_PT = float(_cfg["font_size_pt"])
PHASE_PT = LABEL_PT
HAIR_PT = float(_cfg["line_width_pt"])

ACC = palette.MAIN         # the shipped glyphs are drawn in one accent
INK = palette.INK
MUT = palette.NEUTRAL["600"]
RULE = palette.NEUTRAL["300"]
CHEV = palette.NEUTRAL["400"]
GROUND = palette.NEUTRAL["100"]
BG = palette.SURFACE
# (dry fill, wet fill) per banded variant. Dry neutral, wet theme tint: the
# experimental half is where the campaign's blue lives, so the tint reads as
# "this is the part that happens in cells" rather than as decoration.
BANDS = {"bands": (palette.NEUTRAL["200"], palette.THEME["100"]),
         # COLORS.md: never hand-pick a hex in a figure script. The paler
         # pair is neutral 100 against the same theme 100, so only the grey
         # changes - which is the question this variant asks.
         "bands_soft": (palette.NEUTRAL["100"], palette.THEME["100"]),
         "bands_named": (palette.NEUTRAL["200"], palette.THEME["100"])}
BAND_GAP = 12          # air between the two bands, so the boundary is an edge
                       # on both sides and the chevron crossing it is clear

STEPS = [
    ("binder-rna", ["RNA-protein binder", "generation"]),
    ("funnel", ["Selection of", "favourable designs"]),
    ("oligo-pool", ["Synthetic DNA", "pool ordering"]),
    ("integration", ["Bxb1 integration"]),
    ("tube", ["RNA extraction"]),
    ("long-reads", ["Long-read NGS", "and analysis"]),
    ("cytometry", ["Flow cytometry", "validation"]),
]
PHASES = [("IN SILICO", [0, 1]), ("EXPERIMENTAL", [2, 3, 4, 5, 6])]
SPLIT = 2                      # first step of the second phase

W, H = 1000, 190
U = W * 0.352778 / WIDTH_MM
X0, PAD = 26, 10
EXTRA = 34                     # air on the boundary, 'gap' only. 34 units is
                               # 6 mm printed: a clear break without the halves
                               # drifting into two separate figures.
RULE_Y = 28
ICON_CY = 92
LAB_Y = 150
LINE_H = LABEL_PT * 9.5 / 7 * U
ICON_K = 1.34


def slots(extra):
    """Step centres, with `extra` units of air inserted on the boundary."""
    s = (W - 2 * X0 - extra) / len(STEPS)
    return [X0 + s / 2 + i * s + (extra if i >= SPLIT else 0)
            for i in range(len(STEPS))], s


def T(x, y, s, pt=LABEL_PT, col=INK, w=400, anch="middle", ls=0):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{pt*U:.2f}" '
            f'fill="{col}" font-weight="{w}" text-anchor="{anch}" '
            f'letter-spacing="{ls}">{s}</text>')


def hair(x1, x2, y=RULE_Y):
    return (f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y:.1f}" '
            f'stroke="{RULE}" stroke-width="{HAIR_PT*U:.2f}"/>')


def phase_rule(x1, x2, name):
    sz, ls = PHASE_PT * U, 1.6
    cx = (x1 + x2) / 2
    half = (len(name) * sz * 0.60 + len(name) * ls) / 2 + 9
    return (hair(x1, cx - half) + hair(cx + half, x2)
            + T(cx, RULE_Y + sz * 0.36, name, PHASE_PT, MUT, 600, "middle", ls))


def chevron(x, y, s=7.0):
    return (f'<path d="M{x-s*0.55:.1f},{y-s:.1f} L{x+s*0.55:.1f},{y:.1f} '
            f'L{x-s*0.55:.1f},{y+s:.1f}" fill="none" stroke="{CHEV}" '
            f'stroke-width="{1.1*U:.2f}" stroke-linecap="round" '
            f'stroke-linejoin="round"/>')


def band(x1, x2, fill):
    return (f'<rect x="{x1:.1f}" y="{RULE_Y-15:.1f}" width="{x2-x1:.1f}" '
            f'height="{H-RULE_Y+9:.1f}" rx="10" fill="{fill}"/>')


def build(variant):
    rich = not variant.endswith("_old")
    variant = variant[:-4] if not rich else variant
    extra = EXTRA if variant == "gap" else 0
    cent, slot = slots(extra)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" '
         f'height="{WIDTH_MM*H/W:.2f}mm" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>']

    if variant == "ground":
        # stop short of the chevron that leaves the phase: a band whose edge
        # cuts through the outgoing chevron reads as clipping it
        x1 = cent[0] - slot / 2 - 4
        x2 = (cent[SPLIT - 1] + cent[SPLIT]) / 2 - 16
        o.append(f'<rect x="{x1:.1f}" y="{RULE_Y-14:.1f}" width="{x2-x1:.1f}" '
                 f'height="{H-RULE_Y+6:.1f}" rx="10" fill="{GROUND}"/>')
    if variant in BANDS:
        dry, wet = BANDS[variant]
        xb = (cent[SPLIT - 1] + cent[SPLIT]) / 2
        o.append(band(cent[0] - slot / 2 - 4, xb - BAND_GAP / 2, dry))
        o.append(band(xb + BAND_GAP / 2, cent[-1] + slot / 2 + 4, wet))
    if variant == "divider":
        xb = (cent[SPLIT - 1] + cent[SPLIT]) / 2
        o.append(f'<line x1="{xb:.1f}" y1="{RULE_Y-12:.1f}" x2="{xb:.1f}" '
                 f'y2="{H-8:.1f}" stroke="{RULE}" '
                 f'stroke-width="{HAIR_PT*U:.2f}"/>')

    if variant == "bands_named":
        for name, cols in PHASES:
            o.append(T(cent[cols[0]] - slot / 2 + 4, RULE_Y + PHASE_PT * U * 0.36,
                       name, PHASE_PT, MUT, 600, "start", 1.6))
    else:
        for name, cols in PHASES:
            o.append(phase_rule(cent[cols[0]] - slot / 2 + PAD,
                                cent[cols[-1]] + slot / 2 - PAD, name))

    for i, (glyph, lines) in enumerate(STEPS):
        o.append(rich_glyphs.group(glyph, cent[i], ICON_CY, ICON_K) if rich
                 else glyphs.icon_group(glyph, cent[i], ICON_CY, ACC, ICON_K))
        for j, ln in enumerate(lines):
            o.append(T(cent[i], LAB_Y + j * LINE_H, ln))
        if i < len(STEPS) - 1 and not (variant == "divider" and i == SPLIT - 1):
            o.append(chevron((cent[i] + cent[i + 1]) / 2, ICON_CY))
    o.append("</svg>")
    return "\n".join(o)


def _check(variant, svg):
    variant = variant[:-4] if variant.endswith("_old") else variant
    _, slot = slots(EXTRA if variant == "gap" else 0)
    assert f'font-size="{LABEL_PT*U:.2f}"' in svg, "labels are not the house size"
    widest = max(len(ln) for _, lines in STEPS for ln in lines)
    assert widest * LABEL_PT * U * 0.52 < slot, "a label line is wider than its slot"
    # the divider replaces one chevron; every other gap between steps keeps one
    want = len(STEPS) - 1 - (1 if variant == "divider" else 0)
    assert svg.count('stroke-linejoin="round"/>') >= want
    if variant in BANDS:
        for fill in BANDS[variant]:
            assert svg.count(f'fill="{fill}"/>') >= 1, f"{variant}: {fill} missing"
    print(f"  {variant}: {WIDTH_MM:.0f} x {WIDTH_MM*H/W:.1f} mm, slot {slot*0.178:.1f} mm")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for variant in ("gap", "ground", "divider", "bands", "bands_soft",
                    "bands_named", "bands_soft_old", "bands_old"):
        svg = build(variant)
        _check(variant, svg)
        # 'rich' or 'plain' says which glyph set, so the two never collide
        stem = (f"Fig1a_plain_{variant[:-4]}" if variant.endswith("_old")
                else f"Fig1a_rich_{variant}")
        p = outdir / f"{stem}.svg"
        p.write_text(svg, encoding="utf-8")
        print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "working" / "Fig1a_drafts")
