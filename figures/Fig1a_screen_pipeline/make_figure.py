"""Fig1 panel A: the experimental screening pipeline, in-silico through validation.

Seven steps in one row: line glyphs in the thesis petrol, thin chevrons between
them, a sentence label under each. The phases are marked by a hairline rule above
the row with its name set into the gap, rather than by filled bands - a filled
band carries no information the rule does not, and reads as a slide.

Every size is written in POINTS and converted by U, which is set by the panel
width: the viewBox is 1000 units across 178 mm, so a unit is 0.178 mm and a point
is 1.982 units. Nothing here is a bare viewBox number. Before 2026-11-11 the sizes
were hand-picked viewBox units against W=1600, which put the labels at 4.3 pt -
nothing said so, because the SVG was only ever looked at zoomed in.

No panel letter is drawn: letters are added at assembly, so the artefact can be
placed as any panel. Hand-authored SVG like Fig2a_pipeline; text stays text for
Illustrator/Inkscape.
Glyphs come from code/common/glyphs.py, so they are shared with every other
figure and recoloured by the parent.

Writes to output/Fig1a_screen_pipeline/:
    Fig1a_screen_pipeline.svg   the panel, text kept as text
    .pdf / .png                      only if cairosvg is installed
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import yaml
from common import paths, palette, glyphs

out = paths.fig_output("Fig1a_screen_pipeline")

_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
DPI = int(_cfg["dpi"])
WIDTH_MM = float(_cfg["width_double_mm"])   # 178 mm, double column
LABEL_PT = float(_cfg["font_size_pt"])      # the house label size
PHASE_PT = LABEL_PT                         # phase names, caps and letterspaced. Was
                                            # LABEL_PT - 1: at a 5.5 pt house size that
                                            # is 4.5, under the 5 pt journal floor
HAIR_PT = float(_cfg["line_width_pt"])      # 0.75 pt hairlines

ACC = palette.MAIN              # the one accent; rule 1
INK = palette.INK
MUT = palette.NEUTRAL["600"]    # phase names
RULE = palette.NEUTRAL["300"]   # hairlines
CHEV = palette.NEUTRAL["400"]   # step-to-step chevrons
BG = palette.SURFACE
# Dry/wet grounds, chosen 2026-11-13. A tint per phase, neutral for the
# in-silico half and the theme tint for the experimental half, so the split is
# a change of colour rather than two words on a rule. Sanctioned by COLORS.md:
# rule 1 reserves theme 100-300 for "fills and backgrounds only", rule 3 the
# neutrals for reference bands. No stroke and no label of their own - a
# stroked pale box is the slide-deck tell STYLE_NOTES warns about, a tint
# behind a phase is not.
#
# MEASURED, and the reason to think before darkening either: neutral 100 is
# relative luminance 0.913 against theme 100's 0.863, a contrast ratio of
# 1.055. The two bands are ~5% apart in lightness, so the split is carried by
# hue and will NOT survive a greyscale print. Neutral 200 is 1.066, barely
# better and the other way round. For greyscale, go to neutral 300 or leave
# the dry half white.
BAND_FILLS = (palette.NEUTRAL["100"], palette.THEME["100"])
BAND_GAP = 12      # white between the bands, so the boundary is an edge on
                   # both sides and the chevron crossing it still reads

# (glyph, label lines). One row, seven steps.
STEPS = [
    ("binder-rna", ["RNA-protein binder", "generation"]),
    ("funnel", ["Selection of", "favourable designs"]),
    ("oligo-pool", ["Synthetic DNA", "pool ordering"]),
    ("integration", ["Bxb1 integration"]),
    ("tube", ["RNA extraction"]),
    ("long-reads", ["Long-read NGS", "and analysis"]),
    ("cytometry", ["Flow cytometry", "validation"]),
]

# Phases as (name, step indices). Edit this list to regroup; the rules follow it.
PHASES = [("IN SILICO", [0, 1]), ("EXPERIMENTAL", [2, 3, 4, 5, 6])]

W, H = 1000, 199        # 15 units taller than the unbanded panel: the bands
                        # need air under the labels or the type sits on an edge.
                        # 2026-11-16 (Leo): was 190, which left the SECOND label
                        # line closer to the band's foot than the first line is
                        # to its head - the two-line steps looked crowded while
                        # the one-line ones did not. Only the bands grow: the
                        # labels keep their LAB_Y, so this is air, not leading.
U = W * 0.352778 / WIDTH_MM     # viewBox units per point; 1.982 at 178 mm
X0, PAD = 30, 10
SLOT = (W - 2 * X0) / len(STEPS)
CENT = [X0 + SLOT / 2 + i * SLOT for i in range(len(STEPS))]

RULE_Y = 28        # phase hairline
ICON_CY = 92       # glyph centre
LAB_Y = 150        # first label line
LINE_H = LABEL_PT * 9.5 / 7 * U   # label leading: 9.5 pt on 7, kept in proportion
ICON_K = 1.34      # glyph scale; 48 pt box -> 64 units, ~11.5 mm printed

# --- 2026-11-16 (Leo): step 1 keeps its binder but gets the REAL hairpin.
# icon_binder_rna() is twisted_hairpin() + ribbon_helix(); only the hairpin
# half is swapped, for the same flat ChimeraX render Fig 2a uses, so the two
# pipeline panels show the same molecule. It is drawn INSIDE the glyph's own
# 48-unit box, over the coordinates the drawn hairpin occupied, so the binder
# does not move and the step still reads "a binder packed against the RNA".
# Much smaller than Fig 2a's: there the hairpin is the subject, here it is one
# step of seven. The render is flat two-tone petrol, not a lit cartoon, which
# is why it can sit among line glyphs at all (see Fig2a_pipeline/render_hairpin.py).
# Falls back to the drawn glyph if the render is missing, so run_all.sh never
# breaks on a clean checkout.
HP_PNG = (paths.OUTPUT / "Fig2a_pipeline" / "hairpin"
          / "hairpin_ladder_x-005_y000_r+14.8_s345.png")
HP_BOX = (3.0, 7.0, 23.0, 38.0)   # x, y, w, h in the 48-unit glyph box: where
                                  # twisted_hairpin(cx=15, ytop=19, ybot=42) sat
BIND_R = 5.0        # the binder helix radius, as icon_binder_rna() drew it
BIND_AT = 0.30      # where up the stem it docks, 0 = top of the ink, 1 = bottom
BIND_GAP = 1.4      # clear air between the render's ink and the ribbon


def _dock(im, x0, y0, f):
    """Where the binder sits: parallel to the stem, against its upper right.

    Returns (angle, cx, cy) in glyph units. The angle is the principal axis of
    the hairpin's own ink - a PCA over the alpha mask - so "parallel to the
    side" is read off the structure rather than eyeballed. The x is the ink's
    RIGHT EDGE at that height plus a gap, so the ribbon packs against the stem
    the way icon_binder_rna() always showed it, at whatever width the stem has
    there. Falls back to the glyph's own -38 deg if numpy is unavailable.
    """
    import math
    try:
        import numpy as np
    except ImportError:                      # pragma: no cover
        return -38, x0 + 30, y0 + 16
    a = np.array(im.getchannel("A"))
    ys, xs = np.nonzero(a > 40)
    pts = np.column_stack([xs, ys]).astype(float)
    axis = np.linalg.svd(pts - pts.mean(0), full_matrices=False)[2][0]
    if axis[1] < 0:                          # point it DOWN the page, so the
        axis = -axis                         # reported angle is the stem's
    ang = math.degrees(math.atan2(axis[1], axis[0]))
    # the ink's right edge in a band around BIND_AT, in image pixels
    row = ys.min() + BIND_AT * (ys.max() - ys.min())
    band = np.abs(ys - row) <= 0.09 * (ys.max() - ys.min())
    right = xs[band].max() if band.any() else xs.max()
    cx = x0 + right * f + BIND_R * 0.9 + BIND_GAP
    cy = y0 + row * f
    return ang, cx, cy


def binder_rna_real(cx, cy, k=ICON_K):
    """Step 1's glyph with the drawn hairpin replaced by the render."""
    if not HP_PNG.exists():
        return glyphs.icon_group("binder-rna", cx, cy, ACC, ICON_K)
    import base64
    import io
    from PIL import Image
    im = Image.open(HP_PNG)
    im = im.crop(im.getchannel("A").getbbox())      # trim to ink, as Fig 2a does
    buf = io.BytesIO()
    im.save(buf, "PNG")
    bx, by, bw, bh = HP_BOX
    f = min(bw / im.width, bh / im.height)          # fit, keeping aspect
    w, h = im.width * f, im.height * f
    x0, y0 = bx + (bw - w) / 2, by + (bh - h) / 2
    b64 = base64.b64encode(buf.getvalue()).decode()
    # The binder is docked ON the render, not at a hand-picked spot: its angle
    # and its place are MEASURED off the same alpha mask that is drawn, so a
    # re-render moves the binder with the hairpin instead of leaving it adrift.
    # Leo 2026-11-16: upper right, parallel to the stem.
    ang, bx_, by_ = _dock(im, x0, y0, f)
    inner = (f'<image x="{x0:.2f}" y="{y0:.2f}" width="{w:.2f}" '
             f'height="{h:.2f}" href="data:image/png;base64,{b64}"/>'
             + glyphs.ribbon_helix(bx_, by_, turns=2.0, R=BIND_R, pitch=10.5,
                                   W=6.0, ang=ang))
    half = 24 * k
    return (f'<g transform="translate({cx-half},{cy-half}) scale({k})" '
            f'fill="none" stroke="{ACC}" color="{ACC}" stroke-width="2.4" '
            f'stroke-linecap="round" stroke-linejoin="round">{inner}</g>')


def T(x, y, s, pt=LABEL_PT, col=INK, w=400, anch="middle", ls=0):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{pt*U:.2f}" '
            f'fill="{col}" font-weight="{w}" text-anchor="{anch}" '
            f'letter-spacing="{ls}">{s}</text>')


def hair(x1, x2, y=RULE_Y):
    return (f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y:.1f}" '
            f'stroke="{RULE}" stroke-width="{HAIR_PT*U:.2f}"/>')


def phase_rule(x1, x2, name, pt=PHASE_PT, ls=1.6):
    """Hairline across the phase with its name set into a gap in the middle."""
    sz = pt * U
    cx = (x1 + x2) / 2
    half = (len(name) * sz * 0.60 + len(name) * ls) / 2 + 9
    return (hair(x1, cx - half) + hair(cx + half, x2)
            + T(cx, RULE_Y + sz * 0.36, name, pt, MUT, 600, "middle", ls))


def chevron(x, y, s=7.0):
    return (f'<path d="M{x-s*0.55:.1f},{y-s:.1f} L{x+s*0.55:.1f},{y} '
            f'L{x-s*0.55:.1f},{y+s:.1f}" fill="none" stroke="{CHEV}" '
            f'stroke-width="{1.1*U:.2f}" stroke-linecap="round" '
            f'stroke-linejoin="round"/>')


def bands():
    """One ground per phase, cut back by BAND_GAP/2 where two meet. Driven by
    PHASES, so regrouping the phases regroups the bands with them."""
    assert len(PHASES) == len(BAND_FILLS), "a fill per phase, see PHASES"
    out, edges = [], []
    for i, (_, cols) in enumerate(PHASES):
        x1 = CENT[cols[0]] - SLOT / 2 - 4
        x2 = CENT[cols[-1]] + SLOT / 2 + 4
        if i:
            x1 = (CENT[PHASES[i-1][1][-1]] + CENT[cols[0]]) / 2 + BAND_GAP / 2
        if i < len(PHASES) - 1:
            x2 = (CENT[cols[-1]] + CENT[PHASES[i+1][1][0]]) / 2 - BAND_GAP / 2
        edges.append((x1, x2))
        out.append(f'<rect x="{x1:.1f}" y="{RULE_Y-15:.1f}" width="{x2-x1:.1f}" '
                   f'height="{H-RULE_Y+9:.1f}" rx="10" fill="{BAND_FILLS[i]}"/>')
    for (_, a), (b, _) in zip(edges, edges[1:]):
        assert b - a >= BAND_GAP - 0.01, "bands overlap"
    return out


def build():
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" '
         f'height="{WIDTH_MM*H/W:.2f}mm" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>']
    o += bands()

    for name, cols in PHASES:
        o.append(phase_rule(CENT[cols[0]] - SLOT / 2 + PAD,
                            CENT[cols[-1]] + SLOT / 2 - PAD, name))

    for i, (glyph, lines) in enumerate(STEPS):
        o.append(binder_rna_real(CENT[i], ICON_CY) if glyph == "binder-rna"
                 else glyphs.icon_group(glyph, CENT[i], ICON_CY, ACC, ICON_K))
        for j, ln in enumerate(lines):
            o.append(T(CENT[i], LAB_Y + j * LINE_H, ln))
        if i < len(STEPS) - 1:
            o.append(chevron((CENT[i] + CENT[i + 1]) / 2, ICON_CY))

    o.append('</svg>')
    return "\n".join(o)


def _check():
    """The labels must come out at the house size and stay inside their slot."""
    svg = build()
    assert f'font-size="{LABEL_PT*U:.2f}"' in svg, "label size is not font_size_pt"
    assert abs(LABEL_PT * U * WIDTH_MM / W / 0.352778 - LABEL_PT) < 1e-6
    widest = max(len(ln) for _, lines in STEPS for ln in lines)
    assert widest * LABEL_PT * U * 0.52 < SLOT, "a label line is wider than its slot"
    assert H / W < 6, "panel is a letterbox strip again"
    for fill in BAND_FILLS:
        assert svg.count(f'fill="{fill}"/>') == 1, f"band {fill} missing or doubled"
    # the last label line must sit inside its band, not on its edge
    last = LAB_Y + (max(len(l) for _, l in STEPS) - 1) * LINE_H
    assert last + LABEL_PT * U * 0.25 < H - 6, "a label overruns its band"
    print("checks pass: labels", LABEL_PT, "pt,", f"{WIDTH_MM:.0f} x {WIDTH_MM*H/W:.1f} mm")


if __name__ == "__main__":
    if "--check" in sys.argv:
        _check()
        sys.exit(0)
    svg_path = out / "Fig1a_screen_pipeline.svg"
    svg_path.write_text(build(), encoding="utf-8")
    glyphs.write_icons(out / "icons")
    pdf, png = (out / "Fig1a_screen_pipeline.pdf",
                out / "Fig1a_screen_pipeline.png")
    try:
        import cairosvg
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf))
        cairosvg.svg2png(url=str(svg_path), write_to=str(png),
                         output_width=int(W * DPI / 96))
    except ImportError:
        # Drop any raster left by an earlier run that had cairosvg: a stale
        # PDF/PNG beside a fresh SVG disagrees with it and nothing says so.
        for stale in (pdf, png):
            stale.unlink(missing_ok=True)
        print("SVG only (no cairosvg). The SVG is the deliverable; raster from "
              "the vector editor if a preview is wanted.")
    # Nothing is copied to final/ automatically; final/ is filled by hand.
    print("wrote", out)
