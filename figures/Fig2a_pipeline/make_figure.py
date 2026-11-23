"""Fig02: the design pipeline, from RNA target to ordered oligo library.

No title and no standfirst in the artefact: both belong to the caption, which
is drafted in this figure's README (STYLE_NOTES.md, checklist item 1).

Schematic, no input data. The three counts are the campaign figures recorded in
Methods 2.9 (generation frozen at 6,992,515; the five-criterion gate leaves
55,706; the clustered draw gives the 12,000 ordered panel).

Writes to output/Fig2a_pipeline/:
    Fig2a_pipeline.svg   assembled figure, text kept as text
    Fig2a_pipeline.pdf   same, for LaTeX (needs cairosvg)
    Fig2a_pipeline.png   300 dpi raster (needs cairosvg)
    icons/*.svg             every glyph on its own, reusable

PDF and PNG are skipped with a note if cairosvg is unavailable; the SVG is
always written.
"""
import base64
import io
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import random
import yaml
from common import paths, palette, glyphs

out = paths.fig_output("Fig2a_pipeline")

# Shared font comes from config/style.yaml. The semantic palette below is local:
# config/style.yaml still ships the placeholder matplotlib colours, so adopting
# them here would drop the diffusion / hallucination colour coding. Move these
# into config/style.yaml once the thesis palette is settled.
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
DPI = int(_cfg["dpi"])

# diffusion / hallucination / sequence-design are PROCESS classes, not stacks,
# so they take theme steps rather than stack colours (rule 2). Rule 6 leaves
# schematic elements to be derived from theme and neutrals; this is that.
TEAL = palette.THEME["700"]   # diffusion
AMB = palette.THEME["400"]    # hallucination
IND = palette.THEME["800"]    # sequence design, selection, library
GREY = palette.NEUTRAL["500"]
LGREY = palette.NEUTRAL["300"]
INK = palette.INK
MUT = palette.NEUTRAL["600"]
BG = palette.SURFACE

W, H = 1520, 430        # was 404: 2026-11-13, air for the ground's top label
YM, YT, YB = 206, 102, 310   # was 190, 86, 294: everything down 16


def T(x, y, s, sz=13, col=INK, w=400, anch="middle", ls=0):
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{sz}" '
            f'fill="{col}" font-weight="{w}" text-anchor="{anch}" '
            f'letter-spacing="{ls}">{s}</text>')


def harrow(x1, x2, y, col=GREY, sw=2.6, head=8):
    return (f'<line x1="{x1}" y1="{y}" x2="{x2-head}" y2="{y}" stroke="{col}" '
            f'stroke-width="{sw}"/>'
            f'<path d="M{x2},{y} L{x2-head},{y-head*0.72} '
            f'L{x2-head},{y+head*0.72} Z" fill="{col}"/>')


DOT_R = 1.8   # one dot size for all three clouds: only the count shrinks


def cloud(cx, cy, rx, ry, n, col, seed, r=DOT_R):
    """n dots in an ellipse, rejection-sampled so no two touch."""
    rnd = random.Random(seed)
    pts = []
    min_d2 = (2.4 * r) ** 2
    for _ in range(n * 200):
        if len(pts) == n:
            break
        u, v = rnd.uniform(-1, 1), rnd.uniform(-1, 1)
        if u * u + v * v > 1:
            continue
        x, y = cx + u * rx, cy + v * ry
        if all((x - a) ** 2 + (y - b) ** 2 >= min_d2 for a, b in pts):
            pts.append((x, y))
    return "".join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}"/>'
                   for x, y in pts)


# --- 2026-11-13 (Leo): the whole panel sits on the in-silico ground, and the
# target hairpin at the far left gets a larger slot ahead of a ChimeraX render.
GROUND = palette.NEUTRAL["100"]     # Fig 1a's IN SILICO band: this is all dry lab
HP_K = 4.0                          # hairpin glyph scale (was the default 2.3)
HP_HALF = 24 * HP_K                 # half its 48-unit box: 96
HP_CX = 26 + HP_HALF                # hairpin centre, 26 in from the left edge
# THE HAIRPIN SLOT, for the render that will replace the glyph: a square of
# side 2*HP_HALF centred on (HP_CX, YM), label below it. Keep a render inside it.
SLOT = (HP_CX - HP_HALF, YM - HP_HALF, 2 * HP_HALF, 2 * HP_HALF)
# Everything after the hairpin keeps its size and is only moved: x positions
# from the first arrow on are mapped onto the room that is left. Positions
# only - no scale group - so text keeps one size for the page assembler.
OLD0, NEW0 = 152, HP_CX + HP_HALF + 14
SQ = (W - 22 - NEW0) / (1473 - OLD0)  # 1473: right edge of the last label
# 2026-11-14 (Leo): the slot now holds the real hairpin, rendered flat by
# code/Fig2a_pipeline/render_hairpin.py so it reads as a glyph among glyphs
# rather than as a photograph dropped into a diagram. Rebuild the PNG there,
# not here. The glyph stays as the fallback so run_all.sh never breaks on a
# missing render.
HP_PNG = out / "hairpin" / "hairpin_ladder_x-005_y000_r+14.8_s345.png"
# 2026-11-14 (Leo): a bracket over the three flipped loop bases, in Fig 1d's
# bracket style. WHICH THREE IS MEASURED, not assumed: ring centroids in the
# target put the continuous stack at nt2-7 and nt11-16, with breaks at 7-8
# (8.2 A), 9-10 (8.6 A) and 10-11 (10.4 A), so U8, A9 and U10 are the bases
# that leave the helix. WHERE they land is measured too - render_hairpin.py
# --probe 8-10 renders those three alone on the same camera, and the bracket
# spans the alpha bbox of that probe mapped through this panel's own crop and
# scale. Re-render the hairpin and the bracket follows it; delete the probe
# and the bracket simply does not draw.
HP_PROBE = out / "hairpin" / "hairpin_probe.png"
FLIP_LABEL = "flipped bases"
FLIP_PT = 11                    # the panel's secondary label size
FLIP_GAP = 7.0                  # base ink to the bracket's arm tips
# Fig 1d draws its module brackets with a 2.6-unit arm and a 0.35 stroke
# against a 5.5 px label. This panel's labels are twice that, so the bracket
# is scaled by the same ratio and reads as the same bracket.
FLIP_K = FLIP_PT / 5.5


def X(x):
    return NEW0 + (x - OLD0) * SQ


def hairpin():
    """The ChimeraX render, trimmed to its ink and fitted in the slot.

    Embedded base64 rather than linked, so the SVG stays a single file that
    survives being moved into the assembly folder or opened in Illustrator.
    """
    if not HP_PNG.exists():
        return glyphs.icon_group("rna-hairpin", HP_CX, YM, TEAL, HP_K)
    from PIL import Image
    im = Image.open(HP_PNG)
    bb = im.getchannel("A").getbbox()
    im = im.crop(bb)
    buf = io.BytesIO()
    im.save(buf, "PNG")
    k = min(SLOT[2] / im.width, SLOT[3] / im.height)
    w, h = im.width * k, im.height * k
    x0, y0 = HP_CX - w / 2, YM - h / 2
    b64 = base64.b64encode(buf.getvalue()).decode()
    o = [f'<image x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" '
         f'height="{h:.1f}" href="data:image/png;base64,{b64}"/>']
    o += flipped_bracket(x0, y0, k, bb)
    return "".join(o)


def flipped_bracket(x0, y0, k, bb):
    """Fig 1d's square bracket over U8, A9 and U10, placed off the probe.

    The probe render is the same camera and the same pixel frame as the
    published one, so a pixel in it maps into the placed image by subtracting
    the full render's ink origin and scaling by the fit factor. Arms point
    down at the bases, as the cassette's module brackets point down at the
    strand.
    """
    if not HP_PROBE.exists():
        return []
    from PIL import Image
    pb = Image.open(HP_PROBE).getchannel("A").getbbox()
    if pb is None:
        return []
    bx0, bx1 = (x0 + (pb[i] - bb[0]) * k for i in (0, 2))
    top = y0 + (pb[1] - bb[1]) * k
    arm, sw = 2.6 * FLIP_K, 0.35 * FLIP_K
    y_arm = top - FLIP_GAP                       # arm tips, just clear of the ink
    y_sp = y_arm - arm                           # the spine
    assert y_sp - FLIP_PT > 6, "the bracket label runs off the top of the panel"
    assert bx1 > bx0, "the probe has no width - wrong residues, or a stale render"
    return [f'<path data-name="Flipped bases bracket" '
            f'd="M{bx0:.2f},{y_arm:.2f} L{bx0:.2f},{y_sp:.2f} '
            f'L{bx1:.2f},{y_sp:.2f} L{bx1:.2f},{y_arm:.2f}" fill="none" '
            f'stroke="{MUT}" stroke-width="{sw:.2f}" stroke-linecap="round" '
            f'stroke-linejoin="round"/>',
            T((bx0 + bx1) / 2, y_sp - 5.0, FLIP_LABEL, FLIP_PT, MUT)]


def build():
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">', f'<rect width="{W}" height="{H}" fill="{BG}"/>',
         f'<rect x="6" y="6" width="{W-12}" height="{H-12}" rx="14" fill="{GROUND}"/>']
    o.append(hairpin())
    o.append(T(HP_CX, YM + HP_HALF + 22, "Csy4 crRNA hairpin", 12.5, INK, 600))
    o.append(harrow(X(152), X(210), YM))

    o.append(glyphs.icon_group("diffusion", X(282), YT, TEAL))
    o.append(T(X(282), YT - 66, "DIFFUSION", 11, TEAL, 700, "middle", 1.3))
    o.append(T(X(282), YT + 72, "RFdiffusion3", 12, INK))
    o.append(T(X(282), YT + 89, "BoltzGen", 12, INK))

    o.append(glyphs.icon_group("hallucination", X(282), YB, AMB))
    o.append(T(X(282), YB - 66, "HALLUCINATION", 11, AMB, 700, "middle", 1.3))
    o.append(T(X(282), YB + 72, "BoltzDesign1", 12, INK))
    o.append(T(X(282), YB + 89, "HalluDesign", 12, INK))

    for y0 in (YT, YB):
        o.append(f'<path d="M{X(340):.1f},{y0} C{X(396):.1f},{y0} {X(414):.1f},{YM} '
                 f'{X(462):.1f},{YM}" fill="none" stroke="{GREY}" stroke-width="2.6"/>')
    o.append(f'<path d="M{X(470):.1f},{YM} L{X(462):.1f},{YM-5.8} L{X(462):.1f},{YM+5.8} Z" '
             f'fill="{GREY}"/>')

    o.append(glyphs.icon_group("inverse-folding", X(536), YM, IND))
    o.append(T(X(536), YM + 78, "Inverse folding", 12.5, INK, 600))
    for i, tool in enumerate(("LigandMPNN", "BoltzIF")):   # one tool per line, as in
        o.append(T(X(536), YM + 96 + 17 * i, tool, 11, MUT))   # the generation column
    o.append(harrow(X(598), X(656), YM))

    o.append(glyphs.icon_group("refolding", X(720), YM, IND))
    o.append(T(X(720), YM + 78, "Re-prediction", 12.5, INK, 600))
    for i, tool in enumerate(("AF3", "Boltz-2", "ESMFold2")):   # on one line, these
        o.append(T(X(720), YM + 96 + 17 * i, tool, 11, MUT))     # met LigandMPNN's at ~5.4 pt
    o.append(harrow(X(782), X(840), YM))

    cy = YM
    o.append(f'<path d="M{X(850):.1f},{cy-78} C{X(992):.1f},{cy-78} {X(1092):.1f},{cy-28} '
             f'{X(1272):.1f},{cy-24}" fill="none" stroke="{LGREY}" stroke-width="1.6"/>')
    o.append(f'<path d="M{X(850):.1f},{cy+78} C{X(992):.1f},{cy+78} {X(1092):.1f},{cy+28} '
             f'{X(1272):.1f},{cy+24}" fill="none" stroke="{LGREY}" stroke-width="1.6"/>')
    o.append(cloud(X(914), cy, 60, 70, 330, LGREY, 11))
    o.append(cloud(X(1080), cy, 32, 36, 66, GREY, 12))
    o.append(cloud(X(1212), cy, 19, 19, 26, IND, 13))
    o.append(T(X(1061), cy - 100, "Five-criterion gate, then structural clustering",
               12.5, INK, 600))
    for x, big, small, col, wt in ((914, "6,992,515", "designs", MUT, 600),
                                   (1080, "55,706", "pass the gate", MUT, 600),
                                   (1212, "12,000", "selected", IND, 700)):
        o.append(T(X(x), cy + 100, big, 13, col, wt))
        o.append(T(X(x), cy + 117, small, 11, MUT))
    o.append(harrow(X(1284), X(1338), YM))

    o.append(glyphs.icon_group("oligo-pool", X(1404), YM, IND))
    o.append(T(X(1404), YM + 78, "Pooled oligo library", 12.5, INK, 600))
    o.append('</svg>')
    return "\n".join(o)


if __name__ == "__main__":
    svg_path = out / "Fig2a_pipeline.svg"
    svg_path.write_text(build(), encoding="utf-8")
    glyphs.write_icons(out / "icons")
    try:
        import cairosvg
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(out / "Fig2a_pipeline.pdf"))
        cairosvg.svg2png(url=str(svg_path), write_to=str(out / "Fig2a_pipeline.png"),
                         output_width=int(W * DPI / 96))
    except ImportError:
        # Drop any raster left by an earlier run that had cairosvg. Keeping it
        # would leave a stale PDF/PNG beside a fresh SVG, which is worse than
        # having no raster at all: the two disagree and nothing says so.
        for stale in (out / "Fig2a_pipeline.pdf", out / "Fig2a_pipeline.png"):
            stale.unlink(missing_ok=True)
        print("SVG only (no cairosvg). The SVG is the deliverable; raster from "
              "the vector editor if a preview is wanted.")
    print("wrote", out)
