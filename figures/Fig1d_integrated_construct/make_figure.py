"""Fig 1d, assembled: the landing-pad cell and everything read out of it.

Four pieces, each already built and each drawn at its own scale:

  left   `Fig1d_cassette_integration.svg`   donor -> landing pad -> integrated
  middle the landing-pad cell
  top r. `Fig1d_panel_right_translation.svg` the transcript, spliced and not
  bot r. `Fig1d_panel_pcr.svg`               the readout and the ratio

**Everything is normalised on text size.** Each source sets type at a
different effective size - 24.2 units in the cassette, 9.29 in the
translation panel, 5.5 in the PCR panel - so each is scaled by TEXT over its
own, and all captions land at the same size in the result.

That also matches the drawings, which is the point. In all three panels a
construct block is 1.091x the caption height, so a scale that equalises the
type equalises the blocks, the hairpins and the strands with it. The ratio is
asserted at build time: if a panel is ever redrawn against a different type
size, this stops rather than silently shipping mismatched parts.

Two things are deliberately NOT normalised: the mechanism inside the
translation panel, which has no construct in it to match, and the "Splicing
ratio" heading, which is a heading and is meant to outrank the captions.

    python code/Fig1d_integrated_construct/make_figure.py [outdir]
"""
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))

import make_panel as MP                                      # noqa: E402

WORK = MP.WORK
LIB = _HERE.parents[2] / "output" / "cartoon_library" / "svg"
NEUT = MP.NEUT
THEME7, THEME3 = "#015d82", "#95c9e8"
INK = "#1d2428"

TEXT = 5.5             # every caption in the finished figure
BLOCK_PER_TEXT = 6.0 / 5.5     # a construct block against its caption

# each panel's own effective caption size, measured from the file
SOURCES = {
    "cassette": ("Fig1d_cassette_integration.svg", 5.5 * 4.4),
    "translation": ("Fig1d_panel_right_translation.svg", 32.582 * 0.285),
    "pcr": ("Fig1d_panel_pcr.svg", 5.5),
}

PAD = 7.0              # inside a dashed group box
GAP = 16.0             # between the two right-hand rows
GAP_C = 15.0           # cell to each frame - equal, so it stays centred
CELL_W = 62.0
ARROW_W, EXTRACT_W, EXTRACT_H = 30.0, 64.0, 26.0


def _vb(t):
    return [float(v) for v in
            re.search(r'viewBox="([^"]+)"', t).group(1).replace(",", " ").split()]


def _prefix(t, tag):
    """Namespace a file's CSS classes so two panels cannot repaint each other."""
    t = re.sub(r'\.(?=[A-Za-z_-])([\w-]+)\s*\{', rf'.{tag}-\1{{', t)
    return re.sub(r'class="([^"]*)"',
                  lambda m: 'class="' + " ".join(f"{tag}-{c}" for c in
                                                 m.group(1).split()) + '"', t)


# Each source's INK, not its canvas, in its own units. A panel's own file
# carries margin - the cassette most of all, 59 units of it below the last
# row - and placing a canvas puts that margin inside the frame around it, so
# the frame looks loose and the drawing looks lost in it. Measured by
# rendering each file and trimming; re-measure if a panel is redrawn:
#   rsvg-convert -w 2400 -b white F.svg -o t.png
#   magick t.png -fuzz 0.5% -trim -format "%w %h %X %Y" info:
INK = {"cassette": (15.80, 22.63, 830.39, 642.19),
       "translation": (0.00, 0.68, 389.06, 199.49),
       "pcr": (-3.60, 7.40, 238.94, 112.08)}
BLEED = 0.6            # so a stroke on the edge is not clipped


def _ink(key, k):
    x, y, w, h = INK[key]
    return (w + 2 * BLEED) * k, (h + 2 * BLEED) * k


def _panel(key, x, y, k):
    """One source file placed at (x, y), scaled by k and cropped to its ink."""
    name, _ = SOURCES[key]
    t = (WORK / name).read_text()
    ix, iy, _, _ = INK[key]
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    body = _prefix(body, key[:3])
    return (f'<g data-name="{key}" transform="translate({x:.3f},{y:.3f}) '
            f'scale({k:.6f}) translate({BLEED - ix:.3f},{BLEED - iy:.3f})">'
            f'{body}</g>',) + _ink(key, k)


def _cell(x, y, w):
    """One cell, lifted out of the three in `cell_cluster.svg`.

    The cluster's largest cell is a blob and a filled nucleus and nothing
    else, which is what this wants: `landingpad_cell.svg` draws the cassette
    inside the nucleus, and at this size that reads as a third, tiny copy of
    the construct rather than as a cell.

    Returns the glyph and the two points the panels zoom out of - the nucleus,
    where the cassette integrates, and a patch of cytoplasm on the side the
    transcript panel sits. Both are read off the glyph, not typed in.
    """
    t = (LIB / "cell_cluster.svg").read_text()
    parts = re.findall(r'<(?:path|circle)[^>]*/>', t)
    blob, nucleus = parts[4], parts[5]      # the big cell of the three
    n = [float(v) for v in
         re.findall(r'-?\d+\.?\d*', re.search(r'\sd="([^"]+)"', blob).group(1))]
    x0, x1 = min(n[0::2]), max(n[0::2])
    y0, y1 = min(n[1::2]), max(n[1::2])
    nx, ny, nr = (float(re.search(rf'{a}="([\d.]+)"', nucleus).group(1))
                  for a in ("cx", "cy", "r"))
    k = w / (x1 - x0)
    blob = re.sub(r'fill="#[0-9a-fA-F]+"', f'fill="{NEUT["200"]}"', blob)
    blob = re.sub(r'stroke="#[0-9a-fA-F]+"', f'stroke="{NEUT["500"]}"', blob)
    nucleus = re.sub(r'fill="#[0-9a-fA-F]+"', f'fill="{THEME7}"', nucleus)
    at = lambda px, py: (x + (px - x0) * k, y + (py - y0) * k)
    # the ring sits in the cytoplasm band, nearer the nucleus than the
    # membrane: the blob is irregular, and its extreme x is not at this height
    cyto_x = nx + nr + (x1 - (nx + nr)) * 0.44
    return (f'<g data-name="Landing-pad cell" transform="translate({x:.2f},'
            f'{y:.2f}) scale({k:.5f}) translate({-x0:.2f},{-y0:.2f})">'
            f'{blob}{nucleus}</g>', w, (y1 - y0) * k,
            # just OUTSIDE the nucleus: a dashed ring drawn on top of that
            # solid fill is nearly invisible, against the cytoplasm it reads
            (at(nx, ny), nr * k * 1.16),
            (at(cyto_x, ny), (x1 - (nx + nr)) / 2 * k * 0.55))


def _cell_landingpad(x, y, w):
    """The landing-pad cell, the library glyph unchanged.

    Returns the placed glyph and, in composite coordinates, the two points the
    panels zoom out of: the nucleus, where the cassette is integrated, and a
    patch of cytoplasm, where it is transcribed and translated. Both come from
    the glyph's own paths, so they follow it if it is ever redrawn.
    """
    t = (LIB / "landingpad_cell.svg").read_text()
    vb = _vb(t)
    k = w / vb[2]
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    paths = re.findall(r'<path[^>]*\sd="([^"]+)"', t)

    def box(d):
        n = [float(v) for v in re.findall(r'-?\d+\.?\d*', d)]
        return min(n[0::2]), min(n[1::2]), max(n[0::2]), max(n[1::2])

    mem, nuc = box(paths[0]), box(paths[1])
    at = lambda px, py: (x + px * k, y + py * k)
    nucleus = at((nuc[0] + nuc[2]) / 2, (nuc[1] + nuc[3]) / 2)
    # the cytoplasm is the band between the two contours; take the right of it,
    # which is the side the transcript panel sits on
    # 0.40, not the midpoint: the membrane's extreme x is reached
    # above and below this height, so a ring centred halfway rides it
    cx = nuc[2] + (mem[2] - nuc[2]) * 0.40
    cyto = at(cx, (nuc[1] + nuc[3]) / 2)
    return (f'<g data-name="Landing-pad cell" '
            f'transform="translate({x:.2f},{y:.2f}) scale({k:.5f})">'
            f'{body}</g>', w, vb[3] * k,
            (nucleus, (nuc[2] - nuc[0]) / 2 * k * 0.46),
            (cyto, (mem[2] - nuc[2]) / 2 * k * 0.60))


def _dashed(x, y, w, h):
    return (f'<rect data-name="Zoom frame" x="{x:.2f}" y="{y:.2f}" '
            f'width="{w:.2f}" height="{h:.2f}" rx="2" fill="none" '
            f'stroke="{THEME3}" stroke-width="0.6" stroke-dasharray="1.2 1.8" '
            f'stroke-linecap="round"/>')


def _leaders(src, box, left):
    """A zoom: ring the compartment, then fan to the frame's near corners.

    `src` is ((x, y), r) - where in the cell this panel comes from. The ring
    is what makes it a zoom rather than a line to the cell in general: the
    left panel is the nucleus, the right one a patch of cytoplasm.
    """
    (cx, cy), r = src
    bx = box[0] + box[2] if left else box[0]
    dash = 'stroke-dasharray="1.2 1.8" stroke-linecap="round"'
    o = [f'<circle data-name="Zoom source" cx="{cx:.2f}" cy="{cy:.2f}" '
         f'r="{r:.2f}" fill="none" stroke="{THEME3}" stroke-width="0.5" '
         f'{dash}/>']
    for by in (box[1], box[1] + box[3]):
        dx, dy = bx - cx, by - cy
        n = (dx * dx + dy * dy) ** 0.5 or 1.0
        o.append(f'<path data-name="Leader" d="M{cx + dx / n * r:.2f},'
                 f'{cy + dy / n * r:.2f} L{bx:.2f},{by:.2f}" fill="none" '
                 f'stroke="{THEME3}" stroke-width="0.5" {dash}/>')
    return "".join(o)


def _extract(x, y):
    """The RNA step itself, without its arrows."""
    cx = x + EXTRACT_W / 2
    o = [f'<rect data-name="RNA step" x="{x:.2f}" y="{y:.2f}" '
         f'width="{EXTRACT_W}" height="{EXTRACT_H}" fill="{NEUT["100"]}" '
         f'stroke="{NEUT["500"]}" stroke-width="0.5"/>']
    for i, line in enumerate(("Extraction of RNA", "RT to cDNA")):
        o.append(f'<text text-anchor="middle" x="{cx:.2f}" '
                 f'y="{y + EXTRACT_H / 2 - 1.1 + i * 6.6:.2f}" '
                 f'font-family="Arial" font-size="{TEXT}" fill="{INK}">'
                 f'{line}</text>')
    return "".join(o)


def _right(x0, x1, y, col=None):
    """The same step arrow, turned."""
    col = col or NEUT["500"]
    head = 2.8
    return (f'<path data-name="Step" d="M{x0:.2f},{y:.2f} '
            f'L{x1 - head:.2f},{y:.2f}" fill="none" stroke="{col}" '
            f'stroke-width="0.7" stroke-linecap="round"/>'
            f'<path d="M{x1 - head:.2f},{y - head * 0.42:.2f} '
            f'L{x1 - head:.2f},{y + head * 0.42:.2f} L{x1:.2f},{y:.2f} Z" '
            f'fill="{col}"/>')


def _down(x, y0, y1, col=None):
    """A step arrow, in the weight the panels' own arrows use.

    Grey and hairline, not a heavy coloured arrow: these say "and then", the
    same thing the arrows inside the panels say, so they are drawn the same.
    """
    col = col or NEUT["500"]
    head = 2.8
    return (f'<path data-name="Step" d="M{x:.2f},{y0:.2f} '
            f'L{x:.2f},{y1 - head:.2f}" fill="none" stroke="{col}" '
            f'stroke-width="0.7" stroke-linecap="round"/>'
            f'<path d="M{x - head * 0.42:.2f},{y1 - head:.2f} '
            f'L{x + head * 0.42:.2f},{y1 - head:.2f} L{x:.2f},{y1:.2f} Z" '
            f'fill="{col}"/>')


def build():
    k = {key: TEXT / size for key, (_, size) in SOURCES.items()}

    cas_w, cas_h = _ink("cassette", k["cassette"])
    tr_w, tr_h = _ink("translation", k["translation"])
    pcr_w, pcr_h = _ink("pcr", k["pcr"])

    # The RNA step sits to the LEFT of the readout, not above it, which is
    # what lets the readout climb to just under the transcript panel: the
    # right column's height is then two panels, not two panels and a box.
    # Each frame hugs its own panel - a frame stretched to the figure's height
    # is mostly margin, and makes the drawing inside it look small.
    right_w = max(tr_w, pcr_w) + 2 * PAD
    cas_box_w, cas_box_h = cas_w + 2 * PAD, cas_h + 2 * PAD
    H = tr_h + 2 * PAD + GAP + pcr_h

    # the cell goes on the figure's centre line, so the gap either side of it
    # absorbs the difference between a one-panel half and a two-panel one
    # equal halves about the cell: cas_box_w + gap_l == right_w + gap_r
    gap_r = GAP_C
    gap_l = right_w - cas_box_w + gap_r
    W = cas_box_w + gap_l + CELL_W + gap_r + right_w
    x_cell = cas_box_w + gap_l
    x_right = x_cell + CELL_W + gap_r
    assert abs((x_cell + CELL_W / 2) - W / 2) < 1e-6, "cell off the centre line"

    cas_box = (0.0, (H - cas_box_h) / 2, cas_box_w, cas_box_h)
    tr_box = (x_right, 0.0, right_w, tr_h + 2 * PAD)

    cell, cw, ch, nucleus, cyto = _cell(x_cell, 0.0, CELL_W)
    shift = (tr_box[3] - ch) / 2        # level with the panel it feeds
    cell = cell.replace(",0.00) scale", f",{shift:.2f}) scale", 1)
    nucleus = ((nucleus[0][0], nucleus[0][1] + shift), nucleus[1])
    cyto = ((cyto[0][0], cyto[0][1] + shift), cyto[1])

    o = [_dashed(*cas_box), _dashed(*tr_box),
         _panel("cassette", PAD, cas_box[1] + PAD, k["cassette"])[0],
         cell,
         # after the cell: its body is opaque, and a ring drawn under it
         # simply disappears
         _leaders(nucleus, cas_box, True),
         _leaders(cyto, tr_box, False),
         _panel("translation", x_right + (right_w - tr_w) / 2, PAD,
                k["translation"])[0]]

    y_pcr = tr_box[3] + GAP
    pcr_x = x_right + (right_w - pcr_w) / 2
    o.append(_panel("pcr", pcr_x, y_pcr, k["pcr"])[0])
    # the RNA comes out of the cells, and the cell is directly above the step,
    # so the path bends once: down out of the cell, then right into the readout
    ex_y = y_pcr + (pcr_h - EXTRACT_H) / 2
    ex_x = x_cell + CELL_W / 2 - EXTRACT_W / 2
    o.append(_down(x_cell + CELL_W / 2, shift + ch + 2.0, ex_y - 1.5))
    o.append(_extract(ex_x, ex_y))
    o.append(_right(ex_x + EXTRACT_W + 1.5, pcr_x - 1.5, ex_y + EXTRACT_H / 2))

    m = 3.0
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{(W + 2 * m) * 5:.0f}" height="{(H + 2 * m) * 5:.0f}" '
            f'viewBox="{-m} {-m} {W + 2 * m:.2f} {H + 2 * m:.2f}">'
            + "".join(o) + "</svg>")


# A construct block's height in each source, in that file's own units and
# measured off it: Matteo's 6 inside the cassette's scale(4.4); the binder
# path in the translation panel; 6 in the PCR panel, which has no outer scale.
BLOCK = {"cassette": 6.0 * 4.4, "translation": 10.133, "pcr": 6.0}


def _check():
    """One scale has to serve both the type and the drawing, in every panel.

    It only does so because all three were authored with the same block-to-
    caption proportion. Assert it rather than trust it: if a panel is redrawn
    against a different type size, normalising on text would quietly leave its
    constructs the wrong size, which is exactly the fault this figure exists
    to avoid.
    """
    out = []
    for key, (_, size) in SOURCES.items():
        k = TEXT / size
        block = BLOCK[key] * k
        out.append((key, size, k, block))
        assert abs(BLOCK[key] / size - BLOCK_PER_TEXT) < 0.01, (
            f"{key} sets blocks at {BLOCK[key] / size:.4f} of its caption, "
            f"not {BLOCK_PER_TEXT:.4f}: one scale cannot match both")
    blocks = [b for *_, b in out]
    assert max(blocks) - min(blocks) < 0.02, blocks
    for key, size, k, block in out:
        print(f"  {key:12} caption {size:7.3f} -> x{k:.4f}   "
              f"block {block:.3f}")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    _check()
    p = outdir / "Fig1d_full.svg"
    p.write_text(build(), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else WORK)
