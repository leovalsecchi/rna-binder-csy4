"""Fig 1d, narrowed: Leo's 'Full Panel V2.svg' with its air taken out.

    python code/Fig1d_integrated_construct/compact_full_panel.py

Fig 1d sets the type size of the whole of Fig 1 (code/Fig1_assembly): it has
to fill the 178 mm block, so every unit of width it loses is type gained. The
cheapest width is the column the cell sits in - 503 units of mostly air between
the cassette frame and the right-hand frames. This moves the whole right-hand
column left, shrinks the cell about its centre by CELL_K and centres it in what
is left, and re-aims the four zoom leaders: each keeps its point on its ring
(moved with the cell) and its point on its frame (moved with the frame).

The frames give width too. Each was drawn with ~47-59 units of air between its
content and its sides; PAD_L and PAD_R take the cassette frame and the two
right-hand frames in by that much in total, content kept centred, so all
three end up with the same ~27 units each side.

Nothing else is touched; the input stays the master. Writes
'Full Panel V4 compact.svg' beside it, sized for 5.5 pt captions at 178 mm.
(V3, the 5.17 pt version, was DELTA 130 / CELL_K 0.82 / no padding trim.)
"""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WORK = REPO / "working" / "Fig1d_integrated_construct"
SRC, DST = WORK / "Full Panel V2.svg", WORK / "Full Panel V4 compact.svg"

COLUMN_CUT = 190.3     # the cell column (frame to frame) loses this many units
CELL_K = 0.75          # cell scale, about its own centre
PAD_L = 39.0           # cassette frame: 46.5 -> 27 units of air each side
PAD_R = 64.4           # right-hand frames: 59.2 -> 27 each side (the readout
                       # box shares their edges and keeps more inside)
RIGHT_SIDE = 2700.0    # x above this is a right-hand frame's right side
LEFT_EDGE = 1038.481   # cassette frame, right edge
CELL_BOX = (1156.8, 1466.7, 510.4, 816.0)    # x0 x1 y0 y1 of the cell's drawn extent

# The cell line, named under the cell (Leo, 2026-11-13). Drawn HERE, not in V2:
# V2's cell is shrunk by CELL_K inside its own group, so a label placed in it
# would print at 75% of the caption size, and a new top-level child of V2
# would shift the positions this script reads its parts by. Caption class,
# two lines to fit the cell column, and a white halo so the dashed zoom
# leader that crosses under it breaks cleanly instead of striking through.
CELL_LABEL = ["Landing-pad HEK293", "cell line"]
CAPTION, CAP_H = 27.496, 0.716 * 27.496
LABEL_GAP = 20.0            # cell's drawn bottom to the first line's cap top
HALO = 9.0                  # white stroke behind the text, file units


def children(lines, start, indent):
    """Top-level children of the group opened at lines[start], as line ranges."""
    pad, out, i = " " * indent, [], start + 1
    while not lines[i].startswith(pad[:-2] + "</g>"):
        j = i
        if lines[i].startswith(pad + "<g"):
            while not lines[j].startswith(pad + "</g>"):
                j += 1
        out.append((i, j + 1))
        i = j + 1
    return out


def main():
    lines = SRC.read_text().split("\n")
    top = next(i for i, l in enumerate(lines) if l.startswith('  <g id="Asset_30"')) + 1
    kids = children(lines, top, 6)
    txt = ["\n".join(lines[a:b]) for a, b in kids]

    # who is who, checked by content rather than trusted by position
    assert txt[0].lstrip().startswith('<path d="M1551.709 1.5 '), "right frame moved"
    assert txt[1].lstrip().startswith('<path d="M11.716 1.5 '), "left frame moved"
    assert "Promoterless_donor_plasmid" in txt[2], "cassette group moved"
    assert "M1465.872 672.941" in txt[3], "cell group moved"
    lead = re.compile(r'<path d="M([\d.]+) ([\d.]+) L([\d.]+) ([\d.]+)" class="cls-26"/>')
    assert all(lead.search(txt[i]) for i in (4, 5, 7, 8)), "leaders moved"
    assert txt[6].lstrip().startswith('<path d="M1443.912 678.014'), "cytoplasm ring moved"
    assert 'id="Asset_1"' in txt[9], "right column moved"

    x0, x1, y0, y1 = CELL_BOX
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    shift = PAD_L + COLUMN_CUT                  # how far the right column moves
    left_edge = LEFT_EDGE - PAD_L
    right_edge = 1541.711 - shift
    ncx = (left_edge + right_edge) / 2          # cell centred in the column
    gap = (right_edge - left_edge - (x1 - x0) * CELL_K) / 2
    assert gap > 30, f"cell column too tight: {gap:.1f} units each side"

    def cell(px, py):
        return ncx + CELL_K * (px - cx), cy + CELL_K * (py - cy)
    cell_t = (f'transform="translate({ncx:.3f} {cy:.3f}) scale({CELL_K}) '
              f'translate({-cx:.3f} {-cy:.3f})"')

    def reaim(s, frame_dx):
        m = lead.search(s)
        ax, ay = cell(float(m.group(1)), float(m.group(2)))
        bx, by = float(m.group(3)) + frame_dx, float(m.group(4))
        return lead.sub(f'<path d="M{ax:.3f} {ay:.3f} L{bx:.3f} {by:.3f}" class="cls-26"/>', s)

    def narrow(s, threshold, dx):
        """A frame path with every point right of threshold moved by dx."""
        d = re.search(r' d="([^"]+)"', s).group(1)
        i = iter(range(10 ** 6))
        nd = re.sub(r"-?\d+\.?\d*", lambda m: (
            f"{float(m.group(0)) + dx:.3f}" if next(i) % 2 == 0 and float(m.group(0)) > threshold
            else m.group(0)), d)
        return s.replace(d, nd, 1)

    g = lambda dx, body: [f'<g transform="translate({dx:.3f} 0)">'] + body + ["</g>"]
    READOUT = '<path d="M1556.21 779.391 '
    assert sum(k.lstrip().startswith(READOUT) for k in txt[9:]) == 1, "readout box moved"
    right = [narrow(k, RIGHT_SIDE, -PAD_R) if k.lstrip().startswith(READOUT) else None
             for k in txt[9:]]
    new = (
        g(-shift, [narrow(txt[0], RIGHT_SIDE, -PAD_R)])
        + [narrow(txt[1], 1000.0, -PAD_L)] + g(-PAD_L / 2, [txt[2]])
        + [f"<g {cell_t}>", txt[3], txt[6], "</g>",
           reaim(txt[4], -PAD_L), reaim(txt[5], -PAD_L),
           reaim(txt[7], -shift), reaim(txt[8], -shift)])
    for k, framed in zip(txt[9:], right):     # frames move; contents move and re-centre
        new += g(-shift, [framed]) if framed else g(-shift - PAD_R / 2, [k])

    sys.path.insert(0, str(REPO / "code"))
    from common.style import text_w
    widest = max(text_w(s, CAPTION) for s in CELL_LABEL)
    assert widest + 2 * HALO < right_edge - left_edge - 20, \
        f"cell label {widest:.0f} wide does not fit the {right_edge - left_edge:.0f} column"
    base = cy + CELL_K * (y1 - cy) + LABEL_GAP + CAP_H
    for i, s in enumerate(CELL_LABEL):
        new.append(f'      <text data-name="Cell line label" transform="matrix(1 0 0 1 '
                   f'{ncx:.3f} {base + i * 1.2 * CAPTION:.3f})" x="0" y="0" class="cls-10" '
                   f'style="paint-order:stroke;stroke:#ffffff;stroke-width:{HALO};'
                   f'stroke-linejoin:round"><tspan x="0" y="0">{s}</tspan></text>')
    a, b = kids[0][0], kids[-1][1]
    out = lines[:a] + new + lines[b:]
    t = "\n".join(out)

    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]
    w_new = vb[2] - shift - PAD_R
    h_pt = 504.567 * vb[3] / w_new
    t = re.sub(r'width="[^"]+" height="[^"]+" viewBox="[^"]+"',
               f'width="504.567pt" height="{h_pt:.3f}pt" '
               f'viewBox="{vb[0]:g} {vb[1]:g} {w_new:.3f} {vb[3]:g}"', t, count=1)
    DST.write_text(t, encoding="utf-8")
    t_pt = 27.496 * 504.567 / w_new
    print(f"wrote {DST.name}: {w_new:.1f} units wide (was {vb[2]:.1f}); "
          f"cell gaps {gap:.1f}; captions at 178 mm = {t_pt:.2f} pt")


if __name__ == "__main__":
    main()
