"""Re-set a panel's live text to one printed size at page assembly.

A panel placed at scale s prints a font-size of n units at n * s pt, so the
house size T needs T / s units. Every size goes there except named riders
(deliberate exceptions, as {size in the file: multiple of T}). Stacked labels -
same x, next line close below - keep their leading in proportion, and
letter-spacing keeps its proportion to its own size, so a label block re-set
smaller or larger still reads as set, not stretched.

Keyed on what the file contains, not on what it was drawn at, so a panel
re-run at a new house size needs no change here.

ROTATED labels are re-centred as well (2026-11-19). matplotlib bakes a label's
alignment into its anchor point rather than emitting text-anchor, and a rotated
label's anchor is the END of its baseline, so re-setting the size keeps that end
pinned and slides the whole word along its own length. Fig 2c's "Hallucination"
drifted ~1 mm down its bracket that way - the word ran past the bottom tick with
a gap left at the top, which is what it looked like on the page while the SVG's
numbers all looked right. This assumes a rotated label was CENTRED when drawn,
which is true of every one in this book (axis titles and group labels).
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import style                                  # noqa: E402

SIZE = r'(font-size(?::\s*|="))([0-9.]+)'


def retext(t, T, s, riders=None):
    riders = riders or {}

    def new(old):
        k = next((m for r, m in riders.items() if abs(r - old) < 0.01), 1.0)
        return k * T / s

    t = _lead(t, new)
    t = _recentre_rotated(t, new)
    t = re.sub(r'(font-size="([0-9.]+)"[^>]*?letter-spacing=")([0-9.]+)',
               lambda m: m.group(1) + f"{float(m.group(3)) * new(float(m.group(2))) / float(m.group(2)):.3f}", t)
    return re.sub(SIZE, lambda m: f"{m.group(1)}{new(float(m.group(2))):.4f}", t)


# matplotlib: <text style="font-size: 5.5px; ..." transform="translate(X Y) rotate(-90)">word</text>
ROT = re.compile(r'<text[^>]*?font-size:\s*([0-9.]+)[^>]*?'
                 r'transform="translate\(([-0-9.]+) ([-0-9.]+)\) rotate\((-?90)\)"[^>]*>([^<]*)</text>')


def _recentre_rotated(t, new):
    """Hold a rotated label's CENTRE still while its size changes.

    The anchor sits at one end of the baseline, so a label re-set from f0 to f1
    keeps that end and loses (w0 - w1) off the other. Advance width is linear in
    the size, so the shift is half that difference, signed by which way the text
    runs: rotate(-90) advances UP the page (-y), rotate(90) down it.
    """
    def fix(m):
        f0, word = float(m.group(1)), m.group(5)
        f1 = new(f0)
        if abs(f1 - f0) < 1e-9 or not word.strip():
            return m.group(0)
        w0 = style.text_w(word, f0)                # pt, and 1 SVG unit = 1 pt here
        half = w0 * (f1 / f0 - 1) / 2
        dy = half if m.group(4) == "-90" else -half
        y = float(m.group(3)) + dy
        return m.group(0).replace(f"translate({m.group(2)} {m.group(3)})",
                                  f"translate({m.group(2)} {y:.6f})")
    return ROT.sub(fix, t)


def allowed(T, riders=None):
    """The printed sizes a panel re-set by retext() may contain."""
    return {T} | {m * T for m in (riders or {}).values()}


def _lead(t, new):
    rows = {}
    for m in re.finditer(r'<text x="([0-9.]+)" y="([0-9.]+)"[^>]*?font-size="([0-9.]+)"', t):
        rows.setdefault(m.group(1), []).append(
            (float(m.group(2)), m.start(2), m.end(2), float(m.group(3))))
    moved = {}
    for col in rows.values():
        col.sort()
        for (y0, a0, _, s0), (y1, a1, b1, _) in zip(col, col[1:]):
            if y1 - y0 < 2.0 * s0:                       # the next line of one label
                moved[(a1, b1)] = moved.get((a0, None), y0) + (y1 - y0) * new(s0) / s0
                moved[(a1, None)] = moved[(a1, b1)]
    out, last = [], 0
    for (a, b), y in sorted((k, v) for k, v in moved.items() if k[1] is not None):
        out += [t[last:a], f"{y:.2f}"]
        last = b
    return "".join(out) + t[last:]


if __name__ == "__main__":                               # the smallest check that fails if this breaks
    src = ('<text x="10" y="100" font-size="10">A</text>'
           '<text x="10" y="114" font-size="10" letter-spacing="1">B</text>'
           '<text x="50" y="100" font-size="12">C</text>')
    out = retext(src, 5.0, 0.25)                         # 5 pt at 0.25 pt/unit = 20 units
    assert 'y="128.00"' in out, out                      # 14-unit lead on 10 -> 28 on 20
    assert 'letter-spacing="2.000"' in out, out
    assert out.count('font-size="20.0000"') == 3, out
    assert retext(src, 5.0, 0.25, {12: 1.2}).count('font-size="24.0000"') == 1

    # a rotated label keeps its CENTRE, not its anchor, when it is re-set
    rot = ('<text style="font-size: 10px; font-family: \'Arial\'" '
           'transform="translate(8 100) rotate(-90)">Hallucination</text>')
    w0 = style.text_w("Hallucination", 10.0)
    out = retext(rot, 5.0, 1.0)                          # 10 -> 5, so it halves
    y = float(out.split("translate(8 ")[1].split(")")[0])
    assert abs((y - w0 / 4) - 100 + w0 / 2) < 1e-6, (y, w0)   # centre held: y - w1/2 == 100 - w0/2
    assert abs((100 - w0 / 2) - (y - w0 / 4)) < 1e-6
    # and rotate(90) runs the other way, so it moves the other way
    y90 = float(retext(rot.replace("rotate(-90)", "rotate(90)"), 5.0, 1.0)
                .split("translate(8 ")[1].split(")")[0])
    assert abs((y90 - 100) + (y - 100)) < 1e-6, (y, y90)
    print("retext checks pass")
