"""The SMS receptor, drawn from scratch: binder, kinase lobe, RNase lobe.

Our own drawing - nothing is imported from the shared cartoon set. Flat fills
with one outline weight, as COLORS.md rule 6 asks of a synbio element, each
part in the colour that rule gives it:

    binder   synbio.binder       the designed molecule, so MAIN petrol
    kinase   synbio.enzyme       the lobe that dimerizes
    RNase    synbio.enzyme_alt   the lobe that cuts the RNA

The shapes follow the reference sketch rather than any structure:

* the **binder** is a rounded slab with its near end scalloped into two
  concave scoops - a binding face, deliberately not a fold. A helix bundle was
  tried first and read as "some protein" instead of "the thing that grips";
* the **kinase** is a tapered capsule, two circles of different radius with
  their outer tangents, so it thickens toward the linker;
* the **RNase** is a disc with a round bite out of it, the cleft the RNA sits
  in. The bite is a circle subtracted, not a wedge, so its edge is an arc.

Binder and kinase are one polypeptide - the binder is fused N-terminal to
MgIRE1 - so the linker is chain, a wave rather than a rule, and longer than the
sketch's so it reads as a linker and not a seam.

    python code/cartoon_library/receptor.py <outdir>
    python code/cartoon_library/receptor.py --check
"""
import re
import sys
from math import acos, atan2, cos, hypot, pi, sin
from pathlib import Path
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parent))

from common import palette
from hairpin import spline, wave

S = palette.SYNBIO
LINE = S["outline"]
SW = 3.2           # one outline weight for every part

# binder: a slab, scalloped at the end that does the gripping
SLAB_L, SLAB_T, SLAB_R = 104.0, 54.0, 13.0
SLAB_ANG = -0.42            # its long axis, radians from +x
NOTCH = 21.0                # how deep the two scoops bite into the end

# the linker, from under the slab to the head of the kinase
LINK_LEN, LINK_WAVE, LINK_AMP = 62.0, 31.0, 6.5

# kinase: a tapered capsule
KIN_LEN = 74.0
KIN_R1, KIN_R2 = 27.0, 22.0   # thick at the linker end, thinner at the RNase
KIN_ANG = 2.0                 # radians from +x: down and to the left

# RNase: a disc with a round bite out of it
RNASE_R, BITE_R = 31.0, 17.0
BITE_ANG = 2.55               # where the bite sits, radians from +x
BITE_REACH = 0.80             # bite centre, as a fraction of RNASE_R out

PAD = 10.0
_pts = []


def _p(x, y):
    _pts.append((x, y))
    return x, y


def _at(c, ang, r):
    return _p(c[0] + r * cos(ang), c[1] + r * sin(ang))


def _shape(d, fill):
    return (f'<path d="{d}" fill="{fill}" stroke="{LINE}" stroke-width="{SW}" '
            f'stroke-linejoin="round" stroke-linecap="round"/>')


def slab(cx, cy, ang=SLAB_ANG, L=SLAB_L, T=SLAB_T, r=SLAB_R, notch=NOTCH):
    """A rounded slab whose near end is scalloped into two concave scoops."""
    ca, sa = cos(ang), sin(ang)
    def q(u, v):                       # slab frame -> drawing, origin at centre
        x, y = u - L / 2, v - T / 2
        return _p(cx + x * ca - y * sa, cy + x * sa + y * ca)
    a, b = q(0, 0), q(L - r, 0)
    c, d = q(L, r), q(L, T - r)
    e, f = q(L - r, T), q(0, T)
    m = q(0, T / 2)
    k1, k2 = q(notch, T * 0.74), q(notch, T * 0.26)
    return _shape(
        f"M{a[0]:.1f},{a[1]:.1f} L{b[0]:.1f},{b[1]:.1f} "
        f"A{r},{r} 0 0 1 {c[0]:.1f},{c[1]:.1f} L{d[0]:.1f},{d[1]:.1f} "
        f"A{r},{r} 0 0 1 {e[0]:.1f},{e[1]:.1f} L{f[0]:.1f},{f[1]:.1f} "
        f"Q{k1[0]:.1f},{k1[1]:.1f} {m[0]:.1f},{m[1]:.1f} "
        f"Q{k2[0]:.1f},{k2[1]:.1f} {a[0]:.1f},{a[1]:.1f} Z", S["binder"])


def linker(start, ang, length=LINK_LEN):
    """Chain, not a rule: a wave that starts and ends on its axis."""
    segs = wave(length, LINK_WAVE, LINK_AMP, ang)
    for seg in segs:
        for u, v in seg:
            _p(start[0] + u, start[1] + v)
    d = spline([tuple((start[0] + u, start[1] + v) for u, v in seg)
                for seg in segs])
    return (f'<path d="{d}" fill="none" stroke="{LINE}" stroke-width="{SW}" '
            f'stroke-linecap="round"/>')


def taper(a, ra, b, rb):
    """Two circles of different radius closed by their outer tangents."""
    d = hypot(b[0] - a[0], b[1] - a[1])
    base = atan2(b[1] - a[1], b[0] - a[0])
    t = acos((ra - rb) / d)            # tangent angle off the centre line
    p1, p2 = _at(a, base + t, ra), _at(b, base + t, rb)
    p3, p4 = _at(b, base - t, rb), _at(a, base - t, ra)
    return (f"M{p1[0]:.1f},{p1[1]:.1f} L{p2[0]:.1f},{p2[1]:.1f} "
            f"A{rb},{rb} 0 1 0 {p3[0]:.1f},{p3[1]:.1f} "
            f"L{p4[0]:.1f},{p4[1]:.1f} A{ra},{ra} 0 1 0 {p1[0]:.1f},{p1[1]:.1f} Z")


def bitten_disc(c, R, bite, br):
    """A disc with a round bite: the cut edge is the bite circle's arc."""
    d = hypot(bite[0] - c[0], bite[1] - c[1])
    base = atan2(bite[1] - c[1], bite[0] - c[0])
    t = acos((d * d + R * R - br * br) / (2 * d * R))      # half the bite's span
    p1, p2 = _at(c, base - t, R), _at(c, base + t, R)
    _p(c[0] - R, c[1] - R), _p(c[0] + R, c[1] + R)
    return (f"M{p1[0]:.1f},{p1[1]:.1f} "
            f"A{R},{R} 0 1 0 {p2[0]:.1f},{p2[1]:.1f} "
            f"A{br},{br} 0 0 1 {p1[0]:.1f},{p1[1]:.1f} Z")


def build():
    _pts.clear()
    o = []
    # the slab, and the point on its lower long edge the chain leaves from
    sc = (0.0, 0.0)
    o.append(slab(sc[0], sc[1]))
    foot = _p(sc[0] - (SLAB_L / 2 - 26) * cos(SLAB_ANG)
              - (SLAB_T / 2) * sin(-SLAB_ANG),
              sc[1] - (SLAB_L / 2 - 26) * sin(SLAB_ANG)
              + (SLAB_T / 2) * cos(SLAB_ANG))
    o.append(linker(foot, KIN_ANG))

    head = _p(foot[0] + LINK_LEN * cos(KIN_ANG),
              foot[1] + LINK_LEN * sin(KIN_ANG))
    a = _p(head[0] + KIN_R1 * cos(KIN_ANG), head[1] + KIN_R1 * sin(KIN_ANG))
    b = _p(a[0] + KIN_LEN * cos(KIN_ANG), a[1] + KIN_LEN * sin(KIN_ANG))

    # the RNase first, so the kinase covers where they meet
    # far enough past the kinase that the disc and its bite both read; much
    # less and the lobe swallows it
    reach = RNASE_R * 0.95
    rc = _p(b[0] + reach * cos(KIN_ANG), b[1] + reach * sin(KIN_ANG))
    bite = (rc[0] + RNASE_R * BITE_REACH * cos(BITE_ANG),
            rc[1] + RNASE_R * BITE_REACH * sin(BITE_ANG))
    o.append(_shape(bitten_disc(rc, RNASE_R, bite, BITE_R), S["enzyme_alt"]))
    o.append(_shape(taper(a, KIN_R1, b, KIN_R2), S["enzyme"]))

    m = PAD + SW
    xs, ys = [p[0] for p in _pts], [p[1] for p in _pts]
    x0, y0 = min(xs) - m, min(ys) - m
    w, h = max(xs) - min(xs) + 2 * m, max(ys) - min(ys) + 2 * m
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" '
            f'height="{h:.0f}" viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}">'
            + "".join(o) + "</svg>")


W, H = 230, 300          # nominal only; build() crops to the drawing


def _selfcheck():
    svg = build()
    cols = {c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", svg)}
    assert cols == {c.lower() for c in
                    (S["binder"], S["enzyme"], S["enzyme_alt"], S["outline"])}, \
        f"colours are {sorted(cols)}, not the four synbio tokens"
    # the scoops have to bite without meeting in the middle of the slab
    assert 0 < NOTCH < SLAB_L / 3, "the binder's scoops eat the slab"
    # a taper needs the radius difference to be smaller than the span, or
    # there is no outer tangent and acos() blows up
    assert abs(KIN_R1 - KIN_R2) < KIN_LEN, "no outer tangent for the kinase"
    # the bite must cut the rim and not swallow the disc or miss it
    d = RNASE_R * BITE_REACH
    assert RNASE_R - BITE_R < d < RNASE_R + BITE_R, \
        "the bite either misses the RNase or splits it"
    # the linker arrives on its axis
    assert abs(LINK_LEN / (LINK_WAVE / 2)
               - round(LINK_LEN / (LINK_WAVE / 2))) < 1e-9, \
        "the linker is not a whole number of half waves"
    print(f"receptor: checks pass, {LINK_LEN / LINK_WAVE:.1f} waves on the linker")


if __name__ == "__main__":
    if "--check" in sys.argv:
        _selfcheck()
    else:
        out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
        (out / "sms_receptor.svg").write_text(build(), encoding="utf-8")
        print("wrote sms_receptor.svg to", out)
