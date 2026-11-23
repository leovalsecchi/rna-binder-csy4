"""The AT (XBP1) hairpin, unspliced and spliced - traced off the drawing itself.

Every number below was measured from the hairpins in the existing panel b, at
the pixel level (`inputs/cartoons/` holds five versions of a *different*, compact
hairpin - arms at 45 degrees, two base pairs each - that is not what that panel
shows, so none of them could be lifted). The drawing is wide and flat: two
horizontal arms, a round loop at each end, a small loop over the junction, and a
long stem below it.

The measurement is in the reference image's own pixels, stem centre at x = 0,
and SCALE blows it up; so a constant here can be checked against the picture
with a ruler, which is the point. `_selfcheck()` re-asserts the ones that the
shape depends on.

Two departures from the reference, both asked for:

* **Nothing outside the hairpin.** The reference runs the transcript out to
  either side as a wavy line and flares the stem into it; the rails here stop
  square at the bottom of the stem.
* **Thesis palette.** Exons carry theme 400 - the colour `constructs.py` gives
  the AT module, which COLORS.md ties to this drawing - and the excised intron a
  darker step of the same hue, so it reads as the same molecule. In the
  reference that contrast is black against red.

    python code/cartoon_library/hairpin.py            # into the library
    python code/cartoon_library/hairpin.py --check
"""
import shutil
import subprocess
import sys
from math import cos, pi, radians, sin, sqrt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import palette

# Red is what the hairpin loses, blue what it keeps: the two drawings differ by
# the red, and the loose strand beside them is that red piece. Theme 600 for the
# exons rather than MAIN or theme 400, both of which are taken in this panel -
# MAIN is the binder, theme 400 the RNA target (COLORS.md rule 6).
EXON = palette.THEME["600"]       # the message either side of the cut: it stays
INTRON = palette.SYNBIO["intron"]    # the piece that is cut out and drifts off
TICK = palette.NEUTRAL["400"]     # structure is neutral (COLORS.md rule 6)

SW = 1.30         # rail stroke, in reference pixels
RUNG = 0.90       # a base pair
SCALE = 7.0       # reference pixels to output units

# --- unspliced: measured from the reference ---------------------------------
ARM_Y = 35.50     # the arms' axis
ARM_HALF = 2.00   # half the rail separation along an arm (rails 4.0 apart)
SIDE_CX = 30.00   # the two end loops, left and right of the stem centre
SIDE_R = 7.20
MID_CY = 32.00    # the small loop sitting over the junction
MID_R = 6.00
STEM_HALF = 3.50  # the stem is wider than an arm: it spans the junction gap
CORNER = 2.50     # where an arm's lower rail turns down into the stem
STEM_TOP = 40.00
STEM_BOT = 57.00
ARM_RUNGS, ARM_RUNG_0, ARM_RUNG_D = 4, 4.00, 3.00
STEM_RUNGS, STEM_RUNG_0, STEM_RUNG_D = 8, 41.50, 2.00

# --- spliced: the same, from the hairpin beside it ---------------------------
SP_R = 8.00
SP_TOP = 32.00    # where the loop closes onto the stem rails
SP_BOT = 59.00
SP_RUNGS, SP_RUNG_0, SP_RUNG_D = 9, 34.50, 2.75

# --- the excised intron, drifting free --------------------------------------
# A sine, not a trace. The loose strand in the panel was traced first and it
# came out janky: a long flat plateau, one steep ramp and a flat top, with the
# bases splaying where the ramp turned. Once cut this piece is unpaired and
# free, so there is nothing to trace faithfully - a plain wave is both truer
# and tidier. The footprint is kept close to the traced one so it still sits
# with the hairpins.
#
# Knots every eighth wavelength, Hermite handles at h/3 along the tangent:
# that reproduces a sine to a few hundredths of a percent, and the handles either
# side of a knot are collinear through it by construction, so the curve has no
# kink anywhere. `_selfcheck` asserts both.
DRIFT_LEN = 57.0       # 1.5 wavelengths, about as wide as the traced strand
DRIFT_WAVE = 38.0
DRIFT_AMP = 4.5        # shallow enough that bases on the inside never cross:
                       # the radius at a crest is (WAVE/2pi)^2 / AMP = 8.1
DRIFT_TILT = radians(-12)   # rising left to right, as the panel draws it
DRIFT_BASE = 4.00    # an unpaired base, the length of a pair in the stem
DRIFT_STEP = 3.40    # along the strand, between bases
DRIFT_END = 2.50     # bare strand left at each end, where it was cut


def wave(length, wavelength, amp, tilt):
    """A sine as cubics, knotted every eighth wavelength and then tilted.

    A knot every quarter is the obvious choice and is wrong by about 1% of the
    amplitude; every eighth divides that by 32, which `_selfcheck` holds it to.
    """
    h = wavelength / 8
    n = int(round(length / h))
    k = 2 * pi / wavelength
    rot = lambda q: (q[0] * cos(tilt) - q[1] * sin(tilt),
                     q[0] * sin(tilt) + q[1] * cos(tilt))
    pt = lambda x: (x, amp * sin(k * x))
    slope = lambda x: amp * k * cos(k * x)
    segs = []
    for i in range(n):
        x0, x1 = i * h, (i + 1) * h
        segs.append(tuple(rot(q) for q in (
            pt(x0),
            (x0 + h / 3, amp * sin(k * x0) + h / 3 * slope(x0)),
            (x1 - h / 3, amp * sin(k * x1) - h / 3 * slope(x1)),
            pt(x1))))
    return segs


_DRIFT = wave(DRIFT_LEN, DRIFT_WAVE, DRIFT_AMP, DRIFT_TILT)

# attachment points, from the radii above
_ARM_DX = sqrt(SIDE_R ** 2 - ARM_HALF ** 2)      # 6.917
_MID_DX = sqrt(MID_R ** 2 - (ARM_Y - ARM_HALF - MID_CY) ** 2)   # 5.809
_SP_CY = SP_TOP - sqrt(SP_R ** 2 - STEM_HALF ** 2)


def _path(d, col, w=None):
    return (f'<path d="{d}" fill="none" stroke="{col}" '
            f'stroke-width="{w or SW}" stroke-linecap="round" '
            f'stroke-linejoin="round"/>')


def _arc(r, large, sweep, p):
    return f"A {r} {r} 0 {large} {sweep} {p[0]:.3f},{p[1]:.3f}"


def _tick(x0, y0, x1, y1):
    """A base pair, rail centre to rail centre.

    Butt-ended, and drawn before the strands, so the rails cover both ends: a
    round cap here pokes out past the backbone, which the reference never does.
    """
    return (f'<path d="M {x0:.3f},{y0:.3f} L {x1:.3f},{y1:.3f}" fill="none" '
            f'stroke="{TICK}" stroke-width="{RUNG}" stroke-linecap="butt"/>')


def _unspliced():
    top, bot = ARM_Y - ARM_HALF, ARM_Y + ARM_HALF
    far = SIDE_CX + SIDE_R                       # the outermost point of a loop
    o = []
    for s in (-1, 1):
        for i in range(ARM_RUNGS):
            # an arm runs inward from its loop, so the offsets come off the
            # attachment point toward the junction, not away from it
            x = s * (SIDE_CX - _ARM_DX - ARM_RUNG_0 - i * ARM_RUNG_D)
            o.append(_tick(x, top, x, bot))
    for i in range(STEM_RUNGS):
        y = STEM_RUNG_0 + i * STEM_RUNG_D
        o.append(_tick(-STEM_HALF, y, STEM_HALF, y))
    # the exons: each runs from the far side of its end loop, round the bottom
    # of it, in along the arm's lower rail, and down the stem
    for s in (-1, 1):
        o.append(_path(
            f"M {s*far:.3f},{ARM_Y} "
            + _arc(SIDE_R, 0, 0 if s < 0 else 1, (s * (SIDE_CX - _ARM_DX), bot))
            + f" L {s*(STEM_HALF+CORNER):.3f},{bot:.3f} "
            + _arc(CORNER, 0, 1 if s < 0 else 0, (s * STEM_HALF, STEM_TOP))
            + f" L {s*STEM_HALF:.3f},{STEM_BOT:.3f}", EXON))
    # the intron: one strand over the top, loop to loop, across both arms and
    # the junction loop between them
    o.append(_path(
        f"M {-far:.3f},{ARM_Y} "
        + _arc(SIDE_R, 0, 1, (-(SIDE_CX - _ARM_DX), top))
        + f" L {-_MID_DX:.3f},{top:.3f} "
        + _arc(MID_R, 1, 1, (_MID_DX, top))
        + f" L {SIDE_CX - _ARM_DX:.3f},{top:.3f} "
        + _arc(SIDE_R, 0, 1, (far, ARM_Y)), INTRON))
    return o, (-far, MID_CY - MID_R, far, STEM_BOT)


def spline(segs):
    """The cubics as one path. Public: `receptor.py` draws its linker with it."""
    d = f"M {segs[0][0][0]:.3f},{segs[0][0][1]:.3f}"
    for _, c1, c2, p in segs:
        d += (f" C {c1[0]:.3f},{c1[1]:.3f} {c2[0]:.3f},{c2[1]:.3f} "
              f"{p[0]:.3f},{p[1]:.3f}")
    return d


def _walk(segs, n=60):
    """Dense samples along the same cubics, with the arc length at each."""
    out, total, prev = [], 0.0, None
    for p0, c1, c2, p1 in segs:
        for i in range(n + 1):
            t, u = i / n, 1 - i / n
            q = (u**3 * p0[0] + 3*u*u*t * c1[0] + 3*u*t*t * c2[0] + t**3 * p1[0],
                 u**3 * p0[1] + 3*u*u*t * c1[1] + 3*u*t*t * c2[1] + t**3 * p1[1])
            if prev is not None:
                if q == prev:
                    continue
                total += ((q[0]-prev[0])**2 + (q[1]-prev[1])**2) ** 0.5
            out.append((total, q))
            prev = q
    return out


def _intron():
    """The piece that is cut out, on its own and no longer paired."""
    walk = _walk(_DRIFT)
    o, pts = [], []
    at = DRIFT_END
    for (s0, a), (s1, b) in zip(walk, walk[1:]):
        if not s0 <= at < s1:
            continue
        while s0 <= at < s1 and at <= walk[-1][0] - DRIFT_END:
            tx, ty = b[0] - a[0], b[1] - a[1]
            n = (tx * tx + ty * ty) ** 0.5
            # the lower side of a strand that runs left to right
            nx, ny = -ty / n, tx / n
            q = (a[0] + nx * DRIFT_BASE, a[1] + ny * DRIFT_BASE)
            o.append(_tick(a[0], a[1], q[0], q[1]))
            pts.append(q)
            at += DRIFT_STEP
    o.append(_path(spline(_DRIFT), INTRON))
    pts += [q for _, q in walk]
    xs, ys = [q[0] for q in pts], [q[1] for q in pts]
    return o, (min(xs), min(ys), max(xs), max(ys))


def _spliced():
    o = [_tick(-STEM_HALF, SP_RUNG_0 + i * SP_RUNG_D,
               STEM_HALF, SP_RUNG_0 + i * SP_RUNG_D) for i in range(SP_RUNGS)]
    o.append(_path(f"M {-STEM_HALF:.3f},{SP_BOT:.3f} L {-STEM_HALF:.3f},{SP_TOP:.3f} "
                   + _arc(SP_R, 1, 1, (STEM_HALF, SP_TOP))
                   + f" L {STEM_HALF:.3f},{SP_BOT:.3f}", EXON))
    return o, (-SP_R, _SP_CY - SP_R, SP_R, SP_BOT)


def build(spliced=False, intron=False, scale=SCALE):
    """The hairpin before (default) or after the cut, or the cut-out piece."""
    assert not (spliced and intron), "intron is a third form, not a spliced one"
    body, (x0, y0, x1, y1) = (_intron() if intron else
                              _spliced() if spliced else _unspliced())
    m = SW
    x0, y0, w, h = x0 - m, y0 - m, (x1 - x0) + 2 * m, (y1 - y0) + 2 * m
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w*scale:.0f}" '
            f'height="{h*scale:.0f}" viewBox="0 0 {w*scale:.1f} {h*scale:.1f}">'
            f'<g transform="scale({scale}) translate({-x0:.3f},{-y0:.3f})">'
            + "".join(body) + "</g></svg>")


VARIANTS = [("hairpin_at_unspliced", dict(spliced=False)),
            ("hairpin_at_spliced", dict(spliced=True)),
            ("hairpin_at_intron", dict(intron=True))]


def _selfcheck():
    """The shape the measurements describe, re-derived and compared."""
    # the arms meet their end loops on the loop, not beside it
    for dy in (-ARM_HALF, ARM_HALF):
        r = sqrt(_ARM_DX ** 2 + dy ** 2)
        assert abs(r - SIDE_R) < 1e-9, f"arm leaves the end loop at r={r}"
    # the junction loop is met by the arms' upper rails, likewise
    r = sqrt(_MID_DX ** 2 + (ARM_Y - ARM_HALF - MID_CY) ** 2)
    assert abs(r - MID_R) < 1e-9, f"arm leaves the junction loop at r={r}"
    # the junction loop clears the arms it sits on, or it is not a loop
    assert MID_CY - MID_R < ARM_Y - ARM_HALF, "junction loop is inside the arms"
    # the corner is a quarter turn: it must fit between rail and stem
    assert STEM_HALF + CORNER < SIDE_CX - _ARM_DX, "the corners overlap the loops"
    # the stem is the junction gap, so it is wider than an arm - as drawn
    assert STEM_HALF > ARM_HALF, "stem narrower than an arm; check the trace"
    # rungs stay inside the duplex they cross
    assert ARM_RUNG_0 + (ARM_RUNGS - 1) * ARM_RUNG_D < \
        (SIDE_CX - _ARM_DX) - (STEM_HALF + CORNER), "an arm rung is past the junction"
    assert ARM_RUNG_0 > 0, "an arm rung is inside its end loop"
    assert STEM_RUNG_0 > STEM_TOP and \
        STEM_RUNG_0 + (STEM_RUNGS - 1) * STEM_RUNG_D < STEM_BOT, "stem rung outside"
    assert SP_RUNG_0 > SP_TOP and \
        SP_RUNG_0 + (SP_RUNGS - 1) * SP_RUNG_D < SP_BOT, "spliced rung outside"
    # the two drawings share a helix width, so they read as one molecule
    assert "stroke-width" in build() and "stroke-width" in build(True)
    # the loose intron: bases spaced along it, none running off either cut end
    n = build(intron=True).count(f'stroke="{TICK}"')
    length = _walk(_DRIFT)[-1][0]
    assert n == int((length - 2 * DRIFT_END) // DRIFT_STEP) + 1, \
        f"{n} bases on a strand {length:.1f} long"
    assert 14 <= n <= 22, f"{n} bases on the loose intron; the panel draws ~18"
    # the wave is a sine to within a tenth of a percent of its amplitude
    worst = 0.0
    for i, (p0, c1, c2, p1) in enumerate(_DRIFT):
        back = lambda q: (q[0] * cos(-DRIFT_TILT) - q[1] * sin(-DRIFT_TILT),
                          q[0] * sin(-DRIFT_TILT) + q[1] * cos(-DRIFT_TILT))
        a, b, c, d = (back(q) for q in (p0, c1, c2, p1))
        for j in range(21):
            t, u = j / 20, 1 - j / 20
            x = u**3*a[0] + 3*u*u*t*b[0] + 3*u*t*t*c[0] + t**3*d[0]
            y = u**3*a[1] + 3*u*u*t*b[1] + 3*u*t*t*c[1] + t**3*d[1]
            worst = max(worst, abs(y - DRIFT_AMP * sin(2 * pi * x / DRIFT_WAVE)))
    assert worst < DRIFT_AMP * 1e-3, f"the wave is off a sine by {worst:.4f}"
    # bases on the inside of a crest must not meet: radius there beats length
    assert (DRIFT_WAVE / (2 * pi)) ** 2 / DRIFT_AMP > DRIFT_BASE, \
        "the bases cross on the inside of a crest; flatten the wave"
    # no kink at a joint: the handles either side of a knot are collinear with
    # it, which is the whole reason the strand is written as cubics by hand
    for (_, _, c2, k), (k2, c3, _, _) in zip(_DRIFT, _DRIFT[1:]):
        assert k == k2, f"the strand is broken at {k}"
        cross = ((k[0]-c2[0]) * (c3[1]-k[1])) - ((k[1]-c2[1]) * (c3[0]-k[0]))
        assert abs(cross) < 1e-9, f"kink at {k}: handles are not collinear"
    # proportion: the reference is wide and flat, ~2.2 : 1
    import re
    w, h = (float(v) for v in
            re.search(r'width="(\d+)" height="(\d+)"', build()).groups())
    assert 1.9 < w / h < 2.5, f"unspliced is {w/h:.2f} : 1, the reference is ~2.2"
    print(f"hairpin: checks pass, unspliced {w:.0f} x {h:.0f} ({w/h:.2f} : 1)")


def main(outdir=None):
    from common import paths
    lib = paths.fig_output("cartoon_library")
    svg_dir = Path(outdir) if outdir else lib / "svg"
    svg_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, kw in VARIANTS:
        p = svg_dir / f"{name}.svg"
        p.write_text(build(**kw), encoding="utf-8")
        written.append(p)
    if outdir is None:
        png_dir = lib / "png"
        png_dir.mkdir(parents=True, exist_ok=True)
        ql = shutil.which("qlmanage")
        if ql:
            subprocess.run([ql, "-t", "-s", "1000", "-o", str(png_dir),
                            *map(str, written)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for q in list(png_dir.glob("*.svg.png")):
                if q.name.startswith("._"):
                    q.unlink(missing_ok=True)
                    continue
                # rename fails with EBUSY often enough on the SMB share that a
                # copy is the only reliable move here
                target = q.with_name(q.name[:-len(".svg.png")] + ".png")
                target.write_bytes(q.read_bytes())
                q.unlink(missing_ok=True)
    print(f"wrote {len(written)} to {svg_dir}")


if __name__ == "__main__":
    if "--check" in sys.argv:
        _selfcheck()
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else None)
