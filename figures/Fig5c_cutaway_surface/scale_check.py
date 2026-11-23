"""Are the four panels at one scale?

    ~/.venvs/binder_map/bin/python scale_check.py [render dir]

The panels are at one scale when the camera is set ONCE, before the first
panel, and no panel re-frames or re-zooms: a panel's own `move x/y` slides it
sideways, and translation cannot change size. That is what this asserts, on the
command files render.py wrote - it is a property of the camera, not of a
measurement, so it is checked where it is decided.

It also PRINTS the rendered hairpin's height per panel, which is the thing a
reader actually sees - the hairpin is the same molecule in all four, and each
panel's `turn y` spins about the screen's vertical axis, so height is the one
extent a spin leaves alone. Do not read it as a ruler: the binder hides part of
the hairpin, by a different amount in each panel, which on these four moves the
number by about a fifth on its own. It is here to be eyeballed beside the
render, not to pass or fail.
"""
import sys
from pathlib import Path

import numpy as np
from matplotlib import image as mpimg

TW = ("23", "15", "16", "07")


def one_camera(cxc):
    """No `view` or `zoom` from the first panel onwards, and exactly one of
    each before it."""
    lines = cxc.read_text().splitlines()
    cut = next(i for i, l in enumerate(lines) if l.startswith("save "))
    cam = [l for l in lines[:cut] if l.split()[:1] and l.split()[0] in ("view", "zoom")]
    late = [l for l in lines[cut:] if l.split()[:1] and l.split()[0] in ("view", "zoom")]
    assert not late, f"{cxc.name}: a panel re-frames the camera: {late[:3]}"
    assert len(cam) == 2, f"{cxc.name}: expected one view + one zoom, got {cam}"
    print(f"  {cxc.name}: one camera for all four - {' | '.join(cam)}")


def hairpin_height(png):
    """Vertical extent, in px, of the grey cartoon+nucleotide fill.

    Grey by saturation and LIGHT: the silhouettes are grey too (thesis_grey_700
    against the hairpin's 400) and they outline the binder as well, which would
    report the whole complex instead of the hairpin.
    """
    a = mpimg.imread(png)
    rgb, al = a[..., :3], a[..., 3]
    m = rgb.mean(2)
    g = (rgb.max(2) - rgb.min(2) < 0.06) & (al > 0.5) & (m > 0.50) & (m < 0.95)
    return float(np.ptp(np.where(g)[0]) + 1)


def main(d):
    d = Path(d)
    for name in ("render_whole.cxc", "render_xray.cxc"):
        one_camera(d / name)
    print("visible hairpin (occluded by the binder, see the docstring):")
    h = [hairpin_height(d / f"x{t}.png") for t in TW]
    for t, v in zip(TW, h):
        print(f"  twist {t:>2}: {v:6.0f} px")
    print("OK: one camera, so one scale across the row")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         Path(__file__).resolve().parents[2] / "output" / "Fig5c_cutaway_surface")
