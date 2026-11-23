"""Shared matplotlib style. Sizes from config/style.yaml, colours from config/palette.json."""
from functools import lru_cache as _lru_cache

import yaml
import matplotlib as mpl
from cycler import cycler
from .paths import ROOT
from . import palette

CFG = yaml.safe_load((ROOT / "config" / "style.yaml").read_text())
COLORS = palette.STACKS   # kept for backwards compatibility; prefer palette.* directly


def mm(x: float) -> float:
    """mm to inches."""
    return x / 25.4


def apply():
    mpl.rcParams.update({
        "font.family": CFG["font_family"],
        "font.size": CFG["font_size_pt"],
        "axes.titlesize": CFG["title_size_pt"],
        "axes.linewidth": CFG["line_width_pt"],
        "savefig.dpi": CFG["dpi"],
        "pdf.fonttype": 42,   # editable text in Illustrator/Inkscape
        "svg.fonttype": "none",
        # 2026-11-15 (Leo, font audit): matplotlib's mathtext ignores
        # font.family and ships its own DejaVu Sans, so every $..$ italic
        # (the n, P and q in Fig 5 e/f) came out DejaVu in the SVG and
        # HELVETICA-OBLIQUE in the PDF - the only non-Arial type in the book.
        # "custom" points mathtext at the house family instead, so an italic
        # n is Arial Italic like every other letter on the page.
        "mathtext.fontset": "custom",
        "mathtext.rm": CFG["font_family"],
        "mathtext.it": f"{CFG['font_family']}:italic",
        "mathtext.bf": f"{CFG['font_family']}:bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": cycler(color=[palette.STACKS[k] for k in palette.STACK_ORDER]),
        "text.color": palette.INK,
        "axes.labelcolor": palette.INK,
        "axes.edgecolor": palette.INK,
        "xtick.color": palette.INK,
        "ytick.color": palette.INK,
        # Gridlines OFF by default, stated rather than inherited from
        # matplotlib: STYLE_NOTES.md section 7 allows them only where a
        # value must be read off precisely. A panel that needs them turns
        # them on locally; grid.color is kept for that case.
        "axes.grid": False,
        "grid.color": palette.GRID,
        "image.cmap": "thesis_seq",
    })
    palette.register_colormaps()


def save(fig, outdir, stem, tight=True, svg=True):
    """Save PDF, PNG and (by default) SVG.

    tight=False emits the AUTHORED canvas exactly, for a panel laid out in mm,
    where a tight box would quietly resize what assembly places.

    SVG is written for page assembly, the same way the main-text figures do it
    (added 2026-11-19). It is the format assembly should consume: `svg.fonttype:
    none` is set in apply(), so labels stay live text an editor can restyle,
    where a PDF converted to SVG afterwards arrives as glyph outlines and cannot
    be edited or re-fonted.
    """
    kw = dict(bbox_inches="tight") if tight else {}
    fig.savefig(outdir / f"{stem}.pdf", **kw)
    fig.savefig(outdir / f"{stem}.png", **kw)
    if svg:
        fig.savefig(outdir / f"{stem}.svg", **kw)


# --- measuring type --------------------------------------------------------
# A laid-out row spaces itself off its labels' widths, so those widths have to
# be MEASURED. The usual stand-in, 0.52 em per character, is out by 10-33%
# here and by a different amount for every string - "filtering" by a third,
# "Reverse" by one percent - so the slack it leaves is different in every gap.
# That does not show up as a wrong label; it shows up as a row whose spacing
# is subtly uneven for no visible reason.
_FONT_CHECKED = False


@_lru_cache(maxsize=None)
def text_w(s, size=None, weight="normal"):
    """Ink width of `s` set in the house font, in the SAME unit as `size`.

    Width scales linearly with the font size, so passing a size in points
    gives points and passing one in viewBox units gives units.

    `weight` matters when a label is emphasised: Arial Bold is a few percent
    wider than Arial, which is enough to push a measured-to-fit label over the
    edge of its box. Default "normal", so existing callers are unaffected.
    """
    global _FONT_CHECKED
    from matplotlib.font_manager import FontProperties, findfont
    from matplotlib.textpath import TextPath
    size = CFG["font_size_pt"] if size is None else size
    fp = FontProperties(family=CFG["font_family"], size=size, weight=weight)
    if not _FONT_CHECKED:
        got = findfont(fp)
        assert CFG["font_family"].lower() in got.lower(), (
            f"{CFG['font_family']} not found - matplotlib would measure "
            f"{got}, and every gap laid out from these widths would be wrong")
        _FONT_CHECKED = True
    return TextPath((0, 0), s, prop=fp).get_extents().width if s else 0.0


# --- 2-D density panels ----------------------------------------------------
# The count ramp and the count SCALE for every 2-D density panel in the book,
# shared rather than copied: the main text's Fig3a_backend_agreement and the
# supplementary repredictor_comparison overlap on the AF3-ESMFold2 column, and
# two copies of this would eventually disagree there.
#
# Settled 2026-11-19 (Leo, over four rounds): the house blue is BACK ("if its
# possible we do the original blue"). The detour was THEME 200-800 -> YlGnBu ->
# grey/blue -> all grey -> here. YlGnBu put yellow and green in a figure set
# that has neither; all-grey had to stop at NEUTRAL 700 or the INK threshold
# rules vanished into the dense core, which cost it most of its range.
_DENSITY_STEPS = ("200", "300", "400", "500", "600", "700", "800")


def _lab(h):
    """sRGB hex -> CIELAB (D65). L* is the lightness a reader ranks the ramp by."""
    f = lambda u: (u / 12.92 if u <= .04045 else ((u + .055) / 1.055) ** 2.4)
    r, g, b = [f(int(h.lstrip("#")[i:i + 2], 16) / 255) for i in (0, 2, 4)]
    xyz = [(.4124 * r + .3576 * g + .1805 * b) / .95047,
           (.2126 * r + .7152 * g + .0722 * b),
           (.0193 * r + .1192 * g + .9505 * b) / 1.08883]
    q = [t ** (1 / 3) if t > .008856 else 7.787 * t + 16 / 116 for t in xyz]
    return [116 * q[1] - 16, 500 * (q[0] - q[1]), 200 * (q[1] - q[2])]


def _unlab(lab):
    """CIELAB -> sRGB floats, clipped to gamut (this path needs <0.001 of clip)."""
    fy = (lab[0] + 16) / 116
    fx, fz = fy + lab[1] / 500, fy - lab[2] / 200
    g = lambda u: u ** 3 if u ** 3 > .008856 else (u - 16 / 116) / 7.787
    x, y, z = g(fx) * .95047, g(fy), g(fz) * 1.08883
    lin = (3.2406 * x - 1.5372 * y - .4986 * z,
           -.9689 * x + 1.8758 * y + .0415 * z,
           .0557 * x - .2040 * y + 1.0570 * z)
    f = lambda u: 12.92 * u if u <= .0031308 else 1.055 * u ** (1 / 2.4) - .055
    return tuple(f(min(1.0, max(0.0, c))) for c in lin)


def density_cmap(name="house_blue", n=256):
    """The house theme 200 -> 800, resampled to CONSTANT lightness steps.

    Perceptually uniform in the sense the matplotlib and Kovesi write-ups use:
    L* falls linearly along the colormap, so equal steps in count read as equal
    steps in darkness. This is a refinement, not a rescue - measured in float,
    `from_list`'s sRGB interpolation of these seven stops was already close
    (dL* per colormap step 0.212-0.255, ratio 1.21), because the stops
    themselves are even in lightness (dL* 9.0, then 10.4 five times over).
    Walking the same stops in CIELAB and resampling at constant dL* makes it
    exact: 0.240 per step, L* 87.5 -> 26.3, ratio 1.00. Hue is near-constant
    along the path (245-252 deg) and chroma is a single smooth bulge (15 -> 34
    -> 23, the gamut), so there was no kink to take out and the colour scheme is
    untouched. dE76 per step spans 0.240-0.313 from that chroma bulge, the
    residual non-uniformity plain CIELAB carries; equalising dE instead would
    buy that back at the cost of linear lightness, which is the axis a reader
    actually ranks. Do not measure either ramp through 8-bit hex: 0.24 L* is
    finer than one sRGB code, so rounding alone fakes dL* of 0.00-0.42.
    Empty bins transparent.
    """
    import numpy as np
    from matplotlib.colors import ListedColormap
    lab = np.array([_lab(palette.THEME[k]) for k in _DENSITY_STEPS])
    u = np.linspace(0, 1, 2048)
    path = np.stack([np.interp(u, np.linspace(0, 1, len(lab)), lab[:, j])
                     for j in range(3)], 1)
    assert np.all(np.diff(path[:, 0]) < 0), \
        f"ramp is not monotonically darker: {[round(v, 1) for v in lab[:, 0]]}"
    inc = path[::-1]                        # L* ascending, as np.interp needs
    L = np.linspace(path[0, 0], path[-1, 0], n)
    ramp = np.stack([L] + [np.interp(L, inc[:, 0], inc[:, j]) for j in (1, 2)], 1)
    de = np.linalg.norm(np.diff(ramp, axis=0), axis=1)
    assert de.max() / de.min() < 1.4, \
        f"perceptual steps uneven: dE76 {de.min():.3f}-{de.max():.3f}"
    cm = ListedColormap([_unlab(c) for c in ramp], name)
    cm.set_bad("none")
    return cm


def density_norm(hb, lo_pct=25, hi_pct=99.5):
    """Put `hb` on a log count scale spanning where the data actually are.

    LogNorm(1, max) was the wrong scale twice over, measured on Fig3a's panels:
    (i) the bottom half of occupied bins holds ~1-4% of the points yet took the
    lower HALF of the colours, so the body of the cloud was squeezed into the
    top third and barely varied; (ii) each panel's max differs by up to 48x, so
    panels looked darker or lighter than each other for no reason a reader can
    see. Spanning p25 to p99.5 of the OCCUPIED bins fixes both: below p25
    saturates pale, above p99.5 saturates dark (the one-corner spike, a handful
    of bins), everything between gets the full range. Still per-panel - a shared
    scale would wash the ipTM panel out entirely at DockQ's 770k maximum - so
    the caption has to say each panel is scaled to its own counts.

    Returns (lo, hi, position_of), the last a callable giving a count's place on
    the ramp in 0..1, for the diagnostic line both scripts print.
    """
    import numpy as np
    from matplotlib.colors import LogNorm
    c = hb.get_array()
    c = np.asarray(c[~np.ma.getmaskarray(c)] if np.ma.isMaskedArray(c) else c, float)
    lo, hi = max(1.0, np.percentile(c, lo_pct)), np.percentile(c, hi_pct)
    hb.set_norm(LogNorm(lo, hi))
    return lo, hi, (lambda v: min(1.0, max(0.0, np.log10(max(v, lo) / lo)
                                           / np.log10(hi / lo))))


if __name__ == "__main__":        # the smallest check that fails if the ramp breaks
    import numpy as np
    L = lambda rgb: _lab("#%02x%02x%02x" % tuple(round(255 * c) for c in rgb))[0]
    # Round-trip: _unlab(_lab(x)) == x for every house stop (8-bit exact).
    for k, h in palette.THEME.items():
        assert "#%02x%02x%02x" % tuple(round(255 * c) for c in _unlab(_lab(h))) == h, (k, h)
    # Lightness is linear along the colormap. Measured in FLOAT, not through hex:
    # 0.24 L* per step is finer than one 8-bit code, so rounding fakes the answer.
    cm = density_cmap()
    f = lambda u: u / 12.92 if u <= .04045 else ((u + .055) / 1.055) ** 2.4
    lf = lambda rgb: 116 * sum(w * f(c) for w, c in zip((.2126, .7152, .0722), rgb)) ** (1 / 3) - 16
    d = -np.diff([lf(cm(i / 255)[:3]) for i in range(256)])
    assert d.max() / d.min() < 1.001, f"not lightness-linear: dL* {d.min():.4f}-{d.max():.4f}"
    print(f"ok: round-trip exact on {len(palette.THEME)} stops; "
          f"ramp dL*/step {d.mean():.3f} (ratio {d.max() / d.min():.4f})")
