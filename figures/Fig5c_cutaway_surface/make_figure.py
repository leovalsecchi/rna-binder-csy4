"""Assemble Figure 5 panel c: the four best hits as surfaces, with an X-ray
window cut into each one showing the cartoon through it.

This WAS panel d, a second row under the cartoon close-ups of panel c. Leo,
2026-11-18: "replace the renderings in c by the ones in d" - so it is panel c
now, labels and all, and the row it used to sit under is off the page.
code/Fig5c_best_hits still builds; it is simply no longer promoted.

Two renderings per design, both from render.py and both through THE SAME
CAMERA - w<twist>.png, the solid surfaces, and x<twist>.png, panel c's
see-through cartoon - so they register pixel for pixel and the window is a
straight alpha cut, not an alignment problem. See the X-ray window block below.

ONE BOX SIZE for all four, so the four stay at one scale and one baseline, but
each structure is centred SIDEWAYS in its own copy of it. See lay_out().
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import image as mpimg

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
from common import palette  # noqa: E402

MM = 1 / 25.4
W, H = 178 * MM, 55 * MM          # standalone default, inches: panel c's slot
FONT = 7.5                        # standalone default, pt
INK = palette.NEUTRAL["900"]

PANELS = [("23", "twist 23", "FC 1.56"), ("15", "twist 15", "FC 1.54"),
          ("16", "twist 16", "FC 1.47"), ("07", "twist 7", "FC 1.41")]

# Which set of renderings to assemble. render.py writes both into the same
# output folder: "w" = the whole uncut surfaces in the Fig 3b / Fig 4-grid
# camera (what the panel uses), "c" = the cut-away sections down the hairpin
# axis. One line, so switching the panel between them is one edit, no re-render.
VARIANT = "w"

# On 2026-11-18 this was False: the panel sat directly under Fig5c_best_hits,
# showing the same four designs in the same order, so its labels would have been
# the same four labels twice, two centimetres apart. It IS that panel now, and
# it carries them.
LABELS = True

# 2026-11-18 (Leo): "put figure 5c in the same ... grey background box" as
# Fig 2d. The same tinted card Fig 1d, Fig 2a and Fig 2d sit on: neutral.100,
# with the corner radius Fig 2d derives as 14 units on a 1520-unit viewBox the
# width of the house block - repeated here as a printed pt value so the two
# boxes carry the same radius however this panel is scaled onto the page.
BOX_FILL = palette.NEUTRAL["100"]
W_BLOCK_PT = (210 - 2 * 16) * 72 / 25.4          # 504.567 pt
BOX_RX_PT = 14.0 * W_BLOCK_PT / 1520.0           # 4.65 pt
PAD_PT = 9.0                     # inside the box, top and bottom (Fig 2d's)
LABEL_GAP_PT = 6.0               # image bottom to label top (Fig 2d's)

# --- the X-ray window ------------------------------------------------------
# Leo, 2026-11-18: "you use the surface one as the background, the larger
# background, and then you just crop the figure C style one, have it on top,
# crop it, and then it looks kind of like an X-ray into the figure."
#
# render.py --xray draws panel c's see-through cartoon THROUGH THE SAME CAMERA
# as the surface, so x*.png and w*.png register pixel for pixel and a window
# can simply be cut in the surface with the cartoon showing through it. The
# window is a SQUARE with a hard edge and a hairline rule around it (Leo,
# 2026-11-18: "a square with a clear cut (no blur needed (just thin lines))").
#
# Where the window sits is MEASURED, not chosen: `view` points the camera at the
# centre of what it is framed on - that design's contact residues and the bases
# they read - so the interface is at the centre of the saved frame, and the
# crop offset puts that back in cropped pixels.
XRAY = True
XR_FRAC = 0.24     # HALF the window's side, as a fraction of the panel's
                   # shorter side. A disc of radius 0.42 was the first cut and
                   # showed too much ("the zoom in is too large"); 0.20 was the
                   # answer to that and went a little far ("the squares a bit
                   # larger").
XR_LW = 0.5        # pt, the rule around the window: the book's thinnest
XR_EDGE = "600"    # neutral step for that rule; it marks the cut, it is not
                   # another thing to read
BG_IN = 0.45       # how much of the SURFACE is left inside the window. 0 is a
                   # clean hole; a ghost of the envelope over the cartoon is
                   # what makes it read as an X-ray rather than as a bite taken
                   # out of the surface. It started at 0.25 and went up with
                   # "make the background of the squares less transparent";
                   # 0.65 was rendered too and is where the envelope starts
                   # swallowing the cartoon, worst on twist 23, whose binder
                   # wraps the hairpin. The cartoon was strengthened in the
                   # same pass (render.py XR_T_*), which is what lets the
                   # background carry this much.


def lay_out(raw, pad=10):
    """Where each panel's crop sits. Returns (y0, y1, box width, shifts).

    VERTICALLY the box is shared - the union of the four inks - so the four
    keep one baseline and one scale.

    SIDEWAYS each panel is centred on its OWN ink. The four cells are equally
    spaced, so the structures standing in them have to be too, and one shared
    box does not do that: it holds the four where their cameras left them, so
    the widest binder leaned against its cell's edge while the next cell had a
    hole beside it (Leo, 2026-11-18: "check that they are visually all about
    the same size and spaced correctly"). The box is as wide as the widest ink,
    so nothing is cut, and the narrower panels get transparent margins.

    `shift` is the one number each panel needs: destination x = source x +
    shift, for the surface layer, the cartoon layer and the window alike.
    """
    ink = [np.where(a[..., 3] > 0.02) for a in raw]
    h, w = raw[0].shape[:2]
    y0 = max(min(ys.min() for ys, _ in ink) - pad, 0)
    y1 = min(max(ys.max() for ys, _ in ink) + pad + 1, h)
    bw = max(xs.max() - xs.min() + 1 for _, xs in ink) + 2 * pad
    shifts = [(bw - (xs.max() - xs.min() + 1 + 2 * pad)) // 2 - (xs.min() - pad)
              for _, xs in ink]
    return y0, y1, bw, shifts


def place(a, y0, y1, bw, shift):
    """One panel's crop: rows y0:y1, shifted sideways onto a bw-wide canvas."""
    c = np.zeros((y1 - y0, bw, a.shape[2]), a.dtype)
    src = a[y0:y1]
    j, i = max(shift, 0), max(-shift, 0)
    n = min(bw - j, src.shape[1] - i)
    c[:, j:j + n] = src[:, i:i + n]
    return c


def window(shape, cx, cy, half):
    """1 inside the window, 0 outside. No ramp: the cut is clean and the rule
    drawn on top of it does the work the feather used to do."""
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    return ((np.abs(xx - cx) <= half) & (np.abs(yy - cy) <= half)).astype(float)


def main(in_dir, out_dir, W=W, H=H, FONT=FONT):
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Family from config/style.yaml like every other panel, never a hand-written
    # stack: a local one quietly outranks the house setting and leaks a non-Arial
    # face into the book (font audit, 2026-11-16).
    house = next(l.split(":", 1)[1].strip()
                 for l in (ROOT / "config" / "style.yaml").read_text().splitlines()
                 if l.startswith("font_family:"))
    mpl.rcParams.update({"font.family": house, "svg.fonttype": "none",
                         "pdf.fonttype": 42})

    raw = [mpimg.imread(in_dir / f"{VARIANT}{tw}.png") for tw, _, _ in PANELS]
    full_h, full_w = raw[0].shape[:2]
    y0, y1, bw, shifts = lay_out(raw)
    imgs = [place(a, y0, y1, bw, s) for a, s in zip(raw, shifts)]
    cys = full_h / 2 - y0
    cxs = [full_w / 2 + s for s in shifts]
    if XRAY:
        # The cartoon layer takes the SAME crop and shift, so it stays in
        # register pixel for pixel. The window follows the structure: render.py
        # puts each design's interface at the centre of its frame, so the
        # window is on the interface wherever the structure has been moved to.
        xr = [place(mpimg.imread(in_dir / f"x{tw}.png"), y0, y1, bw, s)
              for (tw, _, _), s in zip(PANELS, shifts)]
        half = XR_FRAC * min(imgs[0].shape[:2])
        for a, cx in zip(imgs, cxs):
            a[..., 3] *= 1.0 - window(a.shape[:2], cx, cys, half) * (1.0 - BG_IN)
        for a, cx in zip(xr, cxs):
            a[..., 3] *= window(a.shape[:2], cx, cys, half)

    fig = plt.figure(figsize=(W, H))
    left = right = gap = 0.012
    pw = (1 - left - right - 3 * gap) / 4
    ar = imgs[0].shape[1] / imgs[0].shape[0]
    # The card is the panel's own ink now, so the layout is driven from its
    # padding rather than from bare fractions: the labels sit PAD_PT above the
    # card's foot, the structures LABEL_GAP_PT above the labels, and the card
    # closes PAD_PT above the structures. Whatever the page's slot has left
    # over falls OUTSIDE the card, where the row gap absorbs it - inside it
    # would print as a dead grey band.
    pad = PAD_PT / 72 / H                              # figure fractions
    lab_h = (FONT * 1.15 / 72 / H) if LABELS else 0.0  # the label's own height
    img_b = pad + (lab_h + LABEL_GAP_PT / 72 / H if LABELS else 0.0)
    img_h = min((pw * W / ar) / H, 1 - pad - img_b)
    fig.add_artist(mpl.patches.FancyBboxPatch(
        (0, 0), W, (img_b + img_h + pad) * H, transform=fig.dpi_scale_trans,
        boxstyle=f"round,pad=0,rounding_size={BOX_RX_PT / 72}",
        facecolor=BOX_FILL, edgecolor="none", zorder=0))

    for i, (tw, name, fc) in enumerate(PANELS):
        x = left + i * (pw + gap)
        ax = fig.add_axes([x, img_b, pw, img_h])
        # NOT lanczos. Its kernel has negative lobes and matplotlib resamples
        # RGBA channel by channel without premultiplying, so at a silhouette the
        # ringing drives alpha below zero and the edge comes out stippled with
        # fully transparent pixels. Same choice as Fig5c_best_hits.
        ax.imshow(imgs[i], interpolation="antialiased")
        if XRAY:
            ax.imshow(xr[i], interpolation="antialiased")
            # Vector, so it stays a hairline in the PDF: matplotlib composites
            # the two image layers into one raster, but not the patches.
            ax.add_patch(mpl.patches.Rectangle(
                (cxs[i] - half, cys - half), 2 * half, 2 * half, fill=False,
                lw=XR_LW, edgecolor=palette.NEUTRAL[XR_EDGE], zorder=5))
        ax.set_axis_off()
        if LABELS:
            # NOT cx/half: those are the window's centre and half-side, in
            # pixels, and rebinding them here moved the next panels' rules onto
            # the label. Name and FC belong to each other, so they are set as
            # one centred pair under their own structure.
            lx, lgap = x + pw / 2, 0.005
            ly = img_b - LABEL_GAP_PT / 72 / H
            fig.text(lx - lgap, ly, name, color=INK, fontsize=FONT,
                     fontweight="bold", ha="right", va="top", zorder=3)
            fig.text(lx + lgap, ly, fc, color=INK, fontsize=FONT,
                     ha="left", va="top", zorder=3)

    for ext in ("pdf", "svg", "png"):
        fig.savefig(out_dir / f"fig5c_cutaway_surface.{ext}", dpi=600,
                    transparent=(ext != "png"),
                    facecolor="white" if ext == "png" else "none")
    print("wrote", out_dir / "fig5c_cutaway_surface.svg")


def _is_num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 4 and all(_is_num(v) for v in a[:3]):
        # House panel CLI used by common.page_assembly.variant:
        #   W_MM H_MM FONT_PT OUTDIR
        # Renderings always come from the standard output dir; only the
        # assembled SVG goes to OUTDIR.
        main(ROOT / "output" / "Fig5c_cutaway_surface", a[3],
             W=float(a[0]) * MM, H=float(a[1]) * MM, FONT=float(a[2]))
    else:
        main(a[0], a[1])
