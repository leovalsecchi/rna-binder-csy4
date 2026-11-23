"""Assemble Figure 5 panel c: a wide four-up strip of the renderings.

The panel carries only the twist number and the flow fold change. The stack,
the model used and the chemistry of each contact live in README.md, not on the
figure. Dashes are drawn on two panels only; see render.py for which and why.
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
W, H = 170 * MM, 48 * MM          # standalone default, inches
FONT = 7.5                        # standalone default, pt
INK, SUB = palette.NEUTRAL["900"], palette.NEUTRAL["600"]

PANELS = [
    ("23", "twist 23", "FC 1.56"),
    ("15", "twist 15", "FC 1.54"),
    ("16", "twist 16", "FC 1.47"),
    ("07", "twist 7", "FC 1.41"),
]


def trim(a, pad=8):
    """Crop to the rendered ink. Returns the crop only; padding to a shared box
    is done afterwards by pad_common, once all four crops are known."""
    alpha = a[..., 3] if a.shape[-1] == 4 else np.ones(a.shape[:2])
    ys, xs = np.where(alpha > 0.02)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, a.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, a.shape[1])
    return a[y0:y1, x0:x1]


def pad_common(crops):
    """Pad every crop to one box, the tallest by the widest.

    All four keep a single scale, so sizes stay comparable between panels, but
    the box is only as tall as the tallest structure needs. Padding each to a
    SQUARE instead (the previous version) baked transparent bands into every
    wide crop, which showed up as a gap between the structures and their labels
    once the panel was stretched into the page's short, wide slot.
    """
    h = max(c.shape[0] for c in crops)
    w = max(c.shape[1] for c in crops)
    out = []
    for c in crops:
        o = np.zeros((h, w, c.shape[2]), dtype=c.dtype)
        oy, ox = (h - c.shape[0]) // 2, (w - c.shape[1]) // 2
        o[oy:oy + c.shape[0], ox:ox + c.shape[1]] = c
        out.append(o)
    return out


def main(in_dir, out_dir, W=W, H=H, FONT=FONT):
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # 2026-11-16 (Leo, font audit): this listed Helvetica FIRST, so every label
    # in this panel - "twist 23", "FC 1.56" - resolved to Helvetica and
    # Helvetica-Bold, and they were the last non-Arial faces in the book. The
    # family now comes from config/style.yaml, like every other panel, instead
    # of a hand-written stack that quietly outranked it.
    _house = next(l.split(":", 1)[1].strip()
                  for l in (ROOT / "config" / "style.yaml").read_text().splitlines()
                  if l.startswith("font_family:"))
    mpl.rcParams.update({"font.family": _house,
                         "svg.fonttype": "none", "pdf.fonttype": 42})

    fig = plt.figure(figsize=(W, H))
    left, right, gap = 0.012, 0.012, 0.012
    pw = (1 - left - right - 3 * gap) / 4

    imgs = pad_common([trim(mpimg.imread(in_dir / f"t{tw}.png"))
                       for tw, _, _ in PANELS])

    # Height the image box NEEDS at this width, from the shared box's own
    # aspect. Without this the axes aspect and the image aspect disagree,
    # imshow fits by width, and the leftover height shows as a gap between the
    # structures and their labels, which is glaring once the page stretches
    # this panel into a short, wide slot.
    ar = imgs[0].shape[1] / imgs[0].shape[0]
    # Anchor the content to the BOTTOM: labels sit in lab_space, structures
    # directly above them, and whatever height is left over falls above the
    # structures, where the page's row gap absorbs it. Anchoring to the top
    # instead left the slack below the labels as a visible dead band.
    lab_space, top_pad = 0.09, 0.02
    img_h = min((pw * W / ar) / H, 1 - top_pad - lab_space)
    img_b = lab_space

    for i, (tw, name, fc) in enumerate(PANELS):
        x = left + i * (pw + gap)
        ax = fig.add_axes([x, img_b, pw, img_h])
        # NOT lanczos. Its kernel has negative lobes, and matplotlib
        # resamples RGBA channel by channel without premultiplying, so at a
        # silhouette the ringing drives alpha below zero and the edge comes
        # out stippled with fully transparent pixels. Same choice as
        # Fig5c_interface_chemistry.
        ax.imshow(imgs[i], interpolation="antialiased")
        ax.set_axis_off()
        # Name and FC belong to each other, so set them as one centred pair
        # under their own structure. Anchoring them to opposite edges of the
        # cell (the previous version) pushed each FC up against the NEXT
        # panel's name, which read as if it labelled that one.
        cx, half_gap = x + pw / 2, 0.005
        fig.text(cx - half_gap, img_b - 0.03, name, color=INK, fontsize=FONT,
                 fontweight="bold", ha="right", va="top")
        fig.text(cx + half_gap, img_b - 0.03, fc, color=INK, fontsize=FONT,
                 ha="left", va="top")

    for ext in ("pdf", "svg", "png"):
        fig.savefig(out_dir / f"fig5c_best_hits.{ext}", dpi=600,
                    transparent=(ext != "png"),
                    facecolor="white" if ext == "png" else "none")
    print("wrote", out_dir / "fig5c_best_hits.pdf")


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
        w_mm, h_mm, font_pt, outdir = float(a[0]), float(a[1]), float(a[2]), a[3]
        main(ROOT / "output" / "Fig5c_best_hits", outdir,
             W=w_mm * MM, H=h_mm * MM, FONT=font_pt)
    else:
        main(a[0], a[1])
