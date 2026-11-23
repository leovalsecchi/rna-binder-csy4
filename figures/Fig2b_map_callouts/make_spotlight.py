#!/usr/bin/env python3
"""Spotlight the interaction inside the full zoom without losing the context.

Takes two renders from the IDENTICAL camera:
  full      the normal zoom, faded context and all
  paironly  the same view with only the stacking pair drawn

and composites `paironly` into a disc centred on the interaction, so inside the
disc the structures in front of and behind the pair disappear, while outside
the disc the context is untouched. The edge is feathered so it reads as a lens
rather than a cut.

The disc centre is MEASURED, not guessed: it is the alpha centroid of the
pair-only render, which contains nothing but the interaction.

Usage: python make_spotlight.py <full> <paironly> <out> [--radius 0.20]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("full"); ap.add_argument("paironly"); ap.add_argument("out")
    ap.add_argument("--radius", type=float, default=0.20,
                    help="disc radius as a fraction of image width")
    ap.add_argument("--feather", type=float, default=0.12,
                    help="soft edge as a fraction of the radius")
    ap.add_argument("--keep", type=float, default=0.0,
                    help="0 removes the context inside the disc, 1 keeps it")
    ap.add_argument("--ring", type=int, default=0, help="ring width in px, 0 for none")
    a = ap.parse_args()

    full = Image.open(a.full).convert("RGBA")
    pair = Image.open(a.paironly).convert("RGBA")
    if full.size != pair.size:
        sys.exit(f"camera mismatch: {full.size} vs {pair.size}")
    w, h = full.size

    al = np.array(pair.getchannel("A"), dtype=float)
    if al.sum() == 0:
        sys.exit("pair-only render is empty")
    ys, xs = np.nonzero(al > 8)
    cx, cy = int(xs.mean()), int(ys.mean())
    r = int(a.radius * w)

    disc = Image.new("L", (w, h), 0)
    ImageDraw.Draw(disc).ellipse((cx - r, cy - r, cx + r, cy + r), fill=255)
    disc = disc.filter(ImageFilter.GaussianBlur(max(1.0, a.feather * r)))

    # inside the disc: the context is dropped to `keep`, the pair is drawn back
    faded = full.copy()
    if a.keep < 1.0:
        fa = np.array(faded.getchannel("A"), dtype=float)
        d = np.array(disc, dtype=float) / 255.0
        fa *= (1.0 - d * (1.0 - a.keep))
        faded.putalpha(Image.fromarray(fa.astype(np.uint8)))

    out = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    out.alpha_composite(faded)
    inside = pair.copy()
    ia = (np.array(inside.getchannel("A"), dtype=float)
          * (np.array(disc, dtype=float) / 255.0))
    inside.putalpha(Image.fromarray(ia.astype(np.uint8)))
    out.alpha_composite(inside)

    if a.ring:
        ring = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        ImageDraw.Draw(ring).ellipse((cx - r, cy - r, cx + r, cy + r),
                                     outline=palette.NEUTRAL["500"], width=a.ring)
        out.alpha_composite(ring)

    out.save(a.out)
    print(f"wrote {a.out}  centre=({cx},{cy}) radius={r}px of {w}")


if __name__ == "__main__":
    main()
