#!/usr/bin/env python3
"""Circular highlight cutout: a round callout showing only the interaction.

Takes the square, pair-centred render and masks it to a disc with everything
outside fully transparent, plus a thin ring in palette neutral 700 so the
callout reads as a lens rather than a cropped photo. Drops straight onto the
fold map next to a leader line.

Usage: python make_cutout.py <square_png> <out_png> [--frac 0.82] [--ring 6]
"""
import argparse
import sys
from pathlib import Path

from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402

SS = 4          # supersample the mask, then downsample: keeps the edge smooth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--frac", type=float, default=0.82,
                    help="disc diameter as a fraction of the shorter side")
    ap.add_argument("--ring", type=int, default=6, help="ring width in px")
    a = ap.parse_args()

    im = Image.open(a.src).convert("RGBA")
    w, h = im.size
    d = int(min(w, h) * a.frac)
    d -= d % 2          # even, so the centred crop is exactly d wide
    cx, cy = w // 2, h // 2
    box = (cx - d // 2, cy - d // 2, cx + d // 2, cy + d // 2)

    mask = Image.new("L", (d * SS, d * SS), 0)
    ImageDraw.Draw(mask).ellipse((0, 0, d * SS - 1, d * SS - 1), fill=255)
    mask = mask.resize((d, d), Image.LANCZOS)

    disc = im.crop(box)
    # keep the render's own transparency: intersect it with the disc
    alpha = disc.getchannel("A").point(lambda v: v)
    disc.putalpha(Image.composite(alpha, Image.new("L", (d, d), 0), mask))

    out = Image.new("RGBA", (d, d), (0, 0, 0, 0))
    out.alpha_composite(disc)

    if a.ring:
        ring = Image.new("RGBA", (d * SS, d * SS), (0, 0, 0, 0))
        r = a.ring * SS
        ImageDraw.Draw(ring).ellipse(
            (r // 2, r // 2, d * SS - 1 - r // 2, d * SS - 1 - r // 2),
            outline=palette.NEUTRAL["700"], width=r)
        out.alpha_composite(ring.resize((d, d), Image.LANCZOS))

    out.save(a.out)
    print(f"wrote {a.out}  {out.size}  disc {a.frac:.2f} of {min(w,h)}px")


if __name__ == "__main__":
    main()
