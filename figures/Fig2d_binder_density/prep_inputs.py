"""Turn the four ChimeraX renders into the panel images Fig 2d embeds.

    python prep_inputs.py [renders_dir] [out_dir]

The renders share one camera, so they are cropped to ONE common box: trimming
each to its own bounding box slides each panel by its own density extent and
the grey hairpins stop lining up across the row.

Downscaled to 1200 px (~44 mm at 600 dpi, the width a quarter of the block
prints at) and quantized to 128 colours, which costs nothing visible on a single
hue over grey and takes each file from ~1.1 MB to ~120 KB. Four full-size RGBA
panels would put 5.5 MB of base64 into the page SVG.
"""
import os
import sys
from pathlib import Path

from PIL import Image

NAMES = ["RFD3", "BoltzGen", "BoltzDesign", "HalluDesign"]
WIDTH = 1200
COLOURS = 128


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "binder_density/renders"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "inputs"
    out.mkdir(parents=True, exist_ok=True)

    ims = {n: Image.open(src / f"density_{n}.png").convert("RGBA") for n in NAMES}
    boxes = [im.getbbox() for im in ims.values()]
    box = (min(b[0] for b in boxes), min(b[1] for b in boxes),
           max(b[2] for b in boxes), max(b[3] for b in boxes))
    print(f"common crop {box} -> {box[2] - box[0]} x {box[3] - box[1]}")

    for n in NAMES:
        im = ims[n].crop(box)
        im = im.resize((WIDTH, round(im.height * WIDTH / im.width)), Image.LANCZOS)
        im = im.quantize(colors=COLOURS, method=Image.FASTOCTREE)
        p = out / f"density_{n}.png"
        im.save(p, optimize=True)
        print(f"{n:<12} {im.width}x{im.height}  {os.path.getsize(p) // 1024} KB")


if __name__ == "__main__":
    main()
