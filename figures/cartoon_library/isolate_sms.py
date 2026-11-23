"""Cut single pieces out of the receptor monomer, each with a viewBox that fits.

Matteo's monomer carries the whole chain plus the target it grips. Dropping a
named group leaves every transform above it untouched, so the geometry is his,
unmoved - only the frame is recomputed, by rendering the file wide and
measuring where the ink actually lands (no SVG geometry engine here, and a
path's control points are not its bounding box).

    python code/cartoon_library/isolate_sms.py [outdir]
"""
import re, sys
import xml.etree.ElementTree as ET
from pathlib import Path

NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)
HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
SRC = REPO / "output/cartoon_library/svg/receptor_target_minibinder_monomer.svg"

# what each cut keeps, by the groups it drops
# "CC v6.3" is the coiled coil that dimerizes the target, not part of our
# chain, so it goes with the target in both cuts.
CUTS = {
    "sms_binder_monomer": ["Target module", "CC v6.3"],
    "sms_only": ["Target module", "CC v6.3", "Minibinder v6.60",
                 "Minibinder to IRE1"],
}


def drop(tree, prefixes):
    """Remove every group whose data-name starts with one of these."""
    for parent in tree.iter():
        for child in list(parent):
            n = child.get("data-name", "")
            if any(n.startswith(p) for p in prefixes):
                parent.remove(child)


def frame(text, margin=4.0):
    """Re-frame on the geometry: every path here is M/C/L/Z under at most one
    translate+scale transform,
    so the control-point hull is the bounding box (a cubic never leaves its
    hull), widened by half a stroke. No renderer in the loop."""
    m = re.search(
        r'transform="translate\(([-\d.]+),([-\d.]+)\)\s*scale\(([\d.]+)\)"', text)
    tx, ty, sc = (float(v) for v in m.groups()) if m else (0.0, 0.0, 1.0)
    pts = [float(v) for d in re.findall(r'\sd="([^"]+)"', text)
           for v in re.findall(r'-?\d+\.?\d*', d)]
    assert pts and len(pts) % 2 == 0, "path data did not read as x,y pairs"
    xs, ys = pts[0::2], pts[1::2]
    # a stroke sits half outside the path, and it may be set either way round
    widths = re.findall(r'stroke-width[:="]+\s*([\d.]+)', text)
    half = (max(float(v) for v in widths) / 2) if widths else 0.0
    x0, y0 = (min(xs) - half) * sc + tx, (min(ys) - half) * sc + ty
    x1, y1 = (max(xs) + half) * sc + tx, (max(ys) + half) * sc + ty
    x0, y0, w, h = x0 - margin, y0 - margin, x1 - x0 + 2 * margin, y1 - y0 + 2 * margin
    out = re.sub(r'viewBox="[^"]+"', f'viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}"',
                 text, count=1)
    if re.search(r'<svg[^>]*\swidth="', out):
        out = re.sub(r'width="[^"]+"', f'width="{w:.1f}"', out, count=1)
        out = re.sub(r'height="[^"]+"', f'height="{h:.1f}"', out, count=1)
    else:
        out = out.replace("<svg", f'<svg width="{w:.1f}" height="{h:.1f}"', 1)
    return out, (w, h)


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for name, gone in CUTS.items():
        tree = ET.fromstring(SRC.read_text())
        drop(tree, gone)
        text, (w, h) = frame(ET.tostring(tree, encoding="unicode"))
        (outdir / f"{name}.svg").write_text(text)
        print(f"{name+'.svg':28} {w:.0f} x {h:.0f}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "output/cartoon_library/svg")
