"""Bold-line variants of the three line-art RNA drawings.

The F30 scaffold and the two AT hairpins are all drawn with one stroke width for
rails and turns (`SW`) and a thinner one for base-pair rungs (`RUNG`). A bold
variant is those two numbers scaled; nothing else about the geometry changes.

    python code/cartoon_library/make_variants.py

Writes `<name>_bold.svg` beside the originals in output/cartoon_library/.

The weights are set on the drawing modules around each call rather than by
editing them, because `f30.py` is being worked on separately and a parameter
added here would collide. If that settles, give both modules a `weight`
argument and delete this file.
"""
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parent))

from common import paths
import f30
import hairpin

BOLD = 1.55          # how much heavier the bold variant is

VARIANTS = [
    ("f30_scaffold", f30, lambda: f30.build()),
    ("f30_scaffold_labelled", f30, lambda: f30.build(True)),
    ("hairpin_at_unspliced", hairpin, lambda: hairpin.build(spliced=False)),
    ("hairpin_at_spliced", hairpin, lambda: hairpin.build(spliced=True)),
    ("hairpin_at_intron", hairpin, lambda: hairpin.build(intron=True)),
]


@contextmanager
def weight(module, factor):
    """Draw this module heavier, then put it back exactly as it was."""
    before = (module.SW, module.RUNG)
    module.SW, module.RUNG = module.SW * factor, module.RUNG * factor
    try:
        yield
    finally:
        module.SW, module.RUNG = before


def main():
    lib = paths.fig_output("cartoon_library")
    svg_dir, png_dir = lib / "svg", lib / "png"
    for d in (svg_dir, png_dir):
        d.mkdir(parents=True, exist_ok=True)

    written = []
    for name, module, build in VARIANTS:
        with weight(module, BOLD):
            svg = build()
        p = svg_dir / f"{name}_bold.svg"
        p.write_text(svg, encoding="utf-8")
        written.append(p)
    # the modules must be back at their normal weight for everything else
    assert (f30.SW, hairpin.SW) == (hairpin.SW * hairpin.SCALE, 1.30), \
        "a weight override leaked"

    ql = shutil.which("qlmanage")
    if ql:
        subprocess.run([ql, "-t", "-s", "1000", "-o", str(png_dir),
                        *map(str, written)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for q in list(png_dir.glob("*.svg.png")):
            if q.name.startswith("._"):
                q.unlink(missing_ok=True)
                continue
            # rename hits EBUSY on the SMB share; copying is what works
            target = q.with_name(q.name[:-len(".svg.png")] + ".png")
            target.write_bytes(q.read_bytes())
            q.unlink(missing_ok=True)
    print(f"wrote {len(written)} bold variants to {svg_dir}")


if __name__ == "__main__":
    main()
