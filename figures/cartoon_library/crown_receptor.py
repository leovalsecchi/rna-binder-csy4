"""The SMS receptor with the crown binder, in our palette.

Matteo's GFP-nanobody assembly carries the shape we want for the binder - the
VHH drawn as a rounded cup with lobes along its rim, a crown rather than a
bundle of helices. His geometry is kept exactly; what changes is

  * the cast: the GFP barrel, the coiled coil above it and his caption text go,
    so what is left is our chain - crown binder, linker, SMS;
  * the colours: his plum/butter SMS becomes our teal kinase and amber RNase,
    and the crown takes `synbio.binder`, the colour rule 6 gives the designed
    molecule.

His fills live in a `<style>` block, not on the elements, so recolouring is a
substitution inside that block - which also catches anything sharing a class.

    python code/cartoon_library/crown_receptor.py [outdir]
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1]))

from common import palette                                   # noqa: E402
from isolate_sms import drop, frame                          # noqa: E402

S = palette.SYNBIO
SRC = REPO / "inputs/cartoons"

# his hex -> ours. Every fill in the file is listed, so a colour that is not
# here is a colour he added since, and the check below will say so.
RECOLOUR = {
    "#a997ca": S["enzyme"],        # kinase lobe
    "#fae6b0": S["enzyme_alt"],    # RNase lobe
    "#87cdc7": S["binder"],        # the crown
    "#96d197": S["reporter"],      # GFP barrel (dropped, but keep it honest)
    "#b0abb5": palette.NEUTRAL["400"],   # coiled coil, back faces
    "#d4d1d7": palette.NEUTRAL["300"],   # coiled coil, front faces
    "#6c6970": S["outline"],       # every contour
}

# what our chain is not: his target, the coiled coil that dimerizes it, and the
# caption he drew into the artboard
GONE = ["GFP target module", "CC v6.3", "CC to GFP"]

CUTS = {"sms_crown_monomer": "crown_monomer.svg",
        "sms_crown_dimer": "crown_dimer.svg"}


def retheme(text):
    for old, new in RECOLOUR.items():
        text = re.sub(old, new, text, flags=re.I)
    left = {h.lower() for h in re.findall(r"#[0-9a-fA-F]{6}", text)}
    allowed = {c.lower() for c in RECOLOUR.values()}
    assert not (left - allowed), f"not in the palette: {sorted(left - allowed)}"
    return text


def strip_text(tree):
    """His artboard caption, and the font rule that styled it."""
    for parent in tree.iter():
        for child in list(parent):
            if child.tag.endswith("text"):
                parent.remove(child)


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for name, src in CUTS.items():
        tree = ET.fromstring(retheme((SRC / src).read_text()))
        drop(tree, GONE)
        strip_text(tree)
        text, (w, h) = frame(ET.tostring(tree, encoding="unicode"))
        (outdir / f"{name}.svg").write_text(text)
        print(f"{name + '.svg':24} {w:.0f} x {h:.0f}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "output/cartoon_library/svg")
