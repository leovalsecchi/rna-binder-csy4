"""Write panel B's pieces as separate SVGs, for hand assembly.

One file per piece, cropped to its content, already in the thesis palette, text
kept as text. Arrange them by hand; nothing here assumes a layout.

Everything is generated from the shared builders, so a piece here is the same
piece panel d uses - edit `common/constructs.py` or the palette and rerun.

    python code/Fig1b_sms_assay/make_elements.py

Molecules are NOT duplicated here: the scaffold, receptors, hairpins and SMS
module live in `output/cartoon_library/`, and copying them would only let the
two folders drift.
"""
import re
import sys
from pathlib import Path
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parents[1] / "cartoon_library"))

import yaml
import f30
import hairpin
from common import paths, palette, cartoon
from common import constructs as K

out = paths.fig_output("Fig1b_sms_assay") / "elements"
out.mkdir(parents=True, exist_ok=True)

_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
MUT = palette.NEUTRAL["600"]
INK = palette.INK
LINE = palette.SYNBIO["outline"]
BG = palette.SURFACE

PAD = 10
Y = 60                 # baseline every bar is drawn on
written = []


def doc(name, w, h, body, note, x0=0.0, y0=0.0):
    (out / f"{name}.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" '
        f'height="{h:.0f}" viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}">'
        f'{body}</svg>', encoding="utf-8")
    written.append((name, note))


def bar_doc(name, items, note):
    """One construct bar, cropped around its strand and anything above it."""
    svg = K.bar(PAD + 14, Y, items)
    w = K.width(items) + 2 * (PAD + 14)
    # ticks and the stop glyph reach above the strand; CDS depth is H/2
    top, bottom = Y - 46, Y + K.H / 2 + PAD
    doc(name, w, bottom - top, svg, note, 0.0, top)


def main():
    # this folder is generated, not kept; the SMB share also drops ._ forks
    for stale in list(out.glob("*.svg")):
        try:
            stale.unlink(missing_ok=True)
        except OSError:
            pass   # held open by a previewer; it gets overwritten below anyway

    for key, note in (("reporter_rna", "construct 1: the target RNA"),
                      ("sensor", "construct 2: the SMS receptor"),
                      ("readout", "construct 3: the reporter"),
                      ("donor", "panel d: the donor cassette"),
                      ("integrated", "panel d: after integration")):
        bar_doc(f"bar_{key}", K.BARS[key], note)

    for key in K.MODULES:
        bar_doc(f"mod_{key}", [key], f"module: {K.MODULES[key][0] or key}")

    # the two reporter states, as whole strands ready to place
    for spliced in (False, True):
        name = "state_spliced" if spliced else "state_unspliced"
        # AUG and the stop are both marks on the strand and take no width, so
        # they cannot share a bar; the stop ships as mod_stop, to be placed
        # between the hairpin and mScarlet
        items = ["aug", "mscarlet" if spliced else "mscarlet_off"]
        bar_doc(name, items,
                "reporter mRNA, cleaved: add the spliced hairpin after AUG"
                if spliced else
                "reporter mRNA, off: add the unspliced hairpin, then mod_stop")

    # the AT hairpin comes from its own builder, which is traced off the panel
    # this one replaces - not from cartoon.fragment("v19.18.svg"), whose hairpin
    # is a different, compact drawing and would disagree with the library
    for name, kw, note in (
            ("hairpin_unspliced", dict(), "the AT hairpin, to sit on a strand"),
            ("hairpin_spliced", dict(spliced=True), "the same after the cut"),
            ("hairpin_intron", dict(intron=True),
             "the excised intron, drifting free")):
        (out / f"{name}.svg").write_text(hairpin.build(**kw), encoding="utf-8")
        written.append((name, note))

    # the F30 scaffold, generated from its own builder rather than copied, so
    # there is still one source for it and the two folders cannot drift
    for name, labelled, note in (
            ("f30_scaffold", False, "the F30 junction carrying two Csy4 hairpins"),
            ("f30_scaffold_labelled", True, "the same, labelled")):
        (out / f"{name}.svg").write_text(f30.build(labelled), encoding="utf-8")
        written.append((name, note))

    lines = ["# Panel B elements", "",
             "One file per piece, cropped, in the thesis palette, text as text.",
             "Generated - do not edit by hand; rerun",
             "`python code/Fig1b_sms_assay/make_elements.py`.", "",
             "Bars and modules come from `common/constructs.py`, so they are the",
             "same shapes panel d uses. **Molecules are not here**, with one",
             "exception: the F30 scaffold is built here from `cartoon_library/",
             "f30.py`, the same builder `output/cartoon_library/` uses, so the",
             "two stay in step. The receptors, SMS module, binder and XBP1",
             "hairpins are in `output/cartoon_library/`.",
             "", "Every bar is drawn on a baseline 60 units from the top of its",
             "own box, so bars line up if you align their boxes.", "",
             "| file | what it is |", "|---|---|"]
    lines += [f"| `{n}.svg` | {note} |" for n, note in written]
    (out / "INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(written)} elements to {out}")


if __name__ == "__main__":
    main()
