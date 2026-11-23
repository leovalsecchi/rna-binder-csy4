"""The H1-null / refolder-reproducibility pair, as a 178 mm block.

    python code/Fig5fg_assembly/make_figure.py [outdir]

One row of two, both panels REDRAWN at the target geometry by their own script via
`variant` (they are y-label heavy, and `pair` cannot move a right-anchored label it
has just resized - see Fig5de's README).

  f  h1null    ten metrics x two refolders, P(binder > non) with bootstrap CIs
  g  refolder  AF3 vs ESMFold2 agreement per interface measurement, 11,999 designs
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common.page_assembly import (GAP_H, MM, PAGE_H, TOP, T_HOUSE, W,  # noqa: E402
                                  floor, page, report, variant, write)

SCRIPT = REPO / "code/Fig5fg_h1null_refolder/make_figure.py"
ROW_H = 52 * MM


def rows(T, workdir):
    half = (W - GAP_H) / 2
    out = [[variant("f", SCRIPT, "h1null.svg", half, ROW_H, T, workdir, sys.executable),
            variant("g", SCRIPT, "refolder.svg", half, ROW_H, T, workdir, sys.executable)]]
    for p in sum(out, []):
        p.retouch(lambda u, s=p.s: floor(u, s))
    return out


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    r = rows(T_HOUSE, outdir / "_variant")
    svg, fit, bot = page(r, "H1 null / refolder reproducibility")
    report(f"H1 null + refolder, T = {T_HOUSE:.2f} pt", r, fit, bot, T_HOUSE)
    assert fit == 1.0, "must fit without shrinking, or the type is no longer T"
    write(outdir, "Fig5fg_h1null_refolder", svg)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "working" / "Fig5fg_assembly")
