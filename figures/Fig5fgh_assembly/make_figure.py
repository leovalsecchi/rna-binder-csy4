"""The H1-null / H2-partial block, assembled as a 178 mm figure.

    python code/Fig5fgh_assembly/make_figure.py [outdir]

Two rows of two. Panels are REDRAWN by their own script at the target geometry via
`variant` rather than re-set in place by `pair`: the panel script takes the
W_MM H_MM FONT_PT OUTDIR CLI, and right-anchored y labels (panels b and c are both
label-heavy) grow leftward out of their own viewBox when `pair` rewrites font-size
without being able to move them. Same route as Fig1e, Fig3e, Fig5de.

Reading order:
  f  scan       the whole feature scan as a uniform QQ plot - H1 is null
  g  h1metrics  the ten metrics the field ranks on, P(binder > non) with CIs
  h  h2legs     each leg of H2 under three stresses - the structural ones collapse
  i  arom       the survivor: aromatic fraction, library / the 32 / natural RBPs

Letters f-i because this block continues the Fig 5 page, whose d|e row is the Csy4
chemistry pair. It also stands alone; the letters are placed at assembly, not here.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common.page_assembly import (GAP_H, GAP_V, MM, PAGE_H, TOP, T_HOUSE, W,  # noqa: E402
                                  floor, page, report, variant, write)

SCRIPT = REPO / "code/Fig5fgh_h1null_h2hold/make_figure.py"
ROW_H = 46 * MM


def rows(T, workdir):
    half = (W - GAP_H) / 2
    mk = lambda letter, stem: variant(letter, SCRIPT, stem, half, ROW_H, T,
                                      workdir, sys.executable)
    out = [[mk("f", "scan.svg"), mk("g", "h1metrics.svg")],
           [mk("h", "h2legs.svg"), mk("i", "arom.svg")]]
    for p in sum(out, []):
        p.retouch(lambda u, s=p.s: floor(u, s))
    return out


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    r = rows(T_HOUSE, outdir / "_variant")
    svg, fit, bot = page(r, "Fig5 f-i  H1 null, H2 partial")
    report(f"H1 null / H2 partial, T = {T_HOUSE:.2f} pt", r, fit, bot, T_HOUSE)
    assert fit == 1.0, "must fit without shrinking, or the type is no longer T"
    write(outdir, "Fig5fgh_h1null_h2hold", svg)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "working" / "Fig5fgh_assembly")
