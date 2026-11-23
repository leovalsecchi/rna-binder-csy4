"""Fig 1 panel b - the SMS reporter assay, in Leo's own arrangement.

Replaces the earlier `make_panel.py`, which drew its parts here from scratch.
**Every element in this one is Matteo's**, imported through `common.cartoon`
or lifted from a file built out of his art:

| slot | element |
|---|---|
| mechanism, left | `Fig1d_mechanism_hairpin_V2.svg` - Leo's Illustrator arrangement of Matteo's VHH v7.9, SMS v6.5, linker, target and hairpin. The FULL four-state version, not the compact one the Fig 1d panel uses |
| hairpins in the two rows | his `v19.18.svg`, unspliced and spliced |
| mScarlet protein | his GFP barrel `Selected_GFP_v1_1x.svg`, recoloured green to red |
| module blocks | his block idiom - rounded box, 6 deep, hairline in `synbio.outline` - the same one the retheme pipeline gives the cassette |

Nothing from `f30.py`, `receptor.py`, `mechanism.py` or `constructs.py` is
used: those are ours, and Leo asked for his.

Layout, which is Leo's and is the point of the rebuild: the mechanism runs
down the left; the uncut reporter sits under it; a Frameshift arrow crosses to
the spliced reporter on the right, which makes mScarlet; a grey box across the
top lists the three transfected components.

    python code/Fig1b_sms_assay/make_panel_b.py
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cartoon, palette                           # noqa: E402
from common.paths import ROOT                                 # noqa: E402

OUT = ROOT / "output/Fig1b_sms_assay"
MECH = ROOT / "working/Fig1d_integrated_construct/Fig1d_mechanism_hairpin_V2.svg"

S, N, TH = palette.SYNBIO, palette.NEUTRAL, palette.THEME
INK, LINE = palette.INK, S["outline"]
OFF = N["300"]                                 # a reporter that is not made

# His cassette band: blocks 6 deep on a strand through their middle, corners
# at 0.8. Kept so this panel and the Fig 1d rows read as the same object.
BAND, R = 6.0, 0.8
KAPPA = 0.4477


def block(x0, x1, ym, fill, label=None, tcol=palette.SURFACE, fs=4.4):
    """One module, in his idiom: rounded box, hairline, label inside."""
    yt, yb, c = ym - BAND / 2, ym + BAND / 2, R * KAPPA
    d = (f"M{x0 + R:g} {yt:g} L{x1 - R:g} {yt:g} "
         f"C{x1 - c:g} {yt:g} {x1:g} {yt + c:g} {x1:g} {yt + R:g} "
         f"L{x1:g} {yb - R:g} C{x1:g} {yb - c:g} {x1 - c:g} {yb:g} "
         f"{x1 - R:g} {yb:g} L{x0 + R:g} {yb:g} "
         f"C{x0 + c:g} {yb:g} {x0:g} {yb - c:g} {x0:g} {yb - R:g} "
         f"L{x0:g} {yt + R:g} C{x0:g} {yt + c:g} {x0 + c:g} {yt:g} "
         f"{x0 + R:g} {yt:g} Z")
    o = (f'<path d="{d}" fill="{fill}" stroke="{LINE}" stroke-width="0.4" '
         f'stroke-linejoin="round"/>')
    if label:
        o += (f'<text x="{(x0 + x1) / 2:g}" y="{ym + fs * 0.35:g}" '
              f'font-size="{fs:g}" fill="{tcol}" class="in">{label}</text>')
    return o


def hairpin(name, foot, ym, k):
    """One of his v19.18 hairpins, stem foot landed on the strand.

    The excised intron imports as neutral - the standard map treats it as
    structure - but COLORS.md rule 4 gives the warm end of the diverging ramp
    to the piece the hairpin loses, and the mechanism beside it already draws
    it that way. Recoloured here only, so the global map is left alone.
    """
    x0, y0, x1, y1 = cartoon.bbox("v19.18.svg", name)
    f = cartoon.fragment("v19.18.svg", name, scale=k,
                         dx=foot - (x0 + x1) / 2 * k, dy=ym - y1 * k)
    return re.sub(r'(<path[^>]*data-name="Excised intron"[^>]*?)'
                  r'stroke="#[0-9a-fA-F]{6}"',
                  lambda m: m.group(1) + f'stroke="{palette.SYNBIO["intron"]}"', f)


def barrel(cx, cy, k, col):
    """His GFP barrel, recoloured. Green is the only fill in that file."""
    f = cartoon.fragment("Selected_GFP_v1_1x.svg", "GFP v1 | 1x reusable",
                         scale=k, at=(cx, cy), anchor="center")
    return f.replace(S["reporter"], col)


def glow(cx, cy, r=3.4, n=8):
    """Rays, the one mark with no source in his set - he has no lightbulb."""
    import math
    o = []
    for i in range(n):
        a = 2 * math.pi * i / n - math.pi / 2
        o.append(f'M{cx + r * math.cos(a):.2f} {cy + r * math.sin(a):.2f} '
                 f'L{cx + (r + 2.2) * math.cos(a):.2f} '
                 f'{cy + (r + 2.2) * math.sin(a):.2f}')
    return (f'<path d="{" ".join(o)}" fill="none" stroke="{S["enzyme_alt"]}" '
            f'stroke-width="0.7" stroke-linecap="round"/>')


def reporter_row(x0, x1, ym, spliced):
    """5' - hairpin - ATG - mScarlet - 3'. The two rows stay in register:
    same ATG, same block, only the hairpin differs."""
    hp = "HP v19.18 Spliced" if spliced else "HP v19.18 Unspliced"
    k = 0.95
    hx = x0 + 26
    mx0, mx1 = x0 + 62, x1 - 12
    o = [f'<path d="M{x0:g} {ym:g} L{x1:g} {ym:g}" fill="none" '
         f'stroke="{INK}" stroke-width="0.4" stroke-linecap="round"/>',
         hairpin(hp, hx, ym, k),
         block(mx0, mx1, ym, S["reporter_alt"] if spliced else OFF,
               "mScarlet", palette.SURFACE if spliced else N["600"]),
         f'<path d="M{mx0:g} {ym + 3.6:g} L{mx0:g} {ym + 5.6:g}" fill="none" '
         f'stroke="{N["500"]}" stroke-width="0.35"/>',
         f'<text x="{mx0:g}" y="{ym + 11:g}" class="lbl">ATG</text>',
         f'<text x="{x0 - 3:g}" y="{ym + 1.9:g}" class="lbl">5′</text>',
         f'<text x="{x1 + 3:g}" y="{ym + 1.9:g}" class="lbl">3′</text>']
    return "".join(o)


def component_bar(x0, ym, parts, w=78):
    """One transfected construct for the grey box."""
    o = [f'<path d="M{x0 - 4:g} {ym:g} L{x0 + w + 4:g} {ym:g}" fill="none" '
         f'stroke="{INK}" stroke-width="0.35" stroke-linecap="round"/>']
    tot = sum(p[0] for p in parts)
    x = x0
    for frac, fill, lab, tcol in parts:
        wd = w * frac / tot
        o.append(block(x, x + wd, ym, fill, lab, tcol, 3.9))
        x += wd
    return "".join(o)


def main():
    o = []
    # ---- the mechanism, down the left -------------------------------------
    t = MECH.read_text()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    # keep <defs>: it carries the <style> block that fills his shapes. Strip
    # only the title. (Dropping defs renders the whole mechanism solid black.)
    body = re.sub(r"<title\b.*?</title>", "", body, flags=re.S)
    mk = 138.0 / vb[3]        # clear of the state label below
    o.append(f'<g transform="translate(6 12) scale({mk:.5f}) '
             f'translate({-vb[0]:g} {-vb[1]:g})">{body}</g>')
    mw = vb[2] * mk

    # ---- grey box of components, across the top ---------------------------
    bx0, bx1, by0 = 118.0, 340.0, 6.0
    o.append(f'<rect x="{bx0:g}" y="{by0:g}" width="{bx1 - bx0:g}" '
             f'height="62" rx="2" fill="{N["100"]}" stroke="{N["300"]}" '
             f'stroke-width="0.5"/>')
    o.append(f'<text x="{bx0 + 6:g}" y="{by0 + 9:g}" class="hd">Components</text>')
    rows = [("1", [(1.0, N["200"], "5′ UTR", INK),
                   (1.5, S["reporter"], "mNeonGreen", INK),
                   (1.4, TH["200"], "RNA scaffold", INK)], "target RNA"),
            ("2", [(1.3, S["binder"], "Binder", palette.SURFACE),
                   (0.25, N["300"], None, INK),
                   (1.5, S["enzyme"], "SMS", INK),
                   (0.55, S["enzyme_alt"], None, INK)], "SMS receptor"),
            ("3", [(0.8, palette.SYNBIO["intron"], "AT", palette.SURFACE),
                   (1.6, S["reporter_alt"], "mScarlet", palette.SURFACE)],
             "reporter")]
    for i, (num, parts, cap) in enumerate(rows):
        ym = by0 + 20 + i * 15
        o.append(f'<circle cx="{bx0 + 10:g}" cy="{ym:g}" r="3.4" fill="none" '
                 f'stroke="{N["500"]}" stroke-width="0.4"/>'
                 f'<text x="{bx0 + 10:g}" y="{ym + 1.4:g}" class="num">{num}</text>')
        o.append(component_bar(bx0 + 20, ym, parts))
        o.append(f'<text x="{bx0 + 106:g}" y="{ym + 1.5:g}" class="cap">{cap}</text>')

    # ---- the two reporter states ------------------------------------------
    ry = 186.0
    o.append(f'<text x="8" y="{ry - 28:g}" class="st">Uncut</text>')
    o.append(reporter_row(10, 146, ry, False))
    o.append(f'<text x="196" y="{ry - 28:g}" class="st">Spliced</text>')
    o.append(reporter_row(198, 334, ry, True))

    # the frameshift arrow between them
    o.append(f'<path d="M156 {ry:g} L186 {ry:g}" fill="none" stroke="{N["500"]}" '
             f'stroke-width="0.9"/>'
             f'<path d="M186 {ry:g} L182 {ry - 1.9:g} L182 {ry + 1.9:g} Z" '
             f'fill="{N["500"]}"/>')
    o.append(f'<text x="171" y="{ry - 4:g}" class="lbl">Frameshift</text>')

    # mScarlet made, over the spliced block
    o.append(barrel(300, ry - 18, 0.95, S["reporter_alt"]))
    o.append(glow(307, ry - 25))
    o.append(f'<path d="M300 {ry - 11:g} L300 {ry - 5:g}" fill="none" '
             f'stroke="{N["500"]}" stroke-width="0.5"/>')

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-6 0 352 206">'
           f'<style>'
           f'text{{font-family:Arial;fill:{INK}}}'
           f'.lbl{{font-size:4.6px;fill:{N["600"]};text-anchor:middle}}'
           f'.in{{text-anchor:middle}}'
           f'.hd{{font-size:4.6px;fill:{N["600"]};letter-spacing:0.4px}}'
           f'.num{{font-size:4.2px;fill:{N["600"]};text-anchor:middle}}'
           f'.cap{{font-size:4.4px;fill:{N["600"]}}}'
           f'.st{{font-size:5.5px;font-weight:bold}}'
           f'</style>' + "".join(o) + "</svg>")

    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "Fig1b_panel_b.svg"
    p.write_text(svg, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-w", "2000", "-b", "white", str(p),
                    "-o", str(p.with_suffix(".png"))], check=False)
    subprocess.run(["inkscape", str(p), "--export-type=pdf",
                    f"--export-filename={p.with_suffix('.pdf')}"], check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    allowed = {c.lower() for c in (*TH.values(), *N.values(), *S.values(),
                                   *palette.DIVERGING, INK, palette.GRID,
                                   palette.SURFACE)}
    stray = sorted({c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", svg)}
                   - allowed)
    print(f"wrote {p}  mechanism {mw:.0f} wide"
          + (f"  LEFTOVER {stray}" if stray else "  palette ok"))


if __name__ == "__main__":
    main()
