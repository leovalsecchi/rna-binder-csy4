"""Shared cartoon library: one drawing per file, as SVG and PNG.

Not tied to a figure. Everything imported from Matteo's snapshots is recoloured
through `common.cartoon`; the landing-pad cell is drawn here because nothing in
either repo had one.

    python code/cartoon_library/make_library.py

Writes output/cartoon_library/{svg,png}/ and an INDEX.md. PNGs are rendered with
macOS `qlmanage`, since cairosvg is deliberately not installed (see the fig01
README); without it the SVGs are still written and the PNG step is skipped.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import math
import shutil
import subprocess
from common import paths, palette, cartoon

out = paths.fig_output("cartoon_library")
svg_dir, png_dir = out / "svg", out / "png"
for d in (svg_dir, png_dir):
    d.mkdir(parents=True, exist_ok=True)

S = palette.SYNBIO
LINE = S["outline"]
PAD = 6

# ---------------------------------------------------------------- imported --
# (name, snapshot, group or None for the whole file, scale, note)
IMPORTED = [
    # --- SMS / IRE1 receptor carrying a binder ---------------------------
    ("receptor_membrane_monomer", "new_pair.svg",
     "Assembly v27.23 | selected minibinder linkers on nanobody | monomer | 1x",
     4.0, "membrane receptor, minibinder + TMD + SMS, one copy"),
    ("receptor_membrane_dimer", "new_pair.svg",
     "Assembly v27.23 | selected minibinder linkers on nanobody | dimer | 1x",
     4.0, "the same receptor, dimerized - the activating state"),
    ("receptor_v20_60_monomer", "v20.60_pair.svg",
     "Assembly v20.60 | monomer | selected v20.57 refinement | 1x", 4.0,
     "receptor v20.60, his selected refinement, one copy"),
    ("receptor_v20_60_dimer", "v20.60_pair.svg",
     "Assembly v20.60 | dimer | selected v20.57 refinement | 1x", 4.0,
     "receptor v20.60, dimerized"),
    ("receptor_cytosolic_monomer", "36_cytosolic.svg",
     "Assembly v27.24 | cytosolic GFP - nanobody | monomer | 1x", 4.0,
     "cytosolic version, GFP target + nanobody"),
    ("receptor_cytosolic_homodimer", "36_cytosolic.svg",
     "Assembly v27.24 | cytosolic GFP - nanobody | homodimer | 1x", 4.0,
     "cytosolic version, dimerized"),
    ("receptor_target_minibinder_monomer", "37_cytosolic.svg",
     "Assembly v20.67 | cytosolic Target - minibinder | monomer | 1x", 4.0,
     "cytosolic, generic target + minibinder - closest to our case"),
    ("receptor_target_minibinder_homodimer", "37_cytosolic.svg",
     "Assembly v20.67 | cytosolic Target - minibinder | homodimer | 1x", 4.0,
     "the same, dimerized"),
    ("receptor_gfp_nanobody_dimer", "source_dimer.svg",
     "Assembly v27.3 | GFP-nanobody | dimer | 1x", 4.0,
     "earlier GFP-nanobody dimer on the membrane"),
    ("sms_module", "Selected_SMS_v6_5_1x.svg", None, 5.0,
     "SMS alone: kinase lobe (enzyme) + RNase lobe (enzyme_alt)"),

    # --- target hairpin, spliced vs unspliced -----------------------------
    ("hairpin_v19_18_unspliced", "v19.18.svg", "HP v19.18 Unspliced", 5.0,
     "XBP1-type hairpin, uncleaved - his latest"),
    ("hairpin_v19_18_spliced", "v19.18.svg", "HP v19.18 Spliced", 5.0,
     "the same, cleaved"),
    ("hairpin_v19_17_unspliced", "v19.17.svg", "HP v19.17 Unspliced", 5.0,
     "variant: thinner ribbons"),
    ("hairpin_v19_17_spliced", "v19.17.svg", "HP v19.17 Spliced", 5.0, "variant"),
    ("hairpin_v19_16_unspliced", "v19.16.svg", "HP v19.16 Unspliced", 5.0,
     "variant: thinner ribbons"),
    ("hairpin_v19_16_spliced", "v19.16.svg", "HP v19.16 Spliced", 5.0, "variant"),
    ("hairpin_v19_7_unspliced", "HP_v19.7.svg", "HP v19.7 Unspliced", 5.0,
     "variant: straight stem base"),
    ("hairpin_v19_7_spliced", "HP_v19.7.svg", "HP v19.7 Spliced", 5.0, "variant"),
    ("hairpin_gapped_v19_9", "HP_v19_9.svg", None, 5.0,
     "variant: gapped pairings"),
    ("hairpin_exon_key", "Assay_colors_v1.5.svg", None, 4.0,
     "5-prime / 3-prime exon and intron colour key"),
]


# ----------------------------------------------------------- landing pad ---
def _blob(cx, cy, r, a1=0.07, a2=0.05, p1=0.0, p2=1.1, k1=3, k2=5, n=84):
    """Irregular closed outline, so a cell never reads as a circle.

    Same two-harmonic radial perturbation as `glyph_cell` in
    Fig1f_scaffold_splicing, kept in step on purpose.
    """
    pts = []
    for i in range(n + 1):
        u = 2 * math.pi * i / n
        rr = r * (1 + a1 * math.sin(k1 * u + p1) + a2 * math.sin(k2 * u + p2))
        pts.append((cx + rr * math.cos(u), cy + rr * math.sin(u)))
    return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"


def landingpad_cell():
    """A cell whose genome carries the integrated cassette.

    Cell outline, nucleus, and inside it the locus: two genome arms, the
    cassette landed between them, att junctions marked. Deliberately carries no
    text - label it in the figure, where the module names are known.
    """
    W, H = 230, 190
    cx, cy = W / 2, H / 2
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">']
    # cell
    o.append(f'<path d="{_blob(cx, cy, 86, p1=0.3, p2=1.7)}" fill="{palette.SURFACE}" '
             f'stroke="{LINE}" stroke-width="2.4" stroke-linejoin="round"/>')
    # nucleus
    o.append(f'<path d="{_blob(cx, cy - 4, 50, a1=0.05, a2=0.04, p1=1.4, p2=0.5)}" '
             f'fill="{palette.NEUTRAL["100"]}" stroke="{palette.NEUTRAL["500"]}" '
             f'stroke-width="1.8" stroke-linejoin="round"/>')
    # the locus inside it: genome arms, cassette, att junctions
    gy = cy - 4
    o.append(f'<path d="M{cx-44},{gy} L{cx-27},{gy}" stroke="{palette.NEUTRAL["500"]}" '
             f'stroke-width="2.0" stroke-linecap="round" fill="none"/>')
    o.append(f'<path d="M{cx+27},{gy} L{cx+44},{gy}" stroke="{palette.NEUTRAL["500"]}" '
             f'stroke-width="2.0" stroke-linecap="round" fill="none"/>')
    mods = [(-27, 18, S["binder"]), (-9, 14, S["enzyme"]), (5, 22, S["rna"])]
    for dx, w, fill in mods:
        o.append(f'<rect x="{cx+dx}" y="{gy-8}" width="{w}" height="16" rx="2.5" '
                 f'fill="{fill}" stroke="{LINE}" stroke-width="1.4"/>')
    for dx in (-27, 27):
        o.append(f'<path d="M{cx+dx},{gy-11} L{cx+dx},{gy+11}" '
                 f'stroke="{palette.NEUTRAL["500"]}" stroke-width="1.6" '
                 f'stroke-linecap="round" fill="none"/>')
    o.append('</svg>')
    return W, H, "".join(o)


def cell_cluster():
    """Three cells, the wet-lab node - an SVG twin of fig04's `glyph_cell`."""
    W, H = 180, 150
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">']
    for cx, cy, r, p1, p2, lw, nr in ((64, 58, 30, 0.9, 2.4, 1.8, 9),
                                      (70, 104, 32, 2.1, 0.4, 1.8, 11),
                                      (116, 80, 46, 0.3, 1.7, 2.2, 19)):
        o.append(f'<path d="{_blob(cx, cy, r, p1=p1, p2=p2)}" fill="{palette.SURFACE}" '
                 f'stroke="{LINE}" stroke-width="{lw}" stroke-linejoin="round"/>')
        o.append(f'<circle cx="{cx}" cy="{cy}" r="{nr}" fill="{LINE}"/>')
    o.append('</svg>')
    return W, H, "".join(o)


def f30_scaffold():
    from f30 import build, W as FW, H as FH
    return FW, FH, build()


def f30_scaffold_labelled():
    from f30 import build, W as FW, H as FH
    return FW, FH, build(True)


def sms_receptor():
    from receptor import build, W as RW, H as RH
    return RW, RH, build()


DRAWN = [("sms_receptor", sms_receptor,
          "SMS receptor drawn from scratch: binder, kinase lobe, RNase lobe"),
         ("f30_scaffold", f30_scaffold,
          "F30 three-way junction carrying two Csy4 hairpins"),
         ("f30_scaffold_labelled", f30_scaffold_labelled,
          "the same, with F30 / Csy4 hairpin / 5' and 3' RNA labelled"),
         ("landingpad_cell", landingpad_cell,
          "cell with the cassette integrated in its genome; no text, label it in the figure"),
         ("cell_cluster", cell_cluster,
          "three cells with nuclei - the SVG twin of fig04's glyph_cell")]


def main():
    written = []
    for name, src, group, k, note in IMPORTED:
        x0, y0, x1, y1 = cartoon.bbox(src, group)
        w, h = (x1 - x0) * k + 2 * PAD, (y1 - y0) * k + 2 * PAD
        body = cartoon.fragment(src, group, scale=k, at=(PAD, PAD))
        (svg_dir / f"{name}.svg").write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.1f}" '
            f'height="{h:.1f}" viewBox="0 0 {w:.1f} {h:.1f}">{body}</svg>',
            encoding="utf-8")
        written.append((name, f"{note} - from `{src}`"))

    for name, fn, note in DRAWN:
        _, _, svg = fn()
        (svg_dir / f"{name}.svg").write_text(svg, encoding="utf-8")
        written.append((name, note))

    # pieces cut out of the imports above, so they run after them
    import crown_receptor
    import isolate_sms
    crown_receptor.main(svg_dir)
    isolate_sms.main(svg_dir)
    written += [
        ("sms_crown_monomer", "crown binder + linker + SMS, one chain - the "
         "binder figures use this one"),
        ("sms_crown_dimer", "the same, dimerized through the kinase lobes"),
        ("sms_binder_monomer", "the older helix-bundle binder + SMS"),
        ("sms_only", "the SMS alone: kinase lobe and RNase lobe"),
    ]

    ql = shutil.which("qlmanage")
    if ql:
        files = [str(svg_dir / f"{n}.svg") for n, _ in written]
        subprocess.run([ql, "-t", "-s", "1000", "-o", str(png_dir), *files],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # qlmanage appends to the name; the SMB share also drops ._ forks
        for p in list(png_dir.glob("*.png")):
            if p.name.startswith("._"):
                p.unlink(missing_ok=True)
            elif p.name.endswith(".svg.png"):
                # rename hits EBUSY on the SMB share; copying is what works
                t = p.with_name(p.name[:-len(".svg.png")] + ".png")
                t.write_bytes(p.read_bytes())
                p.unlink(missing_ok=True)
        n_png = len([f for f in png_dir.glob("*.png")
                     if not f.name.startswith("._")])
    else:
        n_png = 0
        print("no qlmanage; SVG only")

    lines = ["# Cartoon library", "",
             "One drawing per file, cropped, in the thesis palette. SVG is the",
             "master; PNG is a 1000 px preview rendered by `qlmanage`.", "",
             "Rebuild: `python code/cartoon_library/make_library.py`", "",
             "Imported drawings are Matteo's, recoloured on import - provenance in",
             "`inputs/cartoons/SOURCES.md`. His hairpins are XBP1/IRE1, not Csy4:",
             "use them for the reporter hairpin, not for the binder's target.",
             "His two Original_Hairpins files are left out: they carry their",
             "colours in a CSS <style> block, so the recolour cannot reach them.",
             "",
             "| file | what it is |", "|---|---|"]
    lines += [f"| `{n}` | {note} |" for n, note in written]
    (out / "INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(written)} svg, {n_png} png -> {out}")


if __name__ == "__main__":
    main()
