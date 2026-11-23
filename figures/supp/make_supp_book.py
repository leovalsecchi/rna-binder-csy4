"""Every supplementary figure in one SVG, one A4 page per figure.

    ~/.venvs/binder_map/bin/python3 code/supp/make_supp_book.py

Writes to selection/ (README: things arrive there freely; final/ stays manual):
  Supp_figures.svg   one A4 page per figure, side by side, each page's ids
                     prefixed so the pages cannot repaint each other
  Supp_figures.pdf   the same, as real PDF pages

Panels come from output/<slug>/*.svg — the bare panels the figure scripts emit,
since 2026-11-19 in SVG as well as PDF/PNG. SVG is the right source: matplotlib
writes it with `svg.fonttype: none`, so every label stays live text an editor can
restyle or re-font. Converting the PDFs instead would arrive as glyph outlines,
which is exactly what makes an assembled page uneditable.

Nothing is re-rendered here; this only places what is already on disk, so a page
is only as current as its last `make_figure.py` run.

Panel letters are drawn HERE and nowhere else, which is what STYLE_NOTES section 11
asks for: a script emits bare panels, the assembler places and letters them.

Why ids are prefixed: two pages built from the same figure script carry the same
element ids, and an id collision in one SVG lets a later page repaint an earlier
one (the trap recorded in `illustrator-duplicate-ids`). Every page gets `pNN_`.
"""
import re
import sys
from pathlib import Path

import fitz  # PyMuPDF

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette, style  # noqa: E402

PT_MM = 25.4 / 72.0
A4_W, A4_H = 210.0, 297.0
MARGIN = 14.0
GUTTER = 6.0
TITLE_H = 10.0

# Reading order. A slug with no output/ directory is skipped with a warning
# rather than failing the book - the set is still being built.
ORDER = [
    "target_definition",
    "scaffold_design",
    "scaffold_selection",
    "negative_controls",
    "amplicon_design",
    "repredictor_comparison",
    "metric_panel",
    "metric_reliability",
    "backbone_landscape",
    "bif_effect",
    "bif_quasicontrol",
    "library_pcr_cycle_gel",
    "readout_pcr_cycle_test",
    "missing_gel_records",
    "pooled_readout",
    "natural_rbp_benchmark",
    "metric_significance",
    "arrayed_microscopy",
    "arrayed_flow",
]


# Panels a figure script no longer produces but which are still on disk. They are
# kept there deliberately (an Illustrator assembly may still link them) but must
# not reach the book - backbone_landscape's four were superseded by one combined
# panel on 2026-11-19 and would otherwise be placed alongside their replacement.
SUPERSEDED = {
    "backbone_landscape": ("panel_a_campaign", "panel_b_gated",
                           "panel_c_ordered_panel", "panel_d_iptm_counterfactual"),
    # scaffold_design's panel (c) moved to the ORDERED design generation on
    # 2026-11-20; the first-generation drawing is in
    # archive/scaffold_design_panel_c_firstgen_2026-11-20/. Listed here as well
    # as moved, so copying those files back into output/ cannot put the
    # superseded generation into the book beside its replacement.
    "scaffold_design": ("panel_c_target_structures",),
    # repredictor_comparison's gate cascade was removed on 2026-11-20 (Leo:
    # "Remove panel b") and moved to
    # archive/repredictor_comparison_panel_b_gate_cascade_2026-11-20/. Listed
    # here as well as moved: its stem sorts BEFORE panel_b_seed_icc, so a copy
    # back into output/ would not merely add a page panel, it would take the
    # (b) letter off the panel that now holds it.
    "repredictor_comparison": ("panel_b_gate_pass_rate",),
}


def panel_pdfs(slug):
    """The same panels as PDFs, for the PDF book (fitz places PDF, not SVG)."""
    d = REPO / "output" / slug
    if not d.is_dir():
        return []
    drop = SUPERSEDED.get(slug, ())
    return sorted(p for p in d.glob("*.pdf")
                  if not p.name.startswith("._") and p.stem not in drop)


def panels(slug):
    """Panel sources, in filename order, SVG where there is one.

    The PDF set is canonical, because every panel has one. A panel is taken as
    SVG when its sibling .svg exists (live text, editable); otherwise the PDF is
    used and converted. That fallback is not a corner case - a ChimeraX render is
    a raster panel with no SVG, and dropping it would silently lose a panel from
    the page.
    """
    drop = SUPERSEDED.get(slug, ())
    out = []
    for pdf in panel_pdfs(slug):
        svg = pdf.with_suffix(".svg")
        out.append(svg if svg.exists() else pdf)
    return out


def rows_for(sizes, avail):
    """Greedy left-to-right packing into rows no wider than `avail`."""
    rows, cur, w = [], [], 0.0
    for i, (pw, _) in enumerate(sizes):
        add = pw if not cur else pw + GUTTER
        if cur and w + add > avail + 1e-6:
            rows.append(cur)
            cur, w = [i], pw
        else:
            cur.append(i)
            w += add
    if cur:
        rows.append(cur)
    return rows


def svg_size(path):
    """Authored size of a panel, in mm."""
    if path.suffix == ".pdf":
        d = fitz.open(path)
        r = d[0].rect
        d.close()
        return r.width * PT_MM, r.height * PT_MM
    head = path.read_text()[:600]
    w = re.search(r'width="([\d.]+)pt"', head)
    h = re.search(r'height="([\d.]+)pt"', head)
    if not (w and h):                       # matplotlib always writes pt; be loud
        raise ValueError(f"{path.name}: no pt width/height in the SVG header")
    return float(w.group(1)) * PT_MM, float(h.group(1)) * PT_MM


def svg_body(path, prefix):
    """One panel as an SVG fragment, with its ids namespaced."""
    if path.suffix == ".pdf":               # raster/vector panel with no SVG
        doc = fitz.open(path)
        raw = doc[0].get_svg_image(text_as_path=True)
        doc.close()
    else:
        raw = path.read_text()
    inner = raw[raw.index(">", raw.index("<svg")) + 1:raw.rindex("</svg>")]
    # Namespace every id and every reference to one. Without this, page 7's
    # "glyph-0-1" silently becomes page 2's as well.
    inner = re.sub(r'id="([^"]+)"', lambda m: f'id="{prefix}{m.group(1)}"', inner)
    inner = re.sub(r'url\(#([^)]+)\)', lambda m: f'url(#{prefix}{m.group(1)})', inner)
    inner = re.sub(r'(xlink:href|href)="#([^"]+)"',
                   lambda m: f'{m.group(1)}="#{prefix}{m.group(2)}"', inner)
    return inner


def build_page(slug, idx, x0):
    """One A4 page: title, then the figure's panels packed into rows."""
    pdfs = panels(slug)
    if not pdfs:
        return None, 0
    sizes = [svg_size(p) for p in pdfs]

    avail = A4_W - 2 * MARGIN
    # A panel wider than the text block (178 mm double-column) is scaled to fit.
    scale = min(1.0, min(avail / w for w, _ in sizes))
    sizes = [(w * scale, h * scale) for w, h in sizes]
    rows = rows_for(sizes, avail)

    out = [f'<g id="page{idx:02d}">',
           f'<rect x="{x0}" y="0" width="{A4_W}" height="{A4_H}" fill="#ffffff"/>',
           f'<text x="{x0 + MARGIN}" y="{MARGIN}" font-family="Arial" font-size="4.2" '
           f'font-weight="bold" fill="{palette.INK}">{slug}</text>',
           f'<text x="{x0 + A4_W - MARGIN}" y="{MARGIN}" font-family="Arial" '
           f'font-size="3.2" text-anchor="end" fill="{palette.NEUTRAL["500"]}">'
           f'S{idx} · {len(pdfs)} panel{"s" if len(pdfs) > 1 else ""}</text>']

    y = MARGIN + TITLE_H
    for row in rows:
        rh = max(sizes[i][1] for i in row)
        x = x0 + MARGIN
        for i in row:
            pw, ph = sizes[i]
            pre = f"p{idx:02d}_{i}_"
            # The fragment's own user units are pt; scale to mm, then by `scale`
            # if the panel had to be shrunk to the text block.
            s = scale * PT_MM
            out.append(f'<g transform="translate({x},{y}) scale({s})">')
            out.append(svg_body(pdfs[i], pre))
            out.append("</g>")
            # Panel letter, in the left margin above the panel (section 11).
            out.append(
                f'<text x="{x - 3.2}" y="{y + 3.0}" font-family="Arial" '
                f'font-size="4.0" font-weight="bold" fill="{palette.INK}">'
                f'{chr(97 + i)}</text>')
            x += pw + GUTTER
        y += rh + GUTTER
    if y > A4_H - MARGIN:
        print(f"  ! {slug}: content runs {y:.0f} mm, past the {A4_H - MARGIN:.0f} mm "
              f"page body - panels will overflow")
    out.append("</g>")
    return "\n".join(out), len(pdfs)


def main():
    style.apply()
    pages, used = [], []
    for i, slug in enumerate(ORDER, 1):
        body, n = build_page(slug, i, (len(used)) * A4_W)
        if body is None:
            print(f"  - {slug}: no output/ directory, skipped")
            continue
        pages.append(body)
        used.append((slug, n))
        print(f"  {len(used):2d}. {slug:26s} {n} panel(s)")

    total_w = len(used) * A4_W
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'xmlns:xlink="http://www.w3.org/1999/xlink" '
           f'width="{total_w}mm" height="{A4_H}mm" '
           f'viewBox="0 0 {total_w} {A4_H}">\n'
           f'<rect width="{total_w}" height="{A4_H}" fill="#ffffff"/>\n'
           + "\n".join(pages) + "\n</svg>\n")

    outdir = REPO / "selection"
    outdir.mkdir(exist_ok=True)
    svg_path = outdir / "Supp_figures.svg"
    svg_path.write_text(svg)
    print(f"\nwrote {svg_path}  ({len(used)} pages, {total_w:.0f} x {A4_H:.0f} mm)")

    # The same thing as real PDF pages, which is what "one page per figure" means
    # to anything that paginates. Built from the panel PDFs, so it stays vector.
    book = fitz.open()
    for i, (slug, _) in enumerate(used, 1):
        page = book.new_page(width=A4_W / PT_MM, height=A4_H / PT_MM)
        pdfs = panel_pdfs(slug)
        sizes = []
        for p in pdfs:
            d = fitz.open(p)
            sizes.append((d[0].rect.width * PT_MM, d[0].rect.height * PT_MM))
            d.close()
        avail = A4_W - 2 * MARGIN
        sc = min(1.0, min(avail / w for w, _ in sizes))
        sizes = [(w * sc, h * sc) for w, h in sizes]
        page.insert_text((MARGIN / PT_MM, MARGIN / PT_MM), slug,
                         fontname="Helvetica-Bold", fontsize=11)
        y = MARGIN + TITLE_H
        for row in rows_for(sizes, avail):
            rh = max(sizes[j][1] for j in row)
            x = MARGIN
            for j in row:
                pw, ph = sizes[j]
                src = fitz.open(pdfs[j])
                page.show_pdf_page(
                    fitz.Rect(x / PT_MM, y / PT_MM,
                              (x + pw) / PT_MM, (y + ph) / PT_MM), src, 0)
                src.close()
                # The letter, as in the SVG branch. Without this the two books
                # disagree, and the PDF is the one people print.
                page.insert_text(((x - 3.2) / PT_MM, (y + 3.0) / PT_MM),
                                 chr(97 + j), fontname="Helvetica-Bold",
                                 fontsize=10)
                x += pw + GUTTER
            y += rh + GUTTER
    pdf_path = outdir / "Supp_figures.pdf"
    book.save(pdf_path)
    book.close()
    print(f"wrote {pdf_path}  ({len(used)} pages)")


if __name__ == "__main__":
    main()
