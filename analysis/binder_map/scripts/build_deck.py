#!/usr/bin/env python3
"""Assemble the cluster deck from ChimeraX-rendered panels.

One slide per cluster: cluster id as the title, four designs, and the stack
breakdown as percentages with counts.

Every panel was superposed on the same RNA target before rendering, so the
target sits identically in all four images and they differ only in how the
binder approaches it.

Usage:
    python scripts/build_deck.py --structures dist/deck_structures \
        --images dist/deck_images --out dist/cluster_deck.pptx
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

PALETTE = {
    "RFD3+LiMPNN": "4a7ebb", "BG+LiMPNN": "e07b39", "BG+BIF": "4f9d69",
    "BD1+LiMPNN": "b5514f", "HD+LiMPNN": "8264a5",
}
STACKS = list(PALETTE)


def rgb(hexstr: str) -> RGBColor:
    return RGBColor(int(hexstr[0:2], 16), int(hexstr[2:4], 16), int(hexstr[4:6], 16))


def textbox(slide, left, top, width, height, text, size, *, bold=False,
            colour="333333", align_right=False):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = bold
    p.font.color.rgb = rgb(colour)
    if align_right:
        from pptx.enum.text import PP_ALIGN
        p.alignment = PP_ALIGN.RIGHT
    return box


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", type=Path, default=Path("dist/deck_structures"))
    ap.add_argument("--images", type=Path, default=Path("dist/deck_images"))
    ap.add_argument("--out", type=Path, default=Path("dist/cluster_deck.pptx"))
    args = ap.parse_args()

    meta = json.loads((args.structures / "deck_meta.json").read_text())
    panels_by_cluster: dict[int, list[dict]] = defaultdict(list)
    for p in meta["panels"]:
        img = args.images / f"{Path(p['pdb']).stem}.png"
        if img.exists():
            p["img"] = img
            panels_by_cluster[p["cluster"]].append(p)

    clusters = [c for c in meta["clusters"] if panels_by_cluster.get(c["cluster"])]
    n_bb, n_des = meta["n_backbones"], meta["n_designs"]

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]

    # ---------------------------------------------------------- title slide
    s = prs.slides.add_slide(blank)
    textbox(s, Inches(0.7), Inches(0.7), Inches(12.0), Inches(0.9),
            f"Binder fold clusters — {meta['layer']} layer", 34, bold=True,
            colour="1b1b1b")
    textbox(s, Inches(0.7), Inches(1.6), Inches(12.0), Inches(0.5),
            f"{n_bb:,} unique backbones from {n_des:,} passing designs  ·  "
            f"{meta['n_clusters_total']:,} clusters  ·  "
            f"{len(clusters)} shown", 16, colour="444444")

    body = [
        "Method",
        "   Foldseek all-vs-all on the protein backbone; similarity = min(qTM, tTM).",
        "   Edges at TM ≥ 0.40, top-15 neighbours, Leiden clustering (seed 42).",
        "   Backbones deduplicated by coordinate hash, not by design name.",
        "",
        "Reading the panels",
        "   All four structures are superposed on the RNA target, so the target is in",
        "   the same orientation everywhere and the panels differ only in how the",
        "   binder approaches it.",
        "   Protein: helix red, strand amber, coil grey.  RNA target: teal.",
        "   The bar above each panel marks its generation stack.",
        "",
        "Upstream filter",
        "   (rmsd_af3 ≤ 3 OR rmsd_esm ≤ 3) AND (dockq_af3 ≥ 0.7 OR dockq_esm ≥ 0.7)",
        "   AND seq_run ≤ 7 AND mean ipTM ≥ 0.6.  No interface-quality filtering.",
    ]
    tb = s.shapes.add_textbox(Inches(0.7), Inches(2.3), Inches(12.0), Inches(4.6))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, line in enumerate(body):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line
        header = line and not line.startswith(" ")
        p.font.size = Pt(13 if not header else 14)
        p.font.bold = bool(header)
        p.font.color.rgb = rgb("1b1b1b" if header else "555555")

    # ------------------------------------------------------- cluster slides
    POS = [(0.45, 1.5), (3.75, 1.5), (0.45, 4.4), (3.75, 4.4)]
    for rank, info in enumerate(clusters, start=1):
        cid = info["cluster"]
        slide = prs.slides.add_slide(blank)

        textbox(slide, Inches(0.45), Inches(0.28), Inches(8.0), Inches(0.7),
                f"Fold cluster {cid}", 30, bold=True, colour="1b1b1b")
        textbox(slide, Inches(0.45), Inches(0.92), Inches(8.5), Inches(0.45),
                f"{info['backbones']:,} backbones   ·   {info['designs']:,} designs"
                f"   ·   {info['pct']:.2f}% of layer   ·   rank {rank} of "
                f"{len(clusters)}", 14, colour="555555")

        for panel in sorted(panels_by_cluster[cid], key=lambda d: d["slot"])[:4]:
            left, top = POS[panel["slot"]]
            slide.shapes.add_picture(str(panel["img"]), Inches(left), Inches(top),
                                     height=Inches(2.8))
            # stack colour bar + design id under each panel
            bar = slide.shapes.add_textbox(Inches(left), Inches(top + 2.78),
                                           Inches(3.05), Inches(0.28))
            bp = bar.text_frame.paragraphs[0]
            bp.text = panel["design_id"][:44]
            bp.font.size = Pt(7.5)
            bp.font.color.rgb = rgb(PALETTE.get(panel["stack"], "999999"))
            bp.font.bold = True

        # ---- composition panel
        textbox(slide, Inches(7.25), Inches(1.55), Inches(5.6), Inches(0.4),
                "Stack composition", 15, bold=True, colour="1b1b1b")
        y = 2.05
        total = info["backbones"]
        for st in STACKS:
            c = int(info["stacks"].get(st, 0))
            if not c:
                continue
            pct = 100.0 * c / total
            textbox(slide, Inches(7.25), Inches(y), Inches(2.6), Inches(0.32),
                    st, 12, bold=True, colour=PALETTE[st])
            textbox(slide, Inches(9.85), Inches(y), Inches(1.1), Inches(0.32),
                    f"{pct:5.1f}%", 12, colour="333333", align_right=True)
            textbox(slide, Inches(11.0), Inches(y), Inches(1.3), Inches(0.32),
                    f"{c:,}", 12, colour="777777", align_right=True)
            # proportion bar
            barw = max(0.04, 4.9 * pct / 100.0)
            from pptx.enum.shapes import MSO_SHAPE
            shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(7.25),
                                         Inches(y + 0.26), Inches(barw), Inches(0.07))
            shp.fill.solid()
            shp.fill.fore_color.rgb = rgb(PALETTE[st])
            shp.line.fill.background()
            shp.shadow.inherit = False
            y += 0.58

        textbox(slide, Inches(7.25), Inches(6.5), Inches(5.6), Inches(0.75),
                "helix red · strand amber · coil grey · RNA target teal — "
                "all four panels superposed on the RNA target, so they differ "
                "only in how the binder approaches it",
                10, colour="999999")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(args.out)
    print(f"wrote {args.out}  ({len(clusters) + 1} slides, "
          f"{sum(len(v) for v in panels_by_cluster.values())} structures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
