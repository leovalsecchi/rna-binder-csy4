#!/usr/bin/env python3
"""PowerPoint deck: one slide per cluster, four designs each.

Every structure is drawn in the RNA target's own frame, so the target sits in
the same place in all four panels and the panels differ only in where and how
the binder approaches. Rendering each design in its own arbitrary orientation
would make a page of four pictures that cannot be compared.

Rendered with matplotlib rather than ChimeraX because ChimeraX has no OpenGL
offscreen on this machine; the RNA-frame alignment is worth more here than
ray-traced shading.

Usage:
    python scripts/make_cluster_deck.py --layer design --axis fold \
        --resolution cluster_r10 --min-size 4
"""
from __future__ import annotations

import argparse
import gzip
import io as _io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Emu, Inches, Pt

N_RES, N_NT = 80, 16
N_PROT = N_RES * 4
PROT_ATOMS = ["N", "CA", "C", "O"]

PALETTE = {
    "RFD3+LiMPNN": "#4a7ebb", "BG+LiMPNN": "#e07b39", "BG+BIF": "#4f9d69",
    "BD1+LiMPNN": "#b5514f", "HD+LiMPNN": "#8264a5",
}
STACKS = list(PALETTE)
SS_COLOUR = {"H": "#d1495b", "E": "#edae49", "C": "#8d99ae"}
RNA_COLOUR = "#2a9d8f"


# ---------------------------------------------------------------- geometry
def rna_frame(rna_xyz: np.ndarray):
    """Origin + rotation from the RNA itself, matching binder_map.geometry.

    Built from centroid differences (origin->loop, then the 5' half
    orthogonalised) so it is equivariant: rotating the input rotates the frame
    identically, and no sign choice is ever made.
    """
    origin = rna_xyz.mean(axis=0)
    per_nt = rna_xyz.reshape(N_NT, 3, 3).mean(axis=1)      # centroid per nucleotide
    v1 = per_nt[6:10].mean(axis=0) - origin                # toward the loop
    v2 = per_nt[0:4].mean(axis=0) - origin                 # toward the 5' end
    e1 = v1 / (np.linalg.norm(v1) or 1.0)
    v2p = v2 - (v2 @ e1) * e1
    e2 = v2p / (np.linalg.norm(v2p) or 1.0)
    return origin, np.stack([e1, e2, np.cross(e1, e2)])


def decode(blob: np.ndarray, slot: int, grid: float) -> np.ndarray:
    return np.cumsum(blob[slot].astype(np.int32), axis=0) * grid


# ---------------------------------------------------------------- drawing
def draw(ax, xyz: np.ndarray, ss: str):
    """Protein backbone coloured by secondary structure, plus the RNA target."""
    prot, rna = xyz[:N_PROT], xyz[N_PROT:]
    origin, rot = rna_frame(rna)
    prot = (prot - origin) @ rot.T
    rna = (rna - origin) @ rot.T

    ca = prot[1::4]                                        # CA of each residue

    def smooth(p, w):
        """Running mean. A helix CA trace oscillates with ~5.4 A pitch, so a
        narrow window leaves a zig-zag that reads as noise rather than as a
        helix; averaging over a full turn recovers the helical axis."""
        if w <= 0:
            return p.copy()
        out = np.empty_like(p)
        for i in range(len(p)):
            lo, hi = max(0, i - w), min(len(p), i + w + 1)
            out[i] = p[lo:hi].mean(axis=0)
        return out

    # A single box filter cannot null a 3.6-residue-per-turn oscillation - at
    # w=4 the trace still turns 24 deg per step and reads as noise. Cascading
    # two w=3 filters approximates a Gaussian and gets it to 6 deg while
    # staying 2.6 A from the true CA path, so the fold is still faithful.
    helix_axis = smooth(smooth(ca, 3), 3)
    light = smooth(ca, 1)
    is_h = np.array([ss[i] == "H" if i < len(ss) else False for i in range(len(ca))])
    blend = np.convolve(is_h.astype(float), np.ones(5) / 5, mode="same")[:, None]
    ca_s = helix_axis * blend + light * (1 - blend)

    def depth_shade(colour, z, lo, hi):
        """Fade with depth so the fold reads in 3D instead of as a flat tangle."""
        t = 0.35 + 0.65 * (z - lo) / (hi - lo + 1e-9)
        rgb = np.array(matplotlib.colors.to_rgb(colour))
        return tuple(1 - (1 - rgb) * t)

    zs = ca_s[:, 2]
    lo, hi = float(zs.min()), float(zs.max())

    def segs_2d(points):
        return np.stack([points[:-1, :2], points[1:, :2]], axis=1)

    # Painter's algorithm: far segments first, so near ones overlap them.
    mid_z = (ca_s[:-1, 2] + ca_s[1:, 2]) / 2
    order = np.argsort(mid_z)
    seg = segs_2d(ca_s)[order]
    cols, widths = [], []
    for i in order:
        code = ss[i] if i < len(ss) else "C"
        cols.append(depth_shade(SS_COLOUR.get(code, "#8d99ae"), mid_z[i], lo, hi))
        widths.append(5.0 if code == "H" else (4.4 if code == "E" else 2.4))
    ax.add_collection(LineCollection(seg, colors=cols, linewidths=widths,
                                     capstyle="round", joinstyle="round",
                                     zorder=3))

    rna_pts = rna.reshape(N_NT, 3, 3)[:, 1, :]             # C4' per nucleotide
    rna_s = smooth(rna_pts, 1)
    rz = (rna_s[:-1, 2] + rna_s[1:, 2]) / 2
    rorder = np.argsort(rz)
    ax.add_collection(LineCollection(
        segs_2d(rna_s)[rorder],
        colors=[depth_shade(RNA_COLOUR, rz[i], lo, hi) for i in rorder],
        linewidths=7.0, capstyle="round", zorder=2))

    allpts = np.vstack([ca_s[:, :2], rna_s[:, :2]])
    c = allpts.mean(axis=0)
    r = float(np.abs(allpts - c).max()) * 1.10 or 1.0
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_png(xyz, ss, label, stack) -> bytes:
    fig, ax = plt.subplots(figsize=(3.1, 2.7), dpi=190)
    draw(ax, xyz, ss)
    ax.set_title(label, fontsize=6.2, color="#333", pad=3)
    fig.patch.set_facecolor("white")
    # A thin coloured rule ties the panel to its generation stack.
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.05)
    fig.add_artist(plt.Line2D([0.02, 0.98], [0.985, 0.985],
                              color=PALETTE.get(stack, "#999"), linewidth=3,
                              transform=fig.transFigure))
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", facecolor="white")
    plt.close(fig)
    return buf.getvalue()


# ---------------------------------------------------------------- deck
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("local_results"))
    ap.add_argument("--build", type=Path, default=Path("build"))
    ap.add_argument("--layer", default="design")
    ap.add_argument("--axis", default="fold")
    ap.add_argument("--resolution", default="cluster_r10")
    ap.add_argument("--min-size", type=int, default=4)
    ap.add_argument("--per-slide", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    base = args.results / args.layer
    nodes = pd.read_csv(base / "dedup" / "nodes.tsv", sep="\t")
    members = pd.read_csv(base / "dedup" / "members.tsv", sep="\t")
    nodes = nodes.merge(members, on="struct_name", how="left")
    nodes["n_designs"] = nodes["n_designs"].fillna(1).astype(int)
    meta = pd.read_parquet(args.build / "meta" / "designs.parquet").set_index("design_id")
    stack = meta.reindex(nodes["struct_name"])["stack"].to_numpy()
    clusters = pd.read_parquet(
        base / "clusters" / args.axis / "clusters.parquet").sort_values("node_id")
    labels = clusters[args.resolution].to_numpy()

    cx_path = base / "complex" / "complex_g02.bin.gz"
    grid = 0.2
    if not cx_path.exists():
        cx_path = base / "complex" / "complex_int16.bin.gz"
        grid = 0.1
    blob = np.frombuffer(gzip.decompress(cx_path.read_bytes()),
                         dtype="<i2").reshape(-1, N_PROT + N_NT * 3, 3)
    sse_all = np.frombuffer(gzip.decompress((base / "complex" / "sse.bin.gz").read_bytes()),
                            dtype=np.uint8).reshape(-1, N_RES)
    print(f"{args.layer}/{args.axis}/{args.resolution}: "
          f"{len(nodes):,} backbones, blob has {len(blob):,} structures")

    sizes = pd.Series(labels).value_counts()
    chosen = [c for c, n in sizes.items() if n >= args.min_size]
    chosen.sort(key=lambda c: -sizes[c])
    print(f"  {len(chosen)} clusters with >= {args.min_size} backbones "
          f"({sizes[chosen].sum():,} of {len(nodes):,} backbones)")

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]

    # ---- summary slide
    s = prs.slides.add_slide(blank)
    tb = s.shapes.add_textbox(Inches(0.6), Inches(0.5), Inches(12.1), Inches(1.0))
    p = tb.text_frame.paragraphs[0]
    p.text = f"{args.layer} layer — {args.axis} clusters ({args.resolution})"
    p.font.size, p.font.bold = Pt(30), True
    body = s.shapes.add_textbox(Inches(0.6), Inches(1.5), Inches(12.1), Inches(5.4))
    tf = body.text_frame
    tf.word_wrap = True
    total_designs = int(nodes["n_designs"].sum())
    facts = [
        f"{len(nodes):,} unique backbones from {total_designs:,} passing designs",
        f"{len(sizes):,} clusters; {int((sizes == 1).sum()):,} singletons",
        f"{len(chosen)} clusters shown (>= {args.min_size} backbones), "
        f"covering {sizes[chosen].sum():,} backbones "
        f"({100 * sizes[chosen].sum() / len(nodes):.1f}%)",
        "",
        "Method: Foldseek all-vs-all, TM = min(qTM, tTM), edges at TM >= 0.40, "
        "top-15 neighbours, Leiden clustering (seed 42).",
        "Every structure is drawn in the RNA target's own frame, so the target "
        "is in the same orientation in every panel and the panels differ only "
        "in how the binder approaches it.",
        "Protein coloured by secondary structure: helix red, strand amber, "
        "coil grey. RNA target in teal. The bar above each panel is its stack.",
    ]
    for i, line in enumerate(facts):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        para.text = line
        para.font.size = Pt(15)
        para.font.color.rgb = RGBColor(0x33, 0x33, 0x33)

    # ---- one slide per cluster
    for rank, cid in enumerate(chosen, start=1):
        idx = np.flatnonzero(labels == cid)
        n_bb = len(idx)
        n_des = int(nodes["n_designs"].to_numpy()[idx].sum())
        counts = pd.Series(stack[idx]).value_counts()

        slide = prs.slides.add_slide(blank)
        tb = slide.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.7))
        p = tb.text_frame.paragraphs[0]
        p.text = f"{args.axis.capitalize()} cluster {int(cid)}"
        p.font.size, p.font.bold = Pt(28), True

        sub = slide.shapes.add_textbox(Inches(0.5), Inches(0.85), Inches(12.3), Inches(0.45))
        sp = sub.text_frame.paragraphs[0]
        sp.text = (f"#{rank} of {len(chosen)}   |   {n_bb:,} backbones   |   "
                   f"{n_des:,} designs   |   "
                   f"{100 * n_bb / len(nodes):.2f}% of the layer")
        sp.font.size = Pt(14)
        sp.font.color.rgb = RGBColor(0x44, 0x44, 0x44)

        # stack breakdown, percentages of this cluster
        y = Inches(1.35)
        for st in STACKS:
            c = int(counts.get(st, 0))
            if not c:
                continue
            row = slide.shapes.add_textbox(Inches(9.9), y, Inches(3.0), Inches(0.32))
            rp = row.text_frame.paragraphs[0]
            rp.text = f"{st}   {100 * c / n_bb:5.1f}%   ({c:,})"
            rp.font.size = Pt(12)
            hexcol = PALETTE[st].lstrip("#")
            rp.font.color.rgb = RGBColor(int(hexcol[0:2], 16), int(hexcol[2:4], 16),
                                         int(hexcol[4:6], 16))
            rp.font.bold = True
            y = Emu(int(y) + Inches(0.34))

        # up to per-slide structures, evenly spaced through the cluster
        have = [i for i in idx if i < len(blob)]
        pick = ([have[k] for k in np.linspace(0, len(have) - 1,
                                              min(args.per_slide, len(have))).astype(int)]
                if have else [])
        positions = [(Inches(0.55), Inches(1.45)), (Inches(5.0), Inches(1.45)),
                     (Inches(0.55), Inches(4.35)), (Inches(5.0), Inches(4.35))]
        for slot, node in enumerate(pick[:4]):
            xyz = decode(blob, node, grid)
            ss = "".join(chr(c) for c in sse_all[node])
            png = panel_png(xyz, ss, nodes["struct_name"].iloc[node][:52],
                            stack[node])
            left, top = positions[slot]
            slide.shapes.add_picture(_io.BytesIO(png), left, top,
                                     width=Inches(4.25))
        if not pick:
            note = slide.shapes.add_textbox(Inches(0.55), Inches(3.2),
                                            Inches(8.5), Inches(0.5))
            note.text_frame.paragraphs[0].text = (
                "No structures embedded for this cluster.")
        if rank % 10 == 0:
            print(f"    {rank}/{len(chosen)} slides")

    out = args.out or Path(
        f"dist/clusters_{args.layer}_{args.axis}_{args.resolution}.pptx")
    out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(out)
    print(f"  wrote {out} ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
