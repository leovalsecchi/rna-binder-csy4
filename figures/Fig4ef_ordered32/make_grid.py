#!/usr/bin/env python3
"""Fig 3c: the 32 ordered candidates as one grid, all in the Fig 2b shared
RNA frame at the same rotation (y120) as the Fig 3b panels.

Panels come from Fig2b_stack_binders/make_figure.py --set hits32, so the
hairpin is identical in every cell and only the binder moves. Cells are
grouped by stack in the fixed thesis stack order, then by ON/OFF rank, so the
grid reads against the pie in the same panel.

    python3 code/Fig2b_stack_binders/make_figure.py PNGDIR \
        --set hits32 --width 1200 --frames 12 --angle 120
    python3 code/Fig4ef_ordered32/make_grid.py PNGDIR
"""
import argparse
import csv
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402

HITS = Path(__file__).resolve().parents[2] / "working/structures/hits32"
OUT = Path(__file__).resolve().parents[2] / "output/Fig4ef_ordered32"
SCREEN = Path(__file__).resolve().parents[2] / "inputs/screen_points.csv"
FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"


def splicing():
    """design -> splicing ratio, from the same table Fig 4b plots.

    Two joins to get right. `transcript` is a composite
    `POS_A_vs_NEG_A::POS_A::REF_<hash>::<design>`, so the design id is its
    last `::` field; and index.tsv names the comparison `POS_A` where this
    table names it `POS_A_vs_NEG_A`. Keying on (design, comparison) and not on
    design alone matters - a design screened in both comparisons has two rows
    with different ratios, and index.tsv says which one this hit is.
    """
    out = {}
    with SCREEN.open() as fh:
        for r in csv.DictReader(fh):
            out[(r["transcript"].rsplit("::", 1)[-1], r["comparison"])] = \
                float(r["splicing_ratio"])
    return out


def stack_of(name):
    for short in palette.STACKS_SHORT:
        if f"_{short.replace('+', '_')}_" in name:
            return short
    raise ValueError(name)


def rows():
    out = []
    spl = splicing()
    for r in HITS.joinpath("index.tsv").read_text().splitlines()[1:]:
        if not r.strip():
            continue
        f, cmp_, ratio = r.split("\t")
        rank, design = f.split("_", 1)[0], f.split("_", 1)[1][:-4]
        key = (design, f"{cmp_}_vs_NEG_{cmp_[-1]}")
        assert key in spl, f"no splicing ratio for {key}"
        out.append(dict(rank=rank, design=design, stack=stack_of(design),
                        cmp=cmp_, ratio=float(ratio), spl=spl[key]))
    order = {palette.SHORT[k]: i for i, k in enumerate(palette.STACK_ORDER)}
    out.sort(key=lambda d: (order[d["stack"]], -d["ratio"]))
    return out


def hexrgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def tint(h, f):
    """Blend a stack colour f of the way toward white."""
    return tuple(round(c + (255 - c) * f) for c in hexrgb(h))


def fit(text, px, size, font_path=FONT):
    """Largest font at or below `size` whose text fits in `px`."""
    while size > 7:
        f = ImageFont.truetype(font_path, size)
        if f.getlength(text) <= px:
            return f
        size -= 1
    return ImageFont.truetype(font_path, 7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("png_dir")
    # 7 x 5 since 2026-11-14, not 8 x 4: at 178 mm the cells come out 25.4 mm
    # instead of 22.2 and the sheet 123 mm instead of 86, which is most of the
    # white the page had left under it. 35 slots for 32, so the last row runs
    # three short.
    ap.add_argument("--rows", type=int, default=5)
    ap.add_argument("--cols", type=int, default=7)
    ap.add_argument("--angle", type=int, default=120)
    ap.add_argument("--cell", type=int, default=600, help="cell width in px")
    ap.add_argument("--margin", type=float, default=0.06,
                    help="breathing room added around the union ink box, "
                         "as a fraction of its larger side")
    ap.add_argument("--boxes", action="store_true",
                    help="tint each stack's cells as one connected shape")
    ap.add_argument("--names", action="store_true",
                    help="print each design's name under its cell. Off since "
                         "2026-11-14: at 8 columns it set at about 3 pt, half "
                         "the house size, and said nothing a reader could use "
                         "- the identifiers belong in a supplementary table")
    ap.add_argument("--label-pt", type=float, default=5.5,
                    help="printed size of the ON/OFF label, in pt, ASSUMING "
                         "the sheet is placed at --page-mm wide")
    ap.add_argument("--page-mm", type=float, default=178.0)
    a = ap.parse_args()

    recs = rows()
    assert a.rows * a.cols >= len(recs), f"{a.rows}x{a.cols} cannot hold {len(recs)}"
    # 2026-11-18 (Leo, from his supervisor): "fill up to the very bottom with
    # designs, but have this wide gap that we have at the bottom right now at
    # the top left." Every panel is pushed back by the number of free slots,
    # so the hole opens at the START of the grid instead of the end and the
    # key moves into it. Slot index is the only thing that changes - the
    # records stay in stack-then-rank order, so the reading order is intact.
    OFF = a.rows * a.cols - len(recs)
    src = Path(a.png_dir)

    # One window SIZE for all 32, placed per panel on that panel's own ink.
    # A shared window POSITION - the union of all 32 ink boxes, which is what
    # this did until 2026-11-14 - is as wide as the widest binder, because the
    # hairpin sits at a fixed place in every render and the binder extends to
    # one side only: every compact panel then sat left in its cell with dead
    # space to its right, which against the stack tints read as a misplaced
    # rendering (Leo: "not in the center of the squares that color them").
    # Centring each panel costs the shared ink position; the shared camera,
    # rotation and scale, which is what lets the binders be compared, stay.
    src_imgs, inks = {}, []
    for r in recs:
        im = Image.open(src / f"{r['rank']}_y{a.angle:03d}.png").convert("RGBA")
        src_imgs[r["rank"]] = im
        inks.append(im.getchannel("A").getbbox() or im.getbbox())
    # Size: the union of all 32 ink boxes, grown by --margin. Any single
    # panel's ink fits inside a window that size, and keeping the size the
    # union gives keeps the sheet's aspect, so the page geometry under it is
    # unchanged by the centring above.
    u = [f(b[i] for b in inks) for i, f in enumerate((min, min, max, max))]
    src_w, src_h = next(iter(src_imgs.values())).size
    grow = round(a.margin * max(u[2] - u[0], u[3] - u[1]))
    bw = min(src_w, u[2] + grow) - max(0, u[0] - grow)
    bh = min(src_h, u[3] + grow) - max(0, u[1] - grow)
    cell_h = round(a.cell * bh / bw)

    # Strict slot geometry. Every rendering owns one slot of the same size; a
    # slot is the rendering plus its caption plus half a gutter on each side.
    # The stack shapes are slots inset by the SAME amount everywhere, so each
    # boundary falls exactly midway between two renderings and every margin,
    # inner and outer, is equal.
    pad = a.cell // 18                   # full gutter between renderings
    inset = pad // 4                     # slot -> shape inset, uniform
    s1 = max(14, a.cell // 20)           # filename, stack colour
    # The ON/OFF label is sized from where it LANDS, not from the cell: a
    # fraction of the cell printed at ~3 pt once the sheet was scaled to the
    # page. mm per px is page_mm / sheet_px, so the size that prints at
    # label_pt is that, inverted, times the pt in a mm.
    s2 = round(a.label_pt * 25.4 / 72 * a.cols * (a.cell + a.cell // 18)
               / a.page_mm)
    # Corner glyphs, all sized off the cell so they hold at any --cell.
    # 2026-11-18 (Leo, from his supervisor): "for each binder candidate ... in
    # the corner ... one circle indicating the fold change and one kind of bar
    # that adds up to one indicating the splicing ratio".
    # 2026-11-18 (Leo): "align the circle and bar better ... have the bar in
    # the very corner and have the ball basically sitting on the bar ... but
    # then you also add the fold change right next to it." So gm halves (the
    # bar really is in the corner now), every ball SITS on one baseline just
    # above the bar whatever its size, and the fold-change number moves off
    # the caption row to sit on that same baseline beside its ball.
    # 2026-11-18 (Leo): "move everything more into the lower right corner".
    gm = round(pad * 0.18)                      # inset from the cell's edges
    bar_w = round(a.cell * 0.26)
    bar_h = max(2, round(a.cell * 0.022))
    btick = max(1, round(bar_h * 0.75))         # end ticks: the 0 and 1 marks
    tw = max(1, round(a.cell / 400))
    r_max = round(a.cell * 0.072)               # radius of the LARGEST circle
    gapn = round(a.cell * 0.022)                # number to ball
    fc_max = max(r["ratio"] for r in recs)
    fc_min = min(r["ratio"] for r in recs)
    lead = max(4, pad // 6)
    # The ratio used to be set on its own row under each cell; it now rides
    # in the corner next to its ball, so that row is no longer reserved. Five
    # rows x s2 is about 7% of the sheet height given back to the renderings.
    lab = pad // 2 + (s1 + lead if a.names else 0) + pad // 2
    slot_w, slot_h = a.cell + pad, cell_h + lab + pad
    W = a.cols * slot_w + 2 * inset
    H = a.rows * slot_h + 2 * inset
    sheet = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(sheet)

    def slot(k):
        col, row = k % a.cols, k // a.cols
        return inset + col * slot_w, inset + row * slot_h

    def cellbox(k):
        sx, sy = slot(k)
        return sx + pad // 2, sy + pad // 2

    # One tinted shape per stack, drawn first so the panels sit on it. A stack
    # can wrap across rows; the row-runs are bridged over their horizontal
    # overlap so each stack is ONE connected shape, not separate boxes.
    runs = [] if a.boxes else []         # (stack, first_k, last_k)
    for k, r in enumerate(recs if a.boxes else []):
        kk = k + OFF
        if runs and runs[-1][0] == r["stack"] and kk % a.cols != 0:
            runs[-1][2] = kk
        else:
            runs.append([r["stack"], kk, kk])

    def runrect(k0, k1):
        x0, y0 = slot(k0)
        x1, _ = slot(k1)
        return [x0 + inset, y0 + inset, x1 + slot_w - inset, y0 + slot_h - inset]

    mask = Image.new("L", (W, H), 0)
    for i, (stack, k0, k1) in enumerate(runs):
        if stack is None:                          # already drawn with its run group
            continue
        mask.paste(0, (0, 0, W, H))
        md = ImageDraw.Draw(mask)
        rects = [runrect(k0, k1)]
        j = i
        while j + 1 < len(runs) and runs[j + 1][0] == stack:
            j += 1
            rects.append(runrect(runs[j][1], runs[j][2]))
        for rc in rects:
            md.rounded_rectangle(rc, radius=pad // 2, fill=255)
        for ra, rb in zip(rects, rects[1:]):      # bridge consecutive rows
            xa, xb = max(ra[0], rb[0]), min(ra[2], rb[2])
            if xb > xa:
                md.rectangle([xa, ra[3] - pad // 2, xb, rb[1] + pad // 2],
                             fill=255)
        sheet.paste(Image.new("RGB", (W, H), tint(palette.stack_color(stack), 0.80)),
                    (0, 0), mask)
        runs[i][0] = None                          # drawn with its siblings
        for j2 in range(i + 1, len(runs)):
            if runs[j2][0] == stack:
                runs[j2][0] = None

    f2 = ImageFont.truetype(FONT, s2)
    # The widest number plus its ball plus both insets has to sit inside the
    # cell, or a long label would run under the neighbouring candidate.
    _wide = max(f2.getlength(f"{r['ratio']:.0f}\u00d7") for r in recs)
    assert _wide + gapn + 2 * r_max + gm <= a.cell, (
        f"fold-change label needs {_wide + gapn + 2 * r_max + gm:.0f}px of a "
        f"{a.cell}px cell - raise --cell or drop --label-pt")

    def glyphs_corner(d, x, y, r, c):
        """Fold change and splicing ratio in the cell's bottom-right corner.

        Two different encodings on purpose. Fold change spans 12x to 170x, so
        it is a circle whose AREA is proportional to the value (radius goes as
        the square root) - a 14-fold range of radius would not fit the corner,
        and area is the honest encoding for a disc anyway. Splicing ratio is a
        fraction of a whole, so it is a bar against its full 0-1 extent, with
        a tick at each end to say where that extent ends. Without the ticks a
        short bar could be read as a short axis rather than a small fraction.

        The NUMBER and the ball encode the same quantity, so they are one
        mark: same colour, same baseline, side by side. The ball's CENTRE is
        fixed and only its radius moves, which is what lets a reader compare
        balls across candidates; its BOTTOM sits on the baseline, so balls of
        every size rest on the bar instead of floating at different heights.
        """
        gx1 = x + a.cell - gm                   # right edge, shared by both
        by1 = y + cell_h - gm
        by0 = by1 - bar_h
        d.rounded_rectangle([gx1 - bar_w, by0, gx1, by1], radius=bar_h / 2,
                            fill=tint(palette.stack_color(r["stack"]), 0.55))
        # A floor of one bar-height so the smallest ratio is still a mark and
        # not an invisible sliver; at 0.072 the fill is ~2% of the track.
        fill_w = max(bar_h, round(bar_w * r["spl"]))
        d.rounded_rectangle([gx1 - bar_w, by0, gx1 - bar_w + fill_w, by1],
                            radius=bar_h / 2, fill=c)
        for ex in (gx1 - bar_w, gx1):
            d.line([ex, by0 - btick, ex, by1 + btick], fill=c, width=tw)
        # 0.6 of a bar-height of air under the ball, not 0.3: the number is
        # now centred on the ball's middle, so for the SMALLEST ball (12x,
        # whose radius is a quarter of the largest) the text hangs lower than
        # the ball does and would otherwise graze the bar's end ticks.
        base = by0 - btick - max(1, round(bar_h * 0.6))
        rr = r_max * (r["ratio"] / fc_max) ** 0.5
        ccx = gx1 - r_max
        d.ellipse([ccx - rr, base - 2 * rr, ccx + rr, base], fill=c)
        # 2026-11-18 (Leo): "align the number on the horizontal middle height
        # of the ball" - so the number tracks its own ball's centre rather
        # than sitting on the shared baseline. It therefore rides up and down
        # with fold change, which is the point: number and ball are one mark.
        d.text((ccx - r_max - gapn, base - rr), f"{r['ratio']:.0f}\u00d7",
               font=f2, fill=c, anchor="rm")

    for k, r in enumerate(recs):
        # Window centred on this panel's ink. It may run past the source edge;
        # the overhang crops to transparent and composites to the cell.
        b = inks[k]
        x0, y0 = round((b[0] + b[2] - bw) / 2), round((b[1] + b[3] - bh) / 2)
        im = src_imgs[r["rank"]].crop((x0, y0, x0 + bw, y0 + bh))
        ib = im.getchannel("A").getbbox()
        assert ib and (ib[2] - ib[0], ib[3] - ib[1]) == (b[2] - b[0], b[3] - b[1]), \
            f"panel {r['rank']}: the window clips its ink"
        x, y = cellbox(k + OFF)
        bg = tint(palette.stack_color(r["stack"]), 0.80) if a.boxes \
            else (255, 255, 255)
        cell = Image.new("RGBA", im.size, bg + (255,))
        cell.alpha_composite(im)
        cell = cell.convert("RGB").resize((a.cell, cell_h), Image.LANCZOS)
        sheet.paste(cell, (x, y))
        c = hexrgb(palette.stack_color(r["stack"]))
        ly = y + cell_h + pad // 2
        if a.names:
            # drop the constant csy4free_ prefix so the name can be set larger
            name = r["design"].split("_len80_")[1]
            d.text((x, ly), name, font=fit(name, a.cell - 4, s1), fill=c)
            ly += s1 + lead
        glyphs_corner(d, x, y, r, c)

    # Keys, in the free slots the 32 panels leave at the TOP LEFT (see OFF).
    # Three things a reader cannot decode from the sheet itself: the stack
    # tints, and the two corner glyphs. They sit where the eye starts rather
    # than where it ends, which is what the move was for. Frameless and at the
    # label size (STYLE_NOTES section 6); the swatches carry a hairline edge
    # because an 0.80 tint is under 3:1 against white (COLORS rule 9).
    free = OFF
    if a.boxes and free:
        fx, fy = slot(0)
        f3 = ImageFont.truetype(FONT, s2)
        sw_h, gap = round(s2 * 1.3), round(s2 * 0.5)
        sw_w, pitch = round(s2 * 1.3 * 2.2), round(s2 * 1.3 * 1.55)
        names = [palette.SHORT[k] for k in palette.STACK_ORDER]
        block_h = pitch * len(names) - (pitch - sw_h)
        # Left-aligned on the first free slot's cell, so the swatches sit on
        # the same column grid as the renderings above them rather than
        # floating in the middle of the hole.
        x = fx + pad // 2
        y = round(fy + (slot_h - block_h) / 2)
        for i, (long, short) in enumerate(zip(palette.STACK_ORDER, names)):
            c = hexrgb(palette.stack_color(long))
            top = y + i * pitch
            d.rounded_rectangle([x, top, x + sw_w, top + sw_h],
                                radius=pad // 2,
                                fill=tint(palette.stack_color(long), 0.80),
                                outline=c, width=max(1, a.cell // 400))
            d.text((x + sw_w + gap, top + sw_h / 2), short, font=f3, fill=c,
                   anchor="lm")

        # The two corner-glyph keys, in the slots beside the swatches. Both
        # are NEUTRAL, not a stack colour: they explain a mark that every
        # stack carries, so tinting them with one would read as naming it.
        grey, inkc = hexrgb(palette.NEUTRAL["600"]), hexrgb(palette.INK)
        # The shared CENTRELINE, not a baseline: both keys are centred on it, so
        # the balls and the splicing track read as one row across the slots
        # (Leo, 2026-11-18: "align these too"). The bar used to hang below cl by
        # a tick plus a third of its height, which left the two keys on visibly
        # different lines once the balls were centred.
        cl = y + round(block_h * 0.52)
        kb_top = cl - bar_h // 2
        assert abs(kb_top + bar_h / 2 - cl) <= 1, "splicing bar is off the shared centreline"
        if free >= 2:
            kx = slot(1)[0] + pad // 2
            d.text((kx, y), "Fold change", font=f3, fill=inkc, anchor="lt")
            # The two reference circles are the RANGE'S ENDS, not a round
            # pair: with an area encoding the reader needs the extremes to
            # calibrate against, and these are the real 32's extremes.
            #
            # Number then ball, everything CENTRED on cl. Sitting the balls on
            # a common bottom instead puts their centres - and so their
            # numbers, which track their own ball - at two different heights,
            # and the pair reads as misaligned (Leo, 2026-11-18: "alignment").
            # A cell draws one ball and has nothing to line up against; a key
            # draws two side by side, so here a shared centreline wins.
            assert cl + r_max <= fy + slot_h, "the fold-change key overflows its slot"
            cur = kx
            for v in (fc_max, fc_min):
                rv = r_max * (v / fc_max) ** 0.5
                lab_ = f"{v:.0f}\u00d7"
                cur += f3.getlength(lab_)
                d.text((cur, cl), lab_, font=f3, fill=grey, anchor="rm")
                cx_ = cur + gapn + r_max
                d.ellipse([cx_ - rv, cl - rv, cx_ + rv, cl + rv], fill=grey)
                cur = cx_ + r_max + round(s2 * 1.1)
        if free >= 3:
            kx = slot(2)[0] + pad // 2
            d.text((kx, y), "Splicing ratio", font=f3, fill=inkc, anchor="lt")
            kw = round(bar_w * 1.5)
            d.rounded_rectangle([kx, kb_top, kx + kw, kb_top + bar_h],
                                radius=bar_h / 2,
                                fill=hexrgb(palette.NEUTRAL["300"]))
            # Filled to the 32's MEDIAN, so the key shows what a typical bar
            # looks like as well as what the track means. Labelling the ends
            # 0 and 1 is the whole point of the glyph: these designs splice a
            # tenth of their transcript, and a bar that stopped at the data's
            # own maximum would hide that.
            med = sorted(r["spl"] for r in recs)[len(recs) // 2]
            d.rounded_rectangle([kx, kb_top,
                                 kx + max(bar_h, round(kw * med)),
                                 kb_top + bar_h],
                                radius=bar_h / 2, fill=grey)
            for ex in (kx, kx + kw):
                d.line([ex, kb_top - btick, ex, kb_top + bar_h + btick],
                       fill=grey, width=tw)
            for ex, lab in ((kx, "0"), (kx + kw, "1")):
                d.text((ex, kb_top + bar_h + btick + gap), lab, font=f3,
                       fill=grey, anchor="mt")

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / (f"Fig4f_grid_{a.rows}x{a.cols}_y{a.angle:03d}"
                 + ("_boxed" if a.boxes else "")
                 + ("_named" if a.names else "") + ".png")
    sheet.save(out)
    print("wrote", out, sheet.size)


if __name__ == "__main__":
    main()
