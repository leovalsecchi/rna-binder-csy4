"""Bring Leo's `Fluorescence Splicing V1.svg` into house style.

His Illustrator file is the layout; this only corrects it. Geometry is left
alone except the grey container, which is redrawn as a clean rounded rect.

    python code/Fig1b_sms_assay/fix_fluorescence.py

in  working/Fig1d_integrated_construct/Fluorescence Splicing V1.svg
out working/Fig1d_integrated_construct/Fig1b_fluorescence_splicing.{svg,pdf,png}

## What it changes and why

**One stroke weight.** His file carries two: the imported art had no
`stroke-width` at all, so it fell back to 1, while the construct blocks were
1.76. Neither is right. His blocks are 10.95 deep where Matteo's are 6, so the
file sits at 1.825x Matteo's scale and his 0.4 hairline becomes **0.73** here.
Everything strokes at that, so the panel reads as one drawing.

**Palette.** Three fills were off-palette placeholders: `#d4145a` for mScarlet,
`#22b573` for mNeonGreen, `#e6e6e6` for the container. Pure black in two text
classes becomes `ink`.

**Type.** Three sizes in two families, and every text is drawn through a
`matrix(1.193 ...)`, so the size on the page is 1.193x the number in the file.
All of it becomes Arial at one size; the in-block captions are smaller, in the
same ratio to their block that the other panel uses.

**The reporters.** His coloured squares are placeholders. They become Matteo's
GFP barrel - green as mNeonGreen, recoloured as mScarlet.

**Matching the other panel.** `Fig1b_panel_b.svg` draws 6-deep blocks with a
0.4 hairline and 5.5 type. This file is the same drawing at 1.825x, so the
numbers below are those multiplied through. Place the two at the same block
depth and they agree.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cartoon, palette                           # noqa: E402
from common.paths import ROOT                                 # noqa: E402

WORK = ROOT / "working/Fig1d_integrated_construct"
SRC = WORK / "Fluorescence Splicing V1.svg"
DST = WORK / "Fig1b_fluorescence_splicing.svg"

S, N = palette.SYNBIO, palette.NEUTRAL
INK = palette.INK

K = 1.825                       # this file's scale against Matteo's own
STROKE = round(0.4 * K, 2)      # 0.73 - his hairline, brought up to this scale
TYPE = round(5.5 * K, 1)        # 10.0 on the page
MTX = 1.193                     # the scale every one of his texts is drawn through
IN_BLOCK = round(4.4 * K, 1)    # 8.0 - a caption inside a block
RADIUS = round(2.0 * K, 1)

# his placeholder squares -> Matteo's GFP barrel, at the square's centre.
# It is the only fluorescent-protein glyph in his set, so it is the one to
# use; at 1.5 with a neutral-700 outline it read as a tin can, so it is drawn
# nearer his own size with his own lighter rim, and it shines.
BARREL_K = 1.2
BARREL_W, BARREL_H = 9.578, 14.789        # the glyph at 1x
RAY_N = 8
RAY_GAP = round(1.6 * K / 1.825, 2)       # 1.6 - clear of the rim
RAY_LEN = round(4.0 * K / 1.825, 2)       # 4.0
HALO = 0.1                                # a breath of colour behind it
# the red one is lifted 4.5 off his placeholder so his translation arrow
# stops short of the rays instead of running into them
BARRELS = [((11.56 + 22.51) / 2, (13.47 + 24.42) / 2, S["reporter"]),
           ((304.06 + 315.01) / 2, (136.21 + 147.16) / 2 - 4.5, S["reporter_alt"])]
# the lens between his dome and his rim arc, lit: a cylinder with no lit top
# reads flat, which is most of what made it look like a can
CAP_LIT = ("M31.211111 11.605719 C31.281536 9.915523 33.535131 9 36 9 "
           "C38.605719 9 40.788889 10.056373 40.788889 11.605719 "
           "L40.507189 11.18317 C39.802941 13.859314 32.197059 13.859314 "
           "31.49281 11.18317 Z")

# --- what each element IS -------------------------------------------------
# Leo built the bars by duplicating objects, so several ids still name what
# was copied. Three of them are keys in retheme_integration.py's FILL_BY_ID,
# which colours BY ID: left alone, a re-theme paints the mScarlet CDS
# reporter-green and the grey container theme-200. Names follow Matteo's own
# (Promoterless-donor-plasmid.svg) where he has one.
RENAME = {
    # the grey container was a copy of the scaffold block
    "RNA_scaffold": ("Component_box", "Component box"),
    # bar 1: the green block is the reporter, the pale blue one the target
    "RNA_scaffold-3": ("mNeonGreen_replacement_cassette",
                       "mNeonGreen replacement cassette"),
    "RNA_scaffold-2": ("RNA_scaffold", "Target array"),
    # bar 3: the red CDS was copied from the landing pad's green one
    "mNeonGreen_replacement_cassette": ("mScarlet_encoding", "mScarlet | encoding"),
    # the ids the fused receptor ships with (retheme_integration._fuse_receptor);
    # his data-names already match, only the ids differ
    "SMS_kinase_lobe": ("SMS_kinase", "SMS kinase lobe"),
    "SMS_RNase_lobe": ("SMS_rnase", "SMS RNase lobe"),
}
# blocks Illustrator left without a data-name
NAME_ONLY = {"Binder": "Binder", "Linker": "Linker",
             "Kinase": "Kinase", "RNase": "RNase"}

# --- the bars all span the receptor -----------------------------------------
# Bar 2 is one polypeptide at Fig1d's proportions, so its span is fixed; the
# other two are now cut to it. Only the GAPS absorb the difference, so no
# block's width changes except the AT one, which has to hold its caption.
AT_IDS = ("_5_prime_segment", "Intron_segment", "_3_prime_segment",
          "Hairpin_encoding_outline")
AT_DX = -0.892                  # his AT block started right of the other bars
SPAN = 115.842                  # bar 2, the receptor: every bar matches it
BLOCK_W = 49.438                # bar 1's two blocks
MSC_W = 60.230                  # bar 3's mScarlet CDS
GAP = round(SPAN - 2 * BLOCK_W, 3)        # 16.966
AT_W = round(SPAN - MSC_W - GAP, 3)       # 38.646

# --- the component box, set out from its own content ------------------------
# His was 230 across for 152 of content, with the numbers 26 from their bars
# and a 22.0 lead stub on the left against 5.6 on the right. Laid out from the
# content instead: one margin all round, equal stubs, the numbers up against
# the bars. The right edge stays where he put it, level with the strand below.
MARGIN = 15.85                  # = the white between two rows
STUB = 9.0                      # lead either side of the blocks, equal now
NUM_GAP = 9.0                   # ring to the start of its lead
BOX_RIGHT = 391.3
LEAD_X1 = round(BOX_RIGHT - MARGIN, 3)
BAR_X1 = round(LEAD_X1 - STUB, 3)
BAR_X0 = round(BAR_X1 - SPAN, 3)
LEAD_X0 = round(BAR_X0 - STUB, 3)
SHIFT = round(BAR_X0 - 245.640, 3)        # how far his blocks travel
# His AT block was 21.892 across, so its 8.756 intron segment was narrower
# than the "AT" on top of it: the caption ran onto the blue exons either side.
# At 38.646 the intron is 15.46 and the caption sits inside it.
AT_SPLIT = (0.30002, 0.69998)   # his own exon / intron boundaries, kept
AT_R = 1.460                    # Binder's corner radius. His AT block had 0.876
                                # horizontally - it had been squashed in x.
AT_K = 0.5523                   # circular-arc control offset
DX = {"RNA_scaffold": round(BAR_X1 - 357.636, 3),        # bar 1's second block
      "mScarlet_encoding": round(BAR_X0 + AT_W + GAP - 277.348, 3)}
ROW_TOP, ROW_BOT = 21.804, 86.606         # the blocks' own extent

# --- the callouts ----------------------------------------------------------
# The box numbers were orphans: nothing in the panel said what 1, 2 and 3
# referred to. Each now has a twin on the part it names.
CALLOUTS = [(73.5, 17.5, "1"),      # the RNA target in its scaffold
            (74.0, 95.0, "2"),      # clear of the linker curve, which
                                    # reaches x 72.2 down to y 89.2
            (91.6, 178.0, "3")]     # between his cut-site arrowhead (ends
                                    # 84.9) and the mScarlet bracket (98.09),
                                    # level with that bracket
RING = round(2.4 * K, 1)            # 4.4
DIGIT = round(3.6 * K, 1)           # 6.6 - a digit in a ring, not a caption
CAP = 0.716                         # Arial cap height, for centring a digit
NUM_X = round(LEAD_X0 - NUM_GAP - RING, 3)
BOX = (round(NUM_X - RING - MARGIN, 3), round(ROW_TOP - MARGIN, 3),
       BOX_RIGHT, round(ROW_BOT + MARGIN, 3))

# --- start and stop --------------------------------------------------------
# Without them the Frameshift arrow has nothing to act on. Same idiom as
# Fig1d: a tick off the strand, the codon under it.
TICK_W = round(0.35 * K, 2)         # 0.64 - Fig1d's tick against its 0.4 line
STRAND_Y = 189.384
TICK0, TICK1 = STRAND_Y + 1.6, STRAND_Y + 5.2   # 3.6, Fig1d's 2.1 x 1.825
CODON_Y = STRAND_Y + 10.0
# construct_key.svg is the key for these: aug -> "AUG", stop -> "*",
# mscarlet_off -> "mScarlet*". "STOP" is not a name we use anywhere.
STAR_SIZE = round(IN_BLOCK * 1.5, 1)  # 12.0 - an asterisk is drawn small, so
                                      # at the caption size it reads as a speck
STAR_DROP = round(STAR_SIZE * 0.232, 2)  # and it sits high: drop it onto the
                                         # band the other captions use
CODONS = [(29.95, "AUG"), (88.10, "*"),       # unspliced: stop BEFORE the CDS
          (264.30, "AUG"), (390.60, "*")]     # spliced:   stop AFTER it
PAGE_H = 202.5

# The blocks in the component box, captioned. Keyed on the block, so a
# caption follows its block when the block moves - his were absolute, and
# stayed behind when the bars were snapped onto their lead lines.
# "SMS" goes on the kinase lobe, not the outline: centred on the outline it
# straddles the kinase/RNase boundary. "AT" goes on the intron segment, so it
# does not run onto the blue exons either side.
CAPTIONS = [("mNeonGreen_replacement_cassette", "mNeonGreen", INK),
            ("RNA_scaffold", "RNA scaffold", INK),
            ("Binder", "Binder", palette.SURFACE),
            ("SMS_kinase", "SMS", INK),
            ("Intron_segment", "AT", palette.SURFACE),
            ("mScarlet_encoding", "mScarlet", palette.SURFACE)]

COLOURS = {"#d4145a": S["reporter_alt"],      # mScarlet placeholder
           "#22b573": S["reporter"],          # mNeonGreen placeholder
           "#e6e6e6": N["100"],               # container fill
           "#000000": INK}



def _bbox(tag):
    d = re.search(r'\sd="([^"]+)"', tag)
    if not d:
        return None
    n = [float(v) for v in re.findall(r"-?[\d.]+", d.group(1))]
    return (min(n[0::2]), min(n[1::2]), max(n[0::2]), max(n[1::2])) if n else None


def rename(t):
    """Give every block the name of the thing it draws.

    Keyed on the id, and the data-name is rewritten with it - his ids and his
    data-names disagree (id="RNA_scaffold-3" data-name="RNA scaffold" on the
    green mNeonGreen block), so renaming by text alone misses half of them.
    """
    def one(m):
        tag = m.group(0)
        eid = re.search(r'\sid="([^"]+)"', tag)
        if not eid or eid.group(1) not in RENAME:
            return tag
        new_id, new_nm = RENAME[eid.group(1)]
        tag = tag.replace(f'id="{eid.group(1)}"', f'id="{new_id}"', 1)
        return re.sub(r'\sdata-name="[^"]*"', f' data-name="{new_nm}"', tag, count=1)
    t = re.sub(r'<(?:path|rect|g)\b[^>]*>', one, t)
    # Illustrator drops data-name on the first of a duplicated pair
    def fill_in(m):
        tag, eid = m.group(0), m.group(1)
        if "data-name=" in tag or eid not in NAME_ONLY:
            return tag
        return tag.replace(f'id="{eid}"', f'id="{eid}" data-name="{NAME_ONLY[eid]}"', 1)
    t = re.sub(r'<(?:path|rect|g)\b[^>]*\sid="([^"]+)"[^>]*>', fill_in, t)
    # his two strand CDSs carry no name at all
    for tag in re.findall(r'<path[^>]*class="cls-15"[^>]*/>', t):
        nm = "unspliced" if _bbox(tag)[0] < 200 else "spliced"
        t = t.replace(tag, tag.replace(
            "<path", f'<path id="mScarlet_encoding_{nm}" '
                     f'data-name="mScarlet | encoding, {nm}"', 1), 1)
    return t


def snap_rows(t):
    """Put every block back on its own lead line.

    His three leads are an even grid (27.277 / 54.071 / 81.133); the blocks
    had drifted off it by up to 1.13, and bar 1's two blocks were 0.72 apart
    from each other, so the line ran through them at three different heights.
    The leads are the grid, so the blocks move, not the lines.
    """
    leads = sorted({y for tag in re.findall(r'<path[^>]*class="cls-13"[^>]*/>', t)
                    for y in [float(re.findall(r"-?[\d.]+",
                                   re.search(r'\sd="([^"]+)"', tag).group(1))[1])]
                    if y < 100})
    moved = []
    for m in list(re.finditer(r'<(path|rect)\b[^>]*/>', t)):
        tag = m.group(0)
        b = _bbox(tag)
        if not b or b[0] < 220 or b[3] > 100:          # inside the box only
            continue
        cy = (b[1] + b[3]) / 2
        lead = min(leads, key=lambda y: abs(y - cy))
        if 'class="cls-13"' in tag:            # a lead line, redrawn below
            continue
        eid = re.search(r'\sid="([^"]+)"', tag)
        name = eid.group(1) if eid else ""
        dx = AT_DX if name in AT_IDS else DX.get(name, SHIFT)
        dy = lead - cy
        if abs(dx) < 1e-4 and abs(dy) < 1e-4:
            continue
        t = t.replace(tag, tag.replace("<path", f'<path transform="translate({dx:g} {dy:g})"',
                                       1).replace("<rect", f'<rect transform="translate({dx:g} {dy:g})"', 1), 1)
        moved.append((eid.group(1) if eid else "-", dx, dy))
    return t, leads, moved


def relead(t, leads):
    """Equal stubs either side of the blocks, in the box's new position."""
    for tag in re.findall(r'<path[^>]*class="cls-13"[^>]*/>', t):
        y = float(re.findall(r"-?[\d.]+",
                             re.search(r'\sd="([^"]+)"', tag).group(1))[1])
        if y > 100:                                   # the two mRNA strands
            continue
        t = t.replace(tag, re.sub(r'\sd="[^"]+"',
                                  f' d="M{LEAD_X0:g} {y:g} L{LEAD_X1:g} {y:g}"', tag), 1)
    return t


def numbers(t, leads):
    """Centre each number on its row, in his own column, and ring it."""
    rings = []
    for lead, n in zip(leads, "123"):
        t = re.sub(r'(matrix\(1\.193 0 0 1\.193 )[\d.]+( )[\d.]+(\)"[^>]*>\s*<tspan[^>]*>%s<)' % n,
                   lambda m: f"{m.group(1)}{NUM_X:g}{m.group(2)}{lead + TYPE * CAP / 2:.3f}{m.group(3)}",
                   t, count=1)
        rings.append((NUM_X, lead))
    return t, rings


def ring(cx, cy, label=None):
    c = (f'<circle cx="{cx:.3f}" cy="{cy:.3f}" r="{RING:g}" fill="none" '
         f'stroke="{N["500"]}" stroke-width="{TICK_W:g}"/>')
    if label is None:
        return c
    return c + (f'<text x="{cx:.3f}" y="{cy + DIGIT * CAP / 2:.3f}" '
                f'font-family="Arial" font-size="{DIGIT:g}" fill="{INK}" '
                f'text-anchor="middle">{label}</text>')


def codons():
    out = []
    for x, lab in CODONS:
        nm = "stop" if lab == "*" else lab
        star = lab == "*"
        y = CODON_Y + (STAR_DROP if star else 0.0)
        size = STAR_SIZE if star else IN_BLOCK
        out.append(f'<path data-name="{nm} site" d="M{x:g} {TICK0:g} L{x:g} {TICK1:g}" '
                   f'fill="none" stroke="{N["500"]}" stroke-width="{TICK_W:g}" '
                   f'stroke-linecap="round"/>')
        out.append(f'<text data-name="{nm} codon" x="{x:g}" y="{y:g}" '
                   f'font-family="Arial" font-size="{size:g}" fill="{INK}" '
                   f'text-anchor="middle">{lab}</text>')
    return "".join(out)


def shine(cx, cy, col):
    """Rays off a fluorescent protein, so it reads as emitting.

    Matteo has no glow in his set, and the slide used a clip-art lamp. Rays
    are an annotation, not a drawn molecule, and they say the protein itself
    is the light.
    """
    a, b = BARREL_W * BARREL_K / 2, BARREL_H * BARREL_K / 2
    out = [f'<ellipse cx="{cx:.3f}" cy="{cy:.3f}" rx="{a + RAY_GAP * 1.4:.3f}" '
           f'ry="{b + RAY_GAP * 1.4:.3f}" fill="{col}" opacity="{HALO:g}"/>']
    for i in range(RAY_N):
        th = 2 * math.pi * i / RAY_N
        c, s = math.cos(th), math.sin(th)
        r0 = a * b / math.hypot(b * c, a * s) + RAY_GAP
        r1 = r0 + RAY_LEN
        out.append(f'<path d="M{cx + r0 * c:.3f} {cy + r0 * s:.3f} '
                   f'L{cx + r1 * c:.3f} {cy + r1 * s:.3f}" fill="none" '
                   f'stroke="{col}" stroke-width="{TICK_W:g}" '
                   f'stroke-linecap="round"/>')
    return "".join(out)


def _rounded(x0, x1, yt, yb, left, right):
    """A block with Binder's corner radius, rounded on the ends asked for."""
    r, k = AT_R, AT_R * AT_K
    d = [f"M{x0 + (r if left else 0):g} {yt:g}", f"L{x1 - (r if right else 0):g} {yt:g}"]
    if right:
        d += [f"C{x1 - r + k:g} {yt:g} {x1:g} {yt + r - k:g} {x1:g} {yt + r:g}",
              f"L{x1:g} {yb - r:g}",
              f"C{x1:g} {yb - r + k:g} {x1 - r + k:g} {yb:g} {x1 - r:g} {yb:g}"]
    else:
        d += [f"L{x1:g} {yb:g}"]
    d += [f"L{x0 + (r if left else 0):g} {yb:g}"]
    if left:
        d += [f"C{x0 + r - k:g} {yb:g} {x0:g} {yb - r + k:g} {x0:g} {yb - r:g}",
              f"L{x0:g} {yt + r:g}",
              f"C{x0:g} {yt + r - k:g} {x0 + r - k:g} {yt:g} {x0 + r:g} {yt:g}"]
    else:
        d += [f"L{x0:g} {yt:g}"]
    return " ".join(d) + " Z"


def widen_at(t, lead):
    """Redraw the AT block wide enough to hold its caption.

    His was 21.892 across with a 8.756 intron, narrower than the "AT" sitting
    on it, and its corners had been squashed to 0.876 horizontally against
    1.460 everywhere else. Redrawn at 38.646 with circular corners, so the
    bar ends flush with the other two and the caption sits on the intron.
    """
    yt, yb = lead - 5.473, lead + 5.473
    b1 = BAR_X0 + AT_W * AT_SPLIT[0]
    b2 = BAR_X0 + AT_W * AT_SPLIT[1]
    x1 = BAR_X0 + AT_W
    parts = {"_5_prime_segment": _rounded(BAR_X0, b1, yt, yb, True, False),
             "Intron_segment": _rounded(b1, b2, yt, yb, False, False),
             "_3_prime_segment": _rounded(b2, x1, yt, yb, False, True),
             "Hairpin_encoding_outline": _rounded(BAR_X0, x1, yt, yb, True, True)}
    for eid, d in parts.items():
        m = re.search(r'<path[^>]*\sid="%s"[^>]*/>' % eid, t)
        tag = m.group(0)
        new = re.sub(r'\sd="[^"]+"', f' d="{d}"', tag)
        new = re.sub(r'\stransform="[^"]+"', "", new)      # it is drawn in place
        t = t.replace(tag, new, 1)
    return t, (b1 + b2) / 2


def main():
    t = SRC.read_text()

    # --- one stroke weight -------------------------------------------------
    t = re.sub(r"stroke-width:[\d.]+;?", "", t)
    t = re.sub(r"(\.cls-\d+\{[^}]*?stroke:#[0-9a-fA-F]{6})",
               lambda m: m.group(1) + f";stroke-width:{STROKE}", t)

    # Fig1d strokes its brackets and ticks lighter than its blocks and strand
    # (0.35 against 0.4); cls-16 is this file's two brackets.
    t = re.sub(r"(\.cls-16\{[^}]*?stroke-width:)[\d.]+",
               lambda m: m.group(1) + f"{TICK_W:g}", t)

    # --- palette -----------------------------------------------------------
    for bad, good in COLOURS.items():
        t = t.replace(bad, good)
    # the container is a container, not a drawn part: hairline, not outline
    t = re.sub(r"(\.cls-11\{[^}]*?)stroke:#4e5356",
               lambda m: m.group(1) + f'stroke:{N["300"]}', t)

    # --- type: one family, one size ---------------------------------------
    # NB the quoted family is written &apos;...&apos;, and that entity ENDS IN
    # A SEMICOLON. A naive [^;}]+ stops inside it and leaves an unclosed quote,
    # which makes every CSS rule after that point unparseable - the component
    # box renders solid black while the art above it looks fine.
    t = re.sub(r"font-family:\s*(?:&apos;.*?&apos;|'[^']*'|[^;}]+)",
               "font-family:Arial", t)
    t = re.sub(r"font-size:[\d.]+", f"font-size:{TYPE / MTX:.2f}", t)
    t = t.replace("matrix(1.082 0 0 1.193", f"matrix({MTX} 0 0 {MTX}")
    t = t.replace("mScarlett", "mScarlet")
    # construct_key: mscarlet_off -> "mScarlet*". The unspliced reporter is out
    # of frame, and the same asterisk marks the stop that puts it there.
    t = re.sub(r'(matrix\(1\.193 0 0 1\.193 126\.068 [\d.]+\)"[^>]*>\s*<tspan[^>]*>'
               r'mScarlet)(</tspan>)', r"\g<1>*\g<2>", t, count=1)

    # --- the container, redrawn clean --------------------------------------
    x0, y0, x1, y1 = BOX
    box = (f'<rect id="Component_box" data-name="Component box" class="cls-11" x="{x0:g}" y="{y0:g}" width="{x1 - x0:g}" '
           f'height="{y1 - y0:g}" rx="{RADIUS:g}"/>')
    t = re.sub(r'<(path|rect)[^>]*class="cls-11"[^>]*/>', box, t, count=1)

    # --- reporters, and captions in the box --------------------------------
    extra = []
    for cx, cy, col in BARRELS:
        f = cartoon.fragment("Selected_GFP_v1_1x.svg", "GFP v1 | 1x reusable",
                             scale=BARREL_K, at=(cx, cy), anchor="center")
        nm = "mNeonGreen" if col == S["reporter"] else "mScarlet"
        # his own rim is #6C6970; our map darkens it to neutral 700, which is
        # what made the cylinder read as a can
        f = f.replace(S["reporter"], col).replace(N["700"], N["600"])
        f = f.replace('<path data-name="Curved upper cap"',
                      f'<path data-name="Lit cap" d="{CAP_LIT}" '
                      f'fill="{palette.SURFACE}" opacity="0.3"/>'
                      '<path data-name="Curved upper cap"', 1)
        extra.append(f'<g id="{nm}_barrel" data-name="{nm} barrel">'
                     + shine(cx, cy, col) + f + "</g>")
    captions = []       # placed after the blocks have been moved
    # his placeholder squares are replaced, so drop them
    for m in list(re.finditer(r'<(path|rect)[^>]*class="cls-(20|21)"[^>]*/>', t)):
        tag = m.group(0)
        d = re.search(r'\sd="([^"]+)"', tag)
        if not d:
            continue
        n = [float(v) for v in re.findall(r"-?[\d.]+", d.group(1))]
        if max(n[0::2]) - min(n[0::2]) < 14:        # a square, not a bar
            t = t.replace(tag, "")
    # --- names, alignment, callouts, codons --------------------------------
    t = rename(t)
    t, leads, moved = snap_rows(t)
    t = relead(t, leads)
    t, _ = widen_at(t, leads[2])
    for eid, lab, col in CAPTIONS:
        m = re.search(r'<(?:path|rect)[^>]*\sid="%s"[^>]*/>' % eid, t)
        b = _bbox(m.group(0))
        tr = re.search(r'translate\(([-\d.]+) ([-\d.]+)\)', m.group(0))
        dx, dy = (float(tr.group(1)), float(tr.group(2))) if tr else (0.0, 0.0)
        captions.append(
            f'<text data-name="Caption | {lab}" x="{(b[0] + b[2]) / 2 + dx:.3f}" '
            f'y="{(b[1] + b[3]) / 2 + dy + IN_BLOCK * CAP / 2:.3f}" '
            f'font-family="Arial" font-size="{IN_BLOCK:g}" '
            f'fill="{col}" text-anchor="middle">{lab}</text>')
    extra += captions
    t, rings = numbers(t, leads)
    extra += [ring(cx, cy) for cx, cy in rings]
    extra += [ring(cx, cy, lab) for cx, cy, lab in CALLOUTS]
    extra.append(codons())
    t = t.replace('viewBox="0 0 408.66 197.501"', f'viewBox="0 0 408.66 {PAGE_H:g}"')

    t = t.replace("</svg>", "".join(extra) + "</svg>")

    DST.write_text(t, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-w", "2000", "-b", "white", str(DST),
                    "-o", str(DST.with_suffix(".png"))], check=False)
    subprocess.run(["inkscape", str(DST), "--export-type=pdf",
                    f"--export-filename={DST.with_suffix('.pdf')}"],
                   check=False, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    allowed = {c.lower() for c in (*palette.THEME.values(), *N.values(),
                                   *S.values(), *palette.DIVERGING, INK,
                                   palette.GRID, palette.SURFACE)}
    stray = sorted({c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", t)}
                   - allowed)
    css = sorted({float(w) for w in re.findall(r"stroke-width:([\d.]+)", t)})
    assert css == sorted({STROKE, TICK_W}), f"unexpected weights: {css}"
    print(f"wrote {DST}\n  strokes {STROKE} panel / {TICK_W} tick / "
          f"{0.48 * BARREL_K:.2f}+{0.3 * BARREL_K:.2f} his barrel"
          f"   type {TYPE}/{IN_BLOCK}"
          + (f"  LEFTOVER {stray}" if stray else "   palette ok"))


if __name__ == "__main__":
    main()
