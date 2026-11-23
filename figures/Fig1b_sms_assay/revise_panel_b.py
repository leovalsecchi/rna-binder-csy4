#!/usr/bin/env python3
"""Fig 1 panel b: Leo's prio-2 revisions, applied to the live panel SVG.

`output/Fig1b_sms_assay/Fluorescence V1.svg` is the file `Fig1_assembly` places
as panel b. It is a hand-placed descendant of `fix_fluorescence.py`'s output and
no script regenerated it, so this one does - the same idiom as
`fix_fluorescence.py`: Leo's layout is the input, the script only corrects it.

    python code/Fig1b_sms_assay/revise_panel_b.py

in   working/Fig1b_sms_assay/Fluorescence V1 source.svg   (pristine; made on
     first run from the live file, then never written again)
out  output/Fig1b_sms_assay/Fluorescence V1.svg           (what the page uses)

Reading from a pristine copy and writing to the live path makes the script
IDEMPOTENT and leaves the assembly untouched: re-running it re-derives the same
output from the same input rather than re-applying edits to its own output.

## The five changes (Leo, 2026-11-17)

1. **The three component rows become circular plasmids**, drawn as CLOSED
   rounded boxes with every corner rounded (Leo, 2026-11-17: "round edge boxes
   with no edges"). The first attempt drew only three sides, rising from the
   strand and closing across the top, which left two square corners where the
   verticals met the strand.
   The box hangs BELOW its genes (Leo, 2026-11-17: "have the genes at the top of
   the loop so it doesn't clash with the bracket"). It used to rise above them,
   which put each row's bracket and label inside the box. Now the blocks straddle
   the box's TOP edge and the bracket sits clear above it, outside the plasmid
   entirely.
   LOOP_RISE is 9.0, set by the tightest row: below each strand there are 12.40
   units to the next row's label on rows 1-2, and 11.47 to the component box
   floor on row 3, so 9.0 clears by 2.5 at worst. LOOP_R 4.0 leaves 1.0 of
   straight side, so the box reads as a rounded capsule rather than a rectangle.
   They go in **immediately after `Component_box`**, not at the end of the file:
   a closed box has a bottom edge running the length of the row, and appended
   last it would stroke straight through the coloured blocks. Placed there it
   sits above the box background and under everything else, so the bottom edge
   is hidden by the blocks and by the strand that already runs along it.
   Three separate boxes, one per row: these are three separately transfected
   plasmids, not one cassette (Leo's call, 2026-11-17).
   Each row's original STRAND is removed with the box that replaces it. The
   strand ran the full 243.33-377.17 at the row's y, while a rounded box's
   bottom edge only spans x0+r to x1-r, so the strand poked out past both
   rounded corners as a pair of stubs - the "edges" of Leo's note. The box is
   the backbone now. Only the three component rows lose their strand; the two
   reporter strands at y 185.232 are untouched.
2. **The straight arrow above the reporter hairpin is removed** - the vertical
   shaft at x 42.667 and its head. The two CURVED arrows either side and the
   diagonal to the excised intron stay: only the middle, straight one goes.
3. **The stop codons read STOP, not `*`.** `construct_key.svg` names that tick
   `*`, and the panel drew it at 12 pt (an asterisk is a small glyph in its em).
   The word is set in the caption class instead, on AUG's own baseline.
   `mScarlet*` KEEPS its asterisk: there it marks the out-of-frame product, a
   different statement from "a stop codon is here".
4. **"RNA scaffold" -> "Target array"** (the canonical name since 2026-11-23;
   it is narrower than the abbreviation it replaces, so the fit note below is
   historical), at the normal 8 pt, centred on its
   bracket. The full "RNA presenting module" does not fit and CANNOT be made to:
   at 8 pt it is 83.5 wide against 44.9 for the name it replaces, and with
   "mNeonGreen" (46.9) beside it the row needs 130.5 in a 124.2 slot. Setting it
   at 6.6 pt was tried and is worthless here - `common/retext.py` rewrites every
   font-size in a panel to the house size at assembly, so the shrink survives the
   standalone PNG and is undone on the page, which is where it matters. Leo's
   abbreviation is the fix: "RNA pres. module" is 64.0, so the pair comes to
   110.9 and clears by 13.3 at the size the page will actually print.
5. **"Binder" -> "Csy4".**
6. **The component box is 6.0 deeper**, and the mScarlet barrel below it moves
   down with it (Leo, 2026-11-17). Hanging the plasmids below their genes put
   row 3's box foot at y 110.79 against a box floor of 113.27 - a 2.5 margin
   where the top of the box has 5.4. The floor goes to 119.27, which gives row 3
   an 8.5 margin, and the barrel keeps its old clearance because it moves by the
   same 6.0.
7. **The barrel's arrow moves with it.** It is NOT inside the `mScarlet_barrel`
   group - it sits much earlier in the file - so shifting the group alone would
   leave the arrow pointing at where the barrel used to be. Shaft and head are
   both shifted by the same 6.0, which keeps the head's angle on the shaft.

Nothing else moves. Every edit is keyed to the exact coordinate or string it
targets and asserts it hit, so a changed input fails loudly instead of writing a
silently half-edited panel.
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import paths                                         # noqa: E402

LIVE = paths.fig_output("Fig1b_sms_assay") / "Fluorescence V1.svg"
SRC_DIR = paths.ROOT / "working" / "Fig1b_sms_assay"
SRC = SRC_DIR / "Fluorescence V1 source.svg"

STRAND_X0, STRAND_X1 = 243.327, 377.169      # the strand stubs either side
ROW_Y = (33.444, 67.618, 101.792)            # the three component strands
LOOP_RISE, LOOP_R = 9.0, 4.0

# the straight arrow above the reporter hairpin: shaft, then head
ARROW_SHAFT = "M42.667 134.368 L42.667 145.108"
ARROW_HEAD = "M42.667 148.583 L44.12 145.108 L41.214 145.108 Z"

AUG_BASELINE = 201.044                       # the codon row the STOPs join

# 6 + 7. the component box floor, and everything below it that must follow
BOX_DY = 6.0
BOX_BOTTOM_YS = ("109.665", "111.653", "113.265")   # the box path's lower half
ARROW_TO_BARREL = (
    ("M330.162 164.984 L318.637 149.237", ((330.162, 164.984), (318.637, 149.237))),
    ("M317.005 146.871 L320.269 148.751 L317.958 150.514 Z",
     ((317.005, 146.871), (320.269, 148.751), (317.958, 150.514))),
)


def bootstrap() -> str:
    """First run freezes the live file as the pristine input."""
    if not SRC.exists():
        SRC_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(LIVE, SRC)
        print(f"froze pristine source -> {SRC}")
    return SRC.read_text()


def drop_path(svg: str, d: str, what: str) -> str:
    """Remove the one <path> whose d is exactly `d`."""
    pat = re.compile(r'\s*<path[^>]*?d="' + re.escape(d) + r'"[^>]*/?>')
    svg, n = pat.subn("", svg)
    assert n == 1, f"{what}: expected 1 path, removed {n}"
    return svg


def retext(svg: str, old: str, new: str) -> str:
    """Rewrite a <text>'s content, matching the whole element so a substring of
    another label cannot be hit ('Binder' also starts no other string here, but
    'SMS' is a substring of the receptor id, hence matching inside <text>)."""
    pat = re.compile(r'(<text\b[^>]*>(?:(?!</text>).)*?)' + re.escape(old) +
                     r'((?:(?!</text>).)*?</text>)', re.S)
    svg, n = pat.subn(lambda m: m.group(1) + new + m.group(2), svg, count=1)
    assert n == 1, f"retext {old!r} -> {new!r}: matched {n}"
    return svg


def shrink(svg: str, label: str, cls: str) -> str:
    """Re-class the one <text> containing `label` (see change 4 for the sums)."""
    pat = re.compile(r'(<text\b[^>]*?)class="cls-18"((?:(?!</text>).)*?' +
                     re.escape(label) + r'(?:(?!</text>).)*?</text>)', re.S)
    svg, n = pat.subn(lambda m: m.group(1) + f'class="{cls}"' + m.group(2), svg, count=1)
    assert n == 1, f"shrink {label!r}: matched {n}"
    return svg


def move_x(svg: str, label: str, x: float) -> str:
    """Set the x of the one <text> containing `label` (it is middle-anchored)."""
    pat = re.compile(r'(<text transform="matrix\(1 0 0 1 )[\d.]+( [\d.]+\)"'
                     r'(?:(?!</text>).)*?' + re.escape(label) + r'(?:(?!</text>).)*?</text>)',
                     re.S)
    svg, n = pat.subn(lambda m: m.group(1) + f"{x:.3f}" + m.group(2), svg, count=1)
    assert n == 1, f"move_x {label!r}: matched {n}"
    return svg


def deepen_box(svg: str) -> str:
    """Drop the component box floor by BOX_DY. The box is a plain rounded rect,
    so only the three y values in its lower half move."""
    m = re.search(r'(<path id="Component_box"[^>]*?d=")([^"]+)(")', svg)
    assert m, "Component_box not found"
    d = m.group(2)
    for y in BOX_BOTTOM_YS:
        assert y in d, f"box floor value {y} not in the path"
        d = d.replace(y, f"{float(y) + BOX_DY:.3f}")
    return svg[: m.start(2)] + d + svg[m.end(2):]


def move_barrel(svg: str) -> str:
    """Shift the mScarlet barrel group down by BOX_DY."""
    tag = '<g id="mScarlet_barrel"'
    i = svg.find(tag)
    assert i >= 0, "mScarlet_barrel group not found"
    end = svg.index(">", i)
    assert "transform" not in svg[i:end], "barrel group already has a transform"
    return svg[:end] + f' transform="translate(0 {BOX_DY})"' + svg[end:]


def move_arrow(svg: str) -> str:
    """Shift the arrow that points at the barrel, which lives outside its group."""
    for d, pts in ARROW_TO_BARREL:
        assert d in svg, f"arrow path not found: {d[:30]}"
        head = "M" if d.startswith("M") else ""
        parts = [f"{'M' if k == 0 else 'L'}{x} {y + BOX_DY}" for k, (x, y) in enumerate(pts)]
        new = " ".join(parts) + (" Z" if d.rstrip().endswith("Z") else "")
        svg = svg.replace(d, new, 1)
    return svg


def plasmid(y: float) -> str:
    """A closed rounded box hanging BELOW the genes at `y`, corners all rounded.

    `y` is the row's old strand line, which is where the coloured blocks sit, so
    it becomes the box's TOP edge and the blocks straddle it - genes on the
    circle. The top edge is therefore hidden behind the blocks and shows only in
    the stubs either side. Same stroke class as the strand it replaces."""
    bot, r = y + LOOP_RISE, LOOP_R
    x0, x1 = STRAND_X0, STRAND_X1
    return (f'<path class="cls-13" d="M{x0 + r} {y} '
            f'L{x1 - r} {y} '
            f'Q{x1} {y} {x1} {y + r} '
            f'L{x1} {bot - r} '
            f'Q{x1} {bot} {x1 - r} {bot} '
            f'L{x0 + r} {bot} '
            f'Q{x0} {bot} {x0} {bot - r} '
            f'L{x0} {y + r} '
            f'Q{x0} {y} {x0 + r} {y} '
            f'Z"/>')


def stop_codons(svg: str) -> str:
    """The two 12 pt `*` ticks become 8 pt STOP on AUG's baseline.

    Both are already text-anchor:middle (cls-38), so the x is the tick's centre
    and carries over unchanged; only the class and the baseline move."""
    pat = re.compile(
        r'<text transform="matrix\(1 0 0 1 ([\d.]+) [\d.]+\)"([^>]*?)class="cls-38"'
        r'((?:(?!</text>).)*?)</text>', re.S)
    found = pat.findall(svg)
    assert len(found) == 2, f"expected 2 asterisk texts, found {len(found)}"
    for x, _, body in found:
        assert "*" in body, f"text at x={x} is not the asterisk"

    def sub(m):
        x = m.group(1)
        return (f'<text transform="matrix(1 0 0 1 {x} {AUG_BASELINE})" x="0" y="0" '
                f'xml:space="preserve" class="cls-18"><tspan x="0" y="0">STOP</tspan></text>')
    svg, n = pat.subn(sub, svg)
    assert n == 2
    return svg


def main() -> None:
    svg = bootstrap()

    # 2. the straight arrow over the reporter hairpin
    svg = drop_path(svg, ARROW_SHAFT, "straight arrow shaft")
    svg = drop_path(svg, ARROW_HEAD, "straight arrow head")

    # 4 + 5. the two renamed blocks
    svg = retext(svg, "RNA scaffold", "Target array")
    svg = retext(svg, "Binder", "Csy4")

    # 3. stop codons
    svg = stop_codons(svg)

    # 6 + 7. room for the plasmid boxes, and the barrel follows the box floor
    svg = deepen_box(svg)
    svg = move_barrel(svg)
    svg = move_arrow(svg)

    # 1. the plasmid boxes replace each row's strand - see change 1
    for y in ROW_Y:
        svg = drop_path(svg, f"M{STRAND_X0} {y} L{STRAND_X1} {y}", f"row strand y={y}")
    #    inserted right after Component_box so they sit above the box background
    #    and BELOW the blocks
    loops = "".join(plasmid(y) for y in ROW_Y)
    m = re.search(r'<path id="Component_box"[^>]*/>', svg)
    assert m, "Component_box element not found"
    svg = svg[: m.end()] + loops + svg[m.end():]

    LIVE.write_text(svg)
    print(f"wrote {LIVE}")
    print(f"  plasmid loops at y={ROW_Y}, rise {LOOP_RISE}, radius {LOOP_R}")


if __name__ == "__main__":
    main()
