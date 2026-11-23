"""Construct maps: one way of drawing a module bar, shared by every panel.

Panels b and d both draw the same cassette and drew it differently - different
heights, radii and colours. This is the one builder, so a module is the same
shape and the same colour wherever it appears.

    from common import constructs as K
    svg = K.bar(20, 40, ["attB_GT", "at", ("binder", "sms"),
                         "scaffold", "attB_GA"])

Shape follows SBOL Visual
-------------------------
SBOL Visual is the diagram standard for engineering biology, and it carries part
type in the SHAPE of a glyph: its glyphs "purposefully avoided specifying
stylistic features like line width or color, focusing instead on distinctive
shapes", and the standard "does not impose constraints on how diagrams are
generated". So colour is ours to choose, but a plain rectangle for a coding
sequence is the "generic glyph" its compliance review flags as not best
practice. Accordingly:

* **Coding sequences** - binder, SMS, RNase, both reporters - are the CDS
  pentagon: flat back, pointed front, so reading direction is visible.
  Proportions are from `cds.svg` (point length 10.88 against a 22.80 body).
* **attB sites** are the Specific Recombination Site glyph: a right-pointing
  triangle, from `specific-recombination-site.svg`.
* **Untranslated regions** - 5' UTR, AT, the RNA scaffold - have no specific
  glyph in the catalogue, so they stay plain rounded boxes.
* **stop** is SBOL's Translation End glyph - a circled asterisk on a wavy stem -
  drawn from the literal geometry of `translation-end.svg`, scaled so the stem
  foot meets the strand. AUG stays a plain tick: the catalogue has no start-codon
  glyph, and nor does it have one for a hairpin or a 5' UTR, which is why AT and
  the UTR are plain boxes.

Colour is a role, and the cartoons use the same ones
----------------------------------------------------
* **Binder** is the designed protein, so it carries MAIN - the same petrol as
  the binder in the imported cartoons.
* **SMS** carries one label. Its RNase lobe is a short segment at the front of
  the same glyph, in the sand that Matteo's SMS cartoon gives that lobe, so a map
  and the drawing beside it agree without a legend. The whole thing is called
  SMS; the lobe is not labelled separately.
* **One polypeptide is one glyph.** The binder is fused N-terminal to MgIRE1, so
  `("binder", "sms")` draws a single CDS outline with a hairline between the
  two - one coding sequence read in one frame, not two genes in a row.
* **AT** takes the red of the piece the drawn hairpin loses - the warm end of
  the diverging ramp, there being no free red in the synbio block. The cassette
  is there to be cut, so it carries the colour of the cut-out piece. The
  **RNA scaffold** and **F30** are the same object as each other and take a
  lighter step, being bulk context rather than the active element.
* **mNeonGreen** and **mScarlet** are emission-matched, deliberately brighter
  than the rest of the palette, so the maps echo the microscopy panels. A
  reporter that is not expressed goes neutral: on/off is lightness, not hue.
"""
import yaml

from . import palette
from .paths import ROOT

S = palette.SYNBIO
_cfg = yaml.safe_load((ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'

LINE = S["outline"]
INK = palette.INK
MUT = palette.NEUTRAL["600"]
BG = palette.SURFACE

H = 26            # module height
R = 3.0           # corner radius, boxes only
SW = 1.5          # module outline
GAP = 3           # between adjacent modules
PAD = 13          # label padding inside a module
FS = 11.5
POINT = 0.48 * H  # CDS point length, from cds.svg (10.88 / 22.80)
SITE_W = 0.9 * H  # recombination triangle, near-square as in the glyph

# promoter.svg, verbatim: an L from the transcription start, with an arrowhead.
# Its foot is at (7.5, 39.75) in a 45-unit box and it reaches up to y=5.25.
_PROM = ('<path d="m 7.5001114,39.746408 0,-27 28.9999996,0" fill="none" '
         'stroke="{c}" stroke-width="3" stroke-linecap="round" '
         'stroke-linejoin="round"/>'
         '<path d="m 29.000111,5.2464081 8.5,7.4999999 -8.5,7.3333" '
         'fill="none" stroke="{c}" stroke-width="3" stroke-linecap="round" '
         'stroke-linejoin="round"/>')
_PROM_K = H / 34.5        # its height above the strand becomes one module tall
_PROM_FOOT = (7.5, 39.75)
PROM_W = 30.0 * _PROM_K
SPACER_W = 15.0   # bare strand, where two neighbours would otherwise crowd
TAIL_W = 13       # a second domain at the front of one glyph

# name -> (label, fill, text colour, role[, tail fill])
# A tail is a short segment at the front of the same glyph: a second domain of
# one molecule, inside one outline and under one label.
#   cds   coding sequence: pointed front, direction visible
#   box   a region with no specific glyph
#   site  recombination site: a triangle
#   tick  a codon marker: a rule through the strand
MODULES = {
    "utr":          ("5′ UTR",  palette.NEUTRAL["200"], INK, "box"),
    "at":           ("AT",           palette.SYNBIO["intron"], BG,  "box"),
    "scaffold":     ("RNA scaffold", palette.THEME["200"],   INK, "box"),
    "f30":          ("F30",          palette.THEME["200"],   INK, "box"),
    "binder":       ("Binder",       S["binder"],         BG,  "cds"),
    "sms":          ("SMS",          S["enzyme"],         INK, "cds",
                     S["enzyme_alt"]),
    "mneongreen":   ("mNeonGreen",   S["reporter"],       INK, "cds"),
    "mscarlet":     ("mScarlet",     S["reporter_alt"],   BG,  "cds"),
    "mscarlet_off": ("mScarlet",     palette.NEUTRAL["300"], MUT, "cds"),
    "attB_GT":      ("attB-GT",      palette.NEUTRAL["500"], BG, "site"),
    "attB_GA":      ("attB-GA",      palette.NEUTRAL["500"], BG, "site"),
    "attP_GT":      ("attP-GT",      palette.NEUTRAL["300"], INK, "site"),
    "attP_GA":      ("attP-GA",      palette.NEUTRAL["300"], INK, "site"),
    "attL":         ("attL",         palette.NEUTRAL["700"], BG, "site"),
    "attR":         ("attR",         palette.NEUTRAL["700"], BG, "site"),
    "promoter":     ("P",            None,                MUT, "promoter"),
    "spacer":       ("",             None,                MUT, "spacer"),
    "pa":           ("pA",           palette.NEUTRAL["200"], INK, "box"),
    "aug":          ("AUG",          None,                MUT, "tick"),
    "stop":         ("stop",         None,                MUT, "stopsite"),
}


def _w(key):
    """Width one module occupies; a tuple is a fused group, a tick takes none."""
    if isinstance(key, (tuple, list)):
        return sum(_body(k) for k in key) + POINT
    label, _, _, role = MODULES[key][:4]
    if role == "tick":
        return 0.0
    if role == "site":
        return SITE_W
    if role == "promoter":
        return PROM_W
    if role == "spacer":
        return SPACER_W
    w = len(label) * FS * 0.56 + 2 * PAD + (TAIL_W if _tail(key) else 0)
    return w + POINT if role == "cds" else w


def _body(key):
    """Width of a module without its CDS point, so fusions share one point."""
    role = MODULES[key][3]
    return _w(key) - (POINT if role == "cds" else 0)


def _tail(key):
    """The fill of this module's front segment, or None."""
    m = MODULES[key]
    return m[4] if len(m) > 4 else None


def _is_tick(key):
    return (not isinstance(key, (tuple, list))
            and MODULES[key][3] in ("tick", "stopsite"))


# translation-end.svg, verbatim. Its own box is 45x45; the stem foot is at
# (22.5, 39.13) and the circled asterisk sits above it.
_STOP_GLYPH = (
    '<circle cx="22.499998" cy="13.489317" r="7" fill="none" stroke="{c}" '
    'stroke-width="3"/>'
    '<path d="m 22.433046,30.321771 c 4.327115,-4.802042 -4.168807,-4.698615 '
    '0,-8.867421" fill="none" stroke="{c}" stroke-width="3" '
    'stroke-linecap="round" stroke-linejoin="round"/>'
    '<path d="m 22.485815,39.13431 c 4.327115,-4.802042 -4.168806,-4.698614 '
    '0,-8.86742" fill="none" stroke="{c}" stroke-width="3" '
    'stroke-linecap="round" stroke-linejoin="round"/>'
    '<text x="22.592772" y="24.661015" font-family="Arial, Helvetica, sans-serif" '
    'font-size="20" font-weight="bold" fill="{c}" text-anchor="middle">*</text>')
_STOP_K = 0.80            # 45-unit glyph down to about one module height
_STOP_FOOT = (22.5, 39.13)


def stop_site(x, y, col):
    """SBOL Translation End, placed with its stem foot on the strand."""
    dx = x - _STOP_FOOT[0] * _STOP_K
    dy = y - _STOP_FOOT[1] * _STOP_K
    return (f'<g transform="translate({dx:.2f},{dy:.2f}) scale({_STOP_K})">'
            + _STOP_GLYPH.format(c=col) + '</g>')


def width(items):
    """Total width a bar occupies, so a caller can centre it."""
    return max(0.0, sum(_w(k) + GAP for k in items if not _is_tick(k)) - GAP)


def shape(role, x, y, w):
    """Outline path for one module, in SBOL Visual's shape for its part type."""
    t, b = y - H / 2, y + H / 2
    if role == "cds":
        return (f"M{x:.1f},{t:.1f} L{x+w-POINT:.1f},{t:.1f} L{x+w:.1f},{y:.1f} "
                f"L{x+w-POINT:.1f},{b:.1f} L{x:.1f},{b:.1f} Z")
    if role == "site":
        return f"M{x:.1f},{t:.1f} L{x+w:.1f},{y:.1f} L{x:.1f},{b:.1f} Z"
    return (f"M{x+R:.1f},{t:.1f} H{x+w-R:.1f} A{R},{R} 0 0 1 {x+w:.1f},{t+R:.1f} "
            f"V{b-R:.1f} A{R},{R} 0 0 1 {x+w-R:.1f},{b:.1f} H{x+R:.1f} "
            f"A{R},{R} 0 0 1 {x:.1f},{b-R:.1f} V{t+R:.1f} "
            f"A{R},{R} 0 0 1 {x+R:.1f},{t:.1f} Z")


def promoter(x, y, col):
    """SBOL's promoter, placed with its transcription start on the strand."""
    dx = x - _PROM_FOOT[0] * _PROM_K
    dy = y - _PROM_FOOT[1] * _PROM_K
    return (f'<g transform="translate({dx:.2f},{dy:.2f}) scale({_PROM_K:.4f})">'
            + _PROM.format(c=col) + '</g>')


def _text(x, y, s, size, col, weight=500, anchor="middle"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
            f'font-size="{size}" fill="{col}" font-weight="{weight}" '
            f'text-anchor="{anchor}">{s}</text>')


def _outline(d, fill="none"):
    return (f'<path d="{d}" fill="{fill}" stroke="{LINE}" stroke-width="{SW}" '
            f'stroke-linejoin="round"/>')


def _fused(x, y, keys, uid):
    """One coding sequence made of several named parts, inside one outline."""
    total = sum(_body(k) for k in keys) + POINT
    outer = shape("cds", x, y, total)
    clip = f"fuse{uid}"
    o = [f'<clipPath id="{clip}"><path d="{outer}"/></clipPath>',
         f'<g clip-path="url(#{clip})">']
    cursor = x
    for k in keys:                      # square-cut fills; the clip shapes them
        b = _body(k)
        o.append(f'<rect x="{cursor:.1f}" y="{y-H/2:.1f}" width="{b+POINT:.1f}" '
                 f'height="{H}" fill="{MODULES[k][1]}"/>')
        if _tail(k):                    # that part's own second domain
            o.append(f'<rect x="{cursor+b-TAIL_W:.1f}" y="{y-H/2:.1f}" '
                     f'width="{TAIL_W+POINT:.1f}" height="{H}" '
                     f'fill="{_tail(k)}"/>')
        cursor += b
    o.append("</g>")
    cursor = x
    for k in keys:                      # hairlines: between parts, and at a tail
        b = _body(k)
        if _tail(k):
            o.append(f'<path d="M{cursor+b-TAIL_W:.1f},{y-H/2:.1f} '
                     f'L{cursor+b-TAIL_W:.1f},{y+H/2:.1f}" stroke="{LINE}" '
                     f'stroke-width="1.0"/>')
        cursor += b
        if k is not keys[-1]:
            o.append(f'<path d="M{cursor:.1f},{y-H/2:.1f} L{cursor:.1f},'
                     f'{y+H/2:.1f}" stroke="{LINE}" stroke-width="1.0"/>')
    o.append(_outline(outer))
    cursor = x
    for k in keys:
        b = _body(k)
        label, _, text, _ = MODULES[k][:4]
        o.append(_text(cursor + (b - (TAIL_W if _tail(k) else 0)) / 2,
                       y + FS * 0.35, label, FS, text))
        cursor += b
    return "".join(o)


def bar(x, y, items, strand=14):
    """A construct drawn left to right from `x`, centred on baseline `y`.

    """
    for a, b in zip(items, items[1:]):
        if _is_tick(a) and _is_tick(b):
            raise ValueError(
                f"{a!r} and {b!r} are both marks on the strand and take no "
                f"width, so they would be drawn at the same point; put a "
                f"module between them or place one of them in its own bar")
    total = width(items)
    o = [f'<path d="M{x-strand:.1f},{y} L{x+total+strand:.1f},{y}" '
         f'stroke="{LINE}" stroke-width="2.0" stroke-linecap="round" '
         f'fill="none"/>']
    cursor = x
    for i, key in enumerate(items):
        if isinstance(key, (tuple, list)):
            o.append(_fused(cursor, y, list(key), f"{int(x)}_{int(y)}_{i}"))
            cursor += _w(key) + GAP
            continue
        label, fill, text, role = MODULES[key][:4]
        if role == "spacer":
            cursor += _w(key) + GAP
            continue
        if role == "promoter":
            o.append(promoter(cursor, y, LINE))
            o.append(_text(cursor + PROM_W / 2, y + H / 2 + 13, label,
                           FS - 1.5, text, 600))
            cursor += _w(key) + GAP
            continue
        if role == "stopsite":
            o.append(stop_site(cursor, y, LINE))
            continue
        if role == "tick":
            o.append(f'<path d="M{cursor:.1f},{y-9} L{cursor:.1f},{y+9}" '
                     f'stroke="{LINE}" stroke-width="2.0" stroke-linecap="round"/>')
            o.append(_text(cursor, y - 14, label, FS - 1.5, text, 600))
            continue
        w = _w(key)
        outer = shape(role, cursor, y, w)
        o.append(f'<path d="{outer}" fill="{fill}" stroke="none"/>')
        tail = _tail(key)
        if tail:                     # the front segment, clipped to the glyph
            cut = cursor + w - POINT - TAIL_W
            uid = f"t{int(cursor)}_{int(y)}"
            o.append(f'<clipPath id="{uid}"><path d="{outer}"/></clipPath>'
                     f'<rect x="{cut:.1f}" y="{y-H/2:.1f}" width="{w:.1f}" '
                     f'height="{H}" fill="{tail}" clip-path="url(#{uid})"/>'
                     f'<path d="M{cut:.1f},{y-H/2:.1f} L{cut:.1f},{y+H/2:.1f}" '
                     f'stroke="{LINE}" stroke-width="1.0"/>')
        o.append(_outline(outer))
        if role == "site":
            o.append(_text(cursor + w / 2, y - H / 2 - 6, label, FS - 1.5,
                           MUT, 600))
        else:
            body = (w - POINT if role == "cds" else w) - (TAIL_W if _tail(key) else 0)
            o.append(_text(cursor + body / 2, y + FS * 0.35, label, FS, text))
        cursor += w + GAP
    return "".join(o)


# The five bars in figure 1, so panels b and d cannot drift apart.
BARS = {
    "donor":        ["attB_GT", "at", "binder", "sms", "scaffold", "attB_GA"],
    "landing_pad":  ["promoter", "attP_GT", "attP_GA", "pa"],
    "integrated_genome": ["promoter", "attL", "at", "binder", "sms", "f30",
                          "attR", "pa"],
    "reporter_rna": ["utr", "mneongreen", "scaffold"],
    "sensor":       [("binder", "sms")],
    "readout":      ["at", "mscarlet", "stop"],
    "donor":        ["attB_GT", "at", ("binder", "sms"), "scaffold",
                     "attB_GA"],
    "integrated":   ["utr", "at", ("binder", "sms"), "f30"],
}


def key(width_px=980, x=190, step=42):
    """Every module once, then the five bars: the sheet for agreeing colours."""
    o, y = [], 46
    for name in MODULES:
        o.append(_text(20, y + 4, name, 11, MUT, 500, "start"))
        o.append(bar(x, y, [name]))
        y += step
    o.append(_text(20, y + 4, "(binder, sms)", 11, MUT, 500, "start"))
    o.append(bar(x, y, [("binder", "sms")]))
    y += step + 18
    for name, items in BARS.items():
        o.append(_text(20, y + 4, name, 11, MUT, 500, "start"))
        o.append(bar(x, y, items))
        y += 46
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_px}" '
            f'height="{y}" viewBox="0 0 {width_px} {y}">'
            f'<rect width="{width_px}" height="{y}" fill="{BG}"/>'
            + "".join(o) + '</svg>')


def _selfcheck():
    """Every module draws, in its SBOL shape, with palette colours only."""
    import re
    allowed = {c.lower() for c in (*palette.THEME.values(),
                                   *palette.NEUTRAL.values(),
                                   *palette.SYNBIO.values(),
                                   # the exon red AT carries, COLORS.md rule 4
                                   *palette.DIVERGING,
                                   palette.INK, palette.SURFACE)}
    for name, spec in MODULES.items():
        role = spec[3]
        svg = bar(0, 40, [name])
        assert "<path" in svg, name
        for h in re.findall(r"#[0-9A-Fa-f]{6}", svg):
            assert h.lower() in allowed, f"{name}: {h} is not a palette colour"
        if role == "cds":                       # five points, not four
            d = shape("cds", 0, 40, 100)
            assert d.count("L") == 4 and d.endswith("Z"), f"{name}: not a pentagon"
        if role == "site":
            assert shape("site", 0, 40, 20).count("L") == 2, "site is a triangle"
        if role == "stopsite":
            assert "<circle" in svg and ">*<" in svg, "stop site is a circled *"
        if role == "promoter":
            assert svg.count("<path") >= 2 and "scale(" in svg, "promoter is the bent arrow"
    assert width(["binder"]) > 0 and width(["aug"]) == 0, "tick takes no width"
    assert width(["stop"]) == 0, "the stop site sits on the strand"
    try:
        bar(0, 40, ["aug", "stop"])
    except ValueError:
        pass
    else:
        raise AssertionError("two adjacent marks must be refused, not stacked")
    for items in BARS.values():
        flat = [k for it in items for k in (it if isinstance(it, tuple) else [it])]
        assert set(flat) <= set(MODULES), "unknown module in a bar"
    fused = bar(0, 40, [("binder", "sms")])
    assert fused.count(f'stroke-width="{SW}"') == 1, "one outline for one protein"
    assert "Binder" in fused and "SMS" in fused and "RNase" not in fused
    sms = bar(0, 40, ["sms"])
    assert sms.count("SMS") == 1 and "RNase" not in sms, "SMS carries one label"
    assert sms.count(f'stroke-width="{SW}"') == 1, "one outline around the glyph"
    assert S["enzyme_alt"] in sms, "the RNase lobe is drawn"
    # the cartoons and the map must agree on the molecules they share
    assert MODULES["binder"][1] == S["binder"]
    assert (MODULES["sms"][1], _tail("sms")) == (S["enzyme"], S["enzyme_alt"])
    assert MODULES["at"][1] == palette.SYNBIO["intron"], "AT must match the excised piece"
    print(f"ok: {len(MODULES)} modules, {len(BARS)} bars, SBOL shapes")


if __name__ == "__main__":
    _selfcheck()
