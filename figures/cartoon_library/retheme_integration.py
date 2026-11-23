"""Re-theme Matteo's cassette-integration SVGs: our colours and labels, his geometry.

Nothing is redrawn. Each file keeps its paths, its structure and its layout; only
two things change:

* the fills and strokes declared in its `<style>` block, rewritten per CSS class
  according to which named block uses that class, and
* the text inside `<tspan>`s.

Class numbers differ between his files, so the mapping is keyed on the element
ids that use a class (`Kinase_fill`, `_5_prime_segment`, ...), not on the class
number or on the source hex - the same hex means different things in different
files, and `#c54d66` is both a hairpin segment and an att half-site.

His hybrid-site convention is preserved: attL and attR each take one half-site
colour from the donor and one from the landing pad, which is what shows that a
hybrid site is made of both. Donor sites are therefore blues and genomic sites
greys, so the halves stay readable after recombination.

    python retheme_integration.py <indir> <outdir>
"""
import re
import sys
from pathlib import Path

# --- the thesis palette -----------------------------------------------------
# Read from the repo when its share is up, so this tracks whatever COLORS.md
# says; the mirror below only covers the share being down, which it often is.
REPO = next((r for r in (Path(__file__).resolve().parents[2],
        Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_ThesisFigures"),
        Path("/Volumes/bsse/wicky/users/lvalsecchi/Leo_Valsecchi/LV_ThesisFigures"))
        if (r / "config/palette.json").exists()), Path("/nonexistent"))
try:
    sys.path.insert(0, str(REPO / "code"))
    from common import palette as _pal                       # noqa: E402
    THEME, NEUTRAL, SYNBIO = _pal.THEME, _pal.NEUTRAL, _pal.SYNBIO
    DIVERGING_WARM = _pal.DIVERGING[-3:]
    INK, SURFACE = _pal.INK, _pal.SURFACE
    print("palette: from the repo")
except Exception:
    print("palette: mirror (repo share unreachable)")

    THEME = {"200": "#bfe0f5", "300": "#95c9e8", "400": "#65aed8", "500": "#3893c2",
         "600": "#0678a6", "700": "#015d82", "800": "#01435f"}
    NEUTRAL = {"100": "#f0f6fa", "200": "#e3e9ec", "300": "#cfd5d9", "400": "#b3b8bc",
               "500": "#919699", "600": "#6d7276", "700": "#4e5356"}
    SYNBIO = {"rna": "#65aed8", "binder": "#015d82", "reporter": "#3fbf4f",
              "reporter_alt": "#e03127", "enzyme": "#3d8f85",
              "enzyme_alt": "#ffc233", "outline": "#4e5356", "intron": "#a8232a"}
    DIVERGING_WARM = ["#efac8d", "#d97c50", "#ae4e1a"]      # the AT hairpin's family
    INK, SURFACE = "#1d2428", "#ffffff"

# id substring -> our fill. Longest match wins, so "Coding_block-3" beats
# "Coding_block".
# Each module keeps Matteo's geometry; only its fill changes. The receptor is
# the exception: _fuse_receptor redraws it as Binder | linker | SMS.
FILL_BY_ID = {
    # The AT hairpin keeps HIS three segments. They are coloured like the drawn
    # hairpin in cartoon_library/hairpin.py - exons theme 600, excised intron
    # warm red - so the map and the cartoon show the same molecule, and the
    # piece that is cut out is the piece that is red in both.
    "_5_prime_segment": THEME["600"],
    "Intron_segment": SYNBIO["intron"],           # palette v0.5, 2026-11-13
    "_3_prime_segment": THEME["600"],
    # The receptor (his CC, GFP, 2A, VHH, kinase, RNase) is redrawn by
    # _fuse_receptor as Binder | linker | SMS, so no fill is mapped for it here.
    # the landing pad's marker and terminator
    "mNeonGreen_replacement_cassette": SYNBIO["reporter"],
    "pA_block": NEUTRAL["200"],                      # constructs.py "pa"
    # att sites are neutral (constructs.py): donor attB 500, genomic attP 300.
    # Blue is the binder's and must not also mean "donor site".
    "Complementary_left_half-site-2": NEUTRAL["500"],
    "Complementary_right_half-site-2": NEUTRAL["500"],
    "Complementary_left_half-site": NEUTRAL["500"],
    "Complementary_right_half-site": NEUTRAL["500"],
    "Promoter_box___white": SURFACE,
    "Promoter_arrowhead": INK,
}
# the landing pad's own sites are the genomic ones: attP, neutral 300
PAD_SITES = {"Complementary_left_half-site": NEUTRAL["300"],
             "Complementary_right_half-site": NEUTRAL["300"],
             "Complementary_left_half-site-2": NEUTRAL["300"],
             "Complementary_right_half-site-2": NEUTRAL["300"]}
# His hybrid logic, kept: in the integrated locus each site is one half from
# the pad and one from the donor. Only the two hexes change.
HYBRID_SITES = {"Complementary_left_half-site": NEUTRAL["300"],   # pad GT half
                "Complementary_right_half-site": NEUTRAL["500"],  # donor GT half
                "Complementary_left_half-site-2": NEUTRAL["500"],  # donor GA half
                "Complementary_right_half-site-2": NEUTRAL["300"]}  # pad GA half

STROKE_MAP = {"#000000": SYNBIO["outline"], "#8a8790": NEUTRAL["500"]}
TEXT_FILL = {"#000000": INK, "#8a8790": NEUTRAL["500"]}

# Same wording as Fig1d_integrated_construct. The identity entries are listed
# so that _relabel re-centres every bracket label on its own bracket.
LABELS = {
    "Recorder": "Recorder",
    # 2026-11-23: one name per thing with the main text - the four parts are the
    # binder, the sensor, the recorder and the target array. Re-centring on the
    # bracket (below) is what lets the wording change at all.
    "Target module": "Target array",
    "Sensor module": "Sensor",
}


def _classes(text):
    """class -> the set of element ids that carry it."""
    use = {}
    for pat in (r'id="([^"]+)"[^>]*class="(cls-\d+)"',
                r'class="(cls-\d+)"[^>]*id="([^"]+)"'):
        for m in re.finditer(pat, text):
            a, b = m.group(1), m.group(2)
            cls, eid = (b, a) if b.startswith("cls-") else (a, b)
            use.setdefault(cls, set()).add(eid)
    return use


def _fill_for(ids, overrides):
    """Our fill for a class, from the longest id key that matches one of its ids."""
    best = None
    for eid in ids:
        for key in sorted({**FILL_BY_ID, **overrides}, key=len, reverse=True):
            if key in eid:
                if best is None or len(key) > len(best[0]):
                    best = (key, {**FILL_BY_ID, **overrides}[key])
                break
    return best[1] if best else None


HAIRPIN_KEEP = 0.60        # of its original length


def _scale_x(d, x0, k):
    """Scale a path's x coordinates about x0, leaving every y alone."""
    toks = re.findall(r'[A-Za-z]|-?[\d.]+', d)
    out, i, cmd = [], 0, None
    while i < len(toks):
        t = toks[i]
        if t.isalpha():
            cmd, i = t, i + 1
            out.append(t)
            continue
        n = {"M": 2, "L": 2, "C": 6, "Z": 0}.get((cmd or "L").upper(), 2)
        vals = [float(v) for v in toks[i:i + n]]
        for j, v in enumerate(vals):
            out.append(f"{x0 + (v - x0) * k:g}" if j % 2 == 0 else f"{v:g}")
        i += n
    res, k2 = [], 0
    for t in out:
        res.append(t)
    return " ".join(res).replace("M ", "M").replace("L ", "L").replace("C ", "C")


def _shorten_hairpin(text):
    """Bring the AT hairpin down to a short block, and its bracket with it."""
    m = re.search(r'id="Hairpin_encoding_outline"[^>]*\sd="([^"]+)"', text)
    if not m:
        return text
    xs = [float(v) for v in re.findall(r'(?:[ML]\s*|C\s*)(-?[\d.]+)', m.group(1))]
    x0 = min(xs)
    for pid in ("_5_prime_segment", "Intron_segment", "_3_prime_segment",
                "Hairpin_encoding_outline", "Recorder_bracket"):
        for mm in list(re.finditer(rf'(id="{pid}"[^>]*\sd=")([^"]+)(")', text)):
            text = text.replace(mm.group(0),
                                mm.group(1) + _scale_x(mm.group(2), x0,
                                                       HAIRPIN_KEEP) + mm.group(3))
    return text


# The receptor is one polypeptide: Binder - linker - SMS (kinase + RNase lobe).
#
# His cassette has CC | GFP | 2A | VHH | kinase | RNase. The 2A is a SELF-
# CLEAVING peptide: it splits his chain into a separate target protein and a
# separate receptor. Ours has no target to release - the binder is fused
# straight to MgIRE1 - so the 2A goes and a flexible linker takes its place,
# drawn as the shorter grey bar that Fig1d_integrated_construct uses.
#
# Each part is given its own width, because they are not one object being
# scaled: the binder is a short designed protein (~80 aa) against MgIRE1's
# ~450, so it is drawn small, and the room that frees goes to the target.
# The kinase:RNase split inside the SMS is Matteo's.
# 2026-11-13: no linker (Leo) - the binder abuts the SMS. The 4 units it had
# go to the two SMS lobes in Matteo's own 25.42 : 16.08 split, so SMS_END and
# the receptor's span are unchanged.
_OLD_LINKER = 4.0
BINDER_W, LINKER_W = 18.0, 0.0
KINASE_W = 25.42 * (25.42 + 16.08 + _OLD_LINKER) / (25.42 + 16.08)
RNASE_W = 16.08 * (25.42 + 16.08 + _OLD_LINKER) / (25.42 + 16.08)
BINDER_END = BINDER_W
LINKER_END = BINDER_END + LINKER_W
RNASE_START = LINKER_END + KINASE_W
SMS_END = RNASE_START + RNASE_W         # the receptor's width on the strand
LINKER_INSET = 1.5                      # the linker is thinner than a block
ATG_REACH = 7.0                         # how far left the ATG elbow carries, as his
#
# Order on the strand, which is constructs.py's own MAPS["donor"]:
#   attB-GT | AT hairpin | (binder, sms) | RNA scaffold | attB-GA
# so the SENSOR (the protein) sits where his target module was, and the TARGET
# (our RNA scaffold) goes to its right. Shortening the hairpin frees exactly
# the room the scaffold needs, so the receptor slides left by that much and
# keeps its 78.6 span - the span Fig1d uses - and the att sites never move.
# The scaffold is drawn large on purpose: it is an F30 three-way junction
# carrying two Csy4 hairpins, and later panels show it at its real size, so a
# token box here would contradict them. Widening it narrows the receptor by
# the same amount - the strand between the hairpin and attB-GA is fixed.
# 2026-11-13: small gaps only (Leo: "no holes between the genes, small ones are
# okay"). GAP 5 -> 2; the 2 x 3 units saved go to the scaffold, so it still
# ends where it did and the strand up to attB-GA is unchanged.
_OLD_GAP = 5.0
GAP = 2.0                               # between elements on the strand
SCAFFOLD_W = 27.1 + 2 * (_OLD_GAP - GAP)
R = 0.8                                 # his corner radius
_KAPPA = 0.4477                         # circle-to-Bezier, as Illustrator writes
_RECEPTOR = ("Coding_block-4", "Coding_block-3", "Coding_block-2",
             "Coding_block", "Kinase_fill", "RNase_fill", "Single_SMS_exterior")


def _rounded(x0, x1, yt, yb, left_round, right_round):
    """A block, with either end square or rounded at radius R."""
    c = R * _KAPPA
    d = [f"M{x0 + (R if left_round else 0):g} {yt:g}"]
    if right_round:
        d += [f"L{x1 - R:g} {yt:g}",
              f"C{x1 - c:g} {yt:g} {x1:g} {yt + c:g} {x1:g} {yt + R:g}",
              f"L{x1:g} {yb - R:g}",
              f"C{x1:g} {yb - c:g} {x1 - c:g} {yb:g} {x1 - R:g} {yb:g}"]
    else:
        d += [f"L{x1:g} {yt:g}", f"L{x1:g} {yb:g}"]
    if left_round:
        d += [f"L{x0 + R:g} {yb:g}",
              f"C{x0 + c:g} {yb:g} {x0:g} {yb - c:g} {x0:g} {yb - R:g}",
              f"L{x0:g} {yt + R:g}",
              f"C{x0:g} {yt + c:g} {x0 + c:g} {yt:g} {x0 + R:g} {yt:g}", "Z"]
    else:
        d += [f"L{x0:g} {yb:g}", "Z"]
    return " ".join(d)


def _tick_and_label(text, x, yb, eid, label):
    """A label tick under a block, with its caption - his pA/P pattern.

    The two offsets are read off his own ticks: the tick runs from 1 below the
    block to 3.1 below it, and the caption's baseline sits 8.92 below it. That
    holds in each of his rows, which place the strand at a different height.
    """
    cls = re.search(r'<text id="Start_codon"[^>]*class="([^"]+)"', text)
    tick = re.search(r'class="(cls-\d+)"[^>]*/>\s*(?=<text id="Start_codon")',
                     text) or re.search(r'<path id="ATG_site"[^>]*class="(cls-\d+)"',
                                        text)
    return (f'<path id="Label_tick___{eid}" data-name="Label tick | {label}" '
            f'd="M{x:g} {yb + 1:g} L{x:g} {yb + 3.1:g}" '
            f'class="{tick.group(1) if tick else ""}"/>'
            f'<text id="Label___{eid}" data-name="Label | {label}" '
            f'text-anchor="middle" transform="matrix(1 0 0 1 {x:.3f} '
            f'{yb + 8.92:g})" x="0" y="0" xml:space="preserve" '
            f'class="{cls.group(1) if cls else ""}"><tspan x="0" y="0">'
            f'{label}</tspan></text>')


def _fuse_receptor(text):
    """Lay out Binder | linker | SMS, then the RNA scaffold target."""
    m = re.search(r'<path id="Coding_block-4"[^>]*\sd="([^"]+)"[^>]*/>', text)
    hp = re.search(r'id="Hairpin_encoding_outline"[^>]*\sd="([^"]+)"', text)
    if not m or not hp:
        return text
    n = [float(v) for v in re.findall(r'-?[\d.]+', m.group(1))]
    his_x0, yt, yb = min(n[0::2]), min(n[1::2]), max(n[1::2])
    ym = (yt + yb) / 2
    # the receptor starts one gap after the (already shortened) hairpin
    x0 = max(float(v) for v in re.findall(r'[MLC]\s*(-?[\d.]+)', hp.group(1))) + GAP
    sc0 = x0 + SMS_END + GAP
    sc1 = sc0 + SCAFFOLD_W
    O = f' stroke="{SYNBIO["outline"]}" stroke-width="0.4" ' \
        f'stroke-linecap="round" stroke-linejoin="round"'

    def block(i, j, fill, lr, rr, eid, name):
        d = _rounded(x0 + i, x0 + j, yt, yb, lr, rr)
        return f'<path id="{eid}" data-name="{name}" d="{d}" fill="{fill}"/>'

    parts = [
        # binder and SMS are one bar now: square where they meet, one seam
        block(0, BINDER_END, SYNBIO["binder"], True, False, "Binder", "Binder")
        .replace("/>", O + "/>"),
        block(LINKER_END, RNASE_START, SYNBIO["enzyme"], False, False,
              "SMS_kinase", "SMS kinase lobe"),
        block(RNASE_START, SMS_END, SYNBIO["enzyme_alt"], False, True,
              "SMS_rnase", "SMS RNase lobe"),
        f'<path id="SMS_outline" data-name="SMS outline" '
        f'd="{_rounded(x0 + LINKER_END, x0 + SMS_END, yt, yb, False, True)}" '
        f'fill="none"{O}/>',
        # the target: our RNA scaffold, a plain rounded box in theme 200
        # (constructs.py "scaffold"), because a structured UTR has no SBOL glyph
        f'<path id="RNA_scaffold" data-name="RNA scaffold" '
        f'd="{_rounded(sc0, sc1, yt, yb, True, True)}" '
        f'fill="{THEME["200"]}"{O}/>',
        # The scaffold gets no tick: its bracket sits over that one block and
        # already names it, and a caption there collides with the att label.
        _tick_and_label(text, x0 + BINDER_END / 2, yb, "Binder", "Binder"),
        _tick_and_label(text, x0 + (LINKER_END + SMS_END) / 2, yb, "SMS", "SMS"),
    ]
    text = text.replace(m.group(0), "".join(parts))
    for pid in _RECEPTOR[1:]:
        mm = re.search(rf'<path id="{re.escape(pid)}"[^>]*/>', text)
        if mm:
            text = text.replace(mm.group(0), "")
    # The ATG marks the start of the binder, so it travels with it. It keeps
    # Matteo's elbow - down from the block, ATG_REACH to the left, then the
    # tick - because the binder is now too narrow for both captions to sit
    # under it: the elbow is what carries the caption clear while still
    # pointing at the start codon.
    mm = re.search(r'(<path id="ATG_site"[^>]*\sd=")([^"]+)(")', text)
    if mm:
        ys = sorted(float(v) for v in
                    re.findall(r'[ML]\s*-?[\d.]+\s+(-?[\d.]+)', mm.group(2)))
        y0, y1 = ys[0], ys[-1]
        text = text.replace(
            mm.group(0),
            f'{mm.group(1)}M{x0:g} {y0:g} L{x0:g} {(y0 + y1) / 2:g} '
            f'L{x0 - ATG_REACH:g} {(y0 + y1) / 2:g} '
            f'L{x0 - ATG_REACH:g} {y1:g}{mm.group(3)}')
    mm = re.search(r'(<text id="Start_codon"[^>]*transform="matrix\(1 0 0 1 )'
                   r'[\d.-]+( [\d.-]+\)")', text)
    if mm:
        new = f'{mm.group(1)}{x0 - ATG_REACH:g}{mm.group(2)}'
        if 'text-anchor' not in mm.group(0):
            new = new.replace("<text", '<text text-anchor="middle"', 1)
        text = text.replace(mm.group(0), new)
    # one bracket per module: the sensor over the protein, the target over the
    # scaffold. The recorder's bracket belongs to the hairpin and is left alone.
    for bid, (i, j) in (("Sensor_module_bracket", (x0 + 0.2, x0 + SMS_END - 0.2)),
                        ("Target_module_bracket", (sc0 + 0.2, sc1 - 0.2))):
        mm = re.search(rf'(id="{bid}"[^>]*\sd=")([^"]+)(")', text)
        if not mm:
            continue
        ys = [float(v) for v in re.findall(r'[ML]\s*-?[\d.]+\s+(-?[\d.]+)',
                                           mm.group(2))]
        y1, y0 = max(ys), min(ys)
        text = text.replace(
            mm.group(0), f'{mm.group(1)}M{i:g} {y1:g} L{i:g} {y0:g} '
                         f'L{j:g} {y0:g} L{j:g} {y1:g}{mm.group(3)}')
    return text


def retheme(text, overrides):
    use = _classes(text)

    def fix_rule(m):
        cls, body = m.group(1), m.group(2)
        new = _fill_for(use.get(cls, set()), overrides)
        if new:
            body = re.sub(r'fill:#[0-9a-fA-F]{6}', f'fill:{new}', body)
        else:                                   # text and bare strokes
            for old, repl in TEXT_FILL.items():
                body = body.replace(f'fill:{old}', f'fill:{repl}')
        for old, repl in STROKE_MAP.items():
            body = body.replace(f'stroke:{old}', f'stroke:{repl}')
        return f'.{cls}{{{body}}}'

    text = re.sub(r'\.(cls-\d+)\{([^}]*)\}', fix_rule, text)
    # Illustrator writes `font-size:5.5`, which is invalid CSS - a length needs
    # a unit - so a browser throws the declaration away and falls back to 16px.
    # His own renderer tolerates it; everything else draws giant text.
    text = re.sub(r'font-size:([\d.]+)(?![\d.a-z%])', r'font-size:\1px', text)
    text = _shorten_hairpin(text)      # frees the room the scaffold needs
    text = _fuse_receptor(text)
    text = _relabel(text)
    return _prune_styles(text)


def _prune_styles(text):
    """Drop CSS rules no element uses any more.

    The receptor's blocks are gone, so their rules would otherwise sit in the
    <style> block carrying Matteo's hexes - invisible, but enough to make the
    leftover-colour check in main() fire, and enough to mislead anyone reading
    the file.
    """
    body = re.sub(r'<style>.*?</style>', '', text, flags=re.S)
    live = set(re.findall(r'class="(cls-\d+)"', body))
    return re.sub(r'\.(cls-\d+)\{[^}]*\}\s*',
                  lambda m: m.group(0) if m.group(1) in live else "", text)


def _bracket_centre(text, label_id):
    """x midpoint of the bracket a module label belongs to."""
    stem = label_id.replace("_label", "_bracket")
    m = re.search(rf'id="{re.escape(stem)}"[^>]*\sd="([^"]+)"', text)
    if not m:
        return None
    xs = [float(v) for v in re.findall(r'[ML]\s*(-?[\d.]+)', m.group(1))]
    return (min(xs) + max(xs)) / 2 if xs else None


def _relabel(text):
    """Swap label text, and re-centre a module label on its own bracket.

    His labels are placed by a fixed x that suited the old wording; a longer
    word would simply run into its neighbour, so anything sitting on a bracket
    is re-anchored to that bracket's midpoint instead.
    """
    def fix(m):
        head, body, tail = m.group(1), m.group(2), m.group(3)
        new = LABELS.get(body.strip())
        if new is None:
            return m.group(0)
        eid = re.search(r'id="([^"]+)"', head)
        cx = _bracket_centre(text, eid.group(1)) if eid else None
        if cx is not None:
            head = re.sub(r'transform="matrix\(1 0 0 1 [\d.-]+ ',
                          f'transform="matrix(1 0 0 1 {cx:.3f} ', head)
            if 'text-anchor' not in head:
                head = head.replace("<text", '<text text-anchor="middle"', 1)
        return head + new + tail

    return re.sub(r'(<text\b[^>]*>\s*<tspan[^>]*>)([^<]*)(</tspan>)', fix, text)


def main(src, dst):
    dst.mkdir(parents=True, exist_ok=True)
    for f in sorted(src.glob("*.svg")):
        if f.name.startswith("._"):
            continue
        over = (HYBRID_SITES if "Integrated" in f.name else
                PAD_SITES if "landing-pad" in f.name else {})
        out = retheme(f.read_text(), over)
        stray = {h.lower() for h in re.findall(r'#[0-9a-fA-F]{6}', out)}
        allowed = {c.lower() for c in (*THEME.values(), *NEUTRAL.values(),
                                       *SYNBIO.values(), *DIVERGING_WARM,
                                       INK, SURFACE)}
        left = sorted(stray - allowed)
        (dst / f.name).write_text(out, encoding="utf-8")
        print(f"{f.name:52} {'ok' if not left else 'LEFTOVER ' + str(left)}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
