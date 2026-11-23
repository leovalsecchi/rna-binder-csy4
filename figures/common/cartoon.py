"""Lift one drawing out of a snapshot of Matteo's common cartoons, recoloured.

His art is flat vector with every part named (`data-name="Outer contour"`,
`data-name="5-prime exon"`), so a drawing can be pulled out by name and its
colours mapped onto the thesis palette. The snapshots live in `inputs/cartoons/`
with their provenance in SOURCES.md; his own working tree is revised
continuously and must never be read by a figure script.

    from common import cartoon
    svg_fragment = cartoon.fragment("v19.18.svg", "HP v19.18 Spliced",
                                    scale=2.0, dx=120, dy=40)

An unmapped source colour raises. That is deliberate: it is what stops a hex
from Matteo's palette reaching a thesis figure unnoticed (COLORS.md rule 1).
Add it to config/cartoon_colors.json instead.

`bbox()` gives a conservative box (Bezier control points bound their curve), so
`fragment(..., at=(x, y), anchor="center")` places a drawing where you want it
without counting coordinates on Matteo's A4 pages.
"""
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from . import palette
from .paths import ROOT

SVG_NS = "http://www.w3.org/2000/svg"
SOURCES = ROOT / "inputs" / "cartoons"
_CFG = json.loads((ROOT / "config" / "cartoon_colors.json").read_text())
_COLORS = {k.upper(): v for k, v in _CFG["colors"].items()}
_PER_FILE = {f: {k.upper(): v for k, v in m.items()}
             for f, m in _CFG.get("files", {}).items()}

ET.register_namespace("", SVG_NS)

_SINGLE = {"ink": palette.INK, "main": palette.MAIN, "surface": palette.SURFACE}

# Every colour an import is allowed to end up with. Anything else in a finished
# fragment means a source colour slipped past the map - see _assert_in_palette.
_ALLOWED = {c.lower() for c in
            (*palette.THEME.values(), *palette.NEUTRAL.values(),
             *palette.SYNBIO.values(), *palette.DIVERGING,
             palette.INK, palette.GRID, palette.SURFACE, palette.MAIN)}
_HEX = re.compile(r"#[0-9A-Fa-f]{3,6}\b")
_TABLES = {"theme": palette.THEME, "neutral": palette.NEUTRAL,
           "synbio": palette.SYNBIO}


def _resolve(token):
    """'neutral.700' / 'ink' / '#015d82' -> hex."""
    if token.startswith("#"):
        return token
    group, _, step = token.partition(".")
    if group in _TABLES:
        return _TABLES[group][step]
    return _SINGLE[token]


def recolour(value, src=None):
    """One source colour to its thesis colour. Raises on an unmapped hex.

    `src` is the snapshot's file name, so a drawing can override a shared hex.
    """
    key = value.strip().upper()
    over = _PER_FILE.get(Path(src).name) if src else None
    if over and key in over:
        return _resolve(over[key])
    if key not in _COLORS:
        raise KeyError(f"{value} is not in config/cartoon_colors.json; add it "
                       f"rather than letting a foreign hex into a figure")
    return _resolve(_COLORS[key])


def _assert_in_palette(svg, src, name):
    """No hex may survive an import unless the palette owns it.

    Attribute recolouring alone is not enough: Illustrator can export colours in
    a CSS <style> block instead, where nothing sets fill= or stroke= and the
    recolour pass has nothing to rewrite. That export silently kept Matteo's
    own colours until this check was added.
    """
    stray = sorted({h for h in _HEX.findall(svg) if h.lower() not in _ALLOWED})
    if stray:
        raise ValueError(
            f"{src}{' / ' + name if name else ''}: {stray} is not a thesis "
            f"colour. The drawing probably carries a CSS <style> block rather "
            f"than fill/stroke attributes; recolour it by hand or leave it out.")


_CSS_RULE = re.compile(r"\.([A-Za-z0-9_-]+)\s*\{([^}]*)\}")


def _inline_css(root):
    """Move a <style> block's declarations onto the elements that use them.

    Illustrator exports some files with every fill in a CSS block and nothing
    on the elements. Lifting one drawing out of such a file would drop its
    colours on the floor, and the recolour pass, which reads attributes, would
    have nothing to rewrite - that is how an orange of his once reached a
    figure. Inlining first makes both work, whichever way the file was saved.
    """
    style_tag = f"{{{SVG_NS}}}style"
    rules = {}
    for st in root.iter(style_tag):
        for cls, body in _CSS_RULE.findall(st.text or ""):
            rules[cls] = dict(tuple(x.strip() for x in d.split(":", 1))
                              for d in body.split(";") if ":" in d)
    if not rules:
        return root
    for el in root.iter():
        for cls in (el.get("class") or "").split():
            for k, v in rules.get(cls, {}).items():
                el.set(k, v)          # CSS wins over a presentation attribute
        el.attrib.pop("class", None)
    for parent in root.iter():
        for child in list(parent):
            if child.tag == style_tag:
                parent.remove(child)
    return root


def _root(src):
    path = Path(src)
    return _inline_css(
        ET.parse(path if path.is_absolute() else SOURCES / path).getroot())


def _find(root, name):
    """Element whose data-name or id is `name`. Exact first, then prefix."""
    nodes = list(root.iter())
    for attr in ("data-name", "id"):
        for el in nodes:
            if el.get(attr) == name:
                return el
        for el in nodes:
            v = el.get(attr)
            if v and v.startswith(name):
                return el
    raise LookupError(f"no element named {name!r}")


_NUM = re.compile(r"-?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _points(el):
    """Every coordinate pair an element contributes, in its own user units."""
    tag = el.tag.split("}")[-1]
    if tag == "path" and el.get("d"):
        n = [float(v) for v in _NUM.findall(el.get("d"))]
        # path data is always x,y pairs except for arcs, which his art has none of
        return list(zip(n[0::2], n[1::2]))
    g = lambda k, d=0.0: float(el.get(k, d))
    if tag == "rect":
        x, y = g("x"), g("y")
        return [(x, y), (x + g("width"), y + g("height"))]
    if tag in ("circle", "ellipse"):
        cx, cy = g("cx"), g("cy")
        rx = g("r") or g("rx")
        ry = g("r") or g("ry")
        return [(cx - rx, cy - ry), (cx + rx, cy + ry)]
    return []


def bbox(src, name=None):
    """(x0, y0, x1, y1) around a drawing, in the source file's units."""
    root = _root(src)
    el = root if name is None else _find(root, name)
    pts = [p for node in el.iter() for p in _points(node)]
    if not pts:
        raise LookupError(f"nothing to measure in {name!r}")
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


def fragment(src, name=None, scale=1.0, dx=0.0, dy=0.0, at=None,
             anchor="topleft"):
    """SVG string for one drawing, recoloured, optionally scaled and moved.

    `src` is a file in inputs/cartoons/ (or any path). `name` is a data-name or
    id; the whole file is taken when it is None. Give `at=(x, y)` to place the
    drawing's `anchor` ("topleft" or "center") there, instead of dx/dy.
    """
    if at is not None:
        x0, y0, x1, y1 = bbox(src, name)
        if anchor == "center":
            x0, y0 = (x0 + x1) / 2, (y0 + y1) / 2
        elif anchor != "topleft":
            raise ValueError(f"anchor must be topleft or center, not {anchor!r}")
        dx, dy = at[0] - x0 * scale, at[1] - y0 * scale
    root = _root(src)
    el = root if name is None else _find(root, name)

    for node in el.iter():
        for attr in ("fill", "stroke", "stop-color"):
            v = node.get(attr)
            if v and v.startswith("#"):
                node.set(attr, recolour(v, src))

    body = "".join(ET.tostring(child, encoding="unicode") for child in el) \
        if el is root else ET.tostring(el, encoding="unicode")
    body = body.replace(f' xmlns="{SVG_NS}"', "")
    _assert_in_palette(body, Path(src).name, name)
    if (scale, dx, dy) == (1.0, 0.0, 0.0):
        return body
    return (f'<g transform="translate({dx},{dy}) scale({scale})">{body}</g>')


def standalone(src, name=None, scale=1.0, pad=4):
    """The same drawing as a complete SVG file, for a quick look."""
    x0, y0, x1, y1 = bbox(src, name)
    w, h = (x1 - x0) * scale + 2 * pad, (y1 - y0) * scale + 2 * pad
    return (f'<svg xmlns="{SVG_NS}" width="{w:.1f}" height="{h:.1f}" '
            f'viewBox="0 0 {w:.1f} {h:.1f}">'
            f'<rect width="{w:.1f}" height="{h:.1f}" fill="{palette.SURFACE}"/>'
            + fragment(src, name, scale, at=(pad, pad)) + '</svg>')


def _snapshots():
    """The snapshot files, skipping the ._ forks the SMB share leaves behind."""
    return sorted(f for f in SOURCES.glob("*.svg") if not f.name.startswith("._"))


def _selfcheck():
    """Every snapshot maps cleanly, and no source hex survives an import."""
    import re
    seen = set()
    for f in _snapshots():
        # scan the whole file, not just fill/stroke attributes: a CSS <style>
        # block hides colours from an attribute-only regex
        found = {h.upper() for h in _HEX.findall(f.read_text())}
        missing = found - set(_COLORS)
        assert not missing, f"{f.name}: unmapped {sorted(missing)}"
        seen |= found
    assert seen, "no snapshots found in inputs/cartoons/"

    assert recolour("#41B3D8", "SMS_v3.svg") == palette.SYNBIO["enzyme"]
    assert recolour("#41B3D8") != palette.SYNBIO["enzyme"], "override leaked"

    out = fragment("v19.18.svg", "HP v19.18 Spliced", scale=2.0, dx=10, dy=5)
    assert "<path" in out and "scale(2.0)" in out
    for h in seen:
        assert h.lower() not in out.lower(), f"{h} survived the recolour"
    assert recolour("#6c6970") == palette.NEUTRAL["700"]

    x0, y0, x1, y1 = bbox("SMS_v3.svg")
    assert x1 > x0 and y1 > y0, "degenerate bbox"
    placed = fragment("SMS_v3.svg", scale=2.0, at=(100, 50))
    assert "translate(" in placed
    try:
        recolour("#123456")
    except KeyError:
        pass
    else:
        raise AssertionError("an unmapped colour must raise")
    print(f"ok: {len(_snapshots())} snapshots, "
          f"{len(seen)} source colours, all mapped")


if __name__ == "__main__":
    _selfcheck()
