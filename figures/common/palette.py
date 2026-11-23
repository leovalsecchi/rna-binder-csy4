"""Thesis colour palette (v0.1). Single source of truth: config/palette.json.

Guideline: never hand-pick colours in a figure script. Use these names.
See COLORS.md at the repo root.
"""
import json
from .paths import ROOT

_P = json.loads((ROOT / "config" / "palette.json").read_text())

THEME = _P["theme_petrol"]      # "100".."900"; 100-300 tints for fills only, 400-900 for marks
NEUTRAL = _P["neutral"]         # "100".."900"
DIVERGING = _P["diverging"]     # 7 colours, grey midpoint
STACKS = _P["stacks"]           # long stack name -> hex, fixed order
STACK_ORDER = list(STACKS)      # never cycle, never reorder
# Re-predictors (palette v0.6). A refolder is a measurement instrument, never a
# design route, so it has its own colours and must NOT borrow a stack colour.
# Fixed pipeline order, never reordered by rank. See COLORS.md rule 7.
REFOLDERS = _P.get("refolders", {})
REFOLDER_ORDER = list(REFOLDERS)
# Designed RNA scaffold architectures (palette v0.7). Four steps of ONE hue,
# because they are four variants of the same kind of object; F30 is NOT one of
# them - it keeps the main figures' THEME["200"] with its Csy4 sites in the
# diverging warm end. See COLORS.md rule 9.
SCAFFOLDS = _P.get("scaffolds", {})
SCAFFOLD_ORDER = list(SCAFFOLDS)
F30 = THEME["200"]              # the F30 array, as drawn in the Fig 1 construct maps
MAIN = _P["main"]               # default single-series colour (theme 700)
INK = _P["ink"]
GRID = _P["grid"]
SURFACE = _P["surface"]
# Molecular identity in schematics and imported cartoons (rule 6): never a data
# series, never a stack colour. Always drawn with SYNBIO["outline"] around it,
# which is what keeps a pale fill from reading as a data mark.
SYNBIO = _P.get("synbio", {})

# Interaction pseudobonds (contacts, H-bonds). Annotation, not a data
# series: never a stack colour, which means stack identity only.
ANNOTATION = _P.get("annotation", {})

SHORT = {
    "RFD3 + LigandMPNN": "RFD3+LiMPNN",
    "BoltzDesign + LigandMPNN": "BD1+LiMPNN",
    "HalluDesign + LigandMPNN": "HD+LiMPNN",
    "BoltzGen + BoltzIF": "BG+BIF",
    "BoltzGen + LigandMPNN": "BG+LiMPNN",
}
STACKS_SHORT = {SHORT[k]: v for k, v in STACKS.items()}


def stack_color(name: str) -> str:
    """Colour for a stack, by long or short name."""
    return STACKS.get(name) or STACKS_SHORT[name]


def register_colormaps():
    """Register 'thesis_seq' (theme 400-900) and 'thesis_div' with matplotlib."""
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    ramp = [THEME[k] for k in sorted(THEME, key=int) if int(k) >= 400]
    for name, cols in (("thesis_seq", ramp), ("thesis_div", DIVERGING)):
        mpl.colormaps.register(LinearSegmentedColormap.from_list(name, cols), force=True)
