"""config.py — single source of truth for the Csy4 RNA-binder plotting suite.

Thresholds, stack colors, refold conditions, and paths live ONLY here. If any
plot_NN_*.py hardcodes a number/color/path that belongs here, that's a bug.
"""
from pathlib import Path

# ---------------------------------------------------------------- paths
RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
ANALYSES = RB / "analyses"
TOOLS = ANALYSES / "plots_tools"          # this suite's code + sidecars live here

# Which aggregated dataset the suite is run on -> names the output subfolder
# analyses/plots/<DATASET>/{png,html}. Change with the data.
DATASET = "combined_v2_v3"

# Master per-design metrics table (the 3k stratified sample of that dataset).
# Repoint at combined_clean_metrics.csv to generalize to the full ~900k set.
MASTER_CSV = ANALYSES / "cross_campaign" / "combined_sample10k_metrics.csv"

# Optional sidecars merged in by design_id if present (else metrics -> NaN).
DOCKQ_SIDECAR = TOOLS / "sidecars" / "dockq_10k.csv"        # design_id,dockq_af3_msa_off,dockq_b2_msa_off
HBOND_NT_SIDECAR = TOOLS / "sidecars" / "hbond_nts_10k.csv"  # design_id,refold_condition,hbond_nts
# scaffold_fit: does the binder dock the hairpin WITHIN the full 121-nt construct without
# clashing? (scaffold_fit.py — superpose the designed epitope onto a scaffold hairpin, score).
SCAFFOLD_FIT_SIDECAR = TOOLS / "sidecars" / "scaffold_fit_10k.csv"  # design_id,scaffold_clashes,scaffold_interface_frac,scaffold_fits,...

# Outputs: analyses/plots/<DATASET>/{png,html}
OUTDIR = ANALYSES / "plots" / DATASET
PNG_DIR = OUTDIR / "png"
HTML_DIR = OUTDIR / "html"

# Best-designs ranking + export (plots 01_1 / 01_2 + best_binders/). Ranked by consensus
# ipTM = min(af3, b2, esm) (a design ranks high only if ALL THREE refolders agree); only
# designs with all three refolds are eligible.
RANK_COL = "consensus_iptm"
TOP_N_ACROSS = 10       # plot_01_1: global top-N designs table
TOP_N_PER_STACK = 2     # plot_01_2: best-N designs per stack
BEST_BINDERS_DIR = ANALYSES / "best_binders"   # extracted CIFs, uncompressed, rank-named

# ---------------------------------------------------------------- thresholds
# Referenced by every plot that filters. Change here only.
THRESHOLDS = {
    "iptm_min": 0.8,
    "dockq_min": 0.4,
    "ptm_min": 0.65,
    "disorder_max": 0.2,
    "clashes_max": 0,
    "hbonds_min": None,   # TBD — funnels skip this filter while None; trivially set later.
    # tip threshold on tip_weighted_score for the quality plots' "pass all 3" box (plots 14/15).
    # Distinct from the official tip_pass (= tip_bottom_fraction ≤ cutoff); tune to taste.
    "tip_min": 40,
    # scaffold_fit "pass both criteria" box (plot_18): interface predominantly on the hairpin
    # AND ~no clashes with the rest of the construct. Mirror scaffold_fit.py's constants.
    "scaffold_clash_max": 3,
    "scaffold_interface_min": 0.80,
}

# ---------------------------------------------------------------- plot_21 elite funnel
# plot_21_funnel_topx tunables (ALL here). top_x = keep each refolder's top-X designs by its
# OWN ipTM (global, across all stacks); a design passes only if it is in ALL THREE refolders'
# top-X. dockq_min = consensus DockQ = min(AF3,B2,ESM) hard filter. iptm_min = optional
# consensus ipTM = min(AF3,B2,ESM) pre-filter stage (0 = off). ipTM/DockQ thresholds elsewhere
# live in THRESHOLDS above.
TOPX_FUNNEL = {
    "top_x": 1000,
    "dockq_min": THRESHOLDS["dockq_min"],
    "iptm_min": 0.0,
}

# ---------------------------------------------------------------- funnel stages
# STAGGERED cumulative funnel: each (metric, [thresholds], direction) contributes ONE
# step per threshold, so a metric can tighten in stages (e.g. ipTM 0.5 -> 0.65 -> 0.8
# before DockQ). Only stages whose metric has data are shown (others slot in once
# computed). The final ipTM/DockQ stop equals THRESHOLDS[iptm_min]/[dockq_min]. Tune here.
FUNNEL_STAGES = [
    ("iptm", [0.5, 0.65, 0.8], "ge"),
    ("dockq", [0.4], "ge"),
    ("fraction_disordered", [0.2], "le"),
    ("n_clashes", [0], "le"),
]

# ---------------------------------------------------------------- refold conditions
# MSA-OFF ONLY. AF3 MSA-on is filtered out at the data_loader level and must not
# appear in any plot. The two conditions compared throughout (plots 7 & 10):
REFOLD_CONDITIONS = ["af3_msa_off", "b2_msa_off", "esm_msa_off"]
REFOLD_LABELS = {"af3_msa_off": "AF3 (MSA-off)", "b2_msa_off": "Boltz-2 (MSA-off)", "esm_msa_off": "ESMFold2"}
REFOLD_MARKERS = {"af3_msa_off": "o", "b2_msa_off": "s", "esm_msa_off": "^"}  # secondary encoding (relief)

# Cross-refolder ipTM comparison (plot 7 only). ESMFold2 is MSA-off too but is NOT in
# REFOLD_CONDITIONS (kept out of the filtering/funnel scope per spec) -- exposed here
# only so plot 7 can compare which refolder scores highest per stack.
REFOLDER_IPTM_COLS = {"af3": "af3_msaoff_iptm", "b2": "b2_iptm", "esm": "esm_iptm"}
REFOLDER_ORDER = ["af3", "b2", "esm"]
REFOLDER_LABELS = {"af3": "AF3 (MSA-off)", "b2": "Boltz-2 (MSA-off)", "esm": "ESMFold2"}
REFOLDER_COLORS = {"af3": "#0072B2", "b2": "#E69F00", "esm": "#009E73"}  # CVD-safe trio

# ---------------------------------------------------------------- stacks + colors
# Okabe-Ito (colorblind-safe). The two BoltzGen stacks share a generator -> two
# blues (dark + sky), distinguishable by lightness. Fixed order, never cycled.
STACK_COLORS = {
    "rfd3+ligandmpnn":        "#E69F00",   # orange
    "boltzdesign+ligandmpnn": "#009E73",   # bluish-green
    "halludesign+ligandmpnn": "#D55E00",   # vermillion
    "boltzgen+boltzif":       "#0072B2",   # blue  (shade A, darker)
    "boltzgen+ligandmpnn":    "#56B4E9",   # sky   (shade B, lighter)
}
STACK_ORDER = list(STACK_COLORS)

STACK_DISPLAY = {
    "rfd3+ligandmpnn":        "RFD3 + LigandMPNN",
    "boltzdesign+ligandmpnn": "BoltzDesign + LigandMPNN",
    "halludesign+ligandmpnn": "HalluDesign + LigandMPNN",
    "boltzgen+boltzif":       "BoltzGen + BoltzIF",
    "boltzgen+ligandmpnn":    "BoltzGen + LigandMPNN",
}

# Raw stack labels as they appear in the aggregated CSV -> canonical keys above.
STACK_ALIASES = {
    "RFD3+LiMPNN":  "rfd3+ligandmpnn",
    "BD1+LiMPNN":   "boltzdesign+ligandmpnn",
    "HD+LiMPNN":    "halludesign+ligandmpnn",
    "BG+BIF":       "boltzgen+boltzif",
    "BG+LiMPNN":    "boltzgen+ligandmpnn",
}

# ---------------------------------------------------------------- RNA target
RNA_NT_COUNT = 16   # csy4 crRNA 16-nt epitope numbering
RNA_SEQ = "ACUGCCGUAUAGGCAG"              # 16-nt epitope, for heatmap x-tick labels
AA_ORDER = list("RHKDESTNQCGPAVILMFYW")   # charged/polar -> hydrophobic (heatmap y axis)
HEATMAP_CMAP = "magma"                    # sequential, perceptually uniform
INTERFACE_AA_NT_MAP = TOOLS / "sidecars" / "aa_nt_map_af3_10k.csv"  # aa x nt map (hbonds + contacts)

# ---------------------------------------------------------------- style tokens
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#ffffff"
MARK_EDGE = "#33322f"       # thin dark stroke on marks (relief for sub-3:1 fills)
THRESHOLD_LINE = "#d03b3b"  # critical-red for threshold reference lines
FONT_FAMILY = ["system-ui", "-apple-system", "Segoe UI", "DejaVu Sans", "sans-serif"]


def color(stack: str) -> str:
    return STACK_COLORS.get(stack, INK_MUTED)
