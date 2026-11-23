"""Extended metric panel: the campaign-scale metrics the main text does NOT plot, and
the one stored-vs-corrected pair that exists.

One bare artefact. Six small multiples in 2 x 3 are what the audit left: five over the
whole campaign, all 6,992,515 designs, no sampling, and apo RMSD, which exists for the
v2-v7 corpus only (3,934,826 designs, 0 missing) and carries its own n on the canvas.
Under them a two-panel STRIP - sc_DockQ and ligand_RMSD - which the main text DOES plot
which reappear here for one reason: they are the only two metrics in the campaign with a
BoltzGen+BoltzIF true-backbone counterpart IN EXISTENCE, so the strip is the only place
a reader sees the size of the scoring defect against the campaign's other metrics.
Restored 2026-11-19 on Leo's comment, "why did you remove the pre-fix bg+bif. I thought
it was informative". He is right; the rebuild had dropped it with the two metrics.

REBUILT 2026-11-19 on Leo's comment: "most of these metrics we already have in the
main panel. can you show the metrics we have not shown? ... I think we have a lot like
apo_rmsd or delta sasa etc. buns etc." The previous cut drew ipTM, sc_DockQ, ligand_RMSD,
buried dSASA total, Rg and foot clearance; five of those six are in the main text
(Fig 2c, Fig 3a, Fig 3b and the built-but-unpromoted Fig 5d), so the figure said
almost nothing new. README carries the full audit of what is plotted where.

WHAT IS DRAWN, AND WHY THESE SIX - every reason is a measured number, printed by
report() and recorded in the README:
  packstat                 route spread 0.76 x the campaign IQR, and r = -0.328
                           against dsasa_total, so it is not interface size restated.
                           (Both figures are what report() PRINTS, with the exact zeros
                           dropped as the panel drops them. The candidate sweep that kept
                           the zeros said 0.77 and -0.136; this docstring quoted those
                           until 2026-11-19, so it disagreed with its own output.)
  delta_unsat_hbonds       BUNS, which Leo named. The Rosetta InterfaceAnalyzer delta -
                           bound minus unbound - so it IS an interface quantity. It is
                           also the one drawn metric that fails this figure's own
                           independence test (r = +0.628 against dsasa_total), and that
                           is stated rather than dodged: there is no size-free member of
                           this family. The alternative, buns_vsasa, has a larger route
                           separation (1.15 x IQR against 1.00) and a smaller r (+0.368),
                           but ONLY because set_use_ddG_style(True) returns 0.0 for every
                           protein-RNA pose in this PyRosetta build, so buns_vsasa is a
                           WHOLE-POSE absolute count, not a bound-minus-unbound delta
                           (rna-binder-buns-metric). Picking it would have bought the
                           prettier number by quietly swapping the quantity.
  frac_dsasa_base          dSASA, base share. r = -0.305 against dsasa_total.
  frac_dsasa_major         dSASA, major-edge share of the base area. r = +0.382.
  red_total_per_ifaceres   unfavourable contacts; r = +0.022, the one candidate fully
                           independent of interface size, and monotone across routes.
  apo_rmsd_vs_esm          Leo named it. 3,934,826 designs - NOT the 100,000 the 100k
                           table suggests; the 2026-07-29 re-launch finished and its
                           3,835 chunk CSVs were never merged onto the master. n on panel.

NOT DRAWN, each for a measured or file-level reason (README has the table):
  dsasa_polar/_hphobic, dsasa_rna_phosphate/_sugar/_base/_base_major - r = +0.69 to
    +0.98 against dsasa_total, which the main text already shows. They are the
    interface-size panel again. The two FRACTIONS above are the members that are not.
  sc_value          route spread only 0.20 x IQR; nothing to see at campaign scale.
  pi_stack_score    the worst-case consensus is 0 at the median of every route, with
    zero IQR - a violin of it is a line. (code/supp/interface_chemistry drew it per
    route until that figure was archived 2026-11-20; now drawn nowhere.)
  n_h_bonds_interface  medians of 1-2 H-bonds. (Same archived figure drew it per route
    and re-predictor; now drawn nowhere.)
  dockq_fnat, dockq_iRMSD  unplotted anywhere, but they are measured AGAINST THE
    DESIGNED BACKBONE, so they are wrong for BoltzGen+BoltzIF exactly as the stored
    sc_DockQ and ligand_RMSD were - and the true-backbone rescore never re-emitted them
    (see below). Drawing them would put a known-wrong 1.48 M-design series on the page
    with no correction available.
  ipSAE, Rosetta interface energy  do not exist on ANY campaign-scale file. Unchanged
    from the previous cut; README keeps the evidence.

THE STRIP - BOLTZGEN+BOLTZIF, STORED vs CORRECTED, and why exactly two metrics carry it.
That route's stored dockq_* and rmsd_* were scored against BoltzIF's own Boltz
re-prediction instead of against the generated backbone; the true-backbone rescore
(euler:.../analyses/bif_true_bb_rescore/full/bif_true_bb_scores.csv.gz) is the
correction. That file has 16 columns and emits exactly two quantities per re-predictor -
dockq_self and dockq_LRMSD. No other metric on the master has a *_true counterpart, and
that is a file-level fact rather than a judgement: the correction only touches
quantities measured AGAINST THE DESIGNED BACKBONE, while the six panels above are read
off the re-predicted complex (apo, off two ESMFold2 predictions), which the defect never
touched. So the pair is drawable on sc_DockQ and ligand_RMSD and on NOTHING ELSE here - said
plainly in the caption, because it is the reason only those two carry it.

Both ARE in the main text (Fig 2c, Fig 3a, Fig 3b), which breaks this figure's own rule
of drawing what the main text does not. They are re-admitted anyway and flagged as such
in the caption: the main text plots the CORRECTED values only, so nowhere else does a
reader see stored against corrected with four unaffected design routes beside them for
scale. STORED is open, CORRECTED is filled - fill versus outline plus adjacency, never
hatch (palette v0.6) - and both are in the BoltzGen+BoltzIF route colour, since a route
colour means route identity and nothing else. The six panels above keep the rule.

Sources, all local (see README for mtimes):
  inputs/large/metric_panel_v2.csv.gz      the five campaign metrics, per re-predictor,
                                           plus dsasa_total and buns_vsasa, which are
                                           not drawn but ARE reported
  inputs/large/apo_rmsd_all.csv.gz         apo RMSD, both apo runs
  inputs/large/campaign_metrics_full.csv.gz  read ONCE, for the consensus assertion AND
                                           for the strip's STORED sc_DockQ / ligand_RMSD
  inputs/large/bif_true_per_backend.csv.gz the strip's CORRECTED BoltzGen+BoltzIF values

Consensus rule, one rule for the whole figure, the strip included: the worst of the two
gate re-predictors, AF3 (MSA off) and ESMFold2. `min` where high is good, `max` where
low is good, and `mean` for the two composition fractions, which have no better or worse
end - the same exception Rg took in the previous cut. Boltz-2 is excluded, as it was
from the gate. The rule is ASSERTED, not asserted-of: it must reproduce the stored
consensus sc_dockq and sc_rmsd on RFdiffusion3+LigandMPNN to 1e-4, on a design route the
BoltzIF defect never touched, before anything is drawn. The strip applies the same rule
to the *_true columns, so stored and corrected are compared like with like.

No hatch anywhere (palette v0.6, 2026-11-19: the dots and stripes are withdrawn
repo-wide). Nothing in this figure distinguishes a re-predictor - every violin is a
design route, the strip's pair included - so no mark takes a palette.REFOLDERS colour
and every body is drawn in its own route colour from palette.stack_color().

SPLIT VIOLINS, 2026-11-20, in Fig 2c's idiom (Leo: "can you do the split violins for
metric_panel the same way as its done in the main in fig2"). Every body that HAS two
gate re-predictors behind it is now one violin split down its centre line: AF3 (MSA off)
is the LEFT half at full strength, ESMFold2 the RIGHT at alpha 0.5, both in the design
route's own colour, with one key in the grid's row gutter. So the panels draw the two
INPUTS to this figure's consensus rule rather than the consensus itself; the rule is
still asserted against a column somebody else computed, and still what report() tabulates
and what the README's verified medians are, so nothing is lost by it leaving the canvas.

APO RMSD IS NOT SPLIT, and that is the same decision as Fig 2c's first panel being a BAR:
"Designs" is a count with no distribution behind it, so a violin there would invent a
spread, and apo RMSD is ONE ESMFold2 measurement (the binder folded alone) with no AF3
counterpart in existence, so a split there would invent a re-predictor. It keeps one
whole symmetric body, and its n note now names the single re-predictor on the canvas.

Fig 2c's `_half()` is NOT imported, though the idiom is its. It draws horizontally
(`vert=False`), writes to module-global collectors for its own alignment check, and its
density is matplotlib's KDE over a 60,000-design subsample per route. The density here is
a hand-rolled histogram smoothed with a Hann window, computed in the panel's own axis
space (log where the axis is log), over every value - 1.4 M per route, which is what
makes a KDE unaffordable and the subsample necessary there. Importing it would have
changed what six panels draw in order to share nine lines. What IS shared is the
convention, named once in BACKENDS below and nowhere else in this file: which half, which
tone, which name. (The venv also has no usable scipy - _propack dlopen failure - so
every statistic below, median, IQR, Pearson r and the spread ratios, is hand-rolled.
matplotlib's violinplot itself does work: it uses its own mlab.GaussianKDE, not scipy's.)

    ~/.venvs/binder_map/bin/python3 code/supp/metric_panel/make_figure.py
    METRIC_PANEL_CACHE=/tmp/mp2.parquet ...   # optional, caches the merged frame
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
LARGE = paths.INPUTS / "large"
W_MM, H_MM = 178.0, 145.0
# panel box and grid, mm. PH is also the bound on a rotated y label (STYLE_NOTES §3).
# GY is 8 and not the previous cut's 17 because the main grid's rotated route names now
# hang below its bottom row only, into G_STRIP; the height that freed went into PH.
PW, PH, X0, GX, GY = 44.0, 33.0, 13.0, 12.0, 8.0
# The strip row. Two panels at 72 mm put its right edge at 13 + 2*72 + 12 = 169, the
# same as the main grid's 13 + 3*44 + 2*12, so the two blocks share a right margin
# while the strip's six positions get 12 mm each instead of the grid's 8.8.
PW_S = 72.0
YS = 18.5          # bottom margin: holds the strip's longest rotated label, measured
G_STRIP = 18.0     # gap above the strip: holds the MAIN GRID's rotated route names
TOP = 1.5
Y0 = YS + PH + G_STRIP     # bottom of the main grid's lower row
NBIN = 160                 # density bins per panel; ~1.4 M values per route
SMOOTH = 11                # Hann window, in bins
TRIM = 0.005               # draw the body only where the density is at least this
                           # share of its own peak, so a body ends where it still has
                           # width. It no longer carries the clamped-tail needle, which
                           # density() stopped making on 2026-11-19 - that was the clip,
                           # not the tail.
N_CAMPAIGN = 6_992_515

ROUTES = list(palette.STACKS_SHORT)        # fixed order, never sorted by value

# THE SPLIT, declared once. (column suffix, name, which half, tone). AF3 reads first,
# as it does in Fig 2c - there top-to-bottom, here left-to-right, because these violins
# stand up. The tone is a redundant cue, as it is in Fig 2c: the HALF carries the
# re-predictor and the alpha only makes the seam legible where two bodies crowd. Neither
# half takes a palette.REFOLDERS hue - every body here is a design route in its route
# colour, and a re-predictor colour beside them would read as a sixth route (COLORS.md
# rule 7).
BACKENDS = (("af3", "AF3", -1, 1.0), ("esm", "ESMFold2", +1, 0.5))
KEY_SW, KEY_SH = 2.4, 1.8                  # mm, one tone swatch of the key

# key, label, log?, (lo, hi), consensus rule, n annotation (None = the whole campaign)
METRICS = [
    # limits set from the 0.5th-99.5th percentile of each metric's own consensus,
    # measured on the campaign - not guessed. Since density() stopped clamping
    # (2026-11-19) a limit inside that range no longer piles the out-of-range mass
    # against the frame; it drops it off the panel, so a limit set too tight shows as a
    # body running off the frame at full width instead of as a flat top. The percentile
    # rule keeps that at 0.5% a side either way.
    ("packstat",               "Packstat",         False, (0.35, 0.85), "min",  None),
    ("delta_unsat_hbonds",     "Interface BUNS ↓", False, (0, 100),     "max",  None),
    ("red_total_per_ifaceres", "Unfav. contacts ↓", False, (0, 7),      "max",  None),
    # Both shares widened 2026-11-20 by check_limits() over the HALVES, which the mean
    # consensus had been hiding: it averages a 0 against a 0.4 and reports 0.2, so the
    # drawn halves reach further down than the body they replaced. base share 0.1 -> 0.05
    # (RFD3's AF3 half needs 0.0552) and 0.75 -> 0.77 (BG+BIF's ESMFold2 half needs
    # 0.758). major share 0.1 -> 0, where its AF3 halves have a point mass: 0.18-0.96% of
    # each route is EXACTLY 0. Unlike packstat's zeros that is a real composition value -
    # base area with no major-groove-edge share - not Rosetta's "no interface", so it is
    # drawn, at the floor of an axis that honestly starts at zero, rather than dropped.
    ("frac_dsasa_base",        "ΔSASA base share", False, (0.05, 0.77), "mean", None),
    ("frac_dsasa_major",       "ΔSASA major share", False, (0, 1.0),   "mean", None),
    ("apo_rmsd_vs_esm",        "Apo RMSD (Å) ↓",   True,  (0.06, 40),   "one",  "apo"),
]
# per-re-predictor columns exist for everything except apo, which is one ESMFold2
# measurement (binder folded alone against the ESM complex prediction) and so has no
# consensus. The column names on disk keep the historical _af3 / _esm suffixes.
# rule == "one" is therefore also what decides which panel is SPLIT: the five with two
# re-predictors behind them are, apo is not, and the test is the existence of the second
# column rather than a second list that could fall out of step with this one.
PER_BACKEND = [m[0] for m in METRICS if m[4] != "one"]

# --- the strip: the ONE stored-vs-corrected pair that exists --------------------
N_BGBIF = 1_480_140         # BoltzGen+BoltzIF designs; both strip sources must have this
BIF, BIF_ST, BIF_TR = "BG+BIF", "BG+BIF stored", "BG+BIF corrected"
# Six positions, not five: BoltzGen+BoltzIF appears twice, in its own slot in the fixed
# route order (COLORS.md rule 2), stored first. Same colour for both - the pair is
# distinguished by FILL vs OUTLINE and by adjacency, never by hue and never by hatch.
STRIP_X = [r if r != BIF else BIF_ST for r in ROUTES]
STRIP_X.insert(STRIP_X.index(BIF_ST) + 1, BIF_TR)
# key, label, log?, (lo, hi). Limits on the grid's rule - the 0.5th to 99.5th percentile
# - but taken PER SERIES and not over the pooled campaign, which the first render of this
# strip got wrong. Pooled, sc_DockQ's p99.5 is 0.665 and ligand_RMSD's 66.1; but these six
# series are far more unlike each other than the grid's five routes, so a pooled limit is
# set by the crowd and cuts the heavy-tailed series off inside its own body - 1.02% of
# RFdiffusion3's sc_DockQ above 0.7, 1.83% of HalluDesign's ligand_RMSD above 70 would simply
# leave the panel, which draws as a body running off the frame at full width rather than
# as a tail closing in it. Widest series p99.5 is 0.7358
# and 82.11, so the limits are 0.75 and 85, and main() ASSERTS that every drawn series
# fits inside them rather than leaving it to this comment.
# sc_DockQ takes a LOG axis. Its mass is all in the bottom decade (route medians 0.021
# 0.081 under a limit of 0.75), so on a linear axis the whole point of the strip - stored
# 0.058 against corrected 0.020 - is a 5% sliver near the frame; on this log axis it is
# a fifth of the panel. Same argument as the apo panel, and a violin encodes by width,
# not by length from zero, so STYLE_NOTES §8 is not in the way. ligand_RMSD needs no log:
# its medians sit mid-axis and the shift is 30.5 -> 40.2 A.
STRIP = [
    # 0.75 -> 0.85 on 2026-11-20: the halves reach where the min consensus did not look.
    # RFdiffusion3's ESMFold2 half needs 0.8308 against the consensus series' 0.7358.
    ("dockq", "sc_DockQ",      True,  (0.005, 0.85)),
    ("rmsd",  "ligand_RMSD (Å) ↓", False, (0, 85)),
]
# The strip's halves. Both sides of the pair have per-re-predictor columns - STORED in
# the campaign's own dockq_/rmsd_{af3,esm} and CORRECTED in the rescore's *_true twins -
# so the strip splits like the grid, and stored against corrected stays what it was:
# OUTLINE against FILL, in the one route colour, side by side in BoltzGen+BoltzIF's own
# slot. A tone cannot show on an unfilled body, so on the stored one the half alone
# carries the re-predictor; that is the cue the grid leans on too.
STRIP_PB = [f"{k}_{b}" for k in ("dockq", "rmsd") for b, *_ in BACKENDS]


def load():
    cache = os.environ.get("METRIC_PANEL_CACHE")
    if cache and os.path.exists(cache):
        return pd.read_parquet(cache), pd.read_parquet(cache + ".apo")
    # dsasa_total comes along so the "why these five" correlations below are RE-DERIVED
    # on every run rather than quoted from the README. It is not drawn: the main text
    # already shows it (Fig 5d) and the previous cut of this figure did too.
    use = (["stack"]
           + [f"{k}_{b}" for k in PER_BACKEND + ["dsasa_total", "buns_vsasa"]
              for b in ("af3", "esm")])
    d = pd.read_csv(LARGE / "metric_panel_v2.csv.gz", usecols=use,
                    dtype={"stack": "category",
                           **{c: "float32" for c in use if c != "stack"}})
    a = pd.read_csv(LARGE / "apo_rmsd_all.csv.gz",
                    dtype={"stack": "category", "apo_rmsd_vs_esm": "float32"})
    if cache:
        d.to_parquet(cache, index=False)
        a.to_parquet(cache + ".apo", index=False)
    return d, a


def campaign_once():
    """campaign_metrics_full, read ONCE and used twice.

    The consensus assertion needs the RFdiffusion3+LigandMPNN rows; the strip needs all
    five design routes' stored sc_DockQ and ligand_RMSD, plus the BoltzGen+BoltzIF join to
    true-backbone rescore. Returns (those RFD3 rows, the strip's long frame).
    """
    cache = os.environ.get("METRIC_PANEL_CACHE")
    ck_o, ck_s = (cache + ".assert", cache + ".strip") if cache else (None, None)
    if ck_o and os.path.exists(ck_o) and os.path.exists(ck_s):
        return pd.read_parquet(ck_o), pd.read_parquet(ck_s)
    c = pd.read_csv(LARGE / "campaign_metrics_full.csv.gz",
                    usecols=["design_id", "stack", "sc_dockq", "sc_rmsd", "dockq_af3",
                             "dockq_esm", "rmsd_af3", "rmsd_esm", "bif_rescored"],
                    dtype={"stack": "category"})
    assert len(c) == N_CAMPAIGN, f"strip source is not the campaign: {len(c):,}"
    o = c.loc[c["stack"] == "RFD3+LiMPNN",
              ["sc_dockq", "sc_rmsd", "dockq_af3", "dockq_esm",
               "rmsd_af3", "rmsd_esm"]].reset_index(drop=True)
    s = strip_frame(c)
    if ck_o:
        o.to_parquet(ck_o, index=False)
        s.to_parquet(ck_s, index=False)
    return o, s


def strip_frame(c):
    """The strip in long form: one row per drawn value, columns pos / dockq / rmsd /
    paired.

    STORED is the campaign's own per-re-predictor columns under this figure's consensus
    rule; CORRECTED is the SAME rule applied to the *_true columns of the true-backbone
    rescore, so the two sides are compared like with like. They are never mixed inside
    one number - code/supp/bif_effect owns that vocabulary and this figure keeps it.

    `paired` marks, on the two BoltzGen+BoltzIF positions, the designs carrying BOTH a
    stored and a corrected value; report() takes the common-subset medians over it, and
    those answer the obvious objection that the shift is a coverage artefact. It is
    False on the other four design routes, where there is no pair for it to mean.
    """
    c["dockq"] = np.fmin(c["dockq_af3"], c["dockq_esm"])
    c["rmsd"] = np.fmax(c["rmsd_af3"], c["rmsd_esm"])
    t = pd.read_csv(LARGE / "bif_true_per_backend.csv.gz")
    assert len(t) == N_BGBIF, f"true-backbone rescore has {len(t):,}, not {N_BGBIF:,}"
    t["dockq"] = np.fmin(t["dockq_af3_true"], t["dockq_esm_true"])
    t["rmsd"] = np.fmax(t["rmsd_af3_true"], t["rmsd_esm_true"])
    # The two HALVES of each strip body, carried through beside the consensus the strip
    # used to draw whole. Renamed to the plain per-re-predictor names so one drawing call
    # serves both sides of the pair: the suffix that says which side is _st / _tr below,
    # and `_true` never has to be spelled again past this line.
    t = t.rename(columns={f"{c_}_true": c_ for c_ in STRIP_PB})

    b = c.loc[c["stack"] == BIF, ["design_id", "dockq", "rmsd", *STRIP_PB, "bif_rescored"]]
    assert len(b) == N_BGBIF, f"BoltzGen+BoltzIF is {len(b):,}, not {N_BGBIF:,}"
    j = b.merge(t[["design_id", "dockq", "rmsd", *STRIP_PB]], on="design_id", how="inner",
                suffixes=("_st", "_tr"), validate="one_to_one")
    assert len(j) == N_BGBIF, f"the join dropped rows: {len(j):,} of {N_BGBIF:,}"
    # The master's own bif_rescored flag must mark EXACTLY the designs the rescore could
    # score. A cross-file check, and the only thing here that would catch one of these
    # two files being swapped for another vintage of itself.
    assert bool(((j["bif_rescored"] == 1) == j["dockq_tr"].notna()).all()), \
        "bif_rescored disagrees with the true-backbone file on which designs rescored"
    paired = j[["dockq_st", "rmsd_st", "dockq_tr", "rmsd_tr"]].notna().all(1).to_numpy()

    other = c.loc[c["stack"] != BIF, ["stack", "dockq", "rmsd", *STRIP_PB]].rename(
        columns={"stack": "pos"})
    other["pos"] = other["pos"].astype(str)
    other["paired"] = False
    s = pd.concat(
        [other] + [pd.DataFrame({"pos": p, "paired": paired,
                                 "dockq": j[f"dockq_{sfx}"].to_numpy(),
                                 "rmsd": j[f"rmsd_{sfx}"].to_numpy(),
                                 **{c_: j[f"{c_}_{sfx}"].to_numpy() for c_ in STRIP_PB}})
                   for p, sfx in ((BIF_ST, "st"), (BIF_TR, "tr"))],
        ignore_index=True)
    s["pos"] = pd.Categorical(s["pos"], STRIP_X)
    assert s["pos"].notna().all(), "a strip row landed outside the six drawn positions"
    return s


def assert_consensus_rule(o):
    """The figure's own min/max rule must reproduce the stored consensus on a design
    route the BoltzIF defect never touched, or the rule is wrong and every panel using
    it - the strip included - is wrong. Carried over unchanged from the previous cut: it
    is the only place in this repo that checks the worst-of-two convention against a
    column somebody else computed."""
    for stored, a, b, rule in (("sc_dockq", "dockq_af3", "dockq_esm", np.fmin),
                               ("sc_rmsd", "rmsd_af3", "rmsd_esm", np.fmax)):
        bad = int((o[stored] - rule(o[a], o[b])).abs().gt(1e-4).sum())
        assert bad == 0, f"consensus rule does not reproduce {stored} on RFD3: {bad} rows"
    print(f"consensus rule reproduces sc_dockq (min) and sc_rmsd (max) on "
          f"RFdiffusion3+LigandMPNN to 1e-4, {len(o):,} designs")


def derive(d):
    """One rule per metric. min where high is good, max where low is good, mean for the
    two composition fractions, which have no better or worse end."""
    rules = {m[0]: m[4] for m in METRICS}
    rules["dsasa_total"] = "min"          # not drawn; the reference for the r column
    rules["buns_vsasa"] = "max"           # not drawn; the whole-pose BUNS, reported only
    for key, rule in rules.items():
        if rule == "one":
            continue
        af3, esm = d[f"{key}_af3"].to_numpy(), d[f"{key}_esm"].to_numpy()
        d[key] = (np.fmin(af3, esm) if rule == "min" else
                  np.fmax(af3, esm) if rule == "max" else (af3 + esm) / 2)
    # A packstat of EXACTLY 0 is a different population from a badly packed interface:
    # Rosetta returns it where there is no packable interface at all. 45,313/6,992,515
    # (0.65%) of designs, concentrated in BoltzDesign1 (26,808/1,444,614, 1.86%) and
    # HalluDesign (15,244/1,023,007, 1.49%). Between 0 and 0.35 there are only 1,310
    # designs campaign-wide, so it is a point mass with a gap above it, and a violin
    # cannot draw a point mass - drawn from zero it becomes a hairline across the empty
    # gap. It is below this panel's 0.35 limit and so would leave the panel on its own
    # now that density() does not clamp, but it is dropped HERE rather than silently off
    # the frame, because that is what makes report()'s "n with a value" state it on
    # every run.
    # The CONSENSUS column and both HALVES, or the split quietly puts them back: the
    # panels draw the per-re-predictor columns now, and check_limits() caught it as a
    # 0.5th percentile of exactly 0 on BoltzDesign1's and HalluDesign's AF3 halves.
    # Dropped AFTER the consensus above is computed, so every number this figure reports
    # is unchanged - a zero that wins the min is still what made the body drop out.
    # ⚠️ The zeros are AF3's, which only the split could show: 45,310 of the 45,313 are
    # in packstat_af3 and 5 in packstat_esm. So Rosetta finds no packable interface in
    # 1.86% of BoltzDesign1's and 1.49% of HalluDesign's AF3 predictions and in
    # essentially none of anybody's ESMFold2 ones.
    for col in ["packstat"] + [f"packstat_{b}" for b, *_ in BACKENDS]:
        d.loc[d[col] == 0, col] = np.nan
    return d


def series(frame, route, key, col="stack"):
    """One drawn body's values. `col` is "pos" for the strip, whose six positions are
    not five design routes."""
    return frame.loc[frame[col] == route, key].to_numpy(float)


def density(v, lo, hi, log):
    """Smoothed histogram density on the panel's own scale. No scipy.

    Out-of-range values LEAVE THE PANEL. np.histogram already ignores whatever falls
    outside its bins, so nothing is clamped in. Until 2026-11-19 both branches passed
    np.clip(v, lo, hi), which piled every out-of-range value into the end bin; the Hann
    window then spread that pile inwards, and the part of it that cleared TRIM was drawn
    as a hairline running to the frame - a needle made by the clamp, not by data. It was
    visible on the BoltzGen+BoltzIF strip bodies and on unfavourable contacts.
    """
    v = v[np.isfinite(v)]
    if log:
        v = v[v > 0]
        g = np.linspace(np.log10(lo), np.log10(hi), NBIN + 1)
        h, _ = np.histogram(np.log10(v), bins=g)
    else:
        g = np.linspace(lo, hi, NBIN + 1)
        h, _ = np.histogram(v, bins=g)
    w = np.hanning(SMOOTH + 2)[1:-1]
    h = np.convolve(h.astype(float), w / w.sum(), mode="same")
    centres = (g[:-1] + g[1:]) / 2
    if log:
        centres = 10 ** centres
    return centres, h / (h.max() or 1.0)


def body(ax, x, v, lo, hi, log, color, side=0, alpha=1.0, half=0.40, filled=True):
    """One violin body, or ONE HALF of a split one.

    `side` is -1 for the left half, +1 for the right and 0 for a whole symmetric body -
    which only apo RMSD takes, having one re-predictor and so no second half to draw.
    Each half keeps the whole body's 0.40 reach, so the split changes what is inside the
    envelope and not the envelope: the row pitch, the limits and the label geometry are
    untouched by it.

    `filled=False` draws an outline, which is how the strip marks the superseded STORED
    values: fill versus outline, in the one route colour, never hatch (palette v0.6) and
    never a second hue, because a route colour means route identity. The median rule
    follows its own half - SURFACE reads on a filled body and vanishes on an open one, so
    an open body takes the rule in its own colour - and `alpha` is deliberately NOT
    applied to an outline, where it would leave the ESMFold2 half a half-visible
    hairline; there the half itself carries the re-predictor, through the SHAPE of its
    own curve, which is the cue every body on this canvas uses.
    """
    c, dn = density(v, lo, hi, log)
    keep = dn >= TRIM
    if keep.any():                      # contiguous body between the first and last
        i, j = np.argmax(keep), len(keep) - np.argmax(keep[::-1])
        c, dn = c[i:j], dn[i:j]
    reach = half * dn
    if filled:
        ax.fill_betweenx(c, x if side else x - reach, x + (side or 1) * reach,
                         facecolor=color, alpha=alpha, edgecolor=color,
                         linewidth=0.4, zorder=2)
    else:
        # The OUTLINE is the density curve alone, never a closed half. Closing it would
        # run a line down the centre of the body, and on the strip's open stored
        # position that stroke was the loudest mark in the panel - a stem the filled
        # bodies do not have, bisecting the one body whose halves are nearly identical.
        # Two curves meeting at the trim ends read as the hollow violin they did before
        # the split, with each side now traced on its own values. Capped at the ends
        # because at TRIM the curve is 0.5% of the half width from the centre line.
        ax.plot(x + (side or 1) * reach, c, color=color, lw=0.5, zorder=2,
                solid_capstyle="round")
        if not side:
            ax.plot(x - reach, c, color=color, lw=0.5, zorder=2, solid_capstyle="round")
    m = float(np.nanmedian(v))
    ax.plot([x if side else x - half * 0.55, x + (side or 1) * half * 0.55], [m, m],
            color=palette.SURFACE if filled else color,
            lw=0.6, zorder=3, solid_capstyle="butt")


def split(ax, x, get, lo, hi, log, color, filled=True):
    """Both halves of one position: AF3 then ESMFold2, from `get(suffix)`.

    One call site for the convention, so a panel cannot draw the halves in one order and
    the key explain them in another. A half with no finite value is skipped rather than
    drawn as a flat line at the axis floor.
    """
    for suf, _name, side, alpha in BACKENDS:
        v = get(suf)
        if np.isfinite(v).any():
            body(ax, x, v, lo, hi, log, color, side=side, alpha=alpha, filled=filled)


def half_key(fig, fs):
    """The key to the split: which half is which re-predictor.

    Fig 2c's key, turned on its side with its violins. The two tones TOUCH, in the order
    the halves are drawn, each name beside its own side, so the key is a miniature of the
    thing it explains rather than two loose swatches. Neutral INK at the two alphas and
    never a palette.REFOLDERS hue (COLORS.md rule 7, and the reason Fig 2c's key is
    neutral too): every body on this canvas is a design route in its route colour.

    It lives in the gutter BETWEEN the grid's two rows, left-aligned under the top-left
    panel. That band is the only empty full-width space on a canvas whose three rows and
    both blocks of rotated labels are measured to the millimetre - and its right end is
    taken, by the apo panel's n note, which hangs above that panel at axes y = 1.04.
    """
    import matplotlib.patches as mpatches
    PT = 25.4 / 72
    w = {t: style.text_w(t, fs) * PT for t in ("Re-predicted by", "AF3", "ESMFold2")}
    total = (w["Re-predicted by"] + 2.5 + w["AF3"] + 1.0 + 2 * KEY_SW + 1.0
             + w["ESMFold2"])
    assert total <= 2 * PW + GX, f"the half key is {total:.2f} mm, wider than two panels"
    assert KEY_SH <= GY, f"the half key is taller than the {GY} mm row gutter"
    yc = Y0 + PH + GY / 2            # centre of the gap between the grid's two rows
    x = X0
    for t in ("Re-predicted by", "AF3"):
        fig.text(x / W_MM, yc / H_MM, t, ha="left", va="center", color=palette.INK)
        x += w[t] + (2.5 if t == "Re-predicted by" else 1.0)
    for _suf, _name, _side, alpha in BACKENDS:
        fig.add_artist(mpatches.Rectangle(
            (x / W_MM, (yc - KEY_SH / 2) / H_MM), KEY_SW / W_MM, KEY_SH / H_MM,
            transform=fig.transFigure, facecolor=palette.INK, alpha=alpha, lw=0))
        x += KEY_SW
    fig.text((x + 1.0) / W_MM, yc / H_MM, "ESMFold2", ha="left", va="center",
             color=palette.INK)
    print(f"half key: {total:.2f} mm wide, {KEY_SH} mm tall, in the {GY} mm row gutter")


def check_limits(drawn):
    """Every DRAWN series must fit inside its panel's limits - each series on its own.

    `drawn` is (panel label, (lo, hi), series name, values), one entry per body or half
    actually drawn. A limit that lands inside a series' own 0.5-99.5 percentile range cuts
    that body off at the frame instead of letting it close, which no amount of TRIM can
    mend: TRIM only removes thin tails. POOLED percentiles hide exactly this, and the
    first render of the strip was set from them - 1.02% of RFdiffusion3's sc_DockQ and
    1.83% of HalluDesign's ligand_RMSD drew as flat-topped blocks against the frame.
    Since the split, the series are the HALVES, which reach further than the consensus
    they replaced in whichever direction the min/max rule was not looking: packstat's
    halves run above its min-consensus and sc_DockQ's above its, so both panels' upper
    limits come from this check and not from the consensus percentile.
    """
    bad = []
    print("\nper-series axis limits - every DRAWN series, never pooled "
          "(panel: limits | series 0.5-99.5 percentile)")
    for lab in dict.fromkeys(r[0] for r in drawn):
        rows = [r for r in drawn if r[0] == lab]
        (lo, hi) = rows[0][1]
        ps = [(n, *np.percentile(v[np.isfinite(v)], [0.5, 99.5])) for _l, _m, n, v in rows]
        need = (min(p[1] for p in ps), max(p[2] for p in ps))
        print(f"  {lab:18s} {lo:g}-{hi:g} | widest series needs {need[0]:.4g} to "
              f"{need[1]:.4g} ({len(ps)} series)")
        for n, p_lo, p_hi in ps:
            if not (lo <= p_lo and p_hi <= hi):
                bad.append(f"{lab} ({lo:g}, {hi:g}) cuts into {n}: {p_lo:.4g} to {p_hi:.4g}")
    assert not bad, "axis limits cut into a drawn series:\n  " + "\n  ".join(bad)


def main():
    style.apply()
    # Rotated y labels are bounded by the panel HEIGHT (STYLE_NOTES §3), and a label
    # width is measured, never estimated at 0.52 em/char - that stand-in is out by
    # 10-33% per string. Widest here is "ΔSASA major share" at 17.07 mm against 28.0.
    fs = plt.rcParams["font.size"]
    for _k, lab, *_r in METRICS + STRIP:
        w_mm = style.text_w(lab, fs) * 25.4 / 72
        assert w_mm <= PH, f"y label {lab!r} is {w_mm:.2f} mm, taller than the {PH} mm panel"
    # Each block of rotated x labels needs its own measured room, plus the tick and its
    # pad. The main grid's longest is "RFD3+LiMPNN", 15.05 mm all in at the house 5.5 pt,
    # and it now hangs into G_STRIP; the strip's longest is "BG+BIF corrected", 17.24 mm,
    # which is why YS is 18.5 and NOT the 16.5 that held the grid - the pair's labels are
    # longer than any route name. With style.save(tight=False) the authored canvas IS the
    # artefact, so a margin a fraction short clips the label instead of being absorbed by
    # a tight bounding box. Measured, never estimated (figure-measure-dont-estimate-text).
    pad_mm = (2 + plt.rcParams["xtick.major.pad"]) * 25.4 / 72
    need = max(style.text_w(r, fs) for r in ROUTES) * 25.4 / 72 + pad_mm
    need_s = max(style.text_w(r, fs) for r in STRIP_X) * 25.4 / 72 + pad_mm
    assert G_STRIP >= need, f"the gap above the strip ({G_STRIP} mm) cannot hold the " \
                            f"main grid's route labels ({need:.2f} mm)"
    assert YS >= need_s, f"bottom margin {YS} mm cannot hold the strip labels ({need_s:.2f} mm)"
    assert YS + 3 * PH + G_STRIP + GY + TOP <= H_MM, "the grid is taller than the canvas"
    assert X0 + 2 * PW_S + GX == X0 + 3 * PW + 2 * GX, "the strip and the grid disagree " \
                                                       "on the right margin"
    out = paths.fig_output(FIG)
    rfd3, strip = campaign_once()
    assert_consensus_rule(rfd3)
    d, apo = load()
    assert len(d) == N_CAMPAIGN, f"expected the whole campaign, got {len(d):,}"
    d = derive(d)

    # The limits, checked against every series that will be DRAWN - which since the split
    # are the halves, not the consensus. Runs before the first body is laid down.
    drawn = []
    for key, lab, _lg, lim, rule, _n in METRICS:
        frame = apo if rule == "one" else d
        for route in ROUTES:
            if rule == "one":
                drawn.append((lab, lim, route, series(frame, route, key)))
            else:
                drawn += [(lab, lim, f"{route} {nm}", series(frame, route, f"{key}_{suf}"))
                          for suf, nm, _s, _a in BACKENDS]
    for key, lab, _lg, lim in STRIP:
        drawn += [(lab, lim, f"{pos} {nm}", series(strip, pos, f"{key}_{suf}", col="pos"))
                  for pos in STRIP_X for suf, nm, _s, _a in BACKENDS]
    check_limits(drawn)

    fig = plt.figure(figsize=(style.mm(W_MM), style.mm(H_MM)))
    for i, (key, lab, log, (lo, hi), rule, n_note) in enumerate(METRICS):
        r, c = divmod(i, 3)
        frame = apo if rule == "one" else d
        ax = fig.add_axes([(X0 + c * (PW + GX)) / W_MM, (Y0 + (1 - r) * (PH + GY)) / H_MM,
                           PW / W_MM, PH / H_MM])
        for pos, route in enumerate(ROUTES):
            colour = palette.stack_color(route)
            if rule == "one":
                # One ESMFold2 measurement: one whole body. Splitting it would invent a
                # re-predictor, as a violin of "Designs" in Fig 2c would invent a spread.
                body(ax, float(pos), series(frame, route, key), lo, hi, log, colour)
            else:
                split(ax, float(pos),
                      lambda suf, r=route: series(frame, r, f"{key}_{suf}"),
                      lo, hi, log, colour)
        if log:
            ax.set_yscale("log")
        ax.set_ylim(lo, hi)
        ax.set_xlim(-0.7, len(ROUTES) - 0.3)
        ax.set_xticks(range(len(ROUTES)))
        # Only the bottom row is labelled. Every panel now has the same five positions
        # in the same order - the previous cut could not share a row because its two
        # split panels had six - so repeating the names six times is ink, not information.
        ax.set_xticklabels(ROUTES if r == 1 else [""] * len(ROUTES), rotation=90)
        ax.set_ylabel(lab, labelpad=2)
        ax.tick_params(axis="both", length=2, width=0.5)
        ax.minorticks_off()
        if n_note == "apo":
            # The one panel on a different population says so on the canvas: a datum,
            # not a title, and it names the contrast rather than its own n alone.
            # It also names the single re-predictor, because this is the one panel
            # whose body is NOT split and nothing else on the canvas would say why.
            note = (f"ESMFold2 only · n = {int(frame[key].notna().sum()):,} "
                    f"of {N_CAMPAIGN:,}")
            assert style.text_w(note, fs) * 25.4 / 72 <= PW, "n note wider than the panel"
            ax.annotate(note, (1.0, 1.04), xycoords="axes fraction",
                        ha="right", va="bottom", color=palette.INK)

    # The strip. Two panels, six positions: BoltzGen+BoltzIF twice, stored open and
    # corrected filled, in its own slot in the fixed route order. The other four design
    # routes are there as the scale the defect is read against - they are unaffected by
    # it, so they carry one body each and no stored/corrected distinction.
    for i, (key, lab, log, (lo, hi)) in enumerate(STRIP):
        ax = fig.add_axes([(X0 + i * (PW_S + GX)) / W_MM, YS / H_MM,
                           PW_S / W_MM, PH / H_MM])
        for pos, name in enumerate(STRIP_X):
            split(ax, float(pos),
                  lambda suf, n=name: series(strip, n, f"{key}_{suf}", col="pos"),
                  lo, hi, log,
                  palette.stack_color(BIF if name.startswith(BIF) else name),
                  filled=name != BIF_ST)
        if log:
            ax.set_yscale("log")
        ax.set_ylim(lo, hi)
        ax.set_xlim(-0.7, len(STRIP_X) - 0.3)
        ax.set_xticks(range(len(STRIP_X)))
        ax.set_xticklabels(STRIP_X, rotation=90)
        ax.set_ylabel(lab, labelpad=2)
        ax.tick_params(axis="both", length=2, width=0.5)
        ax.minorticks_off()
    half_key(fig, fs)
    # tight=False: every axes here is placed in mm with fig.add_axes, so the
    # authored 178 x 145 canvas IS the artefact. A tight bounding box would crop to
    # the ink and quietly resize what assembly places (style.save, 2026-11-19).
    style.save(fig, out, "metric_panel", tight=False)
    plt.close(fig)
    report(d, apo, strip)


def pearson(a, b):
    """Hand-rolled: scipy is unusable in this venv (_propack dlopen)."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m] - a[m].mean(), b[m] - b[m].mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def report(d, apo, strip):
    print(f"\n{len(d):,} designs, {len(ROUTES)} routes")
    for route in ROUTES:
        print(f"  {route:12s} {int((d['stack'] == route).sum()):>9,}")
    print(f"\napo: {len(apo):,} rows from both apo runs, "
          f"{int(apo['apo_rmsd_vs_esm'].notna().sum()):,} with a value "
          f"({100 * apo['apo_rmsd_vs_esm'].notna().sum() / N_CAMPAIGN:.1f}% of the "
          f"campaign). v2-v7 corpus only - it was never merged onto the master.")
    for route in ROUTES:
        n = int((apo["stack"] == route).sum())
        print(f"  {route:12s} {n:>9,}  "
              f"{100 * n / int((d['stack'] == route).sum()):>5.1f}% of that route")

    # The consensus, which is the rule the figure asserts and the README quotes, AND the
    # two halves that are now what is DRAWN. Both, because the panels no longer show the
    # consensus and a reader checking a median off the canvas needs the halves.
    for key, lab, _lg, _lim, rule, _n in METRICS:
        frame = apo if rule == "one" else d
        print(f"\n{lab}  rule={rule}  (median [IQR], n with a value / n in frame"
              + ("" if rule == "one" else "; then the drawn halves' medians") + ")")
        for route in ROUTES:
            v = pd.Series(series(frame, route, key)).dropna()
            n_route = int((frame["stack"] == route).sum())
            q = v.quantile([0.25, 0.5, 0.75]).to_numpy()
            halves = "" if rule == "one" else "   " + "  ".join(
                f"{nm} {np.nanmedian(series(frame, route, f'{key}_{suf}')):.4g}"
                for suf, nm, _s, _a in BACKENDS)
            print(f"  {route:12s} {q[1]:9.4g} [{q[0]:.4g}-{q[2]:.4g}]  "
                  f"{len(v):>9,}/{n_route:,}{halves}")

    print("\nwhy these five and not the rest of the family, re-derived every run: route "
          "separation as (max-min of route medians) / campaign IQR, and Pearson r "
          "against dsasa_total, which the main text already plots. A member with a high "
          "|r| is the interface-size panel restated.")
    for key, _lab, _lg, _lim, rule, _n in METRICS:
        if rule == "one":
            continue
        v = d[key].to_numpy(float)
        meds = np.array([np.nanmedian(series(d, r, key)) for r in ROUTES])
        q1, q3 = np.nanpercentile(v, [25, 75])
        print(f"  {key:24s} spread {(meds.max() - meds.min()) / (q3 - q1):5.2f} x IQR   "
              f"r vs dsasa_total {pearson(v, d['dsasa_total'].to_numpy(float)):+.3f}   "
              f"{int(np.isfinite(v).sum()):>9,}/{N_CAMPAIGN:,} with a value")

    # The BUNS choice, reported so the swap that was NOT made is on the record.
    for k, what in (("buns_vsasa", "whole-pose VSASA count, NOT an interface delta"),):
        v = d[k].to_numpy(float)
        meds = np.array([np.nanmedian(series(d, r, k)) for r in ROUTES])
        q1, q3 = np.nanpercentile(v, [25, 75])
        print(f"\nnot drawn, reported: {k} ({what})\n"
              f"  spread {(meds.max() - meds.min()) / (q3 - q1):5.2f} x IQR   "
              f"r vs dsasa_total {pearson(v, d['dsasa_total'].to_numpy(float)):+.3f}   "
              f"medians " + "  ".join(f"{m:.4g}" for m in meds))

    print("\nTHE STRIP - BoltzGen+BoltzIF stored vs corrected, the ONE pair that exists."
          "\n  The true-backbone rescore emits dockq_self and dockq_LRMSD per "
          "re-predictor and nothing else, so of the six metrics above NONE has a *_true "
          "counterpart in existence - they are read off the re-predicted complex, which "
          "the defect never touched. sc_DockQ and ligand_RMSD are in the main text (Fig 2c, "
          "Fig 3a, Fig 3b) and are re-admitted here only to carry this pair.")
    for key, lab, _lg, (lo, hi) in STRIP:
        print(f"\n{lab}  (median [IQR], n with a value / n at that position)")
        for name in STRIP_X:
            v = pd.Series(series(strip, name, key, col="pos")).dropna()
            n_pos = int((strip["pos"] == name).sum())
            q = v.quantile([0.25, 0.5, 0.75]).to_numpy()
            halves = "   " + "  ".join(
                f"{nm} {np.nanmedian(series(strip, name, f'{key}_{suf}', col='pos')):.4g}"
                for suf, nm, _s, _a in BACKENDS)
            print(f"  {name:17s} {q[1]:9.4g} [{q[0]:.4g}-{q[2]:.4g}]  "
                  f"{len(v):>9,}/{n_pos:,}{halves}")
        # What a POOLED percentile would have said, which is the trap check_limits()
        # exists for. The per-series numbers it works from are printed there, over the
        # halves, which are what is drawn.
        print(f"  limits drawn {lo}-{hi}; pooled would say "
              f"{np.percentile(strip[key].dropna().to_numpy(), 0.5):.4g} to "
              f"{np.percentile(strip[key].dropna().to_numpy(), 99.5):.4g}, which clips")

    # The common subset answers the obvious objection - that the shift is a coverage
    # artefact, stored and corrected being different populations. It is not: on the
    # designs carrying BOTH values the shift is the same size.
    pair = strip[strip["pos"].isin([BIF_ST, BIF_TR])]
    both = pair[pair["paired"]]
    n_st = int(pair.loc[pair["pos"] == BIF_ST, "dockq"].notna().sum())
    n_tr = int(pair.loc[pair["pos"] == BIF_TR, "dockq"].notna().sum())
    print(f"\n  stored {n_st:,} designs, corrected {n_tr:,}, "
          f"{int((both['pos'] == BIF_ST).sum()):,} carrying both, "
          f"of {N_BGBIF:,} in the design route")
    for key, lab, _lg, _lim in STRIP:
        s_all = pair.loc[pair["pos"] == BIF_ST, key].median()
        t_all = pair.loc[pair["pos"] == BIF_TR, key].median()
        s_c = both.loc[both["pos"] == BIF_ST, key].median()
        t_c = both.loc[both["pos"] == BIF_TR, key].median()
        print(f"  {lab:15s} median {s_all:.4f} -> {t_all:.4f}   "
              f"on the common subset {s_c:.4f} -> {t_c:.4f}")


if __name__ == "__main__":
    main()
