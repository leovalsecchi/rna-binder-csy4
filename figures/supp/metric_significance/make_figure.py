"""metric_significance - (a) metric effect sizes on the two screens, (b) the truncation evidence.

Panel (b), unchanged. Each of the five selection metrics, distributed in the GATED POOL
(55,706 designs, every survivor of the campaign gate) against the ORDERED PANEL (the
12,000 taken forward). Both are drawn as empirical cumulative distributions with the gate
threshold on the axis, because truncation is what the panel is about and an ECDF shows it
directly: every curve runs into its own threshold and stops.

Panel (a), REBUILT 2026-11-20 on Leo's populations: "the subset of the 12k that have >75
reads and dont have fc >10 and dont have SR larger than 0.07 vs the proposed candidates in
our one by one screen". That is the COMPLEMENT of the pooled screen's own hit gate
(code/Fig4c_fc_vs_splicing/make_figure.py, R_GATE/S_GATE/FC_GATE) as the negatives, and the
arrayed flow readout at ON/OFF >= 1.2 as the positives. The two figures read one gate: the
thresholds are PARSED out of Fig 4c's source below, so they cannot drift apart, and the
script asserts that the gate reproduces the 32 ordered hits exactly.

PRECISION-RECALL WAS DROPPED, and that is a result, not a style change. The positives fall
from 32 to 15, so the recall axis advances in steps of 1/15 = 6.7% and the 95% bootstrap
band on the one metric that separates is 41x wide at recall 0.1 and 27x at 0.3 - the curve
would be almost entirely band. The panel is now a rank-AUC effect size per metric with its
stratified bootstrap interval, which is what 15 positives can carry. The no-skill line is
drawn once, at AUC 0.5, with the prevalence stated as n/N on the panel.

TWO contrasts are drawn per metric, because the positives are a SUBSET of the 32 pooled
hits and the negatives are pooled non-hits, so the primary contrast stacks both screens and
cannot by itself say which one a separation comes from:
    A  15 proposed binders vs 8,020 measured non-hits   (pooled AND arrayed)
    B  15 proposed binders vs 17 inactive designs       (arrayed only, inside the 32)

STILL NO TEST IS RUN AND NO P-VALUE IS PRINTED ANYWHERE IN THIS FILE. An effect size with
an interval is not a test, and the AUC bounds below are arithmetic on the two group sizes.
The five gate terms are reported as NOT EVALUABLE, never as null results: panel (b) is the
evidence that each of them is truncated at its own threshold by construction, and that is
unaffected by the change of populations - both arms are still drawn from the gated pool.

scipy is unavailable in ~/.venvs/binder_map (its _propack extension fails to dlopen),
so the ECDFs, the ranks, the medians and the exceedance rates below are hand-rolled over
sorted numpy arrays. None of them needs more than that.

Run from the repo root:
    ~/.venvs/binder_map/bin/python3 code/supp/metric_significance/make_figure.py
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
N_POOL, N_PANEL = 55706, 12000
N_HITS = 32          # the pooled screen's hits, which is also the arrayed screen's panel

# The five ON terms of the campaign gate, read out of
# code/Fig3e_gpu_efficiency/methods/sources/gate/campaign_gate_55706.json.
# `sense` is the gate's own combination rule across the two gate backends:
# rmsd and sc_DockQ were OR (so the better backend decides), footnb was AND (so the
# worse one does). Taking the other end of either would not be the applied gate.
GATES = [
    ("ipTM", "iptm_max", None, ">=", 0.60, (0.55, 1.0)),
    ("ligand_RMSD (Å)", ["rmsd_af3", "rmsd_esm"], "min", "<=", 3.0, (0, 3.2)),
    ("sc_DockQ", ["dockq_af3", "dockq_esm"], "max", ">=", 0.70, (0.66, 1.0)),
    ("seq_run", "seq_run", None, "<=", 5, (0.5, 5.5)),
    ("footnb_010", ["footnb_010_af3", "footnb_010_esm"], "max", "<=", 5, (-0.4, 5.4)),
]

POOL_C = palette.NEUTRAL["500"]
PANEL_C = palette.MAIN

# --- the pooled screen's hit gate, PARSED from Fig 4c so one gate serves both figures ---
# Hard-coding these here would let the two figures drift silently, and the whole point of
# Leo's populations is that the negatives are Fig 4c's own complement. Reading them out of
# its source makes a change there a change here. The read-floor applies to BOTH arms of a
# comparison (reads and neg_reads): that NEG cut is what made the ordered set 32, not 34.
_FIG4C = paths.ROOT / "code" / "Fig4c_fc_vs_splicing" / "make_figure.py"
R_GATE, S_GATE, FC_GATE = (
    float(v) for v in re.search(
        r"^R_GATE, S_GATE, FC_GATE = ([\d.]+), ([\d.]+), ([\d.]+)",
        _FIG4C.read_text(), re.M).groups())

ONOFF_GATE = 1.2     # the standing arrayed threshold: >= is a proposed binder, < inactive

# --- panel (a)'s metric set, Leo's 2026-11-20 list --------------------------------------
# (label, column, direction, flag). Direction is fixed A PRIORI from what the metric means
# and NEVER from whichever way gives the larger AUC: reversing one maps AUC to 1 - AUC, so
# it is part of the result and is drawn on the panel as an arrow.
#
# Flags:
#   "gate"     - a campaign gate term. Scored, as Leo asked, but RANGE-RESTRICTED: every
#                design in both arms is a gate survivor (asserted in gate_truncation()), so
#                the effect is estimated inside the surviving range only.
#   "fixed"    - read from ONE named re-predictor, not from the selected model. ligand_RMSD
#                alone is treated this way, because ligand_RMSD is ALSO the variable the model
#                is chosen by: on the selected model its value is min(af3, esm) BY
#                CONSTRUCTION and the selection manufactures signal (0.602 selected against
#                0.458 fixed-AF3 and 0.553 fixed-ESMFold2 - higher than either backend
#                alone). Scoring it from a fixed backend makes it a measurement again, so
#                it stays a valid metric, as Leo asked, without the artefact. Both fixed
#                values get a row: they are close, they both sit at no effect, and that
#                agreement is itself the answer.
#   "uncertain"- footnb inverts in this project (it has acted as a filter in the opposite
#                sense to its meaning), so its a priori direction is NOT assumed. Drawn
#                without an arrow and scored as given, low-is-better, with the caveat.
#   "noback"   - a sequence metric: it has no re-predictor, so the both-backends rule and
#                the model choice do not apply to it.
#   "norm"     - divided by the chosen model's interface size, see metric_table().
GATE_SET = ["iptm", "sc_rmsd", "dockq", "seq_run", "footnb_010"]
METRICS = [
    ("ipTM", "iptm", +1, "gate"),
    ("ligand_RMSD (Å), AF3", "sc_rmsd_af3", -1, "gate fixed"),
    ("ligand_RMSD (Å), ESMFold2", "sc_rmsd_esm", -1, "gate fixed"),
    ("sc_DockQ", "dockq", +1, "gate"),
    ("seq_run", "seq_run", -1, "gate noback"),
    ("footnb_010", "footnb_010", -1, "gate uncertain"),
    ("packstat", "packstat", +1, ""),
    ("sc_value", "sc_value", +1, ""),
    ("Buried unsat. (optH)", "buns_vsasa_opth", -1, ""),
    ("H-bonds / iface res.", "n_h_bonds_per_ifaceres", +1, "norm"),
    ("π-stacking", "pi_stack_score", +1, ""),
    ("ΔSASA / iface res.", "dsasa_per_ifaceres", +1, "norm"),
    ("ΔSASA polar fraction", "dsasa_polar_frac", +1, ""),
    ("seq_max_aa_frac", "seq_max_aa_frac", -1, "noback"),
    ("Asn at interface", "if_N", +1, ""),
    ("Aromatic at interface", "if_arom", +1, ""),
    ("Acidic at interface", "if_neg", -1, ""),
]
N_BOOT = 1000
# Contrast A is the primary (and the one Leo asked for); B is drawn beside it because A
# alone cannot separate the two screens. A keeps the house single-series colour.
A_C, B_C = palette.MAIN, palette.THEME["400"]
BAND_C = palette.NEUTRAL["100"]

# --- panel (c): precision-recall, restored 2026-11-20 on Leo's "can you do the precision
# recall too?". FOUR curves, not sixteen: sixteen would be unreadable and most of them sit
# on the no-skill line anyway. The four are chosen by the FOREST's own verdicts, written
# down here before any curve was drawn, never by which curve looks best:
#   seq_max_aa_frac - the only metric that clears contrast A's resolvable bound
#   if_N            - the metric that leads contrast B, and this project's standing result
#   sc_value        - the highest-AUC purely structural non-gate metric in contrast A
#   iptm            - the gate term a reader most expects to work, and which does not
# Contrast A ONLY. PR earns its place where prevalence is 0.187% and a ROC would flatter
# everything; contrast B is 15 against 17, prevalence 46.9%, no imbalance to correct for,
# so a PR curve there is a different object with a no-skill line near 0.469 and would
# invite exactly the per-contrast misread the shaded bounds were fixed for.
PR_METRICS = ["seq_max_aa_frac", "if_N", "sc_value", "iptm"]
R_GRID = np.linspace(0.02, 1.0, 50)
PR_BAND_C = palette.THEME["200"]


def load():
    pool = pd.read_csv(paths.INPUTS / "large" / "gated_pool_v6_gate_metrics.csv.gz")
    lib = pd.read_csv(paths.INPUTS / "csy4_binder_library_12k_v6_k30_CORRECTED.csv",
                      usecols=["design_id"])
    assert len(pool) == N_POOL, f"gated pool has {len(pool)} rows, expected {N_POOL:,}"
    assert len(lib) == N_PANEL, f"ordered panel has {len(lib)} rows, expected {N_PANEL:,}"
    sel = pool.design_id.isin(set(lib.design_id))
    assert sel.sum() == N_PANEL, (
        f"only {sel.sum()} of the {N_PANEL:,} ordered designs are in the gated pool; "
        "the panel is supposed to be a subset of it")
    return pool, pool[sel]


def screen_sets():
    """The pooled screen, collapsed to the DESIGN, and the arrayed screen's two classes.

    screen_points.csv is one row per (design, comparison) - POS_A/NEG_A and POS_B/NEG_B are
    two backbone contexts carrying the same designs, not two replicates of a design - so a
    design-level rule is needed and is stated here rather than left to whichever row sorts
    first. The rule:

      * MEASURED: at least one comparison with reads >= 75 in BOTH its POS arm and its
        matched NEG control. Comparisons below the floor are not evidence either way and are
        dropped from the design entirely, not counted as a non-hit.
      * POOLED HIT: at least one measured comparison with splicing_ratio >= 0.07 AND
        fold_change >= 10 in the SAME comparison. This is Fig 4c's own per-point gate; a
        design that clears it in one context is a hit, which is how the 32 were picked.
      * NEGATIVE (Leo's words, taken literally): measured, and EVERY measured comparison has
        fold_change <= 10 AND splicing_ratio <= 0.07. This is stricter than "not a hit": it
        also drops the 843 designs that are high on one axis and low on the other. Those are
        the designs nearest the hit gate, so excluding them makes the negatives
        "comfortably non-hit" rather than merely "non-hit" - which flatters any separation,
        and is why the looser count is printed beside it.

    fold_change_raw is used, not matching_fold_change, because that is what Fig 4c gates on;
    matching_fold_change is the same number capped at 50 (65 of 21,514 rows differ) and gives
    the identical membership at a threshold of 10 - asserted below.
    """
    sp = pd.read_csv(paths.INPUTS / "screen_points.csv")
    sp["design_id"] = sp.transcript.str.split("::").str[-1]
    meas = sp[(sp.reads >= R_GATE) & (sp.neg_reads >= R_GATE)]
    g = meas.groupby("design_id")
    hit = g.apply(lambda d: bool(((d.splicing_ratio >= S_GATE)
                                  & (d.fold_change_raw >= FC_GATE)).any()),
                  include_groups=False)
    hit_capped = g.apply(lambda d: bool(((d.splicing_ratio >= S_GATE)
                                         & (d.matching_fold_change >= FC_GATE)).any()),
                         include_groups=False)
    assert hit.equals(hit_capped), (
        "fold_change_raw and matching_fold_change disagree about who is a pooled hit; "
        "the 50x cap in matching_fold_change is supposed to be irrelevant at FC >= 10")
    neg = g.apply(lambda d: bool(((d.fold_change_raw <= FC_GATE)
                                  & (d.splicing_ratio <= S_GATE)).all()),
                  include_groups=False)

    # The anti-drift check: Fig 4c's gate, applied here, has to give back exactly the 32
    # ordered hits. If Fig 4c's thresholds move, this is where it is noticed.
    ref = pd.read_csv(paths.INPUTS / "candidates32_mapped.csv", usecols=["design_id"])
    assert set(hit[hit].index) == set(ref.design_id), (
        "the pooled gate parsed from Fig 4c does not reproduce candidates32_mapped.csv")

    flow = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv")
    # The two Csy4-H29A controls are not designs: they carry no design_id and must not
    # become positives.
    ctrl = flow.design_id.isna()
    assert ctrl.sum() == 2, f"{ctrl.sum()} control rows in the flow file, expected 2"
    flow = flow[~ctrl]
    assert set(flow.design_id) == set(ref.design_id), (
        "the arrayed flow panel is not the 32 ordered hits")
    return dict(
        pos=set(flow.design_id[flow.ONOFF_avg >= ONOFF_GATE]),
        inactive=set(flow.design_id[flow.ONOFF_avg < ONOFF_GATE]),
        neg=set(neg[neg].index),
        not_hit=set(hit[~hit].index),        # the looser complement, for the README line
        hits=set(hit[hit].index),
        n_measured=int(hit.index.size),
    )


def metric_table():
    """One row per ordered design: a `model` column, then every metric read FROM THAT MODEL.

    Leo's rule, 2026-11-20: a metric is scored only if it exists from BOTH AF3 and
    ESMFold2, and per design every metric is read from **the one model with the lower
    ligand_RMSD**. That is a per-design model choice, not a per-metric best-of-two — the
    difference is the Fig 2b scar, where picking per-metric drew an AF3 pose 29 Å off its
    own call-out, and the fix was exactly this explicit `model` column. It is carried here
    for the same reason, and `report_a()` prints its split. Boltz-2 is dropped: never a
    gate backend, and it fails the both-backends rule for several of these metrics.

    🚨 **ligand_RMSD is the selection variable**, so its value in this table is by construction
    min(af3, esm) and an effect size on it is circular. It is flagged on the panel and the
    two fixed-backend numbers are printed beside it; those are the ones to quote.

    Size normalisation: every ordered binder is exactly 80 aa (asserted), so LENGTH is not
    a confound in this population. Interface SIZE still is — ΔSASA and an interface H-bond
    count both grow with it — so both are divided by the chosen model's own interface
    residue count `n_if`, which is the repo's existing idiom (`red_total_per_ifaceres_*`).
    The raw, unnormalised AUCs are printed beside them so the normalisation is visible.
    `dsasa_polar / dsasa_total` is Leo's ratio and is already size-free.
    """
    be = pd.read_csv(paths.INPUTS / "metric_significance_backends_12k.csv.gz")
    aa = pd.read_csv(paths.INPUTS / "if_interface_aa_12k.csv.gz")
    # that file is written backend-FIRST (af3_if_N); everything else here is
    # metric-first (if_N_af3), so flip it once rather than special-casing the lookup
    aa = aa.rename(columns={c: f"{c[4:]}_{c[:3]}" for c in aa.columns
                            if c.startswith(("af3_", "esm_"))})
    t = be.merge(aa, on="design_id", how="left")
    assert len(t) == N_PANEL, f"{len(t)} rows, expected {N_PANEL:,}"
    assert (t.length == 80).all(), (
        "not every ordered binder is 80 aa, so length IS a confound and the "
        "per-interface-residue normalisation below is no longer the whole story")

    model = np.where(t.sc_rmsd_af3 <= t.sc_rmsd_esm, "af3", "esm")
    out = pd.DataFrame({"design_id": t.design_id, "stack": t["stack"], "model": model})
    pick = lambda stem: np.where(model == "af3", t[f"{stem}_af3"], t[f"{stem}_esm"])

    for stem in ("iptm", "sc_rmsd", "dockq", "footnb_010", "packstat", "sc_value",
                 "buns_vsasa_opth", "pi_stack_score", "dsasa_total", "dsasa_polar",
                 "n_h_bonds_interface", "if_N", "if_arom", "if_neg", "n_if"):
        out[stem] = pick(stem)
    out["seq_run"], out["seq_max_aa_frac"] = t.seq_run, t.seq_max_aa_frac
    out["dsasa_per_ifaceres"] = out.dsasa_total / out.n_if
    out["n_h_bonds_per_ifaceres"] = out.n_h_bonds_interface / out.n_if
    out["dsasa_polar_frac"] = out.dsasa_polar / out.dsasa_total
    # The fixed-backend ligand_RMSD columns, kept only so the circular row can be reported
    # against something that is not circular.
    out["sc_rmsd_af3"], out["sc_rmsd_esm"] = t.sc_rmsd_af3, t.sc_rmsd_esm
    return out


def values(df, cols, how):
    if how is None:
        return df[cols].to_numpy(float)
    return (df[cols].min(axis=1) if how == "min" else df[cols].max(axis=1)).to_numpy(float)


def ecdf(v):
    """Step coordinates for the empirical CDF of v. Hand-rolled; scipy is broken here."""
    s = np.sort(v)
    y = np.arange(1, len(s) + 1) / len(s)
    return s, y


def panel_b(pool, pan, out):
    fig, axes = plt.subplots(1, len(GATES), figsize=(style.mm(178), style.mm(42)))
    for ax, (name, cols, how, cmp_, thr, xlim) in zip(axes, GATES):
        for df, c, lw in ((pool, POOL_C, 0.8), (pan, PANEL_C, 1.0)):
            x, y = ecdf(values(df, cols, how))
            ax.step(x, y, where="post", lw=lw, color=c)

        ax.axvline(thr, color=palette.NEUTRAL["600"], lw=0.5, dashes=(1.5, 1.5))
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 1.04)
        ax.set_xlabel(name)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(["0", "0.5", "1"] if ax is axes[0] else [])
        ax.tick_params(axis="both", length=2, width=0.5)
    axes[0].set_ylabel("Cumulative fraction")

    # The dashed line is named once, over the one sub-panel with clear space at the
    # top left; every other sub-panel uses the same mark for the same thing.
    axes[0].annotate("gate", (0.60, 1.0), textcoords="offset points", xytext=(2, -1),
                     ha="left", va="top", fontsize=style.CFG["font_size_pt"],
                     color=palette.NEUTRAL["600"])
    for i, (lab, c) in enumerate((("Gated pool", POOL_C), ("Ordered panel", PANEL_C))):
        x = 0.36 + i * 0.16
        fig.add_artist(plt.Line2D([x, x + 0.022], [0.035] * 2, lw=1.0, color=c,
                                  transform=fig.transFigure))
        fig.text(x + 0.028, 0.029, lab, fontsize=style.CFG["font_size_pt"],
                 va="baseline", ha="left")

    fig.subplots_adjust(left=0.055, right=0.995, top=0.97, bottom=0.30, wspace=0.18)
    style.save(fig, out, "panel_b_truncation_pool_vs_panel")
    plt.close(fig)


# --- (a) rank-AUC effect sizes ---------------------------------------------
def rank_auc(s, y):
    """P(a proposed binder outranks a comparison design). Mann-Whitney over AVERAGE ranks,
    which is what splits a tie 50/50; the discrete metrics here are mostly ties.
    Hand-rolled: scipy is unavailable in this venv."""
    r = pd.Series(s).rank().to_numpy()
    n1 = int(y.sum())
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1))


def pr_curve(s, y):
    """Precision and recall at every distinct score, and the average precision.

    Tied scores share one threshold: taking them one row at a time would draw a staircase
    inside a tie that no threshold can actually realise, which on a discrete metric is most
    of the curve. AP is the step sum sum((R_n - R_n-1) * P_n), not a trapezoid - trapezoids
    interpolate between thresholds that do not exist. Hand-rolled: no scipy here.
    """
    o = np.argsort(-s, kind="mergesort")
    ys, ss = y[o], s[o]
    tp = np.cumsum(ys)
    k = np.arange(1, len(ys) + 1)
    last = np.r_[ss[1:] != ss[:-1], True]        # last row of each tied block
    tp, k = tp[last], k[last]
    hit = tp > 0                                 # precision is undefined before the 1st
    P, R = (tp / k)[hit], (tp / y.sum())[hit]
    return R, P, float(np.sum(np.diff(np.r_[0, R]) * P))


def average_precision(s, y):
    return pr_curve(s, y)[2]


def resolvable(n1, n0):
    """The smallest rank AUC two groups this size could resolve, from the null SE of the
    Mann-Whitney statistic under a normal approximation. Arithmetic on the group sizes,
    NOT a test and NOT a result."""
    se = np.sqrt((n1 + n0 + 1) / (12 * n1 * n0))
    return se, 0.5 + (1.959964 + 0.841621) * se


def gate_truncation(pool, ids):
    """Are the five gate terms still truncated on BOTH arms of panel (a)?

    The brief for the rebuild assumed they are not — that Leo's negatives, being defined by
    the screen rather than by the gate, escape the truncation that made them unevaluable.
    **They do not.** Both arms are drawn from the 12,000 ordered designs, which are drawn
    from the 55,706 gate survivors, so every design on both sides satisfies all five terms
    by construction and each term's range is still the gated pool's truncated range. This
    function measures it and returns the rows for the README; the terms are therefore
    scored (as asked) but labelled RANGE-RESTRICTED, not reported as ordinary metrics.
    """
    sub = pool[pool.design_id.isin(ids)]
    rows = []
    for name, cols, how, cmp_, thr, _ in GATES:
        a, b = values(pool, cols, how), values(sub, cols, how)
        ok = (b >= thr) if cmp_ == ">=" else (b <= thr)
        rows.append(dict(term=name, gate=f"{cmp_} {thr}", n_pass=int(ok.sum()), n=len(b),
                         pool_lo=a.min(), pool_hi=a.max(), lo=b.min(), hi=b.max(),
                         kept=float(np.ptp(b) / np.ptp(a))))
    return pd.DataFrame(rows)


def contrast(pool, pos, neg, seed=0, pr_metrics=()):
    """AUC and AP per metric with a 95% stratified bootstrap interval.

    Designs are resampled stratified: the positives and the negatives separately, so every
    replicate carries the same n1 and the same prevalence. ONE resample per replicate,
    SHARED by every metric, because the design is the independent unit and all the metrics
    are measured on the same designs; resampling each metric separately would break that.

    What this does NOT cover: there is one biological replicate of each screen, so the
    interval is sampling of designs only, never experimental reproducibility.
    """
    sub = pool[pool.design_id.isin(pos | neg)].reset_index(drop=True)
    assert len(sub) == len(pos | neg), (
        f"{len(pos | neg) - len(sub)} of the contrast's designs are not in the metric table")
    y = sub.design_id.isin(pos).to_numpy()
    scores = {m: sg * sub[m].to_numpy(float) for _, m, sg, _ in METRICS}
    for m, s in scores.items():
        assert np.isfinite(s).all(), f"{m} has {int((~np.isfinite(s)).sum())} missing values"

    rng = np.random.default_rng(seed)
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
    yb = np.r_[np.ones(len(pi), bool), np.zeros(len(ni), bool)]
    boot = {m: np.empty(N_BOOT) for m in scores}
    # The PR band rides the SAME resample as the AUC intervals - one resample per
    # replicate for everything in this contrast, because the design is the independent
    # unit and every metric is measured on the same designs.
    pband = {m: np.empty((N_BOOT, len(R_GRID))) for m in pr_metrics}
    apboot = {m: np.empty(N_BOOT) for m in pr_metrics}
    for b in range(N_BOOT):
        idx = np.r_[rng.choice(pi, len(pi), replace=True),
                    rng.choice(ni, len(ni), replace=True)]
        for m, s in scores.items():
            boot[m][b] = rank_auc(s[idx], yb)
        for m in pr_metrics:
            R, P, ap = pr_curve(scores[m][idx], yb)
            pband[m][b] = P[np.searchsorted(R, R_GRID, side="left").clip(max=len(R) - 1)]
            apboot[m][b] = ap

    se, need = resolvable(int(y.sum()), int((~y).sum()))
    return dict(n1=int(y.sum()), n0=int((~y).sum()), prev=float(y.mean()),
                se=se, need=need, y=y, sub=sub, scores=scores,
                pband={m: np.percentile(pband[m], [2.5, 97.5], axis=0)
                       for m in pr_metrics},
                ap_ci={m: np.percentile(apboot[m], [2.5, 97.5]) for m in pr_metrics},
                stats={m: dict(auc=rank_auc(s, y), ap=average_precision(s, y),
                               auc_ci=np.percentile(boot[m], [2.5, 97.5]))
                       for m, s in scores.items()})


def panel_a(A, B, out):
    """Sixteen rows, two contrasts, rank AUC with its 95% bootstrap interval.

    No PR curve: see the module docstring. The no-skill line is AUC = 0.5, drawn once,
    with the prevalence carried as n/N in the key. The five gate terms are the top block,
    separated by a rule and labelled range-restricted; the circular one (ligand_RMSD) is drawn
    with a hollow marker so it cannot be read as an ordinary result.
    """
    n = len(METRICS)
    fig, ax = plt.subplots(figsize=(style.mm(178), style.mm(17.5 + 4.6 * n)))

    # The shading is PER CONTRAST, in half-row strips, so every marker sits on its own
    # bound. One band for both would be a misread waiting to happen: contrast A resolves
    # at 0.709 and contrast B at 0.791, so if_N's light marker at 0.784 falls OUTSIDE a
    # band drawn at 0.709 while being INSIDE its own contrast's bound - a reader takes it
    # as the one metric that cleared, which is the opposite of the result. The step in the
    # shaded edge falls between the two markers of a row and is the cue that they are
    # judged against different bounds. Group sizes only - not a test - and each band is
    # symmetric because a reversed direction maps AUC to 1 - AUC.
    for i in range(n):
        for d, y0 in ((A, i - 0.5), (B, i)):
            ax.add_patch(plt.Rectangle((1 - d["need"], y0), 2 * d["need"] - 1, 0.5,
                                       facecolor=BAND_C, lw=0, zorder=0))
    ax.axvline(0.5, color=palette.NEUTRAL["600"], lw=0.5, dashes=(1.5, 1.5), zorder=2)

    for i, (label, m, sg, flag) in enumerate(METRICS):
        for d, c, dy in ((A, A_C, -0.17), (B, B_C, +0.17)):
            s = d["stats"][m]
            lo, hi = s["auc_ci"]
            ax.plot([lo, hi], [i + dy] * 2, lw=0.7, color=c, zorder=3,
                    solid_capstyle="butt")
            for x in (lo, hi):          # thin caps, as section 7 asks for
                ax.plot([x] * 2, [i + dy - 0.07, i + dy + 0.07], lw=0.7, color=c, zorder=3)
            ax.plot([s["auc"]], [i + dy], "o", ms=3.0, mfc=c, mec="white", mew=0.4,
                    zorder=4)

    # The rule between the gate block and the rest, drawn in the data's own coordinates so
    # it stays put whatever the canvas height is.
    ng = sum("gate" in f for _, _, _, f in METRICS)
    ax.axhline(ng - 0.5, color=palette.NEUTRAL["300"], lw=0.4, zorder=2)

    # "?" rather than an arrow where the direction is NOT assumed (footnb inverts here).
    # The circular row carries no glyph: Arial has no warning sign and it rendered as
    # tofu; the hollow marker is what flags it, and the note names the marker.
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{lab}  {'?' if 'uncertain' in f else ('↑' if sg > 0 else '↓')}"
                        for lab, _, sg, f in METRICS])
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0.05, 0.95)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.set_xlabel("Rank AUC")
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)

    # Two header lines only. The hairline rule already separates the gate block from the
    # rest, so a second "not selection criteria" label would be the caption's job twice.
    top = 1 - style.mm(2.6) / fig.get_figheight()
    fig.text(0.008, top, "Selection criteria — range-restricted: every design on both "
             "sides is a gate survivor", va="top",
             fontsize=style.CFG["font_size_pt"])
    fig.text(0.992, top, "shaded: below the resolvable effect, EACH CONTRAST AGAINST ITS "
             f"OWN (AUC {A['need']:.3f} dark, {B['need']:.3f} light)", va="top",
             ha="right", fontsize=style.CFG["font_size_pt"],
             color=palette.NEUTRAL["600"])
    fig.text(0.992, top - style.mm(2.3) / fig.get_figheight(),
             "every metric from the re-predictor with the lower ligand_RMSD, except ligand_RMSD "
             "itself — fixed, since it is what the choice is made on",
             va="top", ha="right", fontsize=style.CFG["font_size_pt"],
             color=palette.NEUTRAL["600"])

    # Frameless key along the bottom, shared by the whole panel.
    keys = [(A_C, f"vs measured non-hits, {A['n1']}/{A['n1'] + A['n0']:,} "
                  f"({100 * A['prev']:.3f}%)"),
            # "15/32" beside A's "15/8,035" would invite reading 32 as the same kind of
            # denominator, so B's key names the 17 and says what the 32 is.
            (B_C, f"vs {B['n0']} inactive designs, {B['n1']}/{B['n1'] + B['n0']} of the "
                  f"arrayed panel ({100 * B['prev']:.1f}%)")]
    # Measured, not estimated: a 0.52 em/char guess is out by 10-33% per string and the
    # gaps would be visibly uneven. text_w returns points at the point size given.
    pt_per_fig = 178 / 25.4 * 72
    ky = style.mm(2.2) / fig.get_figheight()
    x = 0.165
    for c, lab in keys:
        fig.add_artist(plt.Line2D([x, x + 0.018], [ky] * 2, lw=0.7, color=c,
                                  transform=fig.transFigure))
        fig.add_artist(plt.Line2D([x + 0.009], [ky], marker="o", ms=3.0, mfc=c,
                                  mec="white", mew=0.4, transform=fig.transFigure))
        fig.text(x + 0.024, ky - style.mm(0.6) / fig.get_figheight(), lab,
                 fontsize=style.CFG["font_size_pt"], va="baseline", ha="left")
        x += 0.024 + style.text_w(lab) / pt_per_fig + 0.035
    fig.add_artist(plt.Line2D([x, x + 0.018], [ky] * 2, lw=0.5, dashes=(1.5, 1.5),
                              color=palette.NEUTRAL["600"], transform=fig.transFigure))
    fig.text(x + 0.024, ky - style.mm(0.6) / fig.get_figheight(), "no skill",
             fontsize=style.CFG["font_size_pt"], va="baseline", ha="left",
             color=palette.NEUTRAL["600"])

    h = fig.get_figheight()
    fig.subplots_adjust(left=0.135, right=0.992, top=1 - style.mm(7.4) / h,
                        bottom=style.mm(10.5) / h)
    style.save(fig, out, "panel_a_auc_arrayed_vs_nonhits")
    plt.close(fig)


def panel_c(A, out):
    """Precision against recall for four metrics, contrast A only.

    Restored 2026-11-20 on Leo's "can you do the precision recall too?" - alongside the
    forest, not instead of it. PR is drawn ONLY for contrast A: prevalence 0.187% is the
    imbalance a PR curve exists to show, while contrast B is 15 against 17 (46.9%) and has
    no imbalance to correct for.

    THE BAND IS DRAWN AT ITS FULL WIDTH AND IS THE POINT OF THE PANEL. At 15 positives
    recall advances in steps of 1/15 = 6.7%, so the left-hand end of every curve rests on
    one or two designs; the 95% band there is tens of times wide. Nothing is clipped or
    smoothed to make it look better - the panel asserts that no curve here can be read, and
    the forest is where the argument lives. The script asserts no band is clipped.
    """
    fig, axes = plt.subplots(1, len(PR_METRICS), figsize=(style.mm(178), style.mm(48)),
                             sharex=True, sharey=True)
    prev = A["prev"]
    # Just below the lowest point any band reaches, NOT the decade below it: flooring to
    # a decade added an empty decade under the data and squashed everything that matters.
    lo_all = min(A["pband"][m][0].min() for m in PR_METRICS)
    ylo = min(lo_all, prev) / 1.6

    for ax, m in zip(axes, PR_METRICS):
        lab, _, sg, _ = next(r for r in METRICS if r[1] == m)
        lo, hi = A["pband"][m]
        ax.fill_between(R_GRID, lo, hi, color=PR_BAND_C, lw=0, zorder=1)
        R, P, ap = pr_curve(A["scores"][m], A["y"])
        ax.step(R, P, where="post", lw=0.9, color=A_C, zorder=3)
        # The no-skill line, in EVERY cell: a PR curve cannot be read without it. One
        # prevalence only, because only contrast A is drawn here.
        ax.axhline(prev, color=palette.NEUTRAL["600"], lw=0.5, dashes=(1.5, 1.5), zorder=2)
        alo, ahi = A["ap_ci"][m]
        ax.annotate(f"{lab} {'↑' if sg > 0 else '↓'}\nAP {100 * ap:.2f}%\n"
                    f"95% {100 * alo:.2f}–{100 * ahi:.2f}%",
                    (0.985, 0.98), xycoords="axes fraction", ha="right", va="top",
                    fontsize=style.CFG["font_size_pt"], linespacing=1.45)
        ax.set_yscale("log")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(ylo, 1.6)
        ax.set_xticks([0, 0.5, 1])
        ax.tick_params(axis="both", length=2, width=0.5)
        assert lo.min() >= ylo, f"{m}: the band is clipped at the bottom of the axis"

    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    dec = [10.0 ** e for e in range(-4, 1) if ylo <= 10.0 ** e <= 1.6]
    axes[0].set_yticks(dec)
    axes[0].set_yticklabels([f"{100 * d:g}%" for d in dec])
    for ax in axes[1:]:
        ax.tick_params(labelleft=False, labelbottom=False)

    h = fig.get_figheight()
    fig.text(0.008, 1 - style.mm(2.6) / h,
             f"15 proposed binders ranked against {A['n0']:,} measured non-hits; band, "
             f"95% over {N_BOOT:,} stratified bootstrap resamples of designs",
             va="top", fontsize=style.CFG["font_size_pt"])
    fig.text(0.992, 1 - style.mm(2.6) / h,
             f"no skill {A['n1']}/{A['n1'] + A['n0']:,} ({100 * prev:.3f}%)", va="top",
             ha="right", fontsize=style.CFG["font_size_pt"],
             color=palette.NEUTRAL["600"])
    fig.subplots_adjust(left=0.052, right=0.992, top=1 - style.mm(6.0) / h,
                        bottom=style.mm(7.0) / h, wspace=0.10)
    style.save(fig, out, "panel_c_precision_recall")
    plt.close(fig)

    print(f"\npanel (c) precision-recall, contrast A only ({A['n1']}/"
          f"{A['n1'] + A['n0']:,} = {prev:.3%}), four metrics chosen from the forest's own "
          f"verdicts before any curve was drawn:")
    for m in PR_METRICS:
        lo, hi = A["pband"][m]
        i10, i30 = int(np.searchsorted(R_GRID, 0.1)), int(np.searchsorted(R_GRID, 0.3))
        w = lambda i: hi[i] / max(lo[i], 1e-12)
        print(f"  {m:<18} AP {A['stats'][m]['ap']:.5f} "
              f"({A['stats'][m]['ap'] / prev:.2f}x no skill)  band width "
              f"{w(i10):.0f}x at recall 0.1, {w(i30):.0f}x at 0.3")
    print("  recall advances in steps of 1/15 = 6.7%; the band IS the result, not a defect")


def report_metrics(A, B, tab, trunc):
    """Everything about the metric set that the panel cannot carry."""
    mdl = A["sub"].model.value_counts()
    print(f"\nmodel choice (lower ligand_RMSD, per design) over the {len(A['sub']):,} designs "
          f"of contrast A: af3 {mdl.get('af3', 0):,} ({mdl.get('af3', 0)/len(A['sub']):.1%}), "
          f"esm {mdl.get('esm', 0):,} ({mdl.get('esm', 0)/len(A['sub']):.1%})")
    print("  carried as an explicit `model` column, not implied: picking per METRIC "
          "instead of per DESIGN is the Fig 2b scar (an AF3 pose 29 A off its own call-out)")

    y = A["y"]
    sel = rank_auc(-A["sub"].sc_rmsd.to_numpy(float), y)
    print("\nligand_RMSD stays a scored metric, but FIXED-BACKEND, because it is also the "
          "variable the model is chosen by. The record of what the choice does:")
    for be in ("af3", "esm"):
        print(f"   fixed {be:<4} (drawn): AUC "
              f"{A['stats'][f'sc_rmsd_{be}']['auc']:.3f}")
    print(f"   selected model (NOT drawn, kept as the demonstration): AUC {sel:.3f} — "
          "higher than either backend alone, which is the selection manufacturing signal, "
          "not the metric carrying it")
    print("   the project has hit this before: 'ligand_RMSD is NOT a predictor — I tested the "
          "selected min = circular', so fixed-backend is the established handling")

    print("\n⚠️ the five gate terms are STILL range-restricted. The rebuild brief assumed "
          "Leo's negatives escape the gate; they do not — both arms are drawn from the "
          "12,000 ordered, which are drawn from the 55,706 gate survivors:")
    print(f"   {'term':<14}{'gate':>10}{'pass BOTH arms':>22}{'gated-pool range':>22}"
          f"{'contrast range':>22}{'range kept':>12}")
    for _, r in trunc.iterrows():
        print(f"   {r.term:<14}{r.gate:>10}{f'{r.n_pass:,}/{r.n:,} ({r.n_pass/r.n:.1%})':>22}"
              f"{f'{r.pool_lo:.3f}-{r.pool_hi:.3f}':>22}"
              f"{f'{r.lo:.3f}-{r.hi:.3f}':>22}{r.kept:>11.1%}")
    print("   so they are SCORED, as asked, but labelled range-restricted on the panel: "
          "the effect is estimated inside the surviving range only. Panel (b) is the "
          "picture of that restriction and stands unchanged.")

    print("\nsize normalisation: every ordered binder is 80 aa (asserted), so length is "
          "not a confound here; interface SIZE is, and the two size-coupled metrics are "
          "divided by the chosen model's own interface residue count:")
    for raw, norm in (("dsasa_total", "dsasa_per_ifaceres"),
                      ("n_h_bonds_interface", "n_h_bonds_per_ifaceres")):
        a = rank_auc(A["sub"][raw].to_numpy(float), y)      # both are scored +1
        print(f"   {raw:<22} raw AUC {a:.3f}  ->  {norm} AUC "
              f"{A['stats'][norm]['auc']:.3f}  (rho with n_if "
              f"{np.corrcoef(A['sub'][raw].rank(), A['sub'].n_if.rank())[0, 1]:+.3f})")
    print("   dsasa_polar_frac = dsasa_polar / dsasa_total is Leo's ratio and is already "
          "size-free")
    print("\n⚠️ `dockq_self` and the gate's `dockq` are the SAME COLUMN (max |d| = 0.0000 "
          "over all 12,000, asserted in export_data.py), so Leo's list has one duplicate: "
          f"{len(METRICS)} rows scoring {len(METRICS) - 1} distinct metrics, ligand_RMSD "
          "counted once and measured on each backend.")
    print("⚠️ footnb_010 is drawn WITHOUT an arrow: footnb has inverted as a filter in "
          "this project, so its a priori direction is not assumed. It is scored "
          "low-is-binder-like and the sign of its result must not be read as confirming "
          "that direction.")


def report_a(A, B, sets, pool, lib_ids):
    print(f"\npanel (a) populations, from Fig 4c's own gate "
          f"(reads >= {R_GATE:.0f} in both arms, splicing >= {S_GATE}, FC >= {FC_GATE:.0f}) "
          f"and the arrayed cut ON/OFF >= {ONOFF_GATE}")
    print(f"  measured designs (read floor, either context): {sets['n_measured']:,}"
          f" of the {N_PANEL:,} ordered ({sets['n_measured'] / N_PANEL:.1%})")
    print(f"  pooled hits: {len(sets['hits'])}  "
          f"= candidates32_mapped.csv exactly (asserted)")
    print(f"  NEGATIVES, literal (FC <= {FC_GATE:.0f} AND SR <= {S_GATE} in every measured "
          f"comparison): {len(sets['neg']):,}")
    print(f"  the looser complement (simply not a pooled hit): {len(sets['not_hit']):,}; "
          f"the {len(sets['not_hit']) - len(sets['neg'])} designs between the two are high "
          f"on one axis and low on the other")
    print(f"  POSITIVES, proposed binders (ON/OFF >= {ONOFF_GATE}): {len(sets['pos'])}"
          f" of {len(sets['pos']) + len(sets['inactive'])}; "
          f"inactive designs: {len(sets['inactive'])}")
    print(f"  OVERLAP: {len(sets['pos'] & sets['neg'])} of the {len(sets['pos'])} proposed "
          f"binders are in the negative set, and {len(sets['pos'] & sets['hits'])}/"
          f"{len(sets['pos'])} are pooled hits - the hit gate excludes every positive from "
          f"the negatives BY CONSTRUCTION, so contrast A stacks both screens")
    print(f"  both arms are inside the ordered panel: "
          f"{len((sets['pos'] | sets['neg']) & lib_ids):,}/"
          f"{len(sets['pos'] | sets['neg']):,}")

    for tag, d, what in (("A", A, "15 proposed binders vs measured non-hits"),
                         ("B", B, "15 proposed binders vs inactive designs, inside the 32")):
        print(f"\ncontrast {tag}: {what}")
        print(f"  n1 = {d['n1']}, n0 = {d['n0']:,}, prevalence "
              f"{d['n1']}/{d['n1'] + d['n0']:,} ({d['prev']:.3%}); null SE of a rank AUC "
              f"{d['se']:.4f}, so two-sided 0.05 at 80% power needs AUC >= {d['need']:.3f} "
              f"(arithmetic on the group sizes, NOT a test)")
        print(f"  {'metric':<24}{'dir':>4}{'AUC':>8}{'AUC 95% boot':>18}"
              f"{'AP':>9}{'AP/no-skill':>13}  resolvable  flag")
        for label, m, sg, flag in METRICS:
            s = d["stats"][m]
            lo, hi = s["auc_ci"]
            res = "yes" if (s["auc"] >= d["need"] or s["auc"] <= 1 - d["need"]) else "-"
            print(f"  {m:<24}{'+' if sg > 0 else '-':>4}{s['auc']:>8.3f}"
                  f"{f'{lo:.3f}-{hi:.3f}':>18}{s['ap']:>9.5f}"
                  f"{s['ap'] / d['prev']:>12.2f}x  {res:>10}  {flag}")

    # Leo's second question: how were the seven chosen? The rule is mechanical - they are
    # whatever was not a gate term - so the honest check is whether they are independent of
    # the gate, on the NEW negative population rather than on the old note's.
    print("\nindependence from the gate, recomputed on the new negatives "
          f"(n = {len(sets['neg']):,}), beside the old 12,000 panel. Spearman |rho| on ranks:")
    # Sensitivity to the negative rule. The literal reading of Leo's sentence gives 8,020
    # and the looser "simply not a pooled hit" gives 8,863; the 843 between are high on one
    # axis and low on the other. Point AUCs only - no second bootstrap, because the
    # question is whether any conclusion moves, not how wide the interval is.
    print(f"\nsensitivity to the negative rule, {len(sets['neg']):,} (literal, built) "
          f"against {len(sets['not_hit']):,} (looser):")
    out = []
    for ids in (sets["neg"], sets["not_hit"]):
        sub = pool[pool.design_id.isin(sets["pos"] | ids)].reset_index(drop=True)
        yy = sub.design_id.isin(sets["pos"]).to_numpy()
        _, need = resolvable(int(yy.sum()), int((~yy).sum()))
        out.append((need, {m: rank_auc(sg * sub[m].to_numpy(float), yy)
                           for _, m, sg, _ in METRICS}))
    (nA, aA), (nB, aB) = out
    worst = max(METRICS, key=lambda r: abs(aB[r[1]] - aA[r[1]]))
    clears = [sorted(m for _, m, _, _ in METRICS if a[m] >= n or a[m] <= 1 - n)
              for n, a in out]
    print(f"  max |dAUC| = {abs(aB[worst[1]] - aA[worst[1]]):.3f} ({worst[1]}); "
          f"resolvable bound {nA:.3f} -> {nB:.3f}; clears its bound: "
          f"{clears[0]} -> {clears[1]}")
    assert clears[0] == clears[1], "the negative rule changes which metrics clear"
    print("  nothing moves: the rule is a presentation choice, not a result")

    hdr = "".join(f"{g:>12}" for g in GATE_SET)
    pops = (("old: 12,000 ordered panel", lib_ids),
            ("new: 8,020 measured non-hits", sets["neg"]))
    for name, ids in pops:
        sub = pool[pool.design_id.isin(ids)]
        print(f"  -- {name} (n = {len(sub):,})")
        print(f"     {'metric':<24}{hdr}{'max':>8}")
        for _, m, _, flag in METRICS:
            if "gate" in flag:
                continue            # a gate term's correlation with itself is 1
            rs = [abs(np.corrcoef(sub[m].rank(), sub[g].rank())[0, 1]) for g in GATE_SET]
            print(f"     {m:<24}" + "".join(f"{r:>12.3f}" for r in rs) + f"{max(rs):>8.3f}")


def report(pool, pan):
    print(f"gated pool n = {len(pool):,}   ordered panel n = {len(pan):,}")
    print(f"{'metric':<12} {'gate':>12}  {'pool median':>11} {'panel median':>12}"
          f"  {'pool at gate':>12} {'panel at gate':>13}")
    for name, cols, how, cmp_, thr, _ in GATES:
        a, b = values(pool, cols, how), values(pan, cols, how)
        ok = (lambda v: v >= thr) if cmp_ == ">=" else (lambda v: v <= thr)
        assert ok(a).all() and ok(b).all(), (
            f"{name}: a design violates its own gate term - this is not the applied gate")
        # "at the gate": within 5% of the term's own passing range of the threshold.
        span = 0.05 * (np.ptp(a) if np.ptp(a) > 0 else 1)
        near = (lambda v: v <= thr + span) if cmp_ == ">=" else (lambda v: v >= thr - span)
        print(f"{name:<12} {cmp_ + ' ' + str(thr):>12}  {np.median(a):>11.3f} "
              f"{np.median(b):>12.3f}  {near(a).mean():>11.1%} {near(b).mean():>12.1%}")

    ip_a, ip_b = pool.iptm_max.to_numpy(), pan.iptm_max.to_numpy()
    print(f"\nipTM is the only term that shifts: median {np.median(ip_a):.3f} -> "
          f"{np.median(ip_b):.3f}; >= 0.8 in {(ip_a >= 0.8).sum():,}/{len(ip_a):,} "
          f"({(ip_a >= 0.8).mean():.1%}) of the pool against "
          f"{(ip_b >= 0.8).sum():,}/{len(ip_b):,} ({(ip_b >= 0.8).mean():.1%}) of the panel")


def main():
    style.apply()
    out = paths.fig_output(FIG)
    pool, pan = load()
    lib_ids = set(pan.design_id)
    sets = screen_sets()
    tab = metric_table()
    A = contrast(tab, sets["pos"], sets["neg"], pr_metrics=PR_METRICS)
    B = contrast(tab, sets["pos"], sets["inactive"], seed=1)
    trunc = gate_truncation(pool, sets["pos"] | sets["neg"] | sets["inactive"])
    panel_b(pool, pan, out)
    panel_a(A, B, out)
    panel_c(A, out)
    report(pool, pan)
    report_a(A, B, sets, tab, lib_ids)
    report_metrics(A, B, tab, trunc)


if __name__ == "__main__":
    main()
