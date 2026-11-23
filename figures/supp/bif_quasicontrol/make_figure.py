"""bif_quasicontrol - the BoltzIF scoring defect used as a quasi-control on the gate.

Two of the five campaign gates, ligand_RMSD and sc_DockQ, were effectively NOT APPLIED to
BoltzGen + BoltzIF: both were scored against BoltzIF's own Boltz re-prediction rather
than against the generated backbone (mechanism and magnitude: code/supp/bif_effect).
One fifth of the ordered panel therefore entered under a materially weaker filter, so
the panel carries an internal comparison the rest of this work has no baseline for.

QUASI-control, never a control. Those designs still passed the other three gates and
were still ipTM-ranked within cluster - a weakened filter, not no filter - and they
were selected for scoring well under the wrong reference, so regression is expected.

Encoding: panel (a) is a route comparison, so route colour from palette.stack_color()
in palette.STACK_ORDER (COLORS.md rule 2), never sorted by value. Panel (b) pools four
routes into one group, which is not a route identity and so takes a neutral (rule 3);
its two BoltzGen + BoltzIF bars keep the route colour and read their gate status off
the tick labels, no second encoding needed (STYLE_NOTES section 0).

scipy is unusable in this venv (_propack dlopen failure), so every statistic below is
hand-rolled: Pearson chi-square with a closed-form upper tail, Fisher's exact test from
math.comb, and a label-permutation p.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from math import comb, erfc, exp, factorial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
STORED = paths.INPUTS / "csy4_binder_library_12k_v6_k30.csv"
CORRECTED = paths.INPUTS / "csy4_binder_library_12k_v6_k30_CORRECTED.csv"
CANDIDATES = paths.INPUTS / "candidates32_mapped.csv"

ROUTE = "BG+BIF"                  # the literal `stack` value in the CSV
RMSD_GATE, DOCKQ_GATE = 3.0, 0.7  # identical to code/supp/bif_effect, do not drift
N_PERM, SEED = 20000, 42

# Short route names for the axis; palette.SHORT is the same mapping the other
# figures use. Order is palette.STACK_ORDER, never the result's order.
ORDER = [palette.SHORT[k] for k in palette.STACK_ORDER]


# --- data -------------------------------------------------------------------
def load():
    """The 12,000-design panel: candidate flag, route, and the CORRECTED gate.

    `cand` marks the 32 designs selected for ordering on ranked splicing ratio and
    fold change. A selection, not a hit rate.

    The CORRECTED library is the repo's own join of the original against
    `large/bif_true_per_backend.csv.gz` (inputs/MANIFEST.md: use it instead of the
    original), so the *_true columns are already in place for BoltzGen + BoltzIF and
    nothing else is touched. Verified element-for-element against the join
    code/supp/bif_effect performs, which is what keeps the two figures from drifting.
    """
    sel = pd.read_csv(CANDIDATES)
    # Provenance check on the right 32: another candidate list would not split 7/25.
    assert len(sel) == 32 and sel.hts_comparison.value_counts().to_dict() == \
        {"POS_B": 25, "POS_A": 7}, sel.hts_comparison.value_counts().to_dict()

    d = pd.read_csv(CORRECTED)
    d["cand"] = d.design_id.isin(set(sel.design_id)).astype(int)
    assert len(d) == 12000 and d.cand.sum() == 32, (len(d), d.cand.sum())

    # The panel is the GATED set: every member passes on the scores the campaign saw.
    # That is what makes this a weakened filter rather than no filter, so it is
    # asserted off the original, pre-correction library rather than assumed.
    st = pd.read_csv(STORED, usecols=["design_id", "sc_rmsd_af3", "sc_rmsd_esm",
                                      "dockq_af3", "dockq_esm"])
    st = d[["design_id"]].merge(st, on="design_id", validate="one_to_one")
    assert gate(st.sc_rmsd_af3, st.sc_rmsd_esm, st.dockq_af3,
                st.dockq_esm).all(), "the panel is not the gated set"

    # Three BoltzGen + BoltzIF members have no true-backbone score (bif_effect: one
    # genuinely unmeasurable, two no_native_interface). None is a candidate, so
    # dropping them cannot move panel (b)'s numerator.
    d["unscored"] = (d["stack"] == ROUTE) & ~d.bif_rescored
    assert d.unscored.sum() == 3 and d.loc[d.unscored, "cand"].sum() == 0
    assert not d.loc[d["stack"] != ROUTE, "bif_rescored"].any(), \
        "a non-BoltzIF route is flagged rescored; the source changed"

    d["corrected_pass"] = gate(d.sc_rmsd_af3, d.sc_rmsd_esm, d.dockq_af3, d.dockq_esm)
    return d


def gate(rmsd_a, rmsd_e, dockq_a, dockq_e):
    """The applied gate's two structural terms; either backend passes either term."""
    return (((rmsd_a <= RMSD_GATE) | (rmsd_e <= RMSD_GATE))
            & ((dockq_a >= DOCKQ_GATE) | (dockq_e >= DOCKQ_GATE)))


def strata(d):
    """Panel (b)'s three groups, in the order they are drawn."""
    s = np.where(d["stack"] != ROUTE, "gated",
                 np.where(d.corrected_pass, "bif_pass", "bif_fail"))
    return pd.Series(s, index=d.index)


def tally(d, by):
    """candidates / members per group, as the panels label them: n/N."""
    t = d.groupby(by, sort=False)["cand"].agg(cand="sum", n="size")
    t["pct"] = 100 * t.cand / t.n
    return t


# --- statistics (hand-rolled; scipy raises on import in this venv) ----------
def chisq(obs):
    """Pearson chi-square of independence with a closed-form upper tail.

    Even df: P(X >= x) = exp(-x/2) * sum_{i<df/2} (x/2)^i / i!.  df = 1: erfc(sqrt(x/2)).
    """
    obs = np.asarray(obs, float)
    rt, ct = obs.sum(1, keepdims=True), obs.sum(0, keepdims=True)
    e = rt @ ct / obs.sum()
    chi2 = float(((obs - e) ** 2 / e).sum())
    df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    h = chi2 / 2
    if df == 1:
        p = erfc(h ** 0.5)
    else:
        assert df % 2 == 0, "the closed form needs df = 1 or an even df"
        p = exp(-h) * sum(h ** i / factorial(i) for i in range(df // 2))
    return chi2, df, p, float(e.min())


def chi2_stat(obs):
    return chisq(obs)[0]


def perm_p(sizes, n_cand, obs_chi2, n=N_PERM, seed=SEED):
    """Permute the candidate LABELS over the panel; p = fraction reaching obs_chi2.

    Warranted here because the asymptotic test's smallest expected candidate count is
    2.14. Monte-Carlo: SE at n = 20,000 is about 0.003, so only two decimals are real.
    """
    labels = np.repeat(np.arange(len(sizes)), sizes)
    cand = np.zeros(labels.size, bool)
    cand[:n_cand] = True
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n):
        c = np.bincount(labels[rng.permutation(cand)], minlength=len(sizes))
        hits += chi2_stat(np.column_stack([np.asarray(sizes) - c, c])) >= obs_chi2 - 1e-9
    return hits / n, hits


def fisher(a, b, c, d):
    """Two-sided Fisher exact p for [[a, b], [c, d]], point-probability convention.

    For a 2x2 with fixed margins this IS the exact permutation test, so no Monte-Carlo
    is needed alongside it.
    """
    n, r1, c1 = a + b + c + d, a + b, a + c
    pr = [comb(r1, k) * comb(n - r1, c1 - k) / comb(n, c1)
          for k in range(max(0, c1 - (n - r1)), min(r1, c1) + 1)]
    p0 = comb(r1, a) * comb(n - r1, c1 - a) / comb(n, c1)
    return sum(x for x in pr if x <= p0 * (1 + 1e-9))


# --- panels -----------------------------------------------------------------
def bars(ax, t, colors, ref, labels, ymax, dy=2):
    """One bar per group, labelled n/N (STYLE_NOTES section 5), reference hairline.

    `dy` lifts the n/N labels clear of the reference hairline where a bar top sits
    within a few points of it, which in panel (b) two of three do.
    """
    x = np.arange(len(t)).astype(float)
    ax.bar(x, t.pct, 0.62, color=colors, edgecolor=palette.INK, linewidth=0.4)
    ax.axhline(ref, color=palette.NEUTRAL["600"], lw=0.5, zorder=0)
    for xi, (_, r) in zip(x, t.iterrows()):
        ax.annotate(f"{int(r.cand)}/{int(r.n):,}", (xi, r.pct),
                    textcoords="offset points", xytext=(0, dy), ha="center",
                    fontsize=style.CFG["font_size_pt"])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Candidates (%)")
    ax.set_ylim(0, ymax)
    ax.tick_params(axis="both", length=2, width=0.5)
    return x


def note(ax, text, xy=(0.0, 1.0), off=(2, -2), ha="left", color=None):
    ax.annotate(text, xy, xycoords="axes fraction", xytext=off,
                textcoords="offset points", va="top", ha=ha,
                fontsize=style.CFG["font_size_pt"], color=color or palette.INK)


def panel_a(d, out):
    """Candidate recovery per design route. The weakly gated route is the highest."""
    t = tally(d, "stack").reindex(ORDER)
    overall = 100 * d.cand.sum() / len(d)
    chi2, df, p, min_exp = chisq(t[["n", "cand"]].assign(
        miss=lambda x: x.n - x.cand)[["miss", "cand"]].to_numpy(float))
    pp, _ = perm_p(t.n.to_numpy(), int(d.cand.sum()), chi2)

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(40)))
    bars(ax, t, [palette.stack_color(k) for k in ORDER], overall, ORDER, 0.52)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.spines["left"].set_bounds(0, 0.4)
    note(ax, f"$\\chi^2$ = {chi2:.2f}, df = {df}, $p$ = {p:.2f}; "
             f"permutation $p$ = {pp:.2f}")
    note(ax, f"all routes {overall:.3f}%", xy=(1.0, 1.0), off=(-2, -2), ha="right",
         color=palette.NEUTRAL["600"])
    fig.subplots_adjust(left=0.11, right=0.99, top=0.97, bottom=0.13)
    style.save(fig, out, "panel_a_candidates_by_route")
    plt.close(fig)
    return t, (chi2, df, p, min_exp, pp)


def panel_b(d, out):
    """The same recovery, split by CORRECTED gate status rather than route.

    The failing stratum is entirely BoltzGen + BoltzIF, so route and gate status are
    confounded by construction; splitting that route in two is as far as the panel can
    separate them, and it is why both BoltzIF bars keep the route colour.
    """
    e = d[~d.unscored]
    t = tally(e.assign(strat=strata(e)), "strat").reindex(
        ["gated", "bif_pass", "bif_fail"])
    overall = 100 * e.cand.sum() / len(e)
    fail, keep = t.loc["bif_fail"], t.loc[["gated", "bif_pass"]].sum()
    pf = fisher(int(fail.n - fail.cand), int(fail.cand),
                int(keep.n - keep.cand), int(keep.cand))

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(40)))
    colors = [palette.NEUTRAL["400"]] + [palette.stack_color(ROUTE)] * 2
    bars(ax, t, colors, overall, ["Four gated\nroutes", "BG+BIF\npasses",
                                  "BG+BIF\nfails"], 0.58, dy=6)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.spines["left"].set_bounds(0, 0.4)
    note(ax, f"fails vs passes, Fisher $p$ = {pf:.2f}")
    note(ax, f"whole panel {overall:.3f}%", xy=(1.0, 1.0), off=(-2, -2), ha="right",
         color=palette.NEUTRAL["600"])
    fig.subplots_adjust(left=0.11, right=0.99, top=0.97, bottom=0.21)
    style.save(fig, out, "panel_b_candidates_by_gate_status")
    plt.close(fig)
    return t, pf


# --- (c), tested and NOT drawn ----------------------------------------------
def panel_c_check(d):
    """Does corrected ligand_RMSD or sc_DockQ track candidate recovery across the panel?

    It does not, so no panel is emitted (the brief made (c) conditional). The whole-panel
    sc_DockQ association is a route-mix artefact: BoltzGen + BoltzIF has both the lowest
    corrected sc_DockQ and the highest candidate rate, and the association dies inside
    either route group. Numbers printed so the README can quote them.
    """
    e = d[~d.unscored].copy()
    # Best-of-two is the gate's own quantity: either backend may pass either term.
    e["rmsd_best"] = e[["sc_rmsd_af3", "sc_rmsd_esm"]].min(axis=1)
    e["dockq_best"] = e[["dockq_af3", "dockq_esm"]].max(axis=1)
    rows = []
    for name, sub in (("whole panel", e),
                      ("four gated routes", e[e["stack"] != ROUTE]),
                      ("BG+BIF only", e[e["stack"] == ROUTE])):
        for col in ("rmsd_best", "dockq_best"):
            rows.append((name, col, *rank_perm(sub, col),
                         sub.loc[sub.cand == 1, col].median(),
                         sub.loc[sub.cand == 0, col].median()))
    return pd.DataFrame(rows, columns=["group", "metric", "p", "k", "n",
                                       "median_cand", "median_rest"])


def rank_perm(sub, col, n=N_PERM, seed=SEED):
    """Two-sided permutation test on the candidates' mean rank of `col`."""
    r = sub[col].rank().to_numpy()
    k = int(sub.cand.sum())
    obs, centre = r[sub.cand == 1].mean(), r.mean()
    rng = np.random.default_rng(seed)
    null = np.array([r[rng.choice(r.size, k, replace=False)].mean()
                     for _ in range(n)])
    return float((np.abs(null - centre) >= abs(obs - centre) - 1e-9).mean()), k, r.size


# --- report ------------------------------------------------------------------
def report(d, a, b, c):
    t_a, (chi2, df, p, min_exp, pp) = a
    t_b, pf = b
    fmt = lambda t: "\n".join(
        f"    {i:<10s} {int(r.cand):>3d}/{int(r.n):<6,d} {r.pct:.3f}%"
        for i, r in t.iterrows())
    print(f"panel (a) candidates per route, n = {len(d):,}, {int(d.cand.sum())} candidates")
    print(fmt(t_a))
    print(f"    chi2 = {chi2:.4f}, df = {df}, p = {p:.4f}; smallest expected "
          f"candidate count {min_exp:.2f}")
    print(f"    {N_PERM:,} label permutations (seed {SEED}): p = {pp:.4f} "
          f"(Monte-Carlo SE ~{(pp * (1 - pp) / N_PERM) ** 0.5:.4f})")
    e = d[~d.unscored]
    fail, keep = t_b.loc["bif_fail"], t_b.loc[["gated", "bif_pass"]].sum()
    bif_n = int(t_b.loc[["bif_pass", "bif_fail"], "n"].sum())
    print(f"\npanel (b) corrected gate, n = {len(e):,} "
          f"({int(d.unscored.sum())} BG+BIF lack a true-backbone score)")
    print(f"    BG+BIF failing the corrected gate: {int(fail.n)}/{bif_n} "
          f"({100 * fail.n / bif_n:.1f}%) - matches bif_effect")
    print(fmt(t_b))
    print(f"    fails {int(fail.cand)}/{int(fail.n):,} vs passes {int(keep.cand)}/"
          f"{int(keep.n):,}, ratio {(fail.cand / fail.n) / (keep.cand / keep.n):.2f}x, "
          f"Fisher p = {pf:.4f}")
    print("    3x2 chi2 = %.4f, df = %d, p = %.4f, smallest expected %.2f"
          % chisq(t_b[["n", "cand"]].assign(miss=lambda x: x.n - x.cand)
                  [["miss", "cand"]].to_numpy(float)))
    bp = t_b.loc["bif_pass"]
    print(f"    within BG+BIF, Fisher p = "
          f"{fisher(int(fail.n - fail.cand), int(fail.cand), int(bp.n - bp.cand), int(bp.cand)):.4f}")
    print("\n(c) tested, NOT drawn - corrected metric vs candidate status:")
    print(c.to_string(index=False,
                      formatters={"p": "{:.4f}".format,
                                  "median_cand": "{:.3f}".format,
                                  "median_rest": "{:.3f}".format}))


def main():
    style.apply()
    out = paths.fig_output(FIG)
    d = load()
    report(d, panel_a(d, out), panel_b(d, out), panel_c_check(d))


if __name__ == "__main__":
    main()
