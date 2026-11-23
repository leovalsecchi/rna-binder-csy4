#!/usr/bin/env python
"""Precision-recall of the Fig 5 metrics, proposed binders vs ONT-negative designs, and
AF3-vs-ESM agreement for the 15 / 32 / 12k.

Classes (chosen by Leo, do not change):
  positives  flow FC (ONOFF_avg) >= 1.2 among the 32 ordered           -> 15
  negatives  ALL HTS non-hits with > 75 reads: gate-eligible designs (inputs/poolcraft_gates/
             gate3_positive.csv, reads > 75 on both the on-target and the negative-control side)
             that never pass the hit gate (reads>=75, splicing>=0.07, FC>=10) in any comparison,
             minus the 32 ordered (the hits themselves, so classes are disjoint).
             Broader than the v1 "splicing < 0.07" rule, which dropped designs that splice
             well but have a low fold change.
No rescaling or subsampling of negatives; baseline precision = prevalence.
Direction is fixed a priori (the direction the field ranks on) and the opposite
direction's AP is stored too. Each backend is scored on its own: never min/max of the
two (selecting the better backend is circular, it gave a spurious P=0.038 once).

numpy only: scipy/sklearn do not import in ~/.venvs/binder_map (propack dlopen error).

    ~/.venvs/binder_map/bin/python code/Fig5pr_precision_recall/export_data.py
"""
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
IN = REPO / "inputs"
# v2 outputs (all ten metrics, gate-correct negatives, BG+BIF true-backbone ligand-RMSD/DockQ); the v1 files stay untouched
OUT_PR, OUT_AP = IN / "fig5pr_v2_pr_curves.csv", IN / "fig5pr_v2_ap_summary.csv"
OUT_AGR, OUT_IDS = IN / "fig5pr_v2_backend_agreement.csv", IN / "fig5pr_v2_class_ids.csv"
OUT_SENS, OUT_COV = IN / "fig5pr_v2_ap_sensitivity.csv", IN / "fig5pr_v2_coverage.csv"
for p in (OUT_PR, OUT_AP, OUT_AGR, OUT_IDS, OUT_SENS, OUT_COV):
    assert not p.exists(), f"refusing to overwrite {p}"
rng = np.random.default_rng(0)


def _prep(v):
    o = np.argsort(-v, kind="stable")
    sv = v[o]
    return o, np.r_[np.where(np.diff(sv) != 0)[0], len(v) - 1]  # last index of each distinct threshold


def AP(y, v):
    """Average precision, sklearn convention (ties grouped at one threshold)."""
    o, ends = _prep(v)
    cs = np.cumsum(y[o])[ends]
    return float(np.sum(np.diff(cs, prepend=0) / y.sum() * cs / (ends + 1)))


def pr_curve(y, v):
    o, ends = _prep(v)
    cs = np.cumsum(y[o])[ends]
    return np.r_[0, cs / y.sum()], np.r_[1, cs / (ends + 1)]


def rho(a, b):
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])



# ---------------- classes
lib = pd.read_csv(IN / "csy4_binder_library_12k_v6_k30.csv")
flow = pd.read_csv(IN / "lv2_34_1_flow_on_off.csv").dropna(subset=["design_id"])
g3 = pd.read_csv(IN / "poolcraft_gates" / "gate3_positive.csv")
sp = pd.read_csv(IN / "screen_points.csv")
sp["design_id"] = sp.transcript.str.split("::").str[-1]
ordered = set(flow.design_id)
pos = set(flow.design_id[flow.ONOFF_avg >= 1.2])
hits = set(g3.design[g3.quadrant.str.startswith("HITS")])
assert len(ordered) == 32 and hits == ordered, "gate3 HITS must be the 32 ordered"
g3p = g3[(g3.on_reads > 75) & (g3.neg_reads > 75)]
pool = set(g3p.design)
neg_broad = (pool - hits) - ordered
neg_narrow = (set(g3p.design[g3p.on_target_splicing_ratio < 0.07])) - ordered
# sensitivity pools: (a) v1 rule on screen_points (on-target reads only, spl < 0.07);
# (b) screen_points, on-target reads only, broad definition (not a hit in the screen)
sp_hit = (sp.reads >= 75) & (sp.neg_reads >= 75) & (sp.splicing_ratio >= .07) & (sp.fold_change_raw >= 10)
v1_neg = set(sp.design_id[(sp.reads > 75) & (sp.splicing_ratio < 0.07)]) - ordered
sp_broad = set(sp.design_id[sp.reads > 75]) - set(sp.design_id[sp_hit]) - ordered
assert not (pos & neg_broad) and (pos | neg_broad) <= set(lib.design_id)
print(f"positives (proposed binders)             n = {len(pos)}")
print(f"gate-eligible designs (reads>75 both sides) {len(pool)}  = 32 hits + {len(pool - hits)} non-hits")
print(f"NEGATIVES, broad (non-hit, reads>75)      n = {len(neg_broad)}   <- used")
print(f"  narrow (splicing<0.07) on same pool     n = {len(neg_narrow)}   (broad adds {len(neg_broad - neg_narrow)})")
print(f"  ordered designs removed from the pools: broad {len((pool - hits) & ordered)}, narrow "
      f"{len(set(g3p.design[g3p.on_target_splicing_ratio < 0.07]) & ordered)}")
print(f"  sensitivity: v1 rule on screen_points n = {len(v1_neg)}; screen_points broad n = {len(sp_broad)}")

# ---------------- metrics: v1 library columns + BG+BIF true-backbone fix + the five recovered ones
lib = lib.merge(pd.read_csv(IN / "large" / "bif_true_per_backend.csv.gz"), on="design_id", how="left")
isbif = lib.inverse_folder.eq("BIF")
for be in ("af3", "esm"):  # stored BG+BIF ligand_RMSD/DockQ are wrong (design CIF = a re-prediction); use true-backbone
    lib.loc[isbif, f"sc_rmsd_{be}"] = lib.loc[isbif, f"rmsd_{be}_true"]
    lib.loc[isbif, f"dockq_{be}"] = lib.loc[isbif, f"dockq_{be}_true"]
rec = pd.read_csv(IN / "fig5pr_recovered_metrics_12k.csv.gz")
lib = lib.merge(rec.drop(columns=[c for c in ("stack",) if c in rec]), on="design_id", how="left", suffixes=("", "_m"))
lib["asn_seq"] = lib.sequence.str.count("N") / lib.sequence.str.len()

# (label, short, lib column stem, +1 higher=binder / -1 lower, primary)
SPEC = [("ipTM", "iptm", 1, 1), ("ligand-RMSD", "sc_rmsd", -1, 1), ("DockQ", "dockq", 1, 1),
        ("Interface ΔSASA", "dsasa", 1, 1), ("Interface H-bonds", "n_hbonds", 1, 1),
        ("Heavy contacts", "n_heavy_contacts_4A", 1, 1), ("Shape compl.", "sc_value", 1, 1),
        ("Packstat", "packstat", 1, 1), ("BUNS", "delta_unsat_hbonds", -1, 1),
        ("π-stacking", "pi_stack_score", 1, 1),
        ("Heavy contacts 4.5 Å", "n_heavy_contacts_4p5A", 1, 0), ("BUNS (optH)", "delta_unsat_hbonds_opth", -1, 0),
        ("BUNS (buried polar VSASA)", "buns_vsasa", -1, 0)]
M = []
for lab, col, s, prim in SPEC:
    for be in ("af3", "esm"):
        M.append((lab, be.upper(), f"{col}_{be}", s, prim))
M.append(("Asn (sequence)", "seq", "asn_seq", 1, 1))
METRICS_PRIMARY = [m[0] for m in SPEC if m[3]] + ["Asn (sequence)"]

def prep(neg):
    D = lib[lib.design_id.isin(pos | neg)].copy()
    D["y"] = D.design_id.isin(pos).astype(int)
    return D

def score(D, boot):
    """AP, permutation P, (optionally) bootstrap interval for every metric; missing values rank last."""
    y = D.y.to_numpy(); base = y.mean(); P = 3000
    ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
    curves, rows = [], []
    for lab, be, col, s, prim in M:
        x = D[col].to_numpy(float)
        miss = np.isnan(x)
        v = np.where(miss, np.nanmin(s * x) - 1, s * x)  # unscorable designs rank below every scored one
        ap, ap_rev = AP(y, v), AP(y, np.where(miss, np.nanmin(-s * x) - 1, -s * x))
        nul = np.array([AP(rng.permutation(y), v) for _ in range(P)])
        p = (np.sum(nul >= ap) + 1) / (P + 1)
        lo = hi = np.nan
        if boot:
            bs = []
            for _ in range(500):  # bootstrap within class
                i = np.r_[rng.choice(ip, len(ip)), rng.choice(ineg, len(ineg))]
                bs.append(AP(y[i], v[i]))
            lo, hi = np.percentile(bs, [2.5, 97.5])
            rc, pr = pr_curve(y, v)
            curves.append(pd.DataFrame(dict(label=lab, backend=be, recall=rc, precision=pr)))
        rows.append(dict(label=lab, backend=be, col=col, primary=prim, direction="higher" if s > 0 else "lower",
                         n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                         n_pos_scored=int((~miss & (y == 1)).sum()), n_neg_scored=int((~miss & (y == 0)).sum()),
                         baseline=base, ap=ap, ap_over_baseline=ap / base, ap_opposite_direction=ap_rev,
                         ap_lo=lo, ap_hi=hi, p_perm=p))
    return pd.DataFrame(rows), curves

D = prep(neg_broad)
y = D.y.to_numpy()
pd.DataFrame({"design_id": D.design_id, "cls": np.where(D.y == 1, "pos", "neg")}).to_csv(OUT_IDS, index=False)
cov = pd.DataFrame([dict(label=l, backend=b, col=c, primary=pr,
                         frac_pos=float(D.loc[D.y == 1, c].notna().mean()), frac_neg=float(D.loc[D.y == 0, c].notna().mean()),
                         frac_32=float(lib.loc[lib.design_id.isin(ordered), c].notna().mean()),
                         frac_12k=float(lib[c].notna().mean())) for l, b, c, s, pr in M])
cov.to_csv(OUT_COV, index=False)
print("\ncoverage (fraction non-missing) of every scored column:\n", cov.round(4).to_string(index=False))

R, curves = score(D, True)
pd.concat(curves).to_csv(OUT_PR, index=False)
R.to_csv(OUT_AP, index=False)
base = y.mean()
print(f"\nclasses: {int(y.sum())} pos / {int((1 - y).sum())} neg, baseline AP = {base:.5f}")
print("ranked by AP (a-priori direction), primary metrics:")
Rp = R[R.primary == 1]
print(Rp.sort_values("ap", ascending=False)[["label", "backend", "direction", "ap", "ap_over_baseline",
      "ap_lo", "ap_hi", "p_perm", "ap_opposite_direction"]].round(4).to_string(index=False))
print(f"P<0.05 expected by chance among {len(Rp)}: {0.05 * len(Rp):.2f}; observed {(Rp.p_perm < 0.05).sum()}")
print("secondary (sensitivity) variants:\n", R[R.primary == 0][["label", "backend", "ap", "ap_over_baseline", "p_perm"]].round(4).to_string(index=False))

# sensitivity: same scoring against the other negative definitions (no bootstrap)
S = []
for name, neg in [("narrow_splicing_lt_0.07", neg_narrow), ("v1_screen_points_onreads_spl", v1_neg),
                  ("screen_points_broad_onreads", sp_broad)]:
    Rs, _ = score(prep(neg), False)
    Rs.insert(0, "neg_definition", name)
    S.append(Rs)
S = pd.concat(S)
S.to_csv(OUT_SENS, index=False)
for name, d in S.groupby("neg_definition"):
    dp = d[d.primary == 1]
    print(f"\n[{name}] n_neg={int(d.n_neg.iloc[0])} baseline={d.baseline.iloc[0]:.5f}  P<0.05: {(dp.p_perm < 0.05).sum()}/{len(dp)}  min P={dp.p_perm.min():.4f} ({dp.loc[dp.p_perm.idxmin(), 'label']} {dp.loc[dp.p_perm.idxmin(), 'backend']})")

# ---------------- deliverable 2: AF3 vs ESM agreement
G = {"15 proposed binders": lib[lib.design_id.isin(pos)],
     "32 ordered": lib[lib.design_id.isin(ordered)], "12k library": lib}
rows = []
for lab, col in [(l, c) for l, c, s, pr in SPEC if pr]:
    L = lib[[f"{col}_af3", f"{col}_esm"]].dropna().to_numpy(float)
    lr = rho(L[:, 0], L[:, 1])
    for g, d in G.items():
        a, b = d[f"{col}_af3"].to_numpy(float), d[f"{col}_esm"].to_numpy(float)
        ok = ~(np.isnan(a) | np.isnan(b))  # pairwise-complete
        a, b = a[ok], b[ok]
        r = rho(a, b)
        rr = []
        for _ in range(4000 if len(d) < 1000 else 300):
            i = rng.integers(0, len(a), len(a))
            rr.append(rho(a[i], b[i]))
        lo, hi = np.nanpercentile(rr, [2.5, 97.5])
        pct = np.nan  # where the subset rho sits among random same-size library subsets
        if g != "12k library":
            sub = [rho(L[i, 0], L[i, 1]) for i in (rng.choice(len(L), len(a), replace=False) for _ in range(2000))]
            pct = float(np.mean(np.array(sub) <= r))
        rows.append(dict(metric=lab, group=g, n=len(a), spearman=r, lo=lo, hi=hi,
                         pearson=float(np.corrcoef(a, b)[0, 1]), lib_spearman=lr,
                         pct_of_random_library_subsets_below=pct))
A = pd.DataFrame(rows)
A.to_csv(OUT_AGR, index=False)
print("\nAF3 vs ESM Spearman (95% bootstrap CI):")
print(A.round(3).to_string(index=False))
