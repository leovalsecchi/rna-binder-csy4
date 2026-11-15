#!/usr/bin/env python3
"""opth_compare.py — what does optimising hydrogens change?

Reads the paired experiment output and answers, in order:
  1. does the _raw pass reproduce the production columns? (sanity)
  2. how big is the optH effect on interface H-bonds and on both BUNS flavours?
  3. does the ranking survive it? (raw vs optH correlation)
  4. do the associations with DockQ / ipTM / interface size change?
  5. does the size-banded BUNS -> DockQ discriminator get stronger with optH?
"""
import numpy as np
import pandas as pd
from scipy import stats

RB = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
d = pd.read_csv(RB + "opth_compare_2k.csv", low_memory=False)
d = d[d["opth_ok"] == 1].copy()
print("rows with paired result: %d / 2000" % len(d))

dq = pd.read_csv(RB + "dockq_rna_v2v7.csv", usecols=["design_id", "dockq_rna_af3"],
                 low_memory=False)
gr = pd.read_csv(RB + "app_metrics_100k_groove.csv",
                 usecols=["design_id", "frac_dsasa_base_af3"], low_memory=False)
d = d.merge(dq, on="design_id", how="left").merge(gr, on="design_id", how="left")
print("joined dockq_rna_af3: %d, frac_dsasa_base: %d"
      % (d.dockq_rna_af3.notna().sum(), d.frac_dsasa_base_af3.notna().sum()))

# ---- 1. sanity: does _raw reproduce production? -----------------------------
print("\n=== 1. _raw vs production columns (must be identical) ===")
for a, b in [("hb_raw", "n_h_bonds_interface_af3"), ("buns_ia_raw", "delta_unsat_hbonds_af3"),
             ("dsasa_total_raw", "dsasa_total_af3")]:
    x, y = d[a].astype(float), d[b].astype(float)
    m = x.notna() & y.notna()
    same = int((np.abs(x[m] - y[m]) < 0.01).sum())
    print("  %-16s vs %-26s identical %d/%d  (max |diff| %.2f)"
          % (a, b, same, int(m.sum()), float(np.abs(x[m] - y[m]).max()) if m.sum() else 0))

# ---- 2. size of the optH effect ---------------------------------------------
print("\n=== 2. effect of optimising hydrogens (paired, same pose) ===")
print("  polar H moved >0.1 A: median %d  mean %.1f  (max %d)"
      % (d.n_polar_h_moved.median(), d.n_polar_h_moved.mean(), d.n_polar_h_moved.max()))
for lab, a, b in [("interface H-bonds", "hb_raw", "hb_opth"),
                  ("BUNS (IA legacy)", "buns_ia_raw", "buns_ia_opth"),
                  ("BUNS (VSASA 1.1)", "buns_vsasa_raw", "buns_vsasa_opth")]:
    x, y = d[a].astype(float), d[b].astype(float)
    m = x.notna() & y.notna()
    x, y = x[m], y[m]
    dl = y - x
    w = stats.wilcoxon(x, y) if (dl != 0).any() else None
    print("\n  %s   n=%d" % (lab, len(x)))
    print("    raw   median %6.1f  mean %7.2f" % (x.median(), x.mean()))
    print("    optH  median %6.1f  mean %7.2f" % (y.median(), y.mean()))
    print("    delta median %6.1f  mean %+7.2f  | changed %.1f%%  up %.1f%%  down %.1f%%"
          % (dl.median(), dl.mean(), 100 * (dl != 0).mean(),
             100 * (dl > 0).mean(), 100 * (dl < 0).mean()))
    print("    mean relative change %+.1f%%   pearson(raw,optH) %.4f  spearman %.4f"
          % (100 * (dl.sum() / x.sum()) if x.sum() else float("nan"),
             stats.pearsonr(x, y)[0], stats.spearmanr(x, y)[0]))
    if w is not None:
        print("    wilcoxon signed-rank p = %.3g" % w.pvalue)

# ---- 3/4. association with quality, raw vs optH -----------------------------
print("\n=== 3. association with other metrics: does optH change the picture? ===")
targets = [("DockQ (RNA-aligned)", "dockq_rna_af3"), ("ipTM", "af3_msaoff_iptm"),
           ("dSASA_total", "dsasa_total_raw"), ("frac_dsasa_base", "frac_dsasa_base_af3")]
metrics = [("hb", "hb_raw", "hb_opth"),
           ("BUNS_ia", "buns_ia_raw", "buns_ia_opth"),
           ("BUNS_vsasa", "buns_vsasa_raw", "buns_vsasa_opth")]
print("  %-12s %-22s %10s %10s %10s" % ("metric", "vs", "r(raw)", "r(optH)", "shift"))
for mlab, a, b in metrics:
    for tlab, t in targets:
        x, y, q = d[a].astype(float), d[b].astype(float), d[t].astype(float)
        m = x.notna() & y.notna() & q.notna()
        if m.sum() < 100:
            continue
        r1 = stats.spearmanr(x[m], q[m])[0]
        r2 = stats.spearmanr(y[m], q[m])[0]
        print("  %-12s %-22s %10.3f %10.3f %10.3f" % (mlab, tlab, r1, r2, r2 - r1))

# ---- 5. size-banded BUNS -> DockQ, raw vs optH ------------------------------
print("\n=== 4. size-banded BUNS -> DockQ pass rate (the discriminator) ===")
d["band"] = pd.cut(d.dsasa_total_raw, [0, 1300, 1600, 1900, 1e9],
                   labels=["<1300", "1300-1600", "1600-1900", ">1900"])
for col, lab in [("buns_ia_raw", "raw"), ("buns_ia_opth", "optH")]:
    rows = []
    for bnd, g in d.groupby("band", observed=True):
        g = g.dropna(subset=[col, "dockq_rna_af3"])
        if len(g) < 120:
            continue
        try:
            q = pd.qcut(g[col], 2, labels=["low", "high"], duplicates="drop")
        except ValueError:
            continue
        for lv, gg in g.groupby(q, observed=True):
            rows.append(dict(band=bnd, BUNS=lv, n=len(gg),
                             med=gg[col].median(),
                             dockq_med=round(gg.dockq_rna_af3.median(), 3),
                             pass04=round(100 * (gg.dockq_rna_af3 >= 0.4).mean(), 2)))
    print("\n  --- %s ---" % lab)
    print(pd.DataFrame(rows).to_string(index=False))

d.to_csv(RB + "opth_compare_2k_joined.csv", index=False)
print("\nwrote " + RB + "opth_compare_2k_joined.csv")
