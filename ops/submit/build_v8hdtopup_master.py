#!/usr/bin/env python3
"""build_v8hdtopup_master.py — fold every sidecar into the v8 HD topup table, reproducing
csy4_master_foot.csv semantics exactly.

Every formula below was verified empirically against real csy4_master_foot.csv rows before
this script was written (300 rows for DockQ, 2000 for ipTM, 6 for foot) — none is guessed:

  dockq_std_LRMSD_<be> = dockq_LRMSD_<be> as emitted by dockq_self.py (STOCK convention)
  dockq_LRMSD_<be>     = dockq_std_LRMSD_<be> * (lrmsd_rna_<be>/lrmsd_prot_<be>)   [300/300]
  dockq_fnat_<be>     *= fnat_factor                                (merge_final.py contract)
  dockq_self_<be>      = dq(fnat', iRMSD, dockq_LRMSD_<be>)                        [300/300]
  dockq_std_self_<be>  = dq(fnat', iRMSD, dockq_std_LRMSD_<be>)                    [300/300]
  iptm_min3 = min3 ; iptm_best2 = MEDIAN of 3 ; mean_iptm3 = mean3                 [2000/2000]
  consensus_iptm = min3 ; *_legacy = min/mean over (af3,b2) only                   [2000/2000]
  foot_<be> = foot_tp_min_dist (TRANSPLANT estimator, user's choice 2026-08-13)    [6/6]

Rows with no fnat_factor keep their DockQ uncorrected (never silently defaulted to 1.0) —
merge_final.py's contract. b2 is absent by instruction, so ==3-gated bare consensus stays
NaN; optH/ipTM bare columns use skipna min/max/mean over the backends present, matching
opth_integrate.py.
"""
import csv, glob, os, sys
import numpy as np, pandas as pd

CC = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
BE = ["af3", "esm"]
FOOT_LO, FOOT_HI = 2138, 2206
MOIETY = ["dsasa_rna_total", "dsasa_rna_backbone", "dsasa_rna_phosphate", "dsasa_rna_sugar",
          "dsasa_rna_o2prime", "dsasa_rna_base", "dsasa_rna_base_major", "dsasa_rna_base_minor",
          "dsasa_rna_base_wc", "dsasa_rna_base_unassigned", "frac_dsasa_base", "frac_dsasa_major"]
OPTH_SRC = {"hb_opth": "n_h_bonds_interface_opth", "buns_ia_opth": "delta_unsat_hbonds_opth",
            "buns_vsasa_raw": "buns_vsasa", "buns_vsasa_opth": "buns_vsasa_opth",
            "n_polar_h_moved": "n_polar_h_moved"}
OPTH_MIN = ["n_h_bonds_interface_opth"]
OPTH_MAX = ["delta_unsat_hbonds_opth", "buns_vsasa", "buns_vsasa_opth"]
dq = lambda f, i, l: (f + 1 / (1 + (i / 1.5) ** 2) + 1 / (1 + (l / 8.5) ** 2)) / 3

m = pd.read_csv(CC + "v8hdtopup_fullmetrics.csv", low_memory=False)
N = len(m)
print(f"base: {N:,} rows x {m.shape[1]} cols", flush=True)
assert m.design_id.is_unique, "duplicate design_id in base"

def cat(pattern, usecols=None):
    fs = sorted(glob.glob(pattern))
    if not fs: sys.exit(f"FATAL: no files match {pattern}")
    d = pd.concat([pd.read_csv(f, usecols=usecols, low_memory=False) for f in fs],
                  ignore_index=True).drop_duplicates("design_id")
    return d

def join(d, label):
    global m
    before = len(m)
    m = m.merge(d, on="design_id", how="left")
    assert len(m) == before, f"{label}: row count changed {before} -> {len(m)}"
    print(f"  + {label}: {d.shape[1]-1} cols, matched "
          f"{m[d.columns[1]].notna().sum():,}/{N:,}", flush=True)

# ---- SASA moiety (per backend) -----------------------------------------------------
for be in BE:
    d = cat(CC + f"sasa_v8hdtopup_{be}/sasa_*.csv", usecols=["design_id"] + MOIETY)
    join(d.rename(columns={c: f"{c}_{be}" for c in MOIETY}), f"sasa_{be}")

# ---- fnat factor -------------------------------------------------------------------
join(cat(CC + "fnatf_v8hdtopup/fnatf_*.csv", usecols=["design_id", "n_rna_res", "fnat_factor"]),
     "fnat")

# ---- RNA-aligned LRMSD ratio -------------------------------------------------------
lr = [f"lrmsd_rna_{b}" for b in BE] + [f"lrmsd_prot_{b}" for b in BE]
join(cat(CC + "lrmsd_rna_v8hdtopup/lrmsd_*.csv", usecols=["design_id"] + lr), "lrmsd_rna")

# ---- optH (with the raw-identity audit opth_integrate.py performs) ------------------
for be in BE:
    d = pd.read_csv(CC + f"opth_v8hdtopup_{be}.csv",
                    usecols=["design_id"] + list(OPTH_SRC) + ["hb_raw"], low_memory=False)
    d = d.drop_duplicates("design_id").rename(
        columns={k: f"{v}_{be}" for k, v in OPTH_SRC.items()} | {"hb_raw": f"__chk_{be}"})
    join(d, f"opth_{be}")
    k = m[f"__chk_{be}"].notna() & m[f"n_h_bonds_interface_{be}"].notna()
    bad = int((m.loc[k, f"__chk_{be}"] != m.loc[k, f"n_h_bonds_interface_{be}"]).sum())
    print(f"    audit hb_raw == n_h_bonds_interface_{be}: {int(k.sum()):,} checked, {bad} mismatches",
          flush=True)
    assert bad == 0, f"optH join is stale/misaligned for {be}"
    m = m.drop(columns=[f"__chk_{be}"])
for c in OPTH_MIN:
    m[c] = m[[f"{c}_{b}" for b in BE]].min(axis=1)
for c in OPTH_MAX:
    m[c] = m[[f"{c}_{b}" for b in BE]].max(axis=1)

# ---- foot (transplant estimator) ---------------------------------------------------
for be in BE:
    fs = [CC + f"foot_clash/out/foot_{be}_{i:05d}.csv" for i in range(FOOT_LO, FOOT_HI + 1)]
    d = pd.concat([pd.read_csv(f, usecols=["design_id", "foot_tp_min_dist"], low_memory=False)
                   for f in fs if os.path.exists(f)], ignore_index=True)
    d = d.drop_duplicates("design_id").rename(columns={"foot_tp_min_dist": f"foot_{be}"})
    join(d, f"foot_{be}")

# ---- DockQ: stock -> RNA-aligned, then the fnat correction --------------------------
ff = pd.to_numeric(m["fnat_factor"], errors="coerce")
have_ff = ff.notna()
print(f"\nfnat_factor present: {int(have_ff.sum()):,}/{N:,} "
      f"({100*have_ff.mean():.2f}%) — the rest keep DockQ uncorrected", flush=True)
for be in BE:
    m[f"dockq_std_LRMSD_{be}"] = m[f"dockq_LRMSD_{be}"]                    # stock, preserved
    ratio = m[f"lrmsd_rna_{be}"] / m[f"lrmsd_prot_{be}"].where(m[f"lrmsd_prot_{be}"] > 0.01)
    m[f"dockq_LRMSD_{be}"] = m[f"dockq_std_LRMSD_{be}"] * ratio            # RNA-aligned
    m.loc[have_ff, f"dockq_fnat_{be}"] = m.loc[have_ff, f"dockq_fnat_{be}"] * ff[have_ff]
    m[f"dockq_self_{be}"] = dq(m[f"dockq_fnat_{be}"], m[f"dockq_iRMSD_{be}"],
                               m[f"dockq_LRMSD_{be}"])
    m[f"dockq_std_self_{be}"] = dq(m[f"dockq_fnat_{be}"], m[f"dockq_iRMSD_{be}"],
                                   m[f"dockq_std_LRMSD_{be}"])
# bare consensus stays ==3-gated => NaN without b2, matching consensus_integrate
for c in ("dockq_fnat", "dockq_self", "dockq_std_self"):
    m[c] = np.nan

# ---- ipTM consensus (b2 ipTM is free from the refold, no b2 metrics needed) ---------
A, B, E = "af3_msaoff_iptm", "b2_iptm", "esm_iptm"
tri = m[[A, B, E]].astype(float)
m["iptm_min3"] = tri.min(axis=1)
m["iptm_best2"] = tri.median(axis=1)
m["mean_iptm3"] = tri.mean(axis=1)
m["consensus_iptm"] = m["iptm_min3"]
m["consensus_iptm_legacy"] = m[[A, B]].astype(float).min(axis=1)
m["mean_iptm_legacy"] = m[[A, B]].astype(float).mean(axis=1)
m["source_table"] = "v8hdtopup"

out = CC + "v8hdtopup_master.csv"
assert len(m) == N and m.design_id.is_unique
m.to_csv(out, index=False)
print(f"\nwrote {out}: {len(m):,} rows x {m.shape[1]} cols", flush=True)
