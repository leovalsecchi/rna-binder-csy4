"""data_loader.py — build the master per-(design × refold_condition) table.

Loads the wide aggregated CSV, keeps MSA-OFF refolds only (AF3 msa-off + Boltz-2),
melts to long (one row per design × refold_condition), maps the real 65-column
schema onto the plotting schema, derives per-nt interface contacts from the
rna_contact_nt_* counts, and merges optional DockQ / per-nt-hbond sidecars.

Metrics not yet computed (dockq without a sidecar, pae, n_clashes,
fraction_disordered, b2 ptm, per-nt hbond_contacts) come through as NaN / empty
so the plots degrade gracefully and auto-populate once the columns exist.
"""
from __future__ import annotations
import ast
import numpy as np
import pandas as pd
from pathlib import Path

import config as cfg

# long-schema columns the plotting layer consumes
SCHEMA = [
    "design_id", "stack", "refold_condition",
    "iptm", "ptm", "dockq", "protein_rmsd", "pae",
    "fraction_disordered", "n_clashes", "n_hbonds",
    "helix_frac", "sheet_frac", "interface_contacts", "hbond_contacts",
]

# per-condition source columns in the wide CSV. None => metric absent for that condition.
_PER_COND = {
    "af3_msa_off": {"iptm": "af3_msaoff_iptm", "ptm": "af3_msaoff_ptm",
                    "has": "has_af3_msaoff", "dockq": "dockq_af3_msa_off"},
    "b2_msa_off":  {"iptm": "b2_iptm", "ptm": None,
                    "has": "has_b2", "dockq": "dockq_b2_msa_off"},
    "esm_msa_off": {"iptm": "esm_iptm", "ptm": None, "has": "has_esm", "dockq": "dockq_esm_msa_off"},
}

# design-level structural metrics (same value for both refold rows of a design).
# Each maps to a tuple of candidate columns tried in order; first that exists AND
# has data wins (DSSP sse_* preferred, plain fraction as fallback on subsets where
# sse_* isn't populated).
_DESIGN_LEVEL = {
    "protein_rmsd": ("ca_rmsd_vs_designed",),
    "n_hbonds": ("n_h_bonds_interface",),
    "helix_frac": ("sse_helix_fraction", "helix_fraction"),
    "sheet_frac": ("sse_strand_fraction", "sheet_fraction"),
}


def _canon_stack(s: str) -> str:
    return cfg.STACK_ALIASES.get(s, s)


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "1.0", "yes"})


def _read_list(v):
    """Parse a stored list-like ('[1, 4, 9]') into a python list; [] on failure."""
    if isinstance(v, (list, tuple)):
        return list(v)
    if not isinstance(v, str) or not v.strip():
        return []
    try:
        out = ast.literal_eval(v)
        return list(out) if isinstance(out, (list, tuple)) else []
    except (ValueError, SyntaxError):
        return []


def _load_sidecar(path: Path) -> pd.DataFrame | None:
    path = Path(path)
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, low_memory=False)
        return df if len(df) else None
    except Exception:  # noqa: BLE001 — a broken sidecar must not sink the whole load
        return None


def load(master_csv=None, stacks=None) -> pd.DataFrame:
    """Return the long per-(design × refold_condition) DataFrame with SCHEMA columns."""
    master_csv = Path(master_csv or cfg.MASTER_CSV)
    w = pd.read_csv(master_csv, low_memory=False)
    w["stack"] = w["stack"].map(_canon_stack)
    w = w[w["stack"].isin(cfg.STACK_ORDER)].copy()
    if stacks:
        keep = {_canon_stack(s) for s in stacks}
        w = w[w["stack"].isin(keep)].copy()

    # derive design-level interface contacts (nt seqids with contact count > 0)
    nt_cols = [f"rna_contact_nt_{i}" for i in range(1, cfg.RNA_NT_COUNT + 1)]
    for c in nt_cols:
        if c not in w.columns:
            w[c] = 0
    nt_mask = w[nt_cols].fillna(0).to_numpy() > 0
    w["interface_contacts"] = [[i + 1 for i in np.nonzero(row)[0]] for row in nt_mask]

    dq = _load_sidecar(cfg.DOCKQ_SIDECAR)
    hb = _load_sidecar(cfg.HBOND_NT_SIDECAR)

    frames = []
    for cond in cfg.REFOLD_CONDITIONS:
        p = _PER_COND[cond]
        sub = w
        if p["has"] in w.columns:
            sub = w[_truthy(w[p["has"]])]
        n = len(sub)
        d = pd.DataFrame({
            "design_id": sub["design_id"].to_numpy(),
            "stack": sub["stack"].to_numpy(),
            "refold_condition": cond,
            "iptm": sub[p["iptm"]].to_numpy() if p["iptm"] in sub.columns else np.nan,
            "ptm": sub[p["ptm"]].to_numpy() if (p["ptm"] and p["ptm"] in sub.columns) else np.nan,
            "pae": np.nan,                 # not computed yet
            "fraction_disordered": np.nan,  # not computed yet
            "n_clashes": np.nan,            # not computed yet
            "interface_contacts": sub["interface_contacts"].tolist(),
        })
        for tgt, srcs in _DESIGN_LEVEL.items():
            col = next((c for c in srcs if c in sub.columns and sub[c].notna().any()), None)
            d[tgt] = sub[col].to_numpy() if col else np.nan

        # DockQ from sidecar (per condition)
        d["dockq"] = np.nan
        if dq is not None and "design_id" in dq.columns and p["dockq"] in dq.columns:
            d = d.merge(dq[["design_id", p["dockq"]]].rename(columns={p["dockq"]: "dockq"}),
                        on="design_id", how="left", suffixes=("", "_sc"))
            if "dockq_sc" in d.columns:
                d["dockq"] = d.pop("dockq_sc")

        # per-nt hbond contacts from sidecar (per design × condition), else []
        d["hbond_contacts"] = [[] for _ in range(len(d))]
        if hb is not None and {"design_id", "refold_condition", "hbond_nts"} <= set(hb.columns):
            hbc = hb[hb["refold_condition"] == cond]
            m = dict(zip(hbc["design_id"], hbc["hbond_nts"].map(_read_list)))
            d["hbond_contacts"] = [m.get(did, []) for did in d["design_id"]]

        frames.append(d)

    long = pd.concat(frames, ignore_index=True)
    for c in SCHEMA:
        if c not in long.columns:
            long[c] = np.nan
    return long[SCHEMA]


def has_metric(df: pd.DataFrame, col: str) -> bool:
    """True if a scalar metric column has any non-NaN value (drives 'pending' panels)."""
    if col not in df.columns:
        return False
    return bool(df[col].notna().any())


def status(df: pd.DataFrame) -> dict:
    """Per-metric availability, for run_all's summary + graceful degradation."""
    scalar = ["iptm", "ptm", "dockq", "protein_rmsd", "pae",
              "fraction_disordered", "n_clashes", "n_hbonds", "helix_frac", "sheet_frac"]
    out = {c: has_metric(df, c) for c in scalar}
    out["interface_contacts"] = bool(df["interface_contacts"].map(len).sum())
    out["hbond_contacts"] = bool(df["hbond_contacts"].map(len).sum())
    return out


def refolder_iptm(master_csv=None, stacks=None) -> pd.DataFrame:
    """Per-design ipTM from each MSA-off refolder (AF3 / Boltz-2 / ESMFold2), long:
    design_id, stack, refolder, iptm. ESM is exposed here only (outside the af3/b2
    filtering scope) so plot 7 can compare which refolder scores highest per stack."""
    w = pd.read_csv(Path(master_csv or cfg.MASTER_CSV), low_memory=False)
    w["stack"] = w["stack"].map(_canon_stack)
    w = w[w["stack"].isin(cfg.STACK_ORDER)]
    if stacks:
        w = w[w["stack"].isin({_canon_stack(s) for s in stacks})]
    frames = []
    for ref, col in cfg.REFOLDER_IPTM_COLS.items():
        if col in w.columns:
            frames.append(pd.DataFrame({"design_id": w["design_id"].to_numpy(),
                                        "stack": w["stack"].to_numpy(),
                                        "refolder": ref, "iptm": w[col].to_numpy()}))
    return (pd.concat(frames, ignore_index=True) if frames
            else pd.DataFrame(columns=["design_id", "stack", "refolder", "iptm"]))


# ---------------------------------------------------------------- thresholds / funnel
# Ordered filter pipeline: (metric, config threshold key, direction). Direction
# "ge" => keep >= threshold, "le" => keep <= threshold.
FILTERS = [
    ("iptm", "iptm_min", "ge"),
    ("dockq", "dockq_min", "ge"),
    ("ptm", "ptm_min", "ge"),
    ("fraction_disordered", "disorder_max", "le"),
    ("n_clashes", "clashes_max", "le"),
    ("n_hbonds", "hbonds_min", "ge"),
]


def available_filters(df):
    """FILTERS whose metric has data AND whose threshold is set (not None)."""
    return [(m, k, d) for (m, k, d) in FILTERS
            if cfg.THRESHOLDS.get(k) is not None and has_metric(df, m)]


def _pass_series(frame, metric, key, direc):
    thr = cfg.THRESHOLDS[key]
    v = frame[metric]
    return (v >= thr) if direc == "ge" else (v <= thr)


def passing_designs(df, mode="avg", filters=None):
    """design_ids passing ALL given (or available) thresholds.

    mode 'avg'    -> the design's mean across refold conditions clears each threshold.
    mode 'strict' -> the design clears each threshold in EVERY refold condition.
    """
    filters = available_filters(df) if filters is None else filters
    if not filters:
        return set(df["design_id"].unique())
    if mode == "avg":
        agg = df.groupby(["design_id"], as_index=False)[[m for m, _, _ in filters]].mean()
        ok = np.ones(len(agg), dtype=bool)
        for m, k, d in filters:
            ok &= _pass_series(agg, m, k, d).fillna(False).to_numpy()
        return set(agg.loc[ok, "design_id"])
    ok_rows = np.ones(len(df), dtype=bool)
    for m, k, d in filters:
        ok_rows &= (_pass_series(df, m, k, d) | df[m].isna()).to_numpy()
    cnt = df.assign(_ok=ok_rows).groupby("design_id")["_ok"].agg(["sum", "count"])
    return set(cnt.index[cnt["sum"] == cnt["count"]])


def _pass_val(frame, metric, thr, direc):
    v = frame[metric]
    return (v >= thr) if direc == "ge" else (v <= thr)


def _funnel_stages(df):
    """Expand cfg.FUNNEL_STAGES into ordered (metric, thr, direc, label) steps, keeping
    only metrics that have data. A STAGGERED metric yields one step per threshold, so the
    funnel tightens in stages (e.g. iptm ≥0.5 → ≥0.65 → ≥0.8)."""
    out = []
    for metric, thrs, direc in cfg.FUNNEL_STAGES:
        if not has_metric(df, metric):
            continue
        sym = "≥" if direc == "ge" else "≤"
        for thr in thrs:
            out.append((metric, thr, direc, f"{metric} {sym} {thr:g}"))
    return out


def _passing_staged(df, mode, stages):
    """design_ids clearing ALL given (metric, thr, direc, label) stages."""
    if not stages:
        return set(df["design_id"].unique())
    metrics = list({m for m, _, _, _ in stages})
    if mode == "avg":
        agg = df.groupby("design_id", as_index=False)[metrics].mean()
        ok = np.ones(len(agg), dtype=bool)
        for m, thr, d, _ in stages:
            ok &= _pass_val(agg, m, thr, d).fillna(False).to_numpy()
        return set(agg.loc[ok, "design_id"])
    ok_rows = np.ones(len(df), dtype=bool)
    for m, thr, d, _ in stages:
        ok_rows &= (_pass_val(df, m, thr, d) | df[m].isna()).to_numpy()
    cnt = df.assign(_ok=ok_rows).groupby("design_id")["_ok"].agg(["sum", "count"])
    return set(cnt.index[cnt["sum"] == cnt["count"]])


def funnel(df, mode="avg"):
    """STAGGERED cumulative filtering funnel (cfg.FUNNEL_STAGES) by absolute design
    counts. Returns (stage_labels, {stack: [count per stage]}, [total per stage], stages)."""
    stages = _funnel_stages(df)
    work = df
    if mode == "strict":
        # Equal-strict with the quality cube (plot_15): only designs refolded by ALL
        # conditions are eligible, and each must clear every threshold in every condition.
        need = set(cfg.REFOLD_CONDITIONS)
        by_design = df.groupby("design_id")["refold_condition"].agg(lambda s: need <= set(s))
        work = df[df["design_id"].isin(set(by_design.index[by_design]))]
    uniq = work.drop_duplicates("design_id")
    stacks = [s for s in cfg.STACK_ORDER if s in set(uniq["stack"])]
    ids_by_stack = {s: set(uniq.loc[uniq["stack"] == s, "design_id"]) for s in stacks}
    labels = ["Total (start)"]
    counts = {s: [len(ids_by_stack[s])] for s in stacks}
    cum = []
    for st in stages:
        cum.append(st)
        keep = _passing_staged(work, mode, cum)
        labels.append(st[3])
        for s in stacks:
            counts[s].append(len(ids_by_stack[s] & keep))
    totals = [sum(counts[s][i] for s in stacks) for i in range(len(labels))]
    return labels, counts, totals, stages


def ranked_wide(master_csv=None, stacks=None):
    """Per-design WIDE table ranked by consensus ipTM = min(AF3, B2, ESM), descending — a design
    ranks high only if ALL THREE refolders agree. Only designs with all three refolds are
    eligible. Carries per-refolder ipTMs, DockQ (merged from the sidecar), length, and the
    designed/refold CIF paths — enough for the ranking tables AND best_binders extraction. Adds
    a 1-based 'rank' + 'rank_in_stack'."""
    w = pd.read_csv(Path(master_csv or cfg.MASTER_CSV), low_memory=False)
    w["stack"] = w["stack"].map(_canon_stack)
    w = w[w["stack"].isin(cfg.STACK_ORDER)].copy()
    if stacks:
        w = w[w["stack"].isin({_canon_stack(s) for s in stacks})].copy()
    # consensus = MIN across all available refolders (af3/b2/esm); every refolder must be high.
    # Recomputed here (overriding any min-of-2 consensus baked into the CSV).
    icols = [cfg.REFOLDER_IPTM_COLS[r] for r in ("af3", "b2", "esm")
             if cfg.REFOLDER_IPTM_COLS[r] in w.columns]
    have = w[icols].notna().all(axis=1)
    w["consensus_iptm"] = w[icols].min(axis=1).where(have)
    w["mean_iptm"] = w[icols].mean(axis=1).where(have)
    w = w[w[cfg.RANK_COL].notna()].copy()
    dq = _load_sidecar(cfg.DOCKQ_SIDECAR)
    if dq is not None and "design_id" in dq.columns:
        keep = ["design_id"] + [c for c in dq.columns if c.startswith("dockq_")]
        w = w.merge(dq[keep].drop_duplicates("design_id"), on="design_id", how="left")
    w = w.sort_values(cfg.RANK_COL, ascending=False, kind="stable").reset_index(drop=True)
    w.insert(0, "rank", range(1, len(w) + 1))
    w["rank_in_stack"] = w.groupby("stack").cumcount() + 1
    return w


def cif_paths(master_csv=None, stacks=None):
    """Per-design refold CIF refs (archive!member, RB-relative) for click-to-open scatters:
    design_id + af3_msaoff_cif_path / b2_cif_path / esm_cif_path. Missing cols -> absent."""
    w = pd.read_csv(Path(master_csv or cfg.MASTER_CSV), low_memory=False)
    if "design_id" not in w.columns:
        return pd.DataFrame(columns=["design_id"])
    cols = ["design_id", "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]
    cols = [c for c in cols if c in w.columns]
    return w[cols].drop_duplicates("design_id")


TIP_COL = "tip_weighted_score"   # the tip-vs-stem weighted contact score (compute_metrics)


def quality_axes(master_csv=None, stacks=None, require=("iptm", "dockq", "tip")):
    """Per-design ipTM / DockQ / tip for the 2-D & 3-D quality plots (+ click refs).
    iptm = consensus min(af3,b2,esm); dockq = consensus = min across the available dockq_*
    sidecar columns; tip = tip_weighted_score. Keeps only rows with every metric in `require`
    present (default all three; plot_15 passes ('iptm','dockq') so tip shows but doesn't gate)."""
    w = ranked_wide(master_csv=master_csv, stacks=stacks)   # consensus_iptm + dockq_* + tip_* + cif paths
    dqcols = [c for c in w.columns if c.startswith("dockq_")]
    w = w.copy()
    w["dockq"] = w[dqcols].min(axis=1) if dqcols else np.nan
    w["tip"] = pd.to_numeric(w[TIP_COL], errors="coerce") if TIP_COL in w.columns else np.nan
    w = w.rename(columns={"consensus_iptm": "iptm"})
    cols = ["design_id", "stack", "iptm", "dockq", "tip",
            "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]
    cols = [c for c in cols if c in w.columns]
    return w[cols].dropna(subset=list(require)).reset_index(drop=True)


def scaffold_fit(master_csv=None, stacks=None):
    """Per-design scaffold-fit metrics + click refs, from the scaffold_fit sidecar merged on
    design_id: does the binder dock a csy4 hairpin WITHIN the full 121-nt construct with its
    interface on the hairpin (scaffold_interface_frac) and no clashes with the rest
    (scaffold_clashes)? Empty if the sidecar is absent."""
    w = ranked_wide(master_csv=master_csv, stacks=stacks)
    sf = _load_sidecar(cfg.SCAFFOLD_FIT_SIDECAR)
    num = ["scaffold_clashes", "scaffold_interface_frac", "scaffold_fits",
           "scaffold_best_hairpin", "scaffold_align_rmsd", "scaffold_n_hairpin_contacts"]
    if sf is None or "design_id" not in sf.columns:
        return pd.DataFrame(columns=["design_id", "stack"] + num)
    keep = ["design_id"] + [c for c in num if c in sf.columns]
    m = w.merge(sf[keep], on="design_id", how="inner")
    for c in num:
        if c in m.columns:
            m[c] = pd.to_numeric(m[c], errors="coerce")
    cols = (["design_id", "stack"] + [c for c in num if c in m.columns]
            + ["af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"])
    cols = [c for c in cols if c in m.columns]
    return m[cols].dropna(subset=["scaffold_clashes", "scaffold_interface_frac"]).reset_index(drop=True)


_TRIPLET_COLS = {
    "iptm": {"af3": "af3_msaoff_iptm", "b2": "b2_iptm", "esm": "esm_iptm"},
    "dockq": {"af3": "dockq_af3_msa_off", "b2": "dockq_b2_msa_off", "esm": "dockq_esm_msa_off"},
}


def refolder_triplet(kind, master_csv=None, stacks=None):
    """Per-design (af3, b2, esm) values of `kind` ('iptm' or 'dockq') + click refs, for the
    per-refolder 3-D cubes. Only designs with all three refolders present. Empty if the metric
    isn't available for all three (e.g. no ESM DockQ column)."""
    src = _TRIPLET_COLS[kind]
    w = ranked_wide(master_csv=master_csv, stacks=stacks)   # carries iptm cols + dockq_* + cif paths
    if not all(c in w.columns for c in src.values()):
        return pd.DataFrame(columns=["design_id", "stack", "af3", "b2", "esm"])
    w = w.rename(columns={src["af3"]: "af3", src["b2"]: "b2", src["esm"]: "esm"})
    cols = ["design_id", "stack", "af3", "b2", "esm",
            "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path"]
    cols = [c for c in cols if c in w.columns]
    return w[cols].dropna(subset=["af3", "b2", "esm"]).reset_index(drop=True)


if __name__ == "__main__":
    d = load()
    print(f"rows={len(d)}  designs={d['design_id'].nunique()}  stacks={d['stack'].nunique()}")
    print("per-condition counts:")
    print(d.groupby("refold_condition").size())
    print("metric availability:")
    for k, v in status(d).items():
        print(f"  {'OK ' if v else '-- '} {k}")
