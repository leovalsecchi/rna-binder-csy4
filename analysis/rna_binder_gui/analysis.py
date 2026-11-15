"""Data layer for the Analysis tab.

Loads an aggregated-metrics CSV (the long, one-row-per-refold-sample output of
`aggregate_metrics.py --structural`), reshapes it to one row per design, and
applies threshold filters. The reshape + composite score mirror the canonical
`analyses/.../full_analysis.py` so the GUI matches the offline analysis.

Pure pandas — no Streamlit import, so it can be unit-tested standalone.

Long CSV: rows carry a `refolder` column (AF3/B2); interface iPTM lives in
`af3_protein_rna_iptm` / `b2_protein_rna_iptm`; structural metrics
(`ca_rmsd_vs_designed`, `interface_bsa_approx`, `n_interface_contacts`,
`n_h_bonds_interface`, ...) are refolder-agnostic (same column either way).
Per-design table: for each (stack, design_key) the best refold sample by iPTM is
kept, separately per refolder, giving `af3_*` and `b2_*` columns.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

# Refolder-agnostic structural columns carried onto each best sample.
STRUCT_COLS = [
    "n_interface_contacts", "n_h_bonds_interface", "interface_bsa_approx",
    "radius_of_gyration_ca", "helix_fraction", "sheet_fraction",
    "loop_fraction", "ca_rmsd_vs_designed", "dockq", "dockq_lrmsd",
]

# Source column -> short per-design suffix (prefixed with af3_/b2_).
_STRUCT_SUFFIX = {
    "n_interface_contacts": "contacts", "n_h_bonds_interface": "hbonds",
    "interface_bsa_approx": "bsa", "ca_rmsd_vs_designed": "rmsd",
    "dockq": "dockq", "dockq_lrmsd": "lrmsd",
    "radius_of_gyration_ca": "rog", "helix_fraction": "helix",
    "sheet_fraction": "sheet", "loop_fraction": "loop",
}

# Identity columns that, together with `stack`, define one design.
_KEY_COMPONENTS = ["design_id", "seq_id", "cycle", "hd_sample", "pack"]


def find_metric_csvs(root: str) -> list[str]:
    """Structural-metrics CSVs under <root>/analyses, newest first.

    Bounded-depth globs only (analyses/<x>/ and analyses/<x>/<y>/) — a recursive
    rglob over the sshfs mount is slow and can wedge it, so we avoid it.
    """
    base = Path(root).expanduser() / "analyses"
    if not base.is_dir():
        return []
    hits: dict[str, float] = {}
    for pattern in ("*structural*.csv", "*/*structural*.csv", "*/*/*structural*.csv"):
        try:
            for p in base.glob(pattern):
                if not p.name.startswith("."):
                    try:
                        hits[str(p)] = p.stat().st_mtime
                    except OSError:
                        continue
        except OSError:          # eulermount hiccup
            continue
    return [p for p, _ in sorted(hits.items(), key=lambda kv: kv[1], reverse=True)]


def _design_key(row) -> str:
    parts = [str(row["stack"])]
    for col in _KEY_COMPONENTS:
        v = row.get(col)
        if pd.notna(v):
            try:
                v = int(v) if float(v).is_integer() else v
            except (TypeError, ValueError):
                pass
            parts.append(f"{col}={v}")
    return "|".join(parts)


def load_metrics(path: str) -> pd.DataFrame:
    """Load the long per-sample CSV; add `design_key` and a unified `iptm`."""
    df = pd.read_csv(path)
    df["design_key"] = df.apply(_design_key, axis=1)
    af3 = df["af3_protein_rna_iptm"] if "af3_protein_rna_iptm" in df else np.nan
    b2 = df["b2_protein_rna_iptm"] if "b2_protein_rna_iptm" in df else np.nan
    esm = df["esm_protein_rna_iptm"] if "esm_protein_rna_iptm" in df else np.nan
    df["iptm"] = np.where(df["refolder"] == "AF3", af3,
                          np.where(df["refolder"] == "B2", b2, esm))
    return df


def _best_per_refolder(df: pd.DataFrame, refolder: str, prefix: str) -> pd.DataFrame:
    """Best (max-iPTM) refold sample per design for one refolder."""
    sub = df[df["refolder"] == refolder].dropna(subset=["iptm"])
    if sub.empty:
        return pd.DataFrame(columns=["stack", "design_key", f"{prefix}_iptm"])
    # ascending sort + keep last == max iptm per (stack, design_key)
    best = sub.sort_values("iptm").drop_duplicates(["stack", "design_key"], keep="last")
    rename = {"iptm": f"{prefix}_iptm"}
    rename.update({c: f"{prefix}_{_STRUCT_SUFFIX[c]}" for c in STRUCT_COLS if c in best})
    return best[["stack", "design_key", *rename]].rename(columns=rename)


def to_per_design(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the long sample table to one row per design (af3_*/b2_*).

    Also carries `designed_source` (path to the designed file) per design when
    present, so the GUI can open a design's 4 variants in ChimeraX.
    """
    af3 = _best_per_refolder(df, "AF3", "af3")
    b2 = _best_per_refolder(df, "B2", "b2")
    esm = _best_per_refolder(df, "ESM", "esm")
    merged = af3.merge(b2, on=["stack", "design_key"], how="outer")
    merged = merged.merge(esm, on=["stack", "design_key"], how="outer")
    if "designed_source" in df.columns:
        src = df[df["designed_source"].notna() & (df["designed_source"].astype(str) != "")]
        src = src.drop_duplicates(["stack", "design_key"])[["stack", "design_key", "designed_source"]]
        merged = merged.merge(src, on=["stack", "design_key"], how="left")
    return merged


def add_composite(designs: pd.DataFrame, prefix: str = "af3") -> pd.DataFrame:
    """composite = iptm * (1-clip(rmsd,0,3)/3) * clip(bsa/1500,0,1) [* dockq if present]."""
    d = designs.copy()
    need = [f"{prefix}_iptm", f"{prefix}_rmsd", f"{prefix}_bsa"]
    if not all(c in d for c in need):
        d[f"{prefix}_composite"] = np.nan
        return d
    rmsd_factor = (1 - d[f"{prefix}_rmsd"].clip(0, 3) / 3).clip(0, 1)
    bsa_factor = (d[f"{prefix}_bsa"] / 1500).clip(0, 1)
    d[f"{prefix}_composite"] = d[f"{prefix}_iptm"] * rmsd_factor * bsa_factor
    # Fold in DockQ pose quality when the column is present (CSVs post-2026-06-05).
    if f"{prefix}_dockq" in d:
        d[f"{prefix}_composite"] *= d[f"{prefix}_dockq"].fillna(0).clip(0, 1)
    return d


def filter_mask(designs: pd.DataFrame, prefix: str, *,
                iptm_min: float = 0.5, rmsd_max: float = 2.0, bsa_min: float = 800.0,
                contacts_min: int = 10, hbonds_min: int = 0,
                dockq_min: float = 0.23) -> pd.Series:
    """Boolean mask of designs passing every active threshold (NaN never passes).

    Strict inequalities for iptm/bsa/contacts and rmsd, matching full_analysis.py;
    hbonds and dockq use >= (set hbonds_min/dockq_min=0 to disable; dockq applied
    only when the dockq column is present). NOTE: dockq_min gates the DockQ *score*
    (0-1, ≥0.23 = CAPRI acceptable) — not the ligand RMSD.
    """
    m = pd.Series(True, index=designs.index)
    checks = [
        (f"{prefix}_iptm", lambda c: c > iptm_min),
        (f"{prefix}_rmsd", lambda c: c < rmsd_max),
        (f"{prefix}_bsa", lambda c: c > bsa_min),
        (f"{prefix}_contacts", lambda c: c > contacts_min),
    ]
    for col, test in checks:
        if col in designs:
            m &= test(designs[col]).fillna(False)
    if hbonds_min and f"{prefix}_hbonds" in designs:
        m &= (designs[f"{prefix}_hbonds"] >= hbonds_min).fillna(False)
    if dockq_min and f"{prefix}_dockq" in designs:
        m &= (designs[f"{prefix}_dockq"] >= dockq_min).fillna(False)
    return m


def pass_rate_table(designs: pd.DataFrame, prefix: str, **thresholds) -> pd.DataFrame:
    """Per-stack pass counts/rates at the given thresholds (+ overall row)."""
    rows = []
    for stack in sorted(designs["stack"].dropna().unique()):
        g = designs[designs["stack"] == stack]
        passed = int(filter_mask(g, prefix, **thresholds).sum())
        n = len(g)
        rows.append({"stack": stack, "designs": n, "pass": passed,
                     "pass %": round(100 * passed / n, 1) if n else 0.0})
    total = len(designs)
    total_pass = int(filter_mask(designs, prefix, **thresholds).sum())
    rows.append({"stack": "ALL", "designs": total, "pass": total_pass,
                 "pass %": round(100 * total_pass / total, 1) if total else 0.0})
    return pd.DataFrame(rows)
