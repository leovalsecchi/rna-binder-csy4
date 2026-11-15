#!/usr/bin/env python3
"""Metric Explorer — build intuition for what each interface metric *does* to designs.

A local Flask app. Top: controls for every metric. Bottom: three turnable 3D panels —
LEFT = a design just *excluded* by the driven metric's threshold, MIDDLE+RIGHT = two
designs just *included*, chosen so every OTHER metric is near-identical (a controlled
comparison: only the driven metric differs). Next/Back walk the matched sets.

Data: the cross-campaign metrics CSV (design_id, stack, per-refolder metrics + the
interface metrics) at data/metrics.csv. Structures: AF3 refold CIFs streamed on demand
from Euler (tar --zstd -xOf over ssh, exactly like open_refolds_server.py) and cached
locally. "Open in ChimeraX" hands off to the existing openrefolds bridge on :45001.

Run:  python3 app.py         # http://127.0.0.1:8000
"""
from __future__ import annotations

import json
import shlex
import subprocess
import threading
from pathlib import Path

import numpy as np
import pandas as pd
from flask import Flask, request, jsonify, send_file, abort
from flask_cors import CORS

import nas_collections as coll     # curated design "folders" on the NAS (see nas_collections.py)
import local_download as dl        # bulk-copy the filtered set to a local folder (local_download.py)
import presets as pre              # named filter presets in presets.json (see presets.py)

APP_DIR = Path(__file__).resolve().parent
DATA_CSV = APP_DIR / "data" / "metrics.csv"
CACHE = APP_DIR / "cache" / "cifs"
CACHE.mkdir(parents=True, exist_ok=True)

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SSH = ["ssh", "-o", "ConnectTimeout=20", "-o", "BatchMode=yes", "euler"]
REFOLD_COL = {"af3": "af3_msaoff_cif_path", "b2": "b2_cif_path", "esm": "esm_cif_path"}
BRIDGE = "http://127.0.0.1:45001"      # open_refolds_server.py

# Metric registry: (key, label, column, higher_is_better). Only those whose column is
# present in the CSV are exposed; consensus_iptm / dockq_consensus are derived on load.
# (key, label, column, higher_is_better, description)
METRICS_SPEC = [
    ("iptm",         "ipTM (consensus)",       "consensus_iptm",       True,
     "Interface predicted TM-score — the folding model's confidence in the binder–RNA "
     "interface pose (0–1), min across AF3 / Boltz-2 / ESMFold2. Higher = more confident, "
     "more agreed interface."),
    ("dockq",        "DockQ (consensus)",      "dockq_consensus",      True,
     "DockQ docking-quality of the predicted complex (0–1), min across refolders. "
     "Higher = the refolders agree more on the docked pose."),
    ("pistack",      "π-stacking (min agreed)","pi_stack_score",       True,
     "Aromatic π-stacking (Trp/Tyr/Phe/His + planar Arg ↔ RNA base), taken as the MIN "
     "across the AF3 / Boltz-2 / ESMFold2 refolds — only stacks all three backends agree "
     "on. Sparse (~0.1%): filter with >0 to isolate the robust stackers. A descriptor, "
     "not a quality ranker (magnitude doesn't track DockQ/ipTM here)."),
    ("dsasa_total",  "ΔSASA total (Å²)",       "dsasa_total",          True,
     "Total solvent-accessible surface area buried when binder and RNA associate "
     "(Rosetta InterfaceAnalyzer). Bigger = larger interface / more contact area. "
     "MIN-consensus across AF3/B2/ESM (worst-of-3); AF3-only value in dsasa_total_af3."),
    ("dsasa_frac",   "Hydrophobic ΔSASA frac", "dsasa_hphobic_frac",   True,
     "Fraction of the buried interface that is hydrophobic (0–1). High = a "
     "hydrophobically-driven interface; low = electrostatic / polar (typical for RNA)."),
    ("sc",           "Shape complementarity",  "sc_value",             True,
     "Rosetta shape complementarity of the interface (0–1) — how snugly the two "
     "surfaces fit together. Higher = tighter geometric fit. "
     "MIN-consensus across AF3/B2/ESM (worst-of-3); AF3-only value in sc_value_af3."),
    ("clashes",      "Interface clashes",      "n_interface_clashes",  False,
     "Number of protein–RNA heavy-atom pairs whose van der Waals spheres overlap by "
     "> 0.4 Å — steric clashes at the interface. Lower = fewer clashes; 0 is ideal. "
     "MAX-consensus across AF3/B2/ESM (worst-of-3 — any refolder's clashes count); "
     "per-backend in n_interface_clashes_af3/_b2/_esm."),
    ("hbonds",       "Interface H-bonds",      "n_h_bonds_interface",  True,
     "Number of hydrogen bonds across the protein–RNA interface (PyRosetta, "
     "hydrogen-aware and angle-checked). MIN-consensus across the AF3/B2/ESM refolds "
     "(fewest of the 3 — a design counts as H-bonding only if all three refolders agree, "
     "so a single-backend H-bond doesn't inflate it). Naturals make ~14 vs designs ~2; "
     "designs badly under-H-bond the RNA. Per-backend in n_h_bonds_interface_af3/_b2/_esm. "
     "Higher = better."),
    ("contacts",     "Heavy-atom contacts",    "n_heavy_contacts_4A",  True,
     "Count of protein–RNA heavy-atom pairs within 4 Å — how much atomic contact / "
     "packing there is at the interface. MIN-consensus across AF3/B2/ESM (worst-of-3); "
     "AF3-only value in n_heavy_contacts_4A_af3."),
    ("frac_base",    "Fraction base contacts", "frac_base_contacts",   True,
     "Share of interface contacts made to the RNA bases vs the phosphate/sugar "
     "backbone. High = reads the bases (sequence-specific); low = grips the backbone. "
     "MIN-consensus across AF3/B2/ESM (worst-of-3); AF3-only in frac_base_contacts_af3. "
     "Blank where a refold makes no base-or-backbone contacts."),
    ("tip",          "Tip-weighted score",     "tip_weighted_score",   True,
     "Geometric score favouring contacts near the RNA hairpin loop tip over the "
     "scaffold-adjacent stem bottom. Higher = the binder grips the tip."),
    ("rog",          "Radius of gyration (Å)", "radius_of_gyration_ca",False,
     "Compactness of the binder's Cα cloud (Å). Lower = a more compact fold."),
    # --- unfavorable ("red") interface interactions (ChimeraX RNAcon port; lower = better).
    # Rates per interface residue; the natural RBP baseline is given in each description.
    ("acid_phos",    "Acid↔phosphate repulsion","red_acid_phos_per_ifaceres", False,
     "Asp/Glu carboxylate O within 4 Å of an RNA phosphate O, per interface residue — "
     "like-charge (negative-on-negative) repulsion. MAX-consensus across the AF3/B2/ESM "
     "refolds (worst of the 3 — a design is clean only if all three agree). Natural RBPs "
     "≈ 0.007; RFD3 and HD are natural-clean, BoltzGen runs 4–5× higher. The cleanest "
     "cross-set red signal — chemistry-driven, robust to the predicted-vs-crystal confound. "
     "Per-backend in red_acid_phos_per_ifaceres_af3/_b2/_esm. Lower = better."),
    ("red_clash",    "Clashes / iface residue","red_clashes_per_ifaceres",   False,
     "Protein–RNA van der Waals overlaps (≥ 0.6 Å) per interface residue "
     "(ChimeraX-calibrated). MAX-consensus across AF3/B2/ESM refolds (a design counts as "
     "clashy if ANY refolder flags it → robust to single-backend artifacts). Natural ≈ "
     "0.000. Large design–natural gap, but CONFOUNDED (predicted models vs clash-free "
     "crystals) — use to rank stacks, not as an absolute score. Per-backend in "
     "red_clashes_per_ifaceres_af3/_b2/_esm. Lower = better."),
    ("red_total",    "Red interactions total", "red_total_per_ifaceres",     False,
     "All unfavorable interface flags (clashes + acid↔phosphate + hydrophobic-C-on-polar-"
     "backbone) per interface residue, MAX-consensus across AF3/B2/ESM refolds. Natural ≈ "
     "0.069. Lower = fewer bad contacts. Blank only for designs undocked in ALL three "
     "refolds (no interface). Per-backend in red_total_per_ifaceres_af3/_b2/_esm."),
    # --- sequence low-complexity (from the binder sequence; lower = better) ---
    ("seq_max_aa",   "Max single-AA fraction", "seq_max_aa_frac",            False,
     "Largest share of the binder sequence taken by a single amino acid (0–1). High = "
     "poly-X / a LigandMPNN sequence collapse. Median 0.20 overall, but RFD3+LiMPNN runs "
     "~0.35. Lower = more complex sequence."),
    ("seq_run",      "Longest same-AA run",    "seq_longest_run",            False,
     "Longest homopolymer stretch (same amino acid in a row) in the binder. Median 3; a "
     "run ≥6 (~19% of designs) flags a repetitive/collapsed sequence. Lower = better."),
    ("seq_kmer",     "Longest repeated k-mer", "seq_max_kmer",               False,
     "Length of the longest substring (k-mer) that repeats within the binder sequence "
     "(occurrences may overlap). Catches PERIODIC / dispersed multi-residue repeats "
     "(e.g. ELELEL, KTKTKT) that the single-AA metrics miss. A random ~80-mer sits ~3–4; "
     "much longer = a repetitive design. Filter out designs whose largest repeated k-mer "
     "exceeds your k. Lower = better."),
    # --- binder core packing (ChimeraX Bbindercore/rnacon port; FOLD quality of the binder
    # itself — its buried monomer core — NOT the RNA interface). SASA on the binder alone. ---
    ("core_hyd",     "Core hydrophobic frac",  "core_hydrophobic_frac",      True,
     "Fraction of the folded binder's BURIED core (residues with monomer SASA < 30 Å²) that is "
     "hydrophobic {A V L I M F W Y P C} — the ChimeraX 'yellow core' as a number. High = a clean "
     "hydrophobic core; low = polar/charged residues buried in the fold (destabilising). A "
     "fold-quality signal independent of the RNA interface."),
    ("core_pack",    "Core packing quality",   "core_packing_score",         True,
     "Geometric core-packing quality (0–1): mean over buried-core residues of how well each is "
     "coordinated — Cβ neighbours within 8 Å, mapped between under-packed (8) and well-packed "
     "(15). High = a densely, void-free packed core; low = loose / under-packed (packing "
     "defects). Complements the hydrophobic fraction — it measures whether the core actually "
     "PACKS, not just its chemistry. Raw coordination is in core_mean_cb_cn."),
    ("core_underpk", "Core under-packed frac", "core_underpacked_frac",      False,
     "Fraction of buried-core residues that are under-coordinated (< 9 Cβ neighbours within "
     "8 Å) — a void / packing-defect proxy. Lower = better packed."),
    # --- apo / free-monomer self-consistency (fold the binder ALONE, no RNA; lower = better) ---
    ("apo_rmsd",     "Apo self-consistency (Å)","apo_rmsd_vs_esm",            False,
     "Fold the binder ALONE with ESMFold2 (no RNA), superpose on the binder chain from the "
     "ESM complex refold → Cα RMSD. Low = the binder holds its designed fold WITHOUT the "
     "target (pre-organised / self-consistent); high = it only takes shape when the RNA is "
     "present. A FOLD-quality signal — distinct from the complex-based ca_rmsd and from "
     "radius_of_gyration; pairs with core-packing. Overall median ~3.9 Å, ~36% fold within "
     "2 Å; BG stacks are the most self-consistent (~70% < 2 Å), HD the least (< 1%). "
     "Lower = better."),
]

_LOCK = threading.Lock()
_DATASETS: dict = {}       # name -> {df, metrics, Z, cols, setcache}

# selectable datasets (name -> CSV); first existing one is the default
DATASET_FILES = [("100k", APP_DIR / "data" / "metrics_100k.csv"),
                 ("100k-2of3", APP_DIR / "data" / "metrics_100k_2of3.csv"),   # median-of-3 (less strict)
                 ("30k",  APP_DIR / "data" / "metrics_30k.csv")]
DEFAULT_DATASET = "100k"


# --------------------------------------------------------------------------- data
def load_one(path):
    df = pd.read_csv(path, low_memory=False)

    # derive consensus (min across refolders) like the plots, if per-refolder present
    def consensus(cols):
        have = [c for c in cols if c in df.columns]
        return df[have].apply(pd.to_numeric, errors="coerce").min(axis=1) if have else None

    ci = consensus(["af3_msaoff_iptm", "b2_iptm", "esm_iptm"])
    if ci is not None and "consensus_iptm" not in df.columns:
        df["consensus_iptm"] = ci
    if "dockq_consensus" not in df.columns:              # respect a stored consensus (e.g. 2of3 median)
        dq = consensus(["dockq_af3_msa_off", "dockq_b2_msa_off", "dockq_esm_msa_off",
                        "af3_dockq", "b2_dockq", "esm_dockq",
                        "af3_msaoff_dockq", "dockq_af3", "dockq_b2", "dockq_esm"])
        if dq is not None:
            df["dockq_consensus"] = dq

    metrics = []
    for key, label, col, hib, desc in METRICS_SPEC:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            if df[col].notna().sum() >= 50:          # enough data to be useful
                metrics.append({"key": key, "label": label, "col": col,
                                "higher_is_better": hib, "desc": desc})
    if not metrics:
        raise SystemExit("no known metric columns found in " + str(path))

    for m in metrics:                                # per-metric stats for sliders
        s = df[m["col"]].dropna()
        m.update(min=float(s.min()), max=float(s.max()), median=float(s.median()),
                 p10=float(s.quantile(.10)), p90=float(s.quantile(.90)))

    cols = [m["col"] for m in metrics]
    Z = df[cols].apply(lambda c: (c - c.mean()) / (c.std(ddof=0) or 1.0)).fillna(0.0)
    return {"df": df, "metrics": metrics, "Z": Z.to_numpy(), "cols": cols, "setcache": {}}


def load_all():
    for name, path in DATASET_FILES:
        if Path(path).exists():
            _DATASETS[name] = load_one(path)
            st = _DATASETS[name]
            print(f"[{name}] loaded {len(st['df']):,} designs, {len(st['metrics'])} metrics: "
                  f"{', '.join(m['key'] for m in st['metrics'])}")
    if not _DATASETS:
        raise SystemExit("no dataset CSVs found in data/")


def _ds(name):
    return (_DATASETS.get(name) or _DATASETS.get(DEFAULT_DATASET)
            or next(iter(_DATASETS.values())))


_OPS = {"gt", "ge", "lt", "le", "eq"}


def _parse_constraints(st, s, skip_key=None):
    """Parse "key:op:value,..." (op in gt/ge/lt/le/eq; legacy "key:value" uses the
    metric's natural direction) into a tuple of (column, op, value). Skips skip_key."""
    by_key = {x["key"]: x for x in st["metrics"]}
    out = []
    for part in (s or "").split(","):
        bits = part.split(":")
        if len(bits) == 3:
            k, op, v = bits
        elif len(bits) == 2:
            k, v, op = bits[0], bits[1], None
        else:
            continue
        mm = by_key.get(k)
        if not mm or k == skip_key:
            continue
        try:
            val = float(v)
        except ValueError:
            continue
        if op is None:
            op = "ge" if mm["higher_is_better"] else "le"
        if op in _OPS:
            out.append((mm["col"], op, val))
    return tuple(out)


def _apply_constraints(df, valid, constraints):
    """AND the hard filters onto a boolean `valid` mask (NaN compares False -> excluded)."""
    for ccol, cdir, cval in constraints:
        cv = pd.to_numeric(df[ccol], errors="coerce").to_numpy()
        if   cdir == "ge": mask = cv >= cval
        elif cdir == "gt": mask = cv >  cval
        elif cdir == "le": mask = cv <= cval
        elif cdir == "lt": mask = cv <  cval
        elif cdir == "eq": mask = cv == cval
        else:              continue
        valid = valid & mask
    return valid


def _selected_idx(st, constrain, metric):
    """Row indices of the designs passing the hard filters (each must have an AF3 refold),
    ordered by `metric` (NaN last). This is THE definition of the 'selected' set — shared by
    /api/browse and /api/download/start so the download copies exactly what Browse shows."""
    df = st["df"]
    constraints = _parse_constraints(st, constrain or "", skip_key=None)
    valid = df[REFOLD_COL["af3"]].notna().to_numpy()
    valid = _apply_constraints(df, valid, constraints)
    idx = np.where(valid)[0]
    m = next((x for x in st["metrics"] if x["key"] == metric), None)
    if m is not None and len(idx):
        vals = pd.to_numeric(df[m["col"]], errors="coerce").to_numpy()[idx]
        fill = -np.inf if m["higher_is_better"] else np.inf   # push NaN to the end
        order = np.argsort(np.where(np.isnan(vals), fill, vals))
        if m["higher_is_better"]:
            order = order[::-1]
        idx = idx[order]
    return idx


# ----------------------------------------------------------------- match selection
def matched_triples(st, driven_key, center, delta, ignore_keys=(), constraints=(),
                    max_sets=60, tol_frac=0.3):
    """Ranked below→center→above comparisons for the driven metric.

    Middle design sits at the driven metric = `center`; left ≈ center-delta, right ≈
    center+delta. All three are matched (nearest neighbours) in the z-scored space of
    every OTHER metric — EXCLUDING the driven one and any in `ignore_keys` — so the
    driven metric is (ideally) the only thing that changes across the row. For each
    middle anchor (nearest to `center`) we grab the nearest below-window and
    above-window design; triples are ranked by other-metric tightness. Returns a list
    of {below, at, above, tightness, m_below, m_at, m_above}.
    """
    from scipy.spatial import cKDTree
    df, metrics, Z, cols = st["df"], st["metrics"], st["Z"], st["cols"]
    m = next(x for x in metrics if x["key"] == driven_key)
    col = m["col"]
    di = cols.index(col)

    valid = (df[col].notna() & df[REFOLD_COL["af3"]].notna()).to_numpy()
    valid = _apply_constraints(df, valid, constraints)       # hard filters
    vidx = np.where(valid)[0]
    vals = df[col].to_numpy()

    # feature space = all metrics EXCEPT the driven one AND any user-ignored ones
    drop = {di}
    for k in ignore_keys:
        mm = next((x for x in metrics if x["key"] == k), None)
        if mm:
            drop.add(cols.index(mm["col"]))
    keep = [j for j in range(Z.shape[1]) if j not in drop]
    Zf = Z[:, keep] if keep else np.zeros((Z.shape[0], 1))   # all ignored -> no constraint

    spread = float(np.nanpercentile(vals[vidx], 90) - np.nanpercentile(vals[vidx], 10)) or 1.0
    tol = max(delta * tol_frac, 0.03 * spread)              # half-width of each target window

    def pool(target):
        return vidx[np.abs(vals[vidx] - target) <= tol]

    mid_pool, lo_pool, hi_pool = pool(center), pool(center - delta), pool(center + delta)
    if len(mid_pool) == 0 or len(lo_pool) == 0 or len(hi_pool) == 0:
        return []

    lo_tree, hi_tree = cKDTree(Zf[lo_pool]), cKDTree(Zf[hi_pool])
    anchors = mid_pool[np.argsort(np.abs(vals[mid_pool] - center))[:max(5, min(len(mid_pool), 200))]]
    sets, seen = [], set()
    for a in anchors:
        dl, li = lo_tree.query(Zf[a]); dh, hi = hi_tree.query(Zf[a])
        lo_i, hi_i = int(lo_pool[li]), int(hi_pool[hi])
        key = (lo_i, int(a), hi_i)
        if lo_i == hi_i or key in seen:
            continue
        seen.add(key)
        sets.append({"below": lo_i, "at": int(a), "above": hi_i,
                     "tightness": float(dl + dh)})
    sets.sort(key=lambda s: s["tightness"])
    return sets[:max_sets]


def _row_payload(st, i):
    df, metrics = st["df"], st["metrics"]
    r = df.iloc[i]
    vals = {}
    for m in metrics:
        v = r[m["col"]]
        vals[m["key"]] = None if pd.isna(v) else round(float(v), 4)
    return {
        "row": int(i),
        "design_id": str(r.get("design_id", "")),
        "stack": str(r.get("stack", "")),
        "metrics": vals,
        "refs": {k: (None if pd.isna(r.get(c)) else str(r.get(c)))
                 for k, c in REFOLD_COL.items()},
    }


def _s(v):
    """NaN/None -> '' , else str (for CSV-sourced fields that may be missing)."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    return str(v)


def _payload_from_record(st, rec):
    """Shape a raw metrics.csv record (dict) from a NAS folder into a browse-mode payload.
    The AF3 CIF is served from that folder's local cif/ dir, so `local_cif` is its basename."""
    def num(col):
        try:
            f = float(rec.get(col))
        except (TypeError, ValueError):
            return None
        return None if pd.isna(f) else round(f, 4)
    af3 = rec.get("af3_msaoff_cif_path")
    local = Path(af3.split("!", 1)[-1]).name if isinstance(af3, str) and af3 else None
    return {
        "design_id": _s(rec.get("design_id")),
        "stack": _s(rec.get("stack")),
        "metrics": {m["key"]: num(m["col"]) for m in st["metrics"]},
        "note": _s(rec.get("note")),
        "local_cif": local,
        "refs": {k: (_s(rec.get(c)) or None) for k, c in REFOLD_COL.items()},
    }


# ------------------------------------------------------------------------- CIF I/O
def fetch_cif(ref: str) -> Path:
    """Stream a CIF member from its Euler tar.zst archive (ssh) and cache it. The ref is
    confined to RB (coll.validate_ref) and options are ended with '--', so a crafted
    ?ref=/etc/passwd or an option-like member can't read outside the repo / inject a tar flag."""
    arch, member = coll.validate_ref(ref)         # raises ValueError on absolute/escape/'-'
    key = Path(member).name if member else Path(arch).name
    out = CACHE / key
    if out.exists() and out.stat().st_size > 0:
        return out
    if member:
        remote = (f"cd {shlex.quote(RB)} && tar --zstd -xOf {shlex.quote(arch)} "
                  f"-- {shlex.quote(member)}")
    else:
        remote = f"cd {shlex.quote(RB)} && cat -- {shlex.quote(arch)}"
    with open(out, "wb") as fh:
        p = subprocess.run(SSH + [remote], stdout=fh, stderr=subprocess.PIPE, timeout=180)
    if p.returncode != 0 or out.stat().st_size == 0:
        if out.exists():
            out.unlink()
        raise RuntimeError(p.stderr.decode()[:300] or "empty CIF")
    return out


# ---------------------------------------------------------------------------- HTTP
app = Flask(__name__, static_folder=str(APP_DIR / "static"))
# Local-only tool: allow same-machine origins (any port) but NOT arbitrary websites — a page
# the user visits must not be able to read /api/cif responses cross-origin. Defense in depth on
# top of ref confinement in fetch_cif.
CORS(app, origins=r"^http://(127\.0\.0\.1|localhost)(:\d+)?$")


@app.route("/")
def index():
    return send_file(APP_DIR / "templates" / "index.html")


@app.route("/api/metrics")
def api_metrics():
    st = _ds(request.args.get("dataset"))
    return jsonify({"metrics": st["metrics"], "n": int(len(st["df"])),
                    "datasets": list(_DATASETS.keys()), "default": DEFAULT_DATASET})


@app.route("/api/triplet")
def api_triplet():
    st = _ds(request.args.get("dataset"))
    driven = request.args.get("metric")
    center = float(request.args.get("threshold"))
    delta = float(request.args.get("delta", 0))
    idx = int(request.args.get("index", 0))
    ignore = request.args.get("ignore", "")
    ignore_keys = tuple(sorted(k for k in ignore.split(",") if k))
    # hard filters (the driven metric is the swept axis, so it is never filtered here)
    constraints = _parse_constraints(st, request.args.get("constrain", ""), skip_key=driven)
    with _LOCK:
        ck = (driven, round(center, 6), round(delta, 6), ignore_keys, constraints)
        cache = st["setcache"]
        if ck not in cache:
            cache[ck] = matched_triples(st, driven, center, delta, ignore_keys, constraints)
        sets = cache[ck]
    if not sets:
        return jsonify({"error": "no matched triple", "n_sets": 0})
    idx = idx % len(sets)
    s = sets[idx]
    return jsonify({
        "index": idx, "n_sets": len(sets), "tightness": round(s["tightness"], 3),
        "below": _row_payload(st, s["below"]),
        "at": _row_payload(st, s["at"]),
        "above": _row_payload(st, s["above"]),
        "driven": driven, "center": center, "delta": delta,
    })


@app.route("/api/browse")
def api_browse():
    """Browse mode: list every design passing the hard filters (no matched-triple
    requirement), sorted by the selected metric, paginated 3 per page. Returns fewer
    than 3 designs when that's all that qualifies (never blanks out on a partial page)."""
    st = _ds(request.args.get("dataset"))
    sortkey = request.args.get("metric")
    page = max(0, int(request.args.get("page", 0)))
    idx = _selected_idx(st, request.args.get("constrain", ""), sortkey)
    total = int(len(idx))
    pg = idx[page * 3: page * 3 + 3]
    return jsonify({"total": total, "page": page, "n_pages": (total + 2) // 3,
                    "sort": sortkey, "designs": [_row_payload(st, int(i)) for i in pg]})


@app.route("/api/cif")
def api_cif():
    ref = request.args.get("ref")
    if not ref:
        abort(400)
    try:
        return send_file(fetch_cif(ref), mimetype="chemical/x-cif")
    except ValueError as e:                       # ref failed the RB-confinement check
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 502


@app.route("/api/chimerax", methods=["POST"])
def api_chimerax():
    """Hand off to the openrefolds bridge (open_refolds_server.py on :45001) to open the
    3 refolds in ChimeraX. Bridge API: GET /open?af3=<ref>&b2=<ref>&esm=<ref>&id=<design_id>."""
    import urllib.parse
    import urllib.request
    d = request.get_json(force=True)
    refs = d.get("refs", {}) or {}
    params = {"id": d.get("design_id", "")}
    for k in ("af3", "b2", "esm"):
        if refs.get(k):
            params[k] = refs[k]
    url = BRIDGE + "/open?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            return jsonify({"ok": True, "msg": r.read().decode()[:400]})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 502


# --------------------------------------------------------------------------- presets API
@app.route("/api/presets")
def api_presets():
    return jsonify({"presets": pre.list_presets()})


@app.route("/api/presets/save", methods=["POST"])
def api_presets_save():
    """Upsert a named filter preset. Body: {name, filters:{key:rawtext}, driven, thresh, delta}."""
    d = request.get_json(force=True) or {}
    try:
        items = pre.save_preset(d.get("name", ""), d.get("filters") or {},
                                d.get("driven"), d.get("thresh"), d.get("delta"))
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    except Exception as e:                             # noqa: BLE001 — corrupt presets.json
        return jsonify({"ok": False, "error": f"presets.json unreadable: {e}"}), 500
    return jsonify({"ok": True, "presets": items})


@app.route("/api/presets/delete", methods=["POST"])
def api_presets_delete():
    d = request.get_json(force=True) or {}
    try:
        items = pre.delete_preset(d.get("name", ""))
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    except Exception as e:                             # noqa: BLE001
        return jsonify({"ok": False, "error": f"presets.json unreadable: {e}"}), 500
    return jsonify({"ok": True, "presets": items})


# --------------------------------------------------------------------- download-selected API
@app.route("/api/download/start", methods=["POST"])
def api_download_start():
    """Start a background copy of the CURRENT filtered set into downloads/<name>/. Body:
    {dataset, constrain, metric, name, confirm}. Returns {need_confirm,count} when the set
    exceeds the cap and confirm is not set; otherwise {job_id, total}."""
    d = request.get_json(force=True) or {}
    st = _ds(d.get("dataset"))
    idx = _selected_idx(st, d.get("constrain", ""), d.get("metric"))
    count = int(len(idx))
    if count == 0:
        return jsonify({"ok": False, "error": "no designs match the current filters"}), 400
    if count > dl.CAP and not d.get("confirm"):
        return jsonify({"ok": True, "need_confirm": True, "count": count})
    df = st["df"]
    designs = []
    for i in idx:
        r = df.iloc[int(i)]
        designs.append({
            "design_id": str(r.get("design_id", "")),
            "refs": {k: (None if pd.isna(r.get(c)) else str(r.get(c))) for k, c in REFOLD_COL.items()},
            "row": r.to_dict(),
        })
    try:
        job_id = dl.start(d.get("name", ""), designs)
    except ValueError as e:                              # bad/empty folder name
        return jsonify({"ok": False, "error": str(e)}), 400
    return jsonify({"ok": True, "job_id": job_id, "total": count})


@app.route("/api/download/status")
def api_download_status():
    """Poll a download job's progress: {status, done, total, cifs, missing, folder, name}."""
    j = dl.job_status(request.args.get("job", ""))
    if j is None:
        return jsonify({"error": "unknown job"}), 404
    return jsonify(j)


# --------------------------------------------------------------------- collections API
@app.route("/api/collections")
def api_collections():
    """List folders (+count); if design_id is given, flag whether it's already in each."""
    did = request.args.get("design_id", "")
    folders = coll.list_folders()
    if did:
        for f in folders:
            f["has"] = did in coll.folder_members(f["name"])
    return jsonify({"available": coll.nas_available(), "folders": folders})


@app.route("/api/collections/add", methods=["POST"])
def api_collections_add():
    """Add the design at (dataset,row) to each named folder: extract its 3 refold CIFs +
    upsert its full metric row. Reconstructs everything from the dataframe (no 230-col POST)."""
    d = request.get_json(force=True) or {}
    st = _ds(d.get("dataset"))
    try:
        i = int(d.get("row"))
    except (TypeError, ValueError):
        return jsonify({"ok": False, "error": "bad row"}), 400
    if i < 0 or i >= len(st["df"]):
        return jsonify({"ok": False, "error": "row out of range"}), 400
    folders = [str(f) for f in (d.get("folders") or []) if str(f).strip()]
    if not folders:
        return jsonify({"ok": False, "error": "no folders given"}), 400
    r = st["df"].iloc[i]
    design_id = str(r.get("design_id", ""))
    refs = {k: (None if pd.isna(r.get(c)) else str(r.get(c))) for k, c in REFOLD_COL.items()}
    results = coll.add_design(folders, design_id, refs, r.to_dict(),
                              note=_s(d.get("note")),
                              dataset=_s(d.get("dataset")),
                              driven_metric=_s(d.get("driven_metric")))
    return jsonify({"ok": all(x["ok"] for x in results),
                    "design_id": design_id, "results": results})


@app.route("/api/collections/remove", methods=["POST"])
def api_collections_remove():
    d = request.get_json(force=True) or {}
    folder, did = str(d.get("folder", "")), str(d.get("design_id", ""))
    if not folder or not did:
        return jsonify({"ok": False, "error": "folder and design_id required"}), 400
    try:
        return jsonify(coll.remove_design(folder, did))
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400


@app.route("/api/collections/browse")
def api_collections_browse():
    st = _ds(request.args.get("dataset"))
    folder = request.args.get("folder", "")
    if not folder:
        return jsonify({"error": "folder required"}), 400
    try:
        recs = coll.browse_folder(folder)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    designs = [_payload_from_record(st, rec) for rec in recs]
    return jsonify({"folder": folder, "total": len(designs), "designs": designs})


@app.route("/api/collections/cif")
def api_collections_cif():
    """Serve a CIF from a folder's local cif/ dir (browse mode; no Euler round-trip)."""
    folder, name = request.args.get("folder", ""), request.args.get("name", "")
    if not folder or not name:
        abort(400)
    try:
        return send_file(coll.serve_cif_path(folder, name), mimetype="chemical/x-cif")
    except Exception as e:
        return jsonify({"error": str(e)}), 404


if __name__ == "__main__":
    load_all()
    app.run(host="127.0.0.1", port=8000, debug=False, threaded=True)
