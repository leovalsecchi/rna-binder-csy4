#!/usr/bin/env python3
"""build_filter_store.py — compact columnar store for the interactive binder-filter tool.

The master CSV is 9.7 GB / 172 columns, which is neither shippable to a laptop nor fast to
filter. The filter tool needs ~20 columns and needs them in RAM, so this streams the master
once and writes one .npy per column into a store directory (~740 MB for 5.34M designs).

⚠️ BOLTZ-2 IS EXCLUDED. Every bare consensus column in the master (consensus_iptm, iptm_min3,
dockq_self, n_h_bonds_interface, pi_stack_score, ...) folds B2 in. The tool must read ONLY the
_af3 / _esm columns and recompute any 2-backend consensus itself. `_assert_no_consensus`
enforces that no bare consensus column can reach the store.

⚠️ STREAM IT. Do NOT pd.read_csv this file — a 9.6 GB pandas read hung a login node for 20
minutes on this project. Run under sbatch, not on a login node.

PATH COMPRESSION. The three CIF path columns are ~150 chars each (450 B/row = 2.4 GB). They
are fully determined by (design_id, uuid8, campaign layout), so each row stores a single
int16 template index instead. A template is the literal path with the design_id replaced by
<ID> and the uuid8 by <U>; reconstruction is literal substitution, and the build ASSERTS a
byte-exact round-trip on every row for all three paths. Any layout this script has not seen
becomes its own template rather than silently corrupting a path.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from rna_binder.analysis.seq_complexity import seq_complexity

csv.field_size_limit(10 ** 7)

ID_WIDTH = 64
UUID_WIDTH = 16
SEQ_CACHE_MAX = 2_000_000          # memoise seq metrics; HD emits many identical sequences

# tool metric -> master column. AF3/ESM only, never a bare consensus column.
FLOAT_COLS = {
    "rmsd_af3":  "dockq_LRMSD_af3",      # RNA-aligned ligand RMSD (pose deviation vs design)
    "rmsd_esm":  "dockq_LRMSD_esm",
    "dockq_af3": "dockq_self_af3",       # RNA-aligned DockQ
    "dockq_esm": "dockq_self_esm",
    "iptm_af3":  "af3_msaoff_iptm",
    "iptm_esm":  "esm_iptm",
    "hb_af3":    "n_h_bonds_interface_af3",
    "hb_esm":    "n_h_bonds_interface_esm",
    "pis_af3":   "pi_stack_score_af3",
    "pis_esm":   "pi_stack_score_esm",
}
# derived in-build from the `sequence` column
SEQ_OUT = {"seq_max_aa_frac": "seq_max_aa_frac", "seq_run": "seq_longest_run",
           "seq_kmer": "seq_max_kmer"}
# no campaign-wide data yet; written as all-NaN so the UI slot exists and a backfill is a
# pure data swap rather than a schema change
PLACEHOLDER_COLS = ("rog_af3", "rog_esm")
CODE_COLS = ("stack", "campaign", "gen", "invfold")
PATH_COLS = ("af3_msaoff_cif_path", "esm_cif_path", "designed_cif_path")

# a bare (B2-contaminated) column reaching the store would silently poison every filter
CONSENSUS_BARE = {
    "consensus_iptm", "consensus_iptm_legacy", "mean_iptm", "mean_iptm_legacy",
    "iptm_min3", "iptm_best2", "mean_iptm3", "dockq_self", "dockq_LRMSD", "dockq_fnat",
    "dockq_iRMSD", "dockq_std_self", "dockq_std_LRMSD", "n_h_bonds_interface",
    "pi_stack_score", "n_h_bonds_interface_opth", "buns_vsasa", "buns_vsasa_opth",
    "delta_unsat_hbonds", "dsasa_total", "red_total_per_ifaceres",
}


def _assert_no_consensus():
    """Every mapped column must be backend-specific, and its backend must match the tool key.
    The backend token is not always a suffix (af3_msaoff_iptm carries it as a prefix), so this
    checks membership, not position. b2 must never appear."""
    bad = sorted(set(FLOAT_COLS.values()) & CONSENSUS_BARE)
    if bad:
        raise SystemExit(f"ABORT: B2-contaminated consensus column(s) in the schema: {bad}")
    for k, v in FLOAT_COLS.items():
        want = "af3" if k.endswith("_af3") else "esm" if k.endswith("_esm") else None
        if want is None:
            raise SystemExit(f"ABORT: tool key {k} has no _af3/_esm backend suffix")
        toks = v.lower().replace("-", "_").split("_")
        if "b2" in toks or "boltz2" in toks:
            raise SystemExit(f"ABORT: {k} -> {v} is a Boltz-2 column; B2 is excluded")
        if want not in toks:
            raise SystemExit(f"ABORT: {k} -> {v} does not carry the '{want}' backend token")


def _f32(v):
    if v is None or v == "":
        return np.float32(np.nan)
    try:
        return np.float32(v)
    except ValueError:
        return np.float32(np.nan)


class Vocab:
    """String -> small int code, in first-seen order."""

    def __init__(self, name):
        self.name, self.items, self.index = name, [], {}

    def code(self, s):
        c = self.index.get(s)
        if c is None:
            c = len(self.items)
            if c > 32000:
                raise SystemExit(f"ABORT: vocab {self.name} exceeded 32000 entries")
            self.items.append(s)
            self.index[s] = c
        return c


def templatize(path: str, did: str, uuid8: str) -> str:
    """Replace design_id then uuid8 with placeholders. design_id CONTAINS the uuid8, so the
    longer token must go first; whatever uuid8 occurrences remain are the standalone ones."""
    return path.replace(did, "<ID>").replace(uuid8, "<U>")


def rebuild(tpl: str, did: str, uuid8: str) -> str:
    return tpl.replace("<ID>", did).replace("<U>", uuid8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="master fullmetrics CSV")
    ap.add_argument("--out", required=True, help="store directory to create")
    ap.add_argument("--n-hint", type=int, default=6_000_000)
    ap.add_argument("--limit", type=int, default=0, help="stop after N rows (smoke test)")
    a = ap.parse_args()

    _assert_no_consensus()
    os.makedirs(a.out, exist_ok=True)

    n_hint = a.n_hint
    ids = np.zeros(n_hint, dtype=f"S{ID_WIDTH}")
    uus = np.zeros(n_hint, dtype=f"S{UUID_WIDTH}")
    tpl_idx = np.zeros(n_hint, dtype=np.int16)
    length = np.zeros(n_hint, dtype=np.int16)
    codes = {c: np.zeros(n_hint, dtype=np.int16) for c in CODE_COLS}
    floats = {k: np.zeros(n_hint, dtype=np.float32) for k in FLOAT_COLS}
    seqf = {k: np.zeros(n_hint, dtype=np.float32) for k in SEQ_OUT}

    vocabs = {c: Vocab(c) for c in CODE_COLS}
    tpl_vocab, tpl_list = {}, []
    seq_cache = {}
    n = 0
    max_id_len = max_uuid_len = 0
    t0 = time.time()

    with open(a.csv, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in list(FLOAT_COLS.values()) + list(PATH_COLS) +
                   ["design_id", "uuid8", "sequence", "length"] + list(CODE_COLS)
                   if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"ABORT: master CSV lacks columns: {missing}")

        for row in reader:
            if n >= len(ids):                       # grow 1.5x rather than guess perfectly
                new = int(len(ids) * 1.5)
                print(f"  growing {len(ids):,} -> {new:,}", flush=True)
                ids = np.resize(ids, new); uus = np.resize(uus, new)
                tpl_idx = np.resize(tpl_idx, new); length = np.resize(length, new)
                codes = {k: np.resize(v, new) for k, v in codes.items()}
                floats = {k: np.resize(v, new) for k, v in floats.items()}
                seqf = {k: np.resize(v, new) for k, v in seqf.items()}

            did, u = row["design_id"], row["uuid8"]
            max_id_len = max(max_id_len, len(did))
            max_uuid_len = max(max_uuid_len, len(u))
            if len(did) > ID_WIDTH or len(u) > UUID_WIDTH:
                raise SystemExit(f"ABORT row {n}: id/uuid exceeds fixed width: {did} {u}")

            key = tuple(templatize(row[p], did, u) for p in PATH_COLS)
            ti = tpl_vocab.get(key)
            if ti is None:
                ti = len(tpl_list)
                tpl_vocab[key] = ti
                tpl_list.append(key)
                print(f"  template {ti}: {key[0]}", flush=True)
            # byte-exact round-trip, EVERY row — a silent path corruption would only surface
            # much later as a failed ChimeraX open
            for tpl, p in zip(key, PATH_COLS):
                if rebuild(tpl, did, u) != row[p]:
                    raise SystemExit(f"ABORT row {n}: path round-trip failed for {p}\n"
                                     f"  orig={row[p]}\n  back={rebuild(tpl, did, u)}")

            ids[n], uus[n], tpl_idx[n] = did.encode(), u.encode(), ti
            try:
                length[n] = int(float(row["length"] or 0))
            except ValueError:
                length[n] = 0
            for c in CODE_COLS:
                codes[c][n] = vocabs[c].code(row.get(c) or "")
            for k, src in FLOAT_COLS.items():
                floats[k][n] = _f32(row.get(src))

            s = row.get("sequence") or ""
            m = seq_cache.get(s)
            if m is None:
                m = seq_complexity(s)
                if len(seq_cache) < SEQ_CACHE_MAX:
                    seq_cache[s] = m
            for k, src in SEQ_OUT.items():
                v = m.get(src)
                seqf[k][n] = np.float32(np.nan) if v is None else np.float32(v)

            n += 1
            if n % 500_000 == 0:
                print(f"  {n:,} rows  {time.time()-t0:.0f}s  "
                      f"{len(tpl_list)} templates  seq_cache={len(seq_cache):,}", flush=True)
            if a.limit and n >= a.limit:
                break

    print(f"read {n:,} rows in {time.time()-t0:.0f}s", flush=True)

    def save(name, arr):
        np.save(os.path.join(a.out, f"{name}.npy"), arr[:n])

    save("design_id", ids); save("uuid8", uus)
    save("tpl_idx", tpl_idx); save("length", length)
    for c in CODE_COLS:
        save(c, codes[c])
    for k in FLOAT_COLS:
        save(k, floats[k])
    for k in SEQ_OUT:
        save(k, seqf[k])
    for k in PLACEHOLDER_COLS:                       # all-NaN until the backfill lands
        save(k, np.full(n, np.nan, dtype=np.float32))

    meta = {
        "n": n,
        "source": os.path.abspath(a.csv),
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "path_cols": list(PATH_COLS),
        "templates": [list(t) for t in tpl_list],
        "vocabs": {c: vocabs[c].items for c in CODE_COLS},
        "float_cols": {k: v for k, v in FLOAT_COLS.items()},
        "seq_cols": list(SEQ_OUT),
        "placeholder_cols": list(PLACEHOLDER_COLS),
        "max_id_len": max_id_len,
        "max_uuid_len": max_uuid_len,
    }
    with open(os.path.join(a.out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"\nwrote {a.out}")
    print(f"  n={n:,}  templates={len(tpl_list)}  max_id_len={max_id_len}")
    for c in CODE_COLS:
        print(f"  {c:10s} {len(vocabs[c].items)} values: {vocabs[c].items[:8]}")
    for k in list(FLOAT_COLS) + list(SEQ_OUT):
        arr = floats.get(k, seqf.get(k))[:n]
        ok = np.isfinite(arr)
        print(f"  {k:16s} {100*ok.mean():6.2f}% non-null   "
              f"min={np.nanmin(arr) if ok.any() else float('nan'):.4g} "
              f"max={np.nanmax(arr) if ok.any() else float('nan'):.4g}")
    total = sum(os.path.getsize(os.path.join(a.out, f)) for f in os.listdir(a.out))
    print(f"  store size: {total/1e6:.0f} MB")


if __name__ == "__main__":
    main()
