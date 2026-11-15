"""export.py — write the passing set out as a CSV, and optionally pull its structures.

Two exports:

  CSV   design_id + every tool metric + the reconstructed CIF refs, plus a header comment
        carrying the exact filter spec that produced the set. Streamed, so a multi-million
        row export does not build a DataFrame in memory.

  CIFs  the design backbone + the AF3 and ESM refolds per design, streamed out of the Euler
        tar.zst archives into exports/<name>/cif/. Background thread, poll for progress,
        idempotent so re-running a name resumes.

The ssh extraction primitive is imported from metric_explorer's nas_collections rather than
reimplemented: `_extract_member` carries a security fix (an earlier version allowed an
arbitrary Euler file read through `?ref=`) and duplicating it would fork that hardening.
Only the primitive is borrowed — the b2-shaped key set and job registry of
metric_explorer's local_download do not fit a tool that excludes Boltz-2.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
EXPORTS_ROOT = APP_DIR / "exports"
CIF_CONFIRM_CAP = 500          # above this many designs the caller must pass confirm=True
CSV_CONFIRM_CAP = 500_000
_MAX_WORKERS = 6
KEYS = ("design", "af3", "esm")           # no b2, by design

sys.path.insert(0, str(Path.home() / "Documents" / "metric_explorer"))
try:
    from nas_collections import _extract_member, _parse_ref, sanitize_folder
    EXTRACT_OK, EXTRACT_ERR = True, ""
except Exception as e:                                        # noqa: BLE001
    EXTRACT_OK, EXTRACT_ERR = False, f"{type(e).__name__}: {e}"

    def sanitize_folder(name):                                # minimal fallback for CSV names
        s = "".join(c for c in str(name) if c.isalnum() or c in "._- ").strip()
        if not s or s in (".", ".."):
            raise ValueError("bad folder name")
        return s

_JOBS: dict = {}
_LOCK = threading.Lock()


# ------------------------------------------------------------------------------ CSV

def export_csv(store, mask, spec, name: str) -> dict:
    """Stream the passing rows to exports/<name>.csv. Returns {path, n, bytes}."""
    import numpy as np

    EXPORTS_ROOT.mkdir(parents=True, exist_ok=True)
    dest = EXPORTS_ROOT / f"{sanitize_folder(name)}.csv"
    idx = np.flatnonzero(mask)
    numeric = store.numeric_cols()
    cats = [c for c in ("stack", "campaign", "gen", "invfold") if c in store]

    tmp = dest.with_suffix(".csv.part")
    with open(tmp, "w", newline="") as f:
        f.write(f"# binder_filter export  n={idx.size}  store_n={store.n}\n")
        f.write(f"# filter={json.dumps(spec, separators=(',', ':'))}\n")
        w = csv.writer(f, lineterminator="\n")     # csv defaults to CRLF; mixed endings once
        w.writerow(["design_id", "uuid8", "length"] + cats + numeric +
                   ["design_cif", "af3_cif", "esm_cif"])
        ids = store.col("design_id")
        uus = store.col("uuid8")
        length = store.col("length")
        ncols = {c: store.col(c) for c in numeric}
        ccols = {c: store.col(c) for c in cats}
        for i in idx.tolist():
            p = store.paths(i)
            w.writerow([ids[i].decode(), uus[i].decode(), int(length[i])] +
                       [store.label(c, int(ccols[c][i])) for c in cats] +
                       [f"{ncols[c][i]:.6g}" for c in numeric] +
                       [p["designed_cif_path"], p["af3_msaoff_cif_path"], p["esm_cif_path"]])
    os.replace(tmp, dest)
    return {"path": str(dest), "n": int(idx.size), "bytes": dest.stat().st_size}


# ----------------------------------------------------------------------------- CIFs

def job_status(job_id: str):
    with _LOCK:
        j = _JOBS.get(job_id)
        return dict(j) if j is not None else None


def _record(job_id, *, cifs_inc=0, done_inc=0, missing=None):
    with _LOCK:
        j = _JOBS[job_id]
        j["cifs"] += cifs_inc
        j["done"] += done_inc
        if missing is not None and len(j["missing"]) < 200:
            j["missing"].append(missing)


def _extract_one(folder: Path, ref: str) -> bool:
    _, _, base = _parse_ref(ref)
    dest = folder / "cif" / base
    if dest.exists() and dest.stat().st_size > 0:
        return True                                   # idempotent: re-running a name resumes
    return bool(_extract_member(ref, dest))


def start_cifs(designs: list[dict], name: str) -> str:
    """designs = [{design_id, refs:{design,af3,esm}}]. Returns a job id immediately."""
    if not EXTRACT_OK:
        raise RuntimeError(f"ssh extraction unavailable: {EXTRACT_ERR}")
    folder = EXPORTS_ROOT / sanitize_folder(name)
    (folder / "cif").mkdir(parents=True, exist_ok=True)
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        _JOBS[job_id] = {"status": "running", "total": len(designs), "done": 0, "cifs": 0,
                         "missing": [], "folder": str(folder), "error": ""}
    threading.Thread(target=_run, args=(job_id, folder, designs), daemon=True).start()
    return job_id


def _run(job_id, folder: Path, designs):
    try:
        tasks, remaining = [], {}
        for di, d in enumerate(designs):
            refs = d.get("refs") or {}
            present = [k for k in KEYS if refs.get(k)]
            remaining[di] = len(present)
            for k in present:
                tasks.append((di, d["design_id"], k, refs[k]))
            for k in KEYS:
                if not refs.get(k):
                    _record(job_id, missing=f"{d['design_id']}:{k}")
            if not present:
                _record(job_id, done_inc=1)

        with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
            fut = {ex.submit(_extract_one, folder, ref): (di, did, k)
                   for (di, did, k, ref) in tasks}
            for f in as_completed(fut):
                di, did, k = fut[f]
                try:
                    ok = bool(f.result())
                except Exception:                     # noqa: BLE001 — a failure is a miss
                    ok = False
                _record(job_id, cifs_inc=1 if ok else 0,
                        missing=None if ok else f"{did}:{k}")
                remaining[di] -= 1
                if remaining[di] == 0:
                    _record(job_id, done_inc=1)
        with _LOCK:
            _JOBS[job_id]["status"] = "done"
    except Exception as e:                            # noqa: BLE001
        with _LOCK:
            _JOBS[job_id]["status"] = "error"
            _JOBS[job_id]["error"] = str(e) or type(e).__name__
