"""local_download.py — bulk-copy the filtered design set to a named LOCAL folder.

The metric explorer's "Download selected" action copies every design passing the current
hard-filter locks into `downloads/<name>/`:

    downloads/<name>/cif/          the 3 refolds (AF3/B2/ESM) per design, original member names
    downloads/<name>/metrics.csv   the full metric rows of the selected set

Structures are streamed READ-ONLY out of the Euler *.tar.zst archives via
`nas_collections._extract_member` (same ssh mechanism as the collections feature and app.py).
The copy runs in a background thread with a small pool; callers poll `job_status()` for progress.
Idempotent: a CIF already present on disk is skipped, so a folder name can be re-used to resume.
"""
from __future__ import annotations

import threading
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import nas_collections as _nc     # reuse the ssh extraction + name-validation helpers

APP_DIR = Path(__file__).resolve().parent
DOWNLOADS_ROOT = APP_DIR / "downloads"
CAP = 500                         # above this many designs the caller must confirm
REFOLD_KEYS = ("af3", "b2", "esm")
_MAX_WORKERS = 6

_JOBS: dict = {}                  # job_id -> progress dict
_LOCK = threading.Lock()


def sanitize_name(name: str) -> str:
    """Validate a user-supplied folder name (reuses the collections rules: letters/digits/
    space/dot/underscore/hyphen; no path separators; no leading dot; rejects ''/'.'/'..')."""
    return _nc.sanitize_folder(name)


def _job_dir(name: str) -> Path:
    return DOWNLOADS_ROOT / sanitize_name(name)


def job_status(job_id: str) -> dict | None:
    """Snapshot (copy) of a job's progress, or None for an unknown id."""
    with _LOCK:
        j = _JOBS.get(job_id)
        return dict(j) if j is not None else None


def _record(job_id: str, *, cifs_inc: int = 0, done_inc: int = 0, missing: str | None = None) -> None:
    with _LOCK:
        j = _JOBS[job_id]
        j["cifs"] += cifs_inc
        j["done"] += done_inc
        if missing is not None:
            j["missing"].append(missing)


def _finish(job_id: str, status: str, error: str = "") -> None:
    with _LOCK:
        _JOBS[job_id]["status"] = status
        _JOBS[job_id]["error"] = error


def _extract_one(folder: Path, ref: str) -> bool:
    """Extract one refold CIF into <folder>/cif/. Idempotent: an existing non-empty file is
    kept (returns True without an ssh round-trip)."""
    _, _, base = _nc._parse_ref(ref)
    dest = folder / "cif" / base
    if dest.exists() and dest.stat().st_size > 0:
        return True
    return bool(_nc._extract_member(ref, dest))


def start(name: str, designs: list[dict]) -> str:
    """Kick off a background copy of `designs` into downloads/<name>/. Returns a job id
    immediately (non-blocking). Each design is {design_id, refs:{af3,b2,esm}, row:dict}.
    Raises ValueError on a bad folder name (before any work)."""
    folder = _job_dir(name)                              # raises ValueError on bad name
    (folder / "cif").mkdir(parents=True, exist_ok=True)
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        _JOBS[job_id] = {"status": "running", "total": len(designs), "done": 0,
                         "cifs": 0, "missing": [], "folder": str(folder),
                         "name": sanitize_name(name), "error": ""}
    threading.Thread(target=_run, args=(job_id, folder, designs), daemon=True).start()
    return job_id


def _run(job_id: str, folder: Path, designs: list[dict]) -> None:
    try:
        tasks = []                                       # (design_idx, design_id, key, ref)
        remaining = {}                                   # design_idx -> refs still pending
        for di, d in enumerate(designs):
            refs = d.get("refs") or {}
            present_keys = [k for k in REFOLD_KEYS if refs.get(k)]
            remaining[di] = len(present_keys)
            for k in REFOLD_KEYS:
                if refs.get(k):
                    tasks.append((di, d["design_id"], k, refs[k]))
                else:
                    _record(job_id, missing=f"{d['design_id']}:{k}")
            if remaining[di] == 0:                       # a design with no refs at all
                _record(job_id, done_inc=1)

        with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as ex:
            futmap = {ex.submit(_extract_one, folder, ref): (di, did, k)
                      for (di, did, k, ref) in tasks}
            for fut in as_completed(futmap):             # consumed in THIS thread only
                di, did, k = futmap[fut]
                ok = False
                try:
                    ok = fut.result()
                except Exception:                        # noqa: BLE001 — extraction error = missing
                    ok = False
                if ok:
                    _record(job_id, cifs_inc=1)
                else:
                    _record(job_id, missing=f"{did}:{k}")
                remaining[di] -= 1
                if remaining[di] == 0:
                    _record(job_id, done_inc=1)

        pd.DataFrame([d["row"] for d in designs]).to_csv(folder / "metrics.csv", index=False)
        _finish(job_id, "done")
    except Exception as e:                               # noqa: BLE001
        _finish(job_id, "error", str(e) or type(e).__name__)
