"""
Read campaign state from the eulermount filesystem.

EULERMOUNT_ROOT is the local mountpoint that mirrors /cluster/work/wicky/lvalsecchi/.
We only do *reads* here — anything that needs to run a command (sbatch, squeue)
goes through euler.ssh().
"""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable
import json
import csv
import re
import os
import time


# Adjust to wherever you mount /cluster/work/wicky/lvalsecchi/ locally.
# Typical: ~/euler/rna_binder/   or   /Volumes/euler/lvalsecchi/rna_binder/
EULERMOUNT_ROOT_DEFAULT = "/Volumes/eulermount/lvalsecchi/rna_binder"


@dataclass
class RunInfo:
    stack: str           # e.g. "rfd3__ligandmpnn"
    name: str            # e.g. "midrun_1ca62a01"
    path: Path
    n_designs: int
    n_af3_refolds: int
    n_b2_refolds: int
    n_esm_refolds: int
    has_meta: bool
    label: str           # launch label shared across a batch, e.g. "midrun"
    submitted: str       # display timestamp "YYYY-MM-DD HH:MM"
    sort_epoch: float    # seconds since epoch, for newest-first sorting


# Run dirs are named "<label>_<uuid8>"; strip the 8-hex-char suffix to recover
# the launch label that a multi-stack batch shares.
_UUID_SUFFIX_RE = re.compile(r"_[0-9a-f]{8}$")


def _read_meta(run_dir: Path) -> dict:
    p = run_dir / "run_meta.json"
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except (ValueError, OSError):
        return {}


def _run_label(run_dir: Path, meta: dict) -> str:
    return meta.get("label") or _UUID_SUFFIX_RE.sub("", run_dir.name)


def _run_time(run_dir: Path, meta: dict) -> tuple[float, str]:
    """(sort_epoch, display) from run_meta 'submitted', falling back to mtime."""
    iso = meta.get("submitted")
    if iso:
        try:
            dt = datetime.fromisoformat(iso)
            return dt.timestamp(), dt.strftime("%Y-%m-%d %H:%M")
        except ValueError:
            pass
    ts = run_dir.stat().st_mtime
    return ts, datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M")


def project_root(root: str | Path = EULERMOUNT_ROOT_DEFAULT) -> Path:
    return Path(root).expanduser()


def list_stacks(root: str | Path = EULERMOUNT_ROOT_DEFAULT,
                target: str = "csy4_crrna_free") -> list[str]:
    p = project_root(root) / "campaigns" / target
    if not p.exists():
        return []
    return sorted([d.name for d in p.iterdir() if d.is_dir()])


def list_runs(root: str | Path = EULERMOUNT_ROOT_DEFAULT,
              target: str = "csy4_crrna_free",
              refold_counts: dict[tuple[str, str], dict[str, int]] | None = None,
              ) -> list[RunInfo]:
    """Walk all stacks under campaigns/<target>/ and summarise each run dir.

    refold_counts: optional {(stack, run): {"af3": n, "b2": n}} from a single
    server-side count (euler.refold_counts). When given, we skip the slow
    recursive rglob over the mount. When None, fall back to local counting.
    """
    base = project_root(root) / "campaigns" / target
    if not base.exists():
        return []

    runs: list[RunInfo] = []
    for stack_dir in sorted(base.iterdir()):
        if not stack_dir.is_dir():
            continue
        for run_dir in sorted(stack_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            # eulermount (sshfs) can drop mid-walk and raise OSError; skip the
            # offending run rather than crashing the whole Recent runs panel.
            try:
                n_des = _count_designs(run_dir)
                if refold_counts is not None:
                    c = refold_counts.get((stack_dir.name, run_dir.name), {})
                    n_af3, n_b2, n_esm = c.get("af3", 0), c.get("b2", 0), c.get("esm", 0)
                else:
                    n_af3 = _count_cifs(run_dir / "refold_af3")
                    n_b2 = _count_cifs(run_dir / "refold_b2")
                    n_esm = _count_cifs(run_dir / "refold_esm")
                meta = _read_meta(run_dir)
                sort_epoch, submitted = _run_time(run_dir, meta)
            except OSError:
                continue
            runs.append(RunInfo(
                stack=stack_dir.name,
                name=run_dir.name,
                path=run_dir,
                n_designs=n_des,
                n_af3_refolds=n_af3,
                n_b2_refolds=n_b2,
                n_esm_refolds=n_esm,
                has_meta=bool(meta),
                label=_run_label(run_dir, meta),
                submitted=submitted,
                sort_epoch=sort_epoch,
            ))
    return runs


def _count_designs(run_dir: Path) -> int:
    """Count a run's design files, archive-aware.

    Phase-2 packs designs to cut inodes, so a plain iterdir reports 0 for packed
    runs. Three layouts:
      * whole-dir archive  <run>/designs.tar.zst(.list)        — single-run / HD
      * per-uuid slices    <run>/designs/<uuid>.tar.zst(.list) — --num-runs>1
      * loose              <run>/designs/*.cif                 — unpacked / legacy
    Count the .list manifest(s) when present; else count loose .pdb/.cif/.gz.
    """
    def _count_list(p: Path) -> int:
        try:
            return sum(1 for ln in p.read_text().splitlines()
                       if ln.rsplit(".", 1)[-1] in ("pdb", "cif", "gz"))
        except OSError:
            return 0

    manifest = run_dir / "designs.tar.zst.list"
    if manifest.exists():
        return _count_list(manifest)
    designs_dir = run_dir / "designs"
    if designs_dir.exists():
        slice_lists = sorted(designs_dir.glob("*.tar.zst.list"))
        if slice_lists:
            return sum(_count_list(p) for p in slice_lists)
        return sum(1 for p in designs_dir.iterdir()
                   if p.suffix in (".pdb", ".cif", ".gz"))
    return 0


def _count_cifs(root: Path) -> int:
    if not root.exists():
        return 0
    n = 0
    try:
        for _ in root.rglob("*.cif"):
            n += 1
    except OSError:
        pass  # eulermount dropped mid-walk; return what we counted
    return n


# ── incremental disk cache ───────────────────────────────────────────────────
# Finished runs never change, so we cache each run's summary locally and only
# recompute runs younger than RECENT_HOURS (and brand-new ones). The slow bits
# (designs iterdir + server-side refold find) then run on a handful of runs
# instead of all of them. Refresh wipes the cache for a full recount.
_CACHE_DIR = Path(os.path.expanduser("~/.cache/rna_binder_gui"))
RECENT_HOURS = 48


def _cache_path(target: str) -> Path:
    return _CACHE_DIR / f"runs_{target}.json"


def clear_runs_cache(target: str) -> None:
    p = _cache_path(target)
    try:
        p.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        pass


def _info_to_dict(r: RunInfo) -> dict:
    return {
        "stack": r.stack, "name": r.name, "n_designs": r.n_designs,
        "n_af3_refolds": r.n_af3_refolds, "n_b2_refolds": r.n_b2_refolds,
        "n_esm_refolds": r.n_esm_refolds,
        "has_meta": r.has_meta, "label": r.label,
        "submitted": r.submitted, "sort_epoch": r.sort_epoch,
    }


def _dict_to_info(d: dict, base: Path) -> RunInfo:
    return RunInfo(
        stack=d["stack"], name=d["name"], path=base / d["stack"] / d["name"],
        n_designs=d["n_designs"], n_af3_refolds=d["n_af3_refolds"],
        n_b2_refolds=d["n_b2_refolds"], n_esm_refolds=d.get("n_esm_refolds", 0),
        has_meta=d["has_meta"],
        label=d["label"], submitted=d["submitted"], sort_epoch=d["sort_epoch"],
    )


def _run_uuid(name: str) -> str:
    m = re.search(r"([0-9a-f]{8})$", name)
    return m.group(1) if m else ""


def list_runs_incremental(
    root: str | Path = EULERMOUNT_ROOT_DEFAULT,
    target: str = "csy4_crrna_free",
    count_fn: Callable[[list[tuple[str, str]]], dict] | None = None,
    active_uuids: set[str] | frozenset[str] | None = None,
    refolds_active: bool = False,
    idle_fallback_hours: float = 6.0,
    refold_fallback_hours: float = 24.0,
) -> list[RunInfo]:
    """list_runs with a local disk cache; recompute only runs that may change.

    A run is recomputed when it is new, OR its uuid8 matches an active SLURM job
    (design/multif jobs carry the uuid8 in their name), OR it falls inside an
    age fallback. The fallback is short when no refold jobs are queued and wider
    while they are (refold jobs aren't uniquely named per run, so we can't map
    them precisely). Everything else is served from the local cache.

    count_fn([(stack, run), ...]) -> {(stack, run): {"af3": n, "b2": n}} does
    the refold count (e.g. euler.refold_counts scoped to those runs). If None,
    refold counts fall back to local rglob for the recomputed runs.
    """
    active_uuids = active_uuids or frozenset()
    base = project_root(root) / "campaigns" / target
    if not base.exists():
        return []

    cache: dict[str, dict] = {}
    cpath = _cache_path(target)
    if cpath.exists():
        try:
            cache = json.loads(cpath.read_text())
        except (ValueError, OSError):
            cache = {}

    # cheap enumeration: stack dirs + their immediate children (run dirs)
    current: list[tuple[str, str, Path]] = []
    try:
        for stack_dir in sorted(base.iterdir()):
            if not stack_dir.is_dir():
                continue
            for run_dir in sorted(stack_dir.iterdir()):
                if run_dir.is_dir():
                    current.append((stack_dir.name, run_dir.name, run_dir))
    except OSError:
        pass

    horizon = refold_fallback_hours if refolds_active else idle_fallback_hours
    cutoff = time.time() - horizon * 3600
    result: dict[str, dict] = {}
    to_recompute: list[tuple[str, str, Path]] = []
    for stack, name, run_dir in current:
        key = f"{stack}/{name}"
        c = cache.get(key)
        recompute = (
            c is None                                   # never seen
            or _run_uuid(name) in active_uuids          # live SLURM job
            or c.get("sort_epoch", 0) >= cutoff         # within age fallback
        )
        if recompute:
            to_recompute.append((stack, name, run_dir))
        else:
            result[key] = c            # finished & idle -> trust cache

    if to_recompute:
        counts: dict = {}
        if count_fn is not None:
            try:
                counts = count_fn([(s, n) for s, n, _ in to_recompute]) or {}
            except Exception:
                counts = {}
        for stack, name, run_dir in to_recompute:
            key = f"{stack}/{name}"
            try:
                n_des = _count_designs(run_dir)
                meta = _read_meta(run_dir)
                sort_epoch, submitted = _run_time(run_dir, meta)
                if count_fn is not None:
                    c = counts.get((stack, name), {})
                    n_af3, n_b2, n_esm = c.get("af3", 0), c.get("b2", 0), c.get("esm", 0)
                else:
                    n_af3 = _count_cifs(run_dir / "refold_af3")
                    n_b2 = _count_cifs(run_dir / "refold_b2")
                    n_esm = _count_cifs(run_dir / "refold_esm")
            except OSError:
                if key in cache:       # mount hiccup -> keep prior value
                    result[key] = cache[key]
                continue
            result[key] = {
                "stack": stack, "name": name, "n_designs": n_des,
                "n_af3_refolds": n_af3, "n_b2_refolds": n_b2, "n_esm_refolds": n_esm,
                "has_meta": bool(meta), "label": _run_label(run_dir, meta),
                "submitted": submitted, "sort_epoch": sort_epoch,
            }

    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cpath.write_text(json.dumps(result))
    except OSError:
        pass

    return [_dict_to_info(d, base) for d in result.values()]


def read_processing_csv(run_dir: Path) -> list[dict]:
    """If the run has a processing_results.csv (HD), return its rows."""
    csv_path = run_dir / "designs" / "processing_results.csv"
    if not csv_path.exists():
        return []
    with csv_path.open() as f:
        return list(csv.DictReader(f))
