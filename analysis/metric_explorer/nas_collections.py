"""collections.py — curated design "folders" (shortlists) materialized on the Wicky NAS.

The metric explorer lets you drop a design into one or more named folders. Each folder is a
directory on the NAS (mounted locally at /Volumes/wicky or /Volumes/bsse/wicky) containing:

    <folder>/metrics.csv   full metric record (all dataset columns) + note/added_utc/…
    <folder>/cif/          the 3 decompressed refolds (AF3/B2/ESM) per design, original names

Only Euler interaction is READ-ONLY: streaming CIF members out of the *.tar.zst archives via
`ssh euler` (same mechanism as app.py:fetch_cif). Nothing is written on Euler, and the
nightly euler->NAS sync is untouched. The folders on the NAS are the single source of truth
(the dropdown is built by listing them; membership/dedup come from each metrics.csv).
"""
from __future__ import annotations

import os
import re
import shlex
import subprocess
import threading
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# --- fixed locations ---------------------------------------------------------
# The Wicky SMB share can land at different local mount points depending on how macOS mounts it
# this session — the share-leaf name (/Volumes/wicky) or the full groups path
# (/Volumes/bsse/wicky). Auto-detect across the known candidates instead of hardcoding one.
NAS_MOUNTS = [Path("/Volumes/wicky"), Path("/Volumes/bsse/wicky")]
_LV_SUBPATH = ("users", "lvalsecchi", "Leo_Valsecchi", "LV_Coding")


def _detect_lv_root() -> Path | None:
    """LV_Coding dir under whichever candidate mount is actually mounted and has the tree;
    None when the share is mounted at no known location."""
    for mnt in NAS_MOUNTS:
        try:
            if os.path.ismount(str(mnt)) and mnt.joinpath(*_LV_SUBPATH).parent.is_dir():
                return mnt.joinpath(*_LV_SUBPATH)
        except OSError:
            continue
    return None


# best-effort at import; nas_available() keeps this current at runtime
LV_ROOT = _detect_lv_root() or NAS_MOUNTS[0].joinpath(*_LV_SUBPATH)
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SSH = ["ssh", "-o", "ConnectTimeout=20", "-o", "BatchMode=yes", "euler"]

ID_COL = "design_id"
EXTRA_COLS = ["note", "added_utc", "dataset", "driven_metric"]   # provenance (spec B + C)
REFOLD_KEYS = ("af3", "b2", "esm")

# serialize every metrics.csv read-modify-write so two quick adds can't clobber the file
_CSV_LOCK = threading.Lock()

_SAFE_NAME = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9._ \-]*$")   # leading '_' ok; leading '.' not


# --------------------------------------------------------------------------- guards
def nas_available() -> bool:
    """True when the Wicky share is mounted (at any known mount point) and the LV_Coding parent
    exists. Refreshes the module-level LV_ROOT so a share mounted at a different path this
    session is picked up without editing the constant."""
    global LV_ROOT
    lv = _detect_lv_root()
    if lv is not None:
        LV_ROOT = lv
        return True
    return False


def sanitize_folder(name: str) -> str:
    """Validate a user-supplied folder name (path-traversal guard). Returns the cleaned name
    or raises ValueError. Allows letters/digits/space/dot/underscore/hyphen; rejects empty,
    '.', '..', and anything with a path separator or leading dot."""
    n = (name or "").strip()
    if not n or n in (".", "..") or "/" in n or "\\" in n or os.sep in n:
        raise ValueError(f"invalid folder name: {name!r}")
    n = re.sub(r"\s+", " ", n)                       # collapse internal whitespace
    if n.startswith(".") or not _SAFE_NAME.match(n):
        raise ValueError(f"invalid folder name: {name!r}")
    return n


def _folder_dir(name: str) -> Path:
    return LV_ROOT / sanitize_folder(name)


# --------------------------------------------------------------------------- CIF I/O
def _parse_ref(ref: str) -> tuple[str, str, str]:
    """'<archive>.tar.zst!<member>' -> (archive, member, member_basename). For a loose path
    (no '!') archive == ref, member == ''."""
    if "!" in ref:
        arch, member = ref.split("!", 1)
        return arch, member, Path(member).name
    return ref, "", Path(ref).name


def validate_ref(ref: str) -> tuple[str, str]:
    """Confine an RB-relative ref to the repo before it reaches the ssh command. shlex.quote
    stops shell-metachar injection but NOT path escape or tar/cat option-injection, so we also
    reject: absolute paths, '~', a leading '-' (option), and any '..' escape. Returns
    (archive, member); raises ValueError otherwise. Callers additionally pass '--' to cat/tar."""
    arch, member, _ = _parse_ref(ref)
    if not arch:
        raise ValueError("empty ref")
    for p in (arch, member):
        if not p:
            continue
        if p.startswith("/") or p.startswith("~") or p.startswith("-") or os.path.isabs(p):
            raise ValueError(f"unsafe ref component: {p!r}")
        if os.path.normpath(p).startswith(".."):
            raise ValueError(f"ref escapes repo: {p!r}")
    return arch, member


def _extract_member(ref: str, dest: Path) -> bool:
    """Stream one CIF member out of its Euler tar.zst (server-side decompress) into `dest`.
    Returns True on a non-empty write, False otherwise (dest removed on failure). Confines the
    ref to RB, ends options with '--', and bounds the ssh call with a timeout."""
    try:
        arch, member = validate_ref(ref)
    except ValueError:
        return False
    if member:
        remote = (f"cd {shlex.quote(RB)} && tar --zstd -xOf {shlex.quote(arch)} "
                  f"-- {shlex.quote(member)}")
    else:
        remote = f"cd {shlex.quote(RB)} && cat -- {shlex.quote(arch)}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with open(tmp, "wb") as fh:
            p = subprocess.run(SSH + [remote], stdout=fh, stderr=subprocess.PIPE, timeout=180)
        if p.returncode == 0 and tmp.stat().st_size > 0:
            os.replace(tmp, dest)                    # atomic
            return True
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        tmp.unlink()
    except OSError:
        pass
    return False


# ------------------------------------------------------------------- metrics.csv upsert
def _read_csv(path: Path) -> pd.DataFrame | None:
    """LENIENT read for display paths (list/members/browse): None on absent OR unreadable —
    a corrupt file there just shows nothing, never destroys data."""
    if path.exists() and path.stat().st_size > 0:
        try:
            return pd.read_csv(path, low_memory=False)
        except Exception:                            # noqa: BLE001 — corrupt/empty -> treat as none
            return None
    return None


def _read_existing(path: Path) -> pd.DataFrame | None:
    """STRICT read for the upsert (write) path: None only when the file is truly absent/empty;
    a parse error on an EXISTING file PROPAGATES so we abort rather than silently overwrite the
    folder's other rows with a single new one (bounded data loss on a flaky SMB mount)."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, low_memory=False)       # parse error -> raises -> upsert aborts


def _folder_count(path: Path) -> int | None:
    """Row count read cheaply (only the design_id column, not all ~230). None if not a readable
    collection metrics.csv."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return int(len(pd.read_csv(path, usecols=[ID_COL])))
    except Exception:                                # noqa: BLE001
        return None


def _write_csv_atomic(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(".csv.tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)                            # atomic: never a half-written metrics.csv


def _upsert_row(csv_path: Path, row: dict) -> None:
    """Insert/replace `row` (keyed on design_id) in metrics.csv, union of columns, atomically.
    Uses the STRICT read so a corrupt existing file raises here instead of being overwritten."""
    new = pd.DataFrame([row])
    existing = _read_existing(csv_path)
    if existing is not None and ID_COL in existing.columns:
        existing = existing[existing[ID_COL].astype(str) != str(row[ID_COL])]
        out = pd.concat([existing, new], ignore_index=True, sort=False)
    else:
        out = new
    _write_csv_atomic(out, csv_path)


# --------------------------------------------------------------------------- public API
def list_folders() -> list[dict]:
    """[{name, count}] for every *collection* folder — i.e. a subdir of LV_ROOT that holds a
    metrics.csv (tool-managed; unrelated sibling dirs without one are skipped). [] if NAS down."""
    if not nas_available() or not LV_ROOT.is_dir():
        return []
    out = []
    for d in sorted(LV_ROOT.iterdir()):
        if not d.is_dir() or d.name.startswith("."):
            continue
        try:
            sanitize_folder(d.name)             # only surface names we can round-trip safely
        except ValueError:
            continue
        n = _folder_count(d / "metrics.csv")
        if n is None:                           # not a (readable) collection -> skip
            continue
        out.append({"name": d.name, "count": n})
    return out


def folder_members(folder: str) -> set[str]:
    """Set of design_ids already in `folder` (for tick-state / dedup)."""
    df = _read_csv(_folder_dir(folder) / "metrics.csv")
    if df is None or ID_COL not in df.columns:
        return set()
    return set(df[ID_COL].astype(str))


def add_design(folders: list[str], design_id: str, refs: dict[str, str],
               row: dict, note: str = "", dataset: str = "",
               driven_metric: str = "") -> list[dict]:
    """Add one design to each named folder: extract its (up to 3) refold CIFs and upsert its
    full metric row. Idempotent per design_id. Returns [{folder, ok, cifs, missing, msg}]."""
    if not nas_available():
        return [{"folder": f, "ok": False, "cifs": 0, "missing": list(REFOLD_KEYS),
                 "msg": "NAS not mounted"} for f in folders]
    extras = {"note": note or "", "dataset": dataset or "",
              "driven_metric": driven_metric or "",
              "added_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    results = []
    for folder in folders:
        try:
            fdir = _folder_dir(folder)                       # raises ValueError on a bad name
            cifdir = fdir / "cif"
            cifdir.mkdir(parents=True, exist_ok=True)
            # CIF extraction is slow (ssh) and touches only per-design files, so it runs
            # OUTSIDE the lock; only the metrics.csv read-modify-write below is serialized.
            got, missing = [], []
            for k in REFOLD_KEYS:
                ref = refs.get(k)
                if not ref:
                    missing.append(k)
                    continue
                _, _, base = _parse_ref(ref)
                (got if _extract_member(ref, cifdir / base) else missing).append(k)
            full = dict(row)
            full[ID_COL] = design_id
            full.update(extras)
            with _CSV_LOCK:                                  # only the CSV RMW is critical
                _upsert_row(fdir / "metrics.csv", full)
            results.append({"folder": folder, "ok": True, "cifs": len(got),
                            "missing": missing, "msg": f"saved ({len(got)}/3 refolds)"})
        except Exception as e:                               # bad name / NAS OSError / corrupt CSV
            results.append({"folder": folder, "ok": False, "cifs": 0,
                            "missing": list(REFOLD_KEYS), "msg": str(e) or type(e).__name__})
    return results


def remove_design(folder: str, design_id: str) -> dict:
    """Drop `design_id` from folder/metrics.csv and delete its cif/<design_id>_*.cif files.
    Returns {ok, msg, emptied} — `emptied` True when that was the folder's last design (the UI
    then leaves folder view, since an empty collection is delisted)."""
    if not nas_available():
        return {"ok": False, "msg": "NAS not mounted"}
    emptied = False
    with _CSV_LOCK:
        fdir = _folder_dir(folder)
        csv_path = fdir / "metrics.csv"
        df = _read_csv(csv_path)
        if df is not None and ID_COL in df.columns:
            kept = df[df[ID_COL].astype(str) != str(design_id)]
            if len(kept):
                _write_csv_atomic(kept, csv_path)
            else:
                emptied = True
                try:
                    csv_path.unlink()                # last one out -> remove the empty table
                except OSError:
                    pass
        # delete this design's cif files (names start with "<design_id>_"; safe re: prefixes
        # because the next char is always '_', so 'seq5' never matches 'seq50'). Also sweep
        # macOS AppleDouble "._<name>" sidecars the SMB mount drops next to each file.
        cifdir = fdir / "cif"
        pfx = str(design_id) + "_"
        removed = 0
        if cifdir.is_dir():
            for fn in os.listdir(cifdir):
                base = fn[2:] if fn.startswith("._") else fn        # strip AppleDouble prefix
                if base.startswith(pfx) and base.endswith(".cif"):
                    try:
                        (cifdir / fn).unlink()
                        if not fn.startswith("._"):
                            removed += 1
                    except OSError:
                        pass
        if emptied:                                  # tidy the now-empty collection dir (best effort)
            _sweep_dotfiles(cifdir)
            for d in (cifdir, fdir):
                try:
                    d.rmdir()                        # only succeeds if truly empty
                except OSError:
                    pass
    return {"ok": True, "emptied": emptied, "msg": f"removed ({removed} cif files)"}


def _sweep_dotfiles(d: Path) -> None:
    """Remove stray macOS AppleDouble/.DS_Store files so an emptied dir can rmdir cleanly."""
    if not d.is_dir():
        return
    for fn in os.listdir(d):
        if fn.startswith("._") or fn == ".DS_Store":
            try:
                (d / fn).unlink()
            except OSError:
                pass


def browse_folder(folder: str) -> list[dict]:
    """Raw metric rows (list of dicts) for `folder`, in file order. app.py shapes these into
    display payloads and serves each AF3 CIF from the local NAS via /api/collections/cif."""
    df = _read_csv(_folder_dir(folder) / "metrics.csv")
    if df is None:
        return []
    return df.to_dict(orient="records")


def serve_cif_path(folder: str, name: str) -> Path:
    """Resolve <folder>/cif/<name> for the browse viewer, refusing any path escape."""
    base = (_folder_dir(folder) / "cif").resolve()
    p = (base / Path(name).name).resolve()           # basename only
    if p.parent != base or not p.is_file():
        raise ValueError(f"no such cif: {name!r}")
    return p
