"""presets.py — named filter presets for the metric explorer.

A preset is a saved filter state you can reapply from a list: the typed hard-filter values
(verbatim, e.g. "<=2" / ">=3" / "=0"), plus which metric was the ● driven/sort metric and its
threshold and Δ. Stored as a JSON list in `presets.json` next to the app, so presets survive a
browser cache clear and are visible from any browser.

    [{"name": "self-consistent", "filters": {"apo_rmsd": "<=2"},
      "driven": "apo_rmsd", "thresh": 2.0, "delta": 1.0,
      "saved_utc": "2026-07-30T09:12:00+00:00"}]

Reads on the DISPLAY path are lenient (a corrupt file just shows no presets); reads on the WRITE
path are strict, so a corrupt presets.json aborts the save instead of being replaced by a single
record (bounded data loss, same principle as nas_collections._read_existing).
"""
from __future__ import annotations

import json
import os
import re
import threading
from datetime import datetime, timezone
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
PRESETS_FILE = APP_DIR / "presets.json"
MAX_NAME = 120

_LOCK = threading.Lock()


def sanitize_preset_name(name: str) -> str:
    """Trim + collapse inner whitespace; require non-empty and <= MAX_NAME chars. A preset name is
    a JSON key (never a filename), so the character set is permissive — no traversal risk."""
    n = re.sub(r"\s+", " ", (name or "").strip())
    if not n:
        raise ValueError("preset name required")
    if len(n) > MAX_NAME:
        raise ValueError(f"preset name too long (max {MAX_NAME})")
    return n


def _load_strict() -> list[dict]:
    """Parse presets.json; raises on a corrupt EXISTING file so writes abort instead of wiping."""
    if not PRESETS_FILE.exists() or PRESETS_FILE.stat().st_size == 0:
        return []
    data = json.loads(PRESETS_FILE.read_text())        # JSONDecodeError propagates
    if not isinstance(data, list):
        raise ValueError("presets.json is not a list")
    return [x for x in data if isinstance(x, dict)]


def _write(items: list[dict]) -> None:
    tmp = PRESETS_FILE.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(items, indent=1))
    os.replace(tmp, PRESETS_FILE)                      # atomic: never a half-written file


def _sorted(items: list[dict]) -> list[dict]:
    return sorted(items, key=lambda x: str(x.get("name", "")))


def list_presets() -> list[dict]:
    """All presets sorted by name. LENIENT: missing/corrupt file -> [] (display must not explode)."""
    try:
        return _sorted(_load_strict())
    except Exception:                                  # noqa: BLE001 — corrupt -> show nothing
        return []


def save_preset(name: str, filters: dict, driven=None, thresh=None, delta=None) -> list[dict]:
    """Upsert a preset keyed on its (sanitised) name. Returns the updated list."""
    n = sanitize_preset_name(name)
    if not isinstance(filters, dict):
        raise ValueError("filters must be an object")
    rec = {
        "name": n,
        "filters": {str(k): str(v) for k, v in filters.items() if str(v).strip() != ""},
        "driven": driven if driven else None,
        "thresh": None if thresh is None else float(thresh),
        "delta": None if delta is None else float(delta),
        "saved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    with _LOCK:
        items = [x for x in _load_strict() if str(x.get("name", "")) != n]   # strict: may raise
        items.append(rec)
        _write(items)
        return _sorted(items)


def delete_preset(name: str) -> list[dict]:
    """Remove a preset by name (no-op if absent). Returns the updated list."""
    n = sanitize_preset_name(name)
    with _LOCK:
        items = _load_strict()                                              # strict: may raise
        kept = [x for x in items if str(x.get("name", "")) != n]
        if len(kept) != len(items):
            _write(kept)
        return _sorted(kept)
