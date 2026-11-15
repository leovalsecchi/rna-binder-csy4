# Filter presets — save a named filter set and reapply it from a list

**Date:** 2026-07-30
**Status:** Approved (design)

## Problem

Setting up a useful filter in the metric explorer means typing several hard filters by hand
(`apo_rmsd ≤2`, `hbonds ≥3`, …) and picking a driven metric. That state is lost on reload and has
to be retyped every time. The user wants to **save a filter under a name and reapply it from a
list**.

## Requirements (locked with the user)

1. **A preset captures the filters + the driven metric** — the typed hard-filter values, plus
   which metric is the ● driven/sort metric and its threshold and Δ.
2. **Stored server-side** in a JSON file in the app dir (survives browser cache clears, works
   across browsers, unit-testable) — consistent with `nas_collections` / `local_download`.
3. **Applying a preset REPLACES the filter state** — fields the preset doesn't mention are
   blanked, so a preset is a complete filter state, not an additive overlay.
4. **Two controls in the header** (next to the dataset / folder / download controls): a dropdown
   to apply, and a button opening a small popover to save/delete.

## Non-goals (YAGNI)

- No preset sharing/export, no per-dataset preset scoping, no ordering/favourites.
- Presets do NOT capture the dataset, the mode (compare/browse), or the 🟢 matchboxes.

## Architecture

### 1. `presets.py` (new module) — storage

- `PRESETS_FILE = APP_DIR / "presets.json"` — a JSON **list** of preset records.
- Record shape:
  ```json
  {"name": "self-consistent + H-bonded",
   "filters": {"apo_rmsd": "<=2", "hbonds": ">=3"},
   "driven": "apo_rmsd", "thresh": 2.0, "delta": 1.0,
   "saved_utc": "2026-07-30T09:12:00+00:00"}
  ```
  `filters` holds the **raw typed text** per metric key (exactly what sits in the `.confield`
  input), so operators (`>`, `>=`, `<`, `<=`, `=`, bare number) round-trip with no re-parsing.
- `sanitize_preset_name(name)` — trim, collapse inner whitespace, require non-empty, cap length
  (120). Raises `ValueError` otherwise. Path traversal is NOT a concern (a JSON key, not a
  filename), so the character set is permissive.
- `list_presets()` — records sorted by name. **Lenient**: a missing or corrupt file yields `[]`
  (display path must never explode).
- `save_preset(name, filters, driven, thresh, delta)` — upsert keyed on name (case-sensitive after
  sanitising), stamps `saved_utc`, returns the updated list. **Strict read**: a parse error on an
  existing `presets.json` PROPAGATES so a corrupt file aborts the write rather than replacing all
  presets with one record.
- `delete_preset(name)` — removes it, returns the updated list; strict read for the same reason.
- All writes atomic (`tmp` + `os.replace`) and serialised under a `threading.Lock`.
- `filters` values are coerced to `str`; non-dict input rejected.

### 2. `app.py` — three endpoints

- `GET  /api/presets` → `{presets: [...]}`
- `POST /api/presets/save` → body `{name, filters, driven, thresh, delta}` → `{ok, presets}`;
  `ValueError` → `{ok:false, error}` 400.
- `POST /api/presets/delete` → body `{name}` → `{ok, presets}`.

Save/delete return the refreshed list so the UI needs no follow-up GET.

### 3. `templates/index.html` — UI

- Header **`⚙ preset…` `<select id="preset-select">`** (mirrors `folder-select`): lists saved
  presets; `onchange` applies immediately, then resets its label back to the placeholder.
- Header **`💾 Presets` button** → popover `#preset-pop` reusing the existing `.folder-pop` CSS:
  a name input + **Save current**, and the saved list, each row with a **✕** delete.
- `applyPreset(p)`:
  1. Every `.confield` gets `p.filters[key] ?? ""` (blanking unmentioned metrics — requirement 3)
     and its `bad` class cleared.
  2. If `p.driven` exists in the current dataset's `METRICS`: clamp `p.thresh` into
     `[m.min, m.max]` → `THRESH[driven]` + set the metric's range input value and `val-<key>`
     label; `DELTA[driven] = p.delta`; check its radio; call the existing **`setDriven(driven)`**
     (which reconfigures the Δ slider and calls `updateRule/updateFilters/scheduleLoad`).
  3. Otherwise keep the current driven metric and just `INDEX=0; updateFilters(); scheduleLoad()`.
  4. Report skipped keys (not in this dataset) in the popover status line.
- `saveCurrentPreset()` collects `{key: input.value}` for all non-empty `.confield`s plus
  `DRIVEN`/`THRESH[DRIVEN]`/`DELTA[DRIVEN]` and POSTs it.

## Edge cases

- **Dataset metric mismatch** (30k exposes 13 metrics, 100k 22): apply only keys present in the
  loaded dataset; surface the skipped ones ("applied 2 · skipped apo_rmsd — not in 30k"). A saved
  `driven` that's absent leaves the current driven metric untouched.
- **Threshold out of range** for the current dataset (2of3 shifts distributions): clamped to the
  metric's min/max.
- **Duplicate name**: overwrites, status says "updated".
- **Empty name**: inline error, nothing saved.
- **Preset with zero filters** is legal — a handy "clear all filters" preset.
- **Corrupt `presets.json`**: the dropdown shows nothing; a save/delete aborts with an error
  rather than destroying the file's other records.

## Testing (`test_presets.py`, plain asserts like `test_nas_collections.py`)

Point `PRESETS_FILE` at a temp path:

- save → list round-trip; `filters` text preserved **verbatim** (e.g. `"<=2"` not `2`).
- upsert: saving the same name twice replaces (list length stays 1) and updates values.
- delete removes only the named preset.
- name validation: rejects empty/whitespace-only; trims; rejects over-long.
- corrupt `presets.json` → `list_presets()` returns `[]` but `save_preset` RAISES (no data loss).
- `driven`/`thresh`/`delta` persist, including `driven=None`.
