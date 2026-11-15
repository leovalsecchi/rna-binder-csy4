# Metric Explorer — Collections ("folders") feature

**Date:** 2026-07-23
**Status:** approved (design), ready to implement
**App:** `~/Documents/metric_explorer` (local Flask + 3Dmol single-page tool)

## Goal

Let the user, while looking at a design in the explorer, drop it into one or more named
**folders** (collections). Each folder is a curated shortlist that materializes on the
**Wicky NAS** as:

- `metrics.csv` — the full metric record (all 230 dataset columns) for every design in the
  folder, one row each, plus provenance columns.
- `cif/` — the three decompressed refold structures (AF3 / Boltz-2 / ESMFold2) for every
  design, under their original archive-member names.

Folders are selectable from a dropdown; multiple folders can be ticked at once. Adding is
one click and feels instant; the heavy work (CIF extraction + NAS write) happens in a
background thread.

## Where things live (decided)

- **NAS is mounted on the Mac** at `/Volumes/wicky` (verified live). The collections root is
  a real local path:

  ```
  NAS_ROOT = /Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_Coding
  ```

  Writing there **is** writing to the NAS. No gobi, no Euler staging, **no change to the
  nightly `sync_euler_to_nas.sh`** — it is untouched.
- **Euler is read-only.** The only Euler interaction is streaming/decompressing CIF members
  out of the `*.tar.zst` archives via `ssh euler` (the exact mechanism `app.py:fetch_cif`
  and `open_refolds_server.py` already use). Nothing is written on Euler.
- **The folders on the NAS are the source of truth.** The dropdown is built by listing
  subdirectories of `LV_ROOT`; membership/dedup come from reading each folder's
  `metrics.csv`. There is **no separate local registry file** to drift.

## Folder layout on disk (NAS)

```
LV_Coding/
  <folder>/
    metrics.csv
    cif/
      <design_id>_AF3_seed1_sample0.cif      # original member basename
      <design_id>_B2_sample0.cif
      <design_id>_ESM_sample0.cif
      ...(3 per design; whichever backends have a ref)
```

- CIF filenames are the **original archive-member basenames** (they already encode
  `design_id` + backend + seed/sample and are unique across the three refolds).
- `metrics.csv` columns = the dataset's columns for that design **plus**: `note` (B),
  `added_utc`, `dataset`, `driven_metric` (C). Union-of-columns across designs; missing
  cells blank.

## Data facts (from the live dataset)

- `design_id` (str) uniquely identifies a design; member basenames start with `design_id`.
- Refs: `af3_msaoff_cif_path`, `b2_cif_path`, `esm_cif_path`, each
  `"<archive>.tar.zst!<member>.cif"` relative to `RB`. Some may be NaN (extract what exists).
- `_row_payload` already carries the dataframe row index (`row`); the frontend passes
  `{dataset, row}` and the backend reconstructs `design_id`, refs, and the full 230-col row
  from `_DATASETS[dataset]["df"].iloc[row]` — no need to ship 230 columns through JS.

## Components

### 1. `collections.py` (new module — all NAS/Euler logic, no Flask)

Pure-ish, independently testable. Key surface:

- `LV_ROOT: Path` and `nas_available() -> bool` — `os.path.ismount("/Volumes/wicky")` and
  `LV_ROOT.parent` exists.
- `sanitize_folder(name) -> str` — allow `[A-Za-z0-9._ -]`, strip, collapse whitespace,
  reject empty, `.`, `..`, and anything with a path separator. Raises `ValueError` on bad
  input. **All folder names from HTTP pass through this** (path-traversal guard).
- `list_folders() -> list[{name, count}]` — subdirs of `LV_ROOT` (skip dotfiles); `count` =
  data rows in that folder's `metrics.csv` (0 if none). Returns [] if NAS down.
- `folder_members(folder) -> set[str]` — `design_id`s in that folder (for tick state/dedup).
- `add_design(folders, design_id, refs, row_dict, extras) -> list[{folder, ok, msg}]`
  - for each folder: `mkdir -p <folder>/cif`;
  - extract each present ref's member → `<folder>/cif/<member_basename>` via
    `ssh euler 'cd RB && tar --zstd -xOf <archive> <member>'` streamed to the dest;
  - upsert the row into `<folder>/metrics.csv` (drop any existing row with same `design_id`,
    append the new one, union columns), **written atomically** (temp file + `os.replace`);
  - idempotent: re-adding overwrites cleanly.
- `remove_design(folder, design_id) -> {ok, msg}` — drop the `design_id` row from
  `metrics.csv` (atomic rewrite) and delete `<folder>/cif/<design_id>_*.cif`.
- `browse_folder(folder) -> list[row-payload]` — read `metrics.csv`; for each design return a
  payload shaped like `_row_payload` but whose AF3 CIF is served from the **local NAS file**
  (no Euler round-trip): `{design_id, stack, metrics{...}, local_cif: "<basename>"}`.
- `serve_cif_path(folder, name) -> Path` — resolve `<folder>/cif/<name>` with sanitization
  (basename only, must stay inside the folder), for the browse-mode viewer.

Concurrency: a module-level `threading.Lock` serializes `metrics.csv` read-modify-write per
process so two quick adds can't clobber the CSV.

### 2. Flask endpoints (in `app.py`)

- `GET  /api/collections?design_id=<id>` → `{available, folders:[{name,count,has}]}` where
  `has` = design already in that folder (drives tick state).
- `POST /api/collections/add` body `{dataset, row, folders:[...], note, driven_metric}` →
  reconstruct design from the df, call `add_design`, return per-folder results.
  **Decision (single, no background job/poll machinery):** the route does the work
  **synchronously** and returns the per-folder result; the *frontend* calls it with a
  non-blocking async `fetch` and shows a spinner → `saved ✓` / error. A single add is 3
  small `ssh euler` extracts (~1–3 s), so no job-id/status endpoint is needed. Flask runs
  `threaded=True`, so a slow add never freezes other requests.
- `POST /api/collections/remove` body `{design_id, folder}` → `remove_design`.
- `GET  /api/collections/browse?folder=<f>` → `{folder, designs:[...]}`.
- `GET  /api/collections/cif?folder=<f>&name=<basename>` → `send_file` the NAS cif.

`collections.py` imports `RB`/`SSH` (or they are passed in) so extraction matches `app.py`.

### 3. Frontend (`templates/index.html`)

- **Per-panel control:** in `.pbtns`, add `＋ Folder ▾`. Clicking opens a small popover:
  a checkable list of existing folders (ticked = design is in it), a `＋ New folder…` input,
  and an optional one-line **note** field. Ticking a folder → `add`; unticking → `remove`.
  A tiny status line shows `saving… / saved ✓ / NAS not mounted / error`.
- **Browse-a-folder (A):** a header control `📁 Folder ▾` listing folders. Selecting one
  enters a **folders view**: the three panels page through that folder's designs (reuse the
  existing panel render + Next/Back), loading each AF3 CIF from
  `/api/collections/cif` (local NAS, no Euler). Leaving the view returns to compare/browse.
- Popover state refreshes when the shown design changes (re-query `has`).

## Error handling

- NAS down (`/Volumes/wicky` not a mount): every collections call returns
  `{available:false}`; UI shows `NAS not mounted` and disables the folder controls. No crash.
- `ssh euler` extract fails for a backend: that CIF is skipped; the row is still written; the
  per-folder result notes which backends were missing. A design with **zero** extractable
  CIFs still records its metrics row (metrics are the point; structures are best-effort).
- Interrupted CSV write can't corrupt the folder: temp-file + `os.replace` is atomic.
- Bad/hostile folder name → `sanitize_folder` raises → 400, no filesystem escape.

## Testing

- **Unit (pure logic, no NAS/Euler):** `sanitize_folder` (valid, empty, `..`, slashes,
  unicode), metrics-row upsert/dedup (same design_id replaces; column union across differing
  rows), remove-by-prefix filename matching (`seq5` must not match `seq50`), ref → member
  basename parsing. Run against a `tmp_path` standing in for `LV_ROOT`.
- **Live smoke (NAS + ssh euler are up now):** create a throwaway folder
  `__smoketest__`, add one real design, assert `metrics.csv` has 1 row with the extra
  columns and `cif/` has the expected files, re-add (assert still 1 row = idempotent),
  remove (assert row + cifs gone), delete the folder. Never touch a real user folder.

## Out of scope (YAGNI)

- No gobi / Euler-staging / nightly-sync changes.
- No push-to-Euler, no on-demand gobi trigger.
- No auth, no multi-user, no server-side rendering of B2/ESM in browse (AF3 only, matching
  the main viewer).

## Non-goals confirmed with user

- Do **not** write anything under `euler_work/` on the NAS.
- Do **not** modify `sync_euler_to_nas.sh`.
- Extras chosen: **A** (browse-a-folder), **B** (note per add), **C** (provenance columns).
  Remove-from-folder is core toggle behavior.

## Post-review hardening (2026-07-23)

A 5-dimension adversarial review (14 confirmed findings) drove these additions beyond the
original design; all covered by tests:

- **Security — ref confinement.** `/api/cif?ref=` was a pre-existing arbitrary-Euler-file
  read (`?ref=/etc/passwd`); `shlex.quote` blocks shell-metachar injection but not path
  escape or tar/cat option-injection. Added `nas_collections.validate_ref` (reject absolute /
  `~` / leading-`-` / `..` escape) + a `--` end-of-options guard, applied in **both**
  `fetch_cif` and `_extract_member`, plus `ssh` timeouts. CORS narrowed from `*` to
  localhost origins.
- **No silent data loss.** `_read_csv` (lenient, None on any failure) is used only for
  display; the upsert path uses `_read_existing` (**strict** — a parse error on an existing
  file propagates), so a corrupt/transiently-unreadable `metrics.csv` aborts the add instead
  of overwriting the folder down to one row.
- **Lock scope + timeouts.** CIF extraction runs outside `_CSV_LOCK`; only the CSV
  read-modify-write is serialized, so a slow `ssh` can't block all collection ops.
- **Robust error surface.** `add_design` catches `Exception` per folder (a NAS `OSError`
  becomes a failed-folder result, not HTTP 500); `list_folders` skips names `sanitize_folder`
  would reject and counts rows cheaply (`usecols=[design_id]`).
- **Empty-folder cleanup + UX.** removing the last design returns `emptied:true`, tidies the
  now-empty dir, and the UI leaves folder view; popover success status is no longer wiped by
  the follow-up refresh; match checkboxes are disabled in folder view.
