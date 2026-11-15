# Download-selected — bulk-copy the filtered design set to a named local folder

**Date:** 2026-07-29
**Status:** Approved (design)

## Problem

The metric explorer lets you filter the design set with hard-filter locks (🔴 conbox).
There is no way to pull that filtered set out of the tool. The user wants a **"Download all
selected"** action that copies every design passing the current filters into a **local folder
they name**, so a filtered shortlist becomes files on disk.

This is the existing NAS "collections" feature (`nas_collections.py`) generalised: same
per-design payload, but written to a **local** directory and acting on a **whole filtered set**
at once instead of one design.

## Requirements (locked with the user)

1. **Selected set = everything passing the current filters.** The same set `/api/browse` pages
   through — i.e. all designs passing the active 🔴 hard-filter locks. Resolved with the SAME
   `_parse_constraints` / `_apply_constraints` code path as Browse so the two never diverge.
2. **Per design, copy all 3 refolds + a metrics.csv.** AF3 + Boltz-2 + ESMFold2 CIFs into
   `<folder>/cif/`, plus a `<folder>/metrics.csv` holding the full metric rows of the selected set.
3. **Destination = `metric_explorer/downloads/<name>/`**, where `<name>` is a user-typed folder
   name. `downloads/` is gitignored (CIFs must never be committed).
4. **Size handling = confirm + cap.** Show the exact passing count first. If `count ≤ CAP` (500),
   download all. If `count > CAP`, return a confirm prompt ("N designs → ~3N CIFs, proceed?");
   only start on explicit confirm. No accidental 10k-fetch run.

## Non-goals (YAGNI)

- No per-design manual checkbox selection (the filter set IS the selection).
- No arbitrary destination path (fixed base `downloads/`, custom leaf name only).
- No zip/archive — a plain folder is the deliverable.
- No writing to Euler or the NAS; Euler is read-only (CIF streaming), as today.

## Architecture

Three units, each independently testable.

### 1. `local_download.py` (new module) — the batch copy engine

- `DOWNLOADS_ROOT = APP_DIR / "downloads"`.
- `CAP = 500` (module constant; the count above which the caller must confirm).
- `sanitize_name(name)` — reuse the `nas_collections.sanitize_folder` rules (letters/digits/
  space/dot/underscore/hyphen; no path separators; no leading dot; reject `.`/`..`).
- **Job registry** `_JOBS: {job_id: {status, total, done, cifs, missing, folder, error, name}}`
  guarded by a `threading.Lock`. `status ∈ {running, done, error}`.
- `start(name, designs) -> job_id` — validates the name, creates `downloads/<name>/cif/`,
  registers a job with `total = len(designs)`, and spawns ONE background `threading.Thread`
  that runs `_run(job_id, folder, designs)`. Returns the job id immediately (non-blocking).
- `_run(...)` — for each design, extract its 3 refolds concurrently via a
  `ThreadPoolExecutor(max_workers=6)` calling `nas_collections._extract_member(ref, dest)`
  (already RB-confined, `--`-guarded, atomic, timeout-bounded). **Idempotent:** skip any CIF
  already present and non-empty on disk (re-running resumes). After each design, increment
  `done`/`cifs`, append any missing-refold notes. At the end, write
  `downloads/<name>/metrics.csv` from the selected rows (full columns via pandas), set
  `status=done`. Any fatal error → `status=error, error=<msg>`.
- `job_status(job_id) -> dict | None` — snapshot of the job dict (copy under lock).

`designs` is a list of `{"design_id": str, "refs": {af3,b2,esm}, "row": dict}` built by the
caller from the dataframe — the module never touches the dataframe (clean boundary).

### 2. `app.py` — two endpoints

- `POST /api/download/start` — body `{dataset, constrain, metric, name, confirm}`.
  1. Resolve passing rows: `_apply_constraints` over `df[REFOLD_COL["af3"]].notna()` (same as
     `/api/browse`), ordered by `metric` for a tidy metrics.csv / fetch order. NOTE: there is no
     top-N truncation — once confirmed, the FULL passing set is copied; `metric` only sets order.
  2. `count = len(idx)`. If `count == 0` → `{ok:false, error:"no designs match"}`.
  3. If `count > local_download.CAP` and not `confirm` → `{ok:true, need_confirm:true, count}`.
  4. Else build `designs` from the rows, call `local_download.start(name, designs)` →
     `{ok:true, job_id, total:count}`.
  Name errors (`ValueError`) → `{ok:false, error}` 400.
- `GET /api/download/status?job=<id>` → the job snapshot, or `{error:"unknown job"}` 404.

### 3. `index.html` — the UI

- A **"⬇ Download selected"** button in the header controls row (near the Browse / filter
  controls).
- Click → small popover: a **folder-name text field**, a **live count** of passing designs
  (from the existing `/api/browse` `total` for the current constraints), and a **Download** button.
- Download → `POST /api/download/start`. If `need_confirm`, show an inline confirm
  ("N designs (> 500) — pulls ~3N CIFs over ssh, can take a while. Proceed?" Proceed/Cancel;
  Proceed re-POSTs with `confirm:true`). Once a `job_id` is returned, poll
  `/api/download/status` every ~1 s and render a progress bar "Downloading X / N…".
- On `status==done`: "Saved to downloads/&lt;name&gt; — N designs, M CIFs (K missing)."
  On `status==error`: show the error.
- Empty passing set → the button is disabled with "0 designs".

## Data flow

```
UI (folder name + current constraints)
  → POST /api/download/start
     → resolve idx via _apply_constraints  (same as /api/browse)
     → count check / cap-confirm
     → build designs[] from df rows (design_id, refs, full row)
     → local_download.start() → background thread
          → ThreadPoolExecutor: nas_collections._extract_member(ref → downloads/<name>/cif/…)
          → write downloads/<name>/metrics.csv
  → UI polls GET /api/download/status?job=<id> → progress bar → summary
```

## Error handling / edge cases

- **Failed CIF fetch** (ssh error / missing member): recorded in the job's `missing`, never
  aborts the batch — mirrors `add_design`.
- **Folder name collision** with an existing `downloads/<name>`: merge/resume into it (idempotent
  skip of present CIFs); the summary still reports the full selected count.
- **Bad/empty/duplicate name**: rejected with a clear message (reused sanitize rules).
- **NAS not required**: this feature reads CIFs from Euler over ssh (like `fetch_cif`), so it
  works whether or not the NAS is mounted.
- **Concurrent downloads**: each `start` gets a fresh `job_id`; the registry is lock-guarded.

## Testing (`test_local_download.py`, plain-asserts like `test_nas_collections.py`)

Stub `nas_collections._extract_member` to drop a placeholder file (no real ssh):

- `sanitize_name` accepts good names, rejects `.`/`..`/path-separators/leading-dot/empty.
- Cap logic: `start` not called when `count > CAP` without confirm (tested at the endpoint layer
  via a small helper, or by asserting the branch).
- Idempotent skip: a second `start` into the same folder does not re-fetch present CIFs.
- `metrics.csv` written with one row per selected design and the expected columns.
- Job lifecycle: `status` goes `running → done`; `done == total`; `missing` captured when a
  ref is absent.

## Reuse / boundaries

- Reuses `nas_collections._extract_member`, `_parse_ref`, `validate_ref` (import) — no
  duplication of the ssh/security logic.
- `local_download.py` has no Flask/pandas-dataframe coupling beyond writing the final CSV; the
  endpoint owns dataframe→designs translation. Each unit testable in isolation.
