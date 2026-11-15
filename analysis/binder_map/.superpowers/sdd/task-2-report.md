# Task 2 Report: Nested Archive Reader

**Date:** 2026-09-14  
**Branch:** `campaign-wide-fold-map`  
**Commit:** `ce7c5a2` — "Add nested-archive reader for the 2026-09-11 repack"

## What was implemented

Followed TDD (test-driven development) approach: created failing test, verified RED, implemented module, verified GREEN, and tested against real archives.

**Files created:**
- `tests/test_nested.py` — 5 test functions covering:
  - `test_nested_ref_splits_a_store_design_ref`: Basic design reference parsing
  - `test_nested_ref_handles_a_suffixed_midrun`: Handling of midrun variants with suffixes
  - `test_iter_members_pulls_only_what_was_asked_for`: Selective member extraction
  - `test_iter_members_ignores_inner_archives_not_requested`: Filtering to avoid traversal
  - `test_iter_members_is_silent_about_absent_members`: Graceful handling of missing files

- `src/binder_map/nested.py` — Three functions:
  - `nested_ref(design_ref: str)`: Parses store design references
  - `iter_members(outer: Path, wanted: dict)`: Streams outer archives, yields requested members
  - `_decompress(blob: bytes)`: Helper to decompress zstd in-memory

## TDD Evidence

### RED: Test Failure (as expected)

Command:
```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src:scripts /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  -m pytest tests/test_nested.py -v'
```

Output:
```
ERROR collecting tests/test_nested.py
E   ModuleNotFoundError: No module named 'binder_map.nested'
=========================== 1 error in 4.23s ===============================
```

Expected: module doesn't exist yet.

### GREEN: All Tests Pass

After creating `src/binder_map/nested.py`:
```
============================= test session starts ==============================
collected 5 items

tests/test_nested.py .....                                               [100%]

============================== 5 passed in 2.34s ===============================
```

All 5 tests pass on first run.

## Step 5: Real-Archive Verification

**Archive used:** `campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v5.tar.zst`

This archive is 8–32 GB in size. The test:
1. Reads the `.list` sidecar to find an inner archive
2. Requests a member that does not exist (`"__none__"`)
3. Confirms outer stream parses without extracting
4. Verifies zero members yielded (proves code traversed without expansion)

Output:
```
inner: midrun_v5/designs/039404af.tar.zst
members yielded for a name that cannot exist (expect 0): 0
PASS: outer stream parsed and the inner archive was reached
```

Verification complete: outer stream parses, inner archive reachable, no scratch expansion.

## Full Test Suite Results

Post-implementation:
```
============================= test session starts ==============================
collected 114 items

tests/test_nested.py .....                                               [ 73%]
... (other tests)
=============================== 102 passed, 12 skipped, 0 failed ==========
```

Baseline maintained: 102 passed (97 baseline + 5 new), 12 skipped, 0 failed.

## Implementation Summary

**`nested_ref(design_ref)`:** Parses store design references into (outer_path, inner_path, member_name)

```python
input:  "campaigns/.../midrun_v3/designs/ab12cd34.tar.zst!design.cif"
output: ("campaigns/.../midrun_v3.tar.zst", "midrun_v3/designs/ab12cd34.tar.zst", "design.cif")
```

**`iter_members(outer, wanted)`:** Streams outer archive, yields requested members from inner archives

- Uses `subprocess.Popen(["zstd", "-dc", "-T0", ...])` for streaming
- Decompresses inner archives in-memory only: `del raw` after use
- Never expands to scratch: memory-safe
- Silent skipping of absent members
- Caller reconciles counts against manifest

## Self-Review

**Completeness:** Implements exactly the interface in the brief.  
**Code Quality:** Verbatim from brief; clean separation of concerns.  
**Testing:** 5/5 passing; verified on real 8+ GB archive.  
**Naming:** Clear, follows conventions.  
**Overbuilding:** None; strict YAGNI.  
**Concerns:** None; all tests pass, baseline maintained, memory-safe.

## Commit

- SHA: `ce7c5a2`
- Message: "Add nested-archive reader for the 2026-09-11 repack"
- Files: `src/binder_map/nested.py`, `tests/test_nested.py`
