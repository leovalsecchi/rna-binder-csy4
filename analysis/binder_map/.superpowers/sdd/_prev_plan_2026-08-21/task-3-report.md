# Task 3 report — Canonical RNA frame and length-variant alignment

Branch: `foldseek-binder-map`
Commit: `f73a86a` — "Add canonical RNA frame mapping with length-variant alignment"

## What was implemented

Followed `.superpowers/sdd/task-3-brief.md` steps 1–5 in order, transcribing
both the test file and the implementation module **verbatim** from the
brief, using true TDD (test written and confirmed red — `ModuleNotFoundError`
— before `src/binder_map/rnaframe.py` existed).

Files created:
- `tests/test_rnaframe.py` — 10 collected test cases: 3 standalone
  (`frame_offset` exact match, closing-pair offset, rejection of an
  unrelated sequence via `FrameError`), 5 from the parametrized
  `test_every_fixture_aligns_to_canonical_frame` (one per fixture:
  `rfd3_design`, `bd1_design`, `bgbif_design`, `af3_refold`, `esm_refold`),
  and 2 covering `canonical_positions` (full 1..16 coverage on
  `af3_refold`, and the out-of-frame path on `bd1_design`, asserted
  **unconditionally** as instructed — no `if len(seq) > RNA_LEN` guard).
- `src/binder_map/rnaframe.py` — `FrameError`, `rna_sequence`,
  `_residue_index` (private helper), `frame_offset`, `canonical_positions`.
  Residues are ordered by first appearance in `(chain_id, res_id, ins_code)`
  coordinate order, not by `res_id` value, per the brief's note that
  `res_id` is not guaranteed to start at 1 across dialects. The canonical
  16-mer is located by exact substring match (`str.find`); flanking
  nucleotides outside the match window are reported as canonical position
  `0` rather than being renumbered or dropped.

## Test commands run, with actual output

### Step 2 — confirm RED (module doesn't exist yet)

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_rnaframe.py -v
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 0 items / 1 error

==================================== ERRORS ====================================
___________________ ERROR collecting tests/test_rnaframe.py ____________________
ImportError while importing test module '/Users/leo/Documents/binder_map/tests/test_rnaframe.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/opt/homebrew/Cellar/python@3.10/3.10.21/Frameworks/Python.framework/Versions/3.10/lib/python3.10/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_rnaframe.py:6: in <module>
    from binder_map.rnaframe import (
E   ModuleNotFoundError: No module named 'binder_map.rnaframe'
=========================== short test summary info ============================
ERROR tests/test_rnaframe.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.34s ===============================
```

Matches the brief's expected failure exactly:
`ModuleNotFoundError: No module named 'binder_map.rnaframe'`.

### Step 4 — confirm GREEN (after writing `src/binder_map/rnaframe.py`)

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_rnaframe.py -v
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 10 items

tests/test_rnaframe.py ..........                                        [100%]

=============================== warnings summary ===============================
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[rfd3_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[bd1_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[bgbif_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[af3_refold]
tests/test_rnaframe.py::test_canonical_positions_cover_all_16
tests/test_rnaframe.py::test_canonical_positions_marks_out_of_frame_as_zero
  /Users/leo/.venvs/binder_map/lib/python3.10/site-packages/biotite/structure/io/pdbx/convert.py:435: UserWarning: Attribute 'auth_comp_id' not found within 'atom_site' category. The fallback attribute 'label_comp_id' will be used instead
    warnings.warn(

tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[rfd3_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[bd1_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[bgbif_design]
tests/test_rnaframe.py::test_every_fixture_aligns_to_canonical_frame[af3_refold]
tests/test_rnaframe.py::test_canonical_positions_cover_all_16
tests/test_rnaframe.py::test_canonical_positions_marks_out_of_frame_as_zero
  /Users/leo/.venvs/binder_map/lib/python3.10/site-packages/biotite/structure/io/pdbx/convert.py:435: UserWarning: Attribute 'auth_atom_id' not found within 'atom_site' category. The fallback attribute 'label_atom_id' will be used instead
    warnings.warn(

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================= 10 passed, 12 warnings in 0.35s ========================
```

10 collected, 10 passed. The `auth_comp_id`/`auth_atom_id` fallback warnings
are the ones flagged as expected/harmless in the task context — not
suppressed, per instruction.

### Full suite (post-implementation, after commit)

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest -v
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
testpaths: tests
collected 20 items

tests/test_chains.py .......                                             [ 35%]
tests/test_constants.py ...                                              [ 50%]
tests/test_rnaframe.py ..........                                        [100%]

======================= 20 passed, 24 warnings in 0.44s ========================
```

All pre-existing Task 1/2 tests (`test_constants.py`, `test_chains.py`)
remain green alongside the new `test_rnaframe.py`.

## Deviation from the brief

**Test count only** — the brief's Step 4 says "Expected: 9 passed"; the
actual, correct collection count is **10**. Counting the test file's
bodies: 3 standalone tests (`test_frame_offset_exact_match`,
`test_frame_offset_with_closing_pair`,
`test_frame_offset_rejects_unrelated_sequence`) + 5 expansions of the
parametrized `test_every_fixture_aligns_to_canonical_frame` (one per
fixture name) + 2 more (`test_canonical_positions_cover_all_16`,
`test_canonical_positions_marks_out_of_frame_as_zero`) = 10. This mirrors
the identical off-by-one miscount noted in the Task 2 report ("7 items
collected" vs. actual parametrize expansion) — a prose-counting slip in
the brief, not a code or test defect. No test file or implementation code
was changed to work around it; the code and tests are transcribed verbatim
as instructed, and all 10 pass.

No other deviations. The implementation module and test file were both
copied verbatim from the brief with no modifications.

## Self-review findings

1. **Ground-truth cross-check.** Ran `rna_sequence`/`frame_offset` directly
   against all 5 fixtures outside the test suite and compared to the
   task's stated ground truth table:
   - `af3_refold`, `esm_refold`, `bgbif_design`, `rfd3_design`: 16 nt,
     `ACUGCCGUAUAGGCAG`, offset 0 — matches exactly.
   - `bd1_design`: 18 nt, `CACUGCCGUAUAGGCAGC`, offset 1 — matches exactly.

2. **`canonical_positions` on `bd1_design`** (the only length-variant
   fixture): 383 atoms total, exactly 40 atoms land at position `0`
   (the two flanking nucleotides — a leading `C` and a trailing `C`), and
   `set(pos[pos > 0]) == {1..16}` — every canonical position covered
   exactly once. No off-by-one at either the offset boundary or the tail
   boundary.

3. **No guard reintroduced.** Confirmed
   `test_canonical_positions_marks_out_of_frame_as_zero` asserts
   `(pos == 0).any()` unconditionally, with no `if len(seq) > RNA_LEN`
   wrapper — matches the brief's explicit instruction not to make this
   assertion silently vacuous.

4. **`res_id`-independent ordering.** `_residue_index` keys residues by
   `(chain_id, res_id, ins_code)` in first-appearance (coordinate) order,
   not by sorting on `res_id`, as required — verified this matters in
   practice since `rfd3_design`'s RNA chain is chain `A` with the protein
   inverted onto chain `B` (per Task 2), and residue numbering dialects
   differ across producers (BG+BIF uses a distinct mmCIF dialect per
   `test_chains.py`).

5. **`frame_offset` ambiguity.** Uses `str.find`, which returns the first
   match. `RNA_CANONICAL = "ACUGCCGUAUAGGCAG"` has no internal repeat that
   would make it self-overlap as a substring of a longer sequence, and
   none of the 5 fixtures produce more than one occurrence, so first-match
   semantics are safe here. Flagged as a latent assumption for future
   fixtures with repetitive flanking sequence, but out of scope for this
   task (brief specifies exact substring match with no disambiguation
   requirement).

6. **Import surface.** `rnaframe.py` consumes `RNA_CANONICAL`/`RNA_LEN`
   from `binder_map.__init__` and `split_chains`/`load_structure` from
   `binder_map.chains`, exactly as scoped in the brief's Interfaces
   section — no reimplementation of chain handling.

No defects found. Task 3 is complete and matches the brief's specification
and verified ground truth.
