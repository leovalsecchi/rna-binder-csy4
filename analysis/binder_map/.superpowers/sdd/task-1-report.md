# Task 1: Backbone Universe Rule — Implementation Report

## Summary

Implemented the per-stack design→backbone mapping module as specified in the task brief. The implementation correctly handles both sequence-variant stacks (RFD3+LiMPNN, BG+LiMPNN, BD1+LiMPNN) and one-backbone-per-design stacks (BG+BIF, HD+LiMPNN).

## What Was Implemented

### Files Created
- `src/binder_map/universe.py` — Core module with:
  - `ONE_BACKBONE_PER_DESIGN` frozenset: identifies stacks where each design is its own backbone
  - `SEQ_VARIANT_STACKS` frozenset: identifies stacks where _seq variants collapse to one backbone
  - `backbone_key(design_id: str, stack: str) -> str` function: maps designs to backbone keys

### Files Modified
- `tests/test_universe.py` — Test suite with 8 test cases covering:
  - Sequence suffix stripping for LiMPNN stacks
  - Sequence variant collapsing behavior
  - One-backbone-per-design behavior
  - Idempotency of the key function
  - Correct stack definitions
  - Unknown stack rejection

## TDD Evidence

### RED Step
```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
PYTHONPATH=src:scripts /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python -m pytest tests/test_universe.py -v'
```

**Output:**
```
ERROR collecting tests/test_universe.py
...
ModuleNotFoundError: No module named 'binder_map.universe'
```

Expected failure confirmed: module does not exist yet.

### GREEN Step
```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
PYTHONPATH=src:scripts /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python -m pytest tests/test_universe.py -v'
```

**Output:**
```
tests/test_universe.py ........                                          [100%]

============================== 8 passed in 2.27s ===============================
```

All 8 tests pass.

## Step 5: Real Data Verification

**Output:**
```
rule keys 31767  true backbones 31763  over-count 4
PASS: over-count is exactly the 4 known cross-prefix backbones
```

Assertions pass. Rule over-counts by exactly 4 as expected.

## Full Test Suite Baseline

**Output (final lines):**
```
tests/test_universe.py ........                                          [ 94%]
tests/test_upstream_parity.py ssssss                                     [100%]

=========== 97 passed, 12 skipped, 112 warnings in 60.07s (0:01:00) ============
```

97 passed (89 baseline + 8 new), 12 skipped, 0 failed. No regressions.

## Commit

```
commit 09dab90
Author: leovalsecchi <leovalsecchi@users.noreply.github.com>

    Add per-stack design->backbone rule

    BG+BIF and HD are one backbone per design; the three LiMPNN stacks dedup by
    stripping _seq<N>. Verified against the v6 coord-hash ground truth: the rule
    over-counts by exactly 4 on the 55,706-design pool.
```

## Files Changed

- `src/binder_map/universe.py` — NEW (54 lines): core mapping logic
- `tests/test_universe.py` — NEW (44 lines): 8 test cases

## Self-Review

### Completeness
- ✓ Both files created and committed
- ✓ All 8 test cases from brief implemented exactly
- ✓ Implementation matches brief specification
- ✓ Real-data verification passes with expected over-count of 4
- ✓ Full test suite passes with no regressions

### Code Quality
- ✓ Module docstring documents the per-stack behavior with ground-truth measurements
- ✓ Clear error message for unknown stacks (raises KeyError with stack name)
- ✓ Regex pattern for sequence suffix is simple and correct
- ✓ Function behavior is idempotent as verified by test
- ✓ Constants are module-level frozensets as specified

### Testing
- ✓ TDD workflow followed: RED → GREEN
- ✓ All test cases pass
- ✓ Real-data verification confirms expected behavior
- ✓ Baseline tests remain passing
- ✓ No overbuilding: implements exactly what was specified

### Naming and Clarity
- ✓ Function and constant names match brief exactly
- ✓ Implementation is minimal (no extra methods, no over-generalization)
- ✓ Docstring explains the why, not just the what

## Concerns

None. The implementation is complete, tested, verified against real data, and commits cleanly without regressions.
