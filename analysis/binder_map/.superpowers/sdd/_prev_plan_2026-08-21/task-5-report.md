# Task 5 report: epitope fingerprint

## What was implemented

- `src/binder_map/fingerprint.py` — verbatim from the brief:
  - `moiety_of(atom_name) -> "phosphate"|"sugar"|"base"`, classifying by a
    fixed set of phosphate atom names (`P, OP1, OP2, OP3, O1P, O2P, O3P, O5', O3'`)
    and sugar atom names (`C1', C2', C3', C4', C5', O2', O4'`); everything else
    (ring/exocyclic atoms) falls through to `"base"`.
  - `Fingerprint` frozen dataclass: `vector` (48, int32), `per_nt` (16, int32),
    `offframe: int`, `total: int`.
  - `fingerprint(protein, rna) -> Fingerprint`: filters both chains to
    heavy atoms (`element != "H"`), builds a `cKDTree` over protein heavy
    atoms, and for every RNA heavy atom counts protein neighbours within
    `RNA_CONTACT_CUTOFF = 4.5` Å with `query_ball_point`. Per-atom neighbour
    counts (not neighbour presence) are summed — this preserves
    protein-atom×RNA-atom **pair** semantics, not "atoms with ≥1 contact".
    Counts are split by `canonical_positions(rna)` into in-frame (binned into
    the 48-dim vector by nucleotide × moiety) and out-of-frame (summed into
    `offframe` only, never touching `vector`).
  - `FINGERPRINT_COLUMNS`: 48 names `nt{i:02d}_{moiety}`, `i` outer loop
    1..16, moiety inner loop in `MOIETIES` order (`phosphate, sugar, base`).
- `tests/test_fingerprint.py` — the 8 test functions (10 cases with
  parametrization) from the brief, verbatim.

No code deviated from the brief's Step 3 listing. `RNA_CONTACT_CUTOFF`,
`RNA_LEN`, `MOIETIES`, `FINGERPRINT_DIM` were consumed from `binder_map`
unchanged; nothing in `src/binder_map/__init__.py` was touched.

## TDD sequence

### Step 2 — red run (before `fingerprint.py` existed)

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
```

Actual output:

```
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 0 items / 1 error

==================================== ERRORS ====================================
__________________ ERROR collecting tests/test_fingerprint.py __________________
ImportError while importing test module '/Users/leo/Documents/binder_map/tests/test_fingerprint.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/opt/homebrew/Cellar/python@3.10/3.10.21/Frameworks/Python.framework/Versions/3.10/lib/python3.10/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_fingerprint.py:6: in <module>
    from binder_map.fingerprint import (
E   ModuleNotFoundError: No module named 'binder_map.fingerprint'
=========================== short test summary info ============================
ERROR tests/test_fingerprint.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.35s ===============================
```

Matches the brief's expected failure exactly (`ModuleNotFoundError`).

### Step 4 — green run, new tests only

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
```

Actual output (tail):

```
tests/test_fingerprint.py ..........                                     [100%]
...
7 warnings (auth_comp_id fallback) + 7 warnings (auth_atom_id fallback) — expected, not suppressed
======================= 10 passed, 14 warnings in 0.69s ========================
```

10 passed (8 test functions, `test_fingerprint_shape_and_consistency` parametrized
over 3 fixtures = 10 cases). Matches the brief's expectation exactly.

### Full suite after implementation

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest -v
```

Actual output:

```
tests/test_backbone.py .......                                           [ 18%]
tests/test_chains.py .......                                             [ 37%]
tests/test_constants.py ...                                              [ 45%]
tests/test_fingerprint.py ..........                                     [ 72%]
tests/test_rnaframe.py ..........                                        [100%]
======================= 37 passed, 52 warnings in 0.54s ========================
```

27 previously-passing tests + 10 new = 37, all green. No regressions.

## Observed `offframe` / `total` values per fixture

Computed directly (`fingerprint(protein, rna)` on `split_chains(load_structure(path))`
for every fixture, not just the 3 the brief parametrizes over):

| fixture | RNA len | `total` | `per_nt.sum()` | `offframe` | `vector.sum()` |
|---|---|---|---|---|---|
| rfd3_design | 16 | 63 | 63 | 0 | 63 |
| rfd3_design_seq7 | 16 | 63 | 63 | 0 | 63 |
| bgbif_design | 16 | 312 | 312 | 0 | 312 |
| bd1_design | 18 | 0 | 0 | 0 | 0 |
| af3_refold | 16 | 250 | 250 | 0 | 250 |
| esm_refold | 16 | 179 | 179 | 0 | 179 |

`test_fingerprint_shape_and_consistency` is parametrized over
`af3_refold` (total 250), `esm_refold` (total 179), `bd1_design` (total 0) —
all pass the consistency identity `total == per_nt.sum() + offframe`
regardless of magnitude, including the degenerate all-zero case.

## Deviation from the brief, and why

**No deviation in the shipped code or tests** — both files are the brief's
Step 1 and Step 3 listings verbatim, byte-for-byte.

One finding changes what the brief's own claim actually demonstrates, so it's
recorded here rather than silently accepted:

**`bd1_design`'s `offframe` is 0, not >0, on the committed fixture.** The
brief states bd1_design "is the ONLY fixture that can produce a non-zero
`offframe`" and its test docstring says contacts to the two flanking
nucleotides "must land in `offframe`". I measured the actual protein–RNA
geometry: the minimum protein-heavy-atom↔RNA-heavy-atom distance in
`bd1_design` is **4.503144 Å** — 0.0031 Å outside the 4.5 Å cutoff — and it is
one of the out-of-frame flanking-nucleotide atoms (`C5` of the offframe
residue closest to protein) that would be the first to make contact if
anything did. Because *no* protein–RNA pair in this structure is within
cutoff, `total = per_nt.sum() = offframe = 0` for this fixture: the identity
holds trivially, so the checked-in test still passes, but it does not
empirically exercise a non-zero `offframe` value. This is a fact about the
fixture's geometry (right at the boundary), not a bug in the implementation
or the test.

I did not alter the cutoff, the fixture, or the test to force a non-zero
value — per instructions the cutoff/binning are pinned and must not be
"improved," and the brief's test code was to be used verbatim.

To confirm the `offframe`-vs-`vector` split logic itself is correct (rather
than just algebraically self-consistent at zero), I ran an out-of-band
manual check (not added to the test suite): nudging `bd1_design`'s RNA
coordinates by ~5.5 Å toward the protein along the vector to its nearest
out-of-frame atom produces `total=373, per_nt.sum()=367, offframe=6,
vector.sum()=367` — i.e. offframe contacts are correctly counted in `total`
and `offframe` and correctly excluded from `vector`/`per_nt`. A second nudge
toward the nearest *in-frame* atom instead produced `total=2, offframe=0,
vector.sum()=2`, landing at canonical position 10 as expected. This
confirms the in-frame/out-of-frame partition (`positions > 0` mask) is wired
correctly; it just isn't exercised as non-zero by the committed fixture set.
Recommendation for whoever owns fixture curation: if a non-zero-`offframe`
regression test is wanted, either pick/generate a fixture with true offframe
contact (bd1_design is right at the edge, so a different bd1-family design or
a slightly relaxed structure would do it), or accept that this path is
covered by code review + the manual check above rather than by CI.

## Self-review findings

1. **Contact = pair, not atom-with-≥1-neighbour, confirmed correct.**
   `counts = [len(n) for n in neighbours]` sums *all* protein atoms within
   cutoff per RNA atom, then those counts (not 0/1 indicators) are what get
   added into `vector`/`offframe`/`total`. This matches the pinned semantics
   from the task instructions.
2. **Atom-name vocabulary coverage checked against all 6 fixtures**, not
   just the 3 the parametrized test touches. The union of RNA atom names
   across `rfd3_design`, `rfd3_design_seq7`, `bgbif_design`, `bd1_design`,
   `af3_refold`, `esm_refold` is exactly 28 names: 6 classify as phosphate,
   7 as sugar, and the remaining 15 (`C2 C4 C5 C6 C8 N1 N2 N3 N4 N6 N7 N9 O2
   O4 O6`) correctly fall through to base. No atom name is silently
   misclassified or missing from either explicit set.
3. **No hydrogens in any fixture** (`"H" not in set(atoms.element)` for both
   chains, all fixtures), so `_heavy()` is a no-op here specifically, but the
   filter is still correct/needed for general applicability (AF3/ESM/other
   producers can emit explicit H).
4. **Determinism**: `query_ball_point`/`np.add.at` give bit-identical
   results across repeated calls on the same input (`test_fingerprint_is_deterministic`
   passes; also reconfirmed manually).
5. **Distant-chain zero case**: translating RNA 500 Å away gives `total=0,
   vector.sum()=0` as required — `cKDTree` with an empty ball query per point
   returns `[]`, `len([])=0`, no exception, no shape bugs.
6. No changes were needed to `binder_map/__init__.py`, `chains.py`, or
   `rnaframe.py` — Task 5 only reads from them, as scoped.

## Process note (not a code issue): commit landed inside an unrelated commit

Per Step 5, I staged exactly `src/binder_map/fingerprint.py` and
`tests/test_fingerprint.py` and went to commit them. By the time I ran
`git add`, `git status` was already clean — the two files had already been
committed as part of `6c88ebe "Remove dead line from Task 7 kNN ranking
code"`, alongside an edit to
`docs/superpowers/plans/2026-08-21-foldseek-binder-map.md` (a `_topk_per_node`
ranking-comment cleanup) that I did not make and is unrelated to Task 5.

This happened because another process is concurrently working in this same
(non-worktree) checkout on later tasks (Task 7's kNN ranking code, per the
plan-file edit and commit message) and, at some point after I wrote
`fingerprint.py`/`test_fingerprint.py` to disk but before I staged them
myself, that process ran a broad `git add` and commit that swept up my
still-untracked files along with its own.

I did not rewrite history to fix this: amending or resetting a commit that
may belong partly to a still-active concurrent process risks colliding with
its next commit (diverged SHA, lost work, or a confusing force-push
situation), and the task's git-safety rules bar amend/reset without explicit
instruction. The commit exists, is on `foldseek-binder-map`, and its content
is byte-identical to what this report describes and to a normal `git show
6c88ebe -- src/binder_map/fingerprint.py tests/test_fingerprint.py` — so
nothing is lost — but the commit message does not mention Task 5/epitope
fingerprint at all, and `git log --oneline -- src/binder_map/fingerprint.py`
returns only that one, misleadingly-titled commit. Flagging this so whoever
integrates the branch is aware and can decide whether to leave it, split it
in a follow-up rebase, or just note it in the PR description.

**Verification commands used:**
```
git status --short
git log --oneline -5
git show --stat 6c88ebe
git show 6c88ebe --name-only
git show 6c88ebe -- docs/superpowers/plans/2026-08-21-foldseek-binder-map.md
```

## Bottom line

Implementation matches the brief exactly; all 10 new tests plus the full
37-test suite pass; the pinned cutoff/binning were not touched. The one
substantive finding is that `bd1_design`'s offframe path is algebraically
tested but not empirically non-zero on the current fixture (verified by
direct geometry measurement and a manual synthetic-contact check), and the
commit ended up folded into an unrelated concurrent commit rather than its
own clean "Add 48-dim epitope fingerprint..." commit.

## Fix: non-vacuous offframe test

`test_offframe_counts_contacts_outside_the_canonical_frame` (flagged above as
algebraically-trivial-at-zero) was replaced with a version that asserts the
real geometry first, then forces the out-of-frame path by moving only the
out-of-frame atoms into contact.

**Before (original test, passing but vacuous):**
```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
============================= test session starts ==============================
collected 10 items

tests/test_fingerprint.py ..........                                     [100%]

======================= 10 passed, 14 warnings in 0.67s ========================
```
(All 10 pass, but `offframe == 0` throughout — the assertions
`fp.total == int(fp.per_nt.sum()) + fp.offframe` and
`fp.vector.reshape(RNA_LEN, 3).sum() == fp.per_nt.sum()` hold trivially at
zero and never exercise the offframe accumulation path.)

**Change applied:**
- Import changed to `from binder_map.rnaframe import canonical_positions, rna_sequence`.
- Test body rewritten to: split chains, assert `len(rna_sequence(rna)) == 18`;
  compute `baseline = fingerprint(protein, rna)` and assert
  `baseline.offframe == 0` (documenting the real 4.503 Å / 4.5 Å near-miss);
  get `positions = canonical_positions(rna)` and assert at least one atom has
  `positions == 0`; copy the RNA, move only the out-of-frame atoms
  (`moved.coord[positions == 0] = protein.coord[0]`) directly onto a protein
  atom; recompute `fp = fingerprint(protein, moved)`; assert `fp.offframe > 0`;
  assert `np.array_equal(fp.vector, baseline.vector)` (the 48-dim vector is
  untouched, proving the new contacts landed in `offframe` and not in-frame);
  assert `fp.total == int(fp.per_nt.sum()) + fp.offframe`.
- `src/binder_map/fingerprint.py` was not modified, per instructions.

**After (corrected test):**
```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
============================= test session starts ==============================
collected 10 items

tests/test_fingerprint.py ..........                                     [100%]

======================= 10 passed, 14 warnings in 0.64s ========================
```
All 10 tests pass, including
`test_offframe_counts_contacts_outside_the_canonical_frame` run in isolation
(`pytest tests/test_fingerprint.py::test_offframe_counts_contacts_outside_the_canonical_frame -v`
→ `1 passed`). `fp.offframe > 0` was genuinely achieved this time — the
implementation correctly routes the forced out-of-frame contacts into
`offframe` while leaving `vector` byte-for-byte equal to the baseline.

**Full suite after the fix:**
```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest -q
37 passed, 52 warnings in 0.82s
```

**Commit:** `891da3d9bae28d4b68d137fcfc8bed78ed7144c5` — "Fix vacuous offframe
test to assert real out-of-frame geometry" — `tests/test_fingerprint.py` only
(staged with `git add tests/test_fingerprint.py`, not `-A`/`.`).

**Bottom line:** the test now genuinely exercises and confirms the offframe
accumulation path; no implementation change was needed — `fingerprint()` in
`src/binder_map/fingerprint.py` already routes out-of-frame contacts
correctly.
