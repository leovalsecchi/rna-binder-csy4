# Foldseek Binder Map — SDD progress ledger

Plan: docs/superpowers/plans/2026-08-21-foldseek-binder-map.md
Branch: foldseek-binder-map
Merge base: 6c72608

## Completed tasks
Task 1: complete (commits 6c72608..c3ae783, review clean — spec ✅, quality approved)
  Minor (deferred to final review): .gitignore scope bleed in c3ae783; Homebrew python@3.10
  is a host-level side effect (justified: llvmlite is arm64-only; contained, not shadowing PATH).
  Controller corrections: 18-nt RNA variant is BD1+LiMPNN not BG+BIF (ce65bae); plan prose
  said "five fixtures" but six are fetched — rfd3_design_seq7 is required by Task 4.
Task 2: complete (commits 1e953c7..f66dc99, review clean — spec ✅, quality approved)
  Reviewer verified the RFD3 guard by monkeypatching a chain-ID revert: it fails as intended.
  Minor (resolved in plan, not code): Task 2 Interfaces falsely listed RNA_LEN as consumed;
  chains scanned twice per structure — Task 13's hot loop now derives the chain id from
  split_chains' result instead of calling protein_chain_id separately (215k calls/layer).
  Minor (deferred): third-chain-by-exclusion — assessed low risk, raises ChainSelectError
  rather than misrouting, and matches the upstream rule we must not drift from.
Controller fix (all tasks): expected test counts in the plan ignored parametrization
  expansion. Recomputed: rnaframe 9->10, fingerprint 9->10, whole-suite 46->67.
Task 3: complete (commits 4ddf8f9..f73a86a, review clean — spec ✅, quality approved)
  Reviewer built adversarial descending-res_id fixtures; ordering logic is genuinely
  res_id-independent. bd1_design out-of-frame assertion confirmed unconditional.
  Minor (deferred): frame_offset uses str.find, so it takes the first match — safe
  because RNA_CANONICAL has no internal self-overlap, but latent if a future target
  has repetitive flanking sequence.
  Accepted debt (deferred to final review): backbone.py's _residue_ordinal (Task 4)
  will duplicate rnaframe.py's _residue_index ordering logic (~10 lines). Not
  refactored: Task 3 is already committed and approved, and the churn exceeds the gain.
Task 4: implemented (commit 7ef3f4e) — Foldseek accepted the generated PDB
  (createdb 1 structure, 0 skipped; convert2pdb round-trip = 80 CA). Review pending.
Task 4: complete (commits b56822b..7ef3f4e, review clean — spec ✅, quality approved)
  Reviewer verified rounding contract empirically (+0.0005 same hash, +0.005 differs),
  sidechain invariance, PDB column layout, and the seq-variant hash collapse.
  Foldseek accepted the PDB: createdb 1 structure / 0 skipped, 80 CA round-tripped.
  Minor (deferred): write_backbone_pdb's `name` param is dead — no live consumer
  anywhere (the tar member name is built from design_id by the caller). Signature is
  spec-mandated verbatim, so left as-is for the final review to triage.
Controller fix: removed a dead `rank = np.ones(...)` line from Task 7's planned knn.py
  (immediately overwritten on the next line) before dispatch.
Task 5: implemented (cf7391d) + fix round (891da3d). Review pending.
  CONTROLLER ERROR: my `git add -A` for a plan fix swept the implementer's in-progress
  fingerprint.py/tests into commit 6c88ebe ("Remove dead line from Task 7..."). Split into
  33c10cf (plan) + cf7391d (code). All later dispatches now mandate path-scoped `git add`.
  CONTROLLER ERROR: the out-of-frame test I added at Task 3-time was STILL vacuous —
  bd1_design's flanking nt sit 4.503 A away, 0.003 A outside the 4.5 A cutoff, so offframe
  is 0 and the assertions were identities true at zero. Rewritten (6acebad) to assert the
  real geometry, then move ONLY out-of-frame atoms into contact and require the 48-dim
  vector to stay byte-identical. Verified: offframe > 0, vector unchanged, impl untouched.
Task 5: complete (commits 33c10cf..891da3d, review clean — spec ✅, quality approved)
  Reviewer verified pair-counting (not atom-presence) by construction, H exclusion on both
  sides, and moiety_of against all 28 atom names in the fixtures (no unintended fallthrough).
  KEY: reviewer located the upstream reference LOCALLY at
  /Users/leo/scratch_scoring/src/rna_binder/analysis/structural.py:213-239 and confirmed our
  cutoff / H-filter / pair-sum semantics match it exactly. That is the property Task 23's
  r>0.9 cross-check depends on, now verified at the source rather than assumed.
  Minor (deferred): 6acebad's message describes a test edit that actually lands in 891da3d
  (blame/bisect noise); test_fingerprint.py's `baseline.offframe == 0` encodes the fixture's
  4.503 A boundary geometry as a precondition — fails loudly not silently if refixtured.
Task 6: fix applied (f9747a0) — CA-only approach. New test failed pre-fix with 0.385 A
  centroid discrepancy, confirming the sidechain artifact was real. Review pending.
Task 6b: complete (dae421a) — NEW task added mid-run. Exact parity with upstream
  rna_contact_counts on all 5 fixtures + constants pinned. 48 tests passing.
  Rationale: Task 23's r>0.9 check against cluster-stored values is late and coarse;
  upstream is importable locally, so exact equality is assertable at edit time.
Controller pre-check: umap-learn 0.5.12 precomputed_knn API VERIFIED working (both
  metric="precomputed" and default). Task 9's riskiest dependency de-risked before dispatch.
Task 6: CRITICAL bug found by reviewer and fixed (5f0109f). reference_frame was NOT
  rigid-motion invariant — SVD singular-vector signs are not equivariant under input
  rotation. Hidden because the test used ONE fixture. Reproduced: 3/6 fixtures failed,
  centroid components sign-flipped. Rebuilt as Gram-Schmidt on molecule-derived vectors
  (no SVD, no sign choice) -> equivariant by construction. Tests now parametrized over
  all 6 fixtures and compare axis DIRECTION (>0.999), not abs(dot) which HID the flip.
  approach() now uses backbone.ca_trace (was reimplementing CA selection inline).
  Re-review dispatched.
Task 6: re-review PASSED spec ✅. Critical fix verified under 120 random rotations
  (20/fixture x 6), worst axis dot 1.00000000, worst centroid drift 2e-6 A.
  Important (fix dispatched, eae8a68 in plan): _frame_atoms now returns positions
  alongside atoms — re-deriving them from the filtered set raised an opaque FrameError
  from rnaframe when a nucleotide lost all frame atoms.
Task 7: complete (commits 5f0109f..f422723, review clean — spec ✅, quality approved)
  Implementer found a REAL BUG in the brief's own code: _topk_per_node double-counted
  mutual top-k picks, eating a rank slot and evicting a distinct neighbour. Fixed by
  deduping directed (src,dst) before ranking; plan synced. Reviewer stress-verified over
  700 seeds; blocked path bit-identical to unblocked; no -inf/NaN leaks.
  MINOR (for final fix wave):
    - knn.py docstring says "~14 MB" but that is only the stored unit vectors; the real
      per-iteration working set at 72k is a (block, n) matrix ~1.2 GB plus a same-size
      negated copy for argpartition (~2-3 GB transient, freed each of ~18 iterations).
      Fits the 64 GB stage-4 job, but the comment misleads whoever tunes `block`.
    - _topk_per_node's `n` parameter is unused (inherited from the brief).
Task 6: COMPLETE (through 2034a0f). Critical + Important both fixed and verified.
  Final Important fix controller-verified directly (4-line change): grep confirms exactly
  one canonical_positions call, in _frame_atoms on the FULL rna — no re-derivation path
  remains. 17 geometry tests, 74 suite, 0 failed. Flagged for the final whole-branch
  review since it did not get its own full re-review agent.
  NOTE: the fix agent BLOCKED rather than weakening a failing test — correct call. The
  failing test was MY error (dropped nt 8 only; centroid_of raises only when the whole
  7-10 group is gone). Corrected + companion graceful-degradation test added.
Task 8: implemented (eb8b265). Review in flight.
Task 8: COMPLETE (eb8b265 + e7d0f48, review spec ✅ quality approved)
  Reviewer verified cross-PROCESS determinism (3 invocations), isolate labels not
  colliding with real clusters, medoid membership on a 300-node/28-cluster graph,
  deterministic tie-breaking, empty-edge and zero-weight cases. Measured Leiden at
  production scale: 22.5s for ~1M edges / 72k nodes, of which graph construction is 0.6s.
  Important FIXED (e7d0f48): igraph silently grows the graph when an edge index >= n,
  returning more labels than nodes with no error at the fault site. Guard added; the test
  genuinely failed first (n=3 grew to 6). Controller-verified directly; flagged for final
  review since it did not get its own re-review agent.
  MINOR (final fix wave): cluster_all() rebuilds the igraph 3x rather than reusing one —
  ~1.4s of ~66s, immaterial, measured not guessed.
Task 10: complete (3602f56) — 9 tests, suite 89. No bugs found; stress-tested both-zero
  fingerprint rows, mixed batches, string labels, ~5000x3000 contingency.
Task 9: review returned ❌ CRITICAL. My ring-placement fix is insufficient: UMAP's
  default spectral init raises IndexError from inside umap-learn once >=3 nodes are
  mutually disconnected, BEFORE the ring code or the isfinite guard runs. Threshold is
  3-4 components, not 1 — and my new test only covered n_lonely=1, hiding it.
  Reviewer verified init="random" avoids the crash and the ring logic then works.
  => I replaced a silent NaN with a hard crash. Fix in progress.
Task 9: COMPLETE (780007e, d9dbc75, 37ee77d). Critical FIXED — disconnected nodes are now
  removed BEFORE UMAP (spectral init crashed inside umap/spectral.py:86 at >=3 components).
  Verified: n_lonely 3/5/12 all crashed pre-fix, n_lonely=1 passed — which is exactly why
  the single-node test hid it. Also fixed: Euclidean (not Chebyshev) ring radius; explicit
  isfinite fallback (NaN is truthy, so `... or 1.0` preserved NaN).
Task 10: COMPLETE (3602f56). No re-review needed — no findings.
=== ALL 11 LIBRARY MODULES DONE. 94 tests passing. ===
Process change from here (user asked about overengineering; agreed to lighten):
  review agents ONLY on stages with real logic (13 extract, 16 graph, 19 payload, 20 UI).
  Pure transcription (sbatch wrappers, stage drivers) goes straight to verification.
Tasks 11+12: COMPLETE (0c3ccc6, 8781a76). 100 tests. Census GATE PASSED on real Euler data.
  EMPIRICAL CORRECTION (3rd revision of this fact): the 18-nt variant CACUGCCGUAUAGGCAGC
  occurs in the DESIGN layer only, and only for BD1+LiMPNN and HD+LiMPNN. All af3/esm
  refolds are the canonical 16-mer. Earlier guesses: BG+BIF (from entity_poly header, wrong),
  then BD1 alone (from one fixture, incomplete). Census over 105 structures governs.
  Also: 8,931 archives, not the ~6,780 I estimated from uuid count. Sharding is by archive
  so it self-balances, but stage-1 task sizing should assume ~140 archives/shard at 64 shards.
Task 13: COMPLETE (5e93acc). RFD3 inversion tally EXACTLY as expected (79/79 chain B;
  all other stacks 75/75 chain A; 154/154 accounted). Inode footprint 19 for both smoke
  shards. TMPDIR cleaned.
  *** MOST CONSEQUENTIAL BUG OF THE RUN, and it was mine. ***
  read_archive_members used `tar -O` + positional zip(members, blobs), assuming tar emits
  members in REQUESTED order. It emits in ARCHIVE order. 84% of archives hold >1 member,
  so this silently cross-assigned CIF content between design_ids for the MAJORITY of the
  dataset: every backbone, coord hash and epitope fingerprint attached to the wrong design.
  Nothing would have crashed. The map would have been complete, plausible, and wrong.
  Fixed by extracting BY NAME into real files. Plan synced.
  (2nd bug: my own step-5 verification command sorted by row_id before cut, breaking join's
  merge-sort precondition — it reported 1 match instead of 154. A broken check that would
  have "passed" the wrong thing.)
=== USER DECISION (post-progress-review): input stays passing.csv (71,866, confirmed NOT
    a 40k subsample — no such export exists in binder_filter/exports/). Ceremony stripped
    for tasks 14-21: controller writes stages directly, no per-task review agents, no new
    unit tests (the library they call is already tested). Verification is running on real
    data. ===

================================================================================
# NEW PLAN — Campaign-Wide Fold Map (three maps)
Plan: docs/superpowers/plans/2026-09-14-campaign-wide-fold-map.md
Spec: docs/superpowers/specs/2026-09-14-campaign-wide-fold-map-design.md
Branch: campaign-wide-fold-map (created off foldseek-binder-map)
Merge base: 4b6a10c
Everything ABOVE this line belongs to the 2026-08-21 plan and is DONE. Do not re-dispatch it.

Controller pre-flight fixes (before Task 1):
  - Authoring loop flipped: author LOCALLY, `bash sync.sh` to Euler, test there, commit
    locally. Laptop has no pytest/pyarrow; Euler venv has pytest 8.3.4 + biotite + zstd.
  - sync.sh written by controller (rsync per-directory; the include/exclude form
    overflowed rsync 3.2.7's rule buffer).
  - PYTHONPATH unified to src:scripts (Tasks 3,5,6,7 import from scripts/).
  - Hoisted a function-local `import io` in Task 5's split_tar.

## Completed tasks
Controller error (caught, no damage): I committed the fence fix (30893bc) while the Task 1
  implementer was live. Explicit-path `git add` meant nothing of the implementer's was swept,
  but it races the index lock and muddies the per-task BASE. RULE FOR THE REST OF THE RUN:
  no controller commits while an implementer is dispatched.
Controller fix (pre-Task-4): the plan had an unbalanced code fence in Task 3, which made
  `task-brief` treat every later heading as fenced text — Task 3's brief ran 1639 lines to EOF
  and Tasks 4-8 "did not exist". Closed the fence (30893bc); all 8 briefs now extract cleanly.
  Task 1's brief predates the fix but its boundaries were already balanced (verified).
Task 1: complete (commits 30893bc..09dab90, review clean — spec ✅, quality approved)
  Reviewer independently re-ran the Step 5 ground-truth check on Euler rather than trusting
  the report: rule keys 31767 / true backbones 31763 / over-count exactly 4. Also md5-verified
  the deployed Euler copy is byte-identical to the diff.
  Controller resolved the reviewer's ⚠️ (unverified full-suite claim): ran it — 97 passed,
  12 skipped, 0 failed, exit 0. Matches the implementer's report exactly.
  Minor (deferred to final review): SEQ_VARIANT_STACKS is public but has no explicit
  membership test, unlike its ONE_BACKBONE_PER_DESIGN counterpart.
Controller (between tasks, no implementer live): committed 17a85c7 — the mutual-kNN work
  (knn.py _mutual_filter, stage4_graph.py --mutual/--k, clustering/htmlbuild/app.js) was
  sitting UNCOMMITTED in the working tree. It is the code that produced
  results_sw_k30_f50_mut, the map the 12k was selected on. A checkout/stash would have
  destroyed it, and Task 6 edits both files so its commit would have absorbed the feature.
  Committed separately, unmodified, labelled as found-not-authored.
Task 2: complete (commits 17a85c7..ce7c5a2, review clean — spec ✅, quality approved)
  Reviewer traced the zstd subprocess lifetime across normal exit, mid-stream exception and
  generator abandonment (GeneratorExit through both `with` blocks into the finally); no leak
  or deadlock path. Confirmed no disk writes anywhere in the module.
  Controller resolved both ⚠️:
    (1) commit body is plain prose, no Co-Authored-By trailer.
    (2) REAL GAP, and it was the plan's fault: Step 5 only tested the NEGATIVE path (asks for
        a member that cannot exist). Positive-path extraction on production data was unproven,
        and cross-assigned tar members were the worst bug of the previous run. Verified
        directly: 3 real members from the SAME inner archive (the multi-member case where
        cross-assignment shows up) each came back matching its expected design_id, 45 KB,
        703 ATOM records. Permanent coverage lands in Task 4's smoke test over 200 real members.
  Minor (deferred to final review): outer zstd Popen has no return-code check or stderr
  capture, so a corrupt archive surfaces as an opaque tarfile.ReadError; inner loop does not
  stop once all wanted members are found; no test for abandoned-generator/corrupt-archive.
Controller hazard caught (pre-Task-4): the 2026-08-21 plan's run used this SAME sdd workspace,
  so stale task-3..task-13 REPORTS from 23 Aug were still present under the exact filenames my
  remaining tasks write. A reviewer handed one would have reviewed a different plan's work with
  no signal anything was wrong. Moved to _prev_plan_2026-08-21/. From here, verify a report's
  mtime is same-day before dispatching its reviewer.
Measured correction to the plan (user prompted: "parallelize extraction on cluster CPUs"):
  extraction ALREADY is a 35-way Slurm array, but measuring it showed my cost model was wrong.
  parse+hash+fingerprint = 0.031 cpu-s/structure (not the 0.55 extrapolated from v6) => 4.4 cpu-h
  for all 508,100. Streaming dominates. AND design members occupy only the first 8-10% of each
  outer archive (.list = tar order: rfd3 v2 8%, bg_bif v8 8%, hd v9 10%) — the old code streamed
  to exhaustion, reading 683 GB to use ~68 GB. Task 4 gains Step 0: stop once every wanted inner
  archive is consumed. cpus-per-task 8 -> 2 (7 sat idle). Stage 1: ~320 cpu-h -> ~5 cpu-h.
Task 3: COMPLETE (commits ce7c5a2..a67469c, review clean — spec ✅, quality approved).
  Tests 108 passed / 12 skipped / 0 failed.
  Reviewer verified by inspection: mask/index alignment between draw_sample and write_build
  (single early reset_index, nodes never resorted after), deterministic representative pick,
  build_shared == map1|map3 with consistent overlap arithmetic (500,000+9,760-1,553=508,207),
  <ID>-before-<U> substitution, manifest column+sort contract.
  OPEN DECISION FOR THE USER (Important, plan-mandated — do not silently skip on resume):
    tests/test_stage0_sample.py covers ONLY draw_sample. build_nodes and write_build — backbone
    dedup, representative selection, in_12k/.any() aggregation, template substitution, manifest
    contract — have NO repeatable tests, only one manual real-data run. That is what the brief
    specified, so it is a plan gap, not an implementer defect. Decide: accept, or add
    fixture-based tests for build_nodes/write_build.
  Controller must resolve on resume (reviewer ⚠️): row_id is assigned BEFORE the
    outer_path/inner_path sort, so it indexes nodes.parquet row order, NOT the on-disk line
    order of design.tsv. Confirm Task 4's stage-1 consumer treats it as a foreign key, not a
    file position. (Task 4 groups by outer_path and never uses row_id -> almost certainly fine,
    but verify before dispatch.)
  Minor (final review): ~7M-row Python loop for backbone_key; read_design_ids opens twice.
  Implementer returned DONE_WITH_CONCERNS; all three concerns were MY SPEC's errors, each
  verified by the controller and fixed in 3fe76bd:
    (1) map3_n is 9,760 not 9,759 — the spec compared a stage-0 NAME-rule count against a
        stage-2 COORD-HASH count. Traced to the exact pair: BG_LiMPNN_..._04bdeaa9_design7 and
        BG_LiMPNN_..._fadec0de_design67 are coordinate-identical under different uuids.
        Both numbers are correct at their own stage. Stage 2 must dedup 9,760 -> 9,759.
    (2) store path: filter_store_v2v9 on Euler is the STALE pre-HD-topup store (6,821,957).
        Correct is filter_store_v2v9_hdtop (6,992,515). Verified both meta.json.
    (3) the library-selection export has its header on LINE 1 (no comment preamble), unlike
        binder_filter exports; hardcoded skiprows=2 was wrong.
DISCOVERY for Task 7: csy4_binder_library_12k_v6_k30.csv already carries backbone_id,
  fold_cluster and epitope_cluster for all 12,000 designs — i.e. the map-2 labels, shipped.
  Use it as the crosswalk source and as a validation check against re-derived labels.

## RESUME HERE
Next: record Task 3's review outcome, then Tasks 4-8 in order.
Task 4 brief already includes the controller-added Step 0 (archive early-exit) and cpus 8->2.
Briefs for all 8 tasks are current in .superpowers/sdd/task-N-brief.md (regenerated after the
fence fix). Stale reports from the 2026-08-21 run are quarantined in _prev_plan_2026-08-21/.
Baseline suite on Euler: 108 passed, 12 skipped, 0 failed.
Task 4: COMPLETE (commits 3fe76bd..f309a02 — 419837b impl, 29f1fa5 fix round 1, f309a02 fix
  round 2). Review clean after two fix rounds. Suite 125 passed / 12 skipped / 0 failed.
  Round 1 reviewer (opus) found 4 Important, each traced experimentally not asserted:
    - by_member keyed on (inner,member) without outer_path: two rows sharing that pair across
      archives collapse, one design vanishes, the other carries a structure that is not its
      own, AND THE ACCOUNTING ASSERT STILL PASSES. Same class as the historical worst bug.
    - corrupt inner archive raised while advancing the generator (outside the per-structure
      try) -> killed the shard, leaked tmp, contradicted the code's own comment.
    - outputs written BEFORE validation (plan-mandated; controller authorised reversing it).
    - loop guard tested `wanted` not `remaining` -> duplicate inner name re-yielded.
  Fix round 1 also found a PRE-EXISTING latent bug: the archive_missing branch double-counted
    in reconciliation; any real shard with one absent outer archive would have crashed.
  Round 2 reviewer found the blocker: the guard fix had NO test — both poison tests put the
    poison LAST, so the early return fired first and they passed with the bug reverted.
    Fixed by ordering [aaaa, poison-aaaa, bbbb]; discrimination proven both ways (fails with
    `wanted`, passes with `remaining`). Controller independently verified the ordering and
    all six round-2 fixes in the committed code, plus a green suite (exit 0, 0 failures).
  *** LESSON, now three occurrences in this codebase: a test written to cover a fix can pass
      against the unfixed code. Always run the new test against the OLD behaviour. ***
  Minor (final review): tarball is the completion marker downstream globs key on; `accounted`
  covers absent+errored archives so the name is broader than "delivered".
Task 5: COMPLETE (commits f309a02..e46a38b — f3f1f72 impl, e46a38b census fixes). Review
  approved; fix round landed. Suite green (exit 0, 0 failures).
  Reviewer confirmed the tar split is genuinely streaming (BytesIO -> dest tar, no loose files)
  and that the output basenames preserve the glob contract stages 2-5 rely on.
  Fix round dispatched for three census gaps:
    (1) REAL DEFECT: the census shortfall check was gated on `name == "lib12k"`, a literal
        match against an arbitrary CLI label. A different spelling silently skips the one
        guarantee protecting against a silent shortfall. Moving to an explicit census flag.
    (2) the check counted only `hashes`, so a design correctly recorded in `missing/` also
        tripped it -> noisy, therefore ignorable. Now distinguishes explained from unexplained.
    (3) manifest_ids and the census path had ZERO tests; adding them plus a discriminate-check.
  Minor (final review): split_tar re-mkdirs a dir main() already created; the implementer's
  report misstated its own file line counts (100/49 actual vs 177/72 claimed).

⚠️ CARRY TO TASK 8 — cannot be verified until Task 4's real extraction runs:
  `split_tar` matches a PDB member stem (`<design_id>.pdb` minus `.pdb`) against design_ids
  read from the manifests. Task 3 writes the manifests and Task 4 writes the PDB names; if
  those two ever disagree on id format by even one character, EVERY shard produces a
  well-formed, EMPTY output tar with no error at all. Task 8 must spot-check a real shard's
  member names against a manifest row BEFORE the 500k search is paid for. Add it beside the
  existing coord-hash identity gate.
  Fix round outcome (agent was killed by an API error AFTER committing, so it never reported;
    controller verified the work directly): census trigger is now the explicit
    should_warn_census(is_census, unexplained) driven by a --map ...:census flag with malformed
    specs rejected, not a string match on a CLI label. reconcile() separates explained from
    unexplained shortfalls. 15 tests in test_split.py now, including 2 for manifest_ids,
    4 for reconcile, 4 for should_warn_census.
    DISCRIMINATION IS STRUCTURAL HERE, not empirical: should_warn_census takes two booleans-ish
    inputs and all four combinations are tested, so breaking it in either direction fails a
    test. That is a stronger guarantee than the temporary-revert experiment.
Task 6: COMPLETE (commits e46a38b..3b005ae — 26786e7 impl, 3b005ae threshold fix).
  Review approved. Suite 147 passed / 12 skipped / 0 failed (159 collected), exit 0.
  TWO AGENTS STALLED on infrastructure faults; the second left load_hits + its tests but had
  not rewired main() or touched knn.py. Controller finished those mechanical edits directly
  and told the reviewer so, since controller-reviewing-own-work is the case review exists for.
  Reviewer independently re-ran the suite on Euler and grepped the LIVE tree (not just diff
  context) to confirm _mutual_filter/_topk_per_node/_assemble/symmetrise_min and the
  --mutual/--k/--tm-floor/--epi-floor surface are untouched. It also verified the knn.py
  docstring's memory arithmetic (1024 x 508k x 4B = 2.08 GB) rather than taking it on trust.
  *** CONTROLLER ERROR, twice over, on the reviewer's Important finding: ***
    (a) I asserted a "145 baseline" I had never measured. Real baseline was 140; 140+6+1=147
        reconciles exactly. Do not state test counts without running them.
    (b) I "confirmed" the float32/float64 threshold flip using the LAPTOP's numpy, which
        predates NEP 50. Euler runs numpy 2.1.3: a float32 ARRAY compared to a bare Python
        float adopts the array dtype, so production already agreed with the old float64 path.
        The reviewer's finding was a false positive on this version — and my first test for it
        was VACUOUS, passing either way. VERIFY IN THE PRODUCTION ENVIRONMENT, not locally.
  Kept the fix anyway as REDUNDANCY_TM = np.float32(0.95): harmless, self-documenting, and it
  makes the behaviour independent of promotion rules. The rewritten test pins agreement with
  the pre-streaming float64 reference and FAILS against a float64 threshold (verified).
  Minor (final review): stage4_graph.py's load_hits docstring says ~12 bytes/hit and ~1.8 GB;
  four 4-byte arrays is 16 bytes/hit and ~2.4 GB at 150M rows (inherited from the brief).
Task 7: COMPLETE (commits 3b005ae..<head> — e1a4c7f impl + folds-column test).
  Review approved. Suite green, exit 0, 0 failures.
  Reviewer independently re-ran the suite on Euler (158/12/0, 170 collected) rather than
  trusting the implementer's dot-counted derivation, and TRACED struct_name back to
  stage2_builddb.py:44 to rule out the "folds silently equals size" failure mode — struct_name
  is the lexicographically-first design_id of a coord-hash group, so it is genuinely
  design-id-shaped text and fold_key's suffix stripping operates on the right thing.
  It also cross-checked fold_key's premises against universe.py's measured split-rate table
  instead of accepting the comment.
  Important finding CLOSED by the controller: the mandatory `folds` column had no test that
  could catch it regressing — it was built inline in an untested main(). Extracted
  cluster_table(labels) as a pure function and added a discriminating test (3 seq variants of
  one backbone = 1 fold / 3 nodes); it FAILS when the aggregation is switched to struct_name.
  Minor (final review): load_labels and design_to_cluster each re-read the same three files
  per map, doubling I/O on the real trees (inherited from the brief).
Task 8: COMPLETE (commit 984ca2c, base 2db16b9 — review gate run 2026-09-24, nine days late;
  the task was committed 09-15 but never review-gated, so the ledger stopped at Task 7).
  Spec ✅ (reviewer verified every frozen parameter END TO END, not by grep: ALN_TYPE=1 reaches
  stage3_search.sbatch:15,39; -c 0.8/--cov-mode 0/-e 10/--max-seqs 300 hardcoded at :38;
  --tm-floor/--epi-floor/--k/--mutual match stage4_graph.py:104-111 and override its non-frozen
  0.5/15 defaults; Leiden SEED is module-level in clustering.py and cluster_r80 hardcoded in
  stage9_compare.py:29, so both are correctly NOT passed; map 2 is read-only at
  stage9_compare.py:34; no b2 anywhere). The implementer's one deviation from the brief
  (`:census` on the lib12k --map) was independently confirmed correct against
  stage1b_split.py:70-97.
  *** IMPORTANT finding, plan-mandated, FIXED BY THE CONTROLLER (uncommitted — see below) ***
    submit_maps.sh:57-58 `J5_LIB=$(submit_map ...)` / `J5_FULL=$(submit_map ...)` had NO
    `|| exit 1`. submit_map runs inside a command substitution, so the `|| exit 1` guards on the
    eight sbatch calls INSIDE it kill only that subshell — the parent script sails on. The
    implementer's report asserts the opposite ("set -u plus || exit 1 after every sbatch call
    means the script aborts immediately"); that claim is FALSE for all eight inner calls.
    Failure mode: a failed lib12k submission does not stop the script, the ~2,100 CPU-h full
    chain is submitted anyway, and stage 9 goes out as `--dependency=afterok:$J5_FULL:` with an
    empty field and no message naming which chain died.
    Verified both ways in a sandbox repro (scratchpad/guard_repro.sh): unguarded reaches the
    stage-9 line with J5='' and exit 0; guarded aborts with exit 1. Fix = append `|| exit 1` to
    both assignments. `bash -n` clean after. The defect is verbatim in task-8-brief.md's Step 2
    block, so it predates the implementation — provenance does not downgrade it.
  Minor (carry to final review):
    - probe_scaling.sh's `for N in 50000 150000` loop has no failure check on its sbatch calls
      either. Low stakes (two cheap diagnostic jobs), byte-identical to the brief.
    - probe_scaling.sh:7 defaults RESULTS=results_full, whose Foldseek DB does not exist until
      full_s2 has run — i.e. the probe CANNOT precede the launch. results_v6 is not a substitute:
      25,359 nodes cannot feed the N=150000 arm. Inherent to the design; document, don't "fix".
    - Brief Step 3 says "seven job ids printed"; the script makes 11 sbatch calls (1+1+4+4+1).
  ⚠️ DOC BUG IN THE BRIEF, controller-found, not the reviewer's: Step 3 expects
    `map3_n` in build_shared/meta/sample.json to be exactly 9,759. The brief's own Notes section
    and the spec's validation §5 both say 9,760 (stage 0 name-rule rows), which stage 2 then
    dedups to 9,759 coord-hash nodes. The Notes call this "an error the spec made once".
    Following Step 3 literally raises a false alarm on a CORRECT run. Use 9,760 at stage 0.
  ⚠️ UNCOMMITTED: the `|| exit 1` fix is in the working tree only. The session that made it was
    pinned to an unrelated git worktree and the harness refuses git against this checkout, so it
    could not commit, could not run `git diff`, and could not run the final whole-branch review.
    The Task 8 review package was assembled by READING the two files, not from git diff — safe
    because 984ca2c adds only those two files, but unproven. NEXT SESSION, from
    ~/Documents/binder_map: commit submit_maps.sh, then run the final whole-branch review
    (merge-base..HEAD) with this Minor list in hand.
