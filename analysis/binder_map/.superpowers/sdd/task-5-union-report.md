# Task 5 — union fold map: colour-by-stack + per-cluster medoid viewer

STATUS: IN PROGRESS (report written early and updated as work lands, per brief)

## What is being added
1. Colour mode: colour points by design stack (payload already carried
   `stack_codes` + `stacks` lookup, unwired).
2. Per-cluster medoid backbone viewer, reusing `web/app.js`'s
   `pdbFromTrace` / `render3D` 3Dmol machinery.

## Findings so far
- `results_union/design/clusters/fold/medoids.parquet`: 160,668 rows across
  5 resolutions (r05 31,886 / r10 32,001 / r20 32,125 / r40 32,262 /
  r80 32,394). **Distinct `struct_name` across all resolutions = 32,707**,
  vs 32,394 for r80 alone — restricting to r80 would save only 313 traces
  (~1%), so the r80-only fallback in the brief buys nothing.
- `results_union/design/shards/` exists (99 shards, 2.5 GB) so stage 7 runs
  against `results_union` unmodified (`--layer design --results .../results_union`).
- Only the fold axis has clusters/medoids, as expected.

(to be filled in)
