# GUARDRAILS — the campaign's institutional memory, executable

Each row is a scar-lesson from the csy4 campaign turned into a `campaign doctor` check that
**refuses to launch** when violated. `doctor.py` is the executable mirror of this table.

| # | Check | Rule enforced | Why (the lesson) |
|---|-------|---------------|------------------|
| 1 | node-local scratch | archive extraction resolves under `$TMPDIR`/`/scratch`/`/tmp`, never `/cluster` | 2026-07-10: atexit-only cleanup leaked 1.35M files onto Lustre → inode brake stalled the whole campaign |
| 2 | mem form | header uses `--mem-per-cpu`, never `--mem` | Euler rejects `--mem` for these partitions |
| 3 | preempt-safe | `--requeue` **and** `--open-mode=append` present | preempt/requeue must resume without clobbering logs; TIMEOUT/scancel are terminal |
| 4 | pack + verify | outputs written as `tar.zst`, `zstd -t` before deleting originals | inode footprint + integrity (verified-before-delete) |
| 5 | per-unit seed | models with `per_unit_seed: true` get a distinct seed per uuid | HalluDesign emits duplicate designs otherwise |
| 6 | omit residues | `constraints.omit_residues` propagated to the inverse-fold flags | cysteine leaked into RFD3-LiMPNN + BG-BIF before configs were fixed |
| 7 | refold routing | a model is submitted only to pools in its `runs_on`; AF3 never on a `<24h` bulk tier; refolds use the 6h bulk walltime | AF3 TIMEOUTs on short bulk; AF3/ESM off scarce pools per compat |
| 8 | blackhole exclude | `--exclude` includes `cluster/euler.yaml: blackhole_excludes` | a wedged-GPU node black-holes (fast-fails) every job routed to it |
| 9 | pool-count truth | HD success is judged by packed-design count, not sacct `State` | HD gen exits FAILED (benign DeepSpeed/Triton atexit) yet packed all 150 |
| 10 | admission control | `hold_governor` is running before a large gen fan-out | Euler evaluates only the top ~250 pending; an aged own-backlog floods that window and starves refolds |
| 11 | requeue semantics | array/bulk jobs set `JobRequeue`-safe flags; refold TIMEOUT treated as terminal, not retried in place | preemption on bulk; NODE_FAIL floods wipe jobs terminally |

Add a row whenever a new failure mode is learned; keep the "why" so the next person understands
the cost, not just the rule.
