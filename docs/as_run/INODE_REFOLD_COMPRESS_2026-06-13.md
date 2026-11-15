# csy4 inode-quota crisis + refold compression (2026-06-13/14)

Single-session as-built. The campaign's jobs were failing with `Disk quota exceeded`;
root-caused to an **inode (file-count) quota**, and reclaimed ~890k inodes by
compressing the loose refold leaves. Complements the concurrent `_backbones`
work (see "Related").

---

## 1. The blocker — INODE quota, not space

Jobs (boltzgen/bd/rfd3 `*_multif`, AF3/B2 refolds) died with
`OSError: [Errno 122] Disk quota exceeded` — including 9/15 of the 2026-06-12 bd
resubmits and ~130 boltzgen refolds.

- `lquota /cluster/work/wicky` → **space 9.5 / 30 TB (fine)**, **files 7.97M→8.12M
  / 7.5M soft / 8.25M hard** — over soft, ~130k under hard.
- **The binding cap is actually the PROJECT quota** (`lfs quota -p /cluster/work/wicky`,
  id 510551) at 8.25M, which `lquota` hides — prefer `lfs quota -p` for true headroom.
- It's the **group/project** quota for the whole Wicky lab on `/cluster/work/wicky`,
  not just lvalsecchi.

## 2. Where the inodes live (csy4_crrna_free ≈ 3.25M of the 8M)

`lfs find` per stack:

| area | inodes |
|---|--:|
| rfd3__ligandmpnn | 1,450,191 |
| boltzdesign__ligandmpnn | 743,786 |
| boltzgen__ligandmpnn | 400,600 |
| boltzgen__boltzif | 292,834 |
| halludesign__ligandmpnn | 254,851 |
| _archive (pre_routeb, ONLY copy — keep) | 49,383 |
| _backbones (already cleaned by the other session) | 18,996 |
| **campaign total** | **3,254,934** |

macOS junk negligible (396 `._*`). The hog is **incomplete refold leaves**: within
a stack, ~92% of inodes are in `midrun_v2/refold_{af3,b2}` (rfd3 midrun_v2 alone:
af3 645k + b2 483k).

## 3. Why refolds are incomplete (the inode source)

A refold leaf only packs when **complete** (gate: `expected.json` `n_source_designs`
vs `outputs/` count). Interrupted leaves keep a full `inputs/` (600) + partial
`outputs/` ≈ ~1700 inodes each and can never pass the gate → stuck loose.

Terminal-failure causes from the leaf `slurm` logs (leaves with ≥1 occurrence;
`--open-mode=append` so they overlap, not terminal-cause):

| cause | refold_af3 (1711) | refold_b2 (857) |
|---|--:|--:|
| **PREEMPTION** (bulk tier) | 759 | 336 |
| **DISK QUOTA** (Errno 122) | 494 | 156 |
| **TIME LIMIT** (timeout) | **160** | **0** |
| node failure / GPU-MPS / OOM | 3 / 7 / 0 | 3 / 0 / 0 |
| no failure sig (completed/in-flight) | 631 | 438 |

- **AF3 leaves fold ~600 designs ×~49 s ≈ 8 h** (needs gpuhe.24h) → exposed to
  preemption + quota + timeout → ~half end loose (rfd3 midrun_v2: 315/600 loose).
- **B2 is ~25 min** → never times out (0), mostly completes & packs.
- Vicious loop: loose leaves inflate inodes → quota fills → more refolds die
  incomplete → more loose leaves.

## 4. The fix — compress loose refold leaves by uuid, all stacks

`~/rna_binder/refold_compress_all.slurm` (job 3261559, normal.24h, --requeue;
log `_refold_compress_2026-06-13.log`). For every `refold_*` leaf it runs
`pack_unit "$leaf/outputs"` and `pack_unit "$leaf/inputs"` (the tested, flock'd,
crash-safe packer: tar→`.tar.zst` + `.list` sidecar, removes dir only on success,
idempotent). Guards: **skip a leaf if its uuid has an active squeue job OR its
slurm/outputs mtime < 30 min** — never tar an in-flight refold.

Result: **COMPLETED 2026-06-14 12:33 (17h13m): 4,285 leaves packed, 0 FAIL, 91 skipped.**
Project inodes **8,119,741 → 5,506,177** (`lfs quota -p 510551`) = **~2.6M reclaimed** —
compress dominant (OK 1,302→4,285 alone dropped ~1.64M with nothing else running) plus
the concurrent `_backbones` cleanup. Now **~2.0M under the 7.5M soft / ~2.74M under the
8.25M hard** — crisis resolved. Idempotent, so re-runnable whenever loose leaves re-accumulate.

- **GOTCHA (cost me 2 restarts):** the enumerating `find` MUST be `-maxdepth 4`
  (refold dirs at depth 3–4). A depthless `find "$CAMP" -type d -name refold_*`
  walks the whole 3.25M-inode tree before emitting anything (ran 6 min, packed 0).
  A `| sort` on it also blocks streaming — removed.

## 5. What got compressed — completeness

Final census over all 2,116 packed leaves — **0 are verifiably complete**, and that's
correct: complete leaves were already packed+cleaned by the pipeline, so they had no
loose `outputs/` for this job to touch.

| bucket | leaves | files (out+in) |
|---|--:|--:|
| COMPLETE (no continuation needed) | 0 | 0 |
| INCOMPLETE (needs unpack + re-fold) | 856 | 1,758,304 |
| no `expected.json` (old/analyzed, ~complete) | 1,260 | 521,081 |

So **to fully continue the interrupted refolds = 856 leaves / ~1.76M files** (~1.23M
partial output predictions + ~0.53M inputs staging). The 1,260 no-`expected.json` leaves
are old-format runs (600_seq / 3kseq / analyzed 33ac0cbd·f2653389), not active partials —
completeness unverifiable but mostly already-complete analyzed data.

Sample archive validated (`zstd -t` OK, 64 MB, 62 output dirs / 1267 files).

## 6. Recovery (how to restore a compressed leaf)

`pack_unit`'s inverse, in `code/scripts/stages/archive_pipeline.sh`:
```bash
source code/scripts/stages/archive_pipeline.sh
export PATH=.../rna_binder/code/bin:$PATH
unpack_unit <leaf>/outputs.tar.zst     # recreates <leaf>/outputs/  (flock'd, idempotent)
unpack_unit <leaf>/inputs.tar.zst
```
Inspect without unpacking: `cat <leaf>/outputs.tar.zst.list`, `count_unit`,
`extract_member <archive> <path>`, `ensure_unit`. The aggregator / GUI / opendesign
are archive-aware (read `.tar.zst` directly). Manual fallback:
`zstd -dc outputs.tar.zst | tar -x -C <leaf>`. **Compression is non-destructive.**

## 7. Follow-ups

- **Completing the 397 partials needs more than a requeue:** each carries a
  `designs.tar.zst` finalize marker → a plain requeue exits early ("already
  complete"). To re-fold: clear the marker + resubmit on a **non-preempt 24 h** slot
  (preemption on bulk is the #1 reason they never finished).
- Consider moving AF3 refolds off the preemptible bulk tier so 8 h jobs can finish.
- Keep an eye on `lfs quota -p` headroom; the compress should be re-runnable (idempotent)
  whenever loose leaves re-accumulate.

## Related
- `_backbones` workdir-scratch leak + source fix (concurrent session): the Tier-B
  `purge "$WORK_*"` enable in the 3 gen multif + 1028-backbone cleanup.
- `boltzgen_patch/ROUTEB_AND_RECOVERY_2026-06-11.md` (finalize-tail / Route B / 6kseq recovery).
- midrun_v2 loose-design recovery (2026-06-12): rfd3 job 3089031, 15 bd resubmits, 11 re-routed.

---

# Resubmission plan — the 856 incomplete refold leaves (2026-06-14)

Target list: `~/rna_binder/_refold_resubmit_targets_2026-06-14.tsv`
(cols: stack, label, rtype, leaf, uuid, ndone, n_source, pct). Wrapper:
`~/rna_binder/resubmit_refolds.sh` (DRY-RUN by default; `--go` to submit).

## Scope (856 leaves = 629 AF3 + 227 B2; 842 are midrun_v2, 14 are 6kseq_4h)

| stack | AF3 | B2 |
|---|--:|--:|
| rfd3__ligandmpnn | 282 | 108 |
| boltzdesign__ligandmpnn | 186 | 67 |
| boltzgen__ligandmpnn | 77 | 30 |
| boltzgen__boltzif | 75 | 22 |
| halludesign__ligandmpnn | 9 | — |

**Reality check — this is a near-full re-fold, not "finishing":** completeness
histogram = **613/856 (72%) <10 % done, 186 at 10–49 %, only 57 ≥50 %.** Most leaves
barely started, so there is no partial work to "continue" — re-folding restarts each
from source. Cost is dominated by the **629 AF3 leaves × ~8 h ≈ thousands of GPU-h.**
**DECIDE SCOPE FIRST** (e.g. B2-only? just the 57 ≥50 %? specific high-value uuids?)
before committing the full AF3 spend.

## Prereqs — the 3 root-cause fixes (else it fails identically)
1. **Non-preempt 24 h for AF3.** Preemption (bulk tier) was the #1 cause (759 AF3 leaves)
   → AF3 on **`gpuhe.24h` + 23:30**, B2 on `gpuhe.4h`. Never bulk.
2. **Clear the finalize marker.** CONFIRMED: both `refold_af3.sh`/`refold_b2.sh` early-exit
   if `$REFOLD_DIR/designs.tar.zst` exists. They read SOURCE designs from the intact
   `$SOURCE_RUN/designs/<uuid>.tar.zst`, so the wrapper deletes the leaf's
   `designs.tar.zst` + my packed `outputs/inputs` and the refold regenerates cleanly.
3. **Throttle inodes.** Re-folding regenerates loose `inputs/`+`outputs/`. Headroom is
   ~2.0 M under soft / ~2.7 M under hard. Run in waves; completed leaves auto-pack (free
   inodes), so steady-state in-flight ≈ concurrency × ~3k. **Watch `lfs quota -p 510551`;
   pause if it climbs past ~7.8 M.**

## Mechanism — `resubmit_refolds.sh <af3|b2|all> [maxN] [--go]`
Per target: skip if uuid has a live job → clear `designs.tar.zst`/`outputs.tar.zst`/
`inputs.tar.zst` (+sidecars) and stale dirs → `sbatch --partition=<part> --time=<t>
<leaf>/run.sbatch`. DRY-RUN prints intended actions; `--go` submits; `maxN` caps a wave.

## Recommended sequence
1. **B2 wave first** — `resubmit_refolds.sh b2 0 --go` (227 leaves, ~25 min each, cheap;
   validates the clear+resubmit path end-to-end).
2. **AF3 pilot** — `resubmit_refolds.sh af3 20 --go`; confirm they run on gpuhe.24h,
   COMPLETE, and auto-pack (leaf gets `designs.tar.zst` back + outputs packed); watch inodes.
3. **AF3 in waves** of ~100–150 only if the pilot is clean and the scope decision says go.
4. Verify per leaf: completeness gate passes (`outputs ≥ n_source_designs`) and the leaf
   re-finalized. Re-run the compress (`refold_compress_all.slurm`, idempotent) afterward
   to repack anything left loose.

---

# Resubmission — execution notes & scheduling (2026-06-14/16)

## Wrapper hardening (`~/rna_binder/resubmit_refolds.sh`) — bugs found during validation
Re-folding a cleared leaf needed more than the first draft. Final wrapper does, per leaf:
- **strip `#SBATCH --dependency=afterok:<gen-job>`** (gen jobs gone → sbatch "Job dependency problem").
- **remove stale `.cache*` dirs** (`.numba_cache`, `.triton_cache`, …) on clean — leftover caches
  caused boltz/lightning `RuntimeError: Missing folder: …/lightning_logs/version_0` (B2 pilot failed 5/5 until added).
- **`--mem-per-cpu=3G` for B2** (=12G vs 32G default; measured peak MaxRSS ~7.5G). Euler rejects
  `--mem=` (by-node) → must use `--mem-per-cpu`. The 32G over-request was the `QOSMaxMemoryPerUser` blocker.
- **active-leaf guard via `scontrol show job … StdOut`** — NOT `squeue %o` (returns `(null)` for
  stdin-submitted jobs); plus **skip leaves already finalized (`designs.tar.zst`)**. Without this it
  re-cleared the 30 running leaves mid-run.
- partition override: B2 → gpuhe.4h, AF3 → gpuhe.24h/23:30 (off `bulk`).

## Refold config mis-tiering (root of the bulk/timeout mess)
`code/configs/refold__{af3,b2}/default.yaml`: **af3 = gpuhe.4h / 04:00:00 (TOO SHORT — af3 needs ~8h
→ the 160 timeouts)**, b2 = gpuhe.24h / 12h. The actual `--partition=bulk` on the run.sbatch came from
a **launch-time override** (GUI preset / `--refold_*_override`), not these defaults. To run the live
pipeline off bulk: fix these defaults (af3→gpuhe.24h, b2→gpuhe.4h) and/or drop the bulk override.

## Scheduling — progress-priority on bulk (`~/rna_binder/refold_priority_babysitter.sh`)
Keep bulk's throughput without the partial-refold bloat by prioritizing near-done refolds. Euler levers:
`scontrol top` DISABLED; users can't RAISE priority; priority ≈ all fairshare (weight 10000 vs age 10,
so age won't auto-help); fairshare ~constant within your own jobs → **`nice` is the intra-user ordering
knob** (verified you can set it). Babysitter nices low-progress refolds DOWN (nice ∝ designs remaining)
so near-done jump ahead; refolds resume per-design, so a prioritized near-done job finishes fast and
drains the backlog. Loop it every 5–10 min in the live pipeline (dry-run default, `--go` to apply).

## Status (2026-06-16)
- **B2:** 30-leaf concurrency batch ran clean (0 fail, cache fix held) → remaining ~197 released.
- **AF3:** 1 pilot launched on gpuhe.24h to measure MaxRSS (→ set lean `AF3_MEM` for QOS) and validate
  the AF3 path before any wave. 629 AF3 = near-full re-fold (72% <10% done) — scope deliberately.
