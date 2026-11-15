# tar.zst archiving — integration sketch

Wrapper: `archive_pipeline.sh` (sources the lab `utils/files/archive.sh`).
Deploy to: `code/scripts/stages/archive_pipeline.sh` on Euler.
**Hold all stage-script edits below until the queue drains** — worker scripts are
read at job run-time, so editing them now changes behaviour for the pending
BD1/RFD3 refolds mid-campaign.

At rest a run dir becomes:

```
<label>_<uuid8>/
  run_meta_<uuid8>.json        # plain (tiny, read constantly) — leave as-is
  resolved_config.yaml         # plain
  slurm_<uuid8>.log / .err     # plain
  designs.tar.zst   + .list    # was designs/        (1 inode + sidecar)
  refold_af3.tar.zst + .list   # was refold_af3/
  refold_b2.tar.zst  + .list   # was refold_b2/
```

Dependency chain launch.py builds (each `pack` is a tiny CPU-only `afterok` job):

```
gen ──▶ pack(designs) ──▶ refold_af3 ──▶ pack(refold_af3) ──┐
                      └──▶ refold_b2  ──▶ pack(refold_b2)  ──┴──▶ aggregate
```

Every compute stage: **unpack inputs to node-local `$TMPDIR`**, compute, **purge
disposables**, and its paired pack job **packs only its own output**. Read-only
inputs are never re-expanded in the shared run dir — no lock contention, no
shared-fs inode churn.

---

## Phase 1 — delete disposables (safe, independent, ship first)

Biggest inode win, no read-path change. Add after the canonicalize call.

**`rfd3_multif.sh` / `boltzgen_multif.sh` / `boltzdesign_multif.sh`** (mpnn working
dirs live at the run root, siblings of designs/ — canonicalize never touches them):
```bash
source "$(dirname "$0")/archive_pipeline.sh"
# ... existing canonicalize_designs.py call ...
purge "$MPNN_RUN_DIR"/mpnn_native_design* "$MPNN_RUN_DIR/_work"
```

**`postprocess_halludesign.sh`** (canonicalize sweeps the nested `recycle_*/` trees
into `_work/`; just delete `_work/`):
```bash
source "$(dirname "$0")/archive_pipeline.sh"
# ... existing canonicalize_designs.py call ...
purge "$RUN_DIR/_work"
```

## Phase 2 — pack/unpack

### Pack jobs (emitted by launch.py)

A pack job is a ~5-line CPU sbatch, no GPU, `afterok` on its compute job:
```bash
#SBATCH --partition=<cpu>           # NOT gpuhe — packing is IO/CPU
#SBATCH --dependency=afterok:<compute_jobid>
#SBATCH --output=<unit_parent>/pack_<unit>.log
set -euo pipefail
source /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/archive_pipeline.sh
pack_unit "<RUN_DIR>/<unit>"        # designs | refold_af3 | refold_b2
```

In `launch.py`, alongside each builder (`build_rfd3_multif_sbatch`,
`build_refold_af3_sbatch`, …), add `build_pack_sbatch(unit_dir, dep_jobid)` and
submit it right after the compute job, capturing nothing:
```python
gen_jobid = submit(gen_sbatch)
pack_designs_jobid = submit(build_pack_sbatch(run_dir/"designs", dep=gen_jobid))
# refolds now depend on the PACK job (so designs.tar.zst exists), not gen:
af3_jobid = submit(build_refold_af3_sbatch(...), dep=pack_designs_jobid)
submit(build_pack_sbatch(source_run/"refold_af3", dep=af3_jobid))
# ...same for b2...
```
For `--num-runs>1` (shared `designs/`): submit ONE `pack(designs)` with
`afterok:<all N gen jobids>` (mirrors how refolds already wait on all N).

### Refold workers read designs from scratch

Top of `refold_af3.sh` / `refold_b2.sh`:
```bash
source "$(dirname "$0")/archive_pipeline.sh"
DESIGNS_IN="$SOURCE_RUN/designs"
if [[ ! -d "$DESIGNS_IN" ]]; then            # at rest it's designs.tar.zst
    unpack_unit "$SOURCE_RUN/designs.tar.zst" "$TMPDIR"
    DESIGNS_IN="$TMPDIR/designs"
fi
# ... use $DESIGNS_IN everywhere the script currently uses $SOURCE_RUN/designs ...
```
Then after refold completes, purge its scratch inputs (the `$TMPDIR` copy is
auto-reaped by SLURM; only purge anything written back into the run dir):
```bash
purge "$REFOLD_DIR"/*/inputs
```

## Phase 3 — read straight from archives (drops the last expansions)

**`aggregate_metrics.py`** — per source run, expand the needed units to scratch
mirroring the run layout, aggregate, discard. Wrap its source-run loop:
```bash
# in run_structural.sbatch, before the python call, for each SRC run:
for u in designs refold_af3 refold_b2; do
    unpack_unit "$SRC/$u.tar.zst" "$SCRATCH/$(basename "$SRC")" 2>/dev/null || true
done
# point aggregate_metrics.py at $SCRATCH/* instead of campaigns/*
```
(tar is sequential → one streaming decompress per unit, which is the *ideal*
access pattern for the aggregator since it reads every member once.)

**`opendesign.py`** (+ make_html_plots click-to-open) — extract the one design +
its refold variants to a `/tmp` cache, open those:
```python
extract_member(f"{run}/designs.tar.zst", member, "/tmp/odcache/")
# build the ChimeraX `open` paths against /tmp/odcache instead of the mount
```

**`inspect_runs.py`** — count/list via the `.list` sidecar (`count_unit`), never
decompress; big perf win over the current dir crawl on the sshfs mount.

---

## Safety / invariants

- `pack_unit` deletes the source dir **only** after `tarz` verifies (`zstd -t`)
  and the archive is atomically renamed into place. A kill mid-pack leaves the
  source dir intact.
- All pack/unpack run under `flock <unit>.lock` — required for shared `designs/`
  (`--num-runs>1`) and parallel refolds.
- Idempotent: `pack_unit` no-ops if already packed; `unpack_unit`/`ensure_unit`
  no-op if already expanded.
- `purge` refuses empty/`/`/`$HOME` and anything not under a `campaigns/` path;
  caller pre-expands globs.
- Per-unit overhead is 2 extra tiny inodes (`.list`, `.lock`) — negligible vs.
  the thousands collapsed.
