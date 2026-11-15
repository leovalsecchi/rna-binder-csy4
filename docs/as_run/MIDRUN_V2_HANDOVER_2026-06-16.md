# midrun_v2 RNA-binder campaign — HANDOVER (2026-06-16)

Campaign dir (Euler): `/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free`
Code root: `/cluster/work/wicky/lvalsecchi/rna_binder/code` (launch.py in `code/scripts`)

## NORTH STAR (the actual goal — we drifted into cleanup, this is the target)
Get **every stack except halludesign to ≥300 completed-generation uuids** (a "completed" uuid = a packed design slice `<stack>/midrun_v2/designs/<uuid>.tar.zst`).
- Submit generation **primarily in `gpuhe.4h`**.
- **boltzdesign `design_samples` halved to 30** (so it finishes inside 4h — the default config is still 60 + 4h, which is the walltime bug; the launch works around it via `--override generator_params.design_samples=30`).
- **auto-refold → `bulk`** (`--refold-af3-override slurm.partition=bulk --refold-af3-override slurm.time=23:30:00 --refold-b2-override slurm.partition=bulk`).
- **all GPUs included** (24G turing pool etc. — already configured; AF3 needs Ampere+/xla, B2 any).

## CURRENT STATE (verified 2026-06-16)
Packed slices (completed gen) + the pending gen wave already queued (gpuhe.4h):

| stack | packed now | gen pending | projected | gap to 300 |
|---|---:|---:|---:|---|
| boltzgen__boltzif | 157 | 143 (bg) | 300 | ✅ covered (no buffer) |
| boltzgen__ligandmpnn | 156 | 143 (same bg jobs feed both) | 299 | ✅ covered (no buffer) |
| boltzdesign__ligandmpnn | 170 | 130 (bd, design_samples=30) | 300 | ✅ covered (no buffer) |
| rfd3__ligandmpnn | 263 | 5 (rfd3) | 268 | ❌ **SHORT ~32** |
| halludesign__ligandmpnn | 29 | (269 hd jobs queued) | — | excluded from target |

**One boltzgen gen job produces designs for BOTH `boltzif` and `ligandmpnn`** (Route B = BIF + LiMPNN from the same run). That's why those two stacks track together (~443 attempted, 157/156 packed).

### #1 NEXT ACTION
**rfd3 is under-queued: only 5 gen runs queued, needs ~37 to reach 300.** Launch ~40–45 more rfd3 (buffer for failures). Also note boltzif/ligandmpnn/bd are queued to *exactly* 300 with **no failure buffer** — once the wave drains, recount and top up any that fell short.

Launch pattern (from `/tmp/gen_launch.sh`, reproduce on Euler `code/scripts`):
```
python3 launch.py ../configs/rfd3__multif/default.yaml --label midrun_v2 --num-runs 45 \
  --auto-refold \
  --override slurm.partition=gpuhe.4h --override slurm.time=04:00:00 \
  --refold-af3-override slurm.partition=bulk --refold-af3-override slurm.time=23:30:00 \
  --refold-b2-override slurm.partition=bulk
```
(boltzdesign adds `--override generator_params.design_samples=30`.)

## ROOT CAUSE of the past gen failures (RESOLVED)
The ~285 boltzgen uuids that produced no designs failed with **`OSError: [Errno 122] Disk quota exceeded` at `output_dir.mkdir()`** — they ran *during* the inode-quota crisis and couldn't create their output dir. **NOT a pipeline bug.** Quota now has headroom, so the re-queued wave should just succeed. (Diagnosed from gen slurm logs: 59 quota / 64 filenotfound / 63 traceback — same cascade; only 11 genuine TIME LIMIT, 3 GPU.)

## INODE STATUS — healthy
`lquota /cluster/work/wicky`: **files 5.94M / 7.5M soft / 8.25M hard** (~1.56M headroom). Crisis resolved (earlier refold compress reclaimed ~2.6M).

## CLEANUP DONE THIS SESSION
1. **Deleted 287 stale failed-gen footprints** (143 boltzif + 144 ligandmpnn) — the Errno-122 failures. Each = `run_<uuid>/`, `run_meta_<uuid>.json`, `slurm_<uuid>.{log,err}`, 3 empty refold leaves. List: `$RB/_stale_delete_2026-06-16.list`. Guard: protected 894 active uuids + 1 boltzdesign uuid that had real refold work.
2. **Deleted 133 empty (no-design) staging dirs** (`<stack>/midrun_v2_<uuid>/` scratch). List: `$RB/_stale_staging_del_2026-06-16.list`.

## OPEN / DEFERRED ITEMS
- **80 staging dirs hold LOOSE designs → must be PACKED, not deleted** (NOT touched). Breakdown: rfd3 46, halludesign 29, boltzdesign 2, boltzgen__boltzif 2, boltzgen__ligandmpnn 1. List: `$RB/_staging_withdesigns_2026-06-16.list`. These are `<stack>/midrun_v2_<uuid>/` (rfd3 holds real .cif/.pdb; many rfd3 are NOT yet packed = possibly only copy). Use the pipeline pack (`pack_unit`/canonicalize via `archive_pipeline.sh`), not `rm`.
- **Refold resume wave:** ~1414 af3 PENDING + 27 RUNNING + 563 b2 PENDING on bulk. Babysitter `refold_prio_loop` (currently job **3576091**, self-chaining on normal.4h) nices low-progress jobs down. STOP it via: `touch $RB/_refold_prio_loop.STOP && scancel --name=refold_prio_loop`.
- **269 hd jobs queued** — halludesign is excluded from the 300-target; decide whether to keep or cancel (they consume scheduler share against the gen wave we care about). Don't cancel without confirming with user.
- **No failure buffer** on boltzif/ligandmpnn/bd (queued to exactly 300). Recount after the wave drains; top up shortfalls.

## KEY GOTCHAS / OPERATIONAL NOTES
- **Euler `/tmp` is node-local** + each SSH may hit a different login node → don't rely on `/tmp` files across SSH calls; write to `$RB` (shared FS). Login banner breaks scp/stdin → use **base64-argv**: `B64=$(base64 <<'EOF' ... EOF | tr -d '\n'); ssh euler "echo $B64 | base64 -d | bash"`.
- **Active-uuid guard** (protect running/pending jobs from deletion): `scontrol show job -o | grep UserId=lvalsecchi`, then extract `[0-9a-f]{8}` from `JobName=` + `StdOut=` tokens. (`squeue %o` StdOut is null for stdin jobs; comma-list `scontrol show job -o <ids>` is "Invalid" — use the all-jobs form + grep.)
- **bd walltime bug:** `code/configs/boltzdesign__multif/default.yaml` defaults to `design_samples: 60` + `time: 04:00:00` (bd@60 needs ~7h → TIMEOUT). Always override `design_samples=30` for the 4h queue, or set time=08:00:00.
- Refold leaf dirs are **pre-created at launch** (run.sbatch + resolved_config.yaml) for every attempted uuid — that's why a stack shows more refold leaves than packed slices; it is NOT a leak.
- Mem: use `--mem-per-cpu=` (af3 4G, b2 3G), NOT `--mem` (Euler lua rejects `--mem` by node).

## STATE-VERIFY ONE-LINER (run first next session)
```
ssh euler 'CAMP=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free; \
for S in boltzgen__boltzif boltzgen__ligandmpnn boltzdesign__ligandmpnn rfd3__ligandmpnn halludesign__ligandmpnn; do \
  printf "%-26s packed=%d\n" "$S" "$(ls $CAMP/$S/midrun_v2/designs/*.tar.zst 2>/dev/null|wc -l)"; done; \
squeue -u lvalsecchi -h -o "%j|%T" | sed -E "s/_[0-9a-f]{8}.*//" | sort | uniq -c; \
lquota /cluster/work/wicky | grep -i file'
```
