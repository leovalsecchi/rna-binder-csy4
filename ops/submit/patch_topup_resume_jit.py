#!/usr/bin/env python3
"""Patch mv9 HD topup run.sbatch to match the original v9 recipe.

Four deltas (topup was generated without them):
  1. do NOT wipe OUTDIR on restart  -> partial recycle_*/ survive requeue
  2. per-task TORCH_EXTENSIONS_DIR  -> no shared DeepSpeed evoformer JIT lock
  3. --resume passed to HalluDesign -> continue instead of restarting
  4. trap-based cleanup of the per-task dir (better than the original's
     trailing rm, which never runs on TIMEOUT/preemption -> 18k leaked dirs)
Idempotent: files already carrying TORCH_EXTENSIONS_DIR are skipped.
"""
import sys, pathlib

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
HD = f"{RB}/campaigns/csy4_crrna_free/halludesign__ligandmpnn"
MAN = f"{RB}/ops/submit/mv9_hd_topup_seed_manifest.txt"
APPLY = "--apply" in sys.argv

OLD_A = 'rm -rf "$OUTDIR"\nmkdir -p "$OUTDIR"\n\napptainer run --nv \\\n'
NEW_A = ('# rm -rf "$OUTDIR"  # resume: keep partial recycle_*/ across requeue\n'
         'mkdir -p "$OUTDIR"\n'
         '\n'
         'HDTE=' + RB + '/containers/torch_ext_pertask/${SLURM_JOB_ID:-$$}\n'
         'trap \'rm -rf "$HDTE" 2>/dev/null || true\' EXIT TERM INT\n'
         'rm -rf "$HDTE"; mkdir -p "$HDTE"\n'
         'cp -a ' + RB + '/containers/runtime_cache/torch_extensions/. "$HDTE"/ 2>/dev/null || true\n'
         '\n'
         'apptainer run --nv \\\n')

OLD_B = '  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \\\n  $CONTAINER \\\n'
NEW_B = ('  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \\\n'
         '  --env TORCH_EXTENSIONS_DIR=$HDTE \\\n'
         '  $CONTAINER \\\n')

OLD_C = '    python ./HalluDesign_run.py \\\n      --HalluDesign_model \\\n'
NEW_C = '    python ./HalluDesign_run.py \\\n      --resume \\\n      --HalluDesign_model \\\n'

uuids = []
for line in open(MAN):
    f = line.split()
    if len(f) >= 2 and len(f[1]) == 8:
        try:
            int(f[1], 16); uuids.append(f[1])
        except ValueError:
            pass
uuids = sorted(set(uuids))

done = skip = missing = fail = 0
for u in uuids:
    p = pathlib.Path(f"{HD}/midrun_v9_{u}/run.sbatch")
    if not p.is_file():
        missing += 1; continue
    s = p.read_text()
    if "TORCH_EXTENSIONS_DIR" in s:
        skip += 1; continue
    if OLD_A not in s or OLD_B not in s or OLD_C not in s:
        fail += 1
        if fail <= 3:
            print(f"  ANCHOR MISS {u}: A={OLD_A in s} B={OLD_B in s} C={OLD_C in s}")
        continue
    s = s.replace(OLD_A, NEW_A, 1).replace(OLD_B, NEW_B, 1).replace(OLD_C, NEW_C, 1)
    if APPLY:
        p.write_text(s)
    done += 1

print(f"uuids={len(uuids)} patched={done} already={skip} no_file={missing} anchor_fail={fail} apply={APPLY}")
