#!/usr/bin/env python3
"""Patch fresh midrun_v9 HD run.sbatch to the PROVEN v7 config (idempotent):
 (1) per-task TORCH_EXTENSIONS_DIR  -> evoformer_attn JIT hang fix (shared torch_ext deadlock)
 (2) --resume + comment `rm -rf $OUTDIR` -> bulk preemption safety (keep partial recycle_*/)
Anchored .replace() (exact once each); verifies every anchor hit or reports the file as FAILED.
"""
import glob, sys

HD = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn"
RC = "/cluster/work/wicky/lvalsecchi/rna_binder/containers"

HDTE_BLOCK = (
    '\nmkdir -p "$OUTDIR"\n'
    "\nHDTE=%s/torch_ext_pertask/${SLURM_JOB_ID:-$$}\n"
    'rm -rf "$HDTE"; mkdir -p "$HDTE"\n'
    'cp -a %s/runtime_cache/torch_extensions/. "$HDTE"/ 2>/dev/null || true\n'
) % (RC, RC)

files = sorted(glob.glob("%s/midrun_v9_*/run.sbatch" % HD))
print("found %d midrun_v9 HD run.sbatch" % len(files))
patched = skipped = failed = 0
fails = []
for f in files:
    txt = open(f).read()
    if "TORCH_EXTENSIONS_DIR=$HDTE" in txt:
        skipped += 1
        continue
    orig = txt
    # 1) comment rm -rf OUTDIR (resume: keep partials)
    txt = txt.replace(
        '\nrm -rf "$OUTDIR"\n',
        '\n# rm -rf "$OUTDIR"  # resume: keep partial recycle_*/ across requeue\n', 1)
    # 2) per-task HDTE block right after mkdir -p OUTDIR
    txt = txt.replace('\nmkdir -p "$OUTDIR"\n', HDTE_BLOCK, 1)
    # 3) TORCH_EXTENSIONS_DIR env just before $CONTAINER
    txt = txt.replace('\n  $CONTAINER \\\n',
                      '\n  --env TORCH_EXTENSIONS_DIR=$HDTE \\\n  $CONTAINER \\\n', 1)
    # 4) --resume on the python call
    txt = txt.replace('\n    python ./HalluDesign_run.py \\\n',
                      '\n    python ./HalluDesign_run.py \\\n      --resume \\\n', 1)
    # 5) HDTE cleanup at end
    txt = txt.rstrip("\n") + '\nrm -rf "$HDTE" 2>/dev/null || true\n'

    ok = ('TORCH_EXTENSIONS_DIR=$HDTE' in txt
          and '--resume \\' in txt
          and '# rm -rf "$OUTDIR"' in txt
          and 'HDTE=%s/torch_ext_pertask' % RC in txt
          and txt != orig)
    if not ok:
        failed += 1; fails.append(f); continue
    open(f, "w").write(txt)
    patched += 1

print("patched=%d skipped(already)=%d FAILED=%d" % (patched, skipped, failed))
for x in fails[:10]:
    print("  FAILED:", x)
sys.exit(1 if failed else 0)
