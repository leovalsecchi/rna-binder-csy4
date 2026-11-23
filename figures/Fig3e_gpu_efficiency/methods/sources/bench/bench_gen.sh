#!/usr/bin/env bash
# gpubench_v1 GEN -- per-stack GPU-efficiency bench.
#
# WHY: the campaign ledger (197,797 GPU-h) cannot rank the stacks by efficiency,
# because (a) 32% of those GPU-h ended preempt/fail/timeout, (b) stacks ran on
# different GPU generations (RFD3/BG/BD1 on Turing eu-g3, HD on Ampere/Ada), and
# (c) BoltzGen gen is ONE job feeding two stacks and was never split. This
# re-runs campaign-identical gen jobs on ONE card model so the only remaining
# variable is the stack.
#
# Model params, seeds policy, containers, cpus and mem are the PRODUCTION
# configs untouched. Only the slurm block is overridden, exactly as the campaign
# drivers did. The four deliberate deviations, each required for a clean clock:
#
#  1. ONE card model: --gpus=rtx_4090:1. All 4090s on Euler are eu-g6 (80 nodes,
#     all in gpuhe.24h), so the card model IS the node pin -- no --nodelist, which
#     also avoids the decommissioned eu-g4-032 baked into the refold configs.
#     Bench numbers will be FASTER than the ledger for RFD3/BG/BD1. That is the
#     point: it removes GPU generation as a confounder.
#  2. Non-preemptible partition only (gpuhe.24h, never .bulk). A preempted job's
#     ElapsedRaw is not a slow clock, it is a wrong clock.
#  3. --requeue / --open-mode=append stripped from EVERY stack (all four
#     templates emit them). A requeued attempt must not silently become a second
#     data point; the collector rejects multi-attempt jobs, so let them fail.
#  4. HD only: --ntasks=4 -> --cpus-per-task=4. launch.py's HD template
#     (launch.py:209) emits --ntasks with no --cpus-per-task; the campaign fixed
#     this in the array bucketer. Left alone, HD's single process is bound to
#     1 core instead of 4 and its number means nothing. This rewrite REPRODUCES
#     campaign conditions, it does not depart from them.
#
# Counts: 2 x RFD3 (640 designs each), 2 x BG (600 BIF + 600 LiMPNN each, shared
# backbone), 2 x BD1 (~600), 4 x HD (150 each = 600). Refolds are a fixed 600
# per stack and live in bench_refold.sh.
#
# Usage:  bash bench_gen.sh            # dry run, submits nothing
#         bash bench_gen.sh --submit
set -uo pipefail

RB=/cluster/work/wicky/lvalsecchi/rna_binder
BENCH="$RB/bench_gpueff"
LABEL=gpubench_v1

PART=gpuhe.24h
GPUS="rtx_4090:1"
# eu-g6-064 is on the campaign refold driver's permanent exclude list
# (mv8_refold_bulk.py:231) and it is in our pinned pool.
EXCLUDE="eu-g6-064"

HD_SEED_BASE=${HD_SEED_BASE:-900001}   # campaign used 90001-92400 (v8) and 100001-102400 (v9)
HD_UUIDS=${HD_UUIDS:-4}            # 4 x 150 = 600, matching the other stacks
REPS=${REPS:-2}                # RFD3 / BG / BD1 replicates

SUBMIT=0
[ "${1:-}" = "--submit" ] && SUBMIT=1

mkdir -p "$BENCH"
MAN="$BENCH/gen_manifest.tsv"
if [ "$SUBMIT" = 1 ] && [ ! -s "$MAN" ]; then
  printf "role\tstack\tbackend\tcfg\tuuid8\tjobid\trun_dir\tn_designs\n" > "$MAN"
fi

cd "$RB/code" || { echo "no $RB/code" >&2; exit 1; }
note(){ printf '>>> %s\n' "$*"; }

# HD ONLY: the shared DeepSpeed evoformer JIT lock.
# launch.py does NOT emit the per-task TORCH_EXTENSIONS_DIR block -- it is a
# POST-generation patch that every new HD batch has to re-apply, and omitting it
# is what wrecked one of this bench's first four HD jobs: it sat 2h at 0% GPU
# waiting on the one shared lock in containers/runtime_cache/torch_extensions/
# evoformer_attn (ds_accelerator at 16:10 for a 14:10 start), then stalled for
# good with 1 file written. The other three raced for the same lock and paid a
# ~7 min kernel rebuild the campaign's pre-seeded jobs never pay.
# Same surgery as ops/submit/patch_topup_resume_jit.py, minus --resume and minus
# the OUTDIR-preserving change: bench jobs have --requeue stripped, so there is
# nothing to resume and a wiped OUTDIR is correct for a single clean run.
# Idempotent, and it ABORTS on a missing anchor rather than submitting unpatched.
# OPT-IN (TS_STAGES=1): timestamp every line the worker prints, so the Stage
# 1 / 2 / 3 boundaries inside ONE allocation become datable.
# Why this is needed: boltzgen_multif.sh prints "=== Stage N ... ===" banners with
# NO timestamps anywhere in its log (verified: 0 timestamped lines), and the two
# inverse folders' designs are both packed at job end with identical mtimes. So
# after the fact the only datable boundary is end-of-backbone-design (the
# _backbones/spec_* mtimes). Splitting the SHARED backbone stage from each
# inverse folder's OWN cost needs the Stage 2 -> Stage 3 boundary, which only a
# timestamped run can give.
# This changes no production file and no GPU work -- it wraps the worker's stdout
# in the bench's own run.sbatch. Do NOT use such a job as a gen clock data point;
# use it only for the stage split.
ts_stages(){
  python3 - "$1" <<'PY'
import pathlib, re, sys
p = pathlib.Path(sys.argv[1]); s = p.read_text()
if "__RC=" in s:
    print("     stage timestamping already present"); raise SystemExit(0)
m = re.search(r'^bash (\S*boltzgen_multif\.sh.*)$', s, re.M)
if not m:
    raise SystemExit("no boltzgen_multif.sh invocation to wrap")
wrapped = ('{ bash ' + m.group(1) + ' ; echo "__RC=$?" ; } 2>&1 \\\n'
           '  | while IFS= read -r __l; do printf \'[%s] %s\\n\' "$(date +%T)" "$__l"; done')
p.write_text(s[:m.start()] + wrapped + s[m.end():])
print("     stage timestamping applied")
PY
}

hd_jit_fix(){
  python3 - "$1" <<'PY'
import pathlib, sys
p = pathlib.Path(sys.argv[1]); s = p.read_text()
if "TORCH_EXTENSIONS_DIR" in s:
    print("     HDTE block already present"); raise SystemExit(0)
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
old_a = 'mkdir -p "$OUTDIR"\n\napptainer run --nv \\\n'
new_a = ('mkdir -p "$OUTDIR"\n\n'
         'HDTE=' + RB + '/containers/torch_ext_pertask/${SLURM_JOB_ID:-$$}\n'
         'trap \'rm -rf "$HDTE" 2>/dev/null || true\' EXIT TERM INT\n'
         'rm -rf "$HDTE"; mkdir -p "$HDTE"\n'
         'cp -a ' + RB + '/containers/runtime_cache/torch_extensions/. "$HDTE"/ 2>/dev/null || true\n'
         '\napptainer run --nv \\\n')
old_b = '  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \\\n  $CONTAINER \\\n'
new_b = ('  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \\\n'
         '  --env TORCH_EXTENSIONS_DIR=$HDTE \\\n'
         '  $CONTAINER \\\n')
for o in (old_a, old_b):
    if s.count(o) != 1:
        raise SystemExit(f"HDTE anchor not found exactly once: {o!r}")
p.write_text(s.replace(old_a, new_a).replace(old_b, new_b))
print("     HDTE JIT fix applied")
PY
}


# launch.py writes the run dir + run.sbatch AND submits immediately. We always
# cancel that native submit and sbatch the script ourselves, so the header edits
# apply and every stack goes through one identical submit path.
launch_one(){  # stack cfg [extra overrides...]
  local stack="$1" cfg="$2"; shift 2
  local out uuid jid rdir sb newjid

  if [ "$SUBMIT" = 0 ]; then
    out=$(python3 scripts/launch.py "$cfg" --label "$LABEL" --num-runs 1 --dry-run \
          --override slurm.partition="$PART" --override slurm.gpus="$GPUS" "$@" 2>&1)
    uuid=$(printf '%s\n' "$out" | grep -oE "UUID8: [0-9a-f]{8}" | grep -oE "[0-9a-f]{8}" | head -1)
    note "DRY $stack uuid=${uuid:-PARSE_FAIL}"
    printf '%s\n' "$out" | sed -n '/^#SBATCH/p' | sed 's/^/     /'
    printf '%s\n' "$out" | grep -qE "^#SBATCH --(requeue|open-mode)" \
      && echo "     [would strip --requeue / --open-mode]"
    printf '%s\n' "$out" | grep -qE "^#SBATCH --ntasks=" \
      && echo "     [would rewrite --ntasks -> --cpus-per-task]"
    return 0
  fi

  out=$(python3 scripts/launch.py "$cfg" --label "$LABEL" --num-runs 1 \
        --override slurm.partition="$PART" --override slurm.gpus="$GPUS" "$@" 2>&1)
  uuid=$(printf '%s\n' "$out" | grep -oE "uuid8=[0-9a-f]{8}|UUID8: [0-9a-f]{8}" \
         | grep -oE "[0-9a-f]{8}" | head -1)
  jid=$(printf '%s\n' "$out" | grep -oE "jobid=[0-9]+" | grep -oE "[0-9]+" | head -1)
  if [ -z "$uuid" ]; then
    echo "   !! no uuid8 parsed for $stack -- launch.py tail:" >&2
    printf '%s\n' "$out" | tail -15 >&2
    return 1
  fi
  [ -n "$jid" ] && scancel "$jid" 2>/dev/null

  # BoltzGen multif creates TWO dirs sharing the uuid (boltzif + ligandmpnn) and
  # only one carries run.sbatch, so pick by the file, not by the first match.
  rdir=""
  while IFS= read -r cand; do
    [ -f "$cand/run.sbatch" ] && { rdir="$cand"; break; }
  done < <(find "$RB/campaigns/csy4_crrna_free" -maxdepth 2 -type d \
           -name "${LABEL}_${uuid}" 2>/dev/null | grep -v _backbones)
  sb="$rdir/run.sbatch"
  if [ -z "$rdir" ] || [ ! -f "$sb" ]; then
    echo "   !! no run.sbatch for $stack $uuid" >&2
    find "$RB/campaigns/csy4_crrna_free" -maxdepth 2 -type d \
         -name "${LABEL}_${uuid}" 2>/dev/null >&2
    return 1
  fi

  # --- the bench-only header edits (see header notes 3 and 4) ---
  sed -i \
    -e 's/^#SBATCH --ntasks=\([0-9]\+\)$/#SBATCH --cpus-per-task=\1/' \
    -e '/^#SBATCH --requeue$/d' \
    -e '/^#SBATCH --open-mode=append$/d' \
    "$sb"
  grep -q '^#SBATCH --exclude=' "$sb" || \
    sed -i "s|^\(#SBATCH --partition=.*\)$|\1\n#SBATCH --exclude=$EXCLUDE|" "$sb"

  # Refuse to submit a job whose header still carries a clock-corrupting flag.
  if grep -qE '^#SBATCH --(requeue|ntasks=)' "$sb"; then
    echo "   !! header edit failed for $stack $uuid -- NOT submitting:" >&2
    grep -E '^#SBATCH --(requeue|ntasks=)' "$sb" >&2
    return 1
  fi

  if [ "${TS_STAGES:-0}" = 1 ] && grep -q "boltzgen_multif.sh" "$sb"; then
    ts_stages "$sb" || { echo "   !! stage-timestamp patch failed -- NOT submitting" >&2; return 1; }
  fi

  if grep -q "HalluDesign_run.py" "$sb"; then
    hd_jit_fix "$sb" || { echo "   !! HD JIT patch failed for $uuid -- NOT submitting" >&2; return 1; }
  fi

  newjid=$(sbatch --parsable "$sb" 2>&1)
  if ! [[ "$newjid" =~ ^[0-9]+$ ]]; then
    echo "   !! sbatch failed for $stack $uuid: $newjid" >&2
    return 1
  fi
  printf "gen\t%s\t-\t%s\t%s\t%s\t%s\t\n" "$stack" "$cfg" "$uuid" "$newjid" "$rdir" >> "$MAN"
  note "$stack uuid=$uuid jobid=$newjid"
}

note "lane: $PART  gpus=$GPUS  exclude=$EXCLUDE  submit=$SUBMIT"

for ((r=0; r<REPS; r++)); do
  launch_one rfd3_limpnn configs/rfd3__multif/default.yaml        --override slurm.time=04:00:00
  launch_one bg_shared   configs/boltzgen__multif/default.yaml    --override slurm.time=23:30:00
  launch_one bd1_limpnn  configs/boltzdesign__multif/default.yaml --override slurm.time=23:30:00
done

for ((k=0; k<HD_UUIDS; k++)); do
  launch_one hd_limpnn configs/halludesign/random_init.yaml \
    --override slurm.time=23:30:00 \
    --override generator_params.params.seed=$((HD_SEED_BASE + k))
done

note "done. manifest: $MAN"
[ "$SUBMIT" = 0 ] && note "DRY RUN -- nothing submitted. Re-run with --submit."
exit 0
