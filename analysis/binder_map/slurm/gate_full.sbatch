#!/bin/bash
# Task 9: hard gate between stage 2 and stage 3 for the `full` map.
#   sbatch --job-name=full_gate --dependency=afterok:$J2 slurm/gate_full.sbatch
#
# Submitted by submit_maps.sh with GATE=1; stage 3's 512-task, ~2,100 CPU-h
# array hangs off THIS job with --dependency=afterok. That is the whole design:
# if anything below exits non-zero, crashes, or the job times out, stage 3 sits
# in DependencyNeverSatisfied and Slurm reaps it. `set -e` is therefore load
# bearing -- do NOT add `|| true`, `|| echo ...`, or any handler that could let
# a failed check leave this job with status 0.
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
# Node-local scratch reservation. The 150k arm builds three sub-DBs and an
# inode-heavy foldseek tmp tree under $TMPDIR; stage1_extract_nested.sbatch
# reserves the same 20 GB for less. A node whose /scratch fills mid-run writes
# blank-but-right-length output (the July 2026 incident) -- exactly the class of
# silent failure this gate exists to stop.
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%j.log
set -euo pipefail
# `set -e` does NOT survive into $( ) subshells by itself. Without this, a
# foldseek failure inside probe_arm falls through to the function's closing
# `echo`, whose status 0 becomes the status of `T50=$(probe_arm 50000)`, and the
# gate releases a 2,100 CPU-h array having measured nothing. Belt and braces:
# every command inside probe_arm also carries an explicit `|| exit 1`, which is
# what the already-proven stage3_search.sbatch does. Requires bash >= 4.4; on an
# older shell `shopt` itself fails and `set -e` kills the job, which is the
# correct (closed) direction.
shopt -s inherit_errexit

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
cd "$FM"

echo "=== gate_full on $(hostname) at $(date) ==="

# --- 1. correctness checks ---------------------------------------------------
# Reads tar members in memory; expands nothing under results/ (hard inode cap).
PYTHONPATH=src "$PY" scripts/gate_checks.py \
  --shared "$FM/results_shared" \
  --ref "$FM/results_v6" \
  --manifest "$FM/build_shared/manifests/design.tsv" \
  --tree "$FM/results_full" \
  --tree "$FM/results_lib12k"

# --- 2. scaling probe --------------------------------------------------------
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960
command -v foldseek >/dev/null

OUT="$FM/results_full/design"
WORK="${TMPDIR:-/tmp}/gate_full_${SLURM_JOB_ID:-$$}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# A short lookup would silently shrink an arm's target set and produce a
# meaningless exponent, so refuse rather than measure the wrong thing.
N_LOOKUP=$(wc -l < "$OUT/db/uniq.lookup")
echo "uniq.lookup entries: $N_LOOKUP"
if [ "$N_LOOKUP" -lt 150000 ]; then
  echo "FATAL: uniq.lookup holds $N_LOOKUP entries, need >= 150000 for the probe"
  exit 1
fi

# Body copied from probe_scaling.sh: sub-DBs over the "" _ss _ca suffix triple,
# 2000 queries, frozen search flags. Everything but the elapsed seconds goes to
# stderr so the command substitution captures only the timing.
#
# Every step is guarded: a step that fails, or succeeds while producing nothing,
# must abort the arm. An unguarded failure here is indistinguishable from a fast
# search, and a fast search is what "the full run is cheap, go ahead" looks like.
# tests/test_gate_shell.py runs this function against a stub foldseek and
# asserts it cannot release the array.
probe_arm () {
  local N=$1
  local ARM="$WORK/n$N"
  # Under this floor an arm did not search: 2000 TM-align queries against 50k+
  # targets take minutes on 16 cores (the reference run: 189 s and 570 s).
  local min_s=5
  mkdir -p "$ARM" || exit 1
  head -n "$N" "$OUT/db/uniq.lookup" | cut -f1 > "$ARM/keys.txt" || exit 1
  head -n 2000 "$ARM/keys.txt" > "$ARM/qkeys.txt" || exit 1
  local suffix
  for suffix in "" _ss _ca; do
    foldseek createsubdb "$ARM/keys.txt"  "$OUT/db/uniq$suffix" "$ARM/sub$suffix" >&2 || exit 1
    foldseek createsubdb "$ARM/qkeys.txt" "$OUT/db/uniq$suffix" "$ARM/q$suffix"   >&2 || exit 1
  done
  # stage3_search.sbatch:34's precedent: createsubdb can exit 0 and still leave
  # nothing usable behind, and a search over an empty query DB is fast and happy.
  [ -s "$ARM/sub_ss.index" ] || { echo "FATAL: probe arm N=$N: target sub-DB not built" >&2; exit 1; }
  [ -s "$ARM/q_ss.index" ]   || { echo "FATAL: probe arm N=$N: query sub-DB not built"  >&2; exit 1; }
  echo "=== probe arm N=$N : 2000 queries against a $N-target DB ===" >&2
  local t0 t1 elapsed
  t0=$(date +%s)
  foldseek search "$ARM/q" "$ARM/sub" "$ARM/aln" "$ARM/tmp" \
    -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16 >&2 || exit 1
  t1=$(date +%s)
  [ -s "$ARM/aln.index" ] || { echo "FATAL: probe arm N=$N: search wrote no alignment db" >&2; exit 1; }
  elapsed=$(( t1 - t0 ))
  if [ "$elapsed" -lt "$min_s" ]; then
    echo "FATAL: probe arm N=$N returned in ${elapsed}s, under the ${min_s}s floor -- it cannot have searched $N targets" >&2
    exit 1
  fi
  echo "$elapsed"
}

T50=$(probe_arm 50000)
T150=$(probe_arm 150000)
echo "probe timings: t50=${T50}s t150=${T150}s"

# An unparseable or zero timing must stop the chain, not release it.
for T in "$T50" "$T150"; do
  case "$T" in
    ''|*[!0-9]*) echo "FATAL: unparseable probe timing '$T'"; exit 1 ;;
  esac
  if [ "$T" -le 0 ]; then
    echo "FATAL: probe timing '$T' is not positive"
    exit 1
  fi
done

# --- 3. projection against the ceiling ---------------------------------------
"$PY" scripts/probe_project.py --t50 "$T50" --t150 "$T150" \
  --ceiling "${GATE_CPUH_CEILING:-6000}"

echo "=== gate_full passed at $(date) ==="
