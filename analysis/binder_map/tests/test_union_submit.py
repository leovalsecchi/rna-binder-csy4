"""Tests for submit_union.sh -- the union build driver (stages 2-5, fold axis).

Nothing here submits a job: `sbatch` is a stub on PATH. The stub records its own
argv and hands back a fresh fake job id per call, so every assertion below is
about what WOULD have been submitted, and about the ids the chain was wired
with.

Three things make these discriminating rather than decorative:

* the dependency assertions compare the ids the stub actually issued, so a chain
  wired to the wrong stage fails -- and so does `--dependency=afterok:` with an
  empty right-hand side, which is exactly what a swallowed `sbatch` failure
  produces;
* the failure-injection tests fail the Nth `sbatch` and assert that nothing
  after it is submitted. `set -e` does not cover every shape of this by itself
  (a failure inside a multi-command $( ) capture never reaches the parent unless
  `inherit_errexit` is set), so the shell preamble is pinned textually too;
* stage 3 is asserted to be submitted as `slurm/stage3_search.sbatch`, not
  re-wrapped, because that file's node-local staging of the target DB is what
  stopped 265 concurrent readers from collapsing the Lustre client on
  2026-09-25.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SUBMIT_UNION = ROOT / "submit_union.sh"
STAGE3_SBATCH = ROOT / "slurm" / "stage3_search.sbatch"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None,
                                reason="bash not available")


# --------------------------------------------------------------------------- #
# harness
# --------------------------------------------------------------------------- #
STUB_SBATCH = """#!/bin/bash
# File-backed counter: each sbatch is its own process, so a shell variable would
# not survive and every fake job id would collide -- which would make the
# dependency assertions vacuous.
n=$(cat "$SBATCH_COUNTER")
n=$((n + 1))
printf '%s' "$n" > "$SBATCH_COUNTER"
if [ -n "${SBATCH_FAIL_AT:-}" ] && [ "$n" = "$SBATCH_FAIL_AT" ]; then
  echo "stub sbatch: refusing submission $n" >&2
  exit 1
fi
id=$((1000 + n))
printf '%s\\t%s\\n' "$id" "$*" >> "$SBATCH_LOG"
printf '%s\\n' "$id"
"""


def _run_with_stub(tmp_path):
    """The plainest possible run: a stub that only echoes a job id.

    Uses the driver's real FM default, so this one also covers the production
    path through the preflight.
    """
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "sbatch").write_text("#!/bin/bash\necho 12345\n")
    (stub / "sbatch").chmod(0o755)
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    return subprocess.run(["bash", str(SUBMIT_UNION)],
                          capture_output=True, text=True, env=env, cwd=tmp_path)


def _fake_fm(tmp_path, *, shards=2, hashes=2):
    """A results_union tree with just enough in it to satisfy the preflight."""
    fm = tmp_path / "fm"
    design = fm / "results_union" / "design"
    (design / "shards").mkdir(parents=True)
    (design / "hashes").mkdir(parents=True)
    for i in range(shards):
        (design / "shards" / f"shard_shared_{i:02d}.tar.gz").write_bytes(b"")
    for i in range(hashes):
        (design / "hashes" / f"shard_shared_{i:02d}.tsv").write_text("")
    return fm


def _run(tmp_path, *, fm=None, fail_at=None):
    """Run the real driver against the recording stub. Returns (proc, rows)."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "sbatch"
    stub.write_text(STUB_SBATCH)
    stub.chmod(0o755)

    log = tmp_path / "sbatch.log"
    log.write_text("")
    counter = tmp_path / "counter"
    counter.write_text("0")

    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["SBATCH_LOG"] = str(log)
    env["SBATCH_COUNTER"] = str(counter)
    env["FM"] = str(fm if fm is not None else _fake_fm(tmp_path))
    if fail_at is not None:
        env["SBATCH_FAIL_AT"] = str(fail_at)

    proc = subprocess.run(["bash", str(SUBMIT_UNION)], env=env,
                          capture_output=True, text=True, cwd=tmp_path, timeout=120)
    rows = [line.split("\t", 1) for line in log.read_text().splitlines() if line]
    return proc, [(job_id, args) for job_id, args in rows]


def _one(rows, needle):
    hits = [row for row in rows if needle in row[1]]
    assert len(hits) == 1, f"expected exactly one {needle!r} submission, got {len(hits)}"
    return hits[0]


# --------------------------------------------------------------------------- #
# the axis, and the frozen parameters
# --------------------------------------------------------------------------- #
def test_submits_fold_axis_only(tmp_path):
    out = _run_with_stub(tmp_path).stdout
    assert "fold" in out
    assert "epitope" not in out


def test_frozen_parameters_present():
    text = SUBMIT_UNION.read_text()
    for token in ("--tm-floor 0.5", "--epi-floor 0.6", "--k 30", "--mutual", "ALN_TYPE=1"):
        assert token in text, token


def test_stage5_runs_once_on_the_fold_axis(tmp_path):
    _, rows = _run(tmp_path)
    stage5 = _one(rows, "scripts/stage5_cluster.py")
    assert "--axis fold" in stage5[1]
    assert not [row for row in rows if "--axis epitope" in row[1]], \
        "the epitope axis is out of scope as of 2026-09-28"


# --------------------------------------------------------------------------- #
# error handling: the reason this driver exists in this shape
# --------------------------------------------------------------------------- #
def test_shell_options_are_the_hardened_ones():
    """`set -e` alone does not survive a $( ) capture; both lines are the fix."""
    text = SUBMIT_UNION.read_text()
    assert "set -euo pipefail" in text
    assert "shopt -s inherit_errexit" in text


def test_every_sbatch_capture_carries_its_own_exit_guard():
    text = SUBMIT_UNION.read_text()
    captures = text.count("$(sbatch")
    assert captures == 4, f"expected 4 sbatch captures, found {captures}"
    assert text.count(") || exit 1") == captures


@pytest.mark.parametrize("fail_at", [1, 2, 3, 4])
def test_a_failed_submission_stops_the_chain(tmp_path, fail_at):
    """Nothing may be submitted after a refused sbatch.

    Without the guards the driver carries on with an empty job id and submits
    the rest of the chain with `--dependency=afterok:` -- jobs that either never
    run or, worse, run immediately against a tree stage 2 never built.
    """
    proc, rows = _run(tmp_path, fail_at=fail_at)
    assert proc.returncode != 0, proc.stdout + proc.stderr
    # Proves the driver really reached sbatch, so fail_at=1 cannot pass merely
    # because the script blew up earlier -- or does not exist at all.
    assert "stub sbatch: refusing submission" in proc.stderr, proc.stderr
    assert len(rows) == fail_at - 1, [row[1] for row in rows]


def test_refuses_to_submit_when_the_union_tree_is_missing(tmp_path):
    fm = tmp_path / "empty_fm"
    fm.mkdir()
    proc, rows = _run(tmp_path, fm=fm)
    assert proc.returncode != 0, proc.stdout + proc.stderr
    assert "FATAL" in proc.stderr
    assert rows == [], "nothing may be submitted without task 1's tree"


def test_refuses_to_submit_when_the_hashes_are_missing(tmp_path):
    """Shards without hashes is the 2026-09-24 failure mode, not an empty tree."""
    fm = _fake_fm(tmp_path, shards=2, hashes=0)
    proc, rows = _run(tmp_path, fm=fm)
    assert proc.returncode != 0, proc.stdout + proc.stderr
    assert "FATAL" in proc.stderr
    assert rows == []


def test_a_complete_union_tree_is_accepted(tmp_path):
    """Positive control: the two refusals above must not pass vacuously."""
    proc, rows = _run(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert len(rows) == 4, [row[1] for row in rows]


def test_the_driver_parses():
    proc = subprocess.run(["bash", "-n", str(SUBMIT_UNION)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# --------------------------------------------------------------------------- #
# the chain
# --------------------------------------------------------------------------- #
def test_only_stages_2_to_5_are_submitted(tmp_path):
    """Task 1 produced the tree, so there is no stage 1 or 1b to run."""
    _, rows = _run(tmp_path)
    assert len(rows) == 4, [row[1] for row in rows]
    for banned in ("stage1_extract", "stage1b_split", "gate_full.sbatch"):
        assert not [row for row in rows if banned in row[1]], banned


def test_the_chain_is_wired_s2_s3_s4_s5(tmp_path):
    proc, rows = _run(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    s2 = _one(rows, "--job-name=union_s2")
    s3 = _one(rows, "--job-name=union_s3")
    s4 = _one(rows, "--job-name=union_s4")
    s5 = _one(rows, "--job-name=union_s5_fold")

    # Stage 2 starts immediately: its input is already on disk.
    assert "--dependency" not in s2[1]
    assert f"--dependency=afterok:{s2[0]}" in s3[1]
    assert f"--dependency=afterok:{s3[0]}" in s4[1]
    assert f"--dependency=afterok:{s4[0]}" in s5[1]
    # ... and no stage hangs off a stage it does not follow.
    assert f"--dependency=afterok:{s2[0]}" not in s4[1]
    assert f"--dependency=afterok:{s3[0]}" not in s5[1]


def test_every_stage_targets_the_union_tree(tmp_path):
    _, rows = _run(tmp_path)
    assert len(rows) == 4, [row[1] for row in rows]
    for _, args in rows:
        assert "results_union" in args, args
        for other in ("results_full", "results_shared", "results_lib12k",
                      "results_v6", "results_sw_"):
            assert other not in args, f"{other} in {args}"


# --------------------------------------------------------------------------- #
# per-stage shapes
# --------------------------------------------------------------------------- #
def test_stage2_builds_the_db_then_splits_the_queries(tmp_path):
    _, rows = _run(tmp_path)
    args = _one(rows, "--job-name=union_s2")[1]
    assert "RESULTS=results_union bash scripts/stage2_builddb.sh design" in args
    assert "scripts/stage3_split_queries.py" in args
    assert "--n-chunks 512" in args
    # The split must run only if the DB build succeeded, or stage 3 searches a
    # stale chunk list.
    assert "&&" in args


def test_stage3_carries_the_measured_shape(tmp_path):
    """56 s elapsed / 927 MB peak RSS, one-chunk test 2026-09-25.

    The previous 04:00:00 / 4G-per-cpu request throttled placement badly for no
    benefit, so these are pinned.
    """
    _, rows = _run(tmp_path)
    args = _one(rows, "--job-name=union_s3")[1]
    assert "--array=0-511%128" in args
    assert "--time=00:20:00" in args
    assert "--mem-per-cpu=1G" in args
    assert "--export=ALL,ALN_TYPE=1,RESULTS=results_union" in args
    assert args.rstrip().endswith("slurm/stage3_search.sbatch design")


def test_stage3_array_width_matches_the_query_split(tmp_path):
    """An array narrower than the split silently drops query chunks."""
    _, rows = _run(tmp_path)
    chunks = int(re.search(r"--n-chunks (\d+)", _one(rows, "--job-name=union_s2")[1]).group(1))
    last = int(re.search(r"--array=0-(\d+)%", _one(rows, "--job-name=union_s3")[1]).group(1))
    assert last == chunks - 1


def test_stage3_is_the_shipped_sbatch_not_a_wrap(tmp_path):
    """Re-wrapping stage 3 would lose its node-local staging of the target DB.

    2026-09-25: 265 array tasks reading the target DB off the shared filesystem
    collapsed the Lustre client and not one of 512 tasks finished in 105 min.
    """
    _, rows = _run(tmp_path)
    assert "--wrap" not in _one(rows, "--job-name=union_s3")[1]
    stage3 = STAGE3_SBATCH.read_text()
    assert 'cp -L "$OUT"/db/uniq*' in stage3
    assert "TMPDIR" in stage3


def test_stage4_carries_the_frozen_graph_parameters(tmp_path):
    _, rows = _run(tmp_path)
    args = _one(rows, "--job-name=union_s4")[1]
    assert "scripts/stage4_graph.py" in args
    assert "--layer design" in args
    for token in ("--tm-floor 0.5", "--epi-floor 0.6", "--k 30", "--mutual"):
        assert token in args, token
