"""Tests for the gate's SHELL logic: slurm/gate_full.sbatch and submit_maps.sh.

tests/test_gate.py covers the two Python modules. The riskiest file in the gate
is neither of them -- it is the sbatch script, where `set -e` does not behave the
way it reads:

    T50=$(probe_arm 50000)

`set -e` is NOT inherited by a command-substitution subshell unless
`shopt -s inherit_errexit` is set, so a `foldseek` that fails inside `probe_arm`
does not abort the function. It falls through to the closing `echo $(( t1 - t0 ))`
whose status 0 becomes the status of the assignment, the parent's `set -e` never
fires, and two ~2 second "timings" project to p = log(1)/log(3) = 0 and ~2 CPU-h
-- three orders under the ceiling. The gate would then release a 512-task,
~2,100 CPU-h array having measured nothing.

So these tests run the real shell functions, extracted verbatim from the shipped
scripts, against stub `foldseek` / `sbatch` executables on PATH. Nothing here
submits a job or calls foldseek for real.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
GATE_SBATCH = REPO / "slurm" / "gate_full.sbatch"
SUBMIT_MAPS = REPO / "submit_maps.sh"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None,
                                reason="bash not available")


# --------------------------------------------------------------------------- #
# extracting the real code under test
# --------------------------------------------------------------------------- #
def _shell_options(text: str) -> list[str]:
    """The script's own `set`/`shopt` preamble.

    The harness must run under exactly the shell options the real job runs
    under, because one of those options (`inherit_errexit`) IS the fix. Taking
    them from the file rather than restating them is what makes these tests
    discriminating: delete the shopt line from the sbatch and the harness loses
    it too.
    """
    lines = []
    for line in text.splitlines():
        if line.startswith("set "):
            lines.append(line)
        elif line.startswith("shopt "):
            # bash < 4.4 (e.g. macOS /bin/bash 3.2) has no inherit_errexit; let
            # the harness still run there. Euler's bash 5.1 sets it for real,
            # and the explicit `|| exit 1` guards cover the rest either way.
            lines.append(f"{line} 2>/dev/null || :")
    return lines


def _extract_function(text: str, name: str) -> str:
    """The verbatim body of a shell function defined at column 0."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines)
              if line.startswith(f"{name} ()") or line.startswith(f"{name}()")]
    assert len(starts) == 1, f"{name}: expected exactly one definition, got {len(starts)}"
    start = starts[0]
    ends = [i for i, line in enumerate(lines) if i > start and line == "}"]
    assert ends, f"{name}: no closing brace at column 0"
    return "\n".join(lines[start:ends[0] + 1])


def _assignments(text: str, names) -> list[str]:
    """Top-level `NAME=...` lines, so the harness uses the script's own values."""
    out = []
    for line in text.splitlines():
        for name in names:
            if line.startswith(f"{name}="):
                out.append(line)
    return out


def _call_sites(text: str, prefix: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith(prefix)]


def test_the_extractor_finds_the_real_probe_arm():
    """Guard the guard: if probe_arm is renamed or inlined these tests must not
    quietly become vacuous."""
    body = _extract_function(GATE_SBATCH.read_text(), "probe_arm")
    assert "foldseek search" in body
    assert "foldseek createsubdb" in body


# --------------------------------------------------------------------------- #
# probe_arm against a stub foldseek
# --------------------------------------------------------------------------- #
STUB_FOLDSEEK = """#!/bin/bash
# $4 is the output DB prefix for both `createsubdb A B OUT` and
# `search Q T OUT TMP`, so one line covers the outputs of both.
if [ "$STUB_MODE" = fail ]; then
  sleep 1          # a failure that takes >= 1s is exactly the one that survives
  exit 1           # an unguarded `set -e` and lands a 2-second "timing"
fi
case "$1" in
  createsubdb) printf 'stub\\n' > "$4.index"; exit 0 ;;
  search)
    case "$STUB_MODE" in
      no_output) sleep 6; exit 0 ;;                            # exits 0, writes nothing
      fast_ok)   printf 'stub\\n' > "$4.index"; exit 0 ;;      # instant "search"
      slow_ok)   sleep 6; printf 'stub\\n' > "$4.index"; exit 0 ;;
    esac
    ;;
esac
exit 0
"""


def _run_probe_arm(tmp_path, mode):
    """Run the real probe_arm under the real shell options, with a stub foldseek.

    The harness prints RELEASED only if the arm returned a timing and the
    surrounding `set -e` did not fire -- i.e. only if the gate would have gone on
    to project a cost and release stage 3.
    """
    text = GATE_SBATCH.read_text()

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "foldseek"
    stub.write_text(STUB_FOLDSEEK)
    stub.chmod(0o755)

    out = tmp_path / "out"
    (out / "db").mkdir(parents=True)
    (out / "db" / "uniq.lookup").write_text(
        "".join(f"{i}\tdesign_{i:06d}\t0\n" for i in range(3000)))

    harness = tmp_path / "harness.sh"
    harness.write_text("\n".join([
        *_shell_options(text),
        f'OUT="{out}"',
        f'WORK="{tmp_path / "work"}"',
        'mkdir -p "$WORK"',
        _extract_function(text, "probe_arm"),
        "T50=$(probe_arm 50000)",
        'echo "RELEASED t50=$T50"',
        "",
    ]))

    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["STUB_MODE"] = mode
    return subprocess.run(["bash", str(harness)], env=env,
                          capture_output=True, text=True, timeout=120)


def test_probe_arm_aborts_when_foldseek_fails(tmp_path):
    """CRITICAL: a failing foldseek must not produce a timing.

    Without `inherit_errexit` (or an explicit `|| exit 1`) this passes: the
    subshell shrugs off the failure, echoes an elapsed second, and the gate
    releases the array.
    """
    proc = _run_probe_arm(tmp_path, "fail")
    assert proc.returncode != 0, f"gate would have released:\n{proc.stdout}{proc.stderr}"
    assert "RELEASED" not in proc.stdout


def test_probe_arm_aborts_when_the_search_writes_no_alignment_db(tmp_path):
    """A search that exits 0 over a broken DB returns a plausible timing and
    nothing else; the arm must notice it produced no output."""
    proc = _run_probe_arm(tmp_path, "no_output")
    assert proc.returncode != 0, f"gate would have released:\n{proc.stdout}{proc.stderr}"
    assert "RELEASED" not in proc.stdout
    assert "no alignment db" in proc.stderr


def test_probe_arm_aborts_on_an_impossibly_fast_arm(tmp_path):
    """2000 TM-align queries against 50k targets take minutes. An instant arm
    measured nothing, whatever its exit status."""
    proc = _run_probe_arm(tmp_path, "fast_ok")
    assert proc.returncode != 0, f"gate would have released:\n{proc.stdout}{proc.stderr}"
    assert "RELEASED" not in proc.stdout
    assert "floor" in proc.stderr


def test_probe_arm_returns_the_elapsed_seconds_when_the_search_really_ran(tmp_path):
    """Positive control: the guards above must not make every arm fail, or the
    three failure tests would pass for the wrong reason."""
    proc = _run_probe_arm(tmp_path, "slow_ok")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    match = re.search(r"RELEASED t50=(\d+)", proc.stdout)
    assert match, proc.stdout + proc.stderr
    assert int(match.group(1)) >= 5


# --------------------------------------------------------------------------- #
# submit_map against a stub sbatch
# --------------------------------------------------------------------------- #
STUB_SBATCH = """#!/bin/bash
# File-backed counter: submit_map runs inside $( ), so a shell variable would
# not survive back to the parent and every fake job id would collide.
n=$(cat "$SBATCH_COUNTER")
n=$((n + 1))
printf '%s' "$n" > "$SBATCH_COUNTER"
id=$((1000 + n))
printf '%s\\t%s\\n' "$id" "$*" >> "$SBATCH_LOG"
printf '%s\\n' "$id"
"""


def _run_submit_map(tmp_path, body):
    """Run the real submit_map (and, optionally, its real call sites)."""
    text = SUBMIT_MAPS.read_text()

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "sbatch"
    stub.write_text(STUB_SBATCH)
    stub.chmod(0o755)

    log = tmp_path / "sbatch.log"
    log.write_text("")
    counter = tmp_path / "counter"
    counter.write_text("0")

    harness = tmp_path / "harness.sh"
    harness.write_text("\n".join([
        *_shell_options(text),
        *_assignments(text, ("FM", "PY", "LAYER")),
        _extract_function(text, "submit_map"),
        "J1B=999",
        body,
        "",
    ]))

    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["SBATCH_LOG"] = str(log)
    env["SBATCH_COUNTER"] = str(counter)
    proc = subprocess.run(["bash", str(harness)], env=env,
                          capture_output=True, text=True, timeout=120)
    rows = [line.split("\t", 1) for line in log.read_text().splitlines() if line]
    return proc, [(job_id, args) for job_id, args in rows]


def _one(rows, needle):
    hits = [row for row in rows if needle in row[1]]
    assert len(hits) == 1, f"expected exactly one {needle!r} submission, got {len(hits)}"
    return hits[0]


REAL_CALL_SITES = "\n".join([
    "J5_LIB=$(submit_map  lib12k results_lib12k 4   4   $J1B 0) || exit 1",
    "J5_FULL=$(submit_map full   results_full   512 256 $J1B 1) || exit 1",
])


def test_the_real_call_sites_are_what_the_script_ships():
    """The harness below runs the call sites verbatim, so pin them to the file.

    This also pins the load-bearing `|| exit 1` guards: submit_map runs inside
    $( ), so its internal guards kill only the subshell.
    """
    lib = _call_sites(SUBMIT_MAPS.read_text(), "J5_LIB=")
    full = _call_sites(SUBMIT_MAPS.read_text(), "J5_FULL=")
    assert lib == [REAL_CALL_SITES.splitlines()[0]]
    assert full == [REAL_CALL_SITES.splitlines()[1]]


def test_only_the_full_map_is_gated_and_its_s3_hangs_off_the_gate(tmp_path):
    proc, rows = _run_submit_map(tmp_path, REAL_CALL_SITES)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    lib_s2 = _one(rows, "--job-name=lib12k_s2")
    lib_s3 = _one(rows, "--job-name=lib12k_s3")
    full_s2 = _one(rows, "--job-name=full_s2")
    gate = _one(rows, "--job-name=full_gate")
    full_s3 = _one(rows, "--job-name=full_s3")

    assert not [row for row in rows if "lib12k_gate" in row[1]], "lib12k must not be gated"
    assert "slurm/gate_full.sbatch" in gate[1]
    assert f"--dependency=afterok:{full_s2[0]}" in gate[1]

    # The point of the whole task: the expensive array waits on the GATE, not on
    # stage 2.
    assert f"--dependency=afterok:{gate[0]}" in full_s3[1]
    assert f"--dependency=afterok:{full_s2[0]}" not in full_s3[1]
    # ... while the cheap one still hangs straight off its own stage 2.
    assert f"--dependency=afterok:{lib_s2[0]}" in lib_s3[1]


def test_the_arrays_are_the_widths_the_call_sites_ask_for(tmp_path):
    """Chunking is a pure work partition; the array geometry is still worth
    pinning, because a throttle typo is silent and expensive."""
    _, rows = _run_submit_map(tmp_path, REAL_CALL_SITES)
    assert "--array=0-511%256" in _one(rows, "--job-name=full_s3")[1]
    assert "--n-chunks 512" in _one(rows, "--job-name=full_s2")[1]
    assert "--array=0-3%4" in _one(rows, "--job-name=lib12k_s3")[1]
    assert "--n-chunks 4" in _one(rows, "--job-name=lib12k_s2")[1]


# NB: an EMPTY sixth argument is not in this list -- `${6:-0}` substitutes
# the default for unset *or* empty, so "" legitimately means "ungated".
@pytest.mark.parametrize("bad", ["true", "yes", "01", "2", "-1"])
def test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one(tmp_path, bad):
    """`[ "$GATE" = "1" ]` alone treats every typo as "ungated" and launches the
    array with no gate and no error."""
    proc, rows = _run_submit_map(
        tmp_path, f'submit_map full results_full 512 256 $J1B "{bad}" || exit 1')
    assert proc.returncode != 0, proc.stdout + proc.stderr
    assert "GATE must be 0 or 1" in proc.stderr
    assert rows == [], "nothing may be submitted once the flag is rejected"


@pytest.mark.parametrize("good", ["0", "1"])
def test_submit_map_accepts_both_gate_flags(tmp_path, good):
    proc, rows = _run_submit_map(
        tmp_path, f'submit_map full results_full 4 4 $J1B {good} || exit 1')
    assert proc.returncode == 0, proc.stdout + proc.stderr
    gated = [row for row in rows if "--job-name=full_gate" in row[1]]
    assert len(gated) == int(good)
