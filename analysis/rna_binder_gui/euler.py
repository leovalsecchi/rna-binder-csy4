"""
Thin wrappers around `ssh euler '<cmd>'` for the rna_binder GUI.

Reads happen via the local eulermount filesystem (fast).
Writes/queries happen via SSH (slower but always current).
"""
from __future__ import annotations
import base64
import subprocess
import shlex
import re
from dataclasses import dataclass


SSH_HOST = "euler"  # matches your ~/.ssh/config alias
SSH_TIMEOUT = 30


# Multiplex SSH: the first call opens a master connection that later calls reuse
# for ~2 min, so squeue + refold_counts don't each pay the login/dotfile cost.
_MUX_OPTS = [
    "-o", "ConnectTimeout=10",
    "-o", "ControlMaster=auto",
    "-o", "ControlPath=/tmp/ssh-mux-rnabinder-%r@%h:%p",
    "-o", "ControlPersist=120",
]


def ssh(cmd: str, *, check: bool = True, timeout: int = SSH_TIMEOUT) -> str:
    """Run a command on Euler via ssh, return stdout."""
    full = ["ssh", *_MUX_OPTS, SSH_HOST, cmd]
    r = subprocess.run(full, capture_output=True, text=True, timeout=timeout)
    if check and r.returncode != 0:
        raise RuntimeError(f"ssh failed ({r.returncode}): {r.stderr.strip()}")
    return r.stdout


@dataclass
class SlurmJob:
    jobid: str
    name: str
    state: str
    elapsed: str
    nodelist: str
    reason: str


def squeue(user: str = "lvalsecchi") -> list[SlurmJob]:
    """Current SLURM queue for the user."""
    raw = ssh(
        f"squeue -u {shlex.quote(user)} -h "
        f"-o '%i|%j|%T|%M|%R'"
    )
    jobs = []
    for line in raw.strip().splitlines():
        if not line:
            continue
        parts = line.split("|", 4)
        if len(parts) < 5:
            continue
        jid, name, state, elapsed, reason = parts
        nodelist = reason if state == "RUNNING" else ""
        jobs.append(SlurmJob(jid.strip(), name.strip(), state.strip(),
                             elapsed.strip(), nodelist.strip(), reason.strip()))
    return jobs


def refold_counts(remote_project: str, target: str,
                  runs: list[tuple[str, str]] | None = None,
                  ) -> dict[tuple[str, str], dict[str, int]]:
    """Count refold .cif files per run with a single server-side find.

    Returns {(stack, run): {"af3": n, "b2": n}}. Far faster than walking the
    eulermount client-side (one SSH call vs thousands of sshfs round-trips).

    runs: optional [(stack, run), ...] to scope the find to just those run dirs
    (fast incremental recompute). When None, scans the whole target tree.
    """
    base = f"{remote_project.rstrip('/')}/campaigns/{target}"
    if runs:
        starts = " ".join(shlex.quote(f"{base}/{s}/{n}") for s, n in runs)
    else:
        starts = shlex.quote(base)
    find_part = r"\( -path '*/refold_af3/*' -o -path '*/refold_b2/*' -o -path '*/refold_esm/*' \)"
    # awk reads full paths and pulls <stack>/<run>/<refold_kind> regardless of
    # the find start point, so it works for both whole-tree and scoped runs.
    awk_part = (
        "awk -F/ '{for(i=1;i<=NF;i++)"
        "if($i==\"refold_af3\"||$i==\"refold_b2\"||$i==\"refold_esm\")"
        "{print $(i-2)\"/\"$(i-1)\"/\"$i; break}}'"
    )
    cmd = (
        f"find {starts} {find_part} -name '*.cif' 2>/dev/null "
        f"| {awk_part} | sort | uniq -c"
    )
    raw = _strip_banner(ssh(cmd, check=False, timeout=120))
    counts: dict[tuple[str, str], dict[str, int]] = {}
    for line in raw.strip().splitlines():
        sp = line.split(None, 1)
        if len(sp) != 2 or not sp[0].isdigit():
            continue
        n, rest = int(sp[0]), sp[1]
        parts = rest.split("/")
        if len(parts) != 3 or not parts[2].startswith("refold_"):
            continue
        stack, run, kind = parts[0], parts[1], parts[2][len("refold_"):]
        counts.setdefault((stack, run), {"af3": 0, "b2": 0, "esm": 0})[kind] = n
    return counts


@dataclass
class LaunchResult:
    ok: bool                    # ssh + launch.py exited 0 (dry-run) / launch started (real)
    returncode: int
    dry_run: bool
    jobid: str | None           # primary design jobid (None for dry-run / parse miss / detached)
    run_dir: str | None         # run dir launch.py reported (or would-create, for dry-run)
    refold_jobids: list[str]    # AF3/B2 dependent jobids (when --auto-refold)
    stdout: str                 # login banner stripped
    stderr: str
    detached: bool = False      # real launches run detached on Euler (see below)
    log_path: str | None = None # remote log file the detached launch.py writes to


# Euler prints a dotfile-seeding banner on every ssh login; it lands on stdout
# ahead of the real command output. Strip it so callers get clean launch.py text.
_BANNER_RE = re.compile(
    r"^Seeding dotfiles.*?Done\. Re-source your shell:[^\n]*\n",
    re.DOTALL | re.MULTILINE,
)


def _strip_banner(text: str) -> str:
    return _BANNER_RE.sub("", text, count=1).lstrip("\n")


def _parse_launch_stdout(stdout: str) -> tuple[str | None, str | None, list[str]]:
    """Extract (design jobid, run dir, refold jobids) from launch.py output."""
    m_job = re.search(r"^Submitted[^\n]*jobid=(\d+)", stdout, re.MULTILINE)
    jobid = m_job.group(1) if m_job else None
    m_dir = (re.search(r"^Run dir:\s*(\S+)", stdout, re.MULTILINE)
             or re.search(r"^BIF:\s*(\S+)", stdout, re.MULTILINE)
             or re.search(r"^# Would create:\s*(\S+)", stdout, re.MULTILINE))
    run_dir = m_dir.group(1) if m_dir else None
    refolds = re.findall(r"dep refold \w+: jobid=(\d+)", stdout)
    return jobid, run_dir, refolds


def build_launch_cmd(remote_cwd: str, cfg_path: str, overrides: list[tuple[str, str]],
                     label: str, auto_refold: bool, *,
                     num_runs: int = 1, dry_run: bool = False,
                     refold_af3_overrides: list[tuple[str, str]] | None = None,
                     refold_b2_overrides: list[tuple[str, str]] | None = None,
                     refold_esm_overrides: list[tuple[str, str]] | None = None) -> str:
    """Build the remote shell command that invokes launch.py.

    Returned string is what we hand to `ssh euler '<cmd>'`. Exposed so the GUI
    can preview the exact command before submitting.
    """
    parts = [
        "module purge && module load stack/2024-06 python/3.12.8 &&",
        f"cd {shlex.quote(remote_cwd)} &&",
        "python3 code/scripts/launch.py",
        shlex.quote(cfg_path),
        f"--label {shlex.quote(label)}",
    ]
    if num_runs and num_runs > 1:
        parts.append(f"--num-runs {int(num_runs)}")
    for path, val in overrides:
        parts.append(f"--override {shlex.quote(f'{path}={val}')}")
    if auto_refold:
        parts.append("--auto-refold")
        for path, val in (refold_af3_overrides or []):
            parts.append(f"--refold-af3-override {shlex.quote(f'{path}={val}')}")
        for path, val in (refold_b2_overrides or []):
            parts.append(f"--refold-b2-override {shlex.quote(f'{path}={val}')}")
        for path, val in (refold_esm_overrides or []):
            parts.append(f"--refold-esm-override {shlex.quote(f'{path}={val}')}")
    if dry_run:
        parts.append("--dry-run")
    return " ".join(parts)


def sbatch_via_launcher(remote_cwd: str, cfg_path: str, overrides: list[tuple[str, str]],
                        label: str, auto_refold: bool, dry_run: bool = False,
                        num_runs: int = 1,
                        refold_af3_overrides: list[tuple[str, str]] | None = None,
                        refold_b2_overrides: list[tuple[str, str]] | None = None,
                        refold_esm_overrides: list[tuple[str, str]] | None = None) -> LaunchResult:
    """
    Invoke launch.py on Euler. Returns a parsed LaunchResult.

    overrides: list of (path, value) e.g. [("generator_params.n_batches", "2")]
    refold_*_overrides: applied to the AF3/B2 refold configs (only when auto_refold).
    num_runs: passed through as --num-runs N (omitted when 1, launch.py's default).
    """
    cmd = build_launch_cmd(
        remote_cwd, cfg_path, overrides, label, auto_refold,
        num_runs=num_runs, dry_run=dry_run,
        refold_af3_overrides=refold_af3_overrides,
        refold_b2_overrides=refold_b2_overrides,
        refold_esm_overrides=refold_esm_overrides,
    )

    if dry_run:
        # Fast + synchronous: launch.py just prints the sbatch it would submit.
        r = subprocess.run(
            ["ssh", SSH_HOST, cmd],
            capture_output=True, text=True, timeout=120,
        )
        stdout = _strip_banner(r.stdout)
        jobid, run_dir, refolds = _parse_launch_stdout(stdout)
        return LaunchResult(
            ok=(r.returncode == 0), returncode=r.returncode, dry_run=True,
            jobid=jobid, run_dir=run_dir, refold_jobids=refolds,
            stdout=stdout, stderr=r.stderr,
        )

    # Real launch: run launch.py DETACHED on Euler (nohup, output to a log file),
    # returning as soon as it has started. A multi-run --auto-refold launch
    # submits dozens of sbatch jobs (~1-2 s each), which used to blow past the
    # 120 s ssh timeout — and the timeout KILLED launch.py mid-submission, leaving
    # a half-submitted campaign. Detaching means the local ssh client returning
    # (or timing out) can't interrupt the submission. Jobs appear in squeue as
    # they go in; the log captures launch.py's stdout/stderr.
    # cmd is base64'd into a temp script to avoid nested-quote hell, then run via
    # `bash -l` so Euler's module system is available.
    b64 = base64.b64encode(cmd.encode()).decode()
    safe_label = re.sub(r"[^A-Za-z0-9_.-]", "_", label) or "launch"
    logdir = f"{remote_cwd.rstrip('/')}/launch_logs"
    # The launch MUST be the only backgrounded statement, on its own line. If the
    # whole `setup && nohup … &` chain is backgrounded on one line, `&` binds to
    # the ENTIRE chain — the subshell then foreground-runs launch.py and holds the
    # ssh channel open until it finishes (minutes), so ssh hangs and times out.
    # Keeping setup+echo on the first line (foreground) and `nohup … &` on the
    # second makes ssh return in ~seconds while launch.py runs detached.
    wrap = (
        f"mkdir -p {shlex.quote(logdir)} && "
        "SCRIPT=$(mktemp) && "
        f'printf %s {shlex.quote(b64)} | base64 -d > "$SCRIPT" && '
        f'LOG={shlex.quote(logdir)}/{safe_label}_$(date +%Y%m%d_%H%M%S)_$$.log && '
        'echo "LAUNCHED log=$LOG"\n'
        'nohup bash -l "$SCRIPT" > "$LOG" 2>&1 < /dev/null &\n'
    )
    r = subprocess.run(
        ["ssh", SSH_HOST, wrap],
        capture_output=True, text=True, timeout=60,
    )
    out = _strip_banner(r.stdout)
    m = re.search(r"LAUNCHED log=(\S+)", out)
    log_path = m.group(1) if m else None
    return LaunchResult(
        ok=(r.returncode == 0 and log_path is not None),
        returncode=r.returncode,
        dry_run=False,
        jobid=None, run_dir=None, refold_jobids=[],
        stdout=out, stderr=r.stderr,
        detached=True, log_path=log_path,
    )


def scancel(jobid: str) -> str:
    return ssh(f"scancel {shlex.quote(jobid)}", check=False)
