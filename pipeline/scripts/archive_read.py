"""Shared structure-path resolver for archived (.tar.zst) runs.

The aggregator records a refold/design path as "<dir>.tar.zst!<member>" when the
unit is packed (analysis/aggregate.py `_stable_ref`). Tools that read those
structures by path — dockq_metrics.py, hbonds_pyrosetta.py — call `resolve()`
to get a readable file: a plain path passes straight through; a packed ref is
expanded ONCE per archive to node-local scratch (cached) and the member path
returned.

zstd: the containers (boltzgen.sif, pyrosetta.sif) ship no zstd, so we prepend a
staged host copy at code/bin/zstd to PATH; on the host DockQ venv the system zstd
is used. Raises a clear error if neither is available.
"""
from __future__ import annotations

import atexit
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

_BUNDLED_BIN = "/cluster/work/wicky/lvalsecchi/rna_binder/code/bin"
_SEP = ".tar.zst!"
_CACHE: dict = {}          # archive path -> extracted scratch dir (insertion-ordered = LRU)
_SCRATCH = None
_SEQ = 0                   # monotonic dest-dir counter (survives cache eviction; no name reuse)
# Bound node-local scratch: keep at most this many extracted archives per process,
# evicting (rmtree) the least-recently-used. Unbounded caching filled /local/scratch
# to ~800GB on 2026-07-23 (every touched archive a full ~12GB expand, never reclaimed
# until process exit). Default 3 leaves headroom for the <=2 archives one work item
# (model+native) holds live, so eviction never rmtrees a path still in use.
_CACHE_MAX = max(2, int(os.environ.get("ARCH_READ_CACHE_MAX", "3") or "3"))


def _node_local_base():
    """Node-local scratch dir for extraction. MUST NOT be under /cluster
    (shared Lustre counts against the project inode quota and is never
    auto-purged). Prefers SLURM's per-job $TMPDIR, then node-local /scratch
    or /tmp; returns None so tempfile falls back to its own default."""
    for cand in (os.environ.get("TMPDIR"), "/scratch", "/tmp"):
        if cand and os.path.isdir(cand) and not os.path.realpath(cand).startswith("/cluster"):
            return cand
    return None


def _scratch_root() -> str:
    global _SCRATCH
    if _SCRATCH is None:
        # Extract to node-local scratch (see _node_local_base): SLURM purges
        # it at job end and it does NOT count against the shared project
        # inode quota, so even a SIGKILL that skips the atexit cleanup below
        # cannot leak onto /cluster (the 2026-07-10 arch_read_* blowup).
        _SCRATCH = tempfile.mkdtemp(prefix="arch_read_", dir=_node_local_base())
        atexit.register(shutil.rmtree, _SCRATCH, ignore_errors=True)
    return _SCRATCH


def _env() -> dict:
    e = dict(os.environ)
    e["PATH"] = _BUNDLED_BIN + os.pathsep + e.get("PATH", "")
    return e


def resolve(ref):
    """Readable file path for `ref`. Plain paths pass through unchanged; a
    "<archive>.tar.zst!<member>" ref expands its archive once (cached) and
    returns the member's extracted path."""
    if ref is None:
        return ref
    ref = str(ref)
    i = ref.find(_SEP)
    if i == -1:
        return ref                                   # normal path
    archive = ref[:i + len(".tar.zst")]
    member = ref[i + len(_SEP):]
    global _SEQ
    cached_dir = _CACHE.get(archive)
    if cached_dir is not None:
        _CACHE.pop(archive); _CACHE[archive] = cached_dir      # LRU touch (mark newest)
        return os.path.join(cached_dir, member)
    env = _env()
    if shutil.which("zstd", path=env["PATH"]) is None:
        raise RuntimeError(
            f"need zstd to read {archive} but none on PATH and none at "
            f"{_BUNDLED_BIN}/zstd")
    _SEQ += 1
    dest = Path(_scratch_root()) / f"a{_SEQ}"
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["bash", "-c",
         f"zstd -d -c {shlex.quote(archive)} | tar -x -C {shlex.quote(str(dest))}"],
        check=True, env=env,
    )
    cached_dir = str(dest)
    _CACHE[archive] = cached_dir
    while len(_CACHE) > _CACHE_MAX:                             # evict LRU -> reclaim scratch
        old_a, old_d = next(iter(_CACHE.items()))
        _CACHE.pop(old_a)
        shutil.rmtree(old_d, ignore_errors=True)
    return os.path.join(cached_dir, member)
