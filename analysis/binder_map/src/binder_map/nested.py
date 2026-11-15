"""Read design structures out of the post-2026-09-11 nested archives.

On 2026-09-11 every midrun directory was collapsed into one solid
`midrun_vN.tar.zst` (683 GB across 35 archives). A design that the filter store
still records as

    <stack>/<midrun>/designs/<uuid>.tar.zst!<design_id>.cif

now lives two levels down:

    <stack>/<midrun>.tar.zst ! <midrun>/designs/<uuid>.tar.zst ! <design_id>.cif

Verified 2026-09-14: all 35 midruns referenced by the store's 37 design
templates are packed, none are live, and every design template has this one
shape. 40 of 40 sampled design archives are absent at their recorded paths.

A solid zstd stream has no random access, so the only sane access pattern is one
sequential pass per outer archive, pulling everything needed from it. That is
why callers group work by outer archive rather than by design.
"""
from __future__ import annotations

import io
import subprocess
import tarfile
from pathlib import Path
from typing import Iterator

_DESIGNS = "/designs/"


def nested_ref(design_ref: str) -> tuple[str, str, str]:
    """Rewrite a store design ref into (outer, inner, member) paths."""
    archive, member = design_ref.split("!", 1)
    head, uuid_tar = archive.split(_DESIGNS, 1)
    midrun = head.rsplit("/", 1)[-1]
    return f"{head}.tar.zst", f"{midrun}{_DESIGNS}{uuid_tar}", member


def _decompress(blob: bytes) -> bytes:
    return subprocess.run(
        ["zstd", "-dc"], input=blob, stdout=subprocess.PIPE, check=True
    ).stdout


def iter_members(outer: Path, wanted: dict[str, set[str]]) -> Iterator[tuple[str, str, bytes]]:
    """Stream `outer` exactly once, yielding (inner_name, member_name, blob).

    `wanted` maps an inner archive's name to the member names needed from it.

    Inner archives are decompressed in memory and dropped immediately. Nothing
    is written to scratch: unbounded expansion of archives filled /local/scratch
    to ~800 GB in July 2026 and produced blank-but-correct-length outputs, which
    is the failure mode a resume check cannot catch.

    Members that are absent are skipped silently; the caller reconciles counts
    against its manifest and records the misses with a reason.

    Design members occupy only the first 8-10% of each outer archive (verified
    against the `.list` sidecars, which record tar member order), so once every
    inner archive named in `wanted` has been consumed this returns immediately
    rather than streaming the remaining 90%+ to exhaustion. `return` from inside
    the `with`/`try` is safe here: it runs the same subprocess teardown as normal
    exit, exceptions and generator abandonment already do (see the module tests).
    """
    remaining = set(wanted)
    if not remaining:
        return
    proc = subprocess.Popen(
        ["zstd", "-dc", "-T0", str(outer)],
        stdout=subprocess.PIPE,
        bufsize=1 << 20,
    )
    try:
        with tarfile.open(fileobj=proc.stdout, mode="r|") as outer_tar:
            for info in outer_tar:
                if not info.isfile() or info.name not in remaining:
                    continue
                need = wanted[info.name]
                raw = _decompress(outer_tar.extractfile(info).read())
                with tarfile.open(fileobj=io.BytesIO(raw), mode="r|") as inner_tar:
                    for member in inner_tar:
                        if member.isfile() and member.name in need:
                            yield info.name, member.name, inner_tar.extractfile(member).read()
                del raw
                remaining.discard(info.name)
                if not remaining:
                    return
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        proc.wait()
