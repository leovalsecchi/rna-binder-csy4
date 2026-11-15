import io
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest

import binder_map.nested as nested
from binder_map.nested import iter_members, nested_ref

pytestmark = pytest.mark.skipif(shutil.which("zstd") is None, reason="needs zstd")


def test_nested_ref_splits_a_store_design_ref():
    ref = ("campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3/designs/"
           "ab12cd34.tar.zst!csy4free_RFD3_LiMPNN_len80_ab12cd34_design0_seq0.cif")
    outer, inner, member = nested_ref(ref)
    assert outer == "campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3.tar.zst"
    assert inner == "midrun_v3/designs/ab12cd34.tar.zst"
    assert member == "csy4free_RFD3_LiMPNN_len80_ab12cd34_design0_seq0.cif"


def test_nested_ref_handles_a_suffixed_midrun():
    ref = ("campaigns/c/boltzgen__boltzif/midrun_v8_50e6da5c/designs/"
           "ff00.tar.zst!d.cif")
    outer, inner, member = nested_ref(ref)
    assert outer == "campaigns/c/boltzgen__boltzif/midrun_v8_50e6da5c.tar.zst"
    assert inner == "midrun_v8_50e6da5c/designs/ff00.tar.zst"
    assert member == "d.cif"


def _zstd(src: Path, dst: Path) -> None:
    subprocess.run(["zstd", "-q", "-f", str(src), "-o", str(dst)], check=True)


def _build_nested(tmp_path: Path) -> Path:
    """Build midrun_v9.tar.zst ! midrun_v9/designs/<u>.tar.zst ! <m>.cif"""
    work = tmp_path / "work"
    (work / "midrun_v9" / "designs").mkdir(parents=True)

    for uuid, members in (("aaaa", {"d0.cif": b"ZERO\n", "d1.cif": b"ONE\n"}),
                          ("bbbb", {"d2.cif": b"TWO\n"})):
        plain = work / f"{uuid}.tar"
        with tarfile.open(plain, "w") as tar:
            for name, data in members.items():
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
        _zstd(plain, work / "midrun_v9" / "designs" / f"{uuid}.tar.zst")
        plain.unlink()

    outer_plain = tmp_path / "midrun_v9.tar"
    with tarfile.open(outer_plain, "w") as tar:
        tar.add(work / "midrun_v9", arcname="midrun_v9")
    outer = tmp_path / "midrun_v9.tar.zst"
    _zstd(outer_plain, outer)
    return outer


def test_iter_members_pulls_only_what_was_asked_for(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d1.cif"},
              "midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    got = {(inner, member): blob for inner, member, blob in iter_members(outer, wanted)}
    assert got == {
        ("midrun_v9/designs/aaaa.tar.zst", "d1.cif"): b"ONE\n",
        ("midrun_v9/designs/bbbb.tar.zst", "d2.cif"): b"TWO\n",
    }


def test_iter_members_ignores_inner_archives_not_requested(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    inners = {inner for inner, _, _ in iter_members(outer, wanted)}
    assert inners == {"midrun_v9/designs/bbbb.tar.zst"}


def test_iter_members_is_silent_about_absent_members(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"nope.cif"}}
    assert list(iter_members(outer, wanted)) == []


def _inner_tar_zst(members: dict) -> bytes:
    """Build one in-memory <uuid>.tar.zst blob: the payload an outer member carries."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return subprocess.run(
        ["zstd", "-q", "-c"], input=buf.getvalue(), stdout=subprocess.PIPE, check=True
    ).stdout


def _build_outer(tmp_path: Path, entries: list) -> Path:
    """Assemble midrun_v9.tar.zst from a literal, ordered list of (name, blob)
    entries. Unlike `_build_nested`, this allows a name to repeat - which is how
    the tests below plant a POISON entry after the real one: a member reusing an
    already-satisfied inner archive's name, but holding bytes that are not valid
    zstd. `_decompress` runs `zstd -dc` with check=True, so opening it raises
    subprocess.CalledProcessError if the duplicate is ever handed to it.

    `iter_members`'s loop guard tests `info.name not in remaining`, and
    `remaining` drops a name the moment it has been fully delivered. The two
    tests below plant their poison entry in different places, and that
    placement changes what each one actually proves:

    - `test_iter_members_skips_a_duplicate_name_with_only_one_archive_wanted`
      puts the poison entry *after* the single wanted archive is already
      satisfied. By then `remaining` is empty and the early `return` fires
      on the real entry, before the loop ever reaches the poison one - so
      this test does not exercise the duplicate-name guard at all. It would
      pass unchanged even if the guard were keyed on `wanted` instead of
      `remaining` (or deleted outright), because the poison entry is never
      looked at either way.

    - `test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found`
      puts the poison entry *between* two still-outstanding wanted archives
      (`aaaa`, poison, `bbbb`), so `remaining` is still non-empty and the
      loop genuinely reaches the guard check on it. This is the one that
      discriminates: with the guard keyed on `remaining`, the poison entry
      is skipped (its name was already discarded once `aaaa` was delivered)
      and both real members still yield; with the guard keyed on `wanted`
      (the pre-fix bug), the name is still present in `wanted` unconditionally,
      so the poison entry is not skipped and `_decompress` raises
      `subprocess.CalledProcessError`.

    Neither test proves the early `return` itself fires (a guard keyed on
    `remaining` skips a duplicate whether or not `return` still exists);
    that is proven separately, by measuring how much of the outer stream
    gets read - see
    `test_iter_members_stops_reading_the_outer_stream_once_satisfied` below.
    """
    plain = tmp_path / "midrun_v9.tar"
    with tarfile.open(plain, "w") as tar:
        for name, blob in entries:
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    outer = tmp_path / "midrun_v9.tar.zst"
    _zstd(plain, outer)
    plain.unlink()
    return outer


def test_iter_members_skips_a_duplicate_name_with_only_one_archive_wanted(tmp_path):
    aaaa = _inner_tar_zst({"d1.cif": b"ONE\n"})
    outer = _build_outer(tmp_path, [
        ("midrun_v9/designs/aaaa.tar.zst", aaaa),
        ("midrun_v9/designs/aaaa.tar.zst", b"not a zstd stream"),  # poison duplicate
    ])
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d1.cif"}}
    got = list(iter_members(outer, wanted))
    assert got == [("midrun_v9/designs/aaaa.tar.zst", "d1.cif", b"ONE\n")]


def test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found(tmp_path):
    aaaa = _inner_tar_zst({"d1.cif": b"ONE\n"})
    bbbb = _inner_tar_zst({"d2.cif": b"TWO\n"})
    outer = _build_outer(tmp_path, [
        ("midrun_v9/designs/aaaa.tar.zst", aaaa),
        ("midrun_v9/designs/aaaa.tar.zst", b"not a zstd stream"),  # poison duplicate,
                                                                    # planted before bbbb
                                                                    # is satisfied
        ("midrun_v9/designs/bbbb.tar.zst", bbbb),
    ])
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d1.cif"},
              "midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    got = {(inner, member): blob for inner, member, blob in iter_members(outer, wanted)}
    assert got == {
        ("midrun_v9/designs/aaaa.tar.zst", "d1.cif"): b"ONE\n",
        ("midrun_v9/designs/bbbb.tar.zst", "d2.cif"): b"TWO\n",
    }


class _CountingStdout:
    """Wraps a Popen.stdout pipe and tallies bytes actually read off it."""

    def __init__(self, inner):
        self._inner = inner
        self.nbytes = 0

    def read(self, size=-1):
        chunk = self._inner.read(size)
        self.nbytes += len(chunk)
        return chunk

    def close(self):
        self._inner.close()


def test_iter_members_stops_reading_the_outer_stream_once_satisfied(tmp_path, monkeypatch):
    """Prove the early `return` actually curtails I/O, not just what gets
    yielded. Neither poison-duplicate test above can show this any more (see
    `_build_outer`'s docstring): once the loop guard is keyed on `remaining`,
    a duplicate is skipped by the guard whether or not `return` still exists.
    So instead measure the outer pipe directly: put the one wanted inner
    archive at the front, follow it with several MB of filler the real
    pipeline never reads (refold output, in production), and assert most of
    that filler is never pulled off `subprocess.Popen`'s stdout pipe.
    """
    aaaa = _inner_tar_zst({"d1.cif": b"ONE\n"})
    filler = b"F" * (4 * 1024 * 1024)

    plain = tmp_path / "midrun_v9.tar"
    with tarfile.open(plain, "w") as tar:
        for name, blob in [("midrun_v9/designs/aaaa.tar.zst", aaaa),
                            ("padding/filler.bin", filler)]:
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    uncompressed_size = plain.stat().st_size
    outer = tmp_path / "midrun_v9.tar.zst"
    _zstd(plain, outer)
    plain.unlink()

    readers: list[_CountingStdout] = []
    real_popen = subprocess.Popen

    def counting_popen(*args, **kwargs):
        proc = real_popen(*args, **kwargs)
        # Only wrap the OUTER `zstd -dc -T0 <outer>` pipe that iter_members
        # itself streams via Popen. `_decompress` also goes through this same
        # patched Popen (subprocess.run calls it internally, and patching
        # nested.subprocess.Popen patches the module-global `subprocess.Popen`
        # subprocess.run resolves too) - wrapping ITS stdout as well would
        # break subprocess.run's own internals (it expects a real pipe object
        # with a `.closed` attribute etc.), and it isn't the stream this test
        # cares about anyway.
        cmd = args[0] if args else kwargs.get("args")
        is_outer_stream = (proc.stdout is not None and isinstance(cmd, (list, tuple))
                            and "-T0" in cmd)
        if is_outer_stream:
            proc.stdout = _CountingStdout(proc.stdout)
            readers.append(proc.stdout)
        return proc

    monkeypatch.setattr(nested.subprocess, "Popen", counting_popen)

    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d1.cif"}}
    got = list(iter_members(outer, wanted))

    assert got == [("midrun_v9/designs/aaaa.tar.zst", "d1.cif", b"ONE\n")]
    assert len(readers) == 1, "iter_members should open exactly one outer pipe"
    assert uncompressed_size > 3 * 1024 * 1024
    assert readers[0].nbytes < uncompressed_size // 4, (
        f"read {readers[0].nbytes} of {uncompressed_size} uncompressed bytes off "
        "the outer pipe; the early return should stop well short of the filler")
