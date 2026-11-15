import io
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import stage1_extract_nested as s1
from stage1_extract_nested import plan_archive_work, shard_outers

pytestmark = pytest.mark.skipif(shutil.which("zstd") is None, reason="needs zstd")


def _manifest() -> pd.DataFrame:
    return pd.DataFrame({
        "design_id": ["d0", "d1", "d2", "d3"],
        "outer_path": ["a.tar.zst", "a.tar.zst", "a.tar.zst", "b.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u0.tar.zst",
                       "m/designs/u1.tar.zst", "m/designs/u2.tar.zst"],
        "member_path": ["d0.cif", "d1.cif", "d2.cif", "d3.cif"],
    })


def test_shard_outers_assigns_whole_archives():
    got = shard_outers(_manifest(), shard=0, n_shards=2)
    assert set(got["outer_path"]) == {"a.tar.zst"}
    assert len(got) == 3


def test_every_archive_lands_in_exactly_one_shard():
    manifest = _manifest()
    seen = []
    for shard in range(2):
        seen += list(shard_outers(manifest, shard, 2)["design_id"])
    assert sorted(seen) == ["d0", "d1", "d2", "d3"]


def test_plan_archive_work_groups_members_by_inner_archive():
    wanted = plan_archive_work(_manifest().query("outer_path == 'a.tar.zst'"))
    assert wanted == {
        "m/designs/u0.tar.zst": {"d0.cif", "d1.cif"},
        "m/designs/u1.tar.zst": {"d2.cif"},
    }


def test_plan_archive_work_on_empty_group_is_empty():
    empty = _manifest().iloc[0:0]
    assert plan_archive_work(empty) == {}


# --------------------------------------------------------------------------
# main()-level integration tests.
#
# main() drives real chemistry (biotite CIF parsing, chain selection,
# fingerprinting) that these tests have no interest in exercising - that is
# covered elsewhere (test_chains.py, test_fingerprint.py, ...) and by the
# brief's own real-archive smoke test. So the structural-biology calls
# main() makes are stubbed out here with simple, content-driven fakes; what
# these tests actually exercise is main()'s own bookkeeping: which manifest
# row a delivered blob gets attributed to, how a corrupt archive is
# recovered from, and what gets written where and when.

def _zstd(src: Path, dst: Path) -> None:
    subprocess.run(["zstd", "-q", "-f", str(src), "-o", str(dst)], check=True)


def _inner_tar_zst(members: dict) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return subprocess.run(
        ["zstd", "-q", "-c"], input=buf.getvalue(), stdout=subprocess.PIPE, check=True
    ).stdout


def _write_outer(path: Path, entries: dict) -> None:
    """`entries` maps an inner archive name to its already-zstd-compressed bytes
    (as produced by `_inner_tar_zst`). Assembles and zstd-compresses the outer
    tar at `path`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.with_name(path.name + ".plain")
    with tarfile.open(plain, "w") as tar:
        for name, blob in entries.items():
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    _zstd(plain, path)
    plain.unlink()


def _fake_load_structure(path):
    """Stands in for binder_map.chains.load_structure: returns the raw blob
    bytes untouched, so downstream fakes can branch on its content."""
    return path.read_bytes()


def _fake_split_chains(blob):
    """Stands in for binder_map.chains.split_chains. A blob starting with
    b"BADLEN" produces a too-short chain (triggers the n_res != 80 gate);
    anything else produces exactly 80 residues."""
    n = 5 if blob.startswith(b"BADLEN") else 80
    protein = SimpleNamespace(chain_id=["A"] * n, res_id=list(range(n)), _blob=blob)
    rna = SimpleNamespace()
    return protein, rna


def _fake_coord_hash(protein):
    return f"hash:{protein._blob!r}"


def _fake_fingerprint(protein, rna):
    n = len(s1.FINGERPRINT_COLUMNS)
    # main() calls fp.vector.tolist(), so this must be a real numpy array,
    # not a plain list.
    return SimpleNamespace(vector=np.zeros(n, dtype=np.int32), per_nt=[0] * 16,
                            offframe=0, total=0)


def _fake_approach(protein, rna):
    return SimpleNamespace(centroid=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0))


def _fake_write_backbone_pdb(protein, design_id):
    # Non-repetitive payload so the resulting shard_NN.tar.gz stays well
    # above the empty-tarball FATAL threshold even when gzip would otherwise
    # crush a tiny, highly-repetitive fake PDB down to a few dozen bytes.
    import os
    return os.urandom(512) + design_id.encode()


def _patch_pipeline(monkeypatch):
    monkeypatch.setattr(s1, "load_structure", _fake_load_structure)
    monkeypatch.setattr(s1, "split_chains", _fake_split_chains)
    monkeypatch.setattr(s1, "coord_hash", _fake_coord_hash)
    monkeypatch.setattr(s1, "fingerprint", _fake_fingerprint)
    monkeypatch.setattr(s1, "approach", _fake_approach)
    monkeypatch.setattr(s1, "write_backbone_pdb", _fake_write_backbone_pdb)


def _run_main(tmp_path, monkeypatch, manifest: pd.DataFrame, *, shard=0, n_shards=1):
    """Write `manifest` and args, invoke main(), return (hashes, fp, missing)
    DataFrames read back from disk. Assumes `s1.RB` has already been
    monkeypatched to `tmp_path / "rb"` and the referenced archives built
    there by the caller.
    """
    build = tmp_path / "build"
    (build / "manifests").mkdir(parents=True)
    manifest.to_csv(build / "manifests" / "design.tsv", sep="\t", index=False)
    results = tmp_path / "results"

    monkeypatch.setattr(sys, "argv", [
        "stage1_extract_nested.py",
        "--shard", str(shard), "--n-shards", str(n_shards),
        "--build", str(build), "--results", str(results),
    ])
    rc = s1.main()
    assert rc == 0

    out = results / "design"
    hashes = pd.read_csv(out / "hashes" / f"shard_{shard:02d}.tsv", sep="\t")
    fp = pd.read_parquet(out / "fp" / f"shard_{shard:02d}.parquet")
    missing = pd.read_csv(out / "missing" / f"shard_{shard:02d}.tsv", sep="\t")
    return hashes, fp, missing, out / "shards" / f"shard_{shard:02d}.tar.gz"


def test_by_member_keeps_rows_distinct_across_colliding_outer_archives(tmp_path, monkeypatch):
    """Two manifest rows share the same (inner_path, member_path) but live in
    two different outer archives - exactly the collision shape the reviewer
    traced: keying by (inner_path, member_path) alone collapses d1 and d2
    onto one dict entry and one of the two designs silently vanishes."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)

    shared_inner = "m/designs/u0.tar.zst"
    shared_member = "shared.cif"
    _write_outer(rb / "a.tar.zst",
                 {shared_inner: _inner_tar_zst({shared_member: b"GOOD:d1"})})
    _write_outer(rb / "b.tar.zst",
                 {shared_inner: _inner_tar_zst({shared_member: b"GOOD:d2"})})

    manifest = pd.DataFrame({
        "design_id": ["d1", "d2"],
        "outer_path": ["a.tar.zst", "b.tar.zst"],
        "inner_path": [shared_inner, shared_inner],
        "member_path": [shared_member, shared_member],
    })

    hashes, fp, missing, _tar = _run_main(tmp_path, monkeypatch, manifest)

    assert sorted(hashes["design_id"]) == ["d1", "d2"]
    assert sorted(fp["design_id"]) == ["d1", "d2"]
    assert missing.empty


def test_by_member_collision_assert_fires_on_a_genuine_duplicate_row(tmp_path, monkeypatch):
    """Sanity check on the guard itself: a manifest with a literal duplicate
    row (same outer/inner/member AND design_id twice) is a real key
    collision the assert must catch, not something Fix 2 silently papers
    over.

    Also covers Minor 4: the assert now runs before any scratch directory is
    created, so a real collision leaves nothing under TMPDIR to clean up."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst",
                 {"m/designs/u0.tar.zst": _inner_tar_zst({"x.cif": b"GOOD:d0"})})

    manifest = pd.DataFrame({
        "design_id": ["d0", "d0"],
        "outer_path": ["a.tar.zst", "a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u0.tar.zst"],
        "member_path": ["x.cif", "x.cif"],
    })
    build = tmp_path / "build"
    (build / "manifests").mkdir(parents=True)
    manifest.to_csv(build / "manifests" / "design.tsv", sep="\t", index=False)
    tmpdir_root = tmp_path / "tmpdir_root"
    monkeypatch.setenv("TMPDIR", str(tmpdir_root))
    monkeypatch.setattr(sys, "argv", [
        "stage1_extract_nested.py", "--shard", "0", "--n-shards", "1",
        "--build", str(build), "--results", str(tmp_path / "results"),
    ])
    with pytest.raises(AssertionError, match="by_member key collision"):
        s1.main()

    assert not (tmpdir_root / "bm_s1n_0").exists(), (
        "the collision assert must fire before the scratch dir is created")


def test_duplicate_member_name_within_one_inner_archive_is_counted_once(tmp_path, monkeypatch):
    """Minor 5: if iter_members ever yields the same (outer, inner, member)
    key twice - e.g. an inner tar carrying two members with the same name -
    the second delivery must be a no-op, not a second accounting hit. Before
    the fix the row landed in hashes/fp twice, and on a fuller manifest that
    trips the len(hashes)+len(missing)==len(mine) accounting assert."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)

    # tarfile permits two members with the same name in one archive;
    # iter_members's streaming scan yields both, since it never deduplicates.
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for data in (b"GOOD:first", b"GOOD:second"):
            info = tarfile.TarInfo("dup.cif")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    inner_zst = subprocess.run(["zstd", "-q", "-c"], input=buf.getvalue(),
                                stdout=subprocess.PIPE, check=True).stdout
    _write_outer(rb / "a.tar.zst", {"m/designs/u0.tar.zst": inner_zst})

    manifest = pd.DataFrame({
        "design_id": ["d0"], "outer_path": ["a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst"], "member_path": ["dup.cif"],
    })
    hashes, fp, missing, _tar = _run_main(tmp_path, monkeypatch, manifest)

    assert list(hashes["design_id"]) == ["d0"]
    assert list(fp["design_id"]) == ["d0"]
    assert missing.empty


def test_tarball_is_promoted_last_so_a_crash_mid_write_leaves_no_tarball(tmp_path, monkeypatch):
    """Minor 1: the tarball must be the LAST thing written, after the three
    metadata files, so a crash between them never leaves a promoted
    shard_NN.tar.gz with no hashes/fp/missing companion - downstream stages
    glob on the tarball as the completion marker."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst",
                 {"m/designs/u0.tar.zst": _inner_tar_zst({"x.cif": b"GOOD:d0"})})
    manifest = pd.DataFrame({
        "design_id": ["d0"], "outer_path": ["a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst"], "member_path": ["x.cif"],
    })
    build = tmp_path / "build"
    (build / "manifests").mkdir(parents=True)
    manifest.to_csv(build / "manifests" / "design.tsv", sep="\t", index=False)
    results = tmp_path / "results"
    monkeypatch.setattr(sys, "argv", [
        "stage1_extract_nested.py", "--shard", "0", "--n-shards", "1",
        "--build", str(build), "--results", str(results),
    ])

    def boom(self, *a, **k):
        raise RuntimeError("simulated crash mid-write")

    # The first metadata write (hashes, a .to_csv call) is made to blow up.
    # If the tarball were still promoted first (the pre-fix order), the
    # promoted file would already exist by the time this fires; the fix
    # requires it to still be sitting under `.partial`.
    monkeypatch.setattr(pd.DataFrame, "to_csv", boom)

    with pytest.raises(RuntimeError, match="simulated crash mid-write"):
        s1.main()

    out = results / "design"
    tar_path = out / "shards" / "shard_00.tar.gz"
    partial = tar_path.with_name(tar_path.name + ".partial")
    assert not tar_path.exists(), "tarball must not be promoted before metadata writes finish"
    assert partial.exists(), "the not-yet-promoted tarball should still be on disk"


def test_corrupt_inner_archive_yields_missing_rows_not_a_crash(tmp_path, monkeypatch):
    """A corrupt inner archive raises CalledProcessError out of iter_members
    while main()'s `for` loop advances it. That must not kill the shard: the
    corrupt archive's rows are recorded as missing and the run moves on to
    the next (good) archive."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)

    _write_outer(rb / "good.tar.zst",
                 {"m/designs/u0.tar.zst": _inner_tar_zst({"ok.cif": b"GOOD:d0"})})
    # The inner "archive" is not a valid zstd stream at all - iter_members's
    # _decompress call raises subprocess.CalledProcessError trying to read it.
    _write_outer(rb / "bad.tar.zst", {"m/designs/u1.tar.zst": b"not a zstd stream"})

    manifest = pd.DataFrame({
        "design_id": ["d0", "d1", "d2"],
        "outer_path": ["good.tar.zst", "bad.tar.zst", "bad.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u1.tar.zst", "m/designs/u1.tar.zst"],
        "member_path": ["ok.cif", "x.cif", "y.cif"],
    })

    hashes, fp, missing, tar_path = _run_main(tmp_path, monkeypatch, manifest)

    assert sorted(hashes["design_id"]) == ["d0"]
    bad = missing.set_index("design_id")["reason"]
    assert bad["d1"] == "archive_read_error:CalledProcessError"
    assert bad["d2"] == "archive_read_error:CalledProcessError"
    assert tar_path.exists() and tar_path.stat().st_size >= 200
    assert not tar_path.with_name(tar_path.name + ".partial").exists()


def test_empty_shard_with_no_outer_archives_assigned_exits_cleanly(tmp_path, monkeypatch):
    """--n-shards greater than the number of outer archives leaves some
    shards with zero rows. That must exit 0 with empty-but-schema'd outputs,
    not a spurious FATAL from the tarball-size check (which used to fire
    unconditionally, unlike the zero-parsed check next to it)."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst",
                 {"m/designs/u0.tar.zst": _inner_tar_zst({"x.cif": b"GOOD:d0"})})

    manifest = pd.DataFrame({
        "design_id": ["d0"], "outer_path": ["a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst"], "member_path": ["x.cif"],
    })
    # shard 1 of 2, with only 1 outer archive total: shard 1 gets nothing.
    hashes, fp, missing, tar_path = _run_main(tmp_path, monkeypatch, manifest,
                                               shard=1, n_shards=2)
    assert hashes.empty and missing.empty
    assert list(fp.columns) == s1.FP_COLUMNS
    assert tar_path.exists()


def test_blank_member_is_member_empty_distinct_from_member_missing(tmp_path, monkeypatch):
    """A member present in the archive but whitespace-only must be tellable
    apart from a member never delivered at all - this project has a history
    of blank-but-right-length data passing undetected. A third, successful
    design (d2) rides along so the shard doesn't hit the separate
    FATAL-zero-parsed check, which is not what this test is about."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst", {
        "m/designs/u0.tar.zst": _inner_tar_zst({"blank.cif": b"   \n",
                                                 "ok.cif": b"GOOD:d2"}),
    })

    manifest = pd.DataFrame({
        "design_id": ["d0", "d1", "d2"],
        "outer_path": ["a.tar.zst", "a.tar.zst", "a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u0.tar.zst", "m/designs/u0.tar.zst"],
        "member_path": ["blank.cif", "absent.cif", "ok.cif"],
    })
    hashes, _, missing, _ = _run_main(tmp_path, monkeypatch, manifest)
    assert list(hashes["design_id"]) == ["d2"]
    reasons = missing.set_index("design_id")["reason"]
    assert reasons["d0"] == "member_empty"
    assert reasons["d1"] == "member_missing"


def test_wrong_length_reason_names_the_length_gate_not_chain_selection(tmp_path, monkeypatch):
    """A second, successful design (d1) rides along so the shard doesn't hit
    the separate FATAL-zero-parsed check, which is not what this test is
    about."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst", {
        "m/designs/u0.tar.zst": _inner_tar_zst({"short.cif": b"BADLEN:d0",
                                                 "ok.cif": b"GOOD:d1"}),
    })
    manifest = pd.DataFrame({
        "design_id": ["d0", "d1"],
        "outer_path": ["a.tar.zst", "a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u0.tar.zst"],
        "member_path": ["short.cif", "ok.cif"],
    })
    hashes, _, missing, _ = _run_main(tmp_path, monkeypatch, manifest)
    assert list(hashes["design_id"]) == ["d1"]
    assert missing.set_index("design_id")["reason"]["d0"] == "wrong_length:n_res=5"


def test_fatal_tarball_empty_message_names_the_partial_path(tmp_path, monkeypatch):
    """Minor 3: the FATAL message for an undersized tarball must name the
    `.partial` file - the one that actually exists on disk at that point -
    not the post-promotion name, which does not exist yet."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    _write_outer(rb / "a.tar.zst",
                 {"m/designs/u0.tar.zst": _inner_tar_zst({"x.cif": b"GOOD:d0"})})
    manifest = pd.DataFrame({
        "design_id": ["d0"], "outer_path": ["a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst"], "member_path": ["x.cif"],
    })
    build = tmp_path / "build"
    (build / "manifests").mkdir(parents=True)
    manifest.to_csv(build / "manifests" / "design.tsv", sep="\t", index=False)
    results = tmp_path / "results"
    monkeypatch.setattr(sys, "argv", [
        "stage1_extract_nested.py", "--shard", "0", "--n-shards", "1",
        "--build", str(build), "--results", str(results),
    ])

    # Force the "tarball too small" branch regardless of what the fake
    # pipeline actually wrote, by reporting a near-zero size for exactly the
    # .partial file main() stats - nothing else is touched.
    real_stat = Path.stat

    def stat_patch(self, *a, **k):
        if self.name.endswith(".tar.gz.partial"):
            return SimpleNamespace(st_size=0)
        return real_stat(self, *a, **k)

    monkeypatch.setattr(Path, "stat", stat_patch)

    with pytest.raises(SystemExit, match=r"FATAL: shard tarball .*\.partial is empty"):
        s1.main()


def test_output_files_are_not_written_when_validation_fails(tmp_path, monkeypatch):
    """Fix 4: a shard that fails the FATAL zero-parsed check must not leave
    hashes/fp/missing files behind, and the shard tarball must stay under
    its `.partial` name rather than being promoted to the final name."""
    _patch_pipeline(monkeypatch)
    rb = tmp_path / "rb"
    monkeypatch.setattr(s1, "RB", rb)
    # The one referenced archive is entirely absent -> 0 parsed, 1 missing,
    # len(mine)=1 and not hashes -> FATAL zero-parsed SystemExit.
    manifest = pd.DataFrame({
        "design_id": ["d0"], "outer_path": ["a.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst"], "member_path": ["x.cif"],
    })
    build = tmp_path / "build"
    (build / "manifests").mkdir(parents=True)
    manifest.to_csv(build / "manifests" / "design.tsv", sep="\t", index=False)
    results = tmp_path / "results"
    monkeypatch.setattr(sys, "argv", [
        "stage1_extract_nested.py", "--shard", "0", "--n-shards", "1",
        "--build", str(build), "--results", str(results),
    ])
    with pytest.raises(SystemExit, match="FATAL"):
        s1.main()

    out = results / "design"
    assert not (out / "hashes" / "shard_00.tsv").exists()
    assert not (out / "fp" / "shard_00.parquet").exists()
    assert not (out / "missing" / "shard_00.tsv").exists()
    assert not (out / "shards" / "shard_00.tar.gz").exists()
