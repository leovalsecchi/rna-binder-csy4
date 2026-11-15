"""Tests for the stage-2 -> stage-3 gate (scripts/gate_checks.py, probe_project.py).

The gate's entire value is whether it FAILS when it should, so every test here
is written to break if the check it covers is inverted or deleted. Fixtures are
real files -- small TSVs and real gzipped tars built with `tarfile` -- so the
readers are exercised, not mocked.
"""
import io
import math
import tarfile

import pytest

from gate_checks import (
    MIN_SHARED_HASHES,
    NO_SHARDS,
    check_hash_identity,
    check_id_format,
    check_split_nonempty,
)
from gate_checks import main as gate_main
from probe_project import main as probe_main
from probe_project import project_cpu_h, valid_timing


# --------------------------------------------------------------------------- #
# fixture builders
# --------------------------------------------------------------------------- #
def _write_hash_shard(tree, name, rows):
    """rows: list of (design_id, coord_hash)."""
    out = tree / "design" / "hashes"
    out.mkdir(parents=True, exist_ok=True)
    lines = ["design_id\tcoord_hash\tprotein_chain\tn_res"]
    lines += [f"{d}\t{h}\tA\t80" for d, h in rows]
    (out / name).write_text("\n".join(lines) + "\n")


def _write_tar(path, members):
    """members: dict name -> bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            tar.addfile(info, io.BytesIO(data))


def _write_shard_tree(tree, members, name="shard_00.tar.gz"):
    _write_tar(tree / "design" / "shards" / name, members)


def _write_manifest(path, ids):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["design_id\tstack"] + [f"{i}\trfd3" for i in ids]
    path.write_text("\n".join(lines) + "\n")


def _ids(n, prefix="design_"):
    return [f"{prefix}{i:06d}" for i in range(n)]


def _rows(ids, hash_of=lambda i: f"hash_{i}"):
    return [(i, hash_of(i)) for i in ids]


N_OK = MIN_SHARED_HASHES + 1000  # comfortably above the floor


def _valid_layout(tmp_path):
    """A layout that passes all three checks; individual tests then break one."""
    ids = _ids(N_OK)
    _write_hash_shard(tmp_path / "results_shared", "shard_00.tsv", _rows(ids))
    _write_hash_shard(tmp_path / "results_v6", "shard_00.tsv", _rows(ids))
    members = {f"{i}.pdb": b"ATOM\n" for i in ids[:10]}
    _write_shard_tree(tmp_path / "results_shared", members)
    _write_shard_tree(tmp_path / "results_full", members)
    _write_shard_tree(tmp_path / "results_lib12k", members)
    _write_manifest(tmp_path / "build_shared" / "manifests" / "design.tsv", ids)
    return ids


def _argv(tmp_path):
    return [
        "--shared", str(tmp_path / "results_shared"),
        "--ref", str(tmp_path / "results_v6"),
        "--manifest", str(tmp_path / "build_shared" / "manifests" / "design.tsv"),
        "--tree", str(tmp_path / "results_full"),
        "--tree", str(tmp_path / "results_lib12k"),
    ]


# --------------------------------------------------------------------------- #
# check 1: hash identity against the pre-repack reference tree
# --------------------------------------------------------------------------- #
def test_hash_identity_reports_zero_mismatches_when_trees_agree(tmp_path):
    ids = _ids(N_OK)
    _write_hash_shard(tmp_path / "shared", "shard_00.tsv", _rows(ids))
    _write_hash_shard(tmp_path / "ref", "shard_00.tsv", _rows(ids))
    n_shared, n_mismatch = check_hash_identity(tmp_path / "shared", tmp_path / "ref")
    assert n_shared == N_OK
    assert n_mismatch == 0


def test_hash_identity_counts_a_single_flipped_hash(tmp_path):
    ids = _ids(N_OK)
    _write_hash_shard(tmp_path / "shared", "shard_00.tsv", _rows(ids))
    ref_rows = _rows(ids)
    ref_rows[7] = (ref_rows[7][0], "0000000000000000deadbeefdeadbeef")
    _write_hash_shard(tmp_path / "ref", "shard_00.tsv", ref_rows)
    n_shared, n_mismatch = check_hash_identity(tmp_path / "shared", tmp_path / "ref")
    assert n_shared == N_OK
    assert n_mismatch == 1


def test_hash_identity_intersects_rather_than_unions(tmp_path):
    """Ids present in only one tree are not compared and do not inflate n_shared."""
    _write_hash_shard(tmp_path / "shared", "shard_00.tsv",
                      _rows(["a", "b", "c"]))
    _write_hash_shard(tmp_path / "ref", "shard_00.tsv",
                      _rows(["b", "c", "d", "e"]))
    n_shared, n_mismatch = check_hash_identity(tmp_path / "shared", tmp_path / "ref")
    assert n_shared == 2
    assert n_mismatch == 0


def test_hash_identity_reads_every_shard(tmp_path):
    _write_hash_shard(tmp_path / "shared", "shard_00.tsv", _rows(["a", "b"]))
    _write_hash_shard(tmp_path / "shared", "shard_01.tsv", _rows(["c"]))
    _write_hash_shard(tmp_path / "ref", "shard_00.tsv", _rows(["a", "b", "c"]))
    assert check_hash_identity(tmp_path / "shared", tmp_path / "ref") == (3, 0)


def test_hash_identity_on_missing_tree_is_empty_not_an_error(tmp_path):
    _write_hash_shard(tmp_path / "shared", "shard_00.tsv", _rows(["a"]))
    assert check_hash_identity(tmp_path / "shared", tmp_path / "nope") == (0, 0)


def test_gate_passes_on_a_clean_layout(tmp_path, capsys):
    _valid_layout(tmp_path)
    assert gate_main(_argv(tmp_path)) == 0
    assert "GATE PASS" in capsys.readouterr().out


def test_gate_fails_on_one_flipped_coord_hash(tmp_path, capsys):
    ids = _valid_layout(tmp_path)
    ref_rows = _rows(ids)
    ref_rows[3] = (ref_rows[3][0], "ffffffffffffffffffffffffffffffff")
    _write_hash_shard(tmp_path / "results_v6", "shard_00.tsv", ref_rows)
    assert gate_main(_argv(tmp_path)) == 1
    out = capsys.readouterr().out
    assert "coord_hash mismatches" in out
    assert "GATE FAIL" in out


def test_gate_fails_when_too_few_ids_are_shared_even_if_all_match(tmp_path, capsys):
    """The floor exists so an almost-empty intersection cannot pass as agreement."""
    _valid_layout(tmp_path)
    few = _ids(MIN_SHARED_HASHES - 1)
    _write_hash_shard(tmp_path / "results_v6", "shard_00.tsv", _rows(few))
    n_shared, n_mismatch = check_hash_identity(tmp_path / "results_shared",
                                               tmp_path / "results_v6")
    assert n_mismatch == 0, "every shared hash still agrees"
    assert n_shared == MIN_SHARED_HASHES - 1
    assert gate_main(_argv(tmp_path)) == 1
    assert "too few" in capsys.readouterr().out


def test_gate_passes_exactly_at_the_shared_id_floor(tmp_path):
    """Boundary: the floor is `< MIN`, not `<= MIN`."""
    _valid_layout(tmp_path)
    exact = _ids(MIN_SHARED_HASHES)
    _write_hash_shard(tmp_path / "results_v6", "shard_00.tsv", _rows(exact))
    assert check_hash_identity(tmp_path / "results_shared",
                               tmp_path / "results_v6")[0] == MIN_SHARED_HASHES
    assert gate_main(_argv(tmp_path)) == 0


# --------------------------------------------------------------------------- #
# check 2: manifest design_ids vs shard .pdb member stems
# --------------------------------------------------------------------------- #
def test_id_format_all_stems_are_manifest_ids(tmp_path):
    _write_manifest(tmp_path / "design.tsv", ["a", "b", "c"])
    _write_tar(tmp_path / "shard_00.tar.gz",
               {"a.pdb": b"X\n", "b.pdb": b"X\n", "c.pdb": b"X\n"})
    assert check_id_format(tmp_path / "design.tsv",
                           tmp_path / "shard_00.tar.gz") == (3, 3)


def test_id_format_catches_a_single_trailing_character(tmp_path):
    """The silent failure: one extra character and the split writes empty tars."""
    _write_manifest(tmp_path / "design.tsv", ["a", "b", "c"])
    _write_tar(tmp_path / "shard_00.tar.gz",
               {"a.pdb": b"X\n", "b_.pdb": b"X\n", "c.pdb": b"X\n"})
    n_members, n_hits = check_id_format(tmp_path / "design.tsv",
                                        tmp_path / "shard_00.tar.gz")
    assert n_members == 3
    assert n_hits == 2


def test_id_format_ignores_non_pdb_members(tmp_path):
    _write_manifest(tmp_path / "design.tsv", ["a"])
    _write_tar(tmp_path / "shard_00.tar.gz", {"a.pdb": b"X\n", "notes.txt": b"X\n"})
    assert check_id_format(tmp_path / "design.tsv",
                           tmp_path / "shard_00.tar.gz") == (1, 1)


def test_id_format_reports_zero_members_for_a_tar_without_pdbs(tmp_path):
    _write_manifest(tmp_path / "design.tsv", ["a"])
    _write_tar(tmp_path / "shard_00.tar.gz", {"notes.txt": b"X\n"})
    assert check_id_format(tmp_path / "design.tsv",
                           tmp_path / "shard_00.tar.gz") == (0, 0)


def test_gate_fails_when_a_member_stem_is_not_a_manifest_id(tmp_path, capsys):
    ids = _valid_layout(tmp_path)
    members = {f"{i}.pdb": b"ATOM\n" for i in ids[:9]}
    members[f"{ids[9]}X.pdb"] = b"ATOM\n"
    _write_shard_tree(tmp_path / "results_shared", members)
    assert gate_main(_argv(tmp_path)) == 1
    assert "not manifest" in capsys.readouterr().out


def test_gate_fails_when_the_shared_shard_has_no_pdb_members(tmp_path, capsys):
    _valid_layout(tmp_path)
    _write_shard_tree(tmp_path / "results_shared", {"readme.txt": b"X\n"})
    assert gate_main(_argv(tmp_path)) == 1
    assert "no .pdb members" in capsys.readouterr().out


def test_gate_fails_when_the_shared_tree_has_no_shards_to_check(tmp_path, capsys):
    _valid_layout(tmp_path)
    (tmp_path / "results_shared" / "design" / "shards" / "shard_00.tar.gz").unlink()
    assert gate_main(_argv(tmp_path)) == 1
    assert "nothing to check the manifest against" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# check 3: the split trees hold structures
# --------------------------------------------------------------------------- #
def test_split_nonempty_counts_members_of_the_first_shard(tmp_path):
    _write_shard_tree(tmp_path / "tree", {"a.pdb": b"X\n", "b.pdb": b"X\n"})
    assert check_split_nonempty(tmp_path / "tree") == 2


def test_split_nonempty_returns_zero_for_an_empty_tar(tmp_path):
    _write_shard_tree(tmp_path / "tree", {})
    assert check_split_nonempty(tmp_path / "tree") == 0


def test_split_nonempty_returns_no_shards_sentinel_when_there_are_none(tmp_path):
    (tmp_path / "tree" / "design" / "shards").mkdir(parents=True)
    assert check_split_nonempty(tmp_path / "tree") == NO_SHARDS
    assert NO_SHARDS != 0, "the sentinel must be distinguishable from an empty shard"


def test_split_nonempty_uses_the_first_shard_by_sorted_name(tmp_path):
    _write_shard_tree(tmp_path / "tree", {"a.pdb": b"X\n"}, name="shard_00.tar.gz")
    _write_shard_tree(tmp_path / "tree", {"b.pdb": b"X\n", "c.pdb": b"X\n"},
                      name="shard_01.tar.gz")
    assert check_split_nonempty(tmp_path / "tree") == 1


def test_gate_fails_when_a_split_tree_shard_is_empty(tmp_path, capsys):
    _valid_layout(tmp_path)
    _write_shard_tree(tmp_path / "results_lib12k", {})
    assert gate_main(_argv(tmp_path)) == 1
    out = capsys.readouterr().out
    assert "first shard is empty" in out
    assert "results_lib12k" in out


def test_gate_fails_when_a_split_tree_has_no_shards(tmp_path, capsys):
    _valid_layout(tmp_path)
    (tmp_path / "results_full" / "design" / "shards" / "shard_00.tar.gz").unlink()
    assert gate_main(_argv(tmp_path)) == 1
    assert "tree has no shards" in capsys.readouterr().out


def test_gate_fails_when_no_tree_was_given(tmp_path, capsys):
    """A gate that checked nothing must not report success."""
    _valid_layout(tmp_path)
    argv = _argv(tmp_path)
    argv = argv[:argv.index("--tree")]
    assert gate_main(argv) == 1
    assert "checked no trees" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# probe projection
# --------------------------------------------------------------------------- #
def test_project_exponent_is_log_ratio_over_log_three():
    p, _ = project_cpu_h(100.0, 300.0)
    assert p == pytest.approx(1.0)
    p9, _ = project_cpu_h(100.0, 900.0)
    assert p9 == pytest.approx(2.0)


def test_project_cpu_h_matches_the_plan_formula():
    p, cpu_h = project_cpu_h(100.0, 300.0)
    expected = 100.0 * (500_000 / 50_000) ** p * (500_000 / 2000) * 16 / 3600
    assert cpu_h == pytest.approx(expected)
    assert cpu_h == pytest.approx(1111.11, rel=1e-4)


def test_project_cpu_h_grows_with_the_exponent():
    _, cheap = project_cpu_h(100.0, 300.0)
    _, dear = project_cpu_h(100.0, 900.0)
    assert dear > cheap


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_project_cpu_h_refuses_garbage_timings(bad):
    with pytest.raises(ValueError):
        project_cpu_h(bad, 300.0)
    with pytest.raises(ValueError):
        project_cpu_h(100.0, bad)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_valid_timing_rejects_garbage(bad):
    assert valid_timing(bad) is False


def test_valid_timing_accepts_a_positive_finite_number():
    assert valid_timing(189.0) is True


def test_probe_main_exits_zero_under_the_ceiling(capsys):
    assert probe_main(["--t50", "100", "--t150", "300", "--ceiling", "6000"]) == 0
    assert "PASS" in capsys.readouterr().out


def test_probe_main_exits_nonzero_over_the_ceiling(capsys):
    assert probe_main(["--t50", "100", "--t150", "900", "--ceiling", "6000"]) == 1
    assert "exceeds" in capsys.readouterr().out


def test_probe_main_ceiling_comparison_is_discriminating():
    """Same projection (~1111 CPU-h): passes under a looser ceiling, fails a tighter."""
    assert probe_main(["--t50", "100", "--t150", "300", "--ceiling", "1200"]) == 0
    assert probe_main(["--t50", "100", "--t150", "300", "--ceiling", "1000"]) == 1


def test_probe_main_default_ceiling_is_6000(capsys):
    assert probe_main(["--t50", "100", "--t150", "300"]) == 0
    assert "6,000" in capsys.readouterr().out


@pytest.mark.parametrize("bad", ["0", "-5", "nan", "inf"])
def test_probe_main_fails_closed_on_a_garbage_timing(bad, capsys):
    """A probe arm that did not run must stop the chain, not release it."""
    assert probe_main(["--t50", bad, "--t150", "300"]) == 1
    assert probe_main(["--t50", "100", "--t150", bad]) == 1


def test_probe_main_rejects_a_non_numeric_timing():
    """An empty or text timing dies at argparse -- still non-zero, never 0."""
    with pytest.raises(SystemExit) as exc:
        probe_main(["--t50", "", "--t150", "300"])
    assert exc.value.code != 0


def test_probe_projection_is_finite_for_realistic_inputs():
    _, cpu_h = project_cpu_h(189.0, 570.0)
    assert math.isfinite(cpu_h)
    assert cpu_h > 0


# --------------------------------------------------------------------------- #
# the shared-id floor, and the degenerate ways check 1 could pass trivially
# --------------------------------------------------------------------------- #
def test_the_shared_id_floor_is_three_thousand():
    """Pinned, because the margin above it is measured and thin.

    Measured on Euler 2026-09-24: build_shared/manifests/design.tsv holds
    508,207 design_ids, results_v6 holds 55,706, and the two intersect in 6,377.
    results_shared's ids are a subset of the manifest, so 6,377 is a hard upper
    bound on n_shared -- the old 5000 floor left 27% of headroom and would have
    hard-blocked a good run after ~22% extraction losses, *after* full_s2 had
    already burned. 3000 still compares thousands of independent hashes, and a
    single mismatch fails the check at any n_shared.
    """
    assert MIN_SHARED_HASHES == 3000


def test_gate_refuses_to_compare_a_tree_with_itself(tmp_path, capsys):
    """A tree always agrees with itself; that must not be a way to pass."""
    _valid_layout(tmp_path)
    argv = _argv(tmp_path)
    argv[argv.index("--ref") + 1] = str(tmp_path / "results_shared")
    assert gate_main(argv) == 1
    out = capsys.readouterr().out
    assert "not evidence of anything" in out
    assert "GATE FAIL" in out


def test_gate_refuses_a_symlink_to_the_tree_it_is_comparing(tmp_path, capsys):
    """The refusal is on the resolved path, not on the spelling."""
    _valid_layout(tmp_path)
    link = tmp_path / "ref_link"
    link.symlink_to(tmp_path / "results_shared")
    argv = _argv(tmp_path)
    argv[argv.index("--ref") + 1] = str(link)
    assert gate_main(argv) == 1
    assert "GATE FAIL" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# the ceiling is the brake, so an unusable ceiling must fail, not widen
# --------------------------------------------------------------------------- #
# ("-inf" is absent on purpose: argparse reads it as a flag, not a value, and
# dies with exit 2 -- non-zero either way, but not through this branch.)
@pytest.mark.parametrize("bad", ["nan", "inf", "0", "-5"])
def test_probe_main_refuses_a_ceiling_that_disables_the_brake(bad, capsys):
    """`GATE_CPUH_CEILING=nan` is forwarded by the sbatch and propagated by
    `sbatch --export=ALL`; with it every `cpu_h > ceiling` is False."""
    assert probe_main(["--t50", "100", "--t150", "300", "--ceiling", bad]) == 1
    assert "ceiling" in capsys.readouterr().out


def test_probe_main_still_accepts_a_finite_positive_ceiling(capsys):
    assert probe_main(["--t50", "100", "--t150", "300", "--ceiling", "1200"]) == 0
    assert "PASS" in capsys.readouterr().out
