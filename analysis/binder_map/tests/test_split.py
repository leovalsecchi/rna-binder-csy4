import io
import sys
import tarfile
from pathlib import Path

import pandas as pd

from stage1b_split import (DUPLICATE_REASON, main, manifest_ids, merge_missing,
                           reconcile, should_warn_census, split_hashes,
                           split_tables, split_tar)


def test_split_tables_keeps_only_requested_designs():
    frame = pd.DataFrame({"design_id": ["a", "b", "c"], "x": [1, 2, 3]})
    out = split_tables(frame, {"a", "c"})
    assert list(out["design_id"]) == ["a", "c"]


def test_split_tables_on_empty_keep_returns_empty_with_columns():
    frame = pd.DataFrame({"design_id": ["a"], "x": [1]})
    out = split_tables(frame, set())
    assert len(out) == 0
    assert list(out.columns) == ["design_id", "x"]


def test_split_tables_ignores_ids_that_are_not_present():
    frame = pd.DataFrame({"design_id": ["a"], "x": [1]})
    assert list(split_tables(frame, {"a", "ghost"})["design_id"]) == ["a"]


def _tar(path, members):
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            tar.addfile(info, io.BytesIO(data))


def test_split_tar_copies_only_matching_members(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"AAA\n", "b.pdb": b"BBB\n", "c.pdb": b"CCC\n"})
    dst = tmp_path / "out.tar.gz"
    n, stems = split_tar(src, dst, {"a", "c"})
    assert n == 2
    assert stems == {"a", "c"}
    with tarfile.open(dst, "r:gz") as tar:
        assert sorted(m.name for m in tar) == ["a.pdb", "c.pdb"]


def test_split_tar_returns_zero_when_nothing_matches(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"AAA\n"})
    assert split_tar(src, tmp_path / "out.tar.gz", {"zz"}) == (0, set())


def test_split_tar_stems_omit_a_wanted_member_the_shard_never_had(tmp_path):
    """The stems are what is on disk, not what was asked for: `b` is in `keep`
    but not in the tar, so it must not come back as written."""
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"AAA\n"})
    n, stems = split_tar(src, tmp_path / "out.tar.gz", {"a", "b"})
    assert n == 1
    assert stems == {"a"}


def test_manifest_ids_reads_design_id_column(tmp_path):
    build_dir = tmp_path / "build"
    (build_dir / "manifests").mkdir(parents=True)
    manifest = build_dir / "manifests" / "design.tsv"
    manifest.write_text("design_id\tother_col\na\tval_a\nb\tval_b\nc\tval_c\n")
    result = manifest_ids(build_dir)
    assert result == {"a", "b", "c"}


def test_manifest_ids_handles_single_entry(tmp_path):
    build_dir = tmp_path / "build"
    (build_dir / "manifests").mkdir(parents=True)
    manifest = build_dir / "manifests" / "design.tsv"
    manifest.write_text("design_id\nonly_one\n")
    result = manifest_ids(build_dir)
    assert result == {"only_one"}


def _hashes(pairs):
    """A stage 1 hashes table: [(design_id, coord_hash), ...]."""
    return pd.DataFrame({
        "design_id": [d for d, _ in pairs],
        "coord_hash": [h for _, h in pairs],
        "protein_chain": ["A"] * len(pairs),
        "n_res": [80] * len(pairs),
    })


def test_split_hashes_keeps_rows_whose_structure_exists():
    frame = _hashes([("a", "h1"), ("b", "h2")])
    kept, dropped = split_hashes(frame, {"a", "b"})
    assert list(kept["design_id"]) == ["a", "b"]
    assert len(dropped) == 0


def test_split_hashes_drops_the_row_with_no_structure():
    """Two designs share a coord_hash; stage 1 wrote only one `.pdb`."""
    frame = _hashes([("ghost", "same"), ("keeper", "same")])
    kept, dropped = split_hashes(frame, {"keeper"})
    assert list(kept["design_id"]) == ["keeper"]
    assert list(dropped["design_id"]) == ["ghost"]
    assert list(dropped["reason"]) == [DUPLICATE_REASON]
    assert list(dropped["coord_hash"]) == ["same"]


def test_split_hashes_dropped_rows_carry_the_missing_schema():
    frame = _hashes([("ghost", "same")])
    _kept, dropped = split_hashes(frame, set())
    assert list(dropped.columns) == ["design_id", "reason", "coord_hash"]


def test_split_hashes_keeps_every_row_when_nothing_was_deduped():
    frame = _hashes([("a", "h1")])
    kept, dropped = split_hashes(frame, {"a"})
    assert len(kept) == 1
    assert dropped.empty


def test_merge_missing_appends_dropped_rows():
    missing = pd.DataFrame({"design_id": ["blank"], "reason": ["member_empty"]})
    dropped = pd.DataFrame({"design_id": ["ghost"], "reason": [DUPLICATE_REASON],
                            "coord_hash": ["same"]})
    out = merge_missing(missing, dropped)
    assert list(out["design_id"]) == ["blank", "ghost"]
    assert list(out["reason"]) == ["member_empty", DUPLICATE_REASON]
    assert pd.isna(out["coord_hash"].iloc[0])
    assert out["coord_hash"].iloc[1] == "same"


def test_merge_missing_with_nothing_dropped_is_the_original_table():
    missing = pd.DataFrame({"design_id": ["blank"], "reason": ["member_empty"]})
    out = merge_missing(missing, pd.DataFrame(
        columns=["design_id", "reason", "coord_hash"]))
    assert list(out["design_id"]) == ["blank"]
    assert list(out.columns) == ["design_id", "reason"]


def test_reconcile_all_accounted_for():
    """All designs either parsed or in missing: no unexplained shortfall."""
    keep = {"a", "b", "c", "d"}
    explained, unexplained = reconcile(keep, hashes_rows=3, missing_rows=1)
    assert explained == 4
    assert unexplained == 0


def test_reconcile_partial_missing():
    """Some designs in missing/: they are accounted for, not unexplained."""
    keep = {"a", "b", "c", "d", "e"}
    explained, unexplained = reconcile(keep, hashes_rows=3, missing_rows=2)
    assert explained == 5
    assert unexplained == 0


def test_reconcile_unexplained_shortfall():
    """Silent gap: designs neither parsed nor in missing/."""
    keep = {"a", "b", "c", "d"}
    explained, unexplained = reconcile(keep, hashes_rows=2, missing_rows=1)
    assert explained == 3
    assert unexplained == 1


def test_reconcile_all_in_missing():
    """All designs are in missing/ (complete failure): still accounted for."""
    keep = {"a", "b"}
    explained, unexplained = reconcile(keep, hashes_rows=0, missing_rows=2)
    assert explained == 2
    assert unexplained == 0


def test_should_warn_census_on_unexplained():
    """Census map with unexplained shortfall: should warn."""
    assert should_warn_census(is_census=True, unexplained=5) is True


def test_should_warn_census_on_zero_unexplained():
    """Census map with no unexplained shortfall: should not warn."""
    assert should_warn_census(is_census=True, unexplained=0) is False


def test_should_warn_census_non_census_with_unexplained():
    """Non-census map with unexplained shortfall: should not warn."""
    assert should_warn_census(is_census=False, unexplained=5) is False


def test_should_warn_census_non_census_zero():
    """Non-census map with no unexplained shortfall: should not warn."""
    assert should_warn_census(is_census=False, unexplained=0) is False


# main()-level tests.
#
# The live failure these pin down: stage 1 appends a hashes row for every
# design_id but skips writing the structure when that coord_hash was already
# seen in the task, so the hashes table is a superset of the shards. Stage 2
# elects the lexicographically FIRST design_id per coord_hash, and when that
# election lands on a design whose structure was deduped away the Foldseek DB
# has no key for it and stage 2 aborts:
#     resolved 499998 / 499999
#     FATAL: 1 canonical structures are not in the Foldseek DB.
# `ghost` sorts before `keeper`, so `ghost` is exactly the id stage 2 would
# elect -- and it is the one with no `.pdb`.


def _shared(root: Path, *, members: dict, hashes: pd.DataFrame,
            missing: pd.DataFrame) -> Path:
    """Build the minimal `results_shared/design` tree stage 1b reads."""
    src = root / "results_shared" / "design"
    for sub in ("shards", "hashes", "fp", "missing"):
        (src / sub).mkdir(parents=True, exist_ok=True)
    _tar(src / "shards" / "shard_00.tar.gz", members)
    hashes.to_csv(src / "hashes" / "shard_00.tsv", sep="\t", index=False)
    missing.to_csv(src / "missing" / "shard_00.tsv", sep="\t", index=False)
    return root / "results_shared"


def _manifest(root: Path, design_ids) -> None:
    build = root / "build" / "manifests"
    build.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"design_id": list(design_ids)}).to_csv(
        build / "design.tsv", sep="\t", index=False)


def _run_main(root: Path, monkeypatch, *, census=True):
    """Invoke main() on `root`, return (hashes, missing) read back from disk."""
    spec = "lib12k:build:results" + (":census" if census else "")
    monkeypatch.setattr(sys, "argv", [
        "stage1b_split.py",
        "--shared", str(root / "results_shared"),
        "--root", str(root),
        "--map", spec,
    ])
    assert main() == 0
    out = root / "results" / "design"
    return (pd.read_csv(out / "hashes" / "shard_00.tsv", sep="\t"),
            pd.read_csv(out / "missing" / "shard_00.tsv", sep="\t"))


def _no_missing():
    return pd.DataFrame(columns=["design_id", "reason"])


def test_main_drops_the_hashes_row_whose_structure_was_deduped_away(
        tmp_path, monkeypatch):
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper"])

    hashes, _missing = _run_main(tmp_path, monkeypatch)

    assert list(hashes["design_id"]) == ["keeper"]


def test_main_records_the_dropped_design_in_missing_with_its_hash(
        tmp_path, monkeypatch):
    """Dropped, not lost: the row keeps a reason and the coord_hash, so the
    design still points at the node whose structure represents it."""
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper"])

    _hashes_out, missing = _run_main(tmp_path, monkeypatch)

    assert list(missing["design_id"]) == ["ghost"]
    assert list(missing["reason"]) == [DUPLICATE_REASON]
    assert list(missing["coord_hash"]) == ["dd7384b2"]


def test_main_keeps_the_reconcile_identity_after_dropping(tmp_path, monkeypatch):
    """hashes_rows + missing_rows == len(keep), still, with stage 1's own
    missing rows preserved alongside the newly dropped one."""
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=pd.DataFrame({"design_id": ["blank"],
                                  "reason": ["member_empty"]}))
    _manifest(tmp_path, ["blank", "ghost", "keeper"])

    hashes, missing = _run_main(tmp_path, monkeypatch)

    assert sorted(missing["design_id"]) == ["blank", "ghost"]
    explained, unexplained = reconcile({"blank", "ghost", "keeper"},
                                       len(hashes), len(missing))
    assert (explained, unexplained) == (3, 0)


def test_main_leaves_no_hashes_row_without_a_structure(tmp_path, monkeypatch):
    """The invariant stage 2 needs: every id it can elect has a DB key."""
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper"])

    hashes, _missing = _run_main(tmp_path, monkeypatch)

    shard = tmp_path / "results" / "design" / "shards" / "shard_00.tar.gz"
    with tarfile.open(shard, "r:gz") as tar:
        stems = {m.name[:-4] for m in tar if m.isfile()}
    assert set(hashes["design_id"]) <= stems


def test_main_census_is_quiet_when_the_drop_is_accounted_for(
        tmp_path, monkeypatch, capsys):
    """A deduped structure is not a data loss event: no census WARNING."""
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper"])

    _run_main(tmp_path, monkeypatch)

    assert "WARNING" not in capsys.readouterr().out


def test_main_census_still_warns_on_a_real_shortfall(tmp_path, monkeypatch, capsys):
    """`lost` is in the manifest but in neither table: the census check must
    still fire, and the dropped `ghost` row must not paper over it."""
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper", "lost"])

    _run_main(tmp_path, monkeypatch)

    assert "WARNING: 1 library designs unaccounted for" in capsys.readouterr().out


def test_main_non_census_does_not_warn_on_a_shortfall(tmp_path, monkeypatch, capsys):
    _shared(tmp_path, members={"keeper.pdb": b"ATOM  keeper\n"},
            hashes=_hashes([("ghost", "dd7384b2"), ("keeper", "dd7384b2")]),
            missing=_no_missing())
    _manifest(tmp_path, ["ghost", "keeper", "lost"])

    _run_main(tmp_path, monkeypatch, census=False)

    assert "WARNING" not in capsys.readouterr().out
