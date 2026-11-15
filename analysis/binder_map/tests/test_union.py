import fnmatch
import io
import sys
import tarfile

import pandas as pd
import pytest
from pathlib import Path

import stageU0_union
from stageU0_union import (DUPLICATE_REASON, MISSING_COLUMNS, dedup_rows,
                           filter_hashes, install_shard, layer_flags, main,
                           present_stems, shard_name, source_tag)


def _shard(path: Path, stems: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as tar:
        for s in stems:
            body = b"ATOM\n"
            info = tarfile.TarInfo(f"{s}.pdb")
            info.size = len(body)
            import io
            tar.addfile(info, io.BytesIO(body))


def test_present_stems_reads_every_shard(tmp_path):
    _shard(tmp_path / "shards" / "shard_00.tar.gz", ["a", "b"])
    _shard(tmp_path / "shards" / "shard_01.tar.gz", ["c"])
    assert present_stems(tmp_path / "shards") == {"a", "b", "c"}


def test_present_stems_ignores_non_pdb(tmp_path):
    p = tmp_path / "shards" / "shard_00.tar.gz"
    p.parent.mkdir(parents=True)
    import io
    with tarfile.open(p, "w:gz") as tar:
        for name in ("a.pdb", "notes.txt"):
            body = b"x"
            info = tarfile.TarInfo(name)
            info.size = len(body)
            tar.addfile(info, io.BytesIO(body))
    assert present_stems(p.parent) == {"a"}


def test_filter_hashes_drops_rows_without_a_structure():
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h1"]})
    kept, dropped = filter_hashes(frame, {"b"})
    assert list(kept["design_id"]) == ["b"]
    assert list(dropped["design_id"]) == ["a"]


def test_filter_hashes_keeps_everything_when_all_present():
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h2"]})
    kept, dropped = filter_hashes(frame, {"a", "b"})
    assert len(kept) == 2 and len(dropped) == 0


def test_filter_hashes_never_drops_a_whole_hash_group():
    # every coord_hash that had a row must still have one afterwards
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h1"]})
    kept, _ = filter_hashes(frame, {"a"})
    assert set(kept["coord_hash"]) == {"h1"}


def test_layer_flags_marks_conserved_nodes():
    hash_of = {"d1": "h1", "d2": "h2", "d3": "h3"}
    layers = {"sample500k": {"d1", "d2"}, "passing60k": {"d2"}, "library12k": {"d3"}}
    out = layer_flags(layers, hash_of).set_index("coord_hash")
    assert bool(out.loc["h2", "in_sample500k"]) and bool(out.loc["h2", "in_passing60k"])
    assert not bool(out.loc["h1", "in_passing60k"])
    assert bool(out.loc["h3", "in_12k"])
    assert len(out) == 3


def test_layer_flags_every_node_carries_at_least_one_flag():
    hash_of = {"d1": "h1"}
    out = layer_flags({"sample500k": {"d1"}, "passing60k": set(), "library12k": set()}, hash_of)
    cols = ["in_sample500k", "in_passing60k", "in_12k"]
    assert bool(out[cols].any(axis=1).all())


# --- present_stems ---------------------------------------------------------


def test_present_stems_on_a_dir_with_no_shards_is_empty(tmp_path):
    (tmp_path / "shards").mkdir()
    assert present_stems(tmp_path / "shards") == set()


def test_present_stems_collapses_a_stem_two_shards_both_carry(tmp_path):
    """Two source trees can hold the same structure; the set says 'present
    once', which is what the hashes filter needs."""
    _shard(tmp_path / "shards" / "shard_shared_00.tar.gz", ["a", "b"])
    _shard(tmp_path / "shards" / "shard_v6_00.tar.gz", ["b", "c"])
    assert present_stems(tmp_path / "shards") == {"a", "b", "c"}


# --- naming ----------------------------------------------------------------


def test_source_tag_strips_the_results_prefix():
    assert source_tag("results_shared") == "shared"
    assert source_tag("results_v6") == "v6"


def test_source_tag_leaves_an_unprefixed_name_alone():
    assert source_tag("scratch_tree") == "scratch_tree"


def test_source_tag_uses_the_last_path_component():
    assert source_tag("/cluster/work/foldseek_map/results_v6") == "v6"


def test_shard_name_inserts_the_tag_after_the_shard_prefix():
    assert shard_name("shared", "shard_00.tar.gz") == "shard_shared_00.tar.gz"
    assert shard_name("v6", "shard_07.tsv") == "shard_v6_07.tsv"


def test_shard_name_still_matches_the_repo_wide_shard_glob():
    """stage2/stage7/validate/gate_checks all glob `shard_*`; a name that
    starts with the source instead would be invisible to every one of them."""
    for pattern, original in (("shard_*.tar.gz", "shard_00.tar.gz"),
                              ("shard_*.tsv", "shard_00.tsv"),
                              ("shard_*.parquet", "shard_00.parquet")):
        assert fnmatch.fnmatch(shard_name("shared", original), pattern)


def test_shard_name_keeps_two_sources_apart():
    assert shard_name("shared", "shard_00.tar.gz") != shard_name("v6", "shard_00.tar.gz")


# --- filter_hashes ---------------------------------------------------------


def _hashes(pairs):
    """A stage 1 hashes table: [(design_id, coord_hash), ...]."""
    return pd.DataFrame({
        "design_id": [d for d, _ in pairs],
        "coord_hash": [h for _, h in pairs],
        "protein_chain": ["A"] * len(pairs),
        "n_res": [80] * len(pairs),
    })


def test_filter_hashes_keeps_the_other_columns():
    kept, _ = filter_hashes(_hashes([("a", "h1")]), {"a"})
    assert list(kept.columns) == ["design_id", "coord_hash", "protein_chain", "n_res"]


def test_filter_hashes_dropped_rows_carry_the_missing_schema():
    _kept, dropped = filter_hashes(_hashes([("ghost", "h1")]), set())
    assert list(dropped.columns) == MISSING_COLUMNS
    assert list(dropped["reason"]) == [DUPLICATE_REASON]
    assert list(dropped["coord_hash"]) == ["h1"]


def test_filter_hashes_on_an_empty_present_set_keeps_nothing():
    kept, dropped = filter_hashes(_hashes([("a", "h1"), ("b", "h2")]), set())
    assert len(kept) == 0 and len(dropped) == 2


# --- dedup_rows ------------------------------------------------------------


def test_dedup_rows_collapses_a_design_both_sources_carry():
    """6,377 design_ids appear in both source hash tables. Left alone they
    would be counted twice in every node's n_designs."""
    frame = _hashes([("a", "h1"), ("a", "h1"), ("b", "h2")])
    out = dedup_rows(frame)
    assert list(out["design_id"]) == ["a", "b"]


def test_dedup_rows_keeps_one_row_per_distinct_design():
    frame = _hashes([("a", "h1"), ("b", "h2"), ("c", "h3")])
    assert len(dedup_rows(frame)) == 3


def test_dedup_rows_keeps_two_designs_that_share_a_hash():
    frame = _hashes([("a", "h1"), ("b", "h1")])
    assert list(dedup_rows(frame)["design_id"]) == ["a", "b"]


def test_dedup_rows_rejects_one_design_with_two_hashes():
    """Same name, different backbone: silently keeping the first would put a
    structure under the wrong node."""
    frame = _hashes([("a", "h1"), ("a", "h2")])
    with pytest.raises(ValueError):
        dedup_rows(frame)


# --- install_shard ---------------------------------------------------------


def _tar(path: Path, members: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            tar.addfile(info, io.BytesIO(data))


def _members(path: Path) -> set[str]:
    with tarfile.open(path, "r:gz") as tar:
        return {m.name for m in tar if m.isfile()}


def test_install_shard_hard_links_when_nothing_collides(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"A\n", "b.pdb": b"B\n"})
    dst = tmp_path / "out" / "shard_shared_00.tar.gz"
    dst.parent.mkdir()
    stems, rewritten = install_shard(src, dst, set())
    assert rewritten is False
    assert stems == {"a", "b"}
    assert dst.stat().st_ino == src.stat().st_ino


def test_install_shard_rewrites_when_a_stem_is_already_claimed(tmp_path):
    """The collision that kills stage 2: two tars carrying `b.pdb` put two
    entries under one name in db/full.lookup, so the canonical name resolves
    to two keys and `resolved N / KEEP` overshoots."""
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"b.pdb": b"B\n", "c.pdb": b"C\n"})
    dst = tmp_path / "out" / "shard_v6_00.tar.gz"
    dst.parent.mkdir()
    stems, rewritten = install_shard(src, dst, {"b"})
    assert rewritten is True
    assert stems == {"c"}
    assert _members(dst) == {"c.pdb"}


def test_install_shard_rewrite_is_a_new_file_not_a_link(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"b.pdb": b"B\n", "c.pdb": b"C\n"})
    dst = tmp_path / "out" / "shard_v6_00.tar.gz"
    dst.parent.mkdir()
    install_shard(src, dst, {"b"})
    assert dst.stat().st_ino != src.stat().st_ino


def test_install_shard_drops_a_whole_shard_that_is_entirely_claimed(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"A\n"})
    dst = tmp_path / "out" / "shard_v6_00.tar.gz"
    dst.parent.mkdir()
    stems, rewritten = install_shard(src, dst, {"a"})
    assert (stems, rewritten) == (set(), True)
    assert _members(dst) == set()


def test_install_shard_reports_only_the_stems_it_installed(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"A\n", "b.pdb": b"B\n", "c.pdb": b"C\n"})
    dst = tmp_path / "out" / "shard_v6_00.tar.gz"
    dst.parent.mkdir()
    stems, _ = install_shard(src, dst, {"a", "zz"})
    assert stems == {"b", "c"}


# --- layer_flags -----------------------------------------------------------


def test_layer_flags_columns_are_in_the_documented_order():
    out = layer_flags({"sample500k": {"d1"}, "passing60k": set(), "library12k": set()},
                      {"d1": "h1"})
    assert list(out.columns) == ["coord_hash", "in_sample500k", "in_passing60k", "in_12k"]


def test_layer_flags_collapses_two_designs_of_one_backbone():
    """Each map elects its own canonical design, so the same backbone reaches
    layer_flags under two names; it is still one node."""
    out = layer_flags({"sample500k": {"d1"}, "passing60k": {"d2"}, "library12k": set()},
                      {"d1": "h1", "d2": "h1"})
    assert len(out) == 1
    assert bool(out.loc[0, "in_sample500k"]) and bool(out.loc[0, "in_passing60k"])


def test_layer_flags_rejects_an_unknown_layer_name():
    with pytest.raises(ValueError):
        layer_flags({"sample500k": {"d1"}, "passing60k": set(),
                     "library12k": set(), "extra": set()}, {"d1": "h1"})


def test_layer_flags_rejects_a_missing_layer_name():
    with pytest.raises(ValueError):
        layer_flags({"sample500k": {"d1"}, "passing60k": set()}, {"d1": "h1"})


def test_layer_flags_rejects_a_node_with_no_coord_hash():
    """A layer name that the union hash table does not carry would silently
    lose its flag, and the conserved counts would come out low."""
    with pytest.raises(ValueError):
        layer_flags({"sample500k": {"d1", "stranger"}, "passing60k": set(),
                     "library12k": set()}, {"d1": "h1"})


# --- main() ----------------------------------------------------------------
#
# The fixture tree below is the live shape in miniature:
#   shared: structures a, b        hashes a->h1, b->h2, ghost->h1
#   v6:     structures b, c        hashes b->h2, c->h3
# `ghost` is the 2026-09-24 killer -- a hashes row whose structure was deduped
# away, and the lexicographically first id for h1, so stage 2 would elect it.
# `b` is the second killer this assembly has to defend against: it has a .pdb
# in BOTH trees, so a straight copy would put `b.pdb` in the Foldseek DB twice.

EXPECT_UNION = 3
EXPECT_LAYER = {"sample500k": 1, "passing60k": 2, "library12k": 1}
EXPECT_CONSERVED = {
    ("sample500k", "passing60k"): 0,
    ("sample500k", "library12k"): 0,
    ("passing60k", "library12k"): 1,
    ("sample500k", "passing60k", "library12k"): 0,
}


def _source(root: Path, name: str, stems: list[str], pairs, missing=None) -> None:
    src = root / name / "design"
    for sub in ("shards", "hashes", "fp", "missing"):
        (src / sub).mkdir(parents=True, exist_ok=True)
    _tar(src / "shards" / "shard_00.tar.gz",
         {f"{s}.pdb": b"ATOM\n" for s in stems})
    frame = _hashes(pairs)
    frame.to_csv(src / "hashes" / "shard_00.tsv", sep="\t", index=False)
    if missing is None:
        missing = pd.DataFrame(columns=["design_id", "reason"])
    missing.to_csv(src / "missing" / "shard_00.tsv", sep="\t", index=False)
    pd.DataFrame({"design_id": list(frame["design_id"]),
                  "fp0": [0.0] * len(frame)}).to_parquet(
        src / "fp" / "shard_00.parquet", index=False)


def _nodes(root: Path, name: str, struct_names: list[str]) -> str:
    dedup = root / name / "design" / "dedup"
    dedup.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"node_id": range(len(struct_names)),
                  "struct_name": struct_names,
                  "n_designs": [1] * len(struct_names)}).to_csv(
        dedup / "nodes.tsv", sep="\t", index=False)
    return f"{name}/design/dedup/nodes.tsv"


def _tree(root: Path, *, layers=None) -> list[str]:
    _source(root, "results_shared", ["a", "b"],
            [("a", "h1"), ("b", "h2"), ("ghost", "h1")])
    _source(root, "results_v6", ["b", "c"], [("b", "h2"), ("c", "h3")])
    layers = layers or {"sample500k": ["a"], "passing60k": ["b", "c"],
                        "library12k": ["c"]}
    return [f"{k}:{_nodes(root, f'results_{k}', v)}" for k, v in layers.items()]


def _run_main(root: Path, monkeypatch, *, layers=None, expect_union=EXPECT_UNION):
    specs = _tree(root, layers=layers)
    monkeypatch.setattr(stageU0_union, "EXPECTED_UNION", expect_union)
    monkeypatch.setattr(stageU0_union, "EXPECTED_LAYER", EXPECT_LAYER)
    monkeypatch.setattr(stageU0_union, "EXPECTED_CONSERVED", EXPECT_CONSERVED)
    argv = ["stageU0_union.py", "--root", str(root),
            "--source", "results_shared", "--source", "results_v6",
            "--out", "results_union"]
    for spec in specs:
        argv += ["--layer-nodes", spec]
    monkeypatch.setattr(sys, "argv", argv)
    return main()


def _union(root: Path) -> Path:
    return root / "results_union" / "design"


def _read_hashes(root: Path) -> pd.DataFrame:
    return pd.concat([pd.read_csv(f, sep="\t") for f in
                      sorted((_union(root) / "hashes").glob("shard_*.tsv"))],
                     ignore_index=True)


def test_main_returns_zero_when_the_counts_match(tmp_path, monkeypatch):
    assert _run_main(tmp_path, monkeypatch) == 0


def test_main_drops_the_hashes_row_whose_structure_is_absent(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    assert sorted(_read_hashes(tmp_path)["design_id"]) == ["a", "b", "c"]


def test_main_records_the_drop_in_missing_with_a_reason(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    missing = pd.concat([pd.read_csv(f, sep="\t") for f in
                         sorted((_union(tmp_path) / "missing").glob("shard_*.tsv"))],
                        ignore_index=True)
    assert list(missing["design_id"]) == ["ghost"]
    assert list(missing["reason"]) == [DUPLICATE_REASON]
    assert list(missing["coord_hash"]) == ["h1"]


def test_main_leaves_no_hashes_row_without_a_structure(tmp_path, monkeypatch):
    """The invariant stage 2 needs: every id it can elect has a DB key."""
    _run_main(tmp_path, monkeypatch)
    assert set(_read_hashes(tmp_path)["design_id"]) <= present_stems(
        _union(tmp_path) / "shards")


def test_main_writes_each_structure_exactly_once(tmp_path, monkeypatch):
    """`b` lives in both source trees. Two `b.pdb` members across the union
    shards means two db/full.lookup rows under one name, and stage 2's
    `resolved N / KEEP` check fails on the overshoot."""
    _run_main(tmp_path, monkeypatch)
    names = []
    for shard in sorted((_union(tmp_path) / "shards").glob("shard_*.tar.gz")):
        with tarfile.open(shard, "r:gz") as tar:
            names += [m.name for m in tar if m.isfile()]
    assert sorted(names) == ["a.pdb", "b.pdb", "c.pdb"]


def test_main_keeps_one_hashes_row_per_design(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    hashes = _read_hashes(tmp_path)
    assert not hashes["design_id"].duplicated().any()


def test_main_shard_names_stay_inside_the_shard_glob(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    shards = sorted(p.name for p in (_union(tmp_path) / "shards").glob("shard_*.tar.gz"))
    assert shards == ["shard_shared_00.tar.gz", "shard_v6_00.tar.gz"]


def test_main_copies_the_fingerprints_from_both_sources(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    fps = sorted(p.name for p in (_union(tmp_path) / "fp").glob("shard_*.parquet"))
    assert fps == ["shard_shared_00.parquet", "shard_v6_00.parquet"]


def test_main_writes_a_layers_table_with_one_row_per_node(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    layers = pd.read_parquet(_union(tmp_path) / "meta" / "layers.parquet")
    assert len(layers) == EXPECT_UNION
    assert sorted(layers["coord_hash"]) == ["h1", "h2", "h3"]


def test_main_layers_table_marks_the_conserved_node(tmp_path, monkeypatch):
    _run_main(tmp_path, monkeypatch)
    layers = pd.read_parquet(_union(tmp_path) / "meta" / "layers.parquet"
                             ).set_index("coord_hash")
    assert bool(layers.loc["h3", "in_passing60k"]) and bool(layers.loc["h3", "in_12k"])
    assert not bool(layers.loc["h1", "in_12k"])


def test_main_fails_when_the_union_count_is_wrong(tmp_path, monkeypatch, capsys):
    """The assert is the whole point of step 7: a tree that does not come out
    at the measured size must not be handed to stage 2."""
    assert _run_main(tmp_path, monkeypatch, expect_union=999) == 1
    assert "FATAL" in capsys.readouterr().out


def test_main_fails_when_a_layer_count_is_wrong(tmp_path, monkeypatch, capsys):
    layers = {"sample500k": ["a", "b"], "passing60k": ["b", "c"], "library12k": ["c"]}
    assert _run_main(tmp_path, monkeypatch, layers=layers) == 1
    assert "FATAL" in capsys.readouterr().out


def test_main_prints_the_union_size(tmp_path, monkeypatch, capsys):
    _run_main(tmp_path, monkeypatch)
    assert "union: 3 nodes" in capsys.readouterr().out
