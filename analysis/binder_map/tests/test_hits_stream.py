import numpy as np
import pandas as pd

from stage4_graph import load_hits

COLUMNS = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
           "lddt", "fident", "alnlen", "evalue", "bits"]


def _write_hits(path, rows):
    frame = pd.DataFrame(rows, columns=COLUMNS)
    frame.to_csv(path, sep="\t", header=False, index=False)


def test_load_hits_maps_names_to_node_ids(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
    ])
    qi, tj, qtm, ttm, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert list(qi) == [0] and list(tj) == [1]
    assert np.isclose(qtm[0], 0.9) and np.isclose(ttm[0], 0.8)
    assert unknown == 0


def test_load_hits_strips_the_pdb_suffix(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
    ])
    qi, _, _, _, _ = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 1


def test_load_hits_drops_and_counts_unknown_names(tmp_path):
    index = pd.Series([0], index=["s0"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s0.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
        ["s0.pdb", "ghost.pdb", 0.7, 0.6, 0.65, 0.6, 0.4, 80, 1e-7, 90],
    ])
    qi, _, _, _, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 1
    assert unknown == 1


def test_load_hits_concatenates_multiple_chunk_files(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100]])
    _write_hits(tmp_path / "hits_01.tsv", [
        ["s1.pdb", "s0.pdb", 0.6, 0.5, 0.55, 0.6, 0.4, 80, 1e-7, 90]])
    qi, tj, _, _, _ = load_hits(
        sorted(tmp_path.glob("hits_*.tsv")), index)
    assert sorted(zip(qi.tolist(), tj.tolist())) == [(0, 1), (1, 0)]


def test_load_hits_returns_compact_dtypes(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100]])
    qi, tj, qtm, ttm, _ = load_hits([tmp_path / "hits_00.tsv"], index)
    assert qi.dtype == np.int32 and tj.dtype == np.int32
    assert qtm.dtype == np.float32 and ttm.dtype == np.float32


def test_load_hits_on_empty_file_returns_empty_arrays(tmp_path):
    index = pd.Series([0], index=["s0"])
    (tmp_path / "hits_00.tsv").write_text("")
    qi, tj, qtm, ttm, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 0 and len(tj) == 0 and unknown == 0


def test_redundancy_threshold_agrees_with_the_float64_reference():
    """An exactly-0.95 hit must be counted redundant, as it was pre-streaming.

    `sim` became float32 when hits started streaming. float32(0.95) is
    0.9499999880..., strictly below the float64 literal 0.95, so a float64
    threshold would silently drop every exactly-0.95 hit the old whole-frame
    float64 path counted. This asserts the new path agrees with that reference.

    It does NOT discriminate against a bare Python float threshold, and should
    not: under NEP 50 a Python float adopts the array's dtype, so that spelling
    is already correct. What it catches is a float64 threshold.
    """
    from stage4_graph import REDUNDANCY_TM

    streamed = np.array([0.95], dtype=np.float32)   # what load_hits now yields
    reference = np.array([0.95], dtype=np.float64)  # what the old path yielded
    assert bool((streamed >= REDUNDANCY_TM)[0]) == bool((reference >= 0.95)[0]), (
        "the redundancy threshold disagrees with the pre-streaming float64 path"
    )
