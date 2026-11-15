import pandas as pd
import pytest

from stage9_compare import (cluster_table, crosswalk_scores, expand_members,
                            fold_key, summarise)


def _labels():
    return pd.DataFrame({
        "cluster":   [0, 0, 0, 1, 1, 2],
        "stack":     ["RFD3+LiMPNN"] * 3 + ["BG+BIF"] * 2 + ["HD+LiMPNN"],
        "in_12k":    [True, False, True, False, False, False],
        "in_passing":[True, True, True, True, False, False],
        "n_designs": [10, 10, 10, 1, 1, 1],
    })


def test_summarise_counts_nodes_and_clusters():
    got = summarise(_labels(), "m1")
    assert got["map"] == "m1"
    assert got["nodes"] == 6
    assert got["clusters"] == 3


def test_summarise_counts_singletons():
    # cluster 2 has exactly one member
    assert summarise(_labels(), "m1")["singletons"] == 1


def test_summarise_reports_membership_totals():
    got = summarise(_labels(), "m1")
    assert got["n_12k"] == 2
    assert got["n_passing"] == 4


def test_summarise_reports_largest_cluster():
    assert summarise(_labels(), "m1")["largest_cluster"] == 3


def test_fold_key_collapses_bg_bif_seq_variants():
    ids = [f"csy4free_BG_BIF_len80_0000df81_design0_seq{i}" for i in range(4)]
    assert len({fold_key(i) for i in ids}) == 1


def test_fold_key_collapses_hd_sample_variants():
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample8_pack1"
    assert fold_key(a) == fold_key(b)


def test_fold_key_keeps_different_hd_recycles_apart():
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle9_sample0_pack1"
    assert fold_key(a) != fold_key(b)


def test_expand_members_splits_the_design_id_list():
    members = pd.DataFrame({"struct_name": ["s0", "s1"],
                            "design_ids": ["d0,d1, d2", "d3"]})
    out = expand_members(members)
    assert list(out["design_id"]) == ["d0", "d1", "d2", "d3"]
    assert list(out["struct_name"]) == ["s0", "s0", "s0", "s1"]


def test_crosswalk_scores_are_one_for_an_identical_partition():
    a = pd.Series([0, 0, 1, 1], index=["a", "b", "c", "d"])
    b = pd.Series([7, 7, 9, 9], index=["a", "b", "c", "d"])
    got = crosswalk_scores(a, b)
    assert got["n_shared"] == 4
    assert got["ari"] == pytest.approx(1.0)
    assert got["ami"] == pytest.approx(1.0)


def test_crosswalk_scores_use_only_shared_designs():
    a = pd.Series([0, 0, 1], index=["a", "b", "c"])
    b = pd.Series([5, 5, 6], index=["b", "c", "zz"])
    assert crosswalk_scores(a, b)["n_shared"] == 2


def test_crosswalk_scores_on_no_overlap_are_nan():
    import math
    a = pd.Series([0], index=["a"])
    b = pd.Series([1], index=["z"])
    got = crosswalk_scores(a, b)
    assert got["n_shared"] == 0
    assert math.isnan(got["ari"])


def test_cluster_table_counts_folds_below_nodes_when_a_backbone_repeats():
    """`folds` must collapse the dedup key's over-splitting; `size` must not.

    Cluster 0 holds three nodes that are really ONE fold (BG+BIF seq variants
    of one backbone, which the coord-hash key splits apart). Cluster 1 holds two
    genuinely different folds. If this ever aggregated struct_name instead of
    fold_key, folds would equal size and the inflation would be invisible.
    """
    labels = pd.DataFrame({
        "cluster": [0, 0, 0, 1, 1],
        "struct_name": [
            "csy4free_BG_BIF_len80_0000df81_design0_seq0",
            "csy4free_BG_BIF_len80_0000df81_design0_seq1",
            "csy4free_BG_BIF_len80_0000df81_design0_seq2",
            "csy4free_RFD3_LiMPNN_len80_aaaaaaaa_design1_seq0",
            "csy4free_RFD3_LiMPNN_len80_bbbbbbbb_design2_seq0",
        ],
        "stack": ["BG+BIF"] * 3 + ["RFD3+LiMPNN"] * 2,
        "in_passing": [True] * 5,
        "in_12k": [False] * 5,
        "n_designs": [1] * 5,
    })
    labels["fold_key"] = [fold_key(s) for s in labels["struct_name"]]
    out = cluster_table(labels).set_index("cluster")

    assert out.loc[0, "size"] == 3, "three nodes in cluster 0"
    assert out.loc[0, "folds"] == 1, (
        "three seq variants of one backbone are ONE fold - folds is tracking "
        "nodes, so the over-splitting this column exists to expose is hidden"
    )
    assert out.loc[1, "size"] == 2 and out.loc[1, "folds"] == 2, (
        "two genuinely distinct backbones must stay two folds"
    )
