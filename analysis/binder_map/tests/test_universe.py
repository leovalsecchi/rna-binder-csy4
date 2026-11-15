import pytest
from binder_map.universe import ONE_BACKBONE_PER_DESIGN, backbone_key


def test_limpnn_stacks_strip_the_seq_suffix():
    assert backbone_key(
        "csy4free_RFD3_LiMPNN_len80_a546ac7c_design31_seq0", "RFD3+LiMPNN"
    ) == "csy4free_RFD3_LiMPNN_len80_a546ac7c_design31"


def test_seq_variants_of_one_backbone_collapse_to_one_key():
    ids = [f"csy4free_BG_LiMPNN_len80_b2f7f57e_design13_seq{i}" for i in range(10)]
    assert len({backbone_key(i, "BG+LiMPNN") for i in ids}) == 1


def test_bd1_also_collapses():
    ids = [f"csy4free_BD1_LiMPNN_len80_000158e9_design10_seq{i}" for i in range(3)]
    assert len({backbone_key(i, "BD1+LiMPNN") for i in ids}) == 1


def test_bg_bif_is_one_backbone_per_design():
    # BoltzIF is generative all-atom: 1,629 of 6,082 name-prefixes in the v6
    # ground truth split into multiple coordinate hashes.
    ids = [f"csy4free_BG_BIF_len80_0000df81_design0_seq{i}" for i in range(4)]
    assert len({backbone_key(i, "BG+BIF") for i in ids}) == 4


def test_hd_is_one_backbone_per_design():
    # All 802 HD designs in the v6 pool carry distinct coord hashes across
    # 447 (uuid, recycle) groups -- sample0..9 are not one backbone.
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample1_pack1"
    assert backbone_key(a, "HD+LiMPNN") != backbone_key(b, "HD+LiMPNN")


def test_key_is_idempotent():
    key = backbone_key("csy4free_RFD3_LiMPNN_len80_a546ac7c_design31_seq0", "RFD3+LiMPNN")
    assert backbone_key(key, "RFD3+LiMPNN") == key


def test_one_backbone_stacks_are_exactly_these_two():
    assert ONE_BACKBONE_PER_DESIGN == frozenset({"BG+BIF", "HD+LiMPNN"})


def test_unknown_stack_raises():
    with pytest.raises(KeyError):
        backbone_key("whatever_seq0", "NOT+A+STACK")
