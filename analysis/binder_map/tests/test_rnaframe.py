import numpy as np
import pytest

from binder_map import RNA_CANONICAL, RNA_LEN
from binder_map.chains import load_structure, split_chains
from binder_map.rnaframe import (
    FrameError, canonical_positions, frame_offset, rna_sequence,
)


def test_frame_offset_exact_match():
    assert frame_offset(RNA_CANONICAL) == 0


def test_frame_offset_with_closing_pair():
    """A BD1+LiMPNN design carries CACUGCCGUAUAGGCAGC: the 16-mer plus one pair.

    Verified against the committed fixtures 2026-08-23.
    """
    assert frame_offset("C" + RNA_CANONICAL + "C") == 1


def test_frame_offset_rejects_unrelated_sequence():
    with pytest.raises(FrameError):
        frame_offset("AAAAAAAAAAAAAAAA")


@pytest.mark.parametrize(
    "name", ["rfd3_design", "bd1_design", "bgbif_design", "af3_refold", "esm_refold"]
)
def test_every_fixture_aligns_to_canonical_frame(name, request):
    _, rna = split_chains(load_structure(request.getfixturevalue(name)))
    seq = rna_sequence(rna)
    offset = frame_offset(seq)
    assert seq[offset : offset + RNA_LEN] == RNA_CANONICAL


def test_canonical_positions_cover_all_16(af3_refold):
    _, rna = split_chains(load_structure(af3_refold))
    pos = canonical_positions(rna)
    assert len(pos) == len(rna)
    assert set(pos[pos > 0]) == set(range(1, RNA_LEN + 1))


def test_canonical_positions_marks_out_of_frame_as_zero(bd1_design):
    """bd1_design carries the 18-nt variant: 2 nucleotides outside the frame.

    Verified 2026-08-23: bd1_design is CACUGCCGUAUAGGCAGC (offset 1), while
    bgbif_design is a plain 16-mer. This is the only fixture that exercises
    the out-of-frame path, so the assertion is unconditional on purpose.
    """
    _, rna = split_chains(load_structure(bd1_design))
    assert len(rna_sequence(rna)) == 18
    pos = canonical_positions(rna)
    assert (pos == 0).any()
    assert set(pos[pos > 0]) == set(range(1, RNA_LEN + 1))
