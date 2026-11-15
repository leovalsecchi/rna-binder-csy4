import numpy as np
import pytest

from binder_map import FINGERPRINT_DIM, RNA_LEN
from binder_map.chains import load_structure, split_chains
from binder_map.fingerprint import (
    FINGERPRINT_COLUMNS, fingerprint, moiety_of,
)
from binder_map.rnaframe import canonical_positions, rna_sequence


def _fp(path):
    protein, rna = split_chains(load_structure(path))
    return fingerprint(protein, rna)


def test_moiety_classification():
    assert moiety_of("P") == "phosphate"
    assert moiety_of("OP1") == "phosphate"
    assert moiety_of("O5'") == "phosphate"
    assert moiety_of("C1'") == "sugar"
    assert moiety_of("O2'") == "sugar"
    assert moiety_of("N1") == "base"
    assert moiety_of("C5") == "base"


def test_columns_match_vector_layout():
    assert len(FINGERPRINT_COLUMNS) == FINGERPRINT_DIM == 48
    assert FINGERPRINT_COLUMNS[0] == "nt01_phosphate"
    assert FINGERPRINT_COLUMNS[-1] == "nt16_base"


@pytest.mark.parametrize("name", ["af3_refold", "esm_refold", "bd1_design"])
def test_fingerprint_shape_and_consistency(name, request):
    fp = _fp(request.getfixturevalue(name))
    assert fp.vector.shape == (48,)
    assert fp.per_nt.shape == (RNA_LEN,)
    # per_nt is the vector summed over the 3 moieties
    assert np.array_equal(fp.per_nt, fp.vector.reshape(RNA_LEN, 3).sum(axis=1))
    # total accounts for everything, in-frame and out
    assert fp.total == int(fp.per_nt.sum()) + fp.offframe


def test_refold_has_contacts(af3_refold):
    """These designs passed a DockQ >= 0.7 gate, so they must have an interface."""
    fp = _fp(af3_refold)
    assert fp.total > 0
    assert fp.per_nt.sum() > 0


def test_offframe_is_zero_for_16nt_target(af3_refold):
    fp = _fp(af3_refold)
    assert fp.offframe == 0


def test_offframe_counts_contacts_outside_the_canonical_frame(bd1_design):
    """Contacts to nucleotides outside the frame must land in `offframe`.

    If they leaked into the 48-dim vector instead, the epitope frame would
    silently shift and every downstream comparison would be wrong while
    looking normal.

    bd1_design is the only 18-nt fixture, but on the real coordinates its two
    flanking nucleotides sit 4.503 A from the nearest protein atom - 0.003 A
    outside the 4.5 A cutoff - so it yields offframe == 0 as it stands. An
    earlier version of this test asserted only algebraic identities, which
    hold trivially at zero and tested nothing. So: assert the real geometry,
    then move ONLY the out-of-frame atoms into contact and assert where the
    new counts land.
    """
    protein, rna = split_chains(load_structure(bd1_design))
    assert len(rna_sequence(rna)) == 18

    baseline = fingerprint(protein, rna)
    assert baseline.offframe == 0

    positions = canonical_positions(rna)
    assert (positions == 0).any(), "fixture must have out-of-frame atoms"

    moved = rna.copy()
    moved.coord[positions == 0] = protein.coord[0]

    fp = fingerprint(protein, moved)
    assert fp.offframe > 0, "out-of-frame contacts were not counted"
    # Nothing in-frame moved, so the 48-dim vector must be untouched: this is
    # what proves the new counts went to `offframe` and not into the frame.
    assert np.array_equal(fp.vector, baseline.vector)
    assert fp.total == int(fp.per_nt.sum()) + fp.offframe


def test_fingerprint_is_deterministic(af3_refold):
    a, b = _fp(af3_refold), _fp(af3_refold)
    assert np.array_equal(a.vector, b.vector)


def test_distant_chains_give_zero_contacts(af3_refold):
    """Translating the RNA 500 A away must produce an empty fingerprint."""
    protein, rna = split_chains(load_structure(af3_refold))
    rna = rna.copy()
    rna.coord = rna.coord + np.array([500.0, 0.0, 0.0], dtype=rna.coord.dtype)
    fp = fingerprint(protein, rna)
    assert fp.total == 0
    assert fp.vector.sum() == 0
