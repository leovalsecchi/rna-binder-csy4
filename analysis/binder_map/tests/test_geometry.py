import numpy as np
import pytest

from binder_map.chains import load_structure, split_chains
from binder_map.geometry import approach, reference_frame
from binder_map.rnaframe import canonical_positions


def _pair(path):
    return split_chains(load_structure(path))


ALL_FIXTURES = [
    "af3_refold", "esm_refold", "bd1_design",
    "rfd3_design", "rfd3_design_seq7", "bgbif_design",
]


@pytest.mark.parametrize("name", ALL_FIXTURES)
def test_reference_frame_is_orthonormal(name, request):
    _, rna = _pair(request.getfixturevalue(name))
    centroid, rot = reference_frame(rna)
    assert centroid.shape == (3,)
    assert rot.shape == (3, 3)
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-6)
    # Right-handed, not merely orthogonal: det must be +1, never -1.
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-6)


def test_approach_shapes(af3_refold):
    protein, rna = _pair(af3_refold)
    a = approach(protein, rna)
    assert a.centroid.shape == (3,)
    assert a.axis.shape == (3,)
    assert np.isclose(np.linalg.norm(a.axis), 1.0, atol=1e-6)


@pytest.mark.parametrize("name", ALL_FIXTURES)
def test_approach_is_invariant_to_rigid_motion(name, request):
    """Rotating and translating the whole complex must not change the result.

    This is the entire point of expressing the binder in the RNA's frame, and
    it MUST be checked on every fixture. An earlier version tested only
    af3_refold and passed while the frame was genuinely not invariant: the
    SVD-based implementation flipped an axis on three other fixtures, giving
    |axis.axis'| of 0.85-0.93 instead of 1.0. One fixture is not a property
    test.
    """
    protein, rna = _pair(request.getfixturevalue(name))
    before = approach(protein, rna)

    theta = 0.7
    rot = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    shift = np.array([13.0, -7.0, 42.0])

    p2, r2 = protein.copy(), rna.copy()
    p2.coord = (p2.coord @ rot.T + shift).astype(p2.coord.dtype)
    r2.coord = (r2.coord @ rot.T + shift).astype(r2.coord.dtype)
    after = approach(p2, r2)

    np.testing.assert_allclose(before.centroid, after.centroid, atol=1e-3)
    # The frame makes no sign choice, so the axis must agree in DIRECTION,
    # not merely as an undirected line. abs() here would hide a flip.
    assert float(before.axis @ after.axis) > 0.999, (
        f"{name}: axis direction changed under rigid motion"
    )


def test_approach_ignores_sidechains(af3_refold):
    """Sidechain presence must not move the centroid or axis.

    Three of five design pipelines ship structures without sidechains while
    all refolds have them. If `approach` used every atom, that asymmetry
    would masquerade as a real difference in binding geometry.
    """
    protein, rna = _pair(af3_refold)
    backbone_only = protein[np.isin(protein.atom_name, ["N", "CA", "C", "O"])]
    full = approach(protein, rna)
    stripped = approach(backbone_only, rna)
    np.testing.assert_allclose(full.centroid, stripped.centroid, atol=1e-6)
    np.testing.assert_allclose(full.axis, stripped.axis, atol=1e-6)


def test_missing_frame_atoms_raise_a_clear_error(af3_refold):
    """A nucleotide with no P/C1'/C4' must fail with THIS module's error.

    Before the frame atoms and their positions were returned together,
    re-deriving positions from the filtered set made such a nucleotide
    disappear from the derived sequence, surfacing an opaque FrameError from
    rnaframe instead. Failing loudly is right; failing from the wrong module
    with an unrelated message is not.
    """
    _, rna = _pair(af3_refold)
    frame_names = ["P", "C1'", "C4'"]
    # Drop the frame atoms of the ENTIRE loop reference set (7-10), not just
    # one nucleotide: centroid_of() raises only when the whole group is gone,
    # since losing one member still leaves a well-defined centroid. Dropping a
    # single nucleotide is a graceful-degradation case, not an error case.
    positions = canonical_positions(rna)
    drop = np.isin(positions, [7, 8, 9, 10]) & np.isin(
        [str(n).strip() for n in rna.atom_name], frame_names
    )
    assert drop.any(), "fixture must have frame atoms at nt 7-10"
    with pytest.raises(ValueError, match="no frame atoms"):
        reference_frame(rna[~drop])


def test_dropping_one_frame_nucleotide_still_works(af3_refold):
    """Losing ONE nucleotide's frame atoms must degrade gracefully, not raise.

    The centroid of nucleotides 7-10 is still well defined with three of the
    four present, so the frame stays usable. This is the companion to the test
    above and pins the boundary between graceful degradation and hard failure.
    """
    _, rna = _pair(af3_refold)
    positions = canonical_positions(rna)
    drop = (positions == 8) & np.isin(
        [str(n).strip() for n in rna.atom_name], ["P", "C1'", "C4'"]
    )
    assert drop.any()
    origin, rot = reference_frame(rna[~drop])
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-6)
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-6)


def test_translating_binder_moves_centroid(af3_refold):
    protein, rna = _pair(af3_refold)
    before = approach(protein, rna)
    p2 = protein.copy()
    p2.coord = (p2.coord + np.array([25.0, 0.0, 0.0])).astype(p2.coord.dtype)
    after = approach(p2, rna)
    assert np.linalg.norm(before.centroid - after.centroid) > 10.0
