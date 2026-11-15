import numpy as np

from binder_map.backbone import (
    BACKBONE_ATOMS, backbone_atoms, ca_trace, coord_hash, write_backbone_pdb,
)
from binder_map.chains import load_structure, split_chains


def _protein(path):
    return split_chains(load_structure(path))[0]


def test_backbone_keeps_only_backbone_atoms(af3_refold):
    bb = backbone_atoms(_protein(af3_refold))
    assert set(bb.atom_name) <= set(BACKBONE_ATOMS)
    assert "CA" in set(bb.atom_name)


def test_coord_hash_is_stable(af3_refold):
    p = _protein(af3_refold)
    assert coord_hash(p) == coord_hash(p)
    assert len(coord_hash(p)) == 32


def test_sequence_variants_share_a_backbone_hash(rfd3_design, rfd3_design_seq7):
    """seq4 and seq7 of the same RFD3 design are coordinate-identical.

    Verified 2026-08-21 by CA-coordinate md5 on Euler. This is what makes
    dedup collapse the design layer from 71,866 to ~32,848.
    """
    assert coord_hash(_protein(rfd3_design)) == coord_hash(_protein(rfd3_design_seq7))


def test_different_designs_have_different_hashes(rfd3_design, bd1_design):
    assert coord_hash(_protein(rfd3_design)) != coord_hash(_protein(bd1_design))


def test_coord_hash_ignores_sidechains(af3_refold):
    """Design CIFs for the LiMPNN stacks have no sidechains; refolds do.

    The hash must depend only on N/CA/C/O so the two are comparable.
    """
    p = _protein(af3_refold)
    keep = np.isin(p.atom_name, ["N", "CA", "C", "O", "CB"])
    assert coord_hash(p) == coord_hash(p[keep])


def test_ca_trace_shape(af3_refold):
    assert ca_trace(_protein(af3_refold)).shape == (80, 3)


def test_write_backbone_pdb_roundtrips(af3_refold):
    data = write_backbone_pdb(_protein(af3_refold), "test_design")
    assert data.startswith(b"ATOM") or b"\nATOM" in data
    assert data.rstrip().endswith(b"END")
    assert data.count(b"\nATOM") + data.startswith(b"ATOM") >= 80
