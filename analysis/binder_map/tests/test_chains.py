import pytest

from binder_map.chains import (
    ChainSelectError, load_structure, protein_chain_id, split_chains,
)

RNA_RESNAMES = {"A", "C", "G", "U"}


def _split(path):
    return split_chains(load_structure(path))


def test_rfd3_design_protein_is_chain_B(rfd3_design):
    """RFD3+LiMPNN design CIFs are inverted: chain A is the RNA target.

    This is 53% of the design layer. If this test ever fails, something
    upstream changed and the pipeline must stop.
    """
    atoms = load_structure(rfd3_design)
    assert protein_chain_id(atoms) == "B"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80
    assert set(rna.res_name) <= RNA_RESNAMES


def test_bd1_design_protein_is_chain_A(bd1_design):
    atoms = load_structure(bd1_design)
    assert protein_chain_id(atoms) == "A"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80


def test_bgbif_dialect_parses(bgbif_design):
    """BG+BIF uses a different mmCIF dialect; a hand-rolled parser failed on it."""
    protein, rna = _split(bgbif_design)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80
    assert set(rna.res_name) <= RNA_RESNAMES


@pytest.mark.parametrize("layer", ["af3_refold", "esm_refold"])
def test_refolds_protein_is_chain_A(layer, request):
    atoms = load_structure(request.getfixturevalue(layer))
    assert protein_chain_id(atoms) == "A"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80


def test_all_fixtures_give_80_residue_protein_and_rna(
    rfd3_design, bd1_design, bgbif_design, af3_refold, esm_refold
):
    for path in (rfd3_design, bd1_design, bgbif_design, af3_refold, esm_refold):
        protein, rna = _split(path)
        assert len(set(zip(protein.chain_id, protein.res_id))) == 80, path
        assert len(rna) > 0, path
        assert set(rna.res_name) <= RNA_RESNAMES, path


def test_protein_only_structure_raises(bd1_design):
    atoms = load_structure(bd1_design)
    protein, _ = split_chains(atoms)
    with pytest.raises(ChainSelectError):
        split_chains(protein)
