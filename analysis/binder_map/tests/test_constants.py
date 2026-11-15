from binder_map import (
    RNA_CANONICAL, RNA_LEN, RNA_CONTACT_CUTOFF, FINGERPRINT_DIM, MOIETIES,
)


def test_canonical_rna_is_16nt():
    assert len(RNA_CANONICAL) == RNA_LEN == 16
    assert set(RNA_CANONICAL) <= {"A", "C", "G", "U"}


def test_contact_cutoff_matches_rna_binder():
    # rna_binder.analysis.structural.RNA_CONTACT_CUTOFF
    assert RNA_CONTACT_CUTOFF == 4.5


def test_fingerprint_dim():
    assert MOIETIES == ("phosphate", "sugar", "base")
    assert FINGERPRINT_DIM == 48
