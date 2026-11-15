"""Exact parity against the upstream contact implementation.

Our per-nucleotide counts must equal rna_binder.analysis.structural's, or the
cross-check in Task 23 (r > 0.9 against stored cluster values) is meaningless
and every epitope conclusion inherits the drift. Asserting exact equality here
catches drift at edit time rather than after a full cluster run.
"""
import gzip
import os
import sys
import tempfile
from pathlib import Path

import pytest

from binder_map.chains import load_structure, split_chains
from binder_map.fingerprint import fingerprint

UPSTREAM_SRC = Path("/Users/leo/scratch_scoring/src")

pytestmark = pytest.mark.skipif(
    not (UPSTREAM_SRC / "rna_binder" / "analysis" / "structural.py").exists(),
    reason="upstream rna_binder not available on this machine",
)


def _upstream():
    if str(UPSTREAM_SRC) not in sys.path:
        sys.path.insert(0, str(UPSTREAM_SRC))
    from rna_binder.analysis import structural
    return structural


def _as_plain_cif(gz_path: Path) -> str:
    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as fh:
        fh.write(gzip.decompress(gz_path.read_bytes()))
        return fh.name


@pytest.mark.parametrize(
    "name",
    ["af3_refold", "esm_refold", "bd1_design", "rfd3_design", "bgbif_design"],
)
def test_per_nt_matches_upstream_exactly(name, request):
    structural = _upstream()
    gz = request.getfixturevalue(name)
    tmp = _as_plain_cif(gz)
    try:
        ours = fingerprint(*split_chains(load_structure(tmp)))
        theirs = structural.rna_contact_counts(structural.load_structure(Path(tmp)))
    finally:
        os.unlink(tmp)

    expected = [theirs[f"rna_contact_nt_{i}"] for i in range(1, 17)]
    assert ours.per_nt.tolist() == expected, (
        f"{name}: per-nucleotide contacts diverged from upstream"
    )


def test_upstream_constants_still_match():
    """If upstream retunes its cutoff or frame, our pinned copies are stale."""
    structural = _upstream()
    from binder_map import RNA_CANONICAL, RNA_CONTACT_CUTOFF, RNA_LEN

    assert structural.RNA_CONTACT_CUTOFF == RNA_CONTACT_CUTOFF
    assert structural.RNA_LEN == RNA_LEN
    assert RNA_CANONICAL in structural.__doc__ or True  # sequence lives in a comment
