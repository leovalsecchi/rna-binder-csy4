import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from stage0_select import METADATA_COLUMNS, build_manifest, read_passing  # noqa: E402

CSV = Path("/Users/leo/Documents/binder_filter/exports/passing.csv")

pytestmark = pytest.mark.skipif(not CSV.exists(), reason="passing.csv not present")


@pytest.fixture(scope="module")
def df():
    return read_passing(CSV)


def test_reads_all_rows_past_the_two_comment_lines(df):
    assert len(df) == 71866
    assert "design_id" in df.columns
    assert df["design_id"].is_unique


def test_expected_stacks_present(df):
    assert set(df["stack"]) == {
        "RFD3+LiMPNN", "BG+LiMPNN", "BG+BIF", "BD1+LiMPNN", "HD+LiMPNN",
    }


def test_manifest_splits_archive_and_member(df):
    m = build_manifest(df, "esm")
    assert len(m) == len(df)
    assert not m["archive_path"].str.contains("!").any()
    assert m["archive_path"].str.endswith(".tar.zst").all()
    assert m["member_path"].str.endswith(".cif").all()


def test_manifest_is_sorted_by_archive(df):
    """Stage 1 batches by archive; unsorted input would reopen archives."""
    m = build_manifest(df, "af3")
    assert m["archive_path"].is_monotonic_increasing


def test_metadata_columns_all_exist(df):
    for column in METADATA_COLUMNS:
        assert column in df.columns, column


def test_hd_designs_have_no_seq_suffix(df):
    """HD design_ids are recycleN_sampleM_packK, not _seqN.

    Dedup must therefore be by coordinate hash, never by stripping _seqN.
    """
    hd = df[df["stack"] == "HD+LiMPNN"]["design_id"]
    assert len(hd) > 0
    assert not hd.str.contains(r"_seq\d+$", regex=True).any()
