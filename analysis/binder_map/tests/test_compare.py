import numpy as np
import pytest

from binder_map.compare import (
    contingency, fingerprint_agreement, partition_scores, stack_profiles,
)


def test_identical_partitions_score_one():
    a = np.array([0, 0, 1, 1, 2, 2])
    s = partition_scores(a, a.copy())
    assert s["ami"] == pytest.approx(1.0)
    assert s["ari"] == pytest.approx(1.0)
    assert s["n"] == 6


def test_relabelled_partition_still_scores_one():
    """Cluster ids are arbitrary; only the grouping matters."""
    a = np.array([0, 0, 1, 1])
    b = np.array([7, 7, 3, 3])
    assert partition_scores(a, b)["ari"] == pytest.approx(1.0)


def test_independent_partitions_score_near_zero():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 5, size=2000)
    b = rng.integers(0, 5, size=2000)
    s = partition_scores(a, b)
    assert abs(s["ami"]) < 0.05
    assert abs(s["ari"]) < 0.05


def test_contingency_counts_co_membership():
    a = np.array([0, 0, 1])
    b = np.array([5, 6, 6])
    mat, la, lb = contingency(a, b)
    assert mat.shape == (len(la), len(lb))
    assert mat.sum() == 3
    assert mat[list(la).index(0), list(lb).index(5)] == 1
    assert mat[list(la).index(0), list(lb).index(6)] == 1


def test_agreement_of_identical_fingerprints_is_one():
    fp = np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 3.0]])
    np.testing.assert_allclose(fingerprint_agreement(fp, fp.copy()), [1.0, 1.0],
                               atol=1e-6)


def test_agreement_of_orthogonal_fingerprints_is_zero():
    a = np.array([[1.0, 0.0, 0.0]])
    b = np.array([[0.0, 1.0, 0.0]])
    assert fingerprint_agreement(a, b)[0] == pytest.approx(0.0)


def test_agreement_is_zero_when_a_row_has_no_contacts():
    """A design with no interface has no epitope; it must not score as
    agreeing with anything."""
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[1.0, 2.0, 3.0]])
    assert fingerprint_agreement(a, b)[0] == 0.0


def test_agreement_ignores_magnitude():
    a = np.array([[1.0, 2.0, 0.0]])
    b = np.array([[10.0, 20.0, 0.0]])
    assert fingerprint_agreement(a, b)[0] == pytest.approx(1.0)


def test_stack_profiles_averages_within_stack():
    per_nt = np.array([
        [2.0] + [0.0] * 15,
        [4.0] + [0.0] * 15,
        [0.0, 6.0] + [0.0] * 14,
    ])
    stacks = np.array(["A", "A", "B"])
    df = stack_profiles(per_nt, stacks)
    assert df.loc["A", "nt01"] == pytest.approx(3.0)
    assert df.loc["B", "nt02"] == pytest.approx(6.0)
    assert list(df.columns) == [f"nt{i:02d}" for i in range(1, 17)]
