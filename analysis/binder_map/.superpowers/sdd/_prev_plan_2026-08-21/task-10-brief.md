### Task 10: Comparison metrics

The point of having two axes. All cheap table operations, but they carry the project's actual conclusions, so they are unit-tested against cases with known answers.

**Files:**
- Create: `src/binder_map/compare.py`
- Create: `tests/test_compare.py`

**Interfaces:**
- Consumes: `binder_map.RNA_LEN`.
- Produces:
  - `partition_scores(a: np.ndarray, b: np.ndarray) -> dict[str, float]` — keys `ami`, `ari`, `n`.
  - `contingency(a, b) -> tuple[np.ndarray, np.ndarray, np.ndarray]` — `(matrix, labels_a, labels_b)`.
  - `fingerprint_agreement(fp_a: np.ndarray, fp_b: np.ndarray) -> np.ndarray` — row-wise cosine similarity, `0.0` where either row is all-zero.
  - `stack_profiles(per_nt: np.ndarray, stacks: np.ndarray) -> pandas.DataFrame` — mean per-nucleotide profile per stack, rows indexed by stack, 16 columns.

- [ ] **Step 1: Write the failing test**

Create `tests/test_compare.py`:

```python
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
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.compare'`

- [ ] **Step 3: Implement `src/binder_map/compare.py`**

```python
"""Compare the fold and epitope partitions, and epitopes across layers."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score, adjusted_rand_score,
)

from binder_map import RNA_LEN


def partition_scores(a, b) -> dict[str, float]:
    """AMI and ARI between two clusterings of the same nodes.

    A low score is a real result — it means backbone fold does not predict
    where a binder lands — not a failure.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if len(a) != len(b):
        raise ValueError(f"partitions differ in length: {len(a)} vs {len(b)}")
    return {
        "ami": float(adjusted_mutual_info_score(a, b)),
        "ari": float(adjusted_rand_score(a, b)),
        "n": int(len(a)),
    }


def contingency(a, b):
    """Co-membership counts between two partitions."""
    a = np.asarray(a)
    b = np.asarray(b)
    labels_a, ia = np.unique(a, return_inverse=True)
    labels_b, ib = np.unique(b, return_inverse=True)
    mat = np.zeros((len(labels_a), len(labels_b)), dtype=np.int64)
    np.add.at(mat, (ia, ib), 1)
    return mat, labels_a, labels_b


def fingerprint_agreement(fp_a, fp_b) -> np.ndarray:
    """Row-wise cosine similarity between two fingerprint matrices.

    Rows with no contacts score 0.0: a design with no interface has no
    epitope, and must not be recorded as agreeing with anything.
    """
    a = np.asarray(fp_a, dtype=np.float64)
    b = np.asarray(fp_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"fingerprint shapes differ: {a.shape} vs {b.shape}")

    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    good = (na > 0) & (nb > 0)
    out = np.zeros(len(a), dtype=np.float32)
    out[good] = ((a[good] * b[good]).sum(axis=1) / (na[good] * nb[good])).astype(
        np.float32
    )
    return out


def stack_profiles(per_nt, stacks) -> pd.DataFrame:
    """Mean per-nucleotide contact profile for each stack.

    Turns "RFD3+LiMPNN prefers nt 3-11" into a directly readable curve rather
    than an inference from cluster colours.
    """
    columns = [f"nt{i:02d}" for i in range(1, RNA_LEN + 1)]
    frame = pd.DataFrame(np.asarray(per_nt, dtype=np.float64), columns=columns)
    frame["stack"] = np.asarray(stacks)
    return frame.groupby("stack", sort=True)[columns].mean()
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

Expected: 9 passed.

- [ ] **Step 5: Run the whole suite**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest -v
```

Expected: 67 passed — constants 3, chains 7, rnaframe 10, backbone 7, fingerprint 10,
geometry 4, knn 8, clustering 6, embedding 3, compare 9. (Counts are collected tests,
so parametrized cases count once per parameter, not once per function.) Every library
module is now covered; the remaining tasks are I/O, orchestration and the front end.

- [ ] **Step 6: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/compare.py tests/test_compare.py
git commit -m "Add partition comparison, epitope agreement, and stack profiles"
```

---

