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
