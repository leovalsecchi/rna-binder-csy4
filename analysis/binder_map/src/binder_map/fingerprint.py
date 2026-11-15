"""48-dimensional epitope fingerprint: 16 canonical nucleotides x 3 moieties.

Counts protein-heavy-atom <-> RNA-heavy-atom pairs within RNA_CONTACT_CUTOFF.
The cutoff and the per-nucleotide binning match
rna_binder.analysis.structural.rna_contact_counts exactly, so the summed
per-nucleotide counts can be cross-checked against the stored
rna_contact_nt_* columns where those exist.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from binder_map import MOIETIES, RNA_CONTACT_CUTOFF, RNA_LEN
from binder_map.rnaframe import canonical_positions

FINGERPRINT_COLUMNS = [
    f"nt{i:02d}_{m}" for i in range(1, RNA_LEN + 1) for m in MOIETIES
]

_PHOSPHATE = frozenset({"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O5'", "O3'"})
_SUGAR = frozenset({"C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O4'"})


def moiety_of(atom_name: str) -> str:
    """Classify an RNA atom as phosphate, sugar, or base."""
    name = atom_name.strip()
    if name in _PHOSPHATE:
        return "phosphate"
    if name in _SUGAR:
        return "sugar"
    return "base"


@dataclass(frozen=True)
class Fingerprint:
    vector: np.ndarray   # (48,) int32, ordered as FINGERPRINT_COLUMNS
    per_nt: np.ndarray   # (16,) int32
    offframe: int        # contacts to nucleotides outside the canonical frame
    total: int           # all protein-RNA contacts


def _heavy(atoms):
    return atoms[atoms.element != "H"]


def fingerprint(protein, rna) -> Fingerprint:
    """Count protein-RNA heavy-atom contacts binned by nucleotide and moiety."""
    p_heavy = _heavy(protein)
    r_heavy = _heavy(rna)

    vector = np.zeros(len(FINGERPRINT_COLUMNS), dtype=np.int32)
    if len(p_heavy) == 0 or len(r_heavy) == 0:
        return Fingerprint(vector, np.zeros(RNA_LEN, dtype=np.int32), 0, 0)

    positions = canonical_positions(r_heavy)
    moiety_idx = np.array(
        [MOIETIES.index(moiety_of(str(n))) for n in r_heavy.atom_name]
    )

    # For every RNA heavy atom, how many protein heavy atoms are within cutoff.
    tree = cKDTree(p_heavy.coord.astype(np.float64))
    neighbours = tree.query_ball_point(
        r_heavy.coord.astype(np.float64), r=RNA_CONTACT_CUTOFF
    )
    counts = np.fromiter((len(n) for n in neighbours), dtype=np.int64,
                         count=len(neighbours))

    total = int(counts.sum())
    in_frame = positions > 0
    offframe = int(counts[~in_frame].sum())

    flat_idx = (positions[in_frame] - 1) * len(MOIETIES) + moiety_idx[in_frame]
    np.add.at(vector, flat_idx, counts[in_frame].astype(np.int32))

    per_nt = vector.reshape(RNA_LEN, len(MOIETIES)).sum(axis=1).astype(np.int32)
    return Fingerprint(vector, per_nt, offframe, total)
