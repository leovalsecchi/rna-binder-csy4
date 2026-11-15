"""Map an RNA chain onto the canonical 16-nt csy4 frame.

Length variants exist: at least one BG+BIF design carries CACUGCCGUAUAGGCAGC,
which is the canonical 16-mer plus one closing base pair. Contacts are only
comparable across designs if every position refers to the same index, so the
canonical 16-mer is located by exact substring match and flanking nucleotides
are reported out-of-frame rather than renumbered.
"""
from __future__ import annotations

import numpy as np

from binder_map import RNA_CANONICAL, RNA_LEN


class FrameError(Exception):
    """The RNA does not contain the canonical csy4 16-mer."""


def rna_sequence(rna_atoms) -> str:
    """One-letter RNA sequence in residue order."""
    _, names = _residue_index(rna_atoms)
    return "".join(names)


def _residue_index(rna_atoms):
    """Return (per-atom residue ordinal, ordered residue names).

    Residues are ordered by first appearance, which is coordinate order in
    every dialect here, rather than by res_id, which is not guaranteed to
    start at 1.
    """
    keys = list(zip(rna_atoms.chain_id, rna_atoms.res_id, rna_atoms.ins_code))
    order: dict[tuple, int] = {}
    names: list[str] = []
    ordinals = np.empty(len(keys), dtype=int)
    for i, key in enumerate(keys):
        if key not in order:
            order[key] = len(order)
            names.append(str(rna_atoms.res_name[i]).strip())
        ordinals[i] = order[key]
    return ordinals, names


def frame_offset(seq: str) -> int:
    """Index into seq where the canonical 16-mer begins."""
    idx = seq.find(RNA_CANONICAL)
    if idx < 0:
        raise FrameError(
            f"canonical frame {RNA_CANONICAL} not found in RNA sequence {seq!r}"
        )
    return idx


def canonical_positions(rna_atoms) -> np.ndarray:
    """Per-atom canonical nucleotide position 1..16, or 0 if out of frame."""
    ordinals, names = _residue_index(rna_atoms)
    offset = frame_offset("".join(names))
    pos = ordinals - offset + 1
    pos[(pos < 1) | (pos > RNA_LEN)] = 0
    return pos
