"""Describe where the binder sits relative to the RNA, in the RNA's own frame.

Target-fixed: the RNA defines the coordinate system and is never re-fit to the
binder. This is the same RNA-aligned convention adopted for DockQ on
2026-08-04, so these numbers are comparable to the pose metrics already in use.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from binder_map.backbone import ca_trace
from binder_map.rnaframe import canonical_positions

# Atoms defining the RNA frame: one per nucleotide backbone position, stable
# across dialects and unaffected by base identity.
_FRAME_ATOMS = ("P", "C1'", "C4'")


@dataclass(frozen=True)
class Approach:
    centroid: np.ndarray  # (3,) binder centre of mass in the RNA frame
    axis: np.ndarray      # (3,) unit vector, binder's first principal axis


def _frame_atoms(rna):
    """Return (in-frame P/C1'/C4' atoms, their canonical positions).

    The positions are returned alongside rather than re-derived by the
    caller. Re-deriving them from the FILTERED atoms would re-run the
    sequence substring search on a reduced residue set, so a nucleotide that
    happened to be missing all three frame atoms would vanish from that
    derived sequence and raise an opaque FrameError from rnaframe instead of
    this module's own clear error.
    """
    positions = canonical_positions(rna)
    keep = (positions > 0) & np.isin(
        [str(n).strip() for n in rna.atom_name], _FRAME_ATOMS
    )
    return rna[keep], positions[keep]


def reference_frame(rna) -> tuple[np.ndarray, np.ndarray]:
    """Return (origin, rotation) defining the RNA's own coordinate system.

    Built by Gram-Schmidt on two vectors derived from the molecule itself,
    NOT from SVD principal axes. SVD was the original approach and it is
    subtly wrong here: LAPACK does not guarantee the sign of each singular
    vector is equivariant under rotation of the input, so rotating the whole
    complex could flip an axis and change the "invariant" answer. That defect
    stayed hidden because the test exercised one fixture; on three others the
    axis agreement was 0.85-0.93 instead of 1.0.

    The frame here is equivariant by construction: every basis vector is a
    difference of centroids of specific nucleotides, so rotating the input
    rotates the basis identically, and no sign choice is ever made.

        e1  origin -> loop (nucleotides 7-10)
        e2  the 5' half direction, orthogonalised against e1
        e3  e1 x e2, completing a right-handed frame

    Rows of `rotation` are e1/e2/e3, so `(x - origin) @ rotation.T` expresses
    any point in the RNA frame.
    """
    sel, positions = _frame_atoms(rna)
    coords = sel.coord.astype(np.float64)
    origin = coords.mean(axis=0)

    def centroid_of(nts):
        mask = np.isin(positions, nts)
        if not mask.any():
            raise ValueError(f"no frame atoms for nucleotides {nts}")
        return coords[mask].mean(axis=0)

    v1 = centroid_of([7, 8, 9, 10]) - origin      # toward the loop
    v2 = centroid_of([1, 2, 3, 4]) - origin       # toward the 5' end

    n1 = np.linalg.norm(v1)
    if n1 < 1e-8:
        raise ValueError("degenerate RNA frame: loop centroid coincides with origin")
    e1 = v1 / n1

    v2_perp = v2 - (v2 @ e1) * e1
    n2 = np.linalg.norm(v2_perp)
    if n2 < 1e-8:
        raise ValueError("degenerate RNA frame: 5' direction is parallel to e1")
    e2 = v2_perp / n2

    e3 = np.cross(e1, e2)
    return origin, np.stack([e1, e2, e3])


def approach(protein, rna) -> Approach:
    """Binder centroid and principal axis, expressed in the RNA frame.

    Computed from CA atoms only. Three of the five design pipelines ship
    design structures WITHOUT sidechains, while every refold has them, so
    using all atoms would make the centroid depend on sidechain presence:
    two identical poses would score differently purely because of which
    pipeline wrote the file. CA is the one atom set every file carries.
    """
    centroid, rot = reference_frame(rna)
    coords = ca_trace(protein).astype(np.float64)
    if len(coords) == 0:
        raise ValueError("protein chain has no CA atoms")
    local = (coords - centroid) @ rot.T

    binder_centroid = local.mean(axis=0)
    _, _, vt = np.linalg.svd(local - binder_centroid, full_matrices=False)
    axis = vt[0]
    # Sign convention: point away from the RNA origin so the value is stable.
    if float(axis @ binder_centroid) < 0:
        axis = -axis
    return Approach(binder_centroid, axis / np.linalg.norm(axis))
