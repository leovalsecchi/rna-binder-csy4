"""Minimal backbone extraction, hashing, and PDB serialisation.

Foldseek reads backbone geometry, so sidechains are irrelevant — which is
fortunate, because the three LiMPNN stacks ship design CIFs without them
(inverse folding is not packing). The coordinate hash therefore uses only
N/CA/C/O, so a design and its refold are hashed on comparable content.
"""
from __future__ import annotations

import hashlib

import numpy as np

BACKBONE_ATOMS = ("N", "CA", "C", "O", "CB")
HASH_ATOMS = ("N", "CA", "C", "O")


def backbone_atoms(protein):
    """Keep only backbone (+CB) atoms."""
    return protein[np.isin(protein.atom_name, BACKBONE_ATOMS)]


def coord_hash(protein) -> str:
    """MD5 over N/CA/C/O coordinates rounded to 3 dp, in residue order."""
    sel = protein[np.isin(protein.atom_name, HASH_ATOMS)]
    order = np.lexsort((
        _atom_rank(sel.atom_name),
        _residue_ordinal(sel),
    ))
    coords = np.round(sel.coord[order].astype(np.float64), 3)
    # -0.0 and 0.0 must hash identically
    coords = coords + 0.0
    return hashlib.md5(coords.tobytes()).hexdigest()


def _atom_rank(names) -> np.ndarray:
    rank = {name: i for i, name in enumerate(HASH_ATOMS)}
    return np.array([rank[str(n)] for n in names])


def _residue_ordinal(atoms) -> np.ndarray:
    keys = list(zip(atoms.chain_id, atoms.res_id, atoms.ins_code))
    order: dict[tuple, int] = {}
    out = np.empty(len(keys), dtype=int)
    for i, key in enumerate(keys):
        if key not in order:
            order[key] = len(order)
        out[i] = order[key]
    return out


def ca_trace(protein) -> np.ndarray:
    """(n_residues, 3) array of CA coordinates in residue order."""
    ca = protein[protein.atom_name == "CA"]
    return ca.coord[np.argsort(_residue_ordinal(ca))].astype(np.float32)


def write_backbone_pdb(protein, name: str) -> bytes:
    """Serialise the backbone as a minimal PDB.

    Written by hand rather than via biotite's writer because this runs 215k
    times and the output only needs to be readable by Foldseek.
    """
    sel = backbone_atoms(protein)
    ordinals = _residue_ordinal(sel)
    lines = []
    for serial, (i, atom) in enumerate(zip(ordinals, sel), start=1):
        x, y, z = atom.coord
        atom_name = str(atom.atom_name)
        # PDB columns: atom names of <4 chars are left-padded by one space
        an = f" {atom_name:<3s}" if len(atom_name) < 4 else atom_name
        lines.append(
            f"ATOM  {serial:>5d} {an}{'':1s}{str(atom.res_name):>3s} A"
            f"{int(i) + 1:>4d}{'':4s}{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}{'':10s}{str(atom.element):>2s}"
        )
    lines.append("END")
    return ("\n".join(lines) + "\n").encode()
