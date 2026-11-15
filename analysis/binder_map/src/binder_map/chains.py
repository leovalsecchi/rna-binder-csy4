"""Composition-based protein/RNA chain selection.

Chain IDs are NOT reliable across this campaign's producers:

    layer   stack                          protein chain
    design  RFD3+LiMPNN                    B   <- inverted, 53% of the layer
    design  BD1/BG+LiMPNN/HD+LiMPNN        A
    design  BG+BIF                         A   (different mmCIF dialect)
    af3     all                            A
    esm     all                            A

So the rule is the one rna_binder.analysis.structural.split_chains uses:
RNA is the chain whose residue names subset {A, C, G, U}; protein is the other.
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

RNA_RESNAMES = frozenset({"A", "C", "G", "U"})


class ChainSelectError(Exception):
    """No unambiguous protein/RNA split exists for this structure."""


def load_structure(path: str | Path):
    """Read a .cif or .cif.gz into an AtomArray (first model, no solvent)."""
    path = Path(path)
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    block = pdbx.CIFFile.read(io.StringIO(raw.decode()))
    atoms = pdbx.get_structure(block, model=1)
    return atoms[~struc.filter_solvent(atoms)]


def _chain_kinds(atoms) -> dict[str, str]:
    """Map chain_id -> 'rna' | 'protein' by residue-name composition."""
    kinds: dict[str, str] = {}
    for chain_id in set(atoms.chain_id):
        names = set(atoms.res_name[atoms.chain_id == chain_id])
        kinds[chain_id] = "rna" if names <= RNA_RESNAMES else "protein"
    return kinds


def protein_chain_id(atoms) -> str:
    """The chain ID holding the protein. Raises if not exactly one."""
    kinds = _chain_kinds(atoms)
    protein = sorted(c for c, k in kinds.items() if k == "protein")
    if len(protein) != 1:
        raise ChainSelectError(
            f"expected exactly one protein chain, got {protein} from {kinds}"
        )
    return protein[0]


def split_chains(atoms):
    """Return (protein, rna) AtomArrays, selected by composition."""
    kinds = _chain_kinds(atoms)
    protein_ids = [c for c, k in kinds.items() if k == "protein"]
    rna_ids = [c for c, k in kinds.items() if k == "rna"]
    if len(protein_ids) != 1 or len(rna_ids) != 1:
        raise ChainSelectError(
            f"need exactly one protein and one RNA chain; got "
            f"protein={sorted(protein_ids)} rna={sorted(rna_ids)}"
        )
    protein = atoms[atoms.chain_id == protein_ids[0]]
    rna = atoms[atoms.chain_id == rna_ids[0]]
    return protein, rna
