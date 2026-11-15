"""Biotite-based structural metrics for refold outputs.

Computes:
    n_interface_contacts        protein residues within INTERFACE_CUTOFF of RNA
    n_h_bonds_interface         protein-RNA interface H-bonds — NOT computed here;
                                filled by the PyRosetta pass (see note below)
    interface_bsa_approx        approximate buried surface area (Å²)
    ca_rmsd_vs_designed         Cα RMSD of refolded vs source-design backbone (Å)
    helix_fraction              fraction of binder residues in helices (annotate_sse, P-SEA)
    sheet_fraction              fraction in beta sheets (annotate_sse)
    loop_fraction               fraction in loops/other (annotate_sse)
    sse_helix_fraction          DSSP (pydssp) fraction of binder residues in H/G/I
    sse_strand_fraction         DSSP (pydssp) fraction in E
    radius_of_gyration          Cα radius of gyration (Å)
    rna_contact_nt_1..N         per-nucleotide protein↔RNA heavy-atom contact counts

All metrics are protein-A and RNA-B centric (chain A = protein, chain B = RNA), but
the protein/RNA chains are found by composition (split_chains), so chain labels need
not literally be A/B (refold CIFs sometimes flip them).
"""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb


INTERFACE_CUTOFF = 5.0   # Å
RNA_CONTACT_CUTOFF = 4.5   # Å, heavy-atom, for per-nucleotide contact counts
RNA_LEN = 16               # csy4 target RNA (ACUGCCGUAUAGGCAG); cols rna_contact_nt_1..RNA_LEN
# NOTE: n_h_bonds_interface is intentionally NOT computed in this module.
# AF3/Boltz CIFs are heavy-atom only and this container (boltzgen.sif) has no
# usable hydrogen-adder, so biotite's struc.hbond can't run and a heavy-atom
# distance criterion overcounts badly (~8x). The column is filled afterwards by
# code/scripts/hbonds_pyrosetta.py (run in pyrosetta.sif) using Rosetta's
# hydrogen-aware, angle-checked detection. See run_structural.sbatch.


def load_structure(path: Path):
    """Load a .cif, .cif.gz, or .pdb file via biotite. Returns the first model."""
    path = Path(path)
    suffix = "".join(path.suffixes)
    if suffix.endswith(".cif.gz"):
        with gzip.open(path, "rt") as f:
            content = f.read()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as tmp:
            tmp.write(content)
            tmp_path = tmp.name
        try:
            cif_file = pdbx.CIFFile.read(tmp_path)
            structure = pdbx.get_structure(cif_file, model=1)
        finally:
            Path(tmp_path).unlink(missing_ok=True)
    elif path.suffix == ".cif":
        cif_file = pdbx.CIFFile.read(path)
        structure = pdbx.get_structure(cif_file, model=1)
    elif path.suffix == ".pdb":
        pdb_file = pdb.PDBFile.read(path)
        structure = pdb.get_structure(pdb_file, model=1)
    else:
        raise ValueError(f"Unsupported file extension: {path}")
    return structure


def split_chains(atoms):
    """Split a structure into (protein chain A, RNA chain B). Returns (protein, rna).

    If chains aren't named exactly A/B, find them by composition.
    """
    chain_ids = np.unique(atoms.chain_id)
    protein, rna = None, None
    for cid in chain_ids:
        chain_atoms = atoms[atoms.chain_id == cid]
        # Detect RNA by residue names (single-letter A/C/G/U)
        unique_residues = set(np.unique(chain_atoms.res_name))
        rna_residues = unique_residues & {"A", "C", "G", "U"}
        if len(rna_residues) > 0 and rna_residues == unique_residues - {"H"}:
            rna = chain_atoms
        else:
            # Treat as protein if 3-letter residue names dominate
            protein = chain_atoms
    return protein, rna


def n_interface_contacts(atoms) -> Optional[int]:
    """Count protein residues with any atom within INTERFACE_CUTOFF of an RNA atom."""
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return None
    # Use biotite's CellList for fast nearby-atom queries
    cell_list = struc.CellList(rna, cell_size=INTERFACE_CUTOFF)
    # For each protein atom, check if there's an RNA atom within cutoff
    contact_mask = np.zeros(protein.array_length(), dtype=bool)
    for i in range(protein.array_length()):
        neighbors = cell_list.get_atoms(protein.coord[i], radius=INTERFACE_CUTOFF)
        if len(neighbors) > 0:
            contact_mask[i] = True
    # Reduce atoms-to-residues
    contact_residues = np.unique(protein.res_id[contact_mask])
    return int(len(contact_residues))


def n_h_bonds_interface(atoms) -> Optional[int]:
    """Always returns None — interface H-bonds are computed by the PyRosetta pass.

    The previous heavy-atom distance criterion (count N/O–N/O pairs within 3.5 Å,
    no angle, no hydrogens) overcounted by ~8x (e.g. ~40 vs a true 5 for an
    80-aa / 16-nt complex). biotite's struc.hbond needs explicit D–H…A triplets,
    which AF3/Boltz heavy-atom CIFs lack, and boltzgen.sif has no hydrogen-adder.
    So this column is filled afterwards by code/scripts/hbonds_pyrosetta.py (run
    in pyrosetta.sif). We return None here so the column exists (as NaN) for that
    pass to overwrite. `atoms` is unused.
    """
    return None



def interface_bsa_approx(atoms) -> Optional[float]:
    """Approximate interface buried surface area: SASA(complex) - SASA(parts isolated).

    Returns total BSA in Å² (sum over both partners' buried surface).
    """
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return None
    try:
        # SASA of the complex
        sasa_complex = struc.sasa(atoms, point_number=200)
        # SASA of protein alone
        sasa_p = struc.sasa(protein, point_number=200)
        # SASA of RNA alone
        sasa_r = struc.sasa(rna, point_number=200)
        # Buried = isolated - in-complex
        # We need protein indices in the complex
        prot_mask = np.isin(atoms.chain_id, np.unique(protein.chain_id))
        rna_mask = np.isin(atoms.chain_id, np.unique(rna.chain_id))
        buried_p = np.nansum(sasa_p) - np.nansum(sasa_complex[prot_mask])
        buried_r = np.nansum(sasa_r) - np.nansum(sasa_complex[rna_mask])
        return float(buried_p + buried_r)
    except Exception:
        return None


def secondary_structure_fractions(atoms) -> dict:
    """Compute fraction of binder residues in helix/sheet/loop using biotite's annotate_sse."""
    protein, _ = split_chains(atoms)
    if protein is None:
        return {"helix_fraction": None, "sheet_fraction": None, "loop_fraction": None}
    try:
        sse = struc.annotate_sse(protein)
    except Exception:
        return {"helix_fraction": None, "sheet_fraction": None, "loop_fraction": None}
    if len(sse) == 0:
        return {"helix_fraction": None, "sheet_fraction": None, "loop_fraction": None}
    # biotite returns 'a' (alpha), 'b' (beta), 'c' (coil)
    n = len(sse)
    return {
        "helix_fraction": float((sse == "a").sum()) / n,
        "sheet_fraction": float((sse == "b").sum()) / n,
        "loop_fraction": float((sse == "c").sum()) / n,
    }


def sse_dssp_fractions(atoms) -> dict:
    """DSSP-based SSE fractions of the binder (protein) chain via pydssp.

    pydssp is a pure-python DSSP (no `mkdssp` binary needed). It emits 3 states
    ('-' loop, 'H' helix, 'E' strand); we report:
        sse_helix_fraction   fraction of binder residues in H/G/I
        sse_strand_fraction  fraction in E
    NOTE: pydssp's single 'H' state already encompasses the helix assignment, so
    G/I are not separately distinguished (they fall under 'H' or loop per pydssp).
    Binder chain found via split_chains (composition), not a hardcoded chain id.
    """
    out = {"sse_helix_fraction": None, "sse_strand_fraction": None}
    protein, _ = split_chains(atoms)
    if protein is None:
        return out
    try:
        import sys as _sys
        _pylibs = str(Path(__file__).resolve().parents[3] / "pylibs")
        if _pylibs not in _sys.path:
            _sys.path.insert(0, _pylibs)
        import pydssp
        order = ["N", "CA", "C", "O"]
        bb = protein[np.isin(protein.atom_name, order)]
        coords = []
        for rid in np.unique(bb.res_id):
            res = bb[bb.res_id == rid]
            d = {nm: res.coord[res.atom_name == nm][0] for nm in order if (res.atom_name == nm).any()}
            if all(k in d for k in order):
                coords.append([d[k] for k in order])
        if len(coords) < 4:
            return out
        arr = np.array(coords, dtype=np.float32)   # (L, 4, 3): N, CA, C, O per residue
        ss = np.asarray(pydssp.assign(arr)).astype(str)
        n = len(ss)
        out["sse_helix_fraction"] = float(np.isin(ss, ["H", "G", "I"]).sum()) / n
        out["sse_strand_fraction"] = float((ss == "E").sum()) / n
    except Exception:
        pass
    return out


def rna_contact_counts(atoms) -> dict:
    """Per-nucleotide protein↔RNA heavy-atom contact counts.

    For each RNA residue (chain order), counts the number of (protein heavy atom,
    RNA-residue heavy atom) pairs within RNA_CONTACT_CUTOFF Å. Returns
    {rna_contact_nt_<i>: count} for positions 1..RNA_LEN. Protein and RNA chains
    are found by composition (split_chains), so chain labels need not be A/B.
    """
    cols = {f"rna_contact_nt_{i}": None for i in range(1, RNA_LEN + 1)}
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return cols
    p_heavy = protein[protein.element != "H"]
    r_heavy = rna[rna.element != "H"]
    if p_heavy.array_length() == 0 or r_heavy.array_length() == 0:
        return cols
    try:
        cell = struc.CellList(p_heavy, cell_size=RNA_CONTACT_CUTOFF)
    except Exception:
        return cols
    for pos, rid in enumerate(np.unique(r_heavy.res_id)[:RNA_LEN], start=1):
        res_atoms = r_heavy[r_heavy.res_id == rid]
        cnt = 0
        for coord in res_atoms.coord:
            cnt += len(cell.get_atoms(coord, radius=RNA_CONTACT_CUTOFF))
        cols[f"rna_contact_nt_{pos}"] = int(cnt)
    return cols


def radius_of_gyration_ca(atoms) -> Optional[float]:
    """Cα radius of gyration of the binder (protein chain) in Å."""
    protein, _ = split_chains(atoms)
    if protein is None:
        return None
    ca_mask = (protein.atom_name == "CA")
    if ca_mask.sum() < 3:
        return None
    coords = protein.coord[ca_mask]
    center = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1))))


def ca_rmsd(atoms_a, atoms_b) -> Optional[float]:
    """Cα RMSD after Kabsch superimposition. Both must have the same N Cα atoms."""
    if atoms_a is None or atoms_b is None:
        return None
    ca_a = atoms_a[(atoms_a.atom_name == "CA")]
    ca_b = atoms_b[(atoms_b.atom_name == "CA")]
    if len(ca_a) != len(ca_b) or len(ca_a) < 3:
        return None
    try:
        fitted, _ = struc.superimpose(ca_a, ca_b)
        return float(struc.rmsd(ca_a, fitted))
    except Exception:
        return None


def ca_rmsd_vs_designed(refold_path: Path, designed_path: Optional[Path]) -> Optional[float]:
    """Cα RMSD between refold structure's protein chain and the designed source's.

    Args:
        refold_path: refold .cif file
        designed_path: path to source design file. If None, returns None.
    """
    if designed_path is None or not designed_path.exists():
        return None
    try:
        refold = load_structure(refold_path)
        designed = load_structure(designed_path)
        refold_prot, _ = split_chains(refold)
        designed_prot, _ = split_chains(designed)
        return ca_rmsd(refold_prot, designed_prot)
    except Exception:
        return None


def protein_sequence(atoms) -> Optional[str]:
    """1-letter amino-acid sequence of the protein (binder) chain, in residue order.

    The refold CIF's protein chain carries the designed binder sequence, so this
    is the sequence that was folded. Non-standard residues map to 'X'.
    """
    from biotite.sequence import ProteinSequence
    protein, _ = split_chains(atoms)
    if protein is None:
        return None
    try:
        _res_ids, res_names = struc.get_residues(protein)
    except Exception:
        return None
    letters = []
    for rn in res_names:
        try:
            letters.append(ProteinSequence.convert_letter_3to1(rn))
        except (KeyError, ValueError):
            letters.append("X")
    return "".join(letters) or None


def compute_all_structural(refold_path: Path, designed_path: Optional[Path] = None) -> dict:
    """Compute all structural metrics for one refold file.

    Args:
        refold_path: path to refold .cif
        designed_path: path to source-designed file for RMSD (optional)

    Returns:
        Flat metric dict with column-friendly names.
    """
    try:
        atoms = load_structure(refold_path)
    except Exception as e:
        return {"structural_error": str(e)}
    out = {
        "n_interface_contacts": n_interface_contacts(atoms),
        "n_h_bonds_interface": n_h_bonds_interface(atoms),
        "interface_bsa_approx": interface_bsa_approx(atoms),
        "radius_of_gyration_ca": radius_of_gyration_ca(atoms),
    }
    out.update(secondary_structure_fractions(atoms))
    out.update(sse_dssp_fractions(atoms))
    out.update(rna_contact_counts(atoms))
    out["ca_rmsd_vs_designed"] = ca_rmsd_vs_designed(refold_path, designed_path)
    out["binder_sequence"] = protein_sequence(atoms)
    return out
