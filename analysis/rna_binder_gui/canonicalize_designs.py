#!/usr/bin/env python3
"""Canonicalize a design-stack designs/ dir (shared post-processing step).

Enforces the unified output convention for every design stack:
  * designs/ contains ONLY canonical design structures, as uncompressed .cif
  * .pdb design files are converted to .cif (then the .pdb is removed)
  * everything else (intermediate subdirs, .fa/.csv/.npz sidecars, raw files)
    is swept into <work_dir> so designs/ is exactly the design set

"Canonical" = a structure file whose name carries the design grammar token
(`_design` or `_recycle`). Directories are always intermediates → swept.

Idempotent and non-destructive: files are MOVED (never deleted) except the
source .pdb after a successful .cif conversion. Runs in any env with gemmi
(e.g. boltzgen.sif).

Usage:
    canonicalize_designs.py DESIGNS_DIR WORK_DIR
"""
import shutil
import sys
from pathlib import Path

import gemmi

_STRUCT_EXT = (".pdb", ".cif", ".cif.gz")


def _is_canonical(name: str) -> bool:
    return "_design" in name or "_recycle" in name


def _pdb_to_cif(pdb: Path) -> Path:
    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    out = pdb.with_suffix(".cif")
    st.make_mmcif_document().write_file(str(out))
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: canonicalize_designs.py DESIGNS_DIR WORK_DIR")
    designs = Path(sys.argv[1])
    work = Path(sys.argv[2])
    if not designs.is_dir():
        sys.exit(f"designs dir not found: {designs}")
    work.mkdir(parents=True, exist_ok=True)

    n_conv = n_keep = n_moved = n_fail = 0
    for entry in sorted(designs.iterdir()):
        name = entry.name
        # Directories are always intermediates.
        if entry.is_dir():
            shutil.move(str(entry), str(work / name))
            n_moved += 1
            continue
        # Canonical structure files stay (as .cif).
        if name.endswith(_STRUCT_EXT) and _is_canonical(name):
            if name.endswith(".pdb"):
                try:
                    _pdb_to_cif(entry)
                    entry.unlink()
                    n_conv += 1
                except Exception as e:
                    print(f"  WARN: pdb->cif failed for {name}: {e}", file=sys.stderr)
                    n_fail += 1
            else:  # .cif / .cif.gz already canonical
                n_keep += 1
            continue
        # Everything else is a sidecar/intermediate.
        shutil.move(str(entry), str(work / name))
        n_moved += 1

    print(f"canonicalize: {n_conv} pdb->cif, {n_keep} cif kept, "
          f"{n_moved} swept to {work}/, {n_fail} failed")


if __name__ == "__main__":
    main()
