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
    canonicalize_designs.py DESIGNS_DIR WORK_DIR [MATCH]

MATCH (optional): only process entries whose name contains this substring. Used
by multi-run launches (--num-runs > 1), where N concurrent jobs share one
designs/ dir: each job passes its own "_<uuid8>_" token so it only touches the
files it wrote. Without it, sibling jobs convert/move each other's files
mid-iteration -> FileNotFoundError, the jobs exit non-zero, and the shared
designs/ is left inconsistent.
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
    if len(sys.argv) not in (3, 4):
        sys.exit("usage: canonicalize_designs.py DESIGNS_DIR WORK_DIR [MATCH]")
    designs = Path(sys.argv[1])
    work = Path(sys.argv[2])
    match = sys.argv[3] if len(sys.argv) == 4 else None
    if not designs.is_dir():
        sys.exit(f"designs dir not found: {designs}")
    work.mkdir(parents=True, exist_ok=True)

    n_conv = n_keep = n_moved = n_fail = 0
    # Scoped (multi-run) mode: glob only this job's own files so a sibling job
    # churning the shared designs/ dir never even enters our listing. iterdir()
    # over the whole shared dir during concurrent create/unlink is the race that
    # made this step transiently fail under a ~60-way launch. Unscoped: walk all.
    entries = sorted(designs.glob(f"*{match}*")) if match is not None else sorted(designs.iterdir())
    for entry in entries:
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
                    print(f"  WARN: pdb->cif failed for {name}: {e} "
                          f"(quarantining to {work}/)", file=sys.stderr)
                    try:  # [Fix 3b] sweep the un-convertible .pdb aside so it is
                        # not packed into the .cif-only slice and silently dropped
                        # by .cif-only readers; n_fail surfaces the count.
                        shutil.move(str(entry), str(work / name))
                    except Exception as e2:
                        print(f"  WARN: quarantine failed for {name}: {e2}", file=sys.stderr)
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
