#!/usr/bin/env python3
"""Rewrite a bench BG+BIF run's design files to carry the TRUE BoltzGen backbone.

THE BUG THIS FIXES
------------------
boltzgen_multif.sh:209-246 harvests BG+BIF's "design" from BoltzGen's
`refold_cif/`, which is a Boltz-2 RE-PREDICTION of the designed sequence, not
the backbone the design step produced. The fold is right (0.82 A median) but the
POSE is not: median 4.4 A RNA-aligned, 25% of designs >10 A from where BoltzGen
actually put the binder. Every BG+BIF ligand_RMSD / sc_DockQ in the campaign was
therefore measured against a different docking than the design step produced.

THE FIX
-------
BG+BIF and BG+LiMPNN share a uuid8 AND the same `spec_<D>.cif` backbones. So the
true design file is: spec_<D> protein backbone (N,CA,C,O) + the RNA, with residue
names set positionally from the BIF sequence.

The spec CANNOT be used directly -- it carries BoltzGen's OWN co-designed
sequence, only ~50% identical to the BIF sequence, which trips
lrmsd_rna_aligned.py's MIN_REGISTER=0.9 gate and makes DockQ's internal
alignment a coin flip. Hence the sequence transplant, reusing the already
validated rescore_lib.transplant() rather than a second implementation.

Backbone-only is correct, not a shortcut: BG+LiMPNN design files are themselves
backbone-only (LigandMPNN writes no sidechains without --pack_side_chains), so
this makes the two inverse folders SYMMETRIC -- which is the whole point of a
stack comparison. ligand_RMSD is backbone-only anyway; only fnat moves on stripping.

This does NOT touch generation. `--steps inverse_folding folding` stays exactly
as production, so the BoltzGen job's GPU cost is unchanged and still comparable
to the campaign. The refolders read only the SEQUENCE, so their cost is
unchanged too. Only the stored coordinates -- the scoring reference -- change.

THE SELF-CHECK
--------------
--validate performs the same transplant with the BG+LiMPNN sequence and compares
it to the REAL BG+LiMPNN design file. Those must agree atom-for-atom, because a
BG+LiMPNN design file IS the spec backbone. Any drift means the transplant, the
spec->design index mapping, or the chain handling is wrong. Run it before
trusting any rewritten design.

WHERE THE REWRITE HAPPENS (corrected 2026-11-11)
------------------------------------------------
Not on Lustre. A finished gen run has no `designs/` -- the finalize tail packs it
into `designs.tar.zst` and deletes the tree (see bench_lib.py). Rewriting in
place would mean unpacking and repacking a production gen artefact, which is
both destructive and pointless, since the only consumer is the refold pool. So
`rewrite_dir()` is called by bench_refold.py on the STAGED copy while a pool is
being built: production output is never touched and the re-prediction stays
exactly where the gen job left it.

Usage:
  python3 bif_true_designs.py --label gpubench_v1 --validate   # self-check only
  (the rewrite itself is driven by bench_refold.py's pool build)
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys
import tempfile

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure

RB = pathlib.Path("/cluster/work/wicky/lvalsecchi/rna_binder")
CAMP = RB / "campaigns" / "csy4_crrna_free"
BENCH = RB / "bench_gpueff"
RESCORE = RB / "analyses" / "bif_true_bb_rescore"

sys.path.insert(0, str(BENCH))
import bench_lib as BL  # noqa: E402  (reads the packed designs.tar.zst)

sys.path.insert(0, str(RESCORE))
from rescore_lib import (  # noqa: E402  (path must be set first)
    prot_backbone, seq_one_letter, split_prot_rna, transplant,
)

COORD_TOL = 0.01  # A. The validated recipe reproduces LiMPNN to 0.0005 A.
DESIGN_IDX = re.compile(r"_design(\d+)_seq(\d+)$")
# A pool mixes several uuids, so the donor spec must be chosen per FILE, not per
# run dir: csy4free_BG_BIF_len80_<uuid8>_design<D>_seq<S>.cif
NAME_RE = re.compile(r"_len\d+_([0-9a-f]{8})_design(\d+)_seq(\d+)$")

_SPEC_CACHE: dict[str, dict[int, str]] = {}


def specs_for(uuid8: str) -> dict[int, str]:
    """spec_map with memoisation -- a 600-design pool would otherwise re-glob
    the same 60-file backbone dir 600 times."""
    if uuid8 not in _SPEC_CACHE:
        _SPEC_CACHE[uuid8] = spec_map(uuid8)
    return _SPEC_CACHE[uuid8]


def load_any(path: str):
    p = str(path)
    if p.endswith(".pdb"):
        return PDBFile.read(p).get_structure(model=1)
    return get_structure(CIFFile.read(p), model=1)


def spec_map(uuid8: str) -> dict[int, str]:
    """design index -> spec cif path. Padding is inconsistent across runs
    (spec_7 vs spec_07), so always match on the PARSED int, never the string."""
    out = {}
    for bb_dir in CAMP.glob(f"_backbones/BG_len80_{uuid8}"):
        for f in bb_dir.glob("spec_*.cif"):
            m = re.match(r"spec_(\d+)$", f.stem)
            if m:
                out[int(m.group(1))] = str(f)
    return out


def build_true(spec_path: str, seq: str, chain_prot: str, chain_rna: str):
    """spec backbone + RNA, residue names transplanted from `seq`, chains stamped
    to match whatever the original design file used (the downstream pipeline
    already keys on those, and the convention differs between stacks)."""
    sp, sr = split_prot_rna(load_any(spec_path))
    bb = prot_backbone(sp)
    if len(bb) == 0 or len(sr) == 0:
        raise ValueError("spec split left an empty protein or RNA side")
    prot = transplant(bb, seq).copy()
    rna = sr.copy()
    prot.chain_id[:] = chain_prot
    rna.chain_id[:] = chain_rna
    st = prot + rna
    # set_annotation, NOT `st.b_factor = ...`. biotite's __setattr__ only routes
    # to an annotation that ALREADY exists; a new name becomes a plain Python
    # attribute, so the CIF writer never sees it and the column is silently
    # omitted. That cost this bench a full refold round: BG+BIF's three refolds
    # died in stage 1 on `KeyError: '_atom_site.B_iso_or_equiv'`, because the
    # worker parses designs with BioPython's MMCIFParser, which hard-requires
    # that column (and occupancy). The four stacks whose designs were untouched
    # refolded fine, which is what pointed here.
    st.set_annotation("b_factor", np.zeros(len(st), dtype=float))
    st.set_annotation("occupancy", np.ones(len(st), dtype=float))
    st.set_annotation("charge", np.zeros(len(st), dtype=int))
    return st


def sort_key(st):
    """Stable atom order for comparison: (res_id, atom_name)."""
    return np.lexsort((np.asarray(st.atom_name), np.asarray(st.res_id).astype(int)))


def compare(a, b) -> tuple[float, int, int]:
    """(max coord deviation, resname mismatches, n_atoms) after sorting both."""
    if len(a) != len(b):
        return (float("inf"), -1, len(a))
    ia, ib = sort_key(a), sort_key(b)
    ca, cb = a.coord[ia], b.coord[ib]
    dev = float(np.abs(ca - cb).max())
    bad = int((np.asarray(a.res_name)[ia] != np.asarray(b.res_name)[ib]).sum())
    return (dev, bad, len(a))


def rewrite_one(path: pathlib.Path) -> None:
    """Replace one staged BG+BIF design file with the TRUE BoltzGen backbone,
    carrying its own BIF sequence. Raises on anything unexpected -- a silently
    skipped file would be scored against the re-prediction again, which is the
    entire bug this exists to remove."""
    m = NAME_RE.search(path.stem)
    if not m:
        raise ValueError(f"cannot parse uuid8/design index from {path.name}")
    uuid8, d = m.group(1), int(m.group(2))
    specs = specs_for(uuid8)
    if d not in specs:
        raise ValueError(f"no spec_{d} in _backbones/BG_len80_{uuid8} for {path.name}")
    op, orr = split_prot_rna(load_any(str(path)))
    seq = seq_one_letter(op)
    cp = str(np.asarray(op.chain_id)[0])
    cr = str(np.asarray(orr.chain_id)[0]) if len(orr) else "B"
    st = build_true(specs[d], seq, cp, cr)
    cf = CIFFile()
    set_structure(cf, st)
    cf.write(str(path))


def rewrite_dir(ddir: pathlib.Path) -> int:
    """Rewrite every BG+BIF design in a STAGED pool directory. Called by
    bench_refold.py during the pool build, never against Lustre."""
    files = sorted(f for f in ddir.iterdir() if f.is_file())
    print(f"=== BG+BIF true-backbone rewrite: {len(files)} staged design(s) ===")
    for f in files:
        rewrite_one(f)
    print(f"  rewrote {len(files)} onto the real spec backbone "
          f"({len(_SPEC_CACHE)} uuid(s) of specs used)")
    return len(files)


def validate(label: str, per_run: int = 20) -> bool:
    """Reproduce BG+LiMPNN designs from their spec. Must be exact, because a
    BG+LiMPNN design file IS the spec backbone -- any drift means the transplant,
    the spec->design index mapping or the chain handling is wrong.

    Reads from the packed designs.tar.zst, since a finished run has no designs/.
    `per_run` designs per uuid is plenty: the failure modes this catches (wrong
    spec index, wrong chain, residue-order mismatch) are systematic, not
    per-file, and the 2026-09-30 probe already covered 20 pairs end to end."""
    n_ok = n_bad = n_skip = 0
    worst = 0.0
    print("=== validation: rebuild BG+LiMPNN designs from spec ===")
    for rdir in BL.finished_runs(CAMP / "boltzgen__ligandmpnn", label):
        uuid8 = rdir.name.split("_")[-1]
        specs = specs_for(uuid8)
        if not specs:
            print(f"  {uuid8}: no _backbones/BG_len80_{uuid8} -- SKIP")
            n_skip += 1
            continue
        want = BL.members(rdir)[:per_run]
        with tempfile.TemporaryDirectory(prefix="bifval_") as td:
            BL.extract(rdir, want, pathlib.Path(td))
            for f in sorted((pathlib.Path(td) / "designs").iterdir()):
                m = DESIGN_IDX.search(f.stem)
                if not m:
                    continue
                d = int(m.group(1))
                if d not in specs:
                    n_skip += 1
                    continue
                try:
                    real = load_any(str(f))
                    rp, rr = split_prot_rna(real)
                    seq = seq_one_letter(rp)
                    cp = str(np.asarray(rp.chain_id)[0])
                    cr = str(np.asarray(rr.chain_id)[0]) if len(rr) else "B"
                    built = build_true(specs[d], seq, cp, cr)
                    bp, _ = split_prot_rna(built)
                    dev, bad, n = compare(prot_backbone(rp), bp)
                except Exception as e:  # noqa: BLE001 - report, never mask
                    print(f"  FAIL {f.name}: {type(e).__name__}: {e}")
                    n_bad += 1
                    continue
                worst = max(worst, dev if dev != float("inf") else worst)
                if dev > COORD_TOL or bad != 0:
                    print(f"  FAIL {f.name}: max_dev={dev:.4f} A  "
                          f"resname_mismatch={bad}/{n}")
                    n_bad += 1
                else:
                    n_ok += 1
    print(f"  reproduced {n_ok}, failed {n_bad}, skipped {n_skip}; "
          f"worst coord dev {worst:.4f} A (tol {COORD_TOL})")
    if n_ok == 0:
        print("  !! nothing validated -- refusing to call this a pass")
        return False
    return n_bad == 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="gpubench_v1")
    ap.add_argument("--validate", action="store_true",
                    help="rebuild BG+LiMPNN designs from spec and require exact agreement")
    ap.add_argument("--per-run", type=int, default=20)
    a = ap.parse_args()
    if not a.validate:
        ap.error("pass --validate (the rewrite itself is driven by bench_refold.py)")
    if not validate(a.label, a.per_run):
        sys.exit("validation FAILED -- pools must not be built")
    print("validation PASSED")


if __name__ == "__main__":
    main()
