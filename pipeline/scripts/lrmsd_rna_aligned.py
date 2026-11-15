#!/usr/bin/env python3
"""lrmsd_rna_aligned.py — LRMSD in the BINDER-CAMPAIGN frame, for a chunk of designs.

DockQ assigns receptor = larger chain, so for an ~80-residue binder vs a 16-nt target the
protein is the receptor: dockq_LRMSD superimposes on the BINDER and reports the RNA's
displacement. Holding the target fixed instead (superimpose on the RNA, measure the binder)
is the quantity a binder campaign actually cares about, and it runs ~1.45x larger.

Per design x backend this emits BOTH directions vs the designed backbone:

  lrmsd_rna_<be>       superimpose on RNA backbone     -> RMSD of protein backbone  (new)
  lrmsd_prot_<be>      superimpose on protein backbone -> RMSD of RNA backbone      (= DockQ's)

lrmsd_prot_* is a free self-check: it should reproduce the stored dockq_LRMSD_<be>.

Residue correspondence is recovered by NAME agreement, not by residue number -- design and
refold do not share a numbering origin (see _match).

Only LRMSD is recomputed — DockQ's fnat and iRMSD terms are frame-independent (fnat is
contact counting; iRMSD superimposes on BOTH chains' interface residues), so the rescored
DockQ is recombined later in pandas from the stored columns.

Work items are ARCHIVE TRIPLES so each packed unit is expanded once per chunk. Runs inside
boltzgen.sif; PYTHONPATH must include code/src and code/scripts.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import biotite.structure as struc

from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains

csv.field_size_limit(10 ** 7)

BACKENDS = {"af3": "af3_msaoff_cif_path", "b2": "b2_cif_path", "esm": "esm_cif_path"}
REF_COL = "designed_cif_path"
# DockQ constants.py BACKBONE_ATOMS, split by polymer type
PROT_BB = ("CA", "C", "N", "O")
RNA_BB = ("P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'",
          "C1'", "C2'", "C3'", "C4'", "C5'")

COLUMNS = ([f"lrmsd_rna_{b}" for b in BACKENDS] +
           [f"lrmsd_prot_{b}" for b in BACKENDS] + ["lrmsd_error"])


def _bb(atoms, names):
    return None if atoms is None else atoms[np.isin(atoms.atom_name, list(names))]


# Design and refold do NOT share a residue-numbering origin. The target PDB chain B is
# residues 4-21 (4 and 21 are atom stubs) and each generator resolved that differently:
# RFD3 renumbers 1-16, BoltzGen keeps 2-17, BoltzDesign/HD emit a complete 18-mer as 1-18 --
# while every refold is folded from a bare 16-mer and comes out 1-16. Keying atoms on
# (res_id, atom_name) therefore paired nucleotide i with nucleotide i for all four non-RFD3
# stacks when the same base is at i-1, sliding the fit one rung along a quasi-periodic
# A-form helix: it converges (~4.9 A instead of ~1.3 A) and silently returns a frame that is
# off by roughly one base-pair screw, inflating every LRMSD measured in it.
# So: recover the integer offset by residue-NAME agreement before matching.
_OFFSET_SEARCH = range(-6, 7)


def _index(a):
    """{(res_id, atom_name) -> row}, built ONCE per chain and probed at every offset."""
    pos = {}
    for i, (r, n) in enumerate(zip(a.res_id, a.atom_name)):
        pos.setdefault((int(r), str(n)), i)
    return pos


def _match_at(a, b, shift, pos=None):
    """Pair b to a with a's numbering shifted by `shift` (probe the unshifted index)."""
    pos = _index(a) if pos is None else pos
    ia, ib = [], []
    for j, (r, n) in enumerate(zip(b.res_id, b.atom_name)):
        i = pos.get((int(r) - shift, str(n)))
        if i is not None:
            ia.append(i)
            ib.append(j)
    return np.asarray(ia, int), np.asarray(ib, int)


def _best_offset(a, b, pos=None):
    """Offset to add to a's res_id so the two chains register. Scores candidates by the
    number of paired atoms whose residue NAMES agree; ties resolve toward 0 (the
    already-correct case, e.g. protein chains and RFD3's RNA)."""
    pos = _index(a) if pos is None else pos
    best, best_score = 0, -1
    for shift in _OFFSET_SEARCH:
        ia, ib = _match_at(a, b, shift, pos)
        if len(ia) < 3:
            continue
        score = int(np.sum(a.res_name[ia] == b.res_name[ib]))
        if score > best_score or (score == best_score and abs(shift) < abs(best)):
            best, best_score = shift, score
    return best, best_score


# Below this, the two chains are not the same polymer in the same register and any
# superposition built from them defines a meaningless frame.
MIN_REGISTER = 0.9


def _match(a, b):
    """Registered atom pairing. Third value is the fraction of paired atoms whose residue
    names agree -- 1.0 for a correct registration, ~0.12 for the old off-by-one."""
    pos = _index(a)
    shift, n_same = _best_offset(a, b, pos)
    ia, ib = _match_at(a, b, shift, pos)
    frac = (n_same / len(ia)) if len(ia) else 0.0
    return ia, ib, frac


def _apply(tr, atoms):
    if hasattr(tr, "apply"):
        return tr.apply(atoms)
    return struc.superimpose_apply(atoms, tr)


def _direction(fix_ref, fix_mod, meas_ref, meas_mod):
    """Superimpose model on reference via fix_*, then RMSD of meas_* with NO second fit."""
    if any(x is None or len(x) < 3 for x in (fix_ref, fix_mod, meas_ref, meas_mod)):
        return None
    fa, fb, f_reg = _match(fix_ref, fix_mod)
    ma, mb, m_reg = _match(meas_ref, meas_mod)
    if len(fa) < 3 or len(ma) < 3:
        return None
    if f_reg < MIN_REGISTER or m_reg < MIN_REGISTER:
        return None   # degrade this direction only, exactly as a short match does
    _, tr = struc.superimpose(fix_ref[fa], fix_mod[fb])
    return float(struc.rmsd(meas_ref[ma], _apply(tr, meas_mod[mb])))


def _parts(path):
    a = load_structure(path)
    a = a[a.element != "H"]
    p, r = split_chains(a)
    return _bb(p, PROT_BB), _bb(r, RNA_BB)


def _score_group(item):
    """One archive group: the designed ref is loaded ONCE per design, shared by the 3 backends."""
    _key, rows = item
    out = []
    for did, ref, refs in rows:
        rec = {}
        try:
            rp, rr = _parts(resolve(ref))
            for be, model_ref in refs.items():
                if not model_ref:
                    continue
                mp, mr = _parts(resolve(model_ref))
                rec[f"lrmsd_rna_{be}"] = _direction(rr, mr, rp, mp)    # target fixed
                rec[f"lrmsd_prot_{be}"] = _direction(rp, mp, rr, mr)   # DockQ's convention
        except Exception as e:
            rec["lrmsd_error"] = f"{type(e).__name__}: {e}"
        out.append((did, rec))
    return out


def _archive(ref):
    return ref.split(".tar.zst!")[0] if ".tar.zst!" in ref else ref


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    with open(a.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if a.limit:
        rows = rows[:a.limit]

    groups = defaultdict(list)
    n_skip = 0
    for r in rows:
        ref = (r.get(REF_COL) or "").strip()
        refs = {be: (r.get(col) or "").strip() for be, col in BACKENDS.items()}
        if not ref or not any(refs.values()):
            n_skip += 1
            continue
        key = (_archive(ref),) + tuple(_archive(refs[b]) for b in BACKENDS)
        groups[key].append((r["design_id"], ref, refs))
    items = sorted(groups.items(), key=lambda kv: kv[0])
    print(f"{len(rows)} rows, {len(items)} archive groups, {n_skip} skipped, jobs={a.jobs}",
          file=sys.stderr, flush=True)

    scored, n_err = {}, 0
    with Pool(a.jobs) as pool:
        for res in pool.imap_unordered(_score_group, items):
            for did, rec in res:
                scored[did] = rec
                n_err += int("lrmsd_error" in rec)

    tmp = a.output + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["design_id"] + COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            did = r["design_id"]
            if did in scored:
                w.writerow({"design_id": did, **scored[did]})
    os.replace(tmp, a.output)
    print(f"wrote {a.output}: {len(scored)} designs, {n_err} errors",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
