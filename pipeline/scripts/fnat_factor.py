#!/usr/bin/env python3
"""fnat_factor.py — per-design correction factor for DockQ's fnat denominator leak.

DockQ computes nat_total on the UNTOUCHED native (the design) but nat_correct only over
residues aligned to the model (the refold, always the 16-mer ACUGCCGUAUAGGCAG). Designs whose
RNA carries extra nucleotides (BD1/HD: 18 nt) therefore have contacts in fnat's denominator
that the numerator can never match. Verified on 30 design x backend pairs per stack:
  - trimming the native leaves iRMSD and LRMSD IDENTICAL (100%)
  - nat_correct is UNCHANGED by trimming (100%)   <- only the denominator leaks
  - fnat_trim == fnat_full * (nat_full/nat_trim)  EXACTLY (100%)
so the fix is arithmetic and needs no DockQ rerun.

nat_total is counted WITHIN the native, so the factor depends on the DESIGN ONLY -- one number
per design, valid for all three backends. This reads one structure per design and computes two
residue-distance matrices; no refold is touched and no trimmed file is written (DockQ's
get_residue_distances takes a residue tuple, so the trim is just a filtered tuple).

Emits: design_id, n_rna_res, nat_full, nat_trim, fnat_factor
"""
import argparse, csv, os, sys
from collections import defaultdict
from multiprocessing import Pool
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/pylibs_dockq")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from DockQ.DockQ import load_PDB, get_residue_distances
from DockQ.constants import FNAT_THRESHOLD
from archive_read import resolve

csv.field_size_limit(10 ** 7)
CANON = "ACUGCCGUAUAGGCAG"
REF_COL = "designed_cif_path"
COLUMNS = ["n_rna_res", "nat_full", "nat_trim", "fnat_factor", "fnat_factor_error"]
NUC = {"A", "C", "G", "U"}


def _letter(resname):
    r = resname.strip()
    return r[-1] if r[-1] in NUC else None


def _split(chains):
    """(protein residues, rna residues) by composition, not chain id."""
    prot, rna = None, None
    for c in chains:
        res = tuple(c.get_residues())
        if not res:
            continue
        lets = [_letter(r.get_resname()) for r in res]
        frac = sum(1 for x in lets if x) / len(lets)
        if frac > 0.8:
            if rna is None or len(res) > len(rna):
                rna = res
        else:
            if prot is None or len(res) > len(prot):
                prot = res
    return prot, rna


def _count(prot, rna):
    d = get_residue_distances(tuple(prot), tuple(rna), "ref")
    return int(np.nonzero(np.asarray(d) < FNAT_THRESHOLD ** 2)[0].shape[0])


def _one(item):
    did, ref = item
    try:
        p = resolve(ref)
        prot, rna = _split(load_PDB(p))
        if prot is None or rna is None:
            return did, {"fnat_factor_error": "missing chain"}
        seq = "".join(_letter(r.get_resname()) or "X" for r in rna)
        i = seq.find(CANON)
        nat_full = _count(prot, rna)
        if i < 0:
            return did, {"n_rna_res": len(rna), "nat_full": nat_full,
                         "fnat_factor_error": f"canonical 16-mer not found in {seq}"}
        trim = rna[i:i + len(CANON)]
        nat_trim = nat_full if len(trim) == len(rna) else _count(prot, trim)
        return did, {"n_rna_res": len(rna), "nat_full": nat_full, "nat_trim": nat_trim,
                     "fnat_factor": round(nat_full / nat_trim, 6) if nat_trim else ""}
    except Exception as e:
        return did, {"fnat_factor_error": f"{type(e).__name__}: {e}"}


def _archive(ref):
    return ref.split(".tar.zst!")[0] if ".tar.zst!" in ref else ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()
    with open(a.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    # group by archive so each packed unit is expanded once
    groups = defaultdict(list)
    for r in rows:
        ref = (r.get(REF_COL) or "").strip()
        if ref:
            groups[_archive(ref)].append((r["design_id"], ref))
    items = [it for _k, v in sorted(groups.items()) for it in v]
    print(f"{len(rows)} rows, {len(items)} with a design ref, "
          f"{len(groups)} archives, jobs={a.jobs}", file=sys.stderr, flush=True)
    scored = {}
    with Pool(a.jobs) as pool:
        for did, rec in pool.imap(_one, items, chunksize=16):
            scored[did] = rec
    tmp = a.output + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["design_id"] + COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            did = r["design_id"]
            if did in scored:
                w.writerow({"design_id": did, **scored[did]})
    os.replace(tmp, a.output)
    nerr = sum(1 for v in scored.values() if "fnat_factor_error" in v)
    print(f"wrote {a.output}: {len(scored)} designs, {nerr} errors", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
