#!/usr/bin/env python3
"""pi_stack_pernt.py — pi-stacking resolved PER RNA NUCLEOTIDE.

The stored pi_stack_* columns give a per-design total and `pi_stack_res_nt` is only a COUNT
of distinct bases, so they cannot answer "which nucleotide is being stacked". This emits a
per-nucleotide breakdown so a specific base can be excluded afterwards.

Why that matters here: the Csy4 target was CROPPED for design. In the refold CIFs it is 16 nt
(ACUGCCGUAUAGGCAG, ids 1-16) with stem pairs 2-16, 3-15, 4-14, 5-13, 6-12, 7-11 and loop
8-9-10. **nt 1 (A) is unpaired in 100% of refolds** — a dangling 5' base whose faces are bare
only because the hairpin was truncated. Stacking on it is a construct artefact, not
recognition of the target, so totals that include it overstate real pi-stacking.
(Designed backbones carry 18 nt, CACUGCCGUAUAGGCAGC, i.e. one extra C at each end; refold
index i == designed index i+1. All pi_stack numbers here are in REFOLD 16-nt numbering.)

Geometry, thresholds and weights are imported from rna_binder.analysis.pi_stack — this file
adds only the per-nucleotide bookkeeping, so it can never drift from the production metric.

Emits per design: pis_nt<i>_count / pis_nt<i>_score for i in 1..NT_MAX, plus the totals and
the same totals with nt 1 removed.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains
from rna_binder.analysis import pi_stack as PS

csv.field_size_limit(10 ** 7)
NT_MAX = 16
DANGLING = 1                     # the cropped, permanently unpaired 5' base


def _pernt(atoms):
    """Same pairing loop as pi_stack.pi_stack, but keyed by RNA residue id."""
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return {"pis_error": "split_chains"}
    aa_rings = PS._rings(protein, PS.AA_RING)
    nt_rings = PS._rings(rna, PS.BASE_RING)
    cnt = {i: 0 for i in range(1, NT_MAX + 1)}
    sco = {i: 0.0 for i in range(1, NT_MAX + 1)}
    for ak, arn, ac, an in aa_rings:
        for nk, nrn, nc, nn in nt_rings:
            v = nc - ac
            d = float(np.linalg.norm(v))
            if d > PS.DIST_MAX:
                continue
            a = PS._plane_angle(an, nn)
            axial = abs(float(np.dot(v, an)))
            lateral = float(np.sqrt(max(0.0, d * d - axial * axial)))
            if a <= PS.ANGLE_PAR and lateral <= PS.OFFSET_MAX:
                s = (float(np.exp(-((d - PS.D0) / PS.DW) ** 2))
                     * float(np.cos(np.radians(a)) ** 2)
                     * PS.W_AA.get(arn, 0.7) * PS.W_NT.get(nrn, 0.7))
                if s >= PS.PAIR_MIN:
                    rid = int(nk[1])
                    if rid in cnt:
                        cnt[rid] += 1
                        sco[rid] += s
    out = {}
    for i in range(1, NT_MAX + 1):
        out[f"pis_nt{i}_count"] = cnt[i]
        out[f"pis_nt{i}_score"] = round(sco[i], 4)
    out["pis_count_all"] = sum(cnt.values())
    out["pis_score_all"] = round(sum(sco.values()), 4)
    out["pis_count_nodangle"] = sum(v for k, v in cnt.items() if k != DANGLING)
    out["pis_score_nodangle"] = round(sum(v for k, v in sco.items() if k != DANGLING), 4)
    return out


COLUMNS = ([f"pis_nt{i}_{k}" for i in range(1, NT_MAX + 1) for k in ("count", "score")] +
           ["pis_count_all", "pis_score_all", "pis_count_nodangle", "pis_score_nodangle",
            "pis_error"])


def _one(ref):
    try:
        return ref, _pernt(load_structure(resolve(ref)))
    except Exception as e:
        return ref, {"pis_error": f"{type(e).__name__}: {e}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--path-col", default="af3_msaoff_cif_path")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()

    with open(a.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    refs = sorted({(r.get(a.path_col) or "").strip() for r in rows} - {""})
    print(f"{len(rows)} rows, {len(refs)} distinct structures, jobs={a.jobs}",
          file=sys.stderr, flush=True)

    got = {}
    with Pool(a.jobs) as pool:
        for ref, m in pool.imap_unordered(_one, refs, chunksize=8):
            got[ref] = m
    n_err = sum(1 for m in got.values() if m.get("pis_error"))

    tmp = a.output + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["design_id"] + COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            ref = (r.get(a.path_col) or "").strip()
            if ref in got:
                w.writerow({"design_id": r["design_id"], **got[ref]})
    os.replace(tmp, a.output)
    print(f"wrote {a.output}: {len(got)} structures, {n_err} errors",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
