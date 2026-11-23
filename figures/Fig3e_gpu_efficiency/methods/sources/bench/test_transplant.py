#!/usr/bin/env python3
"""Check the true-backbone transplant on real campaign data.

Uses the 20 matched BoltzGen pairs in analyses/bif_backbone_probe/bundle/, where
each dir holds the true backbone (00_backbone_spec<D>.cif), the BG+BIF design
(01_BIF_design_*.cif, the re-prediction) and the BG+LiMPNN design
(04_LiMPNN_design_*.cif, which IS the true backbone atom-for-atom).

THE ASSERTION: rebuilding the BG+LiMPNN design from the spec, via the same
build_true() the bench uses, must reproduce the real LiMPNN file exactly. If it
does, the transplant, the chain handling and the spec->design mapping are right.

It also reports how far the BIF re-prediction sits from the true backbone, which
is the bug being fixed -- expected ~0.8 A after a binder fit but metres apart raw.

Run:  python3 test_transplant.py
"""
from __future__ import annotations

import glob
import os
import pathlib
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bif_true_designs import build_true, compare, load_any  # noqa: E402
from rescore_lib import prot_backbone, seq_one_letter, split_prot_rna  # noqa: E402

BUNDLE = pathlib.Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/"
                      "bif_backbone_probe/bundle")
TOL = 0.01  # A


def main() -> int:
    dirs = sorted(d for d in BUNDLE.iterdir() if d.is_dir())
    if not dirs:
        print(f"no pair dirs under {BUNDLE}")
        return 1

    n_ok = n_bad = n_skip = 0
    worst = 0.0
    bif_dev = []

    for d in dirs:
        spec = next(iter(glob.glob(str(d / "00_backbone_spec*.cif"))), None)
        limp = next(iter(glob.glob(str(d / "04_LiMPNN_design_*.cif"))), None)
        bif = next(iter(glob.glob(str(d / "01_BIF_design_*.cif"))), None)
        if not (spec and limp):
            n_skip += 1
            continue

        real = load_any(limp)
        rp, rr = split_prot_rna(real)
        seq = seq_one_letter(rp)
        cp = str(np.asarray(rp.chain_id)[0])
        cr = str(np.asarray(rr.chain_id)[0]) if len(rr) else "B"

        built = build_true(spec, seq, cp, cr)
        bp, _br = split_prot_rna(built)
        dev, bad, n = compare(prot_backbone(rp), bp)

        if dev > TOL or bad != 0:
            print(f"FAIL {d.name}: max_dev={dev:.4f} A  resname_mismatch={bad}/{n}")
            n_bad += 1
        else:
            n_ok += 1
            worst = max(worst, dev)

        # how wrong the stored BIF design is: raw CA offset, NO superposition at
        # all (not the RNA-aligned figure). Expected 7.5-25.6 A from the probe.
        if bif:
            try:
                bpp, _ = split_prot_rna(load_any(bif))
                sp, _ = split_prot_rna(load_any(spec))
                a, b = prot_backbone(bpp), prot_backbone(sp)
                if len(a) == len(b):
                    ca_a = a[a.atom_name == "CA"].coord
                    ca_b = b[b.atom_name == "CA"].coord
                    if len(ca_a) == len(ca_b) and len(ca_a):
                        bif_dev.append(float(np.sqrt(
                            ((ca_a - ca_b) ** 2).sum(axis=1).mean())))
            except Exception:  # noqa: BLE001 - diagnostic only
                pass

    print(f"\nLiMPNN reproduced from spec: {n_ok} ok, {n_bad} failed, "
          f"{n_skip} skipped; worst coord dev {worst:.5f} A (tol {TOL})")
    if bif_dev:
        bif_dev.sort()
        print(f"stored BIF design vs true backbone, RAW CA RMSD (no superposition) "
              f"(n={len(bif_dev)}): median {bif_dev[len(bif_dev)//2]:.2f} A, "
              f"range {bif_dev[0]:.2f}-{bif_dev[-1]:.2f} A  <-- the bug")

    if n_bad or n_ok == 0:
        print("\nRESULT: FAIL")
        return 1
    print("\nRESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
