#!/usr/bin/env python3
"""export_deisolated_pairs.py — CIFs for backbones that were SINGLETONS at k=15 but joined a
cluster at k=30, each paired with the clustermate it attached to most strongly, superposed on
the RNA so the binder placement is directly comparable.

Purpose: eyeball whether the k=30 inclusion was justified — i.e. is the ex-singleton really the
same fold as its new clustermates, or did a laxer neighbour cap glue together things that are
not alike?

Sampling is STRATIFIED over the attachment strength (best TM to a clustermate), because that is
the whole question: 7% attach at TM<0.55 (barely over the floor) and 28% above 0.80.

RNA superposition uses lrmsd_rna_aligned's registered matching, which handles the fact that
RFD3 numbers the target 1-16, BoltzGen 2-17, and BoltzDesign/HD emit an 18-mer — an unregistered
fit silently produces a meaningless frame (the DockQ register bug).
"""
import csv, os, sys, shutil
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "src"))
import biotite.structure as struc
from biotite.structure.io import save_structure
from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains
from lrmsd_rna_aligned import RNA_BB, _bb, _match, _apply, MIN_REGISTER

FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
OUT = os.path.join(FM, "deisolated_pairs")
BANDS = [(0.50, 0.55), (0.55, 0.65), (0.65, 0.75), (0.75, 0.85), (0.85, 1.01)]
PER_BAND = 4
SEED = 42


def full(path):
    a = load_structure(path)
    return a[a.element != "H"]


def main():
    dei = pd.read_csv(os.path.join(FM, "deisolated_k15_k30.csv"))
    paths = pd.read_csv(os.path.join(FM, "passing_v6_final.csv"), skiprows=2,
                        usecols=["design_id", "design_cif"], low_memory=False)
    p = dict(zip(paths.design_id, paths.design_cif))

    rng = np.random.default_rng(SEED)
    picks = []
    for lo, hi in BANDS:
        sub = dei[(dei.best_tm >= lo) & (dei.best_tm < hi)]
        sub = sub[sub.struct.isin(p) & sub.partner_struct.isin(p)]
        if len(sub) == 0:
            print(f"  band {lo}-{hi}: EMPTY"); continue
        take = sub.iloc[rng.choice(len(sub), min(PER_BAND, len(sub)), replace=False)]
        picks.append(take)
        print(f"  band {lo:.2f}-{hi:.2f}: {len(sub):4d} available -> sampled {len(take)}")
    picks = pd.concat(picks).sort_values("best_tm").reset_index(drop=True)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    man = []
    for n, r in picks.iterrows():
        tag = f"case{n:02d}_tm{r.best_tm:.3f}"
        d = os.path.join(OUT, tag); os.makedirs(d, exist_ok=True)
        try:
            A = full(resolve(p[r.struct]))          # the ex-singleton  (reference, unmoved)
            B = full(resolve(p[r.partner_struct]))  # its new clustermate (superposed onto A)
            ra, rb = _bb(split_chains(A)[1], RNA_BB), _bb(split_chains(B)[1], RNA_BB)
            ia, ib, reg = _match(ra, rb)
            if len(ia) < 3 or reg < MIN_REGISTER:
                man.append(dict(case=tag, best_tm=r.best_tm, status=f"SKIP register={reg:.2f}",
                                ex_singleton=r.struct, clustermate=r.partner_struct))
                print(f"  {tag}: register {reg:.2f} — skipped"); continue
            _, tr = struc.superimpose(ra[ia], rb[ib])
            B = _apply(tr, B)
            save_structure(os.path.join(d, f"A_exsingleton_{r.struct}.cif"), A)
            save_structure(os.path.join(d, f"B_clustermate_{r.partner_struct}.cif"), B)
            # residual RNA rmsd after the fit = quality of the shared frame
            rr = float(struc.rmsd(ra[ia], _bb(split_chains(B)[1], RNA_BB)[ib]))
            man.append(dict(case=tag, best_tm=round(float(r.best_tm), 4), status="ok",
                            rna_fit_rmsd=round(rr, 3), register=round(float(reg), 3),
                            ex_singleton=r.struct, clustermate=r.partner_struct,
                            k30_cluster=int(r.c30)))
            print(f"  {tag}: ok  RNA-fit {rr:.2f} A  register {reg:.2f}")
        except Exception as exc:                                    # noqa: BLE001
            man.append(dict(case=tag, best_tm=r.best_tm, status=f"ERROR {type(exc).__name__}: {exc}",
                            ex_singleton=r.struct, clustermate=r.partner_struct))
            print(f"  {tag}: ERROR {exc}")
    pd.DataFrame(man).to_csv(os.path.join(OUT, "manifest.csv"), index=False)
    ok = sum(1 for m in man if m["status"] == "ok")
    print(f"\n{ok}/{len(man)} pairs written -> {OUT}")


if __name__ == "__main__":
    main()
