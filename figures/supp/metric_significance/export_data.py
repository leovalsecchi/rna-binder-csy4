#!/usr/bin/env python
"""The two per-re-predictor tables panel (a) needs, for all 12,000 ordered designs.

Emits, into `inputs/`:

1. `if_interface_aa_12k.csv.gz` — interface amino-acid composition, computed from
   structures (see below). ~6 min; **skipped if the file already exists**, because nothing
   upstream of it changes.
2. `metric_significance_backends_12k.csv.gz` — every other metric, per re-predictor,
   pulled out of the local 12k master and **corrected for the BG+BIF true-backbone
   rescore**. Seconds.

`make_figure.py` reads only `inputs/`; nothing outside the repo is touched at figure time.

### 1. Interface composition

`if_N`, `if_arom` and `if_neg` are NOT columns anywhere: they are computed from a
structure. `code/Fig5def_asn_interface/export_data.py` derives `if_N` from an interface
mask, and this is that computation, unchanged in substance, run on BOTH the AF3 and the
ESMFold2 model of every ordered design instead of on the selected one. Same cutoff (4.5 A
heavy-atom, protein residue to any RNA atom), same `geom.split` / `geom.AA3`, same
"fraction of interface residues" definition - so `if_N` here is comparable with the
numbers in `inputs/asn/asn32.csv`.

Why both backends: panel (a) now scores every metric on ONE model per design, the one with
the lower ligand_RMSD, so each metric has to exist from both or the design has no comparable
value. Boltz-2 is not computed - it was never gated on and it fails the both-backends rule
for several of the other metrics.

The CIFs are local: /Users/leo/Downloads/design29_742db4d3/an32/cif/<design_id>.<be>.cif,
24,032 files = 12,016 designs x 2 re-predictors, the same tree Fig 5d-f read. ~13 ms per
structure, so ~6 min single-core for the 24,000, ~1 min on 8 workers.

VERIFIED against the repo's own export: `af3_if_arom` / `esm_if_arom` here reproduce
`inputs/hits32_geom_12k.csv.gz` to max |d| = 5e-5 over all 12,000 (r = 1.00000), and the
interface SIZE `n_if` matches its `n_iface` exactly, design for design. So the interface
mask is the same one the rest of the book uses; only `if_N` and `if_neg` are new columns.

### 2. The per-re-predictor metric table

From `/Users/leo/Downloads/design29_742db4d3/an32/m32.pkl`, the 12,000 x 281 local master
that Fig 5d-f and Fig 5f-g also read. Its per-backend gate columns are IDENTICAL to
`inputs/large/gated_pool_v6_gate_metrics.csv.gz` on all 12,000 (`dockq_LRMSD_* == rmsd_*`
and `dockq_self_* == dockq_*`, max |d| = 0.0000), so the two sources agree and the gate
columns are not re-derived here.

⚠️ **BG+BIF is corrected.** The stored BG+BIF ligand_RMSD / DockQ are wrong (the route's design
CIF is a re-prediction, so the stored scores are not against the true backbone); the
rescue columns are in `inputs/large/bif_true_per_backend.csv.gz`. They are substituted
for the 2,397 of 2,400 BG+BIF designs that carry them. This is not cosmetic: because the
re-predictor is CHOSEN by minimising ligand_RMSD, the correction **changes which model is
chosen for 794 of the 12,000 designs**, all inside BG+BIF - a third of that route.

Boltz-2 columns are dropped throughout. B2 was never a gate backend, and the rule is that
a metric must exist from both AF3 and ESMFold2 or it is not scored.

    ~/.venvs/binder_map/bin/python3 code/supp/metric_significance/export_data.py
"""
import sys
from pathlib import Path
from multiprocessing import Pool

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
AN = "/Users/leo/Downloads/design29_742db4d3/an32"
sys.path.insert(0, AN)

import numpy as np
import pandas as pd
import biotite.structure as struc
import geom

from common import paths

CUT = 4.5                       # A, as Fig5def_asn_interface and geom.py both use
AROM = list("FWY")              # Phe/Trp/Tyr; His is NOT counted, matching if_arom in
                                # code/Fig5fg_h1null_refolder/export_data.py
ACID = list("DE")               # Asp/Glu, the if_neg set


def interface(path):
    """L, interface size, and the Asn / aromatic / acidic fraction OF THE INTERFACE."""
    try:
        a = geom.load(path)
    except Exception:
        return None
    prot, rna = geom.split(a)
    if prot is None or rna is None or len(prot) == 0:
        return None
    pres = [(r, prot[prot.res_id == r]) for r in np.unique(prot.res_id)]
    s = np.array([geom.AA3.get(x.res_name[0], "X") for _, x in pres])
    cell = struc.CellList(rna, cell_size=CUT)
    hit = np.array([sum(len(cell.get_atoms(c, radius=CUT)) for c in x.coord) > 0
                    for _, x in pres])
    if not hit.any():           # no interface: every fraction would be 0/0
        return None
    si = s[hit]
    return dict(L=len(s), n_if=int(hit.sum()),
                if_N=float((si == "N").mean()),
                if_arom=float(np.isin(si, AROM).mean()),
                if_neg=float(np.isin(si, ACID).mean()))


def one(did):
    row = {"design_id": did}
    for be in ("af3", "esm"):
        r = interface(f"{AN}/cif/{did}.{be}.cif") or {}
        row.update({f"{be}_{k}": v for k, v in r.items()})
    return row


# The per-re-predictor columns panel (a) scores, as {output stem: (af3 column, esm column)}
# in m32.pkl's spelling. `af3_msaoff_iptm` is AF3's MSA-off ipTM, which is what
# code/Fig5fg_h1null_refolder/export_data.py uses and what is comparable with ESMFold2.
PAIRS = {
    "iptm": ("af3_msaoff_iptm", "esm_iptm"),
    "packstat": ("packstat_af3", "packstat_esm"),
    "sc_value": ("sc_value_af3", "sc_value_esm"),
    "buns_vsasa_opth": ("buns_vsasa_opth_af3", "buns_vsasa_opth_esm"),
    "n_h_bonds_interface": ("n_h_bonds_interface_af3", "n_h_bonds_interface_esm"),
    "pi_stack_score": ("pi_stack_score_af3", "pi_stack_score_esm"),
    "dsasa_total": ("dsasa_total_af3", "dsasa_total_esm"),
    "dsasa_polar": ("dsasa_polar_af3", "dsasa_polar_esm"),
}
MASTER = f"{AN}/m32.pkl"


def interface_table(lib):
    dest = paths.INPUTS / "if_interface_aa_12k.csv.gz"
    if dest.exists():
        print(f"{dest.name} exists, not recomputed")
        return
    with Pool(8) as p:
        rows = p.map(one, lib.design_id.tolist(), chunksize=32)
    out = pd.DataFrame(rows)
    out.to_csv(dest, index=False)
    for be in ("af3", "esm"):
        ok = out[f"{be}_if_N"].notna()
        print(f"{be}: {ok.sum():,}/{len(out):,} with an interface; "
              f"if_N {out[f'{be}_if_N'][ok].mean():.4f}  "
              f"if_arom {out[f'{be}_if_arom'][ok].mean():.4f}  "
              f"if_neg {out[f'{be}_if_neg'][ok].mean():.4f}")
    print(f"wrote {dest}")


def backend_table(lib):
    m = pd.read_pickle(MASTER)
    assert len(m) == 12000 and m.design_id.nunique() == 12000, f"{len(m)} rows in m32.pkl"

    # The gate's own per-backend columns, cross-checked against the gated-pool export.
    gp = pd.read_csv(paths.INPUTS / "large" / "gated_pool_v6_gate_metrics.csv.gz")
    g = gp[gp.design_id.isin(set(lib.design_id))]
    chk = m.merge(g, on="design_id", suffixes=("", "_gp"))
    assert len(chk) == 12000
    for be in ("af3", "esm"):
        assert np.allclose(chk[f"dockq_LRMSD_{be}"], chk[f"rmsd_{be}"], atol=1e-9), be
        assert np.allclose(chk[f"dockq_self_{be}"], chk[f"dockq_{be}"], atol=1e-9), be

    out = m[["design_id", "stack", "length", "sequence"]].copy()
    for stem, (ca, ce) in PAIRS.items():
        out[f"{stem}_af3"], out[f"{stem}_esm"] = m[ca].to_numpy(), m[ce].to_numpy()

    # ligand_RMSD and DockQ, with the BG+BIF true-backbone rescue substituted in.
    bif = pd.read_csv(paths.INPUTS / "large" / "bif_true_per_backend.csv.gz")
    t = m[["design_id", "stack"]].merge(bif, on="design_id", how="left")
    is_bif = (t["stack"] == "BG+BIF").to_numpy()
    for be in ("af3", "esm"):
        for stem, src, true in (("sc_rmsd", f"dockq_LRMSD_{be}", f"rmsd_{be}_true"),
                                ("dockq", f"dockq_self_{be}", f"dockq_{be}_true")):
            out[f"{stem}_{be}"] = np.where(is_bif & t[true].notna().to_numpy(),
                                           t[true].to_numpy(), m[src].to_numpy())
    flip = ((m.dockq_LRMSD_af3 <= m.dockq_LRMSD_esm).to_numpy()
            != (out.sc_rmsd_af3 <= out.sc_rmsd_esm).to_numpy())
    print(f"BG+BIF true-backbone rescue applied to "
          f"{int((is_bif & t.rmsd_af3_true.notna().to_numpy()).sum()):,} of "
          f"{int(is_bif.sum()):,} BG+BIF designs; it changes the lower-ligand_RMSD model for "
          f"{int(flip.sum()):,}/{len(out):,} designs ({int((flip & is_bif).sum()):,} "
          f"inside BG+BIF)")

    # The backend-free gate columns, carried so make_figure.py needs one input fewer.
    out = out.merge(g[["design_id", "seq_run", "seq_max_aa_frac",
                       "footnb_010_af3", "footnb_010_esm"]], on="design_id")
    dest = paths.INPUTS / "metric_significance_backends_12k.csv.gz"
    out.to_csv(dest, index=False)
    print(f"wrote {dest}  {len(out):,} rows x {len(out.columns)} columns")


def main():
    lib = pd.read_csv(paths.INPUTS / "csy4_binder_library_12k_v6_k30_CORRECTED.csv",
                      usecols=["design_id"])
    interface_table(lib)
    backend_table(lib)


if __name__ == "__main__":
    main()
