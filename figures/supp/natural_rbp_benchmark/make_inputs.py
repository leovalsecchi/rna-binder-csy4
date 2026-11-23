"""Emit the two tables behind the natural-RBP benchmark: panel (a)'s interface sizes
WITH per-metric n, and panel (d)'s selection funnel.

Why this exists: neither had a machine-readable source.

The size ratios lived as prose in `nat_vs_design_findings/NORMALIZED_FINDINGS.md` and as
literals in the old `fig_F8.py`. `normalize_metrics.py` prints the block but does not save
it, and the printed block carries no n -- which matters, because `dsasa_total` is n = 53,
not 58, so a blanket "n = 58" on panel (a) would be wrong.

The SELECTION of the set was worse off: it existed only as prose in `~/rbp_rna_db/README.md`
and as the SIZE OF TWO DIRECTORIES. Nothing recorded how 58 complexes were arrived at, and
the figure used the set without ever characterising it (Leo 2026-11-20: "Can you add a
general overview of the dataset we used?"). Every count below is RE-DERIVED from the set's
own build artefacts, with the same rule each build script used, and asserted against the
number that script printed -- no count is copied out of the README.

Same two primary files, same filter and same estimator as `normalize_metrics.py`
(median per side, per-metric n). Run once; the figure reads the CSVs.

    ~/.venvs/binder_map/bin/python3 code/supp/natural_rbp_benchmark/make_inputs.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd

from common import paths

DB = Path.home() / "rbp_rna_db"
NAT = DB / "subset59" / "metrics_natural_59.csv"
DES = Path.home() / "Documents" / "metric_explorer" / "data" / "metrics_100k.csv"
OUT = paths.INPUTS / "natural_rbp_benchmark_size_confound.csv"
OUT_SEL = paths.INPUTS / "natural_rbp_benchmark_selection.csv"

SIZE = {
    "n_contact_res_prot": "Protein interface residues",
    "n_contact_res_rna": "RNA interface residues",
    "dsasa_total": "Interface ΔSASA (Å²)",
}


def chains(s):
    """Comma list of chain ids -> list, the way build_subset59.py's n() reads it."""
    return [] if pd.isna(s) else [c for c in str(s).split(",") if c]


def selection():
    """The funnel from the RCSB search to the 58 complexes this figure uses.

    Each step is recomputed from the artefact the step that produced it left behind,
    not read off `~/rbp_rna_db/README.md`:

      search   `cif/` + `cif_dropped/`  -- refine_rna_len.py MOVED the rejects aside
                                           rather than deleting them, so the two
                                           directories together still are its input.
      RNA      `manifest.csv`           -- rewritten by refine_rna_len.py to the kept set.
      chains   `chain_map.tsv`          -- re-filtered with build_subset59.py's own rule.
      metrics  `metrics_natural_59.csv` -- minus the rows whose biotite SASA failed.
    """
    kept = sorted(p.stem for p in (DB / "cif").glob("*.cif"))
    dropped = sorted(p.stem for p in (DB / "cif_dropped").glob("*.cif"))
    man = pd.read_csv(DB / "manifest.csv")
    cm = pd.read_csv(DB / "chain_map.tsv", sep="\t")
    nat = pd.read_csv(NAT)

    assert sorted(man.pdb_id) == kept, "manifest.csv and cif/ disagree on the kept set"
    assert not set(kept) & set(dropped), "a structure is in both cif/ and cif_dropped/"

    # build_subset59.py: exactly one protein chain, exactly one RNA chain, no DNA and
    # no other polymer. Recomputed rather than trusted, then checked against the file
    # that script wrote.
    one = cm[cm.protein_chains.map(lambda s: len(chains(s)) == 1)
             & cm.rna_chains.map(lambda s: len(chains(s)) == 1)
             & cm.dna_chains.map(lambda s: not chains(s))
             & cm.other_chains.map(lambda s: not chains(s))]
    assert sorted(one.pdb_id) == sorted(pd.read_csv(DB / "subset59" / "manifest_59.csv").pdb_id), \
        "the one-protein-one-RNA rule no longer reproduces manifest_59.csv"

    err = nat["dsasa_biotite_error"]
    ok = nat[err.isna() | (err.astype(str).str.strip() == "")]

    rows = [
        ("search", "RCSB, ≤ 2.5 Å, 30% identity", len(kept) + len(dropped),
         "≥ 1 protein and ≥ 1 RNA entity, ≤ 1,500 deposited polymer residues "
         "(excludes ribosomes and spliceosomes), one representative per RCSB 30% "
         "protein-identity cluster"),
        ("rna_length", "Longest RNA ≥ 10 nt", len(kept),
         "refine_rna_len.py MIN_NT = 10; the rejects are in cif_dropped/"),
        ("two_chains", "One protein + one RNA chain", len(one),
         "build_subset59.py: no DNA, no second copy of either polymer, so every "
         "complex has one interface and a uniform A/B chain layout"),
        ("metrics", "Interface metrics computed", len(ok),
         "8JOZ drops: its biotite SASA failed on a modified residue"),
    ]
    out = pd.DataFrame(rows, columns=["step", "label", "n", "note"])
    out["of_search"] = (out.n / out.n.iloc[0] * 100).round(1)
    # The second n this figure carries: ΔSASA is Rosetta InterfaceAnalyzer and is
    # empty for five of the 58, so panel (a)'s ΔSASA bar rests on 53. Not a further
    # funnel step - a per-metric subset of the final one.
    out["n_dsasa"] = ""
    out.loc[out.index[-1], "n_dsasa"] = int(
        pd.to_numeric(ok["dsasa_total"], errors="coerce").notna().sum())
    assert list(out.n) == [329, 180, 59, 58], f"funnel changed: {list(out.n)}"
    return out


def main():
    nat = pd.read_csv(NAT)
    # 8JOZ drops here: its biotite SASA failed, so every SASA-derived value is void.
    err = nat["dsasa_biotite_error"]
    nat = nat[err.isna() | (err.astype(str).str.strip() == "")]
    assert len(nat) == 58, f"expected 58 natural complexes after the SASA filter, got {len(nat)}"

    des = pd.read_csv(DES, low_memory=False,
                      usecols=lambda c: c in SIZE)
    assert len(des) == 100_000, f"expected 100,000 designs, got {len(des)}"

    rows = []
    for col, label in SIZE.items():
        n = pd.to_numeric(nat[col], errors="coerce").dropna()
        d = pd.to_numeric(des[col], errors="coerce").dropna()
        rows.append({
            "metric": col,
            "label": label,
            "nat_median": round(n.median(), 1),
            "nat_n": len(n),
            "des_median": round(d.median(), 1),
            "des_n": len(d),
            "ratio_nat_over_des": round(n.median() / d.median(), 2),
        })

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print("\nwrote", OUT)

    sel = selection()
    sel.to_csv(OUT_SEL, index=False)
    print()
    print(sel.drop(columns="note").to_string(index=False))
    print("\nwrote", OUT_SEL)


if __name__ == "__main__":
    main()
