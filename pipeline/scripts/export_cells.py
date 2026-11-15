"""Export designs sharing a final (fold x epitope) cell, RNA-superposed on the first."""
import os, sys, shutil
import numpy as np, pandas as pd
HERE = "/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts"
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "src"))
import biotite.structure as struc
from biotite.structure.io import save_structure
from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains
from lrmsd_rna_aligned import RNA_BB, _bb, _match, _apply, MIN_REGISTER

FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
RES = FM + "/results_v6/design/"
OUT = FM + "/same_cell_designs"
NCELL, PER = 8, 4

raw = pd.read_csv(FM + "/passing_v6_final.csv", skiprows=2, low_memory=False,
                  usecols=["design_id", "stack", "iptm_af3", "iptm_esm", "design_cif"])
raw["iptm"] = np.maximum(raw.iptm_af3, raw.iptm_esm)
nodes = pd.read_csv(RES + "dedup/nodes.tsv", sep="\t")
mem = pd.read_csv(RES + "dedup/members.tsv", sep="\t")
rep, m = [], []
for sn, ids in zip(mem.struct_name, mem.design_ids):
    for i in str(ids).split(","): rep.append(sn); m.append(i)
nid = pd.Series(nodes.node_id.to_numpy(), index=nodes.struct_name)
raw["node"] = pd.Series(nid.reindex(pd.Index(rep)).to_numpy(), index=m).reindex(raw.design_id).to_numpy()
for ax, col in (("fold", "fold"), ("epitope", "epi")):
    c = pd.read_parquet(RES + f"clusters/{ax}/clusters.parquet").set_index("node_id")
    raw[col] = c["cluster_r80"].reindex(raw["node"]).to_numpy()
raw = raw.dropna(subset=["fold", "epi", "node"])

lib = set(pd.read_parquet(LS + "/hier_k15_max_selections.parquet")
          .query('method == "HIER fold->epitope"').design_id)
L = raw[raw.design_id.isin(lib)]
# cells the library actually populates with several DISTINCT backbones
g = L.groupby(["fold", "epi"]).agg(n=("design_id", "size"), nbb=("node", "nunique"))
cells = g[g.nbb >= 3].sort_values("nbb", ascending=False)
print(f"library cells with >=3 distinct backbones: {len(cells)}")
rng = np.random.default_rng(0)
pick = cells.iloc[rng.choice(len(cells), min(NCELL, len(cells)), replace=False)]

def full(p):
    a = load_structure(resolve(p)); return a[a.element != "H"]

if os.path.isdir(OUT): shutil.rmtree(OUT)
os.makedirs(OUT)
man = []
for (f, e), row in pick.iterrows():
    sub = L[(L.fold == f) & (L.epi == e)].drop_duplicates("node").nlargest(PER, "iptm")
    tag = f"cell_f{int(f):04d}_e{int(e):04d}_n{int(row.nbb)}"
    d = os.path.join(OUT, tag); os.makedirs(d)
    ref = None
    for i, r in enumerate(sub.itertuples()):
        try:
            A = full(r.design_cif)
            rna = _bb(split_chains(A)[1], RNA_BB)
            if ref is None:
                ref = rna
            else:
                ia, ib, reg = _match(ref, rna)
                if len(ia) < 3 or reg < MIN_REGISTER:
                    print(f"  {tag}: skip {r.design_id[-28:]} register={reg:.2f}"); continue
                _, tr = struc.superimpose(ref[ia], rna[ib]); A = _apply(tr, A)
            save_structure(os.path.join(d, f"{i}_{r.stack.replace('+','_')}_{r.design_id}.cif"), A)
            man.append(dict(cell=tag, fold=int(f), epi=int(e), n_backbones=int(row.nbb),
                            design_id=r.design_id, stack=r.stack, iptm=round(r.iptm, 3)))
        except Exception as exc:
            print(f"  {tag}: ERROR {exc}")
    print(f"  {tag}: {len(sub)} designs, stacks={sorted(set(sub["stack"]))}")
pd.DataFrame(man).to_csv(os.path.join(OUT, "manifest.csv"), index=False)
print(f"\nwrote {len(man)} cifs across {len(pick)} cells -> {OUT}")
