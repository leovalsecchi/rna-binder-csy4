import os, sys
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0,os.path.join(HERE,"..","src"))
import biotite.structure as struc
from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains
from lrmsd_rna_aligned import RNA_BB, _bb, _match, _apply
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
man=pd.read_csv(os.path.join(FM,"deisolated_pairs","manifest.csv")); man=man[man.status=="ok"]
cols=["design_id","iptm_mean","dockq_af3","dockq_esm","rmsd_af3","rmsd_esm","design_cif"]
p=pd.read_csv(os.path.join(FM,"passing_v6_final.csv"),skiprows=2,usecols=cols,low_memory=False).set_index("design_id")
def full(path):
    a=load_structure(path); return a[a.element!="H"]
rows=[]
for _,r in man.iterrows():
    A=full(resolve(p.loc[r.ex_singleton,"design_cif"])); B=full(resolve(p.loc[r.clustermate,"design_cif"]))
    pa,ra=split_chains(A); pb,rb=split_chains(B)
    ia,ib,reg=_match(_bb(ra,RNA_BB),_bb(rb,RNA_BB))
    _,tr=struc.superimpose(_bb(ra,RNA_BB)[ia],_bb(rb,RNA_BB)[ib])
    pb2=_apply(tr,pb)
    ca_a=pa[pa.atom_name=="CA"]; ca_b=pb2[pb2.atom_name=="CA"]
    cen=float(np.linalg.norm(ca_a.coord.mean(0)-ca_b.coord.mean(0)))
    n=min(len(ca_a),len(ca_b))
    pose=float(np.sqrt(((ca_a.coord[:n]-ca_b.coord[:n])**2).sum(1).mean()))
    a_=p.loc[r.ex_singleton]; b_=p.loc[r.clustermate]
    da=max(a_.dockq_af3,a_.dockq_esm); db=max(b_.dockq_af3,b_.dockq_esm)
    rows.append(dict(case=r.case[:14], tm=round(r.best_tm,3), cen_A=round(cen,1), poseRMSD=round(pose,1),
        iptm_a=round(a_.iptm_mean,3), iptm_b=round(b_.iptm_mean,3), d_iptm=round(abs(a_.iptm_mean-b_.iptm_mean),3),
        dq_a=round(da,3), dq_b=round(db,3), d_dq=round(abs(da-db),3),
        rmsd_a=round(min(a_.rmsd_af3,a_.rmsd_esm),1), rmsd_b=round(min(b_.rmsd_af3,b_.rmsd_esm),1)))
d=pd.DataFrame(rows).sort_values("tm")
pd.set_option("display.width",250)
print(d.to_string(index=False)); print()
for lab,g in (("TM>=0.80",d[d.tm>=0.80]),("TM 0.60-0.80",d[(d.tm>=0.60)&(d.tm<0.80)]),("TM<0.60",d[d.tm<0.60])):
    print("%-13s n=%2d  centroid %5.1f A  poseRMSD %5.1f A  |d_iptm| %.3f  |d_dockq| %.3f"
          %(lab,len(g),g.cen_A.median(),g.poseRMSD.median(),g.d_iptm.median(),g.d_dq.median()))
print()
print("corr(TM, centroid dist) = %+.3f   corr(TM, |d_iptm|) = %+.3f   corr(TM, |d_dockq|) = %+.3f"
      %(d.tm.corr(d.cen_A), d.tm.corr(d.d_iptm), d.tm.corr(d.d_dq)))
