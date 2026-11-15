#!/usr/bin/env python
"""DockQ over the top-B2 designs (refold vs designed reference = self-consistency).
Reads top150_b2_for_dockq.csv (cols incl abs_path=B2 refold, designed_source=designed ref,
both as '<archive.tar.zst>!<member>'), extracts members over the eulermount, runs DockQ."""
import sys, subprocess, tempfile, shlex
from pathlib import Path
from collections import defaultdict
import pandas as pd, gemmi
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
sys.path.insert(0, "/Users/leo/Documents/rna_binder_gui")
from ligand_rmsd import _kind

EU, MB = "/cluster/work/wicky", "/Volumes/eulermount"
BASE = "/Users/leo/Documents/rna_binder_gui/analyses/midrun_v2_all"
df = pd.read_csv(f"{BASE}/top150_b2_for_dockq.csv")

# 1) collect needed members per archive, extract each archive once
need = defaultdict(set)
for r in df.itertuples():
    for ref in (r.abs_path, r.designed_source):
        a, m = ref.split("!", 1); need[a.replace(EU, MB, 1)].add(m)
cache = {}
print(f"extracting {len(need)} archives over mount ...", flush=True)
for j, (a, members) in enumerate(need.items(), 1):
    scr = Path(tempfile.mkdtemp(prefix="dq_"))
    cmd = "zstd -dc %s | tar -x -C %s %s" % (shlex.quote(a), shlex.quote(str(scr)),
                                             " ".join(shlex.quote(m) for m in members))
    subprocess.run(cmd, shell=True, stderr=subprocess.DEVNULL)
    cache[a] = scr
    if j % 20 == 0: print(f"  {j}/{len(need)} archives", flush=True)

def member_path(ref):
    a, m = ref.split("!", 1); p = cache[a.replace(EU, MB, 1)] / m
    return p if p.exists() and p.stat().st_size > 0 else None

def ctypes(p):
    st = gemmi.read_structure(str(p)); o = {}
    for c in st[0]:
        k = {_kind(x.name) for x in c} - {"other"}
        if "protein" in k: o[c.name] = "protein"
        elif "na" in k: o[c.name] = "na"
    return o

def cmap(m, n):
    mt = ctypes(m); nbt = {t: ch for ch, t in ctypes(n).items()}
    return {mc: nbt[t] for mc, t in mt.items() if t in nbt}

# 2) DockQ each design (refold vs designed reference)
rows = []
for i, r in enumerate(df.itertuples(), 1):
    rec = dict(stack=r.stack, uuid8=r.uuid8, design_id=r.design_id, seq_id=r.seq_id,
               iptm=round(r.b2_protein_rna_iptm, 4), caRMSD=round(r.ca_rmsd_vs_designed, 2),
               DockQ=None, LRMSD=None, iRMSD=None, fnat=None, err="")
    try:
        m, n = member_path(r.abs_path), member_path(r.designed_source)
        if not m or not n: raise FileNotFoundError("member missing")
        res, _ = run_on_all_native_interfaces(load_PDB(str(m)), load_PDB(str(n)), chain_map=cmap(m, n))
        if res:
            v = next(iter(res.values()))
            rec.update(DockQ=round(v["DockQ"], 4), LRMSD=round(v["LRMSD"], 3),
                       iRMSD=round(v["iRMSD"], 3), fnat=round(v["fnat"], 4))
        else:
            rec["err"] = "no interface"
    except Exception as e:
        rec["err"] = repr(e)[:100]
    rows.append(rec)
    if i % 15 == 0 or rec["err"]:
        print(f"[{i}/{len(df)}] {r.uuid8} d{r.design_id}s{r.seq_id} DockQ={rec['DockQ']} {rec['err']}", flush=True)

out = pd.DataFrame(rows).sort_values("DockQ", ascending=False, na_position="last")
out.to_csv(f"{BASE}/dockq_top150_b2.csv", index=False)
ok = out["DockQ"].notna().sum()
print(f"\n=== {ok}/{len(out)} scored | DockQ>=0.8 high:{int((out.DockQ>=0.8).sum())} "
      f">=0.49 medium+:{int((out.DockQ>=0.49).sum())} >=0.23 acceptable+:{int((out.DockQ>=0.23).sum())} ===")
print("\n=== TOP 25 B2 designs by DockQ ===")
print(out.head(25)[["stack","uuid8","design_id","seq_id","DockQ","LRMSD","fnat","iptm","caRMSD"]].to_string(index=False))
