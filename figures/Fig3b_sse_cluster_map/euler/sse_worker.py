"""P-SEA secondary structure for every backbone representative in one bbmap shard.

One PDB per coord-hash node (551,516 total across 29 shards), so this is the whole
backbone population, not a sample -- which is what lets per-cluster H/E means be
exact instead of resampled. annotate_sse is the same P-SEA call the campaign's
sse_* columns use; the stored DSSP columns were no use here (refold structures,
and only ~42% campaign coverage).
"""
import tarfile, sys, os, numpy as np
from biotite.structure import annotate_sse, AtomArray

D = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results_bbmap/design"
OUT = "/cluster/work/wicky/lvalsecchi/sse_bb"
sh = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

def ca_array(b):
    xs, resid = [], []
    for line in b.splitlines():
        if line.startswith(b"ATOM") and line[12:16].strip() == b"CA":
            xs.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            resid.append(int(line[22:26]))
    if not xs:
        return None
    a = AtomArray(len(xs))
    a.coord = np.array(xs, dtype=np.float32)
    a.res_id = np.array(resid)
    a.atom_name = np.full(len(xs), "CA")
    a.res_name = np.full(len(xs), "GLY")
    a.chain_id = np.full(len(xs), "A")
    a.element = np.full(len(xs), "C")
    return a

tmp = f"{OUT}/{sh}.tsv.part"
n = 0
with open(tmp, "w") as fh, tarfile.open(f"{D}/shards/{sh}.tar.gz", "r:gz") as tf:
    fh.write("struct_name\tn_res\tn_helix\tn_strand\tn_helix_elem\tn_strand_elem\n")
    for m in tf:
        if not m.name.endswith(".pdb"):
            continue
        a = ca_array(tf.extractfile(m).read())
        if a is None or len(a) < 10:
            continue
        s = "".join(annotate_sse(a))
        runs, i = [], 0
        while i < len(s):
            j = i
            while j < len(s) and s[j] == s[i]:
                j += 1
            runs.append((s[i], j - i)); i = j
        nh = sum(1 for c, L in runs if c == "a" and L >= 4)
        ne = sum(1 for c, L in runs if c == "b" and L >= 3)
        fh.write(f"{m.name[:-4]}\t{len(s)}\t{s.count('a')}\t{s.count('b')}\t{nh}\t{ne}\n")
        n += 1
os.replace(tmp, f"{OUT}/{sh}.tsv")          # atomic: a partial file is never mistaken for done
print(f"{sh}: {n} structures")
