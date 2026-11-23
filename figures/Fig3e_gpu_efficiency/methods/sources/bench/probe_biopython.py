"""Parse a design CIF exactly as refold_af3_scratch.sh does, inside the AF3
container. The worker uses BioPython's MMCIFParser, which hard-requires
_atom_site.B_iso_or_equiv -- biotite omits that column unless the annotation is
registered with set_annotation(). This is the authoritative check that a
rewritten BG+BIF design is readable by the refolder."""
import pathlib
import sys

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import index_to_one, three_to_index


def seq_of(path):
    st = MMCIFParser(QUIET=True).get_structure("s", str(path))
    best = ""
    for model in st:
        for ch in model:
            letters = []
            for res in ch:
                try:
                    letters.append(index_to_one(three_to_index(res.get_resname())))
                except Exception:
                    continue
            s = "".join(letters)
            if len(s) > len(best):
                best = s
    return best


rc = 0
for a in sys.argv[1:]:
    p = pathlib.Path(a)
    try:
        s = seq_of(p)
        ok = len(s) >= 5
        rc |= 0 if ok else 1
        print(f"{'OK  ' if ok else 'SHORT'} {p.name:<50} len={len(s):>3}  {s[:40]}")
    except Exception as e:
        rc = 1
        print(f"FAIL  {p.name:<50} {type(e).__name__}: {e}")
sys.exit(rc)
