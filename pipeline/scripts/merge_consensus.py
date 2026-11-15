"""Turn per-backend pi_stack columns into min-consensus.
Input CSV has  pi_stack_* (=AF3, unsuffixed), pi_stack_*_b2, pi_stack_*_esm.
Output: adds pi_stack_*_af3 (copy of current), and REWRITES pi_stack_* = min(af3,b2,esm)
so the explorer's existing pi_stack_score/count/best become the consensus."""
import csv, sys
csv.field_size_limit(10 ** 7)
inp, out = sys.argv[1], sys.argv[2]
BASE = ["pi_stack_score", "pi_stack_count", "pi_stack_best", "pi_stack_res_aa",
        "pi_stack_res_nt", "pi_stack_count_tshape",
        "pi_stack_TRP", "pi_stack_TYR", "pi_stack_PHE", "pi_stack_HIS", "pi_stack_ARG"]
FLOAT = {"pi_stack_score", "pi_stack_best"}          # rest are integer counts

rows = list(csv.DictReader(open(inp)))
fn = list(rows[0].keys())
for b in BASE:                                        # add X_af3 right after consensus X
    if b + "_af3" not in fn:
        fn.insert(fn.index(b) + 1, b + "_af3")

def f(x):
    try: return float(x)
    except Exception: return None

n_cons = 0
for r in rows:
    row_has = False
    for b in BASE:
        a_str = r.get(b, "")                          # current value IS the AF3 value
        r[b + "_af3"] = a_str
        vs = [f(a_str), f(r.get(b + "_b2")), f(r.get(b + "_esm"))]
        vs = [v for v in vs if v is not None]
        if not vs:
            r[b] = ""
            continue
        mn = min(vs)
        r[b] = round(mn, 4) if b in FLOAT else int(round(mn))
        if b == "pi_stack_score" and mn > 0:
            row_has = True
    n_cons += int(row_has)

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fn); w.writeheader(); w.writerows(rows)
print(f"wrote {out}  rows={len(rows)}  consensus pi_stack_score>0: {n_cons} ({100*n_cons/len(rows):.1f}%)")
