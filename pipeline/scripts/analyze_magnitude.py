"""Does the MAGNITUDE of pi-stacking (how many / how strong) track DockQ?
Separates 'many pi-pi' from mere backend-agreement."""
import csv, sys, math, statistics as st
from collections import defaultdict, Counter
csv.field_size_limit(10 ** 7)
rows = list(csv.DictReader(open(sys.argv[1])))

def f(x):
    try: return float(x)
    except Exception: return None
DQ = ["dockq_af3_msa_off", "dockq_b2_msa_off", "dockq_esm_msa_off"]
def dq(r):
    vs = [f(r.get(c)) for c in DQ]; vs = [v for v in vs if v is not None]
    return min(vs) if vs else None
def pearson(pairs):
    m = [(x, y) for x, y in pairs if x is not None and y is not None]
    if len(m) < 10: return float("nan"), len(m)
    xs = [a for a, _ in m]; ys = [b for _, b in m]; mx, my = st.mean(xs), st.mean(ys)
    cov = sum((a-mx)*(b-my) for a, b in m)
    sx = math.sqrt(sum((a-mx)**2 for a in xs)); sy = math.sqrt(sum((b-my)**2 for b in ys))
    return (cov/(sx*sy) if sx and sy else float("nan")), len(m)

recs = []
for r in rows:
    ac, bc, ec = f(r.get("pi_stack_count")), f(r.get("pi_stack_count_b2")), f(r.get("pi_stack_count_esm"))
    a, b, e = f(r.get("pi_stack_score")), f(r.get("pi_stack_score_b2")), f(r.get("pi_stack_score_esm"))
    d = dq(r)
    if None in (ac, bc, ec, a, b, e) or d is None: continue
    recs.append(dict(ac=ac, as_=a, min_count=min(ac, bc, ec), max_count=max(ac, bc, ec),
                     min_score=min(a, b, e), max_score=max(a, b, e), dq=d))
n = len(recs); print(f"n={n}")

def bins_by(key, bucket, label):
    g = defaultdict(list)
    for x in recs: g[bucket(x[key])].append(x)
    print(f"\n{label}")
    for bk in sorted(g):
        gg = g[bk]
        print(f"  {bk}: n={len(gg):6d}  mean DockQ={st.mean(x['dq'] for x in gg):.4f}  "
              f"%>=0.23={100*sum(1 for x in gg if x['dq']>=0.23)/len(gg):.1f}%")
    r, _ = pearson([(x[key], x["dq"]) for x in recs]); print(f"  Pearson r({key},DockQ)={r:+.3f}")

bins_by("ac", lambda v: min(int(v), 3), "A. DockQ by AF3 pi_stack_count (many pi-pi in the AF3 pose; 3 = 3+):")
bins_by("max_count", lambda v: min(int(v), 3), "B. DockQ by MAX stack count over backends (best-case magnitude; 3 = 3+):")
bins_by("min_count", lambda v: min(int(v), 2), "C. DockQ by MIN stack count over backends (robust magnitude; 2 = 2+):")

allg = [x for x in recs if x["min_score"] > 0]
print(f"\nD. the all-3-agree group (min_score>0): n={len(allg)}")
if allg:
    print("   how many pi-pi they actually have (robust count = min over backends):",
          dict(sorted(Counter(int(x["min_count"]) for x in allg).items())))
    print(f"   mean DockQ (all-agree) = {st.mean(x['dq'] for x in allg):.4f}")
    r, _ = pearson([(x["min_count"], x["dq"]) for x in allg])
    print(f"   within all-agree: r(robust_count, DockQ)={r:+.3f}")
    one = [x for x in allg if x["min_count"] <= 1]; many = [x for x in allg if x["min_count"] >= 2]
    if one: print(f"   1 robust stack : n={len(one):3d}  mean DockQ={st.mean(x['dq'] for x in one):.4f}")
    if many: print(f"   >=2 robust stacks: n={len(many):3d}  mean DockQ={st.mean(x['dq'] for x in many):.4f}")
