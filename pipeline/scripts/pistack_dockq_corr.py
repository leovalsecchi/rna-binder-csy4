"""Correlate pi-stacking backend-AGREEMENT with DockQ.
Reads a CSV that has per-backend pi_stack (af3 = pi_stack_score[_af3], _b2, _esm) and
per-backend DockQ (dockq_af3_msa_off / _b2 / _esm)."""
import csv, sys, math, statistics as st
csv.field_size_limit(10 ** 7)
rows = list(csv.DictReader(open(sys.argv[1])))

def f(x):
    try: return float(x)
    except Exception: return None

DQ = ["dockq_af3_msa_off", "dockq_b2_msa_off", "dockq_esm_msa_off"]
def dockq_cons(r):
    vs = [f(r.get(c)) for c in DQ]; vs = [v for v in vs if v is not None]
    return min(vs) if vs else None

def pearson(xs, ys):
    m = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(m) < 10: return float("nan"), len(m)
    xs = [x for x, _ in m]; ys = [y for _, y in m]
    mx, my = st.mean(xs), st.mean(ys)
    cov = sum((x-mx)*(y-my) for x, y in m)
    sx = math.sqrt(sum((x-mx)**2 for x in xs)); sy = math.sqrt(sum((y-my)**2 for y in ys))
    return (cov/(sx*sy) if sx and sy else float("nan")), len(m)

recs = []
for r in rows:
    a = f(r.get("pi_stack_score_af3")); a = f(r.get("pi_stack_score")) if a is None else a
    b = f(r.get("pi_stack_score_b2")); e = f(r.get("pi_stack_score_esm"))
    dq = dockq_cons(r)
    if None in (a, b, e) or dq is None:
        continue
    recs.append(dict(af3=a, b2=b, esm=e, nb=int(a>0)+int(b>0)+int(e>0),
                     mean=(a+b+e)/3, mn=min(a,b,e), mx=max(a,b,e), dq=dq,
                     dqa=f(r.get("dockq_af3_msa_off")), dqb=f(r.get("dockq_b2_msa_off")),
                     dqe=f(r.get("dockq_esm_msa_off"))))
n = len(recs)
print(f"n designs with per-backend pi_stack + DockQ: {n}\n")

print("DockQ(consensus=min) by # backends showing a stack:")
for lvl in (0, 1, 2, 3):
    g = [x for x in recs if x["nb"] == lvl]
    if g:
        print(f"  {lvl} backend(s): n={len(g):6d}  mean DockQ={st.mean(x['dq'] for x in g):.4f}  "
              f"%DockQ>=0.23={100*sum(1 for x in g if x['dq']>=0.23)/len(g):.1f}%")
    else:
        print(f"  {lvl} backend(s): n=0")

print("\nPearson r vs consensus DockQ:")
for lab, key in [("n_backends(0-3)", "nb"), ("mean pi_stack", "mean"), ("min pi_stack", "mn"),
                 ("max pi_stack", "mx"), ("af3 pi_stack", "af3")]:
    r, nn = pearson([x[key] for x in recs], [x["dq"] for x in recs])
    print(f"  {lab:16s} r={r:+.3f} (n={nn})")

print("\nWithin-backend Pearson r (pi_stack vs SAME-backend DockQ):")
for lab, ps, dq in [("af3", "af3", "dqa"), ("b2", "b2", "dqb"), ("esm", "esm", "dqe")]:
    r, nn = pearson([x[ps] for x in recs], [x[dq] for x in recs])
    print(f"  {lab:4s} r={r:+.3f} (n={nn})")

recs.sort(key=lambda x: x["dq"]); k = n // 10
bot, top = recs[:k], recs[-k:]
def summ(name, g):
    print(f"  {name}: mean n_backends={st.mean(x['nb'] for x in g):.3f}  "
          f"%any-stack={100*sum(1 for x in g if x['nb']>=1)/len(g):.1f}%  "
          f"mean af3_pistack={st.mean(x['af3'] for x in g):.3f}")
print(f"\nDockQ deciles ({k}/decile):")
summ(f"bottom-10% (DockQ<={bot[-1]['dq']:.3f})", bot)
summ(f"top-10%    (DockQ>={top[0]['dq']:.3f})", top)
