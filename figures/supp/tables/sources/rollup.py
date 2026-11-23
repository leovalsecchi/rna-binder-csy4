import sys, re
H = {}
for l in sys.stdin:
    if l.startswith("== attribution"): break
    m = re.match(r"  (\S.*?)\s{2,}(\S+)\s+([\d.]+)\s+(\d+)$", l)
    if m: H[(m.group(1).strip(), m.group(2))] = float(m.group(3))
S = {"RFD3+LiMPNN": "rfd3__ligandmpnn", "BD1+LiMPNN": "boltzdesign__ligandmpnn", "HD+LiMPNN": "halludesign__ligandmpnn",
     "BG+BIF": "boltzgen__boltzif", "BG+LiMPNN": "boltzgen__ligandmpnn"}
RF = ("refold_AF3", "refold_B2", "refold_ESM")
g = lambda s, r: H.get((s, r), 0.0)
print("stack gen AF3 B2 ESM refold total")
for k, s in S.items():
    rf = [g(s, r) for r in RF]; gen = g(s, "generation")
    print(k, round(gen, 1), *[round(x, 1) for x in rf], round(sum(rf), 1), round(gen + sum(rf), 1))
bg = g("boltzgen (BIF+LiMPNN shared)", "generation")
print("BG gen shared", bg, " BG combined", round(bg + sum(g(S[k], r) for k in ("BG+BIF", "BG+LiMPNN") for r in RF), 1))
stacks = sum(v for (s, r), v in H.items() if s != "-")
print("all 5 stacks", round(stacks, 1))
for r in RF + ("refold_?", "generation"): print(r, "all", round(sum(v for (s, rr), v in H.items() if rr == r), 1))
for r in ("apo_selfconsistency", "seed_selfconsistency", "other"): print(r, H.get(("-", r)))
print("GRAND", round(sum(H.values()), 1))
