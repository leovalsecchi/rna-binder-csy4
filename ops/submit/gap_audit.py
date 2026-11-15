#!/usr/bin/env python3
"""midrun_v4 gap audit: for GEN and REFOLD, show target/done/in-flight(queued)/genuine-gap.
In-flight = covered by a live squeue job (array SCRIPTS or direct run.sbatch)."""
import re, os, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
C = RB + "/campaigns/csy4_crrna_free"
LABEL = "midrun_v4"


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True).stdout


def live_covered():
    """Every run.sbatch path covered by a live job (array SCRIPTS members + direct run.sbatch)."""
    live = set()
    for o in sh("squeue", "-u", "lvalsecchi", "-h", "-o", "%o").splitlines():
        o = o.strip()
        if not o.endswith(".sbatch"):
            continue
        if o.endswith("/run.sbatch"):
            live.add(o); continue
        try:
            txt = open(o).read()
        except OSError:
            continue
        for m in re.finditer(r"^\s*(\S+/run\.sbatch)\s*$", txt, re.M):
            live.add(m.group(1).strip())
    return live


def u8(p):
    m = re.search(r"([0-9a-f]{8})", os.path.basename(os.path.dirname(p)))
    return m.group(1) if m else None


LIVE = live_covered()
# split live coverage into gen vs refold(stack,backend)->uuids
rf_live = collections.defaultdict(set)
for p in LIVE:
    m = re.search(r"/([^/]+)/%s/refold_(af3|b2|esm)/\w*?([0-9a-f]{8})" % LABEL, p)
    if m:
        rf_live[(m.group(1), m.group(2))].add(m.group(3))

print("live-covered run.sbatch total: %d\n" % len(LIVE))

# ---------- GEN ----------
GEN = [("rfd3__ligandmpnn", 141, "multif"),
       ("boltzgen__boltzif", 150, "multif"),
       ("boltzdesign__ligandmpnn", 150, "multif"),
       ("halludesign__ligandmpnn", 600, "hd")]
print("=== GEN gaps (uuids) ===")
print("%-26s %6s %6s %9s %5s %8s" % ("generator", "target", "done", "inflight", "gap", "notmade"))
for sk, tgt, layout in GEN:
    if layout == "multif":
        runs = glob.glob(f"{C}/{sk}/{LABEL}/run_*/run.sbatch")
        packed = {os.path.basename(x)[:8] for x in glob.glob(f"{C}/{sk}/{LABEL}/designs/*.tar.zst")}
    else:  # hd per-uuid, combined into shared
        runs = glob.glob(f"{C}/{sk}/{LABEL}_*/run.sbatch")
        packed = {os.path.basename(x)[:8] for x in glob.glob(f"{C}/{sk}/{LABEL}/designs/*.tar.zst")}
        for r in runs:
            if os.path.exists(os.path.join(os.path.dirname(r), "designs.tar.zst")):
                m = re.search(r"%s_([0-9a-f]{8})" % LABEL, r)
                if m:
                    packed.add(m.group(1))
    ru = {re.search(r"([0-9a-f]{8})", r).group(1): r for r in runs if re.search(r"([0-9a-f]{8})", r)}
    done = sum(1 for u in ru if u in packed)
    infl = sum(1 for u, r in ru.items() if u not in packed and r in LIVE)
    gap = sum(1 for u, r in ru.items() if u not in packed and r not in LIVE)
    notmade = max(0, tgt - len(ru))
    print("%-26s %6d %6d %9d %5d %8d" % (sk, tgt, done, infl, gap, notmade))

# ---------- REFOLD ----------
RF = ["rfd3__ligandmpnn", "boltzdesign__ligandmpnn", "boltzgen__boltzif",
      "boltzgen__ligandmpnn", "halludesign__ligandmpnn"]
print("\n=== REFOLD gaps (uuids of packed designs missing a backend) ===")
print("%-26s %6s | %-18s | %-18s | %-18s" % ("stack", "packed",
      "AF3 done/queue/gap", "ESM done/queue/gap", "B2  done/queue/gap"))
for sk in RF:
    packed = {os.path.basename(x)[:8] for x in glob.glob(f"{C}/{sk}/{LABEL}/designs/*.tar.zst")}
    cells = []
    for be in ("af3", "esm", "b2"):
        done = set()
        for leaf in glob.glob(f"{C}/{sk}/{LABEL}/refold_{be}/*/designs.tar.zst"):
            u = u8(leaf)
            if u:
                done.add(u)
        q = rf_live.get((sk, be), set()) & packed
        gap = packed - done - q
        cells.append("%3d /%4d /%4d" % (len(done & packed), len(q), len(gap)))
    print("%-26s %6d | %-18s | %-18s | %-18s" % (sk, len(packed), cells[0], cells[1], cells[2]))
