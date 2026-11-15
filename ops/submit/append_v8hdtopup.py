#!/usr/bin/env python3
"""append_v8hdtopup.py — append the v8 HD topup rows to csy4_master_foot.csv.

A pure ROW APPEND: every topup design_id was verified disjoint from what the master already
holds, so nothing is replaced. The master half is copied BYTE-EXACT (binary copyfileobj) —
never re-parsed — so no float in 6.8M existing rows is reformatted, which is what a pandas
round-trip would do.

Not using merge_masters.py for two reasons:
  * the master ALREADY carries a duplicate `source_table` column (pre-existing, from the two
    chained merges); merge_masters appends another unconditionally -> a third copy.
  * merge_masters ABORTS on b-only columns, and the topup carries 4 scratch intermediates
    (lrmsd_rna_*/lrmsd_prot_*) that the master deliberately does not keep.
Duplicate header names are handled POSITIONALLY: both `source_table` slots get the same value.
"""
import csv, os, shutil, sys

csv.field_size_limit(10 ** 7)
CC = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
MASTER, TOPUP = CC + "csy4_master_foot.csv", CC + "v8hdtopup_master.csv"
OUT = CC + "csy4_master_foot_v8top.csv"
DROP = {"lrmsd_rna_af3", "lrmsd_rna_esm", "lrmsd_prot_af3", "lrmsd_prot_esm"}
EXPECT_MASTER = 6_821_957
EXPECT_TOPUP = 170_558

def header_of(p):
    with open(p, newline="") as f:
        return [h.strip() for h in next(csv.reader(f))]

ha, hb = header_of(MASTER), header_of(TOPUP)
sa = set(ha)
extra = [c for c in hb if c not in sa and c not in DROP]
if extra:
    sys.exit(f"ABORT: topup has {len(extra)} cols not in master: {extra}")
print(f"master {len(ha)} cols, topup {len(hb)} cols, dropping {sorted(DROP)}", flush=True)
print(f"master cols topup lacks: {len([c for c in ha if c not in set(hb)])} (blank-filled)", flush=True)

# ---- 1. byte-exact copy of the master, counting newlines as we go -------------------
nl = 0
with open(MASTER, "rb") as f, open(OUT + ".tmp", "wb") as g:
    while True:
        buf = f.read(64 << 20)
        if not buf:
            break
        nl += buf.count(b"\n")
        g.write(buf)
        last = buf[-1:]
print(f"copied master: {nl:,} newlines", flush=True)
if last != b"\n":
    with open(OUT + ".tmp", "ab") as g:
        g.write(b"\n")
    nl += 1
    print("master lacked a trailing newline — added one", flush=True)
master_rows = nl - 1
assert master_rows == EXPECT_MASTER, f"master rows {master_rows:,} != {EXPECT_MASTER:,} (embedded newlines?)"

# ---- 2. append topup rows aligned to the master header ------------------------------
idx = {c: i for i, c in enumerate(hb)}          # first occurrence wins
n = 0
with open(TOPUP, newline="") as f, open(OUT + ".tmp", "a", newline="") as g:
    r = csv.reader(f)
    next(r)
    w = csv.writer(g, lineterminator="\n")
    for row in r:
        w.writerow(["v8hdtopup" if c == "source_table"
                    else (row[idx[c]] if c in idx else "") for c in ha])
        n += 1
assert n == EXPECT_TOPUP, f"appended {n:,} != {EXPECT_TOPUP:,}"
os.replace(OUT + ".tmp", OUT)
print(f"appended {n:,} rows -> {OUT}", flush=True)
print(f"TOTAL {master_rows + n:,} rows x {len(ha)} cols", flush=True)
