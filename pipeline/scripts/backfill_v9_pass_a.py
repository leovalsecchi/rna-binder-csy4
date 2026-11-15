#!/usr/bin/env python3
"""backfill_v9_pass_a.py — fill the v9 half of a merged master IN PLACE.

Unlike fix_iptm_master.py / apply_dockq_rna.py (which APPEND their columns), this fills
columns that already exist in the header, because the v2-v9 master carries them for the
v2-v8 rows and leaves them blank on v9.

(1) DockQ convention. v9 rows currently hold STOCK-convention dockq_self/dockq_LRMSD while
    v2-v8 hold the RNA-aligned (target-fixed) values -> the column is mixed. Move the stock
    numbers into the existing dockq_std_* columns and install the RNA-aligned ones from the
    sidecar, rebuilding the bare consensus with the project convention (min for dockq_self,
    max for dockq_LRMSD, blank unless all three backends scored).
(2) ipTM. Fill iptm_min3 / iptm_best2 / mean_iptm3 / consensus_iptm where blank, gated on
    all three per-backend ipTMs.

Rows already carrying a value are never overwritten -- they are audited and counted, so a
re-run is a no-op. v2-v8 rows (absent from the v9 sidecar) are left untouched.
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(10 ** 7)
BE = ["af3", "b2", "esm"]
IPTM = ["af3_msaoff_iptm", "b2_iptm", "esm_iptm"]


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--sidecar", required=True, help="dockq_rna_v9_new.csv (RNA-aligned DockQ)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

# ---- sidecar: design_id -> (dockq_rna_af3,b2,esm, lrmsd_rna_corr_af3,b2,esm) -------------
print("loading sidecar...", flush=True)
side = {}
with open(a.sidecar, newline="") as f:
    for r in csv.DictReader(f):
        side[r["design_id"]] = (r.get("dockq_rna_af3", ""), r.get("dockq_rna_b2", ""),
                                r.get("dockq_rna_esm", ""), r.get("lrmsd_rna_corr_af3", ""),
                                r.get("lrmsd_rna_corr_b2", ""), r.get("lrmsd_rna_corr_esm", ""))
print(f"sidecar: {len(side):,} designs", flush=True)

with open(a.master, newline="") as f:
    rd = csv.reader(f)
    header = [h.strip() for h in next(rd)]
    ix = {c: i for i, c in enumerate(header)}

    need = ([f"dockq_self_{b}" for b in BE] + [f"dockq_LRMSD_{b}" for b in BE]
            + [f"dockq_std_self_{b}" for b in BE] + [f"dockq_std_LRMSD_{b}" for b in BE]
            + ["dockq_self", "dockq_LRMSD", "dockq_std_self", "dockq_std_LRMSD"]
            + IPTM + ["iptm_min3", "iptm_best2", "mean_iptm3", "consensus_iptm"])
    missing = [c for c in need if c not in ix]
    if missing:
        sys.exit(f"ABORT: master is missing expected columns {missing}")
    i_did = ix["design_id"]
    i_camp = ix.get("campaign")

    n = hit = dq_applied = dq_already = iptm_filled = iptm_audit = iptm_mismatch = 0
    v9_no_sidecar = 0
    tmp = a.out + ".tmp"
    with open(tmp, "w", newline="") as g:
        w = csv.writer(g, lineterminator="\n")
        w.writerow(header)
        for row in rd:
            n += 1
            src = side.get(row[i_did])

            # ---------- (1) DockQ ----------
            if src is not None:
                hit += 1
                if row[ix["dockq_std_self_af3"]] != "":
                    dq_already += 1                      # already converted; leave alone
                else:
                    dq_applied += 1
                    per_self, per_lr = [], []
                    for j, b in enumerate(BE):
                        cs, cl = f"dockq_self_{b}", f"dockq_LRMSD_{b}"
                        row[ix[f"dockq_std_self_{b}"]] = row[ix[cs]]
                        row[ix[cs]] = src[j]
                        per_self.append(fnum(src[j]))
                        row[ix[f"dockq_std_LRMSD_{b}"]] = row[ix[cl]]
                        row[ix[cl]] = src[3 + j]
                        per_lr.append(fnum(src[3 + j]))
                    row[ix["dockq_std_self"]] = row[ix["dockq_self"]]
                    ok = [v for v in per_self if v is not None]
                    row[ix["dockq_self"]] = f"{min(ok):.6f}" if len(ok) == 3 else ""
                    row[ix["dockq_std_LRMSD"]] = row[ix["dockq_LRMSD"]]
                    ok = [v for v in per_lr if v is not None]
                    row[ix["dockq_LRMSD"]] = f"{max(ok):.4f}" if len(ok) == 3 else ""
            elif i_camp is not None and row[i_camp] == "v9":
                v9_no_sidecar += 1

            # ---------- (2) ipTM ----------
            if row[ix["iptm_min3"]] == "":
                v = [fnum(row[ix[c]]) for c in IPTM]
                ok = [x for x in v if x is not None]
                if len(ok) == 3:
                    iptm_filled += 1
                    ok.sort()
                    row[ix["iptm_min3"]] = f"{ok[0]:.6g}"
                    row[ix["iptm_best2"]] = f"{ok[1]:.6g}"
                    row[ix["mean_iptm3"]] = f"{sum(ok) / 3.0:.6g}"
                    row[ix["consensus_iptm"]] = f"{ok[0]:.6g}"
            else:
                iptm_audit += 1
                m, c = fnum(row[ix["iptm_min3"]]), fnum(row[ix["consensus_iptm"]])
                if m is None or c is None or abs(m - c) > 1e-6:
                    iptm_mismatch += 1

            w.writerow(row)
            if n % 500000 == 0:
                print(f"  {n:,} rows | dockq {dq_applied:,} | iptm {iptm_filled:,}", flush=True)

os.replace(tmp, a.out)
print(f"""
wrote {a.out}
  rows                       {n:,}
  matched v9 sidecar         {hit:,}
  dockq converted            {dq_applied:,}
  dockq already converted    {dq_already:,}   (skipped, re-run safe)
  v9 rows w/o sidecar row    {v9_no_sidecar:,}
  iptm filled                {iptm_filled:,}
  iptm pre-existing audited  {iptm_audit:,}
  iptm consensus != min3     {iptm_mismatch:,}   (expect 0)
""", flush=True)
