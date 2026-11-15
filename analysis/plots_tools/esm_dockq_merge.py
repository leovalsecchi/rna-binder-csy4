#!/usr/bin/env python3
"""Compute DockQ(ESM refold vs designed) over the 10k sample and merge the
dockq_esm_msa_off column into the existing dockq_10k sidecar. Mirrors
dockq_refolds.run_condition exactly (same dockq_self.py wrapper)."""
import argparse, subprocess, tempfile
from pathlib import Path
import pandas as pd

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
DOCKQ_SELF = RB / "code" / "scripts" / "dockq_self.py"
INP = RB / "analyses" / "cross_campaign" / "combined_sample10k_metrics.csv"
SIDECAR = RB / "analyses" / "plots_tools" / "sidecars" / "dockq_10k.csv"
PATH_COL, REF_COL, OUTCOL = "esm_cif_path", "designed_cif_path", "dockq_esm_msa_off"

ap = argparse.ArgumentParser()
ap.add_argument("--jobs", type=int, default=16)
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--smoke", action="store_true", help="print result, do NOT merge into sidecar")
a = ap.parse_args()

df = pd.read_csv(INP, low_memory=False)[["design_id", PATH_COL, REF_COL]].dropna(subset=[PATH_COL, REF_COL])
print(f"{len(df)} designs with esm_cif_path + designed_cif_path")
with tempfile.TemporaryDirectory(dir=str(SIDECAR.parent)) as tmp:
    tmp = Path(tmp); slim = tmp / "in.csv"; out = tmp / "out.csv"
    df.to_csv(slim, index=False)
    env = (f"PYTHONPATH={RB}/code/pylibs_dockq:{RB}/code/scratch_dockq "
           f"PATH={RB}/code/bin:$PATH")
    cmd = (f"{env} python3 {DOCKQ_SELF} --csv {slim} --output {out} "
           f"--path-col {PATH_COL} --ref-col {REF_COL} --jobs {a.jobs} "
           + (f"--limit {a.limit} " if a.limit else "") + "--progress")
    subprocess.run(cmd, shell=True, check=True, cwd=str(RB))
    r = pd.read_csv(out, low_memory=False)[["design_id", "dockq_self"]].drop_duplicates("design_id")
    r = r.rename(columns={"dockq_self": OUTCOL})

print(f"computed {OUTCOL}: {int(r[OUTCOL].notna().sum())} non-null; "
      f"min={r[OUTCOL].min():.3f} med={r[OUTCOL].median():.3f} max={r[OUTCOL].max():.3f}")
if a.smoke:
    print("SMOKE: not merging."); raise SystemExit(0)

sc = pd.read_csv(SIDECAR, low_memory=False)
if OUTCOL in sc.columns:
    sc = sc.drop(columns=[OUTCOL])
sc = sc.merge(r, on="design_id", how="left")
sc.to_csv(SIDECAR, index=False)
print(f"merged -> {SIDECAR}: esm non-null {int(sc[OUTCOL].notna().sum())}/{len(sc)} "
      f"(af3={int(sc['dockq_af3_msa_off'].notna().sum())}, b2={int(sc['dockq_b2_msa_off'].notna().sum())})")
