#!/usr/bin/env python3
"""dockq_refolds_full.py — DockQ per refold backend, keeping ALL FOUR components.

Same job as analyses/plots_tools/dockq_refolds.py, with the one difference that matters:
that script returns only `dockq_self` and drops dockq_LRMSD / dockq_iRMSD / dockq_fnat into a
TemporaryDirectory it then deletes. Those three are not optional here —

  * dockq_LRMSD_<be>  IS the filter tool's sc_RMSD axis, and
  * dockq_rna_recombine.py needs fnat + iRMSD + LRMSD to rescore DockQ in the RNA-aligned
    frame, which is the convention the master has used since 2026-08-04.

Emits design_id + {dockq_self,dockq_LRMSD,dockq_iRMSD,dockq_fnat}_{af3,b2,esm}. b2 is
computed because dockq_rna_recombine.py gates its consensus on all three backends; the filter
tool simply never maps the b2 columns.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
DOCKQ_SELF = RB / "code" / "scripts" / "dockq_self.py"
REF_COL = "designed_cif_path"
COMPONENTS = ["dockq_self", "dockq_LRMSD", "dockq_iRMSD", "dockq_fnat"]
CONDITIONS = [("af3", "af3_msaoff_cif_path"),
              ("b2", "b2_cif_path"),
              ("esm", "esm_cif_path")]


def run_condition(inp: Path, path_col: str, be: str, jobs: int, workdir: Path) -> pd.DataFrame:
    slim = workdir / f"in_{be}.csv"
    df = pd.read_csv(inp, low_memory=False)[["design_id", path_col, REF_COL]]
    df = df.dropna(subset=[path_col, REF_COL])
    df.to_csv(slim, index=False)
    out = workdir / f"out_{be}.csv"
    env = (f"PYTHONPATH={RB}/code/pylibs_dockq:{RB}/code/scratch_dockq "
           f"PATH={RB}/code/bin:$PATH")
    cmd = (f"{env} python3 {DOCKQ_SELF} --csv {slim} --output {out} "
           f"--path-col {path_col} --ref-col {REF_COL} --jobs {jobs} --progress")
    print(f"[dockq_full] {be}: {cmd}", file=sys.stderr, flush=True)
    subprocess.run(cmd, shell=True, check=True, cwd=str(RB))
    r = pd.read_csv(out, low_memory=False).drop_duplicates("design_id")
    have = [c for c in COMPONENTS if c in r.columns]
    missing = set(COMPONENTS) - set(have)
    if missing:
        sys.exit(f"FATAL: {be} output lacks {sorted(missing)} — got {list(r.columns)}")
    return r[["design_id"] + have].rename(columns={c: f"{c}_{be}" for c in have})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()

    base = pd.read_csv(a.input, low_memory=False)[["design_id"]].drop_duplicates()
    n0 = len(base)
    outdir = Path(a.out).parent
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(outdir)) as tmp:
        for be, path_col in CONDITIONS:
            base = base.merge(run_condition(Path(a.input), path_col, be, a.jobs, Path(tmp)),
                              on="design_id", how="left")
    if len(base) != n0:
        sys.exit(f"FATAL: row count changed {n0} -> {len(base)} (merge fanout)")
    base.to_csv(a.out, index=False)
    filled = {be: int(base[f"dockq_self_{be}"].notna().sum()) for be, _ in CONDITIONS}
    print(f"wrote {len(base)} rows x {len(base.columns)} cols -> {a.out}  filled={filled}",
          flush=True)


if __name__ == "__main__":
    main()
