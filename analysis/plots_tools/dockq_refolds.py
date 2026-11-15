#!/usr/bin/env python3
"""dockq_refolds.py — DockQ(refold vs designed) for BOTH msa-off refold conditions
(AF3 msa-off + Boltz-2) over a metrics CSV → plot sidecar
(design_id, dockq_af3_msa_off, dockq_b2_msa_off).

Wraps the proven archive-aware, parallel code/scripts/dockq_self.py (host py3.12 +
code/pylibs_dockq) once per condition. GENERALIZES to the full ~900k set: point
--input at a chunk of combined_clean_metrics.csv and fan out with the same array
pattern as full_metrics_parallel.sh; here it's run once on the 3k sample.

Usage (Euler host, NOT a container):
    python3 analyses/plots/dockq_refolds.py \\
        --input analyses/cross_campaign/combined_sample10k_metrics.csv \\
        --out   analyses/plots/sidecars/dockq_10k.csv  [--jobs 16] [--limit N]
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
# (refold_condition, path column in input, output DockQ column)
CONDITIONS = [
    ("af3_msa_off", "af3_msaoff_cif_path", "dockq_af3_msa_off"),
    ("b2_msa_off",  "b2_cif_path",         "dockq_b2_msa_off"),
    ("esm_msa_off", "esm_cif_path",        "dockq_esm_msa_off"),
]


def run_condition(inp: Path, path_col: str, jobs: int, limit: int, workdir: Path) -> pd.DataFrame:
    """Run dockq_self.py for one refold path-col; return df[design_id, dockq_self]."""
    slim = workdir / f"in_{path_col}.csv"
    df = pd.read_csv(inp, low_memory=False)[["design_id", path_col, REF_COL]]
    df = df.dropna(subset=[path_col, REF_COL])
    df.to_csv(slim, index=False)
    out = workdir / f"out_{path_col}.csv"
    env = (f"PYTHONPATH={RB}/code/pylibs_dockq:{RB}/code/scratch_dockq "
           f"PATH={RB}/code/bin:$PATH")
    cmd = (f"{env} python3 {DOCKQ_SELF} --csv {slim} --output {out} "
           f"--path-col {path_col} --ref-col {REF_COL} --jobs {jobs} "
           + (f"--limit {limit} " if limit else "") + "--progress")
    print(f"[dockq_refolds] {path_col}: {cmd}", file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True, cwd=str(RB))
    r = pd.read_csv(out, low_memory=False)
    return r[["design_id", "dockq_self"]].drop_duplicates("design_id")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    base = pd.read_csv(a.input, low_memory=False)[["design_id"]].drop_duplicates()
    sidecar_dir = Path(a.out).parent
    sidecar_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(sidecar_dir)) as tmp:
        tmp = Path(tmp)
        for _cond, path_col, outcol in CONDITIONS:
            got = run_condition(Path(a.input), path_col, a.jobs, a.limit, tmp)
            base = base.merge(got.rename(columns={"dockq_self": outcol}),
                              on="design_id", how="left")
    base.to_csv(a.out, index=False)
    n_af3 = int(base["dockq_af3_msa_off"].notna().sum())
    n_b2 = int(base["dockq_b2_msa_off"].notna().sum())
    print(f"[dockq_refolds] wrote {a.out}: {len(base)} designs — af3={n_af3}, b2={n_b2} with DockQ")


if __name__ == "__main__":
    main()
