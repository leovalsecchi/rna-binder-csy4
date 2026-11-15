#!/usr/bin/env python3
"""hbond_nts_refolds.py — per-design RNA-nucleotide H-bond sets for both MSA-off
refolders → plot sidecar (design_id, refold_condition, hbond_nts).

Wraps code/scripts/hbond_map_pyrosetta.py (PyRosetta protein→RNA H-bonds, inside
pyrosetta.sif) over the af3 + b2 refold structures, then collapses its (aa, nt_pos)
map to the SET of RNA nt positions in >=1 H-bond per (design, refolder). nt_pos uses
the same 1..16 chain-rank numbering as plot 5. Generalizes to the full set the same
way as dockq_refolds.py (chunk the input).

Usage (Euler host):
    python3 analyses/plots_tools/hbond_nts_refolds.py \\
        --input analyses/cross_campaign/combined_sample10k_metrics.csv \\
        --out   analyses/plots_tools/sidecars/hbond_nts_10k.csv  [--jobs 16] [--limit N]
"""
from __future__ import annotations
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
SIF = "/cluster/work/wicky/containers/pyrosetta.sif"
HBMAP = RB / "code" / "scripts" / "hbond_map_pyrosetta.py"
CONDS = [("af3_msa_off", "af3_msaoff_cif_path"), ("b2_msa_off", "b2_cif_path")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    m = pd.read_csv(a.input, low_memory=False)
    parts = []
    for cond, col in CONDS:
        sub = m[["design_id", "stack", col]].dropna(subset=[col]).rename(columns={col: "abs_path"})
        sub["refolder"] = cond
        sub["msa_variant"] = "msaoff"
        parts.append(sub)
    inp = pd.concat(parts, ignore_index=True)

    sidecar_dir = Path(a.out).parent
    sidecar_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(sidecar_dir)) as tmp:
        tmp = Path(tmp)
        in_csv = tmp / "hbmap_in.csv"
        inp.to_csv(in_csv, index=False)
        map_csv = tmp / "hbmap_out.csv"
        cmd = (f"apptainer exec --bind /cluster/work/wicky --bind {tmp} --env MPLCONFIGDIR=/tmp "
               f"{SIF} python3 {HBMAP} --csv {in_csv} --output {map_csv} --path-col abs_path "
               f"--jobs {a.jobs} " + (f"--limit {a.limit} " if a.limit else "") + "--progress")
        print(f"[hbond_nts] {cmd}", file=sys.stderr)
        subprocess.run(cmd, shell=True, check=True, cwd=str(RB))
        mp = pd.read_csv(map_csv, low_memory=False)

    p2d = dict(zip(inp["abs_path"], inp["design_id"]))
    mp["design_id"] = mp["abs_path"].map(p2d)
    mp = mp.dropna(subset=["design_id"])
    agg = (mp.groupby(["design_id", "refolder"])["nt_pos"]
             .apply(lambda s: sorted({int(x) for x in s}))
             .reset_index()
             .rename(columns={"refolder": "refold_condition", "nt_pos": "hbond_nts"}))
    agg.to_csv(a.out, index=False)
    print(f"[hbond_nts] wrote {a.out}: {len(agg)} (design,refolder) rows with >=1 H-bond")


if __name__ == "__main__":
    main()
