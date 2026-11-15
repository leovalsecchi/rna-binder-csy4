#!/usr/bin/env python3
"""Add seq_max_aa_frac + seq_longest_run columns to a metrics CSV (from its
`sequence` column). Idempotent; edits in place with a .bak_pre_seq backup.

    python3 add_seq_metrics.py <csv> [--seq-col sequence] [--no-backup]
"""
import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
try:
    from rna_binder.analysis.seq_complexity import seq_complexity, SEQ_COLUMNS
except ModuleNotFoundError:  # deployed layout: src alongside scripts under code/
    sys.path.insert(0, "/cluster/work/wicky/lvalsecchi/rna_binder/code/src")
    from rna_binder.analysis.seq_complexity import seq_complexity, SEQ_COLUMNS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--seq-col", default="sequence")
    ap.add_argument("--no-backup", action="store_true")
    a = ap.parse_args()

    df = pd.read_csv(a.csv, low_memory=False)
    if a.seq_col not in df.columns:
        sys.exit(f"no '{a.seq_col}' column in {a.csv}: {list(df.columns)[:8]}...")
    vals = [seq_complexity(s) for s in df[a.seq_col]]
    for col in SEQ_COLUMNS:
        df[col] = [v[col] for v in vals]  # overwrite if present (idempotent)

    if not a.no_backup:
        shutil.copy(a.csv, str(a.csv) + ".bak_pre_seq")
    df.to_csv(a.csv, index=False)
    n = df["seq_max_aa_frac"].notna().sum()
    flagged = int(((df["seq_max_aa_frac"] > 0.35) | (df["seq_longest_run"] >= 6)).sum())
    print(f"{a.csv.name}: {len(df)} rows, seq scored {n}, "
          f"low-complexity (max_aa>0.35 or run>=6): {flagged} ({100*flagged/max(len(df),1):.1f}%)")


if __name__ == "__main__":
    main()
