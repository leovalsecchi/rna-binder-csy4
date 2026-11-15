#!/usr/bin/env python3
"""Quantify duplication in combined_clean_v2v7.csv: design_id vs sequence,
within-version vs cross-version."""
import sys
import pandas as pd

path = sys.argv[1] if len(sys.argv) > 1 else \
    "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/combined_clean_v2v7.csv"
d = pd.read_csv(path, usecols=["design_id", "sequence", "campaign"], low_memory=False)
n = len(d)
print(f"file: {path}")
print(f"total rows                : {n}")
print(f"unique design_id          : {d['design_id'].nunique()}   (dup design_ids: {n - d['design_id'].nunique()})")

seq = d["sequence"].dropna()
n_seq = len(seq)
uniq_global = seq.nunique()
print(f"rows with a sequence      : {n_seq}   (null-seq rows: {n - n_seq})")
print(f"unique sequences (global) : {uniq_global}")
print(f"sequence-duplicate rows   : {n_seq - uniq_global}   ({100*(n_seq-uniq_global)/n_seq:.2f}% of seq rows)")

# within-version dups (should be ~0 since each version was cleaned/deduped internally)
wv = d.dropna(subset=["sequence"]).duplicated(subset=["campaign", "sequence"]).sum()
print(f"within-version dup rows   : {wv}   (expected ~0; each version deduped internally)")
print(f"=> cross-version dup rows : {(n_seq - uniq_global) - wv}")

# how many distinct sequences appear in >1 version
vc = d.dropna(subset=["sequence"]).groupby("sequence")["campaign"].nunique()
print(f"sequences in >1 version   : {(vc > 1).sum()}")
print("\nper-version rows:")
print(d["campaign"].value_counts().sort_index().to_string())
