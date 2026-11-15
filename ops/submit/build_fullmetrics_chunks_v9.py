#!/usr/bin/env python3
"""build_fullmetrics_chunks_v9.py — chunk midrun_v9 collected_clean.csv for the
fullmetrics suite, in the EXACT format of fullmetrics_chunks/ and fullmetrics_chunks_v8/:
  cols  = design_id,stack,sequence,af3_msaoff_cif_path,b2_cif_path,esm_cif_path,designed_cif_path
  rows  = 2500 per chunk
  names = chunk_<NNNNN>.csv (5-digit)
Streaming (never loads the 1.1 GB master). Refuses to emit a header-only final chunk,
which would IndexError the task and cancel an afterok merge.
"""
import csv, sys
from pathlib import Path

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
SRC = RB / "analyses/midrun_v9/collected_clean.csv"
OUT = RB / "analyses/cross_campaign/fullmetrics_chunks_v9"
COLS = ["design_id", "stack", "sequence",
        "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path", "designed_cif_path"]
PER = 2500

csv.field_size_limit(10 * 1024 * 1024)
OUT.mkdir(parents=True, exist_ok=True)

n = nchunk = 0
w = fh = None
with SRC.open(newline="") as f:
    r = csv.DictReader(f)
    missing = [c for c in COLS if c not in r.fieldnames]
    if missing:
        sys.exit(f"FATAL: source lacks {missing}")
    for row in r:
        if n % PER == 0:
            if fh: fh.close()
            p = OUT / f"chunk_{nchunk:05d}.csv"
            fh = p.open("w", newline=""); w = csv.writer(fh); w.writerow(COLS)
            nchunk += 1
        w.writerow([row[c] for c in COLS])
        n += 1
if fh: fh.close()

tail = n - (nchunk - 1) * PER
print(f"wrote {n} designs -> {nchunk} chunks in {OUT}")
print(f"last chunk chunk_{nchunk-1:05d}.csv holds {tail} rows")
if tail == 0:
    sys.exit("FATAL: header-only final chunk")
