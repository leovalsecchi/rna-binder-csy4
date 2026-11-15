#!/usr/bin/env python3
"""build_fullmetrics_chunks_v8_hdtopup.py

Chunk the v8 designs recovered by the 2026-08-28 re-collect that do NOT already have
fullmetrics, in the EXACT format fullmetrics_chunk.sbatch expects:
  cols  = design_id,stack,sequence,af3_msaoff_cif_path,b2_cif_path,esm_cif_path,designed_cif_path
  rows  = 2500 per chunk
  names = chunk_<NNNNN>.csv (5-digit)

Coverage is computed by DIFFING against the design_ids actually present in the existing
fullmetrics_chunks_v8 OUTPUT chunks (chunk_af3_*.csv) -- NOT by assuming the old manifest
was complete. That assumption is exactly what produced the v9 manifest bug (271,596
designs processed twice, 271,604 never processed, while every row-count check passed).

Streaming; dedupes on design_id; refuses a header-only final chunk.
"""
import csv, sys, glob
from pathlib import Path

RB   = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
SRC  = RB / "analyses/midrun_v8_recollect/collected_clean.csv"
DONE = sorted(glob.glob(str(RB / "analyses/cross_campaign/fullmetrics_chunks_v8/chunk_af3_*.csv")))
OUT  = RB / "analyses/cross_campaign/fullmetrics_chunks_v8_hdtopup"
COLS = ["design_id", "stack", "sequence",
        "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path", "designed_cif_path"]
PER  = 2500

csv.field_size_limit(10 * 1024 * 1024)

# ---- what already has fullmetrics -------------------------------------------------
done = set()
for p in DONE:
    with open(p, newline="") as f:
        r = csv.reader(f)
        hdr = next(r, None)
        if not hdr or hdr[0] != "design_id":
            sys.exit(f"FATAL: {p} first col is {hdr[0] if hdr else None!r}, expected design_id")
        for row in r:
            if row:
                done.add(row[0])
print(f"already computed: {len(done):,} unique design_id across {len(DONE)} chunk_af3_*.csv")

# ---- stream the re-collect, emit only the uncovered --------------------------------
OUT.mkdir(parents=True, exist_ok=True)
for stale in OUT.glob("chunk_*.csv"):
    stale.unlink()

seen = set()
n_src = n_out = n_dup = n_skip = nchunk = 0
by_stack = {}
w = fh = None
with SRC.open(newline="") as f:
    r = csv.DictReader(f)
    missing = [c for c in COLS if c not in r.fieldnames]
    if missing:
        sys.exit(f"FATAL: source lacks {missing}")
    for row in r:
        n_src += 1
        did = row["design_id"]
        if did in done:
            n_skip += 1
            continue
        if did in seen:
            n_dup += 1
            continue
        seen.add(did)
        if n_out % PER == 0:
            if fh: fh.close()
            p = OUT / f"chunk_{nchunk:05d}.csv"
            fh = p.open("w", newline=""); w = csv.writer(fh); w.writerow(COLS)
            nchunk += 1
        w.writerow([row[c] for c in COLS])
        by_stack[row["stack"]] = by_stack.get(row["stack"], 0) + 1
        n_out += 1
if fh: fh.close()

tail = n_out - (nchunk - 1) * PER
print(f"source rows           : {n_src:,}")
print(f"  already computed    : {n_skip:,}")
print(f"  duplicate in source : {n_dup:,}")
print(f"  EMITTED             : {n_out:,}  -> {nchunk} chunks in {OUT}")
print(f"last chunk chunk_{nchunk-1:05d}.csv holds {tail} rows")
print("by stack:")
for s, c in sorted(by_stack.items(), key=lambda kv: -kv[1]):
    print(f"  {s:16} {c:,}")

# ---- invariants --------------------------------------------------------------------
assert n_out == len(seen), "emitted != unique emitted"
assert not (seen & done), "emitted overlaps already-computed"
if tail == 0:
    sys.exit("FATAL: header-only final chunk")
print("OK: emitted set is disjoint from the computed set and internally unique")
