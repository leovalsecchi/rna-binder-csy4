"""Union manifest of every design we have a refold CIF for: v2v7 fullmetrics table
+ the v8 fullmetrics input-chunk series, deduped on design_id."""
import glob
import sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder/"
CC = RB + "analyses/cross_campaign/"
BE = sys.argv[1] if len(sys.argv) > 1 else "af3"
PATHCOL = {"af3": "af3_msaoff_cif_path", "b2": "b2_cif_path",
           "esm": "esm_cif_path"}[BE]
COLS = ["design_id", "stack", PATHCOL]

parts, n = [], 0
for ch in pd.read_csv(CC + "combined_clean_v2v7_fullmetrics.csv", usecols=COLS,
                      chunksize=500_000, low_memory=False):
    parts.append(ch)
    n += len(ch)
print("v2v7 rows: %d" % n, flush=True)

v8 = sorted(glob.glob(CC + "fullmetrics_chunks_v8/chunk_[0-9]*.csv"))
print("v8 af3-series chunk files: %d" % len(v8), flush=True)
m = 0
for f in v8:
    try:
        c = pd.read_csv(f, usecols=COLS, low_memory=False)
    except ValueError:
        continue
    parts.append(c)
    m += len(c)
print("v8 rows: %d" % m, flush=True)

d = pd.concat(parts, ignore_index=True)
print("concat: %d" % len(d), flush=True)
d = d.dropna(subset=[PATHCOL])
print("with %s path: %d" % (BE, len(d)), flush=True)
d = d.drop_duplicates("design_id")
print("unique design_id: %d" % len(d), flush=True)
out = CC + ("opth_manifest_full_%s.csv" % BE)
d.to_csv(out, index=False)
print("wrote %s" % out, flush=True)
print(d["stack"].value_counts().to_string(), flush=True)
