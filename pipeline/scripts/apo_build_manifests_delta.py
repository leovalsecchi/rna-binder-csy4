#!/usr/bin/env python3
"""Build apo-fold chunk manifests (id, sequence, esm_ref, stack) from a source metrics CSV.
Round-robin split so every chunk spans all stacks (source is stack-grouped, so round-robin is
essential for balance). Skips rows missing sequence or esm ref. Optionally excludes ids already
folded in a prior apo run via --done-glob (matches prior chunk_*_apo.csv on their `id` column).

Opens all nchunks output files at once (round-robin) -> raises RLIMIT_NOFILE soft cap as needed."""
import argparse
import csv
import glob
import os
import resource
import sys

csv.field_size_limit(sys.maxsize)  # long sequences / paths
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"

ap = argparse.ArgumentParser()
ap.add_argument("--src", default=f"{RB}/analyses/cross_campaign/combined_clean_v2v7_filled_incl_v2rfd3.csv")
ap.add_argument("--out", default=f"{RB}/apo_chunks_v2v7")
ap.add_argument("-n", "--nchunks", type=int, default=3835)
ap.add_argument("--done-glob", default="", help="glob of prior *_apo.csv whose `id` column marks already-folded designs to skip")
a = ap.parse_args()

# --- raise open-file soft limit so we can hold nchunks writers open at once ---
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
want = a.nchunks + 256
if soft < want:
    new_soft = min(hard, max(want, 8192))
    resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))
    print(f"raised RLIMIT_NOFILE soft {soft} -> {new_soft} (hard {hard}) for {a.nchunks} chunks")

# --- load already-done ids (optional) ---
done = set()
if a.done_glob:
    files = sorted(glob.glob(a.done_glob))
    for fp in files:
        with open(fp, newline="") as fh:
            r = csv.DictReader(fh)
            col = "id" if (r.fieldnames and "id" in r.fieldnames) else ("design_id" if (r.fieldnames and "design_id" in r.fieldnames) else None)
            if not col:
                continue
            for d in r:
                v = (d.get(col) or "").strip()
                if v:
                    done.add(v)
    print(f"done-glob: {len(files)} files -> {len(done)} already-folded ids to skip")

os.makedirs(a.out, exist_ok=True)
writers = [open(f"{a.out}/chunk_{i:04d}.tsv", "w", newline="") for i in range(a.nchunks)]
ws = [csv.writer(w, delimiter="\t") for w in writers]
for w in ws:
    w.writerow(["id", "sequence", "esm_ref", "stack"])

total = kept = skip_missing = skip_done = 0
with open(a.src, newline="") as fh:
    for d in csv.DictReader(fh):
        total += 1
        did = (d.get("design_id") or "").strip()
        seq = (d.get("sequence") or "").strip()
        esm = (d.get("esm_cif_path") or "").strip()
        if not seq or not esm:
            skip_missing += 1
            continue
        if did in done:
            skip_done += 1
            continue
        ws[kept % a.nchunks].writerow([did, seq, esm, d.get("stack", "")])
        kept += 1
for w in writers:
    w.close()
print(f"rows={total} kept={kept} skip_missing={skip_missing} skip_done={skip_done} "
      f"-> {a.nchunks} chunks in {a.out} (~{kept // a.nchunks}/chunk)")
