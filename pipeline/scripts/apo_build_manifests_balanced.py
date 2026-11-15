#!/usr/bin/env python3
"""Build apo-fold chunk manifests with EQUAL NUMBERS PER STACK IN EVERY CHUNK.

Differs from apo_build_manifests{,_delta}.py in one way that matters: those
round-robin by ROW INDEX (`ws[kept % nchunks]`), which on a stack-grouped source
gives every chunk the same stack PROPORTIONS as the campaign (v9: 29% BG+BIF vs
5% HD). This one interleaves ACROSS STACKS and cuts chunks sequentially, so each
chunk holds ~equal COUNTS per stack for as long as every stack still has supply.

Consequence: any set of finished chunks is stack-balanced, regardless of the
order tasks happen to complete in — you do not have to rely on chunk ordering.
Once the smallest stack is exhausted the interleave continues over the stacks
that remain, so the balanced phase is min(stack)*n_stacks designs long.

    python3 apo_build_manifests_balanced.py --src <clean.csv> --out apo_chunks_v9 \
        -n 1193 [--done-glob '...']
"""
import argparse
import collections
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"

ap = argparse.ArgumentParser()
ap.add_argument("--src", default=f"{RB}/analyses/midrun_v9/collected_clean.csv")
ap.add_argument("--out", default=f"{RB}/apo_chunks_v9")
ap.add_argument("-n", "--nchunks", type=int, default=1193)
ap.add_argument("--done-glob", default="")
a = ap.parse_args()

done = set()
if a.done_glob:
    for fp in sorted(glob.glob(a.done_glob)):
        with open(fp, newline="") as fh:
            r = csv.DictReader(fh)
            col = "id" if (r.fieldnames and "id" in r.fieldnames) else None
            if not col:
                continue
            for d in r:
                v = (d.get(col) or "").strip()
                if v:
                    done.add(v)
    print(f"done-glob -> {len(done)} ids to skip", flush=True)

by_stack = collections.OrderedDict()
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
        by_stack.setdefault(d.get("stack", ""), []).append((did, seq, esm))
        kept += 1

print(f"rows={total} kept={kept} skip_missing={skip_missing} skip_done={skip_done}", flush=True)
for s, v in sorted(by_stack.items(), key=lambda kv: -len(kv[1])):
    print(f"  {s:<14} {len(v):>9,}", flush=True)

# interleave across stacks: one row from each stack in rotation, dropping stacks
# as they run out. Largest-first ordering within a round is irrelevant to balance.
order = sorted(by_stack, key=lambda s: -len(by_stack[s]))
pos = {s: 0 for s in order}
stream = []
live = list(order)
while live:
    for s in list(live):
        i = pos[s]
        if i >= len(by_stack[s]):
            live.remove(s)
            continue
        stream.append((s, by_stack[s][i]))
        pos[s] = i + 1
assert len(stream) == kept

nmin = min(len(v) for v in by_stack.values())
print(f"\nbalanced phase: {nmin:,} per stack x {len(by_stack)} stacks = "
      f"{nmin * len(by_stack):,} designs", flush=True)

os.makedirs(a.out, exist_ok=True)
per = -(-kept // a.nchunks)          # ceil
for i in range(a.nchunks):
    rows = stream[i * per:(i + 1) * per]
    with open(f"{a.out}/chunk_{i:04d}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "sequence", "esm_ref", "stack"])
        for s, (did, seq, esm) in rows:
            w.writerow([did, seq, esm, s])

c0 = collections.Counter(s for s, _ in stream[:per])
print(f"\n{a.nchunks} chunks of ~{per} in {a.out}", flush=True)
print("chunk_0000 stack mix: " + "  ".join(f"{k}={v}" for k, v in sorted(c0.items())), flush=True)
print(f"balanced through chunk ~{(nmin * len(by_stack)) // per}", flush=True)
