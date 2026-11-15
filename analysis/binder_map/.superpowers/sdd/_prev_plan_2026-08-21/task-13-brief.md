### Task 13: Stage 1 — extraction and fingerprinting

The heaviest stage, and the one that protects the inode quota. One pass over every archive produces both the Foldseek input and the epitope fingerprints.

**Files:**
- Create: `scripts/stage1_extract.py`
- Create: `slurm/stage1_extract.sbatch`

**Interfaces:**
- Consumes: everything from Tasks 2–6.
- Produces, under `$FM/results/<layer>/`:
  - `shards/shard_%02d.tar.gz` — backbone PDBs, member name = `design_id.pdb`.
  - `hashes/shard_%02d.tsv` — `design_id, coord_hash, protein_chain, n_res`.
  - `fp/shard_%02d.parquet` — `design_id` + 48 fingerprint columns + `per_nt_1..16` + `contacts_offframe` + `contacts_total` + `cx, cy, cz, ax, ay, az`.
  - `missing/shard_%02d.tsv` — `design_id, reason`.

- [ ] **Step 1: Implement `scripts/stage1_extract.py`**

```python
#!/usr/bin/env python3
"""Stage 1: extract backbones and compute epitope fingerprints.

One pass over the archives does both jobs, because opening 6,780 .tar.zst
files twice would double the most expensive part of the pipeline for no gain.

Everything is streamed through $TMPDIR into a single .tar.gz shard. No loose
per-structure files are ever written: the project has a hard inode cap, and
215k loose files would breach it.

Usage (one array task):
    python scripts/stage1_extract.py --layer esm --shard 7 --n-shards 64
"""
from __future__ import annotations

import argparse
import io
import os
import subprocess
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map import RNA_LEN
from binder_map.backbone import coord_hash, write_backbone_pdb
from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.fingerprint import FINGERPRINT_COLUMNS, fingerprint
from binder_map.geometry import approach
from binder_map.rnaframe import FrameError

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
FM = RB / "analyses" / "foldseek_map"


def shard_archives(manifest: pd.DataFrame, shard: int, n_shards: int):
    """Assign whole archives to shards, so each archive is opened exactly once."""
    archives = sorted(manifest["archive_path"].unique())
    mine = set(archives[shard::n_shards])
    return manifest[manifest["archive_path"].isin(mine)]


def read_archive_members(archive: Path, members: list[str]) -> dict[str, bytes]:
    """Extract the needed members from one .tar.zst in a single tar invocation.

    Python's tarfile has no zstd codec here, so tar streams the whole archive
    once and we sieve it in memory.
    """
    proc = subprocess.Popen(
        ["tar", "--zstd", "-xOf", str(archive)] + members,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    out, err = proc.communicate()
    if proc.returncode != 0:
        return {}
    # tar concatenates members in the order given, so split on the mmCIF
    # data_ header, which every dialect here starts with.
    chunks: dict[str, bytes] = {}
    parts = out.split(b"\ndata_")
    if len(parts) == 1:
        blobs = [out]
    else:
        blobs = [parts[0]] + [b"data_" + p for p in parts[1:]]
    for member, blob in zip(members, blobs):
        chunks[member] = blob
    return chunks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, default=64)
    parser.add_argument("--build", type=Path, default=FM / "build")
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    manifest = pd.read_csv(args.build / "manifests" / f"{args.layer}.tsv", sep="\t")
    mine = shard_archives(manifest, args.shard, args.n_shards)
    print(f"layer={args.layer} shard={args.shard}/{args.n_shards}: "
          f"{len(mine):,d} structures from {mine['archive_path'].nunique():,d} archives")

    out_dir = args.results / args.layer
    for sub in ("shards", "hashes", "fp", "missing"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"bm_{args.layer}_{args.shard}"
    tmp.mkdir(parents=True, exist_ok=True)
    scratch = tmp / "one.cif"

    shard_name = f"shard_{args.shard:02d}"
    tar_path = out_dir / "shards" / f"{shard_name}.tar.gz"

    hashes, fps, missing = [], [], []
    seen_local: set[str] = set()
    t0 = time.time()

    with tarfile.open(tar_path, "w:gz") as tar:
        for archive_path, group in mine.groupby("archive_path", sort=True):
            archive = RB / archive_path
            if not archive.exists():
                for design_id in group["design_id"]:
                    missing.append({"design_id": design_id,
                                    "reason": "archive_missing"})
                continue

            members = group["member_path"].tolist()
            blobs = read_archive_members(archive, members)

            for _, row in group.iterrows():
                design_id = row["design_id"]
                blob = blobs.get(row["member_path"])
                if not blob or not blob.strip():
                    missing.append({"design_id": design_id,
                                    "reason": "member_missing"})
                    continue
                try:
                    scratch.write_bytes(blob)
                    atoms = load_structure(scratch)
                    # One composition scan, not two: split_chains already
                    # determined which chain is protein, and this loop runs
                    # ~215k times per layer.
                    protein, rna = split_chains(atoms)
                    chain = str(protein.chain_id[0])
                    n_res = len(set(zip(protein.chain_id, protein.res_id)))
                    if n_res != 80:
                        missing.append({"design_id": design_id,
                                        "reason": f"chain_select_failed:n_res={n_res}"})
                        continue
                    digest = coord_hash(protein)
                    fp = fingerprint(protein, rna)
                    app = approach(protein, rna)
                except (ChainSelectError, FrameError) as exc:
                    missing.append({"design_id": design_id,
                                    "reason": f"chain_select_failed:{type(exc).__name__}"})
                    continue
                except Exception as exc:  # noqa: BLE001 - one bad file must not kill the shard
                    missing.append({"design_id": design_id,
                                    "reason": f"error:{type(exc).__name__}"})
                    continue

                hashes.append({"design_id": design_id, "coord_hash": digest,
                               "protein_chain": chain, "n_res": n_res})

                record = {"design_id": design_id}
                record.update(dict(zip(FINGERPRINT_COLUMNS, fp.vector.tolist())))
                record.update({f"per_nt_{i+1}": int(v)
                               for i, v in enumerate(fp.per_nt)})
                record.update({
                    "contacts_offframe": fp.offframe,
                    "contacts_total": fp.total,
                    "cx": float(app.centroid[0]), "cy": float(app.centroid[1]),
                    "cz": float(app.centroid[2]),
                    "ax": float(app.axis[0]), "ay": float(app.axis[1]),
                    "az": float(app.axis[2]),
                })
                fps.append(record)

                # Local dedup: within a shard, seq variants of one backbone
                # live in the same archive, so most redundancy dies here.
                if digest in seen_local:
                    continue
                seen_local.add(digest)
                pdb = write_backbone_pdb(protein, design_id)
                info = tarfile.TarInfo(name=f"{design_id}.pdb")
                info.size = len(pdb)
                info.mtime = 0
                tar.addfile(info, io.BytesIO(pdb))

    pd.DataFrame(hashes).to_csv(out_dir / "hashes" / f"{shard_name}.tsv",
                                sep="\t", index=False)
    pd.DataFrame(fps).to_parquet(out_dir / "fp" / f"{shard_name}.parquet",
                                 index=False)
    pd.DataFrame(missing, columns=["design_id", "reason"]).to_csv(
        out_dir / "missing" / f"{shard_name}.tsv", sep="\t", index=False)

    for leftover in tmp.glob("*"):
        leftover.unlink(missing_ok=True)

    print(f"done in {time.time() - t0:.0f}s: {len(hashes):,d} parsed, "
          f"{len(seen_local):,d} unique backbones, {len(missing):,d} missing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Write `slurm/stage1_extract.sbatch`**

```bash
#!/bin/bash
# Stage 1: extract backbones + fingerprints. One array task per shard.
#   sbatch --job-name=bm_s1_esm --array=0-63%16 slurm/stage1_extract.sbatch esm
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

LAYER="${1:?usage: stage1_extract.sbatch <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

cd "$FM" || exit 1
echo "=== stage1 layer=$LAYER shard=$SLURM_ARRAY_TASK_ID on $(hostname) $(date) ==="

"$PY" scripts/stage1_extract.py \
  --layer "$LAYER" \
  --shard "$SLURM_ARRAY_TASK_ID" \
  --n-shards "${N_SHARDS:-64}" \
  --build "$FM/${BUILD:-build}" \
  --results "$FM/${RESULTS:-results}"
rc=$?
echo "=== exit $rc $(date) ==="
exit $rc
```

- [ ] **Step 3: Sync and smoke-test one shard interactively**

Never launch 192 array tasks on untested code.

```bash
cd /Users/leo/Documents/binder_map
rsync -av --exclude '.git' --exclude 'tests/fixtures' --exclude 'build' \
  ./ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/

ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  TMPDIR=/tmp PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/stage1_extract.py --layer esm --shard 0 --n-shards 512'
```

Expected: a few hundred structures parsed, `0 missing` or close to it, and a `shard_00.tar.gz` written. `--n-shards 512` keeps this single interactive task small.

- [ ] **Step 4: Verify the shard contents and the inode footprint**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results/esm
echo "=== shard members ==="; tar -tzf shards/shard_00.tar.gz | head -3
echo "=== count ==="; tar -tzf shards/shard_00.tar.gz | wc -l
echo "=== inodes used by results/ ==="; find . | wc -l
echo "=== protein chain tally (must be A for esm) ==="; cut -f3 hashes/shard_00.tsv | sort | uniq -c
echo "=== fingerprint sanity ==="
PYTHONPATH=../../src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python - <<PY
import pandas as pd
df = pd.read_parquet("fp/shard_00.parquet")
print("rows", len(df), "cols", len(df.columns))
print("contacts_total: min", df.contacts_total.min(), "median", df.contacts_total.median())
print("zero-contact rows:", int((df.contacts_total == 0).sum()))
print("offframe nonzero:", int((df.contacts_offframe > 0).sum()))
PY'
```

Expected: members named `<design_id>.pdb`; the inode count in the tens, not thousands; every `protein_chain` is `A` for the esm layer; median `contacts_total` comfortably above zero.

- [ ] **Step 5: Verify the design layer picks up the RFD3 inversion**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  TMPDIR=/tmp PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/stage1_extract.py --layer design --shard 0 --n-shards 512 && \
  join -t"	" -1 1 -2 1 \
    <(sort -k1,1 results/design/hashes/shard_00.tsv) \
    <(sort -k1,1 -t"	" build/manifests/design.tsv | cut -f2,4) \
  2>/dev/null | awk -F"\t" "{print \$3, \$5}" | sort | uniq -c'
```

Expected: `RFD3+LiMPNN` rows report protein chain `B`; all other stacks report `A`. If RFD3 reports `A`, chain selection is broken — stop.

- [ ] **Step 6: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage1_extract.py slurm/stage1_extract.sbatch
git commit -m "Add stage 1: single-pass backbone extraction and epitope fingerprinting"
```

---

