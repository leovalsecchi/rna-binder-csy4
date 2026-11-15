### Task 4: Stage 1 rewrite — one streaming pass per archive

**Files:**
- Create: `scripts/stage1_extract_nested.py`
- Create: `slurm/stage1_extract_nested.sbatch`
- Test: `tests/test_stage1_nested.py`

**Interfaces:**
- Consumes: `iter_members` (Task 2), the SHARED manifest from Task 3
  (`build_shared/manifests/design.tsv`, ~508,100 rows = map 1 ∪ map 3)
- Also MODIFIES `src/binder_map/nested.py` — see Step 0 below.

- [ ] **Step 0: Stop streaming once every wanted inner archive is consumed**

Measured 2026-09-14: **design members occupy only the first 8-10% of each outer
archive**; everything after is refold output this pipeline never reads.
Per the `.list` sidecars (tar member order): rfd3 midrun_v2 designs end at 8% of
11,694 members, boltzgen__boltzif midrun_v8 at 8% of 21,906, halludesign
midrun_v9 at 10% of 34,280. Streaming to the end therefore throws away ~90% of
the I/O — 683 GB read where ~68 GB would do.

`iter_members` currently runs the outer loop to exhaustion. Track what is left
and return as soon as it is empty. Returning triggers the existing `finally`,
which closes the pipe and reaps `zstd`; the lifetime handling is already correct
for early exit (verified in Task 2's review), so this needs no new teardown.

In `src/binder_map/nested.py`, inside `iter_members`, replace the loop body's
guard and add the exit:

```python
    remaining = set(wanted)
    ...
            for info in outer_tar:
                if not info.isfile() or info.name not in remaining:
                    continue
                need = wanted[info.name]
                raw = _decompress(outer_tar.extractfile(info).read())
                with tarfile.open(fileobj=io.BytesIO(raw), mode="r|") as inner_tar:
                    for member in inner_tar:
                        if member.isfile() and member.name in need:
                            yield info.name, member.name, inner_tar.extractfile(member).read()
                del raw
                remaining.discard(info.name)
                if not remaining:
                    return          # closes the pipe via `finally`; ~90% less I/O
```

Add to `tests/test_nested.py`:

```python
def test_iter_members_stops_once_every_wanted_archive_is_consumed(tmp_path):
    outer = _build_nested(tmp_path)
    # ask only for the FIRST inner archive; the second must never be opened
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d0.cif"}}
    seen = list(iter_members(outer, wanted))
    assert [(i, m) for i, m, _ in seen] == [("midrun_v9/designs/aaaa.tar.zst", "d0.cif")]


def test_iter_members_still_returns_everything_when_all_are_wanted(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d0.cif", "d1.cif"},
              "midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    got = {(i, m) for i, m, _ in iter_members(outer, wanted)}
    assert len(got) == 3
```

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v`
Expected: PASS, 7 passed (5 existing + 2 new). The existing five must stay green —
this changes when the loop stops, never what it yields.
- Produces: per-shard `shards/shard_NN.tar.gz`, `hashes/shard_NN.tsv`
  (`design_id, coord_hash, protein_chain, n_res`), `fp/shard_NN.parquet`,
  `missing/shard_NN.tsv` — the same file contract stage 2 already consumes.

- [ ] **Step 1: Write the failing test**

Create `tests/test_stage1_nested.py`:

```python
import pandas as pd
import pytest

from stage1_extract_nested import plan_archive_work, shard_outers


def _manifest() -> pd.DataFrame:
    return pd.DataFrame({
        "design_id": ["d0", "d1", "d2", "d3"],
        "outer_path": ["a.tar.zst", "a.tar.zst", "a.tar.zst", "b.tar.zst"],
        "inner_path": ["m/designs/u0.tar.zst", "m/designs/u0.tar.zst",
                       "m/designs/u1.tar.zst", "m/designs/u2.tar.zst"],
        "member_path": ["d0.cif", "d1.cif", "d2.cif", "d3.cif"],
    })


def test_shard_outers_assigns_whole_archives():
    got = shard_outers(_manifest(), shard=0, n_shards=2)
    assert set(got["outer_path"]) == {"a.tar.zst"}
    assert len(got) == 3


def test_every_archive_lands_in_exactly_one_shard():
    manifest = _manifest()
    seen = []
    for shard in range(2):
        seen += list(shard_outers(manifest, shard, 2)["design_id"])
    assert sorted(seen) == ["d0", "d1", "d2", "d3"]


def test_plan_archive_work_groups_members_by_inner_archive():
    wanted = plan_archive_work(_manifest().query("outer_path == 'a.tar.zst'"))
    assert wanted == {
        "m/designs/u0.tar.zst": {"d0.cif", "d1.cif"},
        "m/designs/u1.tar.zst": {"d2.cif"},
    }


def test_plan_archive_work_on_empty_group_is_empty():
    empty = _manifest().iloc[0:0]
    assert plan_archive_work(empty) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'stage1_extract_nested'`

- [ ] **Step 3: Write the implementation**

Create `scripts/stage1_extract_nested.py`:

```python
#!/usr/bin/env python3
"""Stage 1 (nested): extract backbones and epitope fingerprints.

Rewritten for the 2026-09-11 repack. The old stage 1 opened many small
`<uuid>.tar.zst` archives; those no longer exist as files. Each midrun is now one
solid `midrun_vN.tar.zst` with no random access, so the access pattern inverts:
ONE array task per outer archive, streaming it once and pulling everything the
manifest needs from it.

Usage (one array task):
    python scripts/stage1_extract_nested.py --shard 7 --n-shards 35 \
        --build build_shared --results results_shared
"""
from __future__ import annotations

import argparse
import io
import os
import shutil
import tarfile
import time
from pathlib import Path

import pandas as pd

from binder_map.backbone import coord_hash, write_backbone_pdb
from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.fingerprint import FINGERPRINT_COLUMNS, fingerprint
from binder_map.geometry import approach
from binder_map.nested import iter_members
from binder_map.rnaframe import FrameError

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
FM = RB / "analyses" / "foldseek_map"


def shard_outers(manifest: pd.DataFrame, shard: int, n_shards: int) -> pd.DataFrame:
    """Assign whole outer archives to shards, so each is streamed exactly once."""
    outers = sorted(manifest["outer_path"].unique())
    mine = set(outers[shard::n_shards])
    return manifest[manifest["outer_path"].isin(mine)]


def plan_archive_work(group: pd.DataFrame) -> dict[str, set[str]]:
    """Map inner archive name -> the member names needed from it."""
    wanted: dict[str, set[str]] = {}
    for inner, member in zip(group["inner_path"], group["member_path"]):
        wanted.setdefault(inner, set()).add(member)
    return wanted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, default=35)
    parser.add_argument("--build", type=Path, default=FM / "build_shared")
    parser.add_argument("--results", type=Path, default=FM / "results_shared")
    args = parser.parse_args()

    manifest = pd.read_csv(args.build / "manifests" / "design.tsv", sep="\t")
    mine = shard_outers(manifest, args.shard, args.n_shards)
    print(f"shard {args.shard}/{args.n_shards}: {len(mine):,d} structures from "
          f"{mine['outer_path'].nunique()} outer archives", flush=True)

    out_dir = args.results / "design"
    for sub in ("shards", "hashes", "fp", "missing"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"bm_s1n_{args.shard}"
    tmp.mkdir(parents=True, exist_ok=True)
    scratch = tmp / "one.cif"

    shard_name = f"shard_{args.shard:02d}"
    tar_path = out_dir / "shards" / f"{shard_name}.tar.gz"

    by_member = {(r.inner_path, r.member_path): r.design_id
                 for r in mine.itertuples()}
    hashes, fps, missing = [], [], []
    seen_local: set[str] = set()
    delivered: set[tuple[str, str]] = set()
    t0 = time.time()

    with tarfile.open(tar_path, "w:gz") as tar:
        for outer_path, group in mine.groupby("outer_path", sort=True):
            outer = RB / outer_path
            if not outer.exists():
                for design_id in group["design_id"]:
                    missing.append({"design_id": design_id, "reason": "archive_missing"})
                continue
            wanted = plan_archive_work(group)
            print(f"  streaming {outer_path} "
                  f"({len(wanted)} inner archives, {len(group):,d} members)", flush=True)

            for inner_name, member_name, blob in iter_members(outer, wanted):
                key = (inner_name, member_name)
                design_id = by_member.get(key)
                if design_id is None or not blob.strip():
                    continue
                delivered.add(key)
                try:
                    scratch.write_bytes(blob)
                    atoms = load_structure(scratch)
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
                record.update({f"per_nt_{i+1}": int(v) for i, v in enumerate(fp.per_nt)})
                record.update({
                    "contacts_offframe": fp.offframe, "contacts_total": fp.total,
                    "cx": float(app.centroid[0]), "cy": float(app.centroid[1]),
                    "cz": float(app.centroid[2]),
                    "ax": float(app.axis[0]), "ay": float(app.axis[1]),
                    "az": float(app.axis[2]),
                })
                fps.append(record)

                if digest in seen_local:
                    continue
                seen_local.add(digest)
                pdb = write_backbone_pdb(protein, design_id)
                info = tarfile.TarInfo(name=f"{design_id}.pdb")
                info.size = len(pdb)
                info.mtime = 0
                tar.addfile(info, io.BytesIO(pdb))

    for key, design_id in by_member.items():
        if key not in delivered:
            missing.append({"design_id": design_id, "reason": "member_missing"})

    pd.DataFrame(hashes, columns=["design_id", "coord_hash", "protein_chain", "n_res"]).to_csv(
        out_dir / "hashes" / f"{shard_name}.tsv", sep="\t", index=False)
    pd.DataFrame(fps).to_parquet(out_dir / "fp" / f"{shard_name}.parquet", index=False)
    pd.DataFrame(missing, columns=["design_id", "reason"]).to_csv(
        out_dir / "missing" / f"{shard_name}.tsv", sep="\t", index=False)

    shutil.rmtree(tmp, ignore_errors=True)

    # Fail loudly rather than leave a plausible-looking empty shard. A blank
    # output with the right row count is the failure mode that defeated the
    # resume check in July 2026, so assert on CONTENT, not on file existence.
    assert len(hashes) + len(missing) == len(mine), (
        f"accounting mismatch: {len(hashes)} parsed + {len(missing)} missing "
        f"!= {len(mine)} manifest rows")
    if len(mine) and not hashes:
        raise SystemExit(f"FATAL: shard {args.shard} parsed 0 of {len(mine):,d} structures")
    if tar_path.stat().st_size < 200:
        raise SystemExit(f"FATAL: shard tarball {tar_path} is empty")

    print(f"done in {time.time() - t0:.0f}s: {len(hashes):,d} parsed, "
          f"{len(seen_local):,d} unique backbones, {len(missing):,d} missing", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py -v`
Expected: PASS, 4 passed

- [ ] **Step 5: Write the sbatch wrapper**

Create `slurm/stage1_extract_nested.sbatch`:

```bash
#!/bin/bash
# Stage 1 (nested): one array task per outer midrun archive.
#   sbatch --array=0-34%35 slurm/stage1_extract_nested.sbatch
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --tmp=20000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

cd "$FM" || exit 1
echo "=== stage1-nested shard=$SLURM_ARRAY_TASK_ID on $(hostname) $(date) ==="
command -v zstd >/dev/null || { echo "FATAL: zstd not on PATH"; exit 1; }

"$PY" scripts/stage1_extract_nested.py \
  --shard "$SLURM_ARRAY_TASK_ID" \
  --n-shards "${N_SHARDS:-35}" \
  --build "$FM/${BUILD:-build_shared}" \
  --results "$FM/${RESULTS:-results_shared}"
rc=$?
echo "=== exit $rc $(date) ==="
exit $rc
```

- [ ] **Step 6: Smoke-test on one real archive**

Build a one-archive manifest and run the extractor outside Slurm:

```bash
cd $FM && PYTHONPATH=src:scripts $PY - <<'EOF'
import pandas as pd
from pathlib import Path
FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
m = pd.read_csv(FM / "build_shared/manifests/design.tsv", sep="\t")
one = m[m.outer_path == sorted(m.outer_path.unique())[0]].head(200)
(FM / "build_smoke_nested/manifests").mkdir(parents=True, exist_ok=True)
one.to_csv(FM / "build_smoke_nested/manifests/design.tsv", sep="\t", index=False)
print("smoke manifest:", len(one), "rows from", one.outer_path.iloc[0])
EOF

cd $FM && PYTHONPATH=src $PY scripts/stage1_extract_nested.py \
  --shard 0 --n-shards 1 --build $FM/build_smoke_nested --results $FM/results_smoke_nested
```

Expected: `done in ...s: 200 parsed, N unique backbones, 0 missing` with N ≤ 200,
and no FATAL. Then confirm the RFD3 chain inversion is still detected:

```bash
cd $FM && $PY -c "
import pandas as pd
h = pd.read_csv('results_smoke_nested/design/hashes/shard_00.tsv', sep='\t')
print(h['protein_chain'].value_counts())
print('n_res unique:', h['n_res'].unique())
assert (h['n_res'] == 80).all()
print('PASS')"
```

Expected: `n_res unique: [80]` and `PASS`. For an RFD3 archive the protein chain
must come out as `B`; if it reports `A`, stop — something upstream changed.

- [ ] **Step 7: Commit**

```bash
cd ~/Documents/binder_map
for f in scripts/stage1_extract_nested.py slurm/stage1_extract_nested.sbatch tests/test_stage1_nested.py; do
  done
git add scripts/stage1_extract_nested.py slurm/stage1_extract_nested.sbatch tests/test_stage1_nested.py
git commit -m "Rewrite stage 1 for one streaming pass per outer archive

Each midrun is now a solid tar.zst with no random access, so the task unit
becomes the archive rather than the design. Asserts on parsed content, not file
existence, because a blank shard with the right row count is the failure mode a
resume check cannot catch."
```

---

