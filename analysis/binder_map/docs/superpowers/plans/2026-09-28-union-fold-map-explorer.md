# Union Fold Map Explorer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Co-embed the campaign sample, the passing pool and the ordered 12k library into one fold map of 526,467 backbones, and render it as a single self-contained explorer with toggleable layers and conserved backbones marked.

**Architecture:** Assemble a new `results_union` tree from shards already on disk (no archive extraction), run the existing stages 2-6 over it with frozen parameters on the fold axis, then build a canvas-rendered page whose per-node payload is packed binary rather than JSON.

**Tech Stack:** Foldseek v8-ef4e960, Slurm, pandas/pyarrow, leidenalg, umap-learn, vanilla JS + canvas. Python 3.12 via the Euler venv.

**Spec:** `docs/superpowers/specs/2026-09-28-union-fold-map-explorer-design.md`

## Global Constraints

- **Author locally, test on Euler, commit locally.** The laptop has NO pytest and no pyarrow. Loop: edit in `~/Documents/binder_map`, `bash sync.sh`, run tests over ssh, commit locally. Never `git add -A` — the repo carries large untracked result trees.
- **Euler paths:** `FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map`, `PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python`.
- **Run tests as:** `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/<file> -v"`. Note `pyproject` sets `addopts = "-q"`; pass `-o addopts=` when you need counts.
- **Frozen parameters, unchanged everywhere:** `--alignment-type 1`, `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a`, `TM = min(qTM,tTM)`, `--tm-floor 0.5`, `--epi-floor 0.6`, `--k 30`, `--mutual`, Leiden fixed seed, operative resolution `cluster_r80`.
- **Do not modify stages 2, 3, 4, 5 or 6**, or `dist/binder_map_v6.html`. Stage 3's target-DB staging to node-local scratch (added 2026-09-25) must stay.
- **Measured node counts — assert, never assume:** union 526,467; `sample500k` 499,999; `passing60k` 31,763; `library12k` 9,759; conserved sample∩passing 5,295, sample∩library 1,553, passing∩library 9,759, all three 1,553.
- **Never write loose per-structure files** under `results/` — hard inode cap. **Never expand a nested archive to scratch.**
- **No Boltz-2 / `b2`** anywhere.
- **Never run `sbatch` from a test.** Tests must not submit cluster jobs.

---

### Task 1: Union assembly (stage U0)

Assemble `results_union/design` from two source trees and reconcile the hashes against the structures actually present. This is the task the spec's "Known hazard" section is about: `results_shared` holds 508,207 hash rows but fewer structures, and `results_v6` holds 55,706 rows against 31,763 nodes. Passing unfiltered hashes to stage 2 reproduces the failure that killed the 2026-09-24 run.

**Files:**
- Create: `scripts/stageU0_union.py`
- Test: `tests/test_union.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `present_stems(shard_dir: Path) -> set[str]` — every `.pdb` member stem across `shard_*.tar.gz`.
  - `filter_hashes(frame: pd.DataFrame, present: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]` — `(kept, dropped)`, split on `design_id in present`.
  - `layer_flags(nodes_by_layer: dict[str, set[str]], hash_of: dict[str, str]) -> pd.DataFrame` — columns `coord_hash, in_sample500k, in_passing60k, in_12k`.
  - CLI writing `results_union/design/{shards,hashes,fp,missing}/` plus `meta/layers.parquet`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_union.py
import tarfile
import pandas as pd
import pytest
from pathlib import Path
from stageU0_union import present_stems, filter_hashes, layer_flags


def _shard(path: Path, stems: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as tar:
        for s in stems:
            body = b"ATOM\n"
            info = tarfile.TarInfo(f"{s}.pdb")
            info.size = len(body)
            import io
            tar.addfile(info, io.BytesIO(body))


def test_present_stems_reads_every_shard(tmp_path):
    _shard(tmp_path / "shards" / "shard_00.tar.gz", ["a", "b"])
    _shard(tmp_path / "shards" / "shard_01.tar.gz", ["c"])
    assert present_stems(tmp_path / "shards") == {"a", "b", "c"}


def test_present_stems_ignores_non_pdb(tmp_path):
    p = tmp_path / "shards" / "shard_00.tar.gz"
    p.parent.mkdir(parents=True)
    import io
    with tarfile.open(p, "w:gz") as tar:
        for name in ("a.pdb", "notes.txt"):
            body = b"x"
            info = tarfile.TarInfo(name)
            info.size = len(body)
            tar.addfile(info, io.BytesIO(body))
    assert present_stems(p.parent) == {"a"}


def test_filter_hashes_drops_rows_without_a_structure():
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h1"]})
    kept, dropped = filter_hashes(frame, {"b"})
    assert list(kept["design_id"]) == ["b"]
    assert list(dropped["design_id"]) == ["a"]


def test_filter_hashes_keeps_everything_when_all_present():
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h2"]})
    kept, dropped = filter_hashes(frame, {"a", "b"})
    assert len(kept) == 2 and len(dropped) == 0


def test_filter_hashes_never_drops_a_whole_hash_group():
    # every coord_hash that had a row must still have one afterwards
    frame = pd.DataFrame({"design_id": ["a", "b"], "coord_hash": ["h1", "h1"]})
    kept, _ = filter_hashes(frame, {"a"})
    assert set(kept["coord_hash"]) == {"h1"}


def test_layer_flags_marks_conserved_nodes():
    hash_of = {"d1": "h1", "d2": "h2", "d3": "h3"}
    layers = {"sample500k": {"d1", "d2"}, "passing60k": {"d2"}, "library12k": {"d3"}}
    out = layer_flags(layers, hash_of).set_index("coord_hash")
    assert bool(out.loc["h2", "in_sample500k"]) and bool(out.loc["h2", "in_passing60k"])
    assert not bool(out.loc["h1", "in_passing60k"])
    assert bool(out.loc["h3", "in_12k"])
    assert len(out) == 3


def test_layer_flags_every_node_carries_at_least_one_flag():
    hash_of = {"d1": "h1"}
    out = layer_flags({"sample500k": {"d1"}, "passing60k": set(), "library12k": set()}, hash_of)
    cols = ["in_sample500k", "in_passing60k", "in_12k"]
    assert bool(out[cols].any(axis=1).all())
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union.py -o addopts= -v"`
Expected: FAIL, `ModuleNotFoundError: No module named 'stageU0_union'`.

- [ ] **Step 3: Implement `scripts/stageU0_union.py`**

Write the three pure functions above plus a `main()` that:
1. Takes `--source` twice (repeatable: `results_shared`, `results_v6`), `--out results_union`, and `--layer-nodes NAME:PATH` three times pointing at each source map's `dedup/nodes.tsv`.
2. Copies (or hard-links) every `shard_*.tar.gz` from each source into `results_union/design/shards/`, prefixing the filename with the source name so the two trees cannot collide (`shared_shard_00.tar.gz`, `v6_shard_00.tar.gz`).
3. Computes `present = present_stems(out/"shards")` **once, after all shards are in place**.
4. Concatenates each source's `hashes/shard_*.tsv`, runs `filter_hashes`, writes the kept rows to `results_union/design/hashes/` and the dropped rows to `missing/` with a `reason` column of `duplicate_structure`.
5. Copies `fp/shard_*.parquet` from both sources (needed for the epitope axis).
6. Builds `meta/layers.parquet` via `layer_flags`, mapping each layer's `struct_name` through the concatenated hash table.
7. Prints the counts and **asserts** them against the Global Constraints table; a mismatch exits non-zero.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union.py -o addopts= -v"`
Expected: PASS, 7 tests.

- [ ] **Step 5: Run it for real and check the asserts**

Run:
```
ssh euler "cd $FM && PYTHONPATH=src $PY scripts/stageU0_union.py \
  --source results_shared --source results_v6 --out results_union \
  --layer-nodes sample500k:results_full/design/dedup/nodes.tsv \
  --layer-nodes passing60k:results_sw_k30_f50_mut/design/dedup/nodes.tsv \
  --layer-nodes library12k:results_lib12k/design/dedup/nodes.tsv"
```
Expected: union 526,467; conserved 5,295 / 1,553 / 9,759 / 1,553; exit 0.

- [ ] **Step 6: Commit**

```bash
git add scripts/stageU0_union.py tests/test_union.py
git commit -m "Assemble the union tree and reconcile its hashes against the shards"
```

---

### Task 2: Union build driver (stages 2-5, fold axis)

**Files:**
- Create: `submit_union.sh`
- Test: `tests/test_union_submit.py`

**Interfaces:**
- Consumes: `results_union/design/{shards,hashes,fp}` from Task 1.
- Produces: job ids on stdout; `results_union/design/clusters/fold/clusters.parquet`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_union_submit.py
import os, subprocess, textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run_with_stub(tmp_path):
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "sbatch").write_text("#!/bin/bash\necho 12345\n")
    (stub / "sbatch").chmod(0o755)
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    return subprocess.run(["bash", str(ROOT / "submit_union.sh")],
                          capture_output=True, text=True, env=env, cwd=tmp_path)


def test_submits_fold_axis_only(tmp_path):
    out = _run_with_stub(tmp_path).stdout
    assert "fold" in out
    assert "epitope" not in out


def test_frozen_parameters_present():
    text = (ROOT / "submit_union.sh").read_text()
    for token in ("--tm-floor 0.5", "--epi-floor 0.6", "--k 30", "--mutual", "ALN_TYPE=1"):
        assert token in text, token
```

- [ ] **Step 2: Run it to verify it fails**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union_submit.py -o addopts= -v"`
Expected: FAIL — `submit_union.sh` does not exist.

- [ ] **Step 3: Write `submit_union.sh`**

Model it on `submit_maps.sh`, with these differences:
- No stage 1 or 1b — Task 1 already produced the tree.
- `set -euo pipefail` **and** `shopt -s inherit_errexit`, and `|| exit 1` on every `sbatch`. Any `$(...)` capture of a function needs its own outer `|| exit 1`; the inner guards do not propagate.
- Stage 2: `RESULTS=results_union bash scripts/stage2_builddb.sh design`, then `stage3_split_queries.py --n-chunks 512 --results $FM/results_union`.
- Stage 3: `--array=0-511%128 --time=00:20:00 --mem-per-cpu=1G --export=ALL,ALN_TYPE=1,RESULTS=results_union slurm/stage3_search.sbatch design`. These are the 2026-09-25 measured values (56 s, 927 MB peak), not guesses.
- Stage 4: `stage4_graph.py --layer design --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual --results $FM/results_union`.
- Stage 5 once: `--axis fold`, depending `afterok` on stage 4. The epitope axis is deliberately out of scope (2026-09-28); the fingerprints stay in `fp/` so it can be added later without redoing the search.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union_submit.py -o addopts= -v"`
Expected: PASS, 2 tests.

- [ ] **Step 5: `bash -n` and submit**

Run: `bash -n submit_union.sh && ssh euler "cd $FM && bash submit_union.sh"`
Expected: job ids for s2, s3, s4 and s5-fold, each depending on the previous.

- [ ] **Step 6: Commit**

```bash
git add submit_union.sh tests/test_union_submit.py
git commit -m "Add the union build driver for both axes"
```

---

### Task 3: Fold embedding at 526k, sized from measurement

`stage6_embed.py` is not modified. This task establishes what it costs at this scale and submits it correctly. The kNN block already had to shrink once for memory here (4096 → 1024, an 8.3 GB allocation), so it is sized from a run, not an estimate.

**Files:**
- Create: `submit_union_embed.sh`
- Modify: none.

**Interfaces:**
- Consumes: `results_union/design/clusters/fold/clusters.parquet`.
- Produces: `results_union/design/embed/fold/embed.parquet`.

- [ ] **Step 1: Measure on the fold axis with a generous ceiling**

Submit one job: `--cpus-per-task=16 --mem-per-cpu=16G --time=08:00:00`, running
`stage6_embed.py --layer design --axis fold --k 30 --results $FM/results_union`.

- [ ] **Step 2: Read the real cost**

Run: `ssh euler "sacct -j <id> --format=JobID%18,State%12,Elapsed,MaxRSS"`
Record Elapsed and MaxRSS. If it OOMs or times out, **do not simply raise the ceiling** — report the numbers and stop; the spec names this as the failure point and the response is a decision, not a retry.

- [ ] **Step 3: Record the measured shape**

Write `submit_union_embed.sh` carrying the measured `--mem-per-cpu` and `--time` (roughly 2× measured MaxRSS per cpu and 3× measured Elapsed), each stated in a comment naming the measurement that produced it, so the epitope axis can be run later without re-measuring.

- [ ] **Step 4: Verify the embedding has a row per node**

```
ssh euler "cd $FM && $PY -c \"
import pandas as pd
e = pd.read_parquet('results_union/design/embed/fold/embed.parquet')
print('fold', len(e))
assert len(e) == 526467, len(e)
print('PASS')\""
```
Expected: `fold 526467`, `PASS`.

- [ ] **Step 5: Commit**

```bash
git add submit_union_embed.sh
git commit -m "Size and submit the union fold embedding from a measured run"
```

---

### Task 4: The explorer page

**Files:**
- Create: `scripts/stageU8_build_union_html.py`
- Create: `src/binder_map/unionhtml.py`
- Create: `web_union/template.html`
- Test: `tests/test_unionhtml.py`

**Interfaces:**
- Consumes: `results_union/design/embed/fold/embed.parquet`, `clusters/fold/clusters.parquet`, `dedup/nodes.tsv`, `meta/layers.parquet`.
- Produces: `dist/binder_map_union.html`.
  - `pack_coords(xs, ys) -> str` — base64 of a little-endian `float32` interleaved `[x0,y0,x1,y1,...]`.
  - `pack_uint8(values) -> str` — base64 of a `uint8` array.
  - `layer_mask(flags: pd.DataFrame) -> np.ndarray` — one `uint8` per node, bit 0 sample, bit 1 passing, bit 2 library.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_unionhtml.py
import base64
import numpy as np
import pandas as pd
from binder_map.unionhtml import pack_coords, pack_uint8, layer_mask


def test_pack_coords_roundtrips():
    xs, ys = [1.5, -2.25], [0.0, 3.75]
    raw = np.frombuffer(base64.b64decode(pack_coords(xs, ys)), dtype="<f4")
    assert raw.tolist() == [1.5, 0.0, -2.25, 3.75]


def test_pack_coords_is_float32_not_float64():
    # 4 bytes per value: a float64 payload would be twice the size
    assert len(base64.b64decode(pack_coords([1.0], [2.0]))) == 8


def test_pack_uint8_roundtrips():
    raw = np.frombuffer(base64.b64decode(pack_uint8([0, 7, 255])), dtype=np.uint8)
    assert raw.tolist() == [0, 7, 255]


def test_layer_mask_sets_one_bit_per_layer():
    flags = pd.DataFrame({"in_sample500k": [True, False, True],
                          "in_passing60k": [False, True, True],
                          "in_12k": [False, False, True]})
    assert layer_mask(flags).tolist() == [1, 2, 7]


def test_layer_mask_conserved_nodes_have_multiple_bits():
    flags = pd.DataFrame({"in_sample500k": [True], "in_passing60k": [True], "in_12k": [False]})
    mask = layer_mask(flags)
    assert bin(int(mask[0])).count("1") == 2
```

- [ ] **Step 2: Run them to verify they fail**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_unionhtml.py -o addopts= -v"`
Expected: FAIL, `ModuleNotFoundError: No module named 'binder_map.unionhtml'`.

- [ ] **Step 3: Implement `src/binder_map/unionhtml.py`**

```python
"""Pack the union map's per-node arrays for the browser.

526,467 nodes x 2 axes cannot travel as JSON objects: one dict per node would
put the page in the hundreds of MB. Everything the page needs per node is a
fixed-width array, so it ships as base64 typed arrays and is decoded once.
"""
import base64

import numpy as np
import pandas as pd


def pack_coords(xs, ys) -> str:
    arr = np.empty(2 * len(xs), dtype="<f4")
    arr[0::2] = np.asarray(xs, dtype="<f4")
    arr[1::2] = np.asarray(ys, dtype="<f4")
    return base64.b64encode(arr.tobytes()).decode("ascii")


def pack_uint8(values) -> str:
    return base64.b64encode(np.asarray(values, dtype=np.uint8).tobytes()).decode("ascii")


#: bit 0 = campaign sample, bit 1 = passing pool, bit 2 = ordered library.
#: A node with more than one bit set is conserved -- that is the whole signal.
def layer_mask(flags: pd.DataFrame) -> np.ndarray:
    return (flags["in_sample500k"].to_numpy(bool).astype(np.uint8)
            | (flags["in_passing60k"].to_numpy(bool).astype(np.uint8) << 1)
            | (flags["in_12k"].to_numpy(bool).astype(np.uint8) << 2))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_unionhtml.py -o addopts= -v"`
Expected: PASS, 5 tests.

- [ ] **Step 5: Write `web_union/template.html` and the builder**

The template is a single page with a `<canvas>` and no external assets. Required behaviour:
- Decode the base64 arrays once into `Float32Array` / `Uint8Array`.
- **Three layer checkboxes.** A node draws when `mask & activeMask` is non-zero.
- **Conserved** nodes — `popcount(mask & activeMask) > 1` — drawn last, in a distinct colour, above the single-layer points, with a live count in the legend.
- Canvas rendering with pan/zoom. Redraw on interaction must stay interactive at 526k points: draw into an `ImageData` buffer rather than issuing one `arc()` per node.
- Hover: nearest node within a few pixels via a uniform grid built once; tooltip shows design id, stack, union cluster id, and which layers it is in.
- Click: `POST` to `http://localhost:45000` with the design id, wrapped in try/catch so a missing server changes nothing on screen.
- A visible line stating the union clustering is its own partition and is not the `compare/` result.
- Colours must be defined as CSS custom properties on `:root`, with a dark-mode block; the page must render correctly in both.

`scripts/stageU8_build_union_html.py` reads the parquets, joins on `node_id`, packs the arrays, and substitutes them into the template. It asserts 526,467 nodes before writing and prints the finished file size.

- [ ] **Step 6: Build it and check the size target**

Run: `ssh euler "cd $FM && PYTHONPATH=src $PY scripts/stageU8_build_union_html.py --results results_union --out dist/binder_map_union.html"`
Expected: `526467 nodes`, and a file at or under 25 MB (the fold-only payload should land well under it). If it exceeds the target, decimate the campaign-only layer for display and print the decimation factor on the page — do not silently drop nodes.

- [ ] **Step 7: Open it and look at it**

Copy to the laptop and open in a browser. Confirm: the map renders, all three toggles work, conserved points are visible and counted, hover shows a tooltip, and the page does not stall while panning. A chart that was never opened is not finished.

- [ ] **Step 8: Commit**

```bash
git add scripts/stageU8_build_union_html.py src/binder_map/unionhtml.py web_union/template.html tests/test_unionhtml.py
git commit -m "Build the union explorer with toggleable layers and conserved backbones"
```

---

## Notes for whoever executes this

- **The union clustering is a fourth clustering.** Never quote its cluster counts as a result of the three-map comparison in `compare/`.
- **`$( )` clears `set -e` and swallows `exit`.** This project has been bitten twice — once in `submit_maps.sh`, once in the gate's probe arm. Use `shopt -s inherit_errexit` *and* explicit `|| exit 1`.
- **Stage 3 stages the target DB to node-local scratch.** Do not "simplify" it back to reading from the project filesystem: 265 concurrent readers collapsed the Lustre client and no task finished in 105 minutes.
- **Counting files is not measuring progress.** Stage 1 overwrites in place, so a stale tree looks complete. Use mtime.
