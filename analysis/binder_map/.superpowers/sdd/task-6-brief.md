### Task 6: Stage 4 scale fixes

Two allocation sites break above ~100k nodes. They are fixed together because
stage 4 calls both and a reviewer would judge them as one change.

**Files:**
- Modify: `scripts/stage4_graph.py:71-95` (the hit loading block)
- Modify: `src/binder_map/knn.py:167` (the `block` default)
- Test: `tests/test_hits_stream.py`

**Interfaces:**
- Consumes: nothing new
- Produces: `load_hits(hit_files: list[Path], index: pandas.Series) -> tuple[ndarray, ndarray, ndarray, ndarray, int]`
  returning `(qi, tj, qtm, ttm, n_unknown)` with `qi`/`tj` int32 and `qtm`/`ttm` float32.

- [ ] **Step 1: Write the failing test**

Create `tests/test_hits_stream.py`:

```python
import numpy as np
import pandas as pd

from stage4_graph import load_hits

COLUMNS = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
           "lddt", "fident", "alnlen", "evalue", "bits"]


def _write_hits(path, rows):
    frame = pd.DataFrame(rows, columns=COLUMNS)
    frame.to_csv(path, sep="\t", header=False, index=False)


def test_load_hits_maps_names_to_node_ids(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
    ])
    qi, tj, qtm, ttm, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert list(qi) == [0] and list(tj) == [1]
    assert np.isclose(qtm[0], 0.9) and np.isclose(ttm[0], 0.8)
    assert unknown == 0


def test_load_hits_strips_the_pdb_suffix(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
    ])
    qi, _, _, _, _ = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 1


def test_load_hits_drops_and_counts_unknown_names(tmp_path):
    index = pd.Series([0], index=["s0"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s0.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100],
        ["s0.pdb", "ghost.pdb", 0.7, 0.6, 0.65, 0.6, 0.4, 80, 1e-7, 90],
    ])
    qi, _, _, _, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 1
    assert unknown == 1


def test_load_hits_concatenates_multiple_chunk_files(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100]])
    _write_hits(tmp_path / "hits_01.tsv", [
        ["s1.pdb", "s0.pdb", 0.6, 0.5, 0.55, 0.6, 0.4, 80, 1e-7, 90]])
    qi, tj, _, _, _ = load_hits(
        sorted(tmp_path.glob("hits_*.tsv")), index)
    assert sorted(zip(qi.tolist(), tj.tolist())) == [(0, 1), (1, 0)]


def test_load_hits_returns_compact_dtypes(tmp_path):
    index = pd.Series([0, 1], index=["s0", "s1"])
    _write_hits(tmp_path / "hits_00.tsv", [
        ["s0.pdb", "s1.pdb", 0.9, 0.8, 0.85, 0.7, 0.5, 80, 1e-9, 100]])
    qi, tj, qtm, ttm, _ = load_hits([tmp_path / "hits_00.tsv"], index)
    assert qi.dtype == np.int32 and tj.dtype == np.int32
    assert qtm.dtype == np.float32 and ttm.dtype == np.float32


def test_load_hits_on_empty_file_returns_empty_arrays(tmp_path):
    index = pd.Series([0], index=["s0"])
    (tmp_path / "hits_00.tsv").write_text("")
    qi, tj, qtm, ttm, unknown = load_hits([tmp_path / "hits_00.tsv"], index)
    assert len(qi) == 0 and len(tj) == 0 and unknown == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_hits_stream.py -v`
Expected: FAIL with `ImportError: cannot import name 'load_hits' from 'stage4_graph'`

- [ ] **Step 3: Add `load_hits` to `scripts/stage4_graph.py`**

Insert this function above `def main()`:

```python
HIT_COLUMNS = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
               "lddt", "fident", "alnlen", "evalue", "bits"]
HIT_CHUNK = 5_000_000


def load_hits(hit_files, index):
    """Stream hit TSVs into compact (qi, tj, qtm, ttm) arrays.

    The previous implementation did `pd.concat([pd.read_csv(f) for f in files])`,
    which materialises every row with two Python-string columns. At 31,763 nodes
    that was 1.3 GB of TSV and fine; at ~508k nodes it is ~20 GB and ~150M rows,
    which exhausts the node before the first edge is built.

    Mapping names to node ids per chunk and keeping int32/float32 costs ~12 bytes
    per directed hit, so the same data sits in ~1.8 GB.
    """
    qs, ts, qtms, ttms = [], [], [], []
    unknown = 0
    for path in hit_files:
        if Path(path).stat().st_size == 0:
            continue
        reader = pd.read_csv(
            path, sep="\t", names=HIT_COLUMNS, header=None,
            usecols=["query", "target", "qtmscore", "ttmscore"],
            dtype={"qtmscore": "float32", "ttmscore": "float32"},
            chunksize=HIT_CHUNK,
        )
        for chunk in reader:
            qi = index.reindex(chunk["query"].str.removesuffix(".pdb")).to_numpy()
            tj = index.reindex(chunk["target"].str.removesuffix(".pdb")).to_numpy()
            known = ~(np.isnan(qi) | np.isnan(tj))
            unknown += int((~known).sum())
            qs.append(qi[known].astype(np.int32))
            ts.append(tj[known].astype(np.int32))
            qtms.append(chunk["qtmscore"].to_numpy()[known])
            ttms.append(chunk["ttmscore"].to_numpy()[known])
    if not qs:
        empty_i = np.empty(0, np.int32)
        empty_f = np.empty(0, np.float32)
        return empty_i, empty_i.copy(), empty_f, empty_f.copy(), unknown
    return (np.concatenate(qs), np.concatenate(ts),
            np.concatenate(qtms), np.concatenate(ttms), unknown)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_hits_stream.py -v`
Expected: PASS, 6 passed

- [ ] **Step 5: Replace the concat block in `main()`**

In `scripts/stage4_graph.py`, delete these lines:

```python
    hits = pd.concat(
        [pd.read_csv(f, sep="\t", names=columns, header=None) for f in hit_files],
        ignore_index=True,
    )
    print(f"  {len(hits):,d} raw hits from {len(hit_files)} chunks")

    # Foldseek reports names with the .pdb suffix the tar members carried.
    for column in ("query", "target"):
        hits[column] = hits[column].str.removesuffix(".pdb")

    qi = index.reindex(hits["query"]).to_numpy()
    tj = index.reindex(hits["target"]).to_numpy()
    known = ~(np.isnan(qi) | np.isnan(tj))
    if (~known).any():
        print(f"  WARNING: {int((~known).sum()):,d} hits reference unknown names")
    qi, tj = qi[known].astype(int), tj[known].astype(int)
    qtm = hits["qtmscore"].to_numpy()[known]
    ttm = hits["ttmscore"].to_numpy()[known]
```

and replace them with:

```python
    qi, tj, qtm, ttm, unknown = load_hits(hit_files, index)
    print(f"  {len(qi):,d} usable hits from {len(hit_files)} chunks")
    if unknown:
        print(f"  WARNING: {unknown:,d} hits reference unknown names")
```

Also delete the now-unused local `columns = [...]` list immediately above the
hit-file glob, and add `from pathlib import Path` to the imports if absent.

- [ ] **Step 6: Change the kNN block default**

In `src/binder_map/knn.py:167`, change the signature:

```python
def knn_from_vectors(vectors, k, floor, block: int = 1024,
```

and replace the docstring's second paragraph with:

```python
    """Exact cosine kNN over L1-normalised rows, by blocked matrix product.

    `block` rows are scored against all n at once, so the temporary is
    block x n float32. At n = 508k a block of 4096 allocates 8.3 GB, which is
    why the default is 1024 (2 GB). The product stays exact; only the blocking
    changes.
    """
```

- [ ] **Step 7: Verify the existing suite still passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/ -v`
Expected: PASS, no regressions. `tests/test_knn.py` and `tests/test_clustering.py`
must be green — they pin the graph semantics this plan must not change.

- [ ] **Step 8: Commit**

```bash
cd ~/Documents/binder_map
for f in scripts/stage4_graph.py src/binder_map/knn.py tests/test_hits_stream.py; do
  done
git add scripts/stage4_graph.py src/binder_map/knn.py tests/test_hits_stream.py
git commit -m "Stream hit tables and shrink the kNN block for 500k nodes

pd.concat of every hit TSV holds ~150M rows with string columns (~20 GB at this
scale); streaming to int32/float32 per chunk costs ~1.8 GB. The epitope kNN's
4096-row block allocates 8.3 GB at n=508k, so the default drops to 1024."
```

---

