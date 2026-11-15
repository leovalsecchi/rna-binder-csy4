### Task 11: Stage 0 — manifests and explorer metadata

Runs locally. The census is a gate: if any RNA variant will not align to the canonical frame, the pipeline stops here, before any compute is spent.

**Files:**
- Create: `scripts/stage0_select.py`
- Create: `tests/test_stage0.py`

**Interfaces:**
- Consumes: `binder_map.RNA_LEN`.
- Produces, under `--out` (default `build/`):
  - `manifests/<layer>.tsv` — columns `row_id, design_id, uuid8, stack, campaign, gen, invfold, archive_path, member_path`, sorted by `archive_path` so stage 1 opens each archive once.
  - `meta/designs.parquet` — curated metadata keyed by `design_id`.
  - `meta/counts.json` — per-layer row counts and the input CSV sha256.
  - Importable helpers `read_passing(path) -> pandas.DataFrame` and `build_manifest(df, layer) -> pandas.DataFrame`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_stage0.py`:

```python
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from stage0_select import METADATA_COLUMNS, build_manifest, read_passing  # noqa: E402

CSV = Path("/Users/leo/Documents/binder_filter/exports/passing.csv")

pytestmark = pytest.mark.skipif(not CSV.exists(), reason="passing.csv not present")


@pytest.fixture(scope="module")
def df():
    return read_passing(CSV)


def test_reads_all_rows_past_the_two_comment_lines(df):
    assert len(df) == 71866
    assert "design_id" in df.columns
    assert df["design_id"].is_unique


def test_expected_stacks_present(df):
    assert set(df["stack"]) == {
        "RFD3+LiMPNN", "BG+LiMPNN", "BG+BIF", "BD1+LiMPNN", "HD+LiMPNN",
    }


def test_manifest_splits_archive_and_member(df):
    m = build_manifest(df, "esm")
    assert len(m) == len(df)
    assert not m["archive_path"].str.contains("!").any()
    assert m["archive_path"].str.endswith(".tar.zst").all()
    assert m["member_path"].str.endswith(".cif").all()


def test_manifest_is_sorted_by_archive(df):
    """Stage 1 batches by archive; unsorted input would reopen archives."""
    m = build_manifest(df, "af3")
    assert m["archive_path"].is_monotonic_increasing


def test_metadata_columns_all_exist(df):
    for column in METADATA_COLUMNS:
        assert column in df.columns, column


def test_hd_designs_have_no_seq_suffix(df):
    """HD design_ids are recycleN_sampleM_packK, not _seqN.

    Dedup must therefore be by coordinate hash, never by stripping _seqN.
    """
    hd = df[df["stack"] == "HD+LiMPNN"]["design_id"]
    assert len(hd) > 0
    assert not hd.str.contains(r"_seq\d+$", regex=True).any()
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_stage0.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'stage0_select'`

- [ ] **Step 3: Implement `scripts/stage0_select.py`**

```python
#!/usr/bin/env python3
"""Stage 0: turn passing.csv into per-layer manifests plus explorer metadata.

Also emits the RNA-variant census, which is a gate: if any variant fails to
align to the canonical 16-mer, the pipeline stops here rather than silently
misregistering every epitope downstream.

Runs locally. Usage:
    python scripts/stage0_select.py --csv <passing.csv> --out build/
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

LAYERS = {"design": "design_cif", "af3": "af3_cif", "esm": "esm_cif"}

MANIFEST_COLUMNS = [
    "row_id", "design_id", "uuid8", "stack", "campaign", "gen", "invfold",
    "archive_path", "member_path",
]

# Curated rather than all 65 columns, to keep the HTML payload bounded.
METADATA_COLUMNS = [
    "design_id", "uuid8", "stack", "campaign", "gen", "invfold", "length",
    "dockq_af3", "dockq_esm", "iptm_af3", "iptm_esm", "iptm_mean", "iptm_min",
    "rmsd_af3", "rmsd_esm", "dsasa_af3", "dsasa_esm", "hb_af3", "hb_esm",
    "rog_af3", "rog_esm", "foot_af3", "foot_esm",
    "pose_disagree", "fold_disagree",
]


def read_passing(path: Path) -> pd.DataFrame:
    """Read passing.csv, skipping its two leading comment lines."""
    return pd.read_csv(path, skiprows=2, low_memory=False)


def build_manifest(df: pd.DataFrame, layer: str) -> pd.DataFrame:
    column = LAYERS[layer]
    split = df[column].str.split("!", n=1, expand=True)
    out = pd.DataFrame({
        "row_id": range(len(df)),
        "design_id": df["design_id"].values,
        "uuid8": df["uuid8"].values,
        "stack": df["stack"].values,
        "campaign": df["campaign"].values,
        "gen": df["gen"].values,
        "invfold": df["invfold"].values,
        "archive_path": split[0].values,
        "member_path": split[1].values,
    })
    return out.sort_values("archive_path", kind="stable").reset_index(drop=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path,
                        default=Path("/Users/leo/Documents/binder_filter/exports/passing.csv"))
    parser.add_argument("--out", type=Path, default=Path("build"))
    parser.add_argument("--subsample", type=int, default=0,
                        help="stratified designs per stack, for the smoke test")
    args = parser.parse_args()

    df = read_passing(args.csv)
    print(f"read {len(df):,d} rows from {args.csv}")

    if args.subsample:
        df = (df.groupby("stack", group_keys=False)
                .apply(lambda g: g.sample(min(len(g), args.subsample),
                                          random_state=42))
                .reset_index(drop=True))
        print(f"stratified subsample: {len(df):,d} rows "
              f"({args.subsample} per stack)")

    (args.out / "manifests").mkdir(parents=True, exist_ok=True)
    (args.out / "meta").mkdir(parents=True, exist_ok=True)

    counts = {}
    for layer in LAYERS:
        manifest = build_manifest(df, layer)
        dest = args.out / "manifests" / f"{layer}.tsv"
        manifest.to_csv(dest, sep="\t", index=False)
        counts[layer] = len(manifest)
        print(f"  {layer:7s} {len(manifest):>7,d} rows, "
              f"{manifest['archive_path'].nunique():>6,d} archives -> {dest}")

    df[METADATA_COLUMNS].to_parquet(args.out / "meta" / "designs.parquet",
                                    index=False)

    meta = {
        "input_csv": str(args.csv),
        "input_sha256": sha256(args.csv),
        "n_rows": int(len(df)),
        "subsample_per_stack": args.subsample or None,
        "layer_counts": counts,
        "stack_counts": df["stack"].value_counts().to_dict(),
    }
    (args.out / "meta" / "counts.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote metadata for {len(df):,d} designs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_stage0.py -v
```

Expected: 6 passed.

- [ ] **Step 5: Run stage 0 for real**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python scripts/stage0_select.py --out build
cat build/meta/counts.json
```

Expected: 71,866 rows; each layer 71,866; roughly 6,780 archives for design and similar for the refolds.

- [ ] **Step 6: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage0_select.py tests/test_stage0.py
git commit -m "Add stage 0: per-layer manifests and curated explorer metadata"
```

---

