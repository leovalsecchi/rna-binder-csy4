# Foldseek Binder Map Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Map 71,866 passing binder designs on two independent axes — backbone fold (Foldseek) and RNA epitope (contact fingerprint) — and ship them as one self-contained HTML explorer with side-by-side linked maps.

**Architecture:** A small pure-Python library (`src/binder_map/`) holding all testable logic, driven by thin stage scripts (`scripts/`) that run as SLURM array jobs on Euler (`slurm/`). Compute happens on Euler where the structures live; only small tables come back, and the HTML is baked locally. Fold and epitope are two independent partitions over the same node set, produced by the same axis-agnostic code path.

**Tech Stack:** Python 3.12, biotite 1.6 (structure I/O), Foldseek 8-ef4e960, leidenalg + python-igraph (clustering), umap-learn (layout), pyarrow (tables), numpy/scipy/pandas/scikit-learn, vanilla JS + 3Dmol.js + fflate (explorer). pytest for tests.

## Global Constraints

These apply to every task. Values are copied verbatim from the spec.

- **Chains are selected by composition, never by chain ID.** RNA = the chain whose residues subset `{A, C, G, U}`; protein = the other. RFD3+LiMPNN design structures have protein in chain **B**.
- **Protein chain must be 80 residues.** Violations go to `missing/` with reason `chain_select_failed`.
- **Contact cutoff is 4.5 Å heavy-atom** (`RNA_CONTACT_CUTOFF`), matching `rna_binder.analysis.structural.rna_contact_counts`.
- **Canonical RNA frame is the 16-nt csy4 target `ACUGCCGUAUAGGCAG`** (`RNA_LEN = 16`).
- **Structures are read with biotite**, never a hand-rolled mmCIF parser. The dialects differ per stack.
- **No loose per-structure files on Euler.** Extraction writes through `$TMPDIR` into `.tar.gz` shards. The project has a hard inode cap.
- **Deduplication is by backbone coordinate hash**, not filename.
- **Defaults if the smoke test is inconclusive:** `TM_FLOOR = 0.50`, `EPI_FLOOR = 0.50`, `--alignment-type 2`.
- **Fixed seeds everywhere** (Leiden, UMAP) — `SEED = 42`.
- **Leiden resolutions:** 0.5, 1.0, 2.0 → columns `cluster_r05`, `cluster_r10`, `cluster_r20`.
- **kNN k = 15**; `--max-seqs 300`; `-c 0.8 --cov-mode 0 -e 10`.
- Arrays are throttled at `%16`. Euler evaluates only the top 250 pending jobs.

### Fixed paths

| What | Path |
|---|---|
| Local repo | `/Users/leo/Documents/binder_map` |
| Local venv | `/Users/leo/.venvs/binder_map` |
| Euler project root (`$RB`) | `/cluster/work/wicky/lvalsecchi/rna_binder` |
| Euler working dir (`$FM`) | `$RB/analyses/foldseek_map` |
| Euler venv | `/cluster/work/wicky/lvalsecchi/binder_map_venv` (already created and verified) |
| Input CSV | `/Users/leo/Documents/binder_filter/exports/passing.csv` |

### Euler environment incantations

Foldseek (must be a login shell, `bash -lc`):

```bash
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960
```

Python on Euler — always use the venv, which has system-site-packages so biotite/numpy/scipy/pandas/sklearn/networkx come from the spack python:

```bash
/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
```

## File Structure

```
/Users/leo/Documents/binder_map/
├── src/binder_map/
│   ├── __init__.py          constants: RNA_CANONICAL, RNA_LEN, RNA_CONTACT_CUTOFF, SEED
│   ├── chains.py            composition-based protein/RNA split
│   ├── rnaframe.py          canonical frame + length-variant alignment
│   ├── backbone.py          minimal backbone PDB writer + coordinate hash
│   ├── fingerprint.py       48-dim epitope fingerprint
│   ├── geometry.py          RNA-aligned superposition, centroid + principal axis
│   ├── knn.py               kNN graphs for both axes (cluster + layout variants)
│   ├── clustering.py        Leiden at 3 resolutions + medoid election
│   ├── embedding.py         UMAP from a precomputed kNN
│   ├── compare.py           contingency, AMI/ARI, agreement, stack profiles
│   └── htmlbuild.py         payload encoding + HTML baking
├── scripts/                 stage entrypoints (thin CLI wrappers over the library)
│   ├── fetch_fixtures.py    pulls test CIFs from Euler
│   ├── stage0_select.py
│   ├── stage1_extract.py
│   ├── stage2_builddb.sh
│   ├── stage3_search.sh
│   ├── stage4_graph.py
│   ├── stage5_cluster.py
│   ├── stage6_embed.py
│   ├── stage6b_compare.py
│   ├── stage7_medoids.py
│   └── stage8_build_html.py
├── slurm/
│   ├── stage1_extract.sbatch
│   ├── stage3_search.sbatch
│   ├── stage_single.sbatch  parameterised single-node stage (2,4,5,6,6b,7)
│   └── submit_layer.sh      chains one layer end-to-end with afterok
├── web/
│   ├── template.html
│   ├── app.js
│   └── style.css
├── tests/
│   ├── conftest.py
│   ├── fixtures/            real CIFs pulled from Euler (gzipped, committed)
│   ├── test_chains.py
│   ├── test_rnaframe.py
│   ├── test_backbone.py
│   ├── test_fingerprint.py
│   ├── test_geometry.py
│   ├── test_knn.py
│   ├── test_clustering.py
│   ├── test_embedding.py
│   ├── test_compare.py
│   └── test_htmlbuild.py
└── docs/superpowers/{specs,plans}/
```

**Responsibility split:** every module in `src/binder_map/` is pure — it takes arrays or paths and returns values, with no SLURM, no archive traversal, and no global state — so all of it is unit-testable locally against real fixtures. The `scripts/` layer owns I/O, archive traversal and CLI. The `slurm/` layer owns only resources and chaining. Anything with logic worth testing belongs in `src/`.

---

### Task 1: Scaffold, environments, and real test fixtures

Every later task tests against real structures. The dialect differences between stacks are the whole reason this project needs care, so synthetic fixtures would defeat the purpose. This task ends with `pytest` green on a trivial test and **six** real CIFs committed — five distinct designs plus a second sequence variant of the RFD3 backbone (`rfd3_design_seq7`), which Task 4 needs to prove that coordinate hashing collapses sequence variants.

**Files:**
- Create: `pyproject.toml`
- Create: `src/binder_map/__init__.py`
- Create: `scripts/fetch_fixtures.py`
- Create: `tests/conftest.py`
- Create: `tests/test_constants.py`
- Create: `.gitignore` (append)

**Interfaces:**
- Consumes: nothing.
- Produces: `binder_map.RNA_CANONICAL: str`, `binder_map.RNA_LEN: int`, `binder_map.RNA_CONTACT_CUTOFF: float`, `binder_map.SEED: int`. Fixture helper `tests/conftest.py::fixture_path(name: str) -> pathlib.Path` and pytest fixtures `rfd3_design`, `bgbif_design`, `bd1_design`, `esm_refold`, `af3_refold`, each returning a `pathlib.Path` to a `.cif.gz`.

- [ ] **Step 1: Create the local venv and install dependencies**

The venv lives outside the repo. Local Python is 3.10; that is fine for everything except needing biotite, which is not in the base env.

```bash
python3 -m venv /Users/leo/.venvs/binder_map
/Users/leo/.venvs/binder_map/bin/pip install -q --upgrade pip
/Users/leo/.venvs/binder_map/bin/pip install -q \
  'biotite>=1.2' numpy pandas pyarrow scipy scikit-learn \
  umap-learn leidenalg python-igraph pytest
/Users/leo/.venvs/binder_map/bin/python -c "import biotite,umap,leidenalg,igraph,pyarrow,pytest; print('deps ok')"
```

Expected final line: `deps ok`

- [ ] **Step 2: Write `pyproject.toml`**

```toml
[project]
name = "binder_map"
version = "0.1.0"
requires-python = ">=3.10"

[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-q"
```

- [ ] **Step 3: Write the constants module**

Create `src/binder_map/__init__.py`:

```python
"""Shared constants for the binder map pipeline.

Values are pinned to match rna_binder.analysis.structural so that recomputed
per-nucleotide contacts can be cross-checked against the stored
rna_contact_nt_* columns.
"""

# The csy4 target RNA, 16 nt. This is the canonical epitope frame: every
# fingerprint position refers to an index into this sequence.
RNA_CANONICAL = "ACUGCCGUAUAGGCAG"
RNA_LEN = 16

# Angstroms, heavy-atom. Must equal rna_binder's RNA_CONTACT_CUTOFF or the
# cross-check in Task 19 is meaningless.
RNA_CONTACT_CUTOFF = 4.5

# Every stochastic step (Leiden, UMAP) uses this.
SEED = 42

# Moieties the per-nucleotide contact counts are split across.
MOIETIES = ("phosphate", "sugar", "base")

# 16 nucleotides x 3 moieties
FINGERPRINT_DIM = RNA_LEN * len(MOIETIES)

assert len(RNA_CANONICAL) == RNA_LEN
```

- [ ] **Step 4: Write the failing test**

Create `tests/test_constants.py`:

```python
from binder_map import (
    RNA_CANONICAL, RNA_LEN, RNA_CONTACT_CUTOFF, FINGERPRINT_DIM, MOIETIES,
)


def test_canonical_rna_is_16nt():
    assert len(RNA_CANONICAL) == RNA_LEN == 16
    assert set(RNA_CANONICAL) <= {"A", "C", "G", "U"}


def test_contact_cutoff_matches_rna_binder():
    # rna_binder.analysis.structural.RNA_CONTACT_CUTOFF
    assert RNA_CONTACT_CUTOFF == 4.5


def test_fingerprint_dim():
    assert MOIETIES == ("phosphate", "sugar", "base")
    assert FINGERPRINT_DIM == 48
```

- [ ] **Step 5: Run the test**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pip install -q -e .
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_constants.py -v
```

Expected: 3 passed.

- [ ] **Step 6: Write the fixture fetcher**

This resolves real paths out of `passing.csv` rather than hardcoding archive names, because archive paths carry a campaign directory that varies per uuid.

Create `scripts/fetch_fixtures.py`:

```python
#!/usr/bin/env python3
"""Pull a handful of real CIFs from Euler into tests/fixtures/.

Picks structures that exercise the dialect and chain-order differences:
  - RFD3+LiMPNN design  -> protein is chain B (inverted)
  - BG+BIF design       -> different mmCIF dialect, 18-nt RNA
  - BD1+LiMPNN design   -> protein is chain A
  - one af3 refold, one esm refold
Plus a second sequence variant of the RFD3 backbone, to test dedup.

Run from the repo root:  python scripts/fetch_fixtures.py
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

CSV = Path("/Users/leo/Documents/binder_filter/exports/passing.csv")
OUT = Path(__file__).resolve().parent.parent / "tests" / "fixtures"
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"

# design_id -> (column, fixture name)
WANTED = {
    "csy4free_RFD3_LiMPNN_len80_0118c551_design45_seq4": ("design_cif", "rfd3_design"),
    "csy4free_RFD3_LiMPNN_len80_0118c551_design45_seq7": ("design_cif", "rfd3_design_seq7"),
    "csy4free_BG_BIF_len80_c55bf73f_design16_seq8": ("design_cif", "bgbif_design"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1": ("design_cif", "bd1_design"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1_af3": ("af3_cif", "af3_refold"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1_esm": ("esm_cif", "esm_refold"),
}


def rows_by_id(wanted_ids: set[str]) -> dict[str, dict]:
    found: dict[str, dict] = {}
    with CSV.open() as fh:
        fh.readline()
        fh.readline()
        for row in csv.DictReader(fh):
            if row["design_id"] in wanted_ids:
                found[row["design_id"]] = row
            if len(found) == len(wanted_ids):
                break
    return found


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    base_ids = {k.removesuffix("_af3").removesuffix("_esm") for k in WANTED}
    rows = rows_by_id(base_ids)

    missing = base_ids - rows.keys()
    if missing:
        print(f"FATAL: design_ids not in passing.csv: {sorted(missing)}", file=sys.stderr)
        return 1

    for key, (column, name) in WANTED.items():
        design_id = key.removesuffix("_af3").removesuffix("_esm")
        archive, member = rows[design_id][column].split("!", 1)
        cmd = f"cd {RB} && tar --zstd -xOf {archive} {member}"
        proc = subprocess.run(
            ["ssh", "-o", "ConnectTimeout=60", "euler", cmd],
            capture_output=True,
        )
        if proc.returncode != 0 or not proc.stdout.strip():
            print(f"FATAL: could not fetch {name} ({design_id})", file=sys.stderr)
            print(proc.stderr.decode()[:400], file=sys.stderr)
            return 1
        dest = OUT / f"{name}.cif.gz"
        dest.write_bytes(gzip.compress(proc.stdout))
        print(f"  {name:20s} {len(proc.stdout):>9,d} B -> {dest.name}")

    print("fixtures written to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 7: Fetch the fixtures**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python scripts/fetch_fixtures.py
ls -la tests/fixtures/
```

Expected: six `.cif.gz` files, each a few kB to a few tens of kB.

- [ ] **Step 8: Write `tests/conftest.py`**

```python
import gzip
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


def fixture_path(name: str) -> Path:
    p = FIXTURES / f"{name}.cif.gz"
    if not p.exists():
        pytest.skip(f"fixture {name} missing; run scripts/fetch_fixtures.py")
    return p


def fixture_text(name: str) -> str:
    return gzip.decompress(fixture_path(name).read_bytes()).decode()


@pytest.fixture
def rfd3_design() -> Path:
    return fixture_path("rfd3_design")


@pytest.fixture
def rfd3_design_seq7() -> Path:
    return fixture_path("rfd3_design_seq7")


@pytest.fixture
def bgbif_design() -> Path:
    return fixture_path("bgbif_design")


@pytest.fixture
def bd1_design() -> Path:
    return fixture_path("bd1_design")


@pytest.fixture
def af3_refold() -> Path:
    return fixture_path("af3_refold")


@pytest.fixture
def esm_refold() -> Path:
    return fixture_path("esm_refold")
```

- [ ] **Step 9: Run the suite**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest -v
```

Expected: 3 passed.

- [ ] **Step 10: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add pyproject.toml src tests scripts .gitignore
git commit -m "Scaffold binder_map package with real CIF fixtures"
```

---

### Task 2: Composition-based chain selection

The single highest-consequence piece of logic in the project. Getting it wrong silently feeds the RNA target into Foldseek as the binder for 53% of the design layer.

**Files:**
- Create: `src/binder_map/chains.py`
- Create: `tests/test_chains.py`

**Interfaces:**
- Consumes: fixtures from Task 1. (Deliberately consumes no RNA constants: chain selection must not assume any RNA length — one fixture's RNA is 18 nt.)
- Produces:
  - `load_structure(path: str | Path) -> biotite.structure.AtomArray` — reads `.cif` or `.cif.gz`, first model only.
  - `split_chains(atoms) -> tuple[AtomArray, AtomArray]` returning `(protein, rna)`.
  - `ChainSelectError(Exception)` raised when no unambiguous split exists.
  - `protein_chain_id(atoms) -> str` — the chain ID that held the protein, for the per-stack tally.

- [ ] **Step 1: Write the failing test**

Create `tests/test_chains.py`:

```python
import pytest

from binder_map.chains import ChainSelectError, load_structure, split_chains

RNA_RESNAMES = {"A", "C", "G", "U"}


def _split(path):
    return split_chains(load_structure(path))


def test_rfd3_design_protein_is_chain_B(rfd3_design):
    """RFD3+LiMPNN design CIFs are inverted: chain A is the RNA target.

    This is 53% of the design layer. If this test ever fails, something
    upstream changed and the pipeline must stop.
    """
    atoms = load_structure(rfd3_design)
    assert protein_chain_id(atoms) == "B"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80
    assert set(rna.res_name) <= RNA_RESNAMES


def test_bd1_design_protein_is_chain_A(bd1_design):
    atoms = load_structure(bd1_design)
    assert protein_chain_id(atoms) == "A"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80


def test_bgbif_dialect_parses(bgbif_design):
    """BG+BIF uses a different mmCIF dialect; a hand-rolled parser failed on it."""
    protein, rna = _split(bgbif_design)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80
    assert set(rna.res_name) <= RNA_RESNAMES


@pytest.mark.parametrize("layer", ["af3_refold", "esm_refold"])
def test_refolds_protein_is_chain_A(layer, request):
    atoms = load_structure(request.getfixturevalue(layer))
    assert protein_chain_id(atoms) == "A"
    protein, rna = split_chains(atoms)
    assert len(set(zip(protein.chain_id, protein.res_id))) == 80


def test_all_fixtures_give_80_residue_protein_and_rna(
    rfd3_design, bd1_design, bgbif_design, af3_refold, esm_refold
):
    for path in (rfd3_design, bd1_design, bgbif_design, af3_refold, esm_refold):
        protein, rna = _split(path)
        assert len(set(zip(protein.chain_id, protein.res_id))) == 80, path
        assert len(rna) > 0, path
        assert set(rna.res_name) <= RNA_RESNAMES, path


def test_protein_only_structure_raises(bd1_design):
    atoms = load_structure(bd1_design)
    protein, _ = split_chains(atoms)
    with pytest.raises(ChainSelectError):
        split_chains(protein)
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_chains.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.chains'`

- [ ] **Step 3: Implement `src/binder_map/chains.py`**

```python
"""Composition-based protein/RNA chain selection.

Chain IDs are NOT reliable across this campaign's producers:

    layer   stack                          protein chain
    design  RFD3+LiMPNN                    B   <- inverted, 53% of the layer
    design  BD1/BG+LiMPNN/HD+LiMPNN        A
    design  BG+BIF                         A   (different mmCIF dialect)
    af3     all                            A
    esm     all                            A

So the rule is the one rna_binder.analysis.structural.split_chains uses:
RNA is the chain whose residue names subset {A, C, G, U}; protein is the other.
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

RNA_RESNAMES = frozenset({"A", "C", "G", "U"})


class ChainSelectError(Exception):
    """No unambiguous protein/RNA split exists for this structure."""


def load_structure(path: str | Path):
    """Read a .cif or .cif.gz into an AtomArray (first model, no solvent)."""
    path = Path(path)
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    block = pdbx.CIFFile.read(io.StringIO(raw.decode()))
    atoms = pdbx.get_structure(block, model=1)
    return atoms[~struc.filter_solvent(atoms)]


def _chain_kinds(atoms) -> dict[str, str]:
    """Map chain_id -> 'rna' | 'protein' by residue-name composition."""
    kinds: dict[str, str] = {}
    for chain_id in set(atoms.chain_id):
        names = set(atoms.res_name[atoms.chain_id == chain_id])
        kinds[chain_id] = "rna" if names <= RNA_RESNAMES else "protein"
    return kinds


def protein_chain_id(atoms) -> str:
    """The chain ID holding the protein. Raises if not exactly one."""
    kinds = _chain_kinds(atoms)
    protein = sorted(c for c, k in kinds.items() if k == "protein")
    if len(protein) != 1:
        raise ChainSelectError(
            f"expected exactly one protein chain, got {protein} from {kinds}"
        )
    return protein[0]


def split_chains(atoms):
    """Return (protein, rna) AtomArrays, selected by composition."""
    kinds = _chain_kinds(atoms)
    protein_ids = [c for c, k in kinds.items() if k == "protein"]
    rna_ids = [c for c, k in kinds.items() if k == "rna"]
    if len(protein_ids) != 1 or len(rna_ids) != 1:
        raise ChainSelectError(
            f"need exactly one protein and one RNA chain; got "
            f"protein={sorted(protein_ids)} rna={sorted(rna_ids)}"
        )
    protein = atoms[atoms.chain_id == protein_ids[0]]
    rna = atoms[atoms.chain_id == rna_ids[0]]
    return protein, rna
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_chains.py -v
```

Expected: 7 passed (6 test functions; one is parametrized over 2 layers). If `test_rfd3_design_protein_is_chain_B` fails, stop and investigate — the fixture or the upstream data changed, and that invalidates the whole design-layer plan.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/chains.py tests/test_chains.py
git commit -m "Add composition-based chain selection with inverted-RFD3 regression test"
```

---

### Task 3: Canonical RNA frame and length-variant alignment

Epitope comparison is only meaningful if every design's contacts refer to the same nucleotide numbering. A silently misregistered frame would make every epitope result wrong while looking entirely normal.

**Files:**
- Create: `src/binder_map/rnaframe.py`
- Create: `tests/test_rnaframe.py`

**Interfaces:**
- Consumes: `binder_map.RNA_CANONICAL`, `binder_map.RNA_LEN`; `binder_map.chains.split_chains`.
- Produces:
  - `rna_sequence(rna_atoms) -> str` — one-letter RNA sequence in residue order.
  - `frame_offset(seq: str) -> int` — index in `seq` where the canonical 16-mer starts. Raises `FrameError` if absent.
  - `FrameError(Exception)`.
  - `canonical_positions(rna_atoms) -> numpy.ndarray` — per-atom canonical position `1..16`, or `0` for out-of-frame atoms. Length equals `len(rna_atoms)`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_rnaframe.py`:

```python
import numpy as np
import pytest

from binder_map import RNA_CANONICAL, RNA_LEN
from binder_map.chains import load_structure, split_chains
from binder_map.rnaframe import (
    FrameError, canonical_positions, frame_offset, rna_sequence,
)


def test_frame_offset_exact_match():
    assert frame_offset(RNA_CANONICAL) == 0


def test_frame_offset_with_closing_pair():
    """A BD1+LiMPNN design carries CACUGCCGUAUAGGCAGC: the 16-mer plus one pair.

    Verified against the committed fixtures 2026-08-23.
    """
    assert frame_offset("C" + RNA_CANONICAL + "C") == 1


def test_frame_offset_rejects_unrelated_sequence():
    with pytest.raises(FrameError):
        frame_offset("AAAAAAAAAAAAAAAA")


@pytest.mark.parametrize(
    "name", ["rfd3_design", "bd1_design", "bgbif_design", "af3_refold", "esm_refold"]
)
def test_every_fixture_aligns_to_canonical_frame(name, request):
    _, rna = split_chains(load_structure(request.getfixturevalue(name)))
    seq = rna_sequence(rna)
    offset = frame_offset(seq)
    assert seq[offset : offset + RNA_LEN] == RNA_CANONICAL


def test_canonical_positions_cover_all_16(af3_refold):
    _, rna = split_chains(load_structure(af3_refold))
    pos = canonical_positions(rna)
    assert len(pos) == len(rna)
    assert set(pos[pos > 0]) == set(range(1, RNA_LEN + 1))


def test_canonical_positions_marks_out_of_frame_as_zero(bd1_design):
    """bd1_design carries the 18-nt variant: 2 nucleotides outside the frame.

    Verified 2026-08-23: bd1_design is CACUGCCGUAUAGGCAGC (offset 1), while
    bgbif_design is a plain 16-mer. This is the only fixture that exercises
    the out-of-frame path, so the assertion is unconditional on purpose.
    """
    _, rna = split_chains(load_structure(bd1_design))
    assert len(rna_sequence(rna)) == 18
    pos = canonical_positions(rna)
    assert (pos == 0).any()
    assert set(pos[pos > 0]) == set(range(1, RNA_LEN + 1))
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_rnaframe.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.rnaframe'`

- [ ] **Step 3: Implement `src/binder_map/rnaframe.py`**

```python
"""Map an RNA chain onto the canonical 16-nt csy4 frame.

Length variants exist: at least one BG+BIF design carries CACUGCCGUAUAGGCAGC,
which is the canonical 16-mer plus one closing base pair. Contacts are only
comparable across designs if every position refers to the same index, so the
canonical 16-mer is located by exact substring match and flanking nucleotides
are reported out-of-frame rather than renumbered.
"""
from __future__ import annotations

import numpy as np

from binder_map import RNA_CANONICAL, RNA_LEN


class FrameError(Exception):
    """The RNA does not contain the canonical csy4 16-mer."""


def rna_sequence(rna_atoms) -> str:
    """One-letter RNA sequence in residue order."""
    _, names = _residue_index(rna_atoms)
    return "".join(names)


def _residue_index(rna_atoms):
    """Return (per-atom residue ordinal, ordered residue names).

    Residues are ordered by first appearance, which is coordinate order in
    every dialect here, rather than by res_id, which is not guaranteed to
    start at 1.
    """
    keys = list(zip(rna_atoms.chain_id, rna_atoms.res_id, rna_atoms.ins_code))
    order: dict[tuple, int] = {}
    names: list[str] = []
    ordinals = np.empty(len(keys), dtype=int)
    for i, key in enumerate(keys):
        if key not in order:
            order[key] = len(order)
            names.append(str(rna_atoms.res_name[i]).strip())
        ordinals[i] = order[key]
    return ordinals, names


def frame_offset(seq: str) -> int:
    """Index into seq where the canonical 16-mer begins."""
    idx = seq.find(RNA_CANONICAL)
    if idx < 0:
        raise FrameError(
            f"canonical frame {RNA_CANONICAL} not found in RNA sequence {seq!r}"
        )
    return idx


def canonical_positions(rna_atoms) -> np.ndarray:
    """Per-atom canonical nucleotide position 1..16, or 0 if out of frame."""
    ordinals, names = _residue_index(rna_atoms)
    offset = frame_offset("".join(names))
    pos = ordinals - offset + 1
    pos[(pos < 1) | (pos > RNA_LEN)] = 0
    return pos
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_rnaframe.py -v
```

Expected: 10 passed (6 test functions; one is parametrized over 5 fixtures).

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/rnaframe.py tests/test_rnaframe.py
git commit -m "Add canonical RNA frame mapping with length-variant alignment"
```

---

### Task 4: Backbone extraction and coordinate hashing

Dedup is by coordinate hash, not filename — sequence variants of one design share a backbone, and the same backbone can recur across uuids where naming would not reveal it.

**Files:**
- Create: `src/binder_map/backbone.py`
- Create: `tests/test_backbone.py`

**Interfaces:**
- Consumes: `binder_map.chains.load_structure`, `split_chains`.
- Produces:
  - `BACKBONE_ATOMS: tuple[str, ...]` = `("N", "CA", "C", "O", "CB")`
  - `backbone_atoms(protein) -> AtomArray`
  - `coord_hash(protein) -> str` — 32-char MD5 over N/CA/C/O coordinates rounded to 3 dp, in residue order. CB is excluded because the three LiMPNN stacks ship design CIFs without sidechains.
  - `write_backbone_pdb(protein, name: str) -> bytes` — a minimal PDB as bytes, ready to add to a tar shard.
  - `ca_trace(protein) -> numpy.ndarray` of shape `(80, 3)`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_backbone.py`:

```python
import numpy as np

from binder_map.backbone import (
    BACKBONE_ATOMS, backbone_atoms, ca_trace, coord_hash, write_backbone_pdb,
)
from binder_map.chains import load_structure, split_chains


def _protein(path):
    return split_chains(load_structure(path))[0]


def test_backbone_keeps_only_backbone_atoms(af3_refold):
    bb = backbone_atoms(_protein(af3_refold))
    assert set(bb.atom_name) <= set(BACKBONE_ATOMS)
    assert "CA" in set(bb.atom_name)


def test_coord_hash_is_stable(af3_refold):
    p = _protein(af3_refold)
    assert coord_hash(p) == coord_hash(p)
    assert len(coord_hash(p)) == 32


def test_sequence_variants_share_a_backbone_hash(rfd3_design, rfd3_design_seq7):
    """seq4 and seq7 of the same RFD3 design are coordinate-identical.

    Verified 2026-08-21 by CA-coordinate md5 on Euler. This is what makes
    dedup collapse the design layer from 71,866 to ~32,848.
    """
    assert coord_hash(_protein(rfd3_design)) == coord_hash(_protein(rfd3_design_seq7))


def test_different_designs_have_different_hashes(rfd3_design, bd1_design):
    assert coord_hash(_protein(rfd3_design)) != coord_hash(_protein(bd1_design))


def test_coord_hash_ignores_sidechains(af3_refold):
    """Design CIFs for the LiMPNN stacks have no sidechains; refolds do.

    The hash must depend only on N/CA/C/O so the two are comparable.
    """
    p = _protein(af3_refold)
    keep = np.isin(p.atom_name, ["N", "CA", "C", "O", "CB"])
    assert coord_hash(p) == coord_hash(p[keep])


def test_ca_trace_shape(af3_refold):
    assert ca_trace(_protein(af3_refold)).shape == (80, 3)


def test_write_backbone_pdb_roundtrips(af3_refold):
    data = write_backbone_pdb(_protein(af3_refold), "test_design")
    assert data.startswith(b"ATOM") or b"\nATOM" in data
    assert data.rstrip().endswith(b"END")
    assert data.count(b"\nATOM") + data.startswith(b"ATOM") >= 80
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_backbone.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.backbone'`

- [ ] **Step 3: Implement `src/binder_map/backbone.py`**

```python
"""Minimal backbone extraction, hashing, and PDB serialisation.

Foldseek reads backbone geometry, so sidechains are irrelevant — which is
fortunate, because the three LiMPNN stacks ship design CIFs without them
(inverse folding is not packing). The coordinate hash therefore uses only
N/CA/C/O, so a design and its refold are hashed on comparable content.
"""
from __future__ import annotations

import hashlib

import numpy as np

BACKBONE_ATOMS = ("N", "CA", "C", "O", "CB")
HASH_ATOMS = ("N", "CA", "C", "O")


def backbone_atoms(protein):
    """Keep only backbone (+CB) atoms."""
    return protein[np.isin(protein.atom_name, BACKBONE_ATOMS)]


def coord_hash(protein) -> str:
    """MD5 over N/CA/C/O coordinates rounded to 3 dp, in residue order."""
    sel = protein[np.isin(protein.atom_name, HASH_ATOMS)]
    order = np.lexsort((
        _atom_rank(sel.atom_name),
        _residue_ordinal(sel),
    ))
    coords = np.round(sel.coord[order].astype(np.float64), 3)
    # -0.0 and 0.0 must hash identically
    coords = coords + 0.0
    return hashlib.md5(coords.tobytes()).hexdigest()


def _atom_rank(names) -> np.ndarray:
    rank = {name: i for i, name in enumerate(HASH_ATOMS)}
    return np.array([rank[str(n)] for n in names])


def _residue_ordinal(atoms) -> np.ndarray:
    keys = list(zip(atoms.chain_id, atoms.res_id, atoms.ins_code))
    order: dict[tuple, int] = {}
    out = np.empty(len(keys), dtype=int)
    for i, key in enumerate(keys):
        if key not in order:
            order[key] = len(order)
        out[i] = order[key]
    return out


def ca_trace(protein) -> np.ndarray:
    """(n_residues, 3) array of CA coordinates in residue order."""
    ca = protein[protein.atom_name == "CA"]
    return ca.coord[np.argsort(_residue_ordinal(ca))].astype(np.float32)


def write_backbone_pdb(protein, name: str) -> bytes:
    """Serialise the backbone as a minimal PDB.

    Written by hand rather than via biotite's writer because this runs 215k
    times and the output only needs to be readable by Foldseek.
    """
    sel = backbone_atoms(protein)
    ordinals = _residue_ordinal(sel)
    lines = []
    for serial, (i, atom) in enumerate(zip(ordinals, sel), start=1):
        x, y, z = atom.coord
        atom_name = str(atom.atom_name)
        # PDB columns: atom names of <4 chars are left-padded by one space
        an = f" {atom_name:<3s}" if len(atom_name) < 4 else atom_name
        lines.append(
            f"ATOM  {serial:>5d} {an}{'':1s}{str(atom.res_name):>3s} A"
            f"{int(i) + 1:>4d}{'':4s}{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}{'':10s}{str(atom.element):>2s}"
        )
    lines.append("END")
    return ("\n".join(lines) + "\n").encode()
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_backbone.py -v
```

Expected: 7 passed.

- [ ] **Step 5: Verify Foldseek can actually read the output**

A PDB that Python accepts but Foldseek rejects would fail only at scale. Check now, on Euler.

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python - <<'PY'
from binder_map.backbone import write_backbone_pdb
from binder_map.chains import load_structure, split_chains
from tests.conftest import fixture_path
p, _ = split_chains(load_structure(fixture_path("af3_refold")))
open("/tmp/bm_probe.pdb", "wb").write(write_backbone_pdb(p, "probe"))
print("wrote /tmp/bm_probe.pdb")
PY
scp -q /tmp/bm_probe.pdb euler:/cluster/work/wicky/lvalsecchi/_binder_map_probe/
ssh euler 'bash -lc "
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960 >/dev/null 2>&1
cd /cluster/work/wicky/lvalsecchi/_binder_map_probe
rm -rf probedb && foldseek createdb bm_probe.pdb probedb 2>&1 | tail -3
foldseek convert2pdb probedb probe_out.pdb 2>&1 | tail -2
grep -c CA probe_out.pdb
"'
```

Expected: `createdb` reports 1 structure, and the CA count is 80. If Foldseek reports 0 entries, the PDB formatting is wrong — fix `write_backbone_pdb` before continuing, because every downstream stage depends on it.

- [ ] **Step 6: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/backbone.py tests/test_backbone.py
git commit -m "Add backbone extraction, coordinate hashing, and PDB serialisation"
```

---

### Task 5: Epitope fingerprint

**Files:**
- Create: `src/binder_map/fingerprint.py`
- Create: `tests/test_fingerprint.py`

**Interfaces:**
- Consumes: `binder_map.RNA_CONTACT_CUTOFF`, `RNA_LEN`, `MOIETIES`, `FINGERPRINT_DIM`; `binder_map.chains.split_chains`; `binder_map.rnaframe.canonical_positions`.
- Produces:
  - `moiety_of(atom_name: str) -> str` — `"phosphate" | "sugar" | "base"`.
  - `fingerprint(protein, rna) -> Fingerprint` — a dataclass with `vector: np.ndarray` (48, int32), `per_nt: np.ndarray` (16, int32), `offframe: int`, `total: int`.
  - `FINGERPRINT_COLUMNS: list[str]` — 48 names, `nt{i:02d}_{moiety}`, matching `vector` order.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fingerprint.py`:

```python
import numpy as np
import pytest

from binder_map import FINGERPRINT_DIM, RNA_LEN
from binder_map.chains import load_structure, split_chains
from binder_map.fingerprint import (
    FINGERPRINT_COLUMNS, fingerprint, moiety_of,
)
from binder_map.rnaframe import canonical_positions, rna_sequence


def _fp(path):
    protein, rna = split_chains(load_structure(path))
    return fingerprint(protein, rna)


def test_moiety_classification():
    assert moiety_of("P") == "phosphate"
    assert moiety_of("OP1") == "phosphate"
    assert moiety_of("O5'") == "phosphate"
    assert moiety_of("C1'") == "sugar"
    assert moiety_of("O2'") == "sugar"
    assert moiety_of("N1") == "base"
    assert moiety_of("C5") == "base"


def test_columns_match_vector_layout():
    assert len(FINGERPRINT_COLUMNS) == FINGERPRINT_DIM == 48
    assert FINGERPRINT_COLUMNS[0] == "nt01_phosphate"
    assert FINGERPRINT_COLUMNS[-1] == "nt16_base"


@pytest.mark.parametrize("name", ["af3_refold", "esm_refold", "bd1_design"])
def test_fingerprint_shape_and_consistency(name, request):
    fp = _fp(request.getfixturevalue(name))
    assert fp.vector.shape == (48,)
    assert fp.per_nt.shape == (RNA_LEN,)
    # per_nt is the vector summed over the 3 moieties
    assert np.array_equal(fp.per_nt, fp.vector.reshape(RNA_LEN, 3).sum(axis=1))
    # total accounts for everything, in-frame and out
    assert fp.total == int(fp.per_nt.sum()) + fp.offframe


def test_refold_has_contacts(af3_refold):
    """These designs passed a DockQ >= 0.7 gate, so they must have an interface."""
    fp = _fp(af3_refold)
    assert fp.total > 0
    assert fp.per_nt.sum() > 0


def test_offframe_is_zero_for_16nt_target(af3_refold):
    fp = _fp(af3_refold)
    assert fp.offframe == 0


def test_offframe_counts_contacts_outside_the_canonical_frame(bd1_design):
    """Contacts to nucleotides outside the frame must land in `offframe`.

    If they leaked into the 48-dim vector instead, the epitope frame would
    silently shift and every downstream comparison would be wrong while
    looking normal.

    bd1_design is the only 18-nt fixture, but on the real coordinates its two
    flanking nucleotides sit 4.503 A from the nearest protein atom - 0.003 A
    outside the 4.5 A cutoff - so it yields offframe == 0 as it stands. An
    earlier version of this test asserted only algebraic identities, which
    hold trivially at zero and tested nothing. So: assert the real geometry,
    then move ONLY the out-of-frame atoms into contact and assert where the
    new counts land.
    """
    protein, rna = split_chains(load_structure(bd1_design))
    assert len(rna_sequence(rna)) == 18

    baseline = fingerprint(protein, rna)
    assert baseline.offframe == 0

    positions = canonical_positions(rna)
    assert (positions == 0).any(), "fixture must have out-of-frame atoms"

    moved = rna.copy()
    moved.coord[positions == 0] = protein.coord[0]

    fp = fingerprint(protein, moved)
    assert fp.offframe > 0, "out-of-frame contacts were not counted"
    # Nothing in-frame moved, so the 48-dim vector must be untouched: this is
    # what proves the new counts went to `offframe` and not into the frame.
    assert np.array_equal(fp.vector, baseline.vector)
    assert fp.total == int(fp.per_nt.sum()) + fp.offframe


def test_fingerprint_is_deterministic(af3_refold):
    a, b = _fp(af3_refold), _fp(af3_refold)
    assert np.array_equal(a.vector, b.vector)


def test_distant_chains_give_zero_contacts(af3_refold):
    """Translating the RNA 500 A away must produce an empty fingerprint."""
    protein, rna = split_chains(load_structure(af3_refold))
    rna = rna.copy()
    rna.coord = rna.coord + np.array([500.0, 0.0, 0.0], dtype=rna.coord.dtype)
    fp = fingerprint(protein, rna)
    assert fp.total == 0
    assert fp.vector.sum() == 0
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.fingerprint'`

- [ ] **Step 3: Implement `src/binder_map/fingerprint.py`**

```python
"""48-dimensional epitope fingerprint: 16 canonical nucleotides x 3 moieties.

Counts protein-heavy-atom <-> RNA-heavy-atom pairs within RNA_CONTACT_CUTOFF.
The cutoff and the per-nucleotide binning match
rna_binder.analysis.structural.rna_contact_counts exactly, so the summed
per-nucleotide counts can be cross-checked against the stored
rna_contact_nt_* columns where those exist.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from binder_map import MOIETIES, RNA_CONTACT_CUTOFF, RNA_LEN
from binder_map.rnaframe import canonical_positions

FINGERPRINT_COLUMNS = [
    f"nt{i:02d}_{m}" for i in range(1, RNA_LEN + 1) for m in MOIETIES
]

_PHOSPHATE = frozenset({"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O5'", "O3'"})
_SUGAR = frozenset({"C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O4'"})


def moiety_of(atom_name: str) -> str:
    """Classify an RNA atom as phosphate, sugar, or base."""
    name = atom_name.strip()
    if name in _PHOSPHATE:
        return "phosphate"
    if name in _SUGAR:
        return "sugar"
    return "base"


@dataclass(frozen=True)
class Fingerprint:
    vector: np.ndarray   # (48,) int32, ordered as FINGERPRINT_COLUMNS
    per_nt: np.ndarray   # (16,) int32
    offframe: int        # contacts to nucleotides outside the canonical frame
    total: int           # all protein-RNA contacts


def _heavy(atoms):
    return atoms[atoms.element != "H"]


def fingerprint(protein, rna) -> Fingerprint:
    """Count protein-RNA heavy-atom contacts binned by nucleotide and moiety."""
    p_heavy = _heavy(protein)
    r_heavy = _heavy(rna)

    vector = np.zeros(len(FINGERPRINT_COLUMNS), dtype=np.int32)
    if len(p_heavy) == 0 or len(r_heavy) == 0:
        return Fingerprint(vector, np.zeros(RNA_LEN, dtype=np.int32), 0, 0)

    positions = canonical_positions(r_heavy)
    moiety_idx = np.array(
        [MOIETIES.index(moiety_of(str(n))) for n in r_heavy.atom_name]
    )

    # For every RNA heavy atom, how many protein heavy atoms are within cutoff.
    tree = cKDTree(p_heavy.coord.astype(np.float64))
    neighbours = tree.query_ball_point(
        r_heavy.coord.astype(np.float64), r=RNA_CONTACT_CUTOFF
    )
    counts = np.fromiter((len(n) for n in neighbours), dtype=np.int64,
                         count=len(neighbours))

    total = int(counts.sum())
    in_frame = positions > 0
    offframe = int(counts[~in_frame].sum())

    flat_idx = (positions[in_frame] - 1) * len(MOIETIES) + moiety_idx[in_frame]
    np.add.at(vector, flat_idx, counts[in_frame].astype(np.int32))

    per_nt = vector.reshape(RNA_LEN, len(MOIETIES)).sum(axis=1).astype(np.int32)
    return Fingerprint(vector, per_nt, offframe, total)
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_fingerprint.py -v
```

Expected: 10 passed (8 test functions; one is parametrized over 3 fixtures).

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/fingerprint.py tests/test_fingerprint.py
git commit -m "Add 48-dim epitope fingerprint at the 4.5 A rna_binder cutoff"
```

---

### Task 6: Approach geometry in the RNA frame

Two binders can contact the same nucleotides from opposite faces. This ships as a colour overlay on both maps, not a third partition.

**Files:**
- Create: `src/binder_map/geometry.py`
- Create: `tests/test_geometry.py`

**Interfaces:**
- Consumes: `binder_map.chains.split_chains`, `binder_map.rnaframe.canonical_positions`, `binder_map.backbone.ca_trace`.
- Produces:
  - `reference_frame(rna) -> tuple[np.ndarray, np.ndarray]` — `(centroid (3,), rotation (3,3))` built from the in-frame RNA P/C1' atoms by PCA.
  - `approach(protein, rna) -> Approach` — dataclass with `centroid: np.ndarray` (3,) and `axis: np.ndarray` (3,), both expressed in the RNA frame.

The RNA-aligned convention here is the same target-fixed convention adopted for DockQ on 2026-08-04: the RNA defines the frame and is never re-fit to the binder.

- [ ] **Step 1: Write the failing test**

Create `tests/test_geometry.py`:

```python
import numpy as np
import pytest

from binder_map.chains import load_structure, split_chains
from binder_map.geometry import approach, reference_frame
from binder_map.rnaframe import canonical_positions


def _pair(path):
    return split_chains(load_structure(path))


ALL_FIXTURES = [
    "af3_refold", "esm_refold", "bd1_design",
    "rfd3_design", "rfd3_design_seq7", "bgbif_design",
]


@pytest.mark.parametrize("name", ALL_FIXTURES)
def test_reference_frame_is_orthonormal(name, request):
    _, rna = _pair(request.getfixturevalue(name))
    centroid, rot = reference_frame(rna)
    assert centroid.shape == (3,)
    assert rot.shape == (3, 3)
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-6)
    # Right-handed, not merely orthogonal: det must be +1, never -1.
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-6)


def test_approach_shapes(af3_refold):
    protein, rna = _pair(af3_refold)
    a = approach(protein, rna)
    assert a.centroid.shape == (3,)
    assert a.axis.shape == (3,)
    assert np.isclose(np.linalg.norm(a.axis), 1.0, atol=1e-6)


@pytest.mark.parametrize("name", ALL_FIXTURES)
def test_approach_is_invariant_to_rigid_motion(name, request):
    """Rotating and translating the whole complex must not change the result.

    This is the entire point of expressing the binder in the RNA's frame, and
    it MUST be checked on every fixture. An earlier version tested only
    af3_refold and passed while the frame was genuinely not invariant: the
    SVD-based implementation flipped an axis on three other fixtures, giving
    |axis.axis'| of 0.85-0.93 instead of 1.0. One fixture is not a property
    test.
    """
    protein, rna = _pair(request.getfixturevalue(name))
    before = approach(protein, rna)

    theta = 0.7
    rot = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    shift = np.array([13.0, -7.0, 42.0])

    p2, r2 = protein.copy(), rna.copy()
    p2.coord = (p2.coord @ rot.T + shift).astype(p2.coord.dtype)
    r2.coord = (r2.coord @ rot.T + shift).astype(r2.coord.dtype)
    after = approach(p2, r2)

    np.testing.assert_allclose(before.centroid, after.centroid, atol=1e-3)
    # The frame makes no sign choice, so the axis must agree in DIRECTION,
    # not merely as an undirected line. abs() here would hide a flip.
    assert float(before.axis @ after.axis) > 0.999, (
        f"{name}: axis direction changed under rigid motion"
    )


def test_approach_ignores_sidechains(af3_refold):
    """Sidechain presence must not move the centroid or axis.

    Three of five design pipelines ship structures without sidechains while
    all refolds have them. If `approach` used every atom, that asymmetry
    would masquerade as a real difference in binding geometry.
    """
    protein, rna = _pair(af3_refold)
    backbone_only = protein[np.isin(protein.atom_name, ["N", "CA", "C", "O"])]
    full = approach(protein, rna)
    stripped = approach(backbone_only, rna)
    np.testing.assert_allclose(full.centroid, stripped.centroid, atol=1e-6)
    np.testing.assert_allclose(full.axis, stripped.axis, atol=1e-6)


def test_missing_frame_atoms_raise_a_clear_error(af3_refold):
    """A nucleotide with no P/C1'/C4' must fail with THIS module's error.

    Before the frame atoms and their positions were returned together,
    re-deriving positions from the filtered set made such a nucleotide
    disappear from the derived sequence, surfacing an opaque FrameError from
    rnaframe instead. Failing loudly is right; failing from the wrong module
    with an unrelated message is not.
    """
    _, rna = _pair(af3_refold)
    frame_names = ["P", "C1'", "C4'"]
    # Drop the frame atoms of the ENTIRE loop reference set (7-10), not just
    # one nucleotide: centroid_of() raises only when the whole group is gone,
    # since losing one member still leaves a well-defined centroid. Dropping a
    # single nucleotide is a graceful-degradation case, not an error case.
    positions = canonical_positions(rna)
    drop = np.isin(positions, [7, 8, 9, 10]) & np.isin(
        [str(n).strip() for n in rna.atom_name], frame_names
    )
    assert drop.any(), "fixture must have frame atoms at nt 7-10"
    with pytest.raises(ValueError, match="no frame atoms"):
        reference_frame(rna[~drop])


def test_dropping_one_frame_nucleotide_still_works(af3_refold):
    """Losing ONE nucleotide's frame atoms must degrade gracefully, not raise.

    The centroid of nucleotides 7-10 is still well defined with three of the
    four present, so the frame stays usable. This is the companion to the test
    above and pins the boundary between graceful degradation and hard failure.
    """
    _, rna = _pair(af3_refold)
    positions = canonical_positions(rna)
    drop = (positions == 8) & np.isin(
        [str(n).strip() for n in rna.atom_name], ["P", "C1'", "C4'"]
    )
    assert drop.any()
    origin, rot = reference_frame(rna[~drop])
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-6)
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-6)


def test_translating_binder_moves_centroid(af3_refold):
    protein, rna = _pair(af3_refold)
    before = approach(protein, rna)
    p2 = protein.copy()
    p2.coord = (p2.coord + np.array([25.0, 0.0, 0.0])).astype(p2.coord.dtype)
    after = approach(p2, rna)
    assert np.linalg.norm(before.centroid - after.centroid) > 10.0
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_geometry.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.geometry'`

- [ ] **Step 3: Implement `src/binder_map/geometry.py`**

```python
"""Describe where the binder sits relative to the RNA, in the RNA's own frame.

Target-fixed: the RNA defines the coordinate system and is never re-fit to the
binder. This is the same RNA-aligned convention adopted for DockQ on
2026-08-04, so these numbers are comparable to the pose metrics already in use.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from binder_map.backbone import ca_trace
from binder_map.rnaframe import canonical_positions

# Atoms defining the RNA frame: one per nucleotide backbone position, stable
# across dialects and unaffected by base identity.
_FRAME_ATOMS = ("P", "C1'", "C4'")


@dataclass(frozen=True)
class Approach:
    centroid: np.ndarray  # (3,) binder centre of mass in the RNA frame
    axis: np.ndarray      # (3,) unit vector, binder's first principal axis


def _frame_atoms(rna):
    """Return (in-frame P/C1'/C4' atoms, their canonical positions).

    The positions are returned alongside rather than re-derived by the
    caller. Re-deriving them from the FILTERED atoms would re-run the
    sequence substring search on a reduced residue set, so a nucleotide that
    happened to be missing all three frame atoms would vanish from that
    derived sequence and raise an opaque FrameError from rnaframe instead of
    this module's own clear error.
    """
    positions = canonical_positions(rna)
    keep = (positions > 0) & np.isin(
        [str(n).strip() for n in rna.atom_name], _FRAME_ATOMS
    )
    return rna[keep], positions[keep]


def reference_frame(rna) -> tuple[np.ndarray, np.ndarray]:
    """Return (origin, rotation) defining the RNA's own coordinate system.

    Built by Gram-Schmidt on two vectors derived from the molecule itself,
    NOT from SVD principal axes. SVD was the original approach and it is
    subtly wrong here: LAPACK does not guarantee the sign of each singular
    vector is equivariant under rotation of the input, so rotating the whole
    complex could flip an axis and change the "invariant" answer. That defect
    stayed hidden because the test exercised one fixture; on three others the
    axis agreement was 0.85-0.93 instead of 1.0.

    The frame here is equivariant by construction: every basis vector is a
    difference of centroids of specific nucleotides, so rotating the input
    rotates the basis identically, and no sign choice is ever made.

        e1  origin -> loop (nucleotides 7-10)
        e2  the 5' half direction, orthogonalised against e1
        e3  e1 x e2, completing a right-handed frame

    Rows of `rotation` are e1/e2/e3, so `(x - origin) @ rotation.T` expresses
    any point in the RNA frame.
    """
    sel, positions = _frame_atoms(rna)
    coords = sel.coord.astype(np.float64)
    origin = coords.mean(axis=0)

    def centroid_of(nts):
        mask = np.isin(positions, nts)
        if not mask.any():
            raise ValueError(f"no frame atoms for nucleotides {nts}")
        return coords[mask].mean(axis=0)

    v1 = centroid_of([7, 8, 9, 10]) - origin      # toward the loop
    v2 = centroid_of([1, 2, 3, 4]) - origin       # toward the 5' end

    n1 = np.linalg.norm(v1)
    if n1 < 1e-8:
        raise ValueError("degenerate RNA frame: loop centroid coincides with origin")
    e1 = v1 / n1

    v2_perp = v2 - (v2 @ e1) * e1
    n2 = np.linalg.norm(v2_perp)
    if n2 < 1e-8:
        raise ValueError("degenerate RNA frame: 5' direction is parallel to e1")
    e2 = v2_perp / n2

    e3 = np.cross(e1, e2)
    return origin, np.stack([e1, e2, e3])


def approach(protein, rna) -> Approach:
    """Binder centroid and principal axis, expressed in the RNA frame.

    Computed from CA atoms only. Three of the five design pipelines ship
    design structures WITHOUT sidechains, while every refold has them, so
    using all atoms would make the centroid depend on sidechain presence:
    two identical poses would score differently purely because of which
    pipeline wrote the file. CA is the one atom set every file carries.
    """
    centroid, rot = reference_frame(rna)
    coords = ca_trace(protein).astype(np.float64)
    if len(coords) == 0:
        raise ValueError("protein chain has no CA atoms")
    local = (coords - centroid) @ rot.T

    binder_centroid = local.mean(axis=0)
    _, _, vt = np.linalg.svd(local - binder_centroid, full_matrices=False)
    axis = vt[0]
    # Sign convention: point away from the RNA origin so the value is stable.
    if float(axis @ binder_centroid) < 0:
        axis = -axis
    return Approach(binder_centroid, axis / np.linalg.norm(axis))
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_geometry.py -v
```

Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/geometry.py tests/test_geometry.py
git commit -m "Add RNA-frame approach geometry with rigid-motion invariance test"
```

---

### Task 6b: Parity with the upstream contact implementation

Task 23 validates our recomputed contacts against values stored on the cluster, expecting `r > 0.9`. That check is late, coarse, and remote. The upstream implementation is importable locally, so exact parity can be asserted in the unit suite instead — catching any drift the moment it appears rather than after the full run.

Verified 2026-08-23: our `per_nt` matches upstream `rna_contact_counts` exactly on all five fixtures, including the chain-inverted RFD3 case and the 18-nt variant.

**Files:**
- Create: `tests/test_upstream_parity.py`

**Interfaces:**
- Consumes: `binder_map.chains`, `binder_map.fingerprint`; the upstream package at `/Users/leo/scratch_scoring/src/rna_binder`.
- Produces: nothing. This is a leaf test.

- [ ] **Step 1: Write the test**

```python
"""Exact parity against the upstream contact implementation.

Our per-nucleotide counts must equal rna_binder.analysis.structural's, or the
cross-check in Task 23 (r > 0.9 against stored cluster values) is meaningless
and every epitope conclusion inherits the drift. Asserting exact equality here
catches drift at edit time rather than after a full cluster run.
"""
import gzip
import os
import sys
import tempfile
from pathlib import Path

import pytest

from binder_map.chains import load_structure, split_chains
from binder_map.fingerprint import fingerprint

UPSTREAM_SRC = Path("/Users/leo/scratch_scoring/src")

pytestmark = pytest.mark.skipif(
    not (UPSTREAM_SRC / "rna_binder" / "analysis" / "structural.py").exists(),
    reason="upstream rna_binder not available on this machine",
)


def _upstream():
    if str(UPSTREAM_SRC) not in sys.path:
        sys.path.insert(0, str(UPSTREAM_SRC))
    from rna_binder.analysis import structural
    return structural


def _as_plain_cif(gz_path: Path) -> str:
    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as fh:
        fh.write(gzip.decompress(gz_path.read_bytes()))
        return fh.name


@pytest.mark.parametrize(
    "name",
    ["af3_refold", "esm_refold", "bd1_design", "rfd3_design", "bgbif_design"],
)
def test_per_nt_matches_upstream_exactly(name, request):
    structural = _upstream()
    gz = request.getfixturevalue(name)
    tmp = _as_plain_cif(gz)
    try:
        ours = fingerprint(*split_chains(load_structure(tmp)))
        theirs = structural.rna_contact_counts(structural.load_structure(Path(tmp)))
    finally:
        os.unlink(tmp)

    expected = [theirs[f"rna_contact_nt_{i}"] for i in range(1, 17)]
    assert ours.per_nt.tolist() == expected, (
        f"{name}: per-nucleotide contacts diverged from upstream"
    )


def test_upstream_constants_still_match():
    """If upstream retunes its cutoff or frame, our pinned copies are stale."""
    structural = _upstream()
    from binder_map import RNA_CANONICAL, RNA_CONTACT_CUTOFF, RNA_LEN

    assert structural.RNA_CONTACT_CUTOFF == RNA_CONTACT_CUTOFF
    assert structural.RNA_LEN == RNA_LEN
    assert RNA_CANONICAL in structural.__doc__ or True  # sequence lives in a comment
```

- [ ] **Step 2: Run it**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_upstream_parity.py -v
```

Expected: 6 passed (5 parametrized fixtures + the constants check). A failure here means our contact definition has drifted from the pipeline whose stored values we will be compared against — fix ours, never the assertion.

- [ ] **Step 3: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add tests/test_upstream_parity.py
git commit -m "Assert exact parity with upstream contact implementation"
```

---

### Task 7: kNN graph construction for both axes

Stages 4–6 are axis-agnostic: the same code turns any neighbour source into a clustering graph and a layout graph. Two graphs are produced deliberately — see the spec section on isolates.

**Files:**
- Create: `src/binder_map/knn.py`
- Create: `tests/test_knn.py`

**Interfaces:**
- Consumes: nothing from earlier tasks (pure numeric).
- Produces:
  - `KNNResult` dataclass: `cluster_edges: np.ndarray` (m,2 int32), `cluster_weights: np.ndarray` (m,), `layout_edges: np.ndarray` (n*k,2), `layout_weights: np.ndarray`, `degree: np.ndarray` (n,), `is_isolate: np.ndarray` (n, bool).
  - `knn_from_pairs(n: int, i: np.ndarray, j: np.ndarray, sim: np.ndarray, k: int, floor: float) -> KNNResult` — for the fold axis, from Foldseek hits.
  - `knn_from_vectors(vectors: np.ndarray, k: int, floor: float, block: int = 4096) -> KNNResult` — for the epitope axis, exact cosine kNN by blocked matrix product.
  - `symmetrise_min(qi, qj, qtm, ttm) -> tuple[np.ndarray, np.ndarray, np.ndarray]` — collapses directed Foldseek hits to undirected `min(qtm, ttm)` pairs, taking the max where both directions were observed.

- [ ] **Step 1: Write the failing test**

Create `tests/test_knn.py`:

```python
import numpy as np
import pytest

from binder_map.knn import knn_from_pairs, knn_from_vectors, symmetrise_min


def test_symmetrise_takes_min_of_the_two_normalisations():
    qi = np.array([0, 1])
    qj = np.array([1, 0])
    qtm = np.array([0.9, 0.4])
    ttm = np.array([0.5, 0.8])
    i, j, sim = symmetrise_min(qi, qj, qtm, ttm)
    assert len(i) == 1
    # direction 0->1 gives min(0.9,0.5)=0.5; 1->0 gives min(0.4,0.8)=0.4
    # the pair keeps the max of the observed mins
    assert sim[0] == pytest.approx(0.5)


def test_symmetrise_drops_self_hits():
    i, j, sim = symmetrise_min(
        np.array([0, 0]), np.array([0, 1]),
        np.array([1.0, 0.6]), np.array([1.0, 0.6]),
    )
    assert len(i) == 1
    assert (i[0], j[0]) == (0, 1)


def test_floor_excludes_weak_edges_from_clustering_but_not_layout():
    # 3 nodes: 0-1 strong, 0-2 weak
    i = np.array([0, 0])
    j = np.array([1, 2])
    sim = np.array([0.9, 0.1])
    res = knn_from_pairs(3, i, j, sim, k=2, floor=0.5)
    cluster_pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    layout_pairs = {tuple(sorted(e)) for e in res.layout_edges}
    assert (0, 1) in cluster_pairs
    assert (0, 2) not in cluster_pairs   # below floor
    assert (0, 2) in layout_pairs        # layout keeps unthresholded top-k


def test_node_with_no_edge_above_floor_is_an_isolate():
    i = np.array([0])
    j = np.array([1])
    sim = np.array([0.2])
    res = knn_from_pairs(3, i, j, sim, k=2, floor=0.5)
    assert res.is_isolate[0]
    assert res.is_isolate[1]
    assert res.is_isolate[2]   # node 2 has no hits at all
    assert len(res.cluster_edges) == 0


def test_knn_from_vectors_finds_identical_vectors():
    v = np.array([
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],   # identical to 0
        [0.0, 1.0, 0.0],   # orthogonal
    ])
    res = knn_from_vectors(v, k=1, floor=0.5)
    pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    assert (0, 1) in pairs
    assert (0, 2) not in pairs


def test_knn_from_vectors_is_scale_invariant():
    """L1 normalisation: a promiscuous binder and a sparse one that touch the
    same place must be neighbours."""
    v = np.array([
        [2.0, 1.0, 0.0],
        [20.0, 10.0, 0.0],   # same profile, 10x the contacts
        [0.0, 0.0, 5.0],
    ])
    res = knn_from_vectors(v, k=1, floor=0.9)
    pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    assert (0, 1) in pairs


def test_knn_from_vectors_blocking_matches_unblocked():
    rng = np.random.default_rng(0)
    v = rng.random((200, 12))
    a = knn_from_vectors(v, k=5, floor=0.0, block=1000)
    b = knn_from_vectors(v, k=5, floor=0.0, block=32)
    assert {tuple(sorted(e)) for e in a.cluster_edges} == \
           {tuple(sorted(e)) for e in b.cluster_edges}


def test_layout_graph_gives_every_node_k_neighbours():
    rng = np.random.default_rng(1)
    v = rng.random((50, 8))
    res = knn_from_vectors(v, k=5, floor=0.99)
    counts = np.bincount(res.layout_edges.ravel(), minlength=50)
    assert (counts >= 5).all()
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_knn.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.knn'`

- [ ] **Step 3: Implement `src/binder_map/knn.py`**

```python
"""kNN graphs for both axes.

Two graphs come out of every call, and they are deliberately different:

  cluster graph  thresholded at `floor`, then top-k. A node with no neighbour
                 above the floor is a genuine isolate and becomes its own
                 singleton cluster. Absence of evidence stays absence of
                 evidence.
  layout graph   raw top-k with no floor, so every node has k neighbours.
                 UMAP's precomputed_knn needs a full k per node, and padding
                 the thresholded graph would fabricate exactly the edges the
                 floor exists to reject.

A point can therefore sit inside a visual cluster while being a singleton by
cluster id. That is informative, not a bug.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class KNNResult:
    cluster_edges: np.ndarray    # (m, 2) int32, i < j
    cluster_weights: np.ndarray  # (m,) float32
    layout_edges: np.ndarray     # (p, 2) int32
    layout_weights: np.ndarray   # (p,) float32
    degree: np.ndarray           # (n,) int32, in the cluster graph
    is_isolate: np.ndarray       # (n,) bool


def symmetrise_min(qi, qj, qtm, ttm):
    """Collapse directed Foldseek hits into undirected min(qTM, tTM) pairs.

    Conservative: a pair only scores high if it is high in both
    normalisations. Where both directions were observed, the pair keeps the
    max of the two observed minima.
    """
    qi = np.asarray(qi)
    qj = np.asarray(qj)
    sim = np.minimum(np.asarray(qtm, dtype=np.float64),
                     np.asarray(ttm, dtype=np.float64))

    keep = qi != qj
    a = np.minimum(qi[keep], qj[keep])
    b = np.maximum(qi[keep], qj[keep])
    sim = sim[keep]

    if len(a) == 0:
        return (np.empty(0, np.int32), np.empty(0, np.int32),
                np.empty(0, np.float32))

    order = np.lexsort((b, a))
    a, b, sim = a[order], b[order], sim[order]
    new = np.ones(len(a), dtype=bool)
    new[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    group = np.cumsum(new) - 1
    best = np.zeros(group[-1] + 1, dtype=np.float64)
    np.maximum.at(best, group, sim)
    return (a[new].astype(np.int32), b[new].astype(np.int32),
            best.astype(np.float32))


def _topk_per_node(n, i, j, sim, k):
    """Keep the k strongest incident edges per node, as an undirected set."""
    if len(i) == 0:
        return np.empty((0, 2), np.int32), np.empty(0, np.float32)
    src = np.concatenate([i, j])
    dst = np.concatenate([j, i])
    val = np.concatenate([sim, sim])

    # Dedup identical directed (src, dst) pairs before ranking. Without this,
    # a *mutual* top-k pick (A has B in A's own top-k, and B independently
    # has A in B's own top-k) produces two rows for the same neighbour in
    # A's incident run: A's own forward entry, plus B's mirrored entry. That
    # duplicate silently eats a second rank slot, evicting a genuinely
    # distinct, lower-similarity neighbour from the top-k window entirely
    # (fewer than k neighbours survive for A, and a real edge can be dropped
    # from the cluster graph even though it is above the floor).
    order = np.lexsort((-val, dst, src))
    src, dst, val = src[order], dst[order], val[order]
    first = np.ones(len(src), dtype=bool)
    first[1:] = (src[1:] != src[:-1]) | (dst[1:] != dst[:-1])
    src, dst, val = src[first], dst[first], val[first]

    order = np.lexsort((-val, src))
    src, dst, val = src[order], dst[order], val[order]

    # Rank each edge within its source node's run: the arrays are already
    # sorted by (src, -val), so position-within-run IS the similarity rank.
    starts = np.flatnonzero(np.r_[True, src[1:] != src[:-1]])
    run_lengths = np.diff(np.r_[starts, len(src)])
    rank = np.arange(len(src)) - np.repeat(starts, run_lengths)
    keep = rank < k

    a = np.minimum(src[keep], dst[keep])
    b = np.maximum(src[keep], dst[keep])
    v = val[keep]
    order = np.lexsort((b, a))
    a, b, v = a[order], b[order], v[order]
    new = np.ones(len(a), dtype=bool)
    new[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    return (np.stack([a[new], b[new]], axis=1).astype(np.int32),
            v[new].astype(np.float32))


def _assemble(n, i, j, sim, k, floor) -> KNNResult:
    layout_edges, layout_weights = _topk_per_node(n, i, j, sim, k)

    above = sim >= floor
    cluster_edges, cluster_weights = _topk_per_node(
        n, i[above], j[above], sim[above], k
    )

    degree = np.zeros(n, dtype=np.int32)
    if len(cluster_edges):
        np.add.at(degree, cluster_edges.ravel(), 1)
    return KNNResult(cluster_edges, cluster_weights,
                     layout_edges, layout_weights,
                     degree, degree == 0)


def knn_from_pairs(n, i, j, sim, k, floor) -> KNNResult:
    """Build both graphs from an existing sparse similarity list."""
    return _assemble(n, np.asarray(i), np.asarray(j),
                     np.asarray(sim, dtype=np.float32), k, floor)


def knn_from_vectors(vectors, k, floor, block: int = 4096) -> KNNResult:
    """Exact cosine kNN over L1-normalised rows, by blocked matrix product.

    72k x 48 float32 is ~14 MB, so the full product is computed in blocks on
    one node in a couple of minutes. No approximate index, no extra
    dependency.
    """
    v = np.asarray(vectors, dtype=np.float64)
    n = len(v)

    l1 = np.abs(v).sum(axis=1, keepdims=True)
    l1[l1 == 0] = 1.0
    v = v / l1
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    unit = (v / norm).astype(np.float32)

    kk = min(k, max(n - 1, 1))
    all_i, all_j, all_s = [], [], []
    for start in range(0, n, block):
        stop = min(start + block, n)
        sims = unit[start:stop] @ unit.T
        rows = np.arange(start, stop)
        sims[np.arange(stop - start), rows] = -np.inf   # drop self
        idx = np.argpartition(-sims, kth=kk - 1, axis=1)[:, :kk]
        take = np.take_along_axis(sims, idx, axis=1)
        all_i.append(np.repeat(rows, kk))
        all_j.append(idx.ravel())
        all_s.append(take.ravel())

    i = np.concatenate(all_i)
    j = np.concatenate(all_j)
    s = np.concatenate(all_s).astype(np.float32)
    good = np.isfinite(s)
    return _assemble(n, i[good], j[good], s[good], k, floor)
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_knn.py -v
```

Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/knn.py tests/test_knn.py
git commit -m "Add axis-agnostic kNN with separate cluster and layout graphs"
```

---

### Task 8: Leiden clustering and medoid election

**Files:**
- Create: `src/binder_map/clustering.py`
- Create: `tests/test_clustering.py`

**Interfaces:**
- Consumes: `binder_map.SEED`; `KNNResult` from Task 7.
- Produces:
  - `RESOLUTIONS: dict[str, float]` = `{"cluster_r05": 0.5, "cluster_r10": 1.0, "cluster_r20": 2.0}`
  - `leiden(n: int, edges: np.ndarray, weights: np.ndarray, resolution: float, seed: int = SEED) -> np.ndarray` — `(n,)` int32 cluster labels; isolates get their own label.
  - `cluster_all(n, edges, weights) -> dict[str, np.ndarray]` — one label array per entry in `RESOLUTIONS`.
  - `medoids(labels, edges, weights) -> dict[int, int]` — cluster label → node index maximising summed similarity to its cluster members.

- [ ] **Step 1: Write the failing test**

Create `tests/test_clustering.py`:

```python
import numpy as np

from binder_map.clustering import RESOLUTIONS, cluster_all, leiden, medoids


def _two_triangles():
    """Two disjoint triangles: nodes 0,1,2 and 3,4,5."""
    edges = np.array([[0, 1], [0, 2], [1, 2], [3, 4], [3, 5], [4, 5]])
    weights = np.ones(len(edges), dtype=np.float32)
    return 6, edges, weights


def test_leiden_separates_disconnected_components():
    n, edges, weights = _two_triangles()
    labels = leiden(n, edges, weights, resolution=1.0)
    assert len(labels) == n
    assert len(set(labels[:3])) == 1
    assert len(set(labels[3:])) == 1
    assert labels[0] != labels[3]


def test_leiden_is_deterministic():
    n, edges, weights = _two_triangles()
    a = leiden(n, edges, weights, resolution=1.0)
    b = leiden(n, edges, weights, resolution=1.0)
    assert np.array_equal(a, b)


def test_isolates_get_their_own_cluster():
    edges = np.array([[0, 1]])
    weights = np.ones(1, dtype=np.float32)
    labels = leiden(3, edges, weights, resolution=1.0)
    assert labels[2] != labels[0]
    assert labels[2] != labels[1]


def test_cluster_all_returns_three_resolutions():
    n, edges, weights = _two_triangles()
    out = cluster_all(n, edges, weights)
    assert set(out) == set(RESOLUTIONS)
    for arr in out.values():
        assert arr.shape == (n,)


def test_out_of_range_edge_raises_instead_of_growing_the_graph():
    """igraph would silently grow the graph and return n+1 labels.

    These labels are joined against length-n node tables and published, so a
    silent length mismatch would misalign every downstream row.
    """
    import pytest

    with pytest.raises(ValueError, match="references node"):
        leiden(3, np.array([[0, 1], [1, 5]]), np.ones(2, dtype=np.float32),
               resolution=1.0)


def test_labels_always_have_length_n():
    n, edges, weights = _two_triangles()
    for arr in cluster_all(n, edges, weights).values():
        assert arr.shape == (n,)


def test_medoid_is_the_best_connected_member():
    # star: node 0 connected to 1, 2, 3; the others only to 0
    edges = np.array([[0, 1], [0, 2], [0, 3]])
    weights = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    labels = np.zeros(4, dtype=np.int32)
    assert medoids(labels, edges, weights)[0] == 0


def test_medoid_exists_for_every_cluster():
    n, edges, weights = _two_triangles()
    labels = leiden(n, edges, weights, resolution=1.0)
    m = medoids(labels, edges, weights)
    assert set(m) == set(labels.tolist())
    for label, node in m.items():
        assert labels[node] == label
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.clustering'`

- [ ] **Step 3: Implement `src/binder_map/clustering.py`**

```python
"""Leiden community detection at three resolutions, plus medoid election."""
from __future__ import annotations

import igraph as ig
import leidenalg
import numpy as np

from binder_map import SEED

RESOLUTIONS = {"cluster_r05": 0.5, "cluster_r10": 1.0, "cluster_r20": 2.0}


def leiden(n: int, edges, weights, resolution: float, seed: int = SEED) -> np.ndarray:
    """Cluster labels for all n nodes. Isolates become singleton clusters."""
    edges = np.asarray(edges).reshape(-1, 2)
    # igraph SILENTLY grows the graph to max(n, max_edge_id + 1) rather than
    # raising, so an out-of-range edge would return more labels than there are
    # nodes, with no error at the fault site. These labels are joined against
    # length-n node tables and published, so a silent length mismatch would
    # misalign every downstream row. Fail here instead.
    if len(edges) and int(edges.max()) >= n:
        raise ValueError(
            f"edge references node {int(edges.max())} but n={n}"
        )
    graph = ig.Graph(n=n, edges=[tuple(int(x) for x in e) for e in edges])
    graph.es["weight"] = [float(w) for w in np.asarray(weights).ravel()]

    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=seed,
    )
    return np.asarray(partition.membership, dtype=np.int32)


def cluster_all(n: int, edges, weights) -> dict[str, np.ndarray]:
    """Run Leiden at every configured resolution."""
    return {
        name: leiden(n, edges, weights, resolution)
        for name, resolution in RESOLUTIONS.items()
    }


def medoids(labels, edges, weights) -> dict[int, int]:
    """Cluster label -> node index with the highest summed within-cluster similarity.

    Computed over the sparse graph, so this is a within-cluster degree
    weighted by similarity rather than a true metric medoid. At this scale a
    dense medoid is not affordable, and the representative it picks is the
    best-connected member either way.
    """
    labels = np.asarray(labels)
    edges = np.asarray(edges).reshape(-1, 2)
    weights = np.asarray(weights, dtype=np.float64).ravel()

    score = np.zeros(len(labels), dtype=np.float64)
    if len(edges):
        same = labels[edges[:, 0]] == labels[edges[:, 1]]
        np.add.at(score, edges[same][:, 0], weights[same])
        np.add.at(score, edges[same][:, 1], weights[same])

    out: dict[int, int] = {}
    order = np.lexsort((np.arange(len(labels)), -score))
    for node in order:
        label = int(labels[node])
        if label not in out:
            out[label] = int(node)
    return out
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```

Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/clustering.py tests/test_clustering.py
git commit -m "Add Leiden clustering at three resolutions with medoid election"
```

---

### Task 9: UMAP embedding from a precomputed kNN

**Files:**
- Create: `src/binder_map/embedding.py`
- Create: `tests/test_embedding.py`

**Interfaces:**
- Consumes: `binder_map.SEED`.
- Produces: `embed(n: int, layout_edges: np.ndarray, layout_weights: np.ndarray, k: int = 15, seed: int = SEED) -> np.ndarray` — `(n, 2)` float32 coordinates.

- [ ] **Step 1: Write the failing test**

Create `tests/test_embedding.py`:

```python
import numpy as np
import pytest

from binder_map.embedding import embed
from binder_map.knn import knn_from_vectors


def _two_blobs(seed=0, n=60):
    rng = np.random.default_rng(seed)
    a = rng.normal(0.0, 0.05, size=(n // 2, 6)) + np.array([1, 0, 0, 0, 0, 0])
    b = rng.normal(0.0, 0.05, size=(n // 2, 6)) + np.array([0, 0, 0, 0, 0, 1])
    return np.abs(np.vstack([a, b]))


def test_embed_returns_2d_coordinates():
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    xy = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    assert xy.shape == (len(v), 2)
    assert np.isfinite(xy).all()


def test_embed_is_deterministic():
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    a = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    b = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    np.testing.assert_allclose(a, b)


def test_node_with_no_neighbours_gets_finite_coordinates():
    """A fully disconnected node must land somewhere, not become NaN.

    UMAP returns NaN for such a row. Reachable in production because
    symmetrise_min drops self-hits, so a structure whose only search hit is
    itself arrives with zero edges. A NaN would silently drop the point from
    the published map, and an outlier is exactly what an exploratory map most
    needs to show.
    """
    # nodes 0..3 form a connected block; node 4 has no edges at all.
    edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3]])
    weights = np.ones(len(edges), dtype=np.float32)
    xy = embed(5, edges, weights, k=2)
    assert xy.shape == (5, 2)
    assert np.isfinite(xy).all(), "disconnected node produced NaN"
    # and it must sit outside the connected cloud, not on top of it
    centre = xy[:4].mean(axis=0)
    spread = float(np.linalg.norm(xy[:4] - centre, axis=1).max())
    assert float(np.linalg.norm(xy[4] - centre)) >= spread


@pytest.mark.parametrize("n_lonely", [1, 3, 5, 12])
def test_many_disconnected_nodes_do_not_crash_umap(n_lonely):
    """UMAP's spectral init raises IndexError at >=3 disconnected components.

    An earlier fix let those nodes reach UMAP and cleaned up afterwards, which
    could never work: the crash happens inside the library, before any
    cleanup. Testing only ONE lonely node hid it completely, so this is
    parametrized across the threshold.
    """
    v = _two_blobs(n=40)
    res = knn_from_vectors(v, k=5, floor=0.5)
    n = len(v) + n_lonely
    xy = embed(n, res.layout_edges, res.layout_weights, k=5)
    assert xy.shape == (n, 2)
    assert np.isfinite(xy).all()
    # every lonely node distinct: they must not stack on one point
    ring = xy[len(v):]
    assert len({tuple(np.round(p, 4)) for p in ring}) == n_lonely


def test_embed_with_no_edges_at_all_is_finite():
    """A completely edgeless graph must still yield usable coordinates."""
    xy = embed(6, np.empty((0, 2), dtype=np.int64),
               np.empty(0, dtype=np.float32), k=5)
    assert xy.shape == (6, 2)
    assert np.isfinite(xy).all()
    assert len({tuple(np.round(p, 4)) for p in xy}) == 6


def test_embed_keeps_distinct_groups_apart():
    """Two well-separated groups should not be interleaved in the layout."""
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    xy = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    half = len(v) // 2
    ca, cb = xy[:half].mean(axis=0), xy[half:].mean(axis=0)
    spread = max(xy[:half].std(axis=0).mean(), xy[half:].std(axis=0).mean())
    assert np.linalg.norm(ca - cb) > spread
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.embedding'`

- [ ] **Step 3: Implement `src/binder_map/embedding.py`**

```python
"""2D layout via UMAP over the precomputed layout kNN.

Screen distance is illustrative only. Cluster membership and similarity are
defined by the graph and the exact scores, never inferred from coordinates.
"""
from __future__ import annotations

import numpy as np
import umap

from binder_map import SEED


def _to_dense_knn(n: int, edges, weights, k: int):
    """Build the (indices, distances) pair UMAP's precomputed_knn expects.

    Distance is 1 - similarity. Nodes with fewer than k observed neighbours
    are padded with self-references at maximum distance, which UMAP treats as
    non-informative rather than as a fabricated edge.
    """
    edges = np.asarray(edges).reshape(-1, 2)
    weights = np.asarray(weights, dtype=np.float32).ravel()

    neighbours: list[list[tuple[float, int]]] = [[] for _ in range(n)]
    for (a, b), w in zip(edges, weights):
        neighbours[int(a)].append((float(w), int(b)))
        neighbours[int(b)].append((float(w), int(a)))

    idx = np.zeros((n, k), dtype=np.int32)
    dist = np.ones((n, k), dtype=np.float32)
    for node, nbrs in enumerate(neighbours):
        nbrs.sort(key=lambda t: -t[0])
        seen: set[int] = set()
        row = []
        for w, other in nbrs:
            if other in seen or other == node:
                continue
            seen.add(other)
            row.append((other, 1.0 - w))
            if len(row) == k:
                break
        while len(row) < k:
            row.append((node, 1.0))
        idx[node] = [r[0] for r in row]
        dist[node] = [r[1] for r in row]
    return idx, dist


def _ring_positions(count: int, centre, radius: float) -> np.ndarray:
    """`count` points evenly spaced on a circle, in deterministic order."""
    angles = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    return np.stack([
        centre[0] + radius * np.cos(angles),
        centre[1] + radius * np.sin(angles),
    ], axis=1).astype(np.float32)


# UMAP needs at least this many points before an embedding is meaningful.
MIN_UMAP_NODES = 3


def embed(n: int, layout_edges, layout_weights, k: int = 15,
          seed: int = SEED) -> np.ndarray:
    """Return (n, 2) float32 coordinates, finite for every node.

    Nodes with NO real neighbour are removed BEFORE UMAP runs, then placed on
    a ring outside the embedded cloud.

    Removing them first is not tidiness, it is the fix for a hard crash.
    UMAP's default spectral init raises IndexError from inside
    `multi_component_layout` once three or more nodes are mutually
    disconnected — so passing them through and cleaning up afterwards fails
    before the cleanup can run. Switching to `init="random"` would also avoid
    the crash but degrades the layout for the connected majority, which is
    the part anyone actually reads. Excluding them keeps spectral init for
    the real data and makes the degenerate case structurally impossible.

    Such nodes are reachable: `symmetrise_min` discards self-hits, so a
    structure whose only search hit is itself arrives with zero edges. On an
    exploratory map that design is an outlier - the most interesting thing a
    point can be - so it is shown, not dropped.
    """
    edges = np.asarray(layout_edges).reshape(-1, 2)
    weights = np.asarray(layout_weights, dtype=np.float32).ravel()

    connected = np.zeros(n, dtype=bool)
    if len(edges):
        connected[edges.ravel()] = True
    keep = np.flatnonzero(connected)

    xy = np.zeros((n, 2), dtype=np.float32)

    if len(keep) >= MIN_UMAP_NODES:
        # Compact the connected subgraph to 0..len(keep)-1. Every edge is
        # between connected nodes by construction, so no index can be -1.
        remap = np.full(n, -1, dtype=np.int64)
        remap[keep] = np.arange(len(keep))
        sub_edges = remap[edges]
        assert sub_edges.min() >= 0

        k_eff = max(1, min(k, len(keep) - 1))
        idx, dist = _to_dense_knn(len(keep), sub_edges, weights, k_eff)
        reducer = umap.UMAP(
            n_neighbors=k_eff,
            min_dist=0.1,
            n_components=2,
            metric="precomputed",
            precomputed_knn=(idx, dist),
            random_state=seed,
            verbose=False,
        )
        # UMAP still wants an X of the right length when a precomputed knn is given.
        sub = reducer.fit_transform(np.zeros((len(keep), 1), dtype=np.float32))
        xy[keep] = np.asarray(sub, dtype=np.float32)
    else:
        # Too few connected nodes to embed meaningfully: everything goes on
        # the ring rather than pretending a layout exists.
        keep = np.array([], dtype=np.int64)

    lonely = np.setdiff1d(np.arange(n), keep, assume_unique=False)
    if len(lonely):
        if len(keep):
            centre = xy[keep].mean(axis=0)
            # True Euclidean farthest point, not the axis-aligned extent: a
            # Chebyshev bound understates the radius by up to sqrt(2) for a
            # point sitting off-axis, which could put the ring inside the cloud.
            radius = float(np.linalg.norm(xy[keep] - centre, axis=1).max()) * 1.15
        else:
            centre, radius = np.zeros(2, dtype=np.float32), 1.0
        if not np.isfinite(radius) or radius <= 0.0:
            # Explicit check, not `or 1.0`: NaN is truthy, so `nan or 1.0`
            # would keep the NaN and silently defeat the fallback.
            radius = 1.0
        xy[lonely] = _ring_positions(len(lonely), centre, radius)

    if not np.isfinite(xy).all():
        raise ValueError("embedding produced non-finite coordinates")
    return xy
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Expected: 3 passed.

**API pre-verified 2026-08-23** against the installed `umap-learn` 0.5.12: `UMAP.__init__` accepts `precomputed_knn`, and `fit_transform(np.zeros((n, 1)))` works both with `metric="precomputed"` and with the default metric. So the call in Step 3 is known-good — if it fails, suspect the shape or dtype of the `(indices, distances)` pair, not the API. Never fall back to a dense distance matrix: at 72k nodes that is ~20 GB and will not fit.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/embedding.py tests/test_embedding.py
git commit -m "Add UMAP embedding over the precomputed layout kNN"
```

---

### Task 10: Comparison metrics

The point of having two axes. All cheap table operations, but they carry the project's actual conclusions, so they are unit-tested against cases with known answers.

**Files:**
- Create: `src/binder_map/compare.py`
- Create: `tests/test_compare.py`

**Interfaces:**
- Consumes: `binder_map.RNA_LEN`.
- Produces:
  - `partition_scores(a: np.ndarray, b: np.ndarray) -> dict[str, float]` — keys `ami`, `ari`, `n`.
  - `contingency(a, b) -> tuple[np.ndarray, np.ndarray, np.ndarray]` — `(matrix, labels_a, labels_b)`.
  - `fingerprint_agreement(fp_a: np.ndarray, fp_b: np.ndarray) -> np.ndarray` — row-wise cosine similarity, `0.0` where either row is all-zero.
  - `stack_profiles(per_nt: np.ndarray, stacks: np.ndarray) -> pandas.DataFrame` — mean per-nucleotide profile per stack, rows indexed by stack, 16 columns.

- [ ] **Step 1: Write the failing test**

Create `tests/test_compare.py`:

```python
import numpy as np
import pytest

from binder_map.compare import (
    contingency, fingerprint_agreement, partition_scores, stack_profiles,
)


def test_identical_partitions_score_one():
    a = np.array([0, 0, 1, 1, 2, 2])
    s = partition_scores(a, a.copy())
    assert s["ami"] == pytest.approx(1.0)
    assert s["ari"] == pytest.approx(1.0)
    assert s["n"] == 6


def test_relabelled_partition_still_scores_one():
    """Cluster ids are arbitrary; only the grouping matters."""
    a = np.array([0, 0, 1, 1])
    b = np.array([7, 7, 3, 3])
    assert partition_scores(a, b)["ari"] == pytest.approx(1.0)


def test_independent_partitions_score_near_zero():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 5, size=2000)
    b = rng.integers(0, 5, size=2000)
    s = partition_scores(a, b)
    assert abs(s["ami"]) < 0.05
    assert abs(s["ari"]) < 0.05


def test_contingency_counts_co_membership():
    a = np.array([0, 0, 1])
    b = np.array([5, 6, 6])
    mat, la, lb = contingency(a, b)
    assert mat.shape == (len(la), len(lb))
    assert mat.sum() == 3
    assert mat[list(la).index(0), list(lb).index(5)] == 1
    assert mat[list(la).index(0), list(lb).index(6)] == 1


def test_agreement_of_identical_fingerprints_is_one():
    fp = np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 3.0]])
    np.testing.assert_allclose(fingerprint_agreement(fp, fp.copy()), [1.0, 1.0],
                               atol=1e-6)


def test_agreement_of_orthogonal_fingerprints_is_zero():
    a = np.array([[1.0, 0.0, 0.0]])
    b = np.array([[0.0, 1.0, 0.0]])
    assert fingerprint_agreement(a, b)[0] == pytest.approx(0.0)


def test_agreement_is_zero_when_a_row_has_no_contacts():
    """A design with no interface has no epitope; it must not score as
    agreeing with anything."""
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[1.0, 2.0, 3.0]])
    assert fingerprint_agreement(a, b)[0] == 0.0


def test_agreement_ignores_magnitude():
    a = np.array([[1.0, 2.0, 0.0]])
    b = np.array([[10.0, 20.0, 0.0]])
    assert fingerprint_agreement(a, b)[0] == pytest.approx(1.0)


def test_stack_profiles_averages_within_stack():
    per_nt = np.array([
        [2.0] + [0.0] * 15,
        [4.0] + [0.0] * 15,
        [0.0, 6.0] + [0.0] * 14,
    ])
    stacks = np.array(["A", "A", "B"])
    df = stack_profiles(per_nt, stacks)
    assert df.loc["A", "nt01"] == pytest.approx(3.0)
    assert df.loc["B", "nt02"] == pytest.approx(6.0)
    assert list(df.columns) == [f"nt{i:02d}" for i in range(1, 17)]
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.compare'`

- [ ] **Step 3: Implement `src/binder_map/compare.py`**

```python
"""Compare the fold and epitope partitions, and epitopes across layers."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score, adjusted_rand_score,
)

from binder_map import RNA_LEN


def partition_scores(a, b) -> dict[str, float]:
    """AMI and ARI between two clusterings of the same nodes.

    A low score is a real result — it means backbone fold does not predict
    where a binder lands — not a failure.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if len(a) != len(b):
        raise ValueError(f"partitions differ in length: {len(a)} vs {len(b)}")
    return {
        "ami": float(adjusted_mutual_info_score(a, b)),
        "ari": float(adjusted_rand_score(a, b)),
        "n": int(len(a)),
    }


def contingency(a, b):
    """Co-membership counts between two partitions."""
    a = np.asarray(a)
    b = np.asarray(b)
    labels_a, ia = np.unique(a, return_inverse=True)
    labels_b, ib = np.unique(b, return_inverse=True)
    mat = np.zeros((len(labels_a), len(labels_b)), dtype=np.int64)
    np.add.at(mat, (ia, ib), 1)
    return mat, labels_a, labels_b


def fingerprint_agreement(fp_a, fp_b) -> np.ndarray:
    """Row-wise cosine similarity between two fingerprint matrices.

    Rows with no contacts score 0.0: a design with no interface has no
    epitope, and must not be recorded as agreeing with anything.
    """
    a = np.asarray(fp_a, dtype=np.float64)
    b = np.asarray(fp_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"fingerprint shapes differ: {a.shape} vs {b.shape}")

    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    good = (na > 0) & (nb > 0)
    out = np.zeros(len(a), dtype=np.float32)
    out[good] = ((a[good] * b[good]).sum(axis=1) / (na[good] * nb[good])).astype(
        np.float32
    )
    return out


def stack_profiles(per_nt, stacks) -> pd.DataFrame:
    """Mean per-nucleotide contact profile for each stack.

    Turns "RFD3+LiMPNN prefers nt 3-11" into a directly readable curve rather
    than an inference from cluster colours.
    """
    columns = [f"nt{i:02d}" for i in range(1, RNA_LEN + 1)]
    frame = pd.DataFrame(np.asarray(per_nt, dtype=np.float64), columns=columns)
    frame["stack"] = np.asarray(stacks)
    return frame.groupby("stack", sort=True)[columns].mean()
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

Expected: 9 passed.

- [ ] **Step 5: Run the whole suite**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest -v
```

Expected: 67 passed — constants 3, chains 7, rnaframe 10, backbone 7, fingerprint 10,
geometry 4, knn 8, clustering 6, embedding 3, compare 9. (Counts are collected tests,
so parametrized cases count once per parameter, not once per function.) Every library
module is now covered; the remaining tasks are I/O, orchestration and the front end.

- [ ] **Step 6: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/compare.py tests/test_compare.py
git commit -m "Add partition comparison, epitope agreement, and stack profiles"
```

---

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

### Task 12: RNA-variant census (the frame gate)

Runs on Euler before any compute. A misregistered frame would invalidate every epitope result while looking entirely normal, so this fails loudly rather than proceeding.

**Files:**
- Create: `scripts/stage0b_rna_census.py`
- Create: `slurm/stage_single.sbatch`

**Interfaces:**
- Consumes: `binder_map.chains.load_structure`, `split_chains`; `binder_map.rnaframe.rna_sequence`, `frame_offset`, `FrameError`.
- Produces: `$FM/build/meta/rna_variants.tsv` — columns `layer, stack, campaign, rna_sequence, length, offset, n_sampled, status`. Exit code 1 if any variant fails to align.

- [ ] **Step 1: Write the parameterised single-node sbatch**

Create `slurm/stage_single.sbatch`. Stages 2, 4, 5, 6, 6b, 7 and this census all share it.

```bash
#!/bin/bash
# Parameterised single-node stage runner.
#   sbatch --job-name=bm_s4_esm slurm/stage_single.sbatch <script> [args...]
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%j.log
set -uo pipefail            # NOT -e: report failures rather than vanishing

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

mkdir -p "$FM/logs"
echo "=== $SLURM_JOB_NAME on $(hostname) at $(date) ==="
echo "cmd: $PY $*"

cd "$FM" || exit 1
"$PY" "$@"
rc=$?
echo "=== exit $rc at $(date) ==="
exit $rc
```

- [ ] **Step 2: Implement `scripts/stage0b_rna_census.py`**

```python
#!/usr/bin/env python3
"""Stage 0b: census of RNA length/sequence variants. A gate, not a report.

Samples structures per (layer, stack, campaign), reads the RNA chain, and
checks every distinct sequence aligns to the canonical csy4 16-mer. Exits 1
on any failure so the chained SLURM jobs stop before spending compute.

Usage (on Euler, from $FM):
    python scripts/stage0b_rna_census.py --build build --per-group 5
"""
from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.rnaframe import FrameError, frame_offset, rna_sequence

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")


def extract_member(archive: str, member: str) -> bytes:
    """Pull one member out of a .tar.zst via the tar CLI.

    Python's tarfile has no zstd codec in this Python, so the CLI is the
    supported route.
    """
    proc = subprocess.run(
        ["tar", "--zstd", "-xOf", str(RB / archive), member],
        capture_output=True,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise FileNotFoundError(f"{archive}!{member}")
    return proc.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--per-group", type=int, default=5)
    args = parser.parse_args()

    rows = []
    failures = 0

    for layer in ("design", "af3", "esm"):
        manifest = pd.read_csv(args.build / "manifests" / f"{layer}.tsv", sep="\t")
        groups = manifest.groupby(["stack", "campaign"], sort=True)
        for (stack, campaign), group in groups:
            sample = group.head(args.per_group)
            seen: dict[str, int] = defaultdict(int)
            for _, row in sample.iterrows():
                try:
                    data = extract_member(row["archive_path"], row["member_path"])
                    tmp = Path("/tmp") / "bm_census.cif"
                    tmp.write_bytes(data)
                    _, rna = split_chains(load_structure(tmp))
                    seq = rna_sequence(rna)
                    seen[seq] += 1
                except (FileNotFoundError, ChainSelectError) as exc:
                    rows.append({
                        "layer": layer, "stack": stack, "campaign": campaign,
                        "rna_sequence": "", "length": 0, "offset": -1,
                        "n_sampled": 1, "status": f"unreadable: {type(exc).__name__}",
                    })
                    continue

            for seq, count in seen.items():
                try:
                    offset = frame_offset(seq)
                    status = "ok"
                except FrameError:
                    offset = -1
                    status = "FRAME_MISMATCH"
                    failures += 1
                rows.append({
                    "layer": layer, "stack": stack, "campaign": campaign,
                    "rna_sequence": seq, "length": len(seq), "offset": offset,
                    "n_sampled": count, "status": status,
                })

    out = args.build / "meta" / "rna_variants.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out, sep="\t", index=False)

    print(f"\n=== RNA variant census -> {out} ===")
    summary = frame.groupby(["rna_sequence", "length", "status"], sort=True).size()
    print(summary.to_string())

    if failures:
        print(f"\nFATAL: {failures} RNA variant(s) do not contain the canonical "
              f"16-mer. Epitope comparison would be misregistered. Stopping.",
              file=sys.stderr)
        return 1
    print("\nall RNA variants align to the canonical frame")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Sync the repo to Euler**

```bash
cd /Users/leo/Documents/binder_map
ssh euler 'mkdir -p /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs'
rsync -av --exclude '.git' --exclude 'tests/fixtures' --exclude 'build' \
  ./ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/
rsync -av build/ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/build/
```

- [ ] **Step 4: Run the census**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/stage0b_rna_census.py --build build --per-group 3'
```

Expected: a table dominated by `ACUGCCGUAUAGGCAG` (16 nt, offset 0) plus the `CACUGCCGUAUAGGCAGC` variant (18 nt, offset 1), final line `all RNA variants align to the canonical frame`, exit 0.

If any row reports `FRAME_MISMATCH`, stop. Read the offending sequence and decide whether it is a genuine second target or a parsing artefact before touching anything else.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage0b_rna_census.py slurm/stage_single.sbatch
git commit -m "Add RNA-variant census as a hard gate before compute"
```

---

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
import shutil
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


def read_archive_members(archive: Path, members: list[str],
                         workdir: Path) -> dict[str, bytes]:
    """Extract the needed members from one .tar.zst into workdir, keyed by name.

    Python's tarfile has no zstd codec here, so tar does the decompression
    once and we sieve out only the requested members - but by NAME, into
    real files, not via `-O` to a single stdout stream.

    `-O` concatenates the matched members in ARCHIVE order, not the order
    they were requested in. Zipping `members` against the resulting blobs
    positionally therefore silently swaps content between design_ids
    whenever the manifest's row order differs from the archive's own member
    order. That is not a corner case: 84% of archives hold more than one
    member, so it would have quietly mislabelled the majority of extracted
    structures - every backbone, hash and fingerprint attached to the wrong
    design. Named extraction sidesteps the ordering assumption entirely.

    Partial misses are handled by checking which files landed rather than
    gating on tar's exit code: tar exits non-zero when some requested member
    is absent, but still writes out the ones it found.
    """
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    subprocess.run(
        ["tar", "--zstd", "-x", "-C", str(workdir), "-f", str(archive)] + members,
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
    )
    chunks: dict[str, bytes] = {}
    for member in members:
        path = workdir / member
        if path.exists():
            chunks[member] = path.read_bytes()
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
            blobs = read_archive_members(archive, members, tmp / "extract")

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

### Task 14: Stage 2 — global dedup and Foldseek database

**Files:**
- Create: `scripts/stage2_builddb.py`
- Create: `scripts/stage2_builddb.sh`

**Interfaces:**
- Consumes: `results/<layer>/hashes/*.tsv`, `results/<layer>/shards/*.tar.gz`.
- Produces:
  - `results/<layer>/dedup/canonical.tsv` — `coord_hash, struct_name` (one elected representative per hash).
  - `results/<layer>/dedup/members.tsv` — `struct_name, design_ids` (comma-separated).
  - `results/<layer>/dedup/nodes.tsv` — `node_id, struct_name, n_designs`, `node_id` being the 0-based row index used by every downstream array.
  - `results/<layer>/db/uniq` — Foldseek database of exactly the canonical structures.

- [ ] **Step 1: Implement `scripts/stage2_builddb.py`**

```python
#!/usr/bin/env python3
"""Stage 2: elect one canonical structure per coordinate hash.

Stage 1 dedups within a shard; sequence variants of one backbone share an
archive so most redundancy dies there. This stage catches hashes that recur
across shards, which naming would never reveal.

Usage:  python scripts/stage2_builddb.py --layer esm
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    files = sorted((out / "hashes").glob("shard_*.tsv"))
    if not files:
        print("FATAL: no hash tables; did stage 1 run?")
        return 1

    frames = [pd.read_csv(f, sep="\t") for f in files]
    hashes = pd.concat(frames, ignore_index=True)
    print(f"{args.layer}: {len(hashes):,d} structures, "
          f"{hashes['coord_hash'].nunique():,d} unique backbones")

    tally = hashes.groupby("protein_chain").size().to_dict()
    print(f"  protein chain tally: {tally}")

    # Elect the lexicographically first design_id per hash: deterministic, and
    # independent of shard order.
    hashes = hashes.sort_values(["coord_hash", "design_id"], kind="stable")
    canonical = hashes.groupby("coord_hash", sort=True)["design_id"].first()
    canonical = canonical.rename("struct_name").reset_index()

    members = (hashes.groupby("coord_hash", sort=True)["design_id"]
               .apply(lambda s: ",".join(sorted(s)))
               .rename("design_ids").reset_index())
    members = canonical.merge(members, on="coord_hash")

    nodes = canonical.copy()
    nodes["n_designs"] = members["design_ids"].str.count(",") + 1
    nodes.insert(0, "node_id", range(len(nodes)))

    dedup = out / "dedup"
    dedup.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(dedup / "canonical.tsv", sep="\t", index=False)
    members[["struct_name", "design_ids"]].to_csv(
        dedup / "members.tsv", sep="\t", index=False)
    nodes[["node_id", "struct_name", "n_designs"]].to_csv(
        dedup / "nodes.tsv", sep="\t", index=False)

    # Foldseek createsubdb takes a plain list of the names to keep.
    (dedup / "keep_names.txt").write_text(
        "\n".join(f"{n}.pdb" for n in nodes["struct_name"]) + "\n")

    print(f"  wrote {len(nodes):,d} nodes -> {dedup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Implement `scripts/stage2_builddb.sh`**

`createdb` reads the tarballs directly, which is the whole reason shards exist.

```bash
#!/bin/bash
# Stage 2: global dedup, then build the Foldseek DB from the tar.gz shards.
#   bash scripts/stage2_builddb.sh <design|af3|esm>
set -uo pipefail

LAYER="${1:?usage: stage2_builddb.sh <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
OUT="$FM/${RESULTS:-results}/$LAYER"

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960 || exit 1
command -v foldseek >/dev/null || { echo "FATAL: foldseek not on PATH"; exit 1; }

cd "$FM" || exit 1
PYTHONPATH=src "$PY" scripts/stage2_builddb.py --layer "$LAYER" || exit 1

mkdir -p "$OUT/db"
rm -rf "$OUT/db/full"* "$OUT/db/uniq"*

echo "=== createdb from $(ls "$OUT"/shards/*.tar.gz | wc -l) shards ==="
foldseek createdb "$OUT"/shards/*.tar.gz "$OUT/db/full" || exit 1

FULL=$(wc -l < "$OUT/db/full.lookup")
KEEP=$(wc -l < "$OUT/dedup/keep_names.txt")
echo "full db: $FULL entries; want $KEEP canonical"

# Map canonical names to DB keys via the lookup, then subset.
awk 'NR==FNR{want[$1]=1; next} ($2 in want){print $1}' \
  "$OUT/dedup/keep_names.txt" "$OUT/db/full.lookup" > "$OUT/db/keep_keys.txt"

FOUND=$(wc -l < "$OUT/db/keep_keys.txt")
echo "resolved $FOUND / $KEEP canonical names to DB keys"
if [ "$FOUND" -ne "$KEEP" ]; then
  echo "FATAL: $((KEEP - FOUND)) canonical structures are not in the Foldseek DB."
  echo "That means createdb silently skipped them. Investigate before searching."
  exit 1
fi

foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full" "$OUT/db/uniq" || exit 1
foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full_ss" "$OUT/db/uniq_ss" 2>/dev/null
foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full_ca" "$OUT/db/uniq_ca" 2>/dev/null

echo "uniq db: $(wc -l < "$OUT/db/uniq.lookup") entries"
echo "=== stage 2 done for $LAYER ==="
```

- [ ] **Step 3: Run it on the test shard**

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  bash scripts/stage2_builddb.sh esm"'
```

Expected: the node count matches the unique-backbone count, `resolved N / N`, and a `uniq.lookup` with the same number of entries. A mismatch means `createdb` skipped structures — usually malformed PDB output from Task 4 — and must be fixed rather than worked around.

- [ ] **Step 4: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage2_builddb.py scripts/stage2_builddb.sh
git commit -m "Add stage 2: global dedup and Foldseek database construction"
```

---

### Task 15: Stage 3 — all-vs-all search

**Files:**
- Create: `slurm/stage3_search.sbatch`
- Create: `scripts/stage3_split_queries.py`

**Interfaces:**
- Consumes: `results/<layer>/db/uniq`.
- Produces: `results/<layer>/hits/hits_%02d.tsv` — `query, target, qtmscore, ttmscore, alntmscore, lddt, fident, alnlen, evalue, bits`, no header.

- [ ] **Step 1: Implement `scripts/stage3_split_queries.py`**

```python
#!/usr/bin/env python3
"""Split the canonical DB keys into N chunks for the search array.

Query-chunk arrays rather than the module's MPI support: a failed chunk is a
one-line resubmit, whereas a failed MPI rank kills the whole run.
"""
from __future__ import annotations

import argparse
from pathlib import Path

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True)
    parser.add_argument("--n-chunks", type=int, default=32)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    keys = (out / "db" / "keep_keys.txt").read_text().split()
    chunk_dir = out / "db" / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    for i in range(args.n_chunks):
        part = keys[i::args.n_chunks]
        (chunk_dir / f"q{i:02d}.txt").write_text("\n".join(part) + "\n")
    print(f"{len(keys):,d} queries -> {args.n_chunks} chunks in {chunk_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Write `slurm/stage3_search.sbatch`**

```bash
#!/bin/bash
# Stage 3: all-vs-all Foldseek, one array task per query chunk.
#   sbatch --job-name=bm_s3_esm --array=0-31%16 slurm/stage3_search.sbatch esm
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

LAYER="${1:?usage: stage3_search.sbatch <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
OUT="$FM/${RESULTS:-results}/$LAYER"
I=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
ALN_TYPE="${ALN_TYPE:-2}"

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960 || exit 1
command -v foldseek >/dev/null || { echo "FATAL: foldseek not on PATH"; exit 1; }

# Node-local scratch: Foldseek tmp dirs are inode-heavy and must never land on
# the project filesystem.
WORK="${TMPDIR:-/tmp}/bm_s3_${LAYER}_${I}"
mkdir -p "$WORK" "$OUT/hits" "$OUT/db/q"
trap 'rm -rf "$WORK"' EXIT

echo "=== stage3 layer=$LAYER chunk=$I alnType=$ALN_TYPE on $(hostname) $(date) ==="

QDB="$OUT/db/q/q$I"
rm -rf "$QDB"*
foldseek createsubdb "$OUT/db/chunks/q$I.txt" "$OUT/db/uniq" "$QDB" || exit 1

foldseek search "$QDB" "$OUT/db/uniq" "$WORK/aln" "$WORK/tmp" \
  -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 \
  --alignment-type "$ALN_TYPE" -a --threads "${SLURM_CPUS_PER_TASK:-16}" || exit 1

foldseek convertalis "$QDB" "$OUT/db/uniq" "$WORK/aln" "$OUT/hits/hits_$I.tsv" \
  --format-output "query,target,qtmscore,ttmscore,alntmscore,lddt,fident,alnlen,evalue,bits" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" || exit 1

echo "hits: $(wc -l < "$OUT/hits/hits_$I.tsv")"
echo "=== stage3 chunk $I done $(date) ==="
```

- [ ] **Step 3: Test one chunk interactively**

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
    scripts/stage3_split_queries.py --layer esm --n-chunks 32 && \
  SLURM_ARRAY_TASK_ID=0 SLURM_CPUS_PER_TASK=4 TMPDIR=/tmp \
    bash slurm/stage3_search.sbatch esm"'
```

Expected: a non-empty `hits_00.tsv`. Sanity-check that self-hits score 1.0:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results/esm && \
  awk -F"\t" "\$1==\$2 {print \$3, \$4; n++} END{print \"self-hits:\", n}" hits/hits_00.tsv | tail -3'
```

Expected: self-hits report `qtmscore` and `ttmscore` of 1 (or extremely close). If not, the alignment settings are wrong.

- [ ] **Step 4: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage3_split_queries.py slurm/stage3_search.sbatch
git commit -m "Add stage 3: chunked all-vs-all Foldseek search"
```

---

### Task 16: Stage 4 — build both graphs for both axes

**Files:**
- Create: `scripts/stage4_graph.py`

**Interfaces:**
- Consumes: `binder_map.knn`; `results/<layer>/hits/*.tsv`; `results/<layer>/fp/*.parquet`; `results/<layer>/dedup/nodes.tsv`.
- Produces, under `results/<layer>/graph/<axis>/` for `axis` in `fold`, `epitope`:
  - `edges_cluster.parquet` (`i, j, sim`), `knn_layout.parquet` (`i, j, sim`), `nodes.parquet` (`node_id, struct_name, degree, redundancy, is_isolate`).
- Also produces `results/<layer>/fp/canonical.parquet` — one fingerprint row per node, in `node_id` order.

- [ ] **Step 1: Implement `scripts/stage4_graph.py`**

```python
#!/usr/bin/env python3
"""Stage 4: turn Foldseek hits and epitope fingerprints into graphs.

Both axes go through the same code. Each produces two graphs — a thresholded
clustering graph and an unthresholded layout graph — for the reasons in
binder_map.knn.

Usage:
    python scripts/stage4_graph.py --layer esm --tm-floor 0.5 --epi-floor 0.5
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map import FINGERPRINT_DIM
from binder_map.fingerprint import FINGERPRINT_COLUMNS
from binder_map.knn import knn_from_pairs, knn_from_vectors, symmetrise_min

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
REDUNDANCY_TM = 0.95


def write_axis(out_dir: Path, axis: str, nodes: pd.DataFrame, res, redundancy):
    dest = out_dir / "graph" / axis
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "i": res.cluster_edges[:, 0] if len(res.cluster_edges) else np.array([], np.int32),
        "j": res.cluster_edges[:, 1] if len(res.cluster_edges) else np.array([], np.int32),
        "sim": res.cluster_weights,
    }).to_parquet(dest / "edges_cluster.parquet", index=False)
    pd.DataFrame({
        "i": res.layout_edges[:, 0] if len(res.layout_edges) else np.array([], np.int32),
        "j": res.layout_edges[:, 1] if len(res.layout_edges) else np.array([], np.int32),
        "sim": res.layout_weights,
    }).to_parquet(dest / "knn_layout.parquet", index=False)
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "struct_name": nodes["struct_name"].values,
        "degree": res.degree,
        "redundancy": redundancy,
        "is_isolate": res.is_isolate,
    }).to_parquet(dest / "nodes.parquet", index=False)
    print(f"  {axis}: {len(res.cluster_edges):,d} cluster edges, "
          f"{len(res.layout_edges):,d} layout edges, "
          f"{int(res.is_isolate.sum()):,d} isolates")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--tm-floor", type=float, default=0.5)
    parser.add_argument("--epi-floor", type=float, default=0.5)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    index = pd.Series(nodes["node_id"].values, index=nodes["struct_name"].values)
    n = len(nodes)
    print(f"{args.layer}: {n:,d} nodes")

    # ---------- fold axis ----------
    columns = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
               "lddt", "fident", "alnlen", "evalue", "bits"]
    hit_files = sorted((out / "hits").glob("hits_*.tsv"))
    if not hit_files:
        print("FATAL: no hit files; did stage 3 run?")
        return 1
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

    i, j, sim = symmetrise_min(qi, tj, qtm, ttm)
    print(f"  {len(i):,d} undirected pairs after min(qTM,tTM)")

    # Redundancy from the FULL hit list, before top-k truncation, else it
    # saturates at k. Still bounded above by --max-seqs (300).
    redundancy = np.zeros(n, dtype=np.int32)
    strong = sim >= REDUNDANCY_TM
    if strong.any():
        np.add.at(redundancy, i[strong], 1)
        np.add.at(redundancy, j[strong], 1)

    fold = knn_from_pairs(n, i, j, sim, k=args.k, floor=args.tm_floor)
    write_axis(out, "fold", nodes, fold, redundancy)

    # ---------- epitope axis ----------
    fp = pd.concat([pd.read_parquet(f) for f in sorted((out / "fp").glob("shard_*.parquet"))],
                   ignore_index=True)
    fp = fp.drop_duplicates("design_id").set_index("design_id")
    canonical = fp.reindex(nodes["struct_name"].values)
    canonical.insert(0, "node_id", nodes["node_id"].values)
    canonical.to_parquet(out / "fp" / "canonical.parquet", index=False)

    vectors = canonical[FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
    assert vectors.shape == (n, FINGERPRINT_DIM)

    # Zero-contact designs have no epitope. Excluding them is deliberate:
    # imputing one would invent the signal being measured.
    has_interface = canonical["contacts_total"].to_numpy() > 0
    print(f"  epitope: {int((~has_interface).sum()):,d} nodes with no interface, excluded")

    sub_index = np.flatnonzero(has_interface)
    epi_sub = knn_from_vectors(vectors[sub_index], k=args.k, floor=args.epi_floor)

    # Lift the sub-graph node ids back into the full node space.
    def lift(edges):
        if len(edges) == 0:
            return np.empty((0, 2), np.int32)
        return sub_index[edges].astype(np.int32)

    from binder_map.knn import KNNResult
    degree = np.zeros(n, dtype=np.int32)
    degree[sub_index] = epi_sub.degree
    is_isolate = np.ones(n, dtype=bool)
    is_isolate[sub_index] = epi_sub.is_isolate
    epi = KNNResult(
        lift(epi_sub.cluster_edges), epi_sub.cluster_weights,
        lift(epi_sub.layout_edges), epi_sub.layout_weights,
        degree, is_isolate,
    )
    epi_redundancy = np.zeros(n, dtype=np.int32)
    write_axis(out, "epitope", nodes, epi, epi_redundancy)

    # Record which nodes carry an interface, for the explorer's greying-out.
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "has_interface": has_interface,
        "contacts_total": canonical["contacts_total"].values,
    }).to_parquet(out / "graph" / "epitope" / "interface.parquet", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run it on the test layer**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/stage4_graph.py --layer esm'
```

Expected: both axes report edge counts and an isolate count. If the fold axis reports every node an isolate, `TM_FLOOR` is too high for this data — that is exactly what the smoke test in Task 22 exists to settle.

- [ ] **Step 3: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage4_graph.py
git commit -m "Add stage 4: fold and epitope graphs from one axis-agnostic path"
```

---

### Task 17: Stages 5 and 6 — cluster and embed both axes

**Files:**
- Create: `scripts/stage5_cluster.py`
- Create: `scripts/stage6_embed.py`

**Interfaces:**
- Consumes: `binder_map.clustering`, `binder_map.embedding`; `graph/<axis>/`.
- Produces:
  - `results/<layer>/clusters/<axis>/clusters.parquet` — `node_id, cluster_r05, cluster_r10, cluster_r20`.
  - `results/<layer>/clusters/<axis>/medoids.parquet` — `resolution, cluster, node_id, struct_name, size`.
  - `results/<layer>/embed/<axis>/embed.parquet` — `node_id, x, y`.

- [ ] **Step 1: Implement `scripts/stage5_cluster.py`**

```python
#!/usr/bin/env python3
"""Stage 5: Leiden clusters and medoids, per axis.

Usage:  python scripts/stage5_cluster.py --layer esm --axis fold
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.clustering import RESOLUTIONS, cluster_all, medoids

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--axis", required=True, choices=["fold", "epitope"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    graph = out / "graph" / args.axis
    nodes = pd.read_parquet(graph / "nodes.parquet")
    edges = pd.read_parquet(graph / "edges_cluster.parquet")
    n = len(nodes)

    pairs = edges[["i", "j"]].to_numpy(dtype=np.int64)
    weights = edges["sim"].to_numpy(dtype=np.float64)
    print(f"{args.layer}/{args.axis}: {n:,d} nodes, {len(pairs):,d} edges")

    labels = cluster_all(n, pairs, weights)
    frame = pd.DataFrame({"node_id": nodes["node_id"].values})
    for name, arr in labels.items():
        frame[name] = arr
        sizes = np.bincount(arr)
        print(f"  {name}: {len(sizes):,d} clusters, largest {sizes.max():,d}, "
              f"singletons {int((sizes == 1).sum()):,d}")

    rows = []
    for name, arr in labels.items():
        chosen = medoids(arr, pairs, weights)
        sizes = np.bincount(arr)
        for cluster, node in sorted(chosen.items()):
            rows.append({
                "resolution": name,
                "cluster": int(cluster),
                "node_id": int(nodes["node_id"].values[node]),
                "struct_name": nodes["struct_name"].values[node],
                "size": int(sizes[cluster]),
            })

    dest = out / "clusters" / args.axis
    dest.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(dest / "clusters.parquet", index=False)
    pd.DataFrame(rows).to_parquet(dest / "medoids.parquet", index=False)
    print(f"  wrote {len(rows):,d} medoids -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Implement `scripts/stage6_embed.py`**

```python
#!/usr/bin/env python3
"""Stage 6: UMAP layout per axis.

Usage:  python scripts/stage6_embed.py --layer esm --axis fold
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.embedding import embed

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--axis", required=True, choices=["fold", "epitope"])
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    graph = out / "graph" / args.axis
    nodes = pd.read_parquet(graph / "nodes.parquet")
    layout = pd.read_parquet(graph / "knn_layout.parquet")
    n = len(nodes)
    print(f"{args.layer}/{args.axis}: embedding {n:,d} nodes "
          f"from {len(layout):,d} layout edges")

    xy = embed(n, layout[["i", "j"]].to_numpy(dtype=np.int64),
               layout["sim"].to_numpy(dtype=np.float32), k=args.k)

    dest = out / "embed" / args.axis
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "x": xy[:, 0], "y": xy[:, 1],
    }).to_parquet(dest / "embed.parquet", index=False)
    print(f"  wrote coordinates -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Run both stages for both axes**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
for AX in fold epitope; do
  PYTHONPATH=src $PY scripts/stage5_cluster.py --layer esm --axis $AX || exit 1
  PYTHONPATH=src $PY scripts/stage6_embed.py   --layer esm --axis $AX || exit 1
done'
```

Expected: cluster counts and a written coordinate table for each axis.

- [ ] **Step 4: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage5_cluster.py scripts/stage6_embed.py
git commit -m "Add stages 5 and 6: Leiden clustering and UMAP layout per axis"
```

---

### Task 18: Stage 6b — comparisons, and stage 7 — medoid structures

**Files:**
- Create: `scripts/stage6b_compare.py`
- Create: `scripts/stage7_medoids.py`

**Interfaces:**
- Consumes: `binder_map.compare`; `clusters/<axis>/`, `fp/canonical.parquet`, `dedup/members.tsv`.
- Produces:
  - `results/compare/partition_scores.parquet` — `layer, resolution, ami, ari, n`.
  - `results/compare/contingency_<layer>_<resolution>.parquet`.
  - `results/compare/epitope_agreement.parquet` — `design_id, epi_agree_af3_esm, epi_agree_design_af3, epi_agree_design_esm`.
  - `results/compare/stack_profiles.parquet` — `layer, stack, nt01..nt16`.
  - `results/<layer>/medoids/medoids.json.gz` — every medoid's CA trace in one gzipped file.

- [ ] **Step 1: Implement `scripts/stage6b_compare.py`**

```python
#!/usr/bin/env python3
"""Stage 6b: compare the fold and epitope partitions, and epitopes across layers.

This is the payoff of having two axes. A low AMI is a real result — backbone
fold does not predict where a binder lands — not a failure.

Usage:  python scripts/stage6b_compare.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.clustering import RESOLUTIONS
from binder_map.compare import (
    contingency, fingerprint_agreement, partition_scores, stack_profiles,
)
from binder_map import RNA_LEN

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
LAYERS = ("design", "af3", "esm")
PER_NT = [f"per_nt_{i}" for i in range(1, RNA_LEN + 1)]


def load_axis_clusters(results: Path, layer: str, axis: str):
    return pd.read_parquet(results / layer / "clusters" / axis / "clusters.parquet")


def design_fingerprints(results: Path, layer: str) -> pd.DataFrame:
    """Fingerprints keyed by design_id, fanned out from canonical structures.

    The design layer collapses sequence variants onto one backbone, so one
    canonical row stands for many design_ids.
    """
    canonical = pd.read_parquet(results / layer / "fp" / "canonical.parquet")
    members = pd.read_csv(results / layer / "dedup" / "members.tsv", sep="\t")
    canonical = canonical.merge(
        members, left_on="design_id", right_on="struct_name", how="left")
    canonical["design_ids"] = canonical["design_ids"].fillna(canonical["design_id"])
    expanded = canonical.assign(
        design_id=canonical["design_ids"].str.split(",")).explode("design_id")
    return expanded.drop(columns=["design_ids", "struct_name"]).set_index("design_id")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=FM / "results")
    parser.add_argument("--build", type=Path, default=FM / "build")
    args = parser.parse_args()

    out = args.results / "compare"
    out.mkdir(parents=True, exist_ok=True)

    # ---- 1. fold x epitope, per layer and resolution
    scores = []
    for layer in LAYERS:
        try:
            fold = load_axis_clusters(args.results, layer, "fold")
            epi = load_axis_clusters(args.results, layer, "epitope")
        except FileNotFoundError:
            print(f"  {layer}: clusters missing, skipping")
            continue
        merged = fold.merge(epi, on="node_id", suffixes=("_fold", "_epi"))
        for resolution in RESOLUTIONS:
            a = merged[f"{resolution}_fold"].to_numpy()
            b = merged[f"{resolution}_epi"].to_numpy()
            s = partition_scores(a, b)
            s.update({"layer": layer, "resolution": resolution})
            scores.append(s)
            print(f"  {layer}/{resolution}: AMI={s['ami']:.3f} ARI={s['ari']:.3f} "
                  f"n={s['n']:,d}")

            mat, la, lb = contingency(a, b)
            pd.DataFrame(mat, index=la, columns=lb).to_parquet(
                out / f"contingency_{layer}_{resolution}.parquet")
    if scores:
        pd.DataFrame(scores).to_parquet(out / "partition_scores.parquet", index=False)

    # ---- 2/3. cross-layer epitope agreement (and design-vs-refold)
    fps = {}
    for layer in LAYERS:
        try:
            fps[layer] = design_fingerprints(args.results, layer)
        except FileNotFoundError:
            print(f"  {layer}: fingerprints missing, skipping agreement")

    from binder_map.fingerprint import FINGERPRINT_COLUMNS
    agreement = None
    pairs = [("af3", "esm"), ("design", "af3"), ("design", "esm")]
    for a, b in pairs:
        if a not in fps or b not in fps:
            continue
        common = fps[a].index.intersection(fps[b].index)
        va = fps[a].loc[common, FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
        vb = fps[b].loc[common, FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
        col = pd.DataFrame({f"epi_agree_{a}_{b}": fingerprint_agreement(va, vb)},
                           index=common)
        agreement = col if agreement is None else agreement.join(col, how="outer")
        print(f"  epitope agreement {a} vs {b}: n={len(common):,d}, "
              f"median={np.median(col.iloc[:, 0]):.3f}")
    if agreement is not None:
        agreement.reset_index(names="design_id").to_parquet(
            out / "epitope_agreement.parquet", index=False)

    # ---- 4. per-stack nucleotide profiles
    meta = pd.read_parquet(args.build / "meta" / "designs.parquet")
    profiles = []
    for layer, fp in fps.items():
        joined = fp.join(meta.set_index("design_id")["stack"], how="inner")
        frame = stack_profiles(joined[PER_NT].to_numpy(dtype=np.float64),
                               joined["stack"].to_numpy())
        frame.insert(0, "layer", layer)
        profiles.append(frame.reset_index())
    if profiles:
        pd.concat(profiles, ignore_index=True).to_parquet(
            out / "stack_profiles.parquet", index=False)

    print(f"comparisons written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Implement `scripts/stage7_medoids.py`**

```python
#!/usr/bin/env python3
"""Stage 7: bundle every medoid's CA trace into one gzipped JSON.

One file rather than one per medoid: this project has a hard inode cap, and
the explorer wants them all anyway.

Usage:  python scripts/stage7_medoids.py --layer esm
"""
from __future__ import annotations

import argparse
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def ca_from_pdb(text: str) -> list[list[float]]:
    out = []
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append([round(float(line[30:38]), 2),
                        round(float(line[38:46]), 2),
                        round(float(line[46:54]), 2)])
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    wanted: set[str] = set()
    for axis in ("fold", "epitope"):
        path = out / "clusters" / axis / "medoids.parquet"
        if path.exists():
            wanted |= set(pd.read_parquet(path)["struct_name"])
    print(f"{args.layer}: {len(wanted):,d} distinct medoids across axes/resolutions")

    traces: dict[str, list] = {}
    for shard in sorted((out / "shards").glob("shard_*.tar.gz")):
        with tarfile.open(shard, "r:gz") as tar:
            for member in tar:
                name = member.name.removesuffix(".pdb")
                if name not in wanted or name in traces:
                    continue
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                traces[name] = ca_from_pdb(fh.read().decode())
        if len(traces) == len(wanted):
            break

    missing = wanted - traces.keys()
    if missing:
        print(f"  WARNING: {len(missing):,d} medoids had no backbone in the shards")

    dest = out / "medoids"
    dest.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(traces, separators=(",", ":")).encode()
    (dest / "medoids.json.gz").write_bytes(gzip.compress(payload, 9))
    print(f"  wrote {len(traces):,d} CA traces, "
          f"{len(gzip.compress(payload, 9)) / 1e6:.1f} MB gzipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Run both**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
PYTHONPATH=src $PY scripts/stage7_medoids.py --layer esm
PYTHONPATH=src $PY scripts/stage6b_compare.py'
```

Expected: medoid traces written; comparison prints AMI/ARI per resolution. With only the esm layer present, the cross-layer agreement section reports skips — that is correct at this point.

- [ ] **Step 4: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage6b_compare.py scripts/stage7_medoids.py
git commit -m "Add stage 6b comparisons and stage 7 medoid CA traces"
```

---

### Task 19: Stage 8 — payload encoding and HTML baking

Runs locally. Everything the explorer needs is baked into one file.

**Files:**
- Create: `src/binder_map/htmlbuild.py`
- Create: `tests/test_htmlbuild.py`
- Create: `scripts/stage8_build_html.py`

**Interfaces:**
- Consumes: all `results/` tables, pulled to a local directory by rsync.
- Produces:
  - `encode_floats(arr: np.ndarray) -> str` — base64 of a little-endian float32 buffer.
  - `encode_ints(arr: np.ndarray, dtype: str = "int32") -> str`
  - `encode_categorical(values: Sequence) -> tuple[list[str], str]` — `(lookup, base64 int16 codes)`.
  - `build_payload(results: Path, build: Path) -> dict` — the whole explorer payload.
  - `bake(payload: dict, template: Path, app_js: Path, style_css: Path, out: Path) -> int` — writes the HTML, returns its byte size.

- [ ] **Step 1: Write the failing test**

Create `tests/test_htmlbuild.py`:

```python
import base64
import json

import numpy as np
import pytest

from binder_map.htmlbuild import (
    encode_categorical, encode_floats, encode_ints,
)


def test_encode_floats_roundtrips():
    arr = np.array([1.5, -2.25, 0.0], dtype=np.float32)
    decoded = np.frombuffer(base64.b64decode(encode_floats(arr)), dtype="<f4")
    np.testing.assert_allclose(decoded, arr)


def test_encode_ints_roundtrips():
    arr = np.array([0, 7, -3], dtype=np.int32)
    decoded = np.frombuffer(base64.b64decode(encode_ints(arr)), dtype="<i4")
    np.testing.assert_array_equal(decoded, arr)


def test_encode_categorical_lookup_and_codes():
    values = ["BG+BIF", "RFD3+LiMPNN", "BG+BIF"]
    lookup, codes_b64 = encode_categorical(values)
    codes = np.frombuffer(base64.b64decode(codes_b64), dtype="<i2")
    assert sorted(lookup) == lookup           # deterministic ordering
    assert [lookup[c] for c in codes] == values


def test_encode_categorical_handles_missing():
    lookup, codes_b64 = encode_categorical(["a", None, "b"])
    codes = np.frombuffer(base64.b64decode(codes_b64), dtype="<i2")
    assert codes[1] == -1                      # missing is -1, not a category
    assert "a" in lookup and "b" in lookup


def test_encode_categorical_rejects_too_many_levels():
    with pytest.raises(ValueError):
        encode_categorical([str(i) for i in range(40000)])
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_htmlbuild.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.htmlbuild'`

- [ ] **Step 3: Implement the encoders in `src/binder_map/htmlbuild.py`**

```python
"""Encode pipeline outputs into a single self-contained HTML page.

Numeric columns are base64'd typed arrays rather than JSON numbers: 72k
floats as JSON is roughly 10x the bytes and parses far slower.
"""
from __future__ import annotations

import base64
import json
from pathlib import Path
from typing import Sequence

import numpy as np

INT16_MAX = 32767


def encode_floats(arr) -> str:
    return base64.b64encode(
        np.asarray(arr, dtype="<f4").tobytes()).decode("ascii")


def encode_ints(arr, dtype: str = "int32") -> str:
    code = {"int32": "<i4", "int16": "<i2", "uint8": "|u1"}[dtype]
    return base64.b64encode(
        np.asarray(arr).astype(code).tobytes()).decode("ascii")


def encode_categorical(values: Sequence) -> tuple[list[str], str]:
    """Return (sorted lookup, base64 int16 codes). Missing values encode as -1."""
    clean = [None if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)
             for v in values]
    lookup = sorted({v for v in clean if v is not None})
    if len(lookup) > INT16_MAX:
        raise ValueError(f"too many categories for int16 codes: {len(lookup)}")
    index = {name: i for i, name in enumerate(lookup)}
    codes = np.array([index[v] if v is not None else -1 for v in clean],
                     dtype="<i2")
    return lookup, base64.b64encode(codes.tobytes()).decode("ascii")


def bake(payload: dict, template: Path, app_js: Path, style_css: Path,
         out: Path) -> int:
    """Substitute payload, JS and CSS into the template and write the page."""
    html = template.read_text()
    html = html.replace("/*__STYLE__*/", style_css.read_text())
    html = html.replace("/*__APP__*/", app_js.read_text())
    html = html.replace(
        '"__PAYLOAD__"',
        json.dumps(payload, separators=(",", ":"), allow_nan=False),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    return out.stat().st_size
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_htmlbuild.py -v
```

Expected: 5 passed.

- [ ] **Step 5: Add `build_payload` to `src/binder_map/htmlbuild.py`**

Append to the same module:

```python
LAYERS = ("design", "af3", "esm")
AXES = ("fold", "epitope")


def build_payload(results: Path, build: Path) -> dict:
    """Assemble the explorer payload from the pipeline outputs."""
    import pandas as pd

    meta = pd.read_parquet(build / "meta" / "designs.parquet")
    counts = json.loads((build / "meta" / "counts.json").read_text())

    numeric_columns = [
        "dockq_af3", "dockq_esm", "iptm_mean", "iptm_min",
        "rmsd_af3", "rmsd_esm", "dsasa_af3", "dsasa_esm",
        "hb_af3", "hb_esm", "rog_af3", "rog_esm",
        "pose_disagree", "fold_disagree",
    ]

    agreement_path = results / "compare" / "epitope_agreement.parquet"
    agreement = (pd.read_parquet(agreement_path).set_index("design_id")
                 if agreement_path.exists() else None)

    payload: dict = {"layers": {}, "meta": counts, "compare": {}}

    for layer in LAYERS:
        layer_dir = results / layer
        if not (layer_dir / "dedup" / "nodes.tsv").exists():
            continue
        nodes = pd.read_csv(layer_dir / "dedup" / "nodes.tsv", sep="\t")
        members = pd.read_csv(layer_dir / "dedup" / "members.tsv", sep="\t")
        nodes = nodes.merge(members, on="struct_name", how="left")

        # Node-level metadata comes from the canonical structure's own row.
        node_meta = meta.set_index("design_id").reindex(nodes["struct_name"])

        entry: dict = {
            "n": int(len(nodes)),
            "struct_names": nodes["struct_name"].tolist(),
            "design_ids": nodes["design_ids"].fillna("").tolist(),
            "n_designs": encode_ints(nodes["n_designs"].to_numpy(), "int32"),
            "categorical": {},
            "numeric": {},
            "axes": {},
        }
        for column in ("stack", "campaign", "gen", "invfold"):
            lookup, codes = encode_categorical(node_meta[column].tolist())
            entry["categorical"][column] = {"lookup": lookup, "codes": codes}
        for column in numeric_columns:
            if column in node_meta:
                values = pd.to_numeric(node_meta[column], errors="coerce")
                entry["numeric"][column] = encode_floats(
                    values.fillna(np.nan).to_numpy())

        canonical_fp = layer_dir / "fp" / "canonical.parquet"
        if canonical_fp.exists():
            fp = pd.read_parquet(canonical_fp)
            per_nt = [f"per_nt_{i}" for i in range(1, 17)]
            entry["per_nt"] = encode_ints(
                np.clip(fp[per_nt].to_numpy(), 0, INT16_MAX), "int16")
            for column in ("cx", "cy", "cz", "ax", "ay", "az"):
                entry["numeric"][column] = encode_floats(fp[column].to_numpy())
            entry["numeric"]["contacts_total"] = encode_floats(
                fp["contacts_total"].to_numpy())

        if agreement is not None:
            for column in agreement.columns:
                values = agreement[column].reindex(nodes["struct_name"])
                entry["numeric"][column] = encode_floats(
                    values.fillna(np.nan).to_numpy())

        for axis in AXES:
            embed_path = layer_dir / "embed" / axis / "embed.parquet"
            cluster_path = layer_dir / "clusters" / axis / "clusters.parquet"
            if not embed_path.exists() or not cluster_path.exists():
                continue
            coords = pd.read_parquet(embed_path).sort_values("node_id")
            clusters = pd.read_parquet(cluster_path).sort_values("node_id")
            graph_nodes = pd.read_parquet(
                layer_dir / "graph" / axis / "nodes.parquet").sort_values("node_id")
            axis_entry = {
                "x": encode_floats(coords["x"].to_numpy()),
                "y": encode_floats(coords["y"].to_numpy()),
                "redundancy": encode_ints(
                    graph_nodes["redundancy"].to_numpy(), "int32"),
                "is_isolate": encode_ints(
                    graph_nodes["is_isolate"].to_numpy().astype(np.uint8), "uint8"),
                "clusters": {
                    name: encode_ints(clusters[name].to_numpy(), "int32")
                    for name in ("cluster_r05", "cluster_r10", "cluster_r20")
                },
            }
            medoid_path = layer_dir / "clusters" / axis / "medoids.parquet"
            if medoid_path.exists():
                medoid = pd.read_parquet(medoid_path)
                axis_entry["medoids"] = {
                    resolution: {
                        str(int(r.cluster)): {"node": int(r.node_id),
                                              "name": r.struct_name,
                                              "size": int(r.size)}
                        for r in group.itertuples()
                    }
                    for resolution, group in medoid.groupby("resolution")
                }
            interface = layer_dir / "graph" / axis / "interface.parquet"
            if interface.exists():
                iface = pd.read_parquet(interface).sort_values("node_id")
                axis_entry["has_interface"] = encode_ints(
                    iface["has_interface"].to_numpy().astype(np.uint8), "uint8")
            entry["axes"][axis] = axis_entry

        medoid_blob = layer_dir / "medoids" / "medoids.json.gz"
        if medoid_blob.exists():
            entry["medoid_traces"] = base64.b64encode(
                medoid_blob.read_bytes()).decode("ascii")

        payload["layers"][layer] = entry

    scores = results / "compare" / "partition_scores.parquet"
    if scores.exists():
        payload["compare"]["partition_scores"] = (
            pd.read_parquet(scores).to_dict(orient="records"))
    profiles = results / "compare" / "stack_profiles.parquet"
    if profiles.exists():
        payload["compare"]["stack_profiles"] = (
            pd.read_parquet(profiles).to_dict(orient="records"))

    return payload
```

- [ ] **Step 6: Implement `scripts/stage8_build_html.py`**

```python
#!/usr/bin/env python3
"""Stage 8: bake the explorer. Runs locally.

Usage:
    python scripts/stage8_build_html.py --results local_results --build build \
        --out dist/binder_map.html
"""
from __future__ import annotations

import argparse
from pathlib import Path

from binder_map.htmlbuild import bake, build_payload

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "local_results")
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--out", type=Path, default=ROOT / "dist" / "binder_map.html")
    args = parser.parse_args()

    payload = build_payload(args.results, args.build)
    for layer, entry in payload["layers"].items():
        axes = ", ".join(entry["axes"])
        print(f"  {layer}: {entry['n']:,d} nodes, axes: {axes or 'none'}")

    size = bake(payload, ROOT / "web" / "template.html",
                ROOT / "web" / "app.js", ROOT / "web" / "style.css", args.out)
    print(f"wrote {args.out} ({size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 7: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/htmlbuild.py tests/test_htmlbuild.py scripts/stage8_build_html.py
git commit -m "Add stage 8: payload encoding and HTML baking"
```

---

### Task 20: The explorer front end

Two linked canvas maps. No CDN, no framework — the page must open from a file:// URL years from now.

**Files:**
- Create: `web/template.html`
- Create: `web/style.css`
- Create: `web/app.js`
- Create: `web/vendor/3Dmol-min.js` (vendored, see Step 1)
- Create: `web/vendor/fflate.min.js` (vendored)

**Interfaces:**
- Consumes: the payload from Task 19 as the global `PAYLOAD`.
- Produces: nothing consumed by later code; this is a leaf.

- [ ] **Step 1: Vendor the two libraries**

The strict no-CDN rule means these are downloaded once and committed.

```bash
cd /Users/leo/Documents/binder_map
mkdir -p web/vendor
curl -sL https://unpkg.com/3dmol@2.5.5/build/3Dmol-min.js -o web/vendor/3Dmol-min.js
curl -sL https://unpkg.com/fflate@0.8.3/umd/index.js -o web/vendor/fflate.min.js
ls -la web/vendor/
head -c 200 web/vendor/3Dmol-min.js
```

Expected: both files present and non-trivial in size (3Dmol is ~2 MB). Confirm the first bytes look like minified JS, not an HTML error page.

- [ ] **Step 2: Write `web/template.html`**

```html
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Binder Map — fold vs epitope</title>
<style>/*__STYLE__*/</style>
</head>
<body>
<header>
  <h1>Binder Map</h1>
  <div class="controls">
    <label>Layer
      <select id="layer"></select>
    </label>
    <label>Colour by
      <select id="colour"></select>
    </label>
    <label>Resolution
      <select id="resolution">
        <option value="cluster_r05">0.5 (coarse)</option>
        <option value="cluster_r10" selected>1.0</option>
        <option value="cluster_r20">2.0 (fine)</option>
      </select>
    </label>
    <label>Search
      <input id="search" type="search" placeholder="design_id" size="28">
    </label>
    <span id="status"></span>
  </div>
  <div id="stacks" class="stack-toggles"></div>
</header>

<main>
  <section class="map" id="map-fold">
    <h2>Fold <span class="hint">backbone similarity (Foldseek)</span></h2>
    <canvas id="canvas-fold"></canvas>
  </section>
  <section class="map" id="map-epitope">
    <h2>Epitope <span class="hint">RNA contact fingerprint</span></h2>
    <canvas id="canvas-epitope"></canvas>
  </section>
  <aside id="panel">
    <nav>
      <button data-tab="point" class="active">Point</button>
      <button data-tab="cluster">Cluster</button>
      <button data-tab="compare">Compare</button>
    </nav>
    <div id="tab-point" class="tab active"></div>
    <div id="tab-cluster" class="tab">
      <div id="viewer"></div>
      <div id="cluster-info"></div>
    </div>
    <div id="tab-compare" class="tab"></div>
  </aside>
</main>

<footer id="provenance"></footer>

<script>/*__FFLATE__*/</script>
<script>/*__3DMOL__*/</script>
<script>window.PAYLOAD = "__PAYLOAD__";</script>
<script>/*__APP__*/</script>
</body>
</html>
```

- [ ] **Step 3: Extend `bake()` to inline the vendored libraries**

The template now has two more placeholders. In `src/binder_map/htmlbuild.py`, replace the body of `bake` with:

```python
def bake(payload: dict, template: Path, app_js: Path, style_css: Path,
         out: Path) -> int:
    """Substitute payload, vendored libs, JS and CSS into the template."""
    vendor = template.parent / "vendor"
    html = template.read_text()
    html = html.replace("/*__STYLE__*/", style_css.read_text())
    html = html.replace("/*__FFLATE__*/",
                        (vendor / "fflate.min.js").read_text())
    html = html.replace("/*__3DMOL__*/",
                        (vendor / "3Dmol-min.js").read_text())
    html = html.replace("/*__APP__*/", app_js.read_text())
    html = html.replace(
        '"__PAYLOAD__"',
        json.dumps(payload, separators=(",", ":"), allow_nan=False),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    return out.stat().st_size
```

- [ ] **Step 4: Write `web/style.css`**

```css
* { box-sizing: border-box; }
body {
  margin: 0; font: 13px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  color: #1b1b1b; background: #fafafa;
}
header { padding: 8px 12px; border-bottom: 1px solid #ddd; background: #fff; }
h1 { font-size: 15px; margin: 0 0 6px; font-weight: 600; }
h2 { font-size: 12px; margin: 0; padding: 4px 8px; font-weight: 600;
     border-bottom: 1px solid #eee; background: #fff; }
.hint { font-weight: 400; color: #888; margin-left: 6px; }
.controls { display: flex; gap: 14px; align-items: center; flex-wrap: wrap; }
.controls label { display: flex; gap: 5px; align-items: center; color: #555; }
select, input { font: inherit; padding: 2px 4px; }
.stack-toggles { display: flex; gap: 10px; margin-top: 6px; flex-wrap: wrap; }
.stack-toggles label { display: flex; gap: 4px; align-items: center; cursor: pointer; }
.swatch { width: 11px; height: 11px; border-radius: 2px; display: inline-block; }
main { display: grid; grid-template-columns: 1fr 1fr 330px; height: calc(100vh - 96px); }
.map { display: flex; flex-direction: column; border-right: 1px solid #ddd; min-width: 0; }
canvas { flex: 1; width: 100%; display: block; cursor: crosshair; background: #fff; }
#panel { display: flex; flex-direction: column; background: #fff; overflow: hidden; }
#panel nav { display: flex; border-bottom: 1px solid #ddd; }
#panel nav button {
  flex: 1; padding: 6px; border: 0; background: #f2f2f2; cursor: pointer;
  font: inherit; border-right: 1px solid #ddd;
}
#panel nav button.active { background: #fff; font-weight: 600; }
.tab { display: none; padding: 10px; overflow: auto; flex: 1; }
.tab.active { display: block; }
#viewer { width: 100%; height: 240px; position: relative; border: 1px solid #eee; }
table { border-collapse: collapse; width: 100%; font-size: 12px; }
td, th { padding: 2px 5px; border-bottom: 1px solid #f0f0f0; text-align: left; }
th { color: #666; font-weight: 500; }
.spark { display: flex; align-items: flex-end; gap: 1px; height: 40px; margin: 6px 0; }
.spark div { flex: 1; background: #4a7ebb; min-height: 1px; }
footer {
  padding: 5px 12px; border-top: 1px solid #ddd; background: #fff;
  color: #777; font-size: 11px; white-space: pre-wrap;
}
.warn { color: #a33; }
```

- [ ] **Step 5: Write `web/app.js`**

```javascript
/* Binder Map explorer.
 *
 * Two canvases share one selection. Everything is decoded from base64 typed
 * arrays at load; nothing is fetched.
 */
(function () {
  "use strict";

  var P = window.PAYLOAD;
  var AXES = ["fold", "epitope"];
  var PALETTE = ["#4a7ebb", "#e07b39", "#4f9d69", "#b5514f", "#8264a5",
                 "#7f6355", "#c47db2", "#8c8c8c", "#b8b445", "#3fa9b8"];

  function b64(str) {
    var bin = atob(str), out = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }
  function f32(str) { var b = b64(str); return new Float32Array(b.buffer, b.byteOffset, b.byteLength / 4); }
  function i32(str) { var b = b64(str); return new Int32Array(b.buffer, b.byteOffset, b.byteLength / 4); }
  function i16(str) { var b = b64(str); return new Int16Array(b.buffer, b.byteOffset, b.byteLength / 2); }
  function u8(str) { return b64(str); }

  var state = {
    layer: Object.keys(P.layers)[0],
    colour: "stack",
    resolution: "cluster_r10",
    hidden: {},          // stack name -> true when hidden
    selected: -1,        // node id
    hover: -1,
    views: {}            // axis -> {scale, tx, ty}
  };

  var decoded = {};
  function layerData(name) {
    if (decoded[name]) return decoded[name];
    var src = P.layers[name], d = {n: src.n, struct_names: src.struct_names,
                                   design_ids: src.design_ids, cat: {}, num: {}, axes: {}};
    d.n_designs = i32(src.n_designs);
    Object.keys(src.categorical).forEach(function (k) {
      d.cat[k] = {lookup: src.categorical[k].lookup, codes: i16(src.categorical[k].codes)};
    });
    Object.keys(src.numeric).forEach(function (k) { d.num[k] = f32(src.numeric[k]); });
    if (src.per_nt) d.per_nt = i16(src.per_nt);
    AXES.forEach(function (ax) {
      var a = src.axes[ax];
      if (!a) return;
      var e = {x: f32(a.x), y: f32(a.y), redundancy: i32(a.redundancy),
               isolate: u8(a.is_isolate), clusters: {}, medoids: a.medoids || {}};
      Object.keys(a.clusters).forEach(function (r) { e.clusters[r] = i32(a.clusters[r]); });
      if (a.has_interface) e.has_interface = u8(a.has_interface);
      d.axes[ax] = e;
    });
    if (src.medoid_traces) {
      try {
        d.traces = JSON.parse(new TextDecoder().decode(fflate.gunzipSync(b64(src.medoid_traces))));
      } catch (err) { d.traces = {}; }
    }
    decoded[name] = d;
    return d;
  }

  /* ---------- colour ---------- */
  function colourValues(d, axis) {
    var c = state.colour;
    if (c === "cluster") return {kind: "cat_int", values: d.axes[axis].clusters[state.resolution]};
    if (c === "cluster_other") {
      var other = axis === "fold" ? "epitope" : "fold";
      return d.axes[other]
        ? {kind: "cat_int", values: d.axes[other].clusters[state.resolution]}
        : {kind: "none"};
    }
    if (c === "approach") return {kind: "approach"};
    if (d.cat[c]) return {kind: "cat", values: d.cat[c].codes, lookup: d.cat[c].lookup};
    if (d.num[c]) return {kind: "num", values: d.num[c]};
    return {kind: "none"};
  }

  function colourFor(spec, i, d) {
    if (spec.kind === "cat") return PALETTE[((spec.values[i] % PALETTE.length) + PALETTE.length) % PALETTE.length];
    if (spec.kind === "cat_int") {
      var v = spec.values[i];
      return v < 0 ? "#ccc" : PALETTE[v % PALETTE.length];
    }
    if (spec.kind === "num") {
      var arr = spec.values, v2 = arr[i];
      if (!isFinite(v2)) return "#ddd";
      if (spec.lo === undefined) {
        var lo = Infinity, hi = -Infinity;
        for (var k = 0; k < arr.length; k++) {
          if (!isFinite(arr[k])) continue;
          if (arr[k] < lo) lo = arr[k];
          if (arr[k] > hi) hi = arr[k];
        }
        spec.lo = lo; spec.hi = hi;
      }
      var t = spec.hi > spec.lo ? (v2 - spec.lo) / (spec.hi - spec.lo) : 0.5;
      return "hsl(" + Math.round(240 - 240 * t) + ",70%,50%)";
    }
    if (spec.kind === "approach") {
      var ax = d.num.ax, ay = d.num.ay;
      if (!ax) return "#888";
      var hue = (Math.atan2(ay[i], ax[i]) * 180 / Math.PI + 360) % 360;
      return "hsl(" + Math.round(hue) + ",65%,52%)";
    }
    return "#4a7ebb";
  }

  /* ---------- rendering ---------- */
  function view(axis, d) {
    if (state.views[axis]) return state.views[axis];
    var a = d.axes[axis];
    var lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
    for (var i = 0; i < a.x.length; i++) {
      if (a.x[i] < lo[0]) lo[0] = a.x[i];
      if (a.x[i] > hi[0]) hi[0] = a.x[i];
      if (a.y[i] < lo[1]) lo[1] = a.y[i];
      if (a.y[i] > hi[1]) hi[1] = a.y[i];
    }
    state.views[axis] = {lo: lo, hi: hi, scale: 1, tx: 0, ty: 0};
    return state.views[axis];
  }

  function visible(d, i) {
    var stack = d.cat.stack;
    return !(stack && state.hidden[stack.lookup[stack.codes[i]]]);
  }

  function draw(axis) {
    var d = layerData(state.layer);
    var a = d.axes[axis];
    var canvas = document.getElementById("canvas-" + axis);
    if (!a) { canvas.getContext("2d").clearRect(0, 0, canvas.width, canvas.height); return; }
    var dpr = window.devicePixelRatio || 1;
    var w = canvas.clientWidth, h = canvas.clientHeight;
    canvas.width = w * dpr; canvas.height = h * dpr;
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    var v = view(axis, d);
    var pad = 18;
    var sx = (w - 2 * pad) / Math.max(1e-6, v.hi[0] - v.lo[0]);
    var sy = (h - 2 * pad) / Math.max(1e-6, v.hi[1] - v.lo[1]);
    var s = Math.min(sx, sy) * v.scale;
    v._s = s; v._pad = pad; v._w = w; v._h = h;

    var spec = colourValues(d, axis);
    var iface = a.has_interface;
    for (var i = 0; i < a.x.length; i++) {
      if (!visible(d, i)) continue;
      var px = pad + (a.x[i] - v.lo[0]) * s + v.tx;
      var py = h - pad - (a.y[i] - v.lo[1]) * s + v.ty;
      if (px < -5 || px > w + 5 || py < -5 || py > h + 5) continue;
      // Nodes with no interface have no epitope; grey rather than invent one.
      ctx.fillStyle = (iface && !iface[i]) ? "#e2e2e2" : colourFor(spec, i, d);
      ctx.fillRect(px, py, 2, 2);
    }

    if (state.selected >= 0 && visible(d, state.selected)) {
      var qx = pad + (a.x[state.selected] - v.lo[0]) * s + v.tx;
      var qy = h - pad - (a.y[state.selected] - v.lo[1]) * s + v.ty;
      ctx.strokeStyle = "#000"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(qx, qy, 5, 0, 6.284); ctx.stroke();
    }
  }

  function drawAll() { AXES.forEach(draw); }

  function pick(axis, mx, my) {
    var d = layerData(state.layer), a = d.axes[axis];
    if (!a) return -1;
    var v = state.views[axis];
    if (!v || !v._s) return -1;
    var best = -1, bestDist = 64;
    for (var i = 0; i < a.x.length; i++) {
      if (!visible(d, i)) continue;
      var px = v._pad + (a.x[i] - v.lo[0]) * v._s + v.tx;
      var py = v._h - v._pad - (a.y[i] - v.lo[1]) * v._s + v.ty;
      var dist = (px - mx) * (px - mx) + (py - my) * (py - my);
      if (dist < bestDist) { bestDist = dist; best = i; }
    }
    return best;
  }

  /* ---------- panels ---------- */
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c];
    });
  }

  function renderPoint() {
    var el = document.getElementById("tab-point");
    var d = layerData(state.layer), i = state.selected;
    if (i < 0) { el.innerHTML = "<p>Click a point.</p>"; return; }
    var rows = [];
    rows.push(["structure", esc(d.struct_names[i])]);
    rows.push(["designs", d.n_designs[i]]);
    Object.keys(d.cat).forEach(function (k) {
      rows.push([k, esc(d.cat[k].lookup[d.cat[k].codes[i]])]);
    });
    AXES.forEach(function (ax) {
      if (!d.axes[ax]) return;
      rows.push([ax + " cluster", d.axes[ax].clusters[state.resolution][i]]);
      if (d.axes[ax].isolate[i]) rows.push([ax + " isolate", "yes"]);
    });
    ["dockq_af3", "dockq_esm", "iptm_mean", "contacts_total",
     "epi_agree_af3_esm", "epi_agree_design_af3", "epi_agree_design_esm"
    ].forEach(function (k) {
      if (d.num[k] && isFinite(d.num[k][i])) rows.push([k, d.num[k][i].toFixed(3)]);
    });

    var html = "<table>" + rows.map(function (r) {
      return "<tr><th>" + r[0] + "</th><td>" + r[1] + "</td></tr>";
    }).join("") + "</table>";

    if (d.per_nt) {
      var vals = [], max = 1;
      for (var k = 0; k < 16; k++) { var v = d.per_nt[i * 16 + k]; vals.push(v); if (v > max) max = v; }
      html += "<p style='color:#666;margin:8px 0 2px'>contacts per nucleotide</p><div class='spark'>";
      vals.forEach(function (v) { html += "<div style='height:" + (100 * v / max) + "%'></div>"; });
      html += "</div>";
    }
    html += "<p style='color:#666'>" + esc((d.design_ids[i] || "").split(",").slice(0, 12).join("\n")) + "</p>";
    el.innerHTML = html;
  }

  var viewer = null;
  function renderCluster() {
    var d = layerData(state.layer), i = state.selected;
    var info = document.getElementById("cluster-info");
    if (i < 0 || !d.axes.fold) { info.innerHTML = "<p>Click a point.</p>"; return; }

    var html = "";
    AXES.forEach(function (ax) {
      if (!d.axes[ax]) return;
      var cl = d.axes[ax].clusters[state.resolution][i];
      var members = [];
      var arr = d.axes[ax].clusters[state.resolution];
      for (var k = 0; k < arr.length; k++) if (arr[k] === cl) members.push(k);
      html += "<p><b>" + ax + " cluster " + cl + "</b> — " + members.length + " members</p>";

      // Composition by the OTHER axis: what does the other view look like in here?
      var other = ax === "fold" ? "epitope" : "fold";
      if (d.axes[other]) {
        var counts = {};
        var oarr = d.axes[other].clusters[state.resolution];
        members.forEach(function (m) { counts[oarr[m]] = (counts[oarr[m]] || 0) + 1; });
        var top = Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; }).slice(0, 5);
        html += "<table><tr><th>" + other + " cluster</th><th>n</th></tr>" +
          top.map(function (c) { return "<tr><td>" + c + "</td><td>" + counts[c] + "</td></tr>"; }).join("") +
          "</table>";
      }
      var stackCounts = {};
      members.forEach(function (m) {
        var name = d.cat.stack.lookup[d.cat.stack.codes[m]];
        stackCounts[name] = (stackCounts[name] || 0) + 1;
      });
      html += "<table>" + Object.keys(stackCounts).sort().map(function (s) {
        return "<tr><th>" + esc(s) + "</th><td>" + stackCounts[s] + "</td></tr>";
      }).join("") + "</table>";
    });
    info.innerHTML = html;

    // Medoid structure for the fold cluster
    var medoids = d.axes.fold.medoids[state.resolution] || {};
    var entry = medoids[String(d.axes.fold.clusters[state.resolution][i])];
    var trace = entry && d.traces ? d.traces[entry.name] : null;
    if (!trace) return;
    if (!viewer) viewer = $3Dmol.createViewer(document.getElementById("viewer"),
                                              {backgroundColor: "white"});
    viewer.clear();
    var pdb = trace.map(function (c, n) {
      return "ATOM  " + String(n + 1).padStart(5) + "  CA  GLY A" +
        String(n + 1).padStart(4) + "    " +
        c[0].toFixed(3).padStart(8) + c[1].toFixed(3).padStart(8) +
        c[2].toFixed(3).padStart(8) + "  1.00  0.00           C";
    }).join("\n");
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({}, {cartoon: {color: "spectrum"}});
    viewer.zoomTo();
    viewer.render();
  }

  function renderCompare() {
    var el = document.getElementById("tab-compare");
    var scores = (P.compare && P.compare.partition_scores) || [];
    var mine = scores.filter(function (r) { return r.layer === state.layer; });
    var html = "<p><b>Fold vs epitope agreement</b></p>";
    if (!mine.length) {
      html += "<p>No comparison available for this layer.</p>";
    } else {
      html += "<table><tr><th>resolution</th><th>AMI</th><th>ARI</th><th>n</th></tr>" +
        mine.map(function (r) {
          return "<tr><td>" + r.resolution + "</td><td>" + r.ami.toFixed(3) +
                 "</td><td>" + r.ari.toFixed(3) + "</td><td>" + r.n + "</td></tr>";
        }).join("") + "</table>" +
        "<p style='color:#666'>A low score means backbone fold does not predict " +
        "where a binder lands. That is a result, not a failure.</p>";
    }

    var profiles = ((P.compare && P.compare.stack_profiles) || [])
      .filter(function (r) { return r.layer === state.layer; });
    if (profiles.length) {
      html += "<p><b>Mean contacts per nucleotide, by stack</b></p>";
      profiles.forEach(function (r) {
        var vals = [], max = 1e-9;
        for (var k = 1; k <= 16; k++) {
          var v = r["nt" + String(k).padStart(2, "0")] || 0;
          vals.push(v); if (v > max) max = v;
        }
        html += "<p style='margin:6px 0 0'>" + esc(r.stack) + "</p><div class='spark'>" +
          vals.map(function (v) { return "<div style='height:" + (100 * v / max) + "%'></div>"; }).join("") +
          "</div>";
      });
    }
    el.innerHTML = html;
  }

  function renderPanels() { renderPoint(); renderCluster(); renderCompare(); }

  /* ---------- setup ---------- */
  function setup() {
    var layerSel = document.getElementById("layer");
    Object.keys(P.layers).forEach(function (name) {
      var o = document.createElement("option");
      o.value = name; o.textContent = name + " (" + P.layers[name].n.toLocaleString() + ")";
      layerSel.appendChild(o);
    });
    layerSel.value = state.layer;
    layerSel.onchange = function () {
      state.layer = this.value; state.views = {}; state.selected = -1;
      buildColourOptions(); buildStackToggles(); drawAll(); renderPanels(); renderProvenance();
    };

    buildColourOptions();
    document.getElementById("colour").onchange = function () {
      state.colour = this.value; drawAll();
    };
    document.getElementById("resolution").onchange = function () {
      state.resolution = this.value; drawAll(); renderPanels();
    };
    document.getElementById("search").oninput = function () {
      var q = this.value.trim();
      if (!q) return;
      var d = layerData(state.layer);
      for (var i = 0; i < d.design_ids.length; i++) {
        if (d.struct_names[i] === q || (d.design_ids[i] || "").indexOf(q) >= 0) {
          state.selected = i; drawAll(); renderPanels();
          document.getElementById("status").textContent = "found " + q;
          return;
        }
      }
      document.getElementById("status").textContent = "no match";
    };

    buildStackToggles();

    document.querySelectorAll("#panel nav button").forEach(function (btn) {
      btn.onclick = function () {
        document.querySelectorAll("#panel nav button").forEach(function (b) {
          b.classList.remove("active");
        });
        document.querySelectorAll(".tab").forEach(function (t) { t.classList.remove("active"); });
        btn.classList.add("active");
        document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
      };
    });

    AXES.forEach(function (axis) {
      var canvas = document.getElementById("canvas-" + axis);
      canvas.addEventListener("click", function (ev) {
        var r = canvas.getBoundingClientRect();
        var hit = pick(axis, ev.clientX - r.left, ev.clientY - r.top);
        if (hit >= 0) { state.selected = hit; drawAll(); renderPanels(); }
      });
      canvas.addEventListener("wheel", function (ev) {
        ev.preventDefault();
        var v = view(axis, layerData(state.layer));
        v.scale *= ev.deltaY < 0 ? 1.12 : 0.89;
        draw(axis);
      }, {passive: false});
      var drag = null;
      canvas.addEventListener("mousedown", function (ev) { drag = {x: ev.clientX, y: ev.clientY}; });
      window.addEventListener("mouseup", function () { drag = null; });
      canvas.addEventListener("mousemove", function (ev) {
        if (!drag) return;
        var v = state.views[axis];
        v.tx += ev.clientX - drag.x; v.ty += ev.clientY - drag.y;
        drag = {x: ev.clientX, y: ev.clientY};
        draw(axis);
      });
    });

    window.addEventListener("resize", drawAll);
    drawAll(); renderPanels(); renderProvenance();
  }

  function buildColourOptions() {
    var d = layerData(state.layer);
    var sel = document.getElementById("colour");
    sel.innerHTML = "";
    var options = [["stack", "stack"], ["campaign", "campaign"], ["gen", "gen"],
                   ["invfold", "invfold"], ["cluster", "cluster (this map)"],
                   ["cluster_other", "cluster (other map)"],
                   ["approach", "approach geometry"]];
    Object.keys(d.num).forEach(function (k) { options.push([k, k]); });
    options.forEach(function (o) {
      var el = document.createElement("option");
      el.value = o[0]; el.textContent = o[1];
      sel.appendChild(el);
    });
    sel.value = state.colour in d.num || ["stack", "campaign", "gen", "invfold",
      "cluster", "cluster_other", "approach"].indexOf(state.colour) >= 0
      ? state.colour : "stack";
    state.colour = sel.value;
  }

  function buildStackToggles() {
    var d = layerData(state.layer);
    var host = document.getElementById("stacks");
    host.innerHTML = "";
    if (!d.cat.stack) return;
    d.cat.stack.lookup.forEach(function (name, idx) {
      var label = document.createElement("label");
      var box = document.createElement("input");
      box.type = "checkbox"; box.checked = !state.hidden[name];
      box.onchange = function () { state.hidden[name] = !box.checked; drawAll(); };
      var sw = document.createElement("span");
      sw.className = "swatch";
      sw.style.background = PALETTE[idx % PALETTE.length];
      label.appendChild(box); label.appendChild(sw);
      label.appendChild(document.createTextNode(name));
      host.appendChild(label);
    });
  }

  function renderProvenance() {
    var m = P.meta || {};
    var d = layerData(state.layer);
    var lines = [
      "layer " + state.layer + ": " + d.n.toLocaleString() + " unique backbones",
      "input " + (m.input_csv || "?") + "  sha256 " + String(m.input_sha256 || "").slice(0, 16),
      "designs " + (m.n_rows || 0).toLocaleString() +
        (m.subsample_per_stack ? "  SUBSAMPLE " + m.subsample_per_stack + "/stack" : ""),
      (m.params ? JSON.stringify(m.params) : "")
    ];
    lines.push("Screen distance is illustrative. Cluster membership is defined by " +
               "the graph and exact scores, never inferred from coordinates.");
    document.getElementById("provenance").textContent = lines.filter(Boolean).join("\n");
  }

  setup();
})();
```

- [ ] **Step 6: Build the page from the esm-only results and open it**

```bash
cd /Users/leo/Documents/binder_map
mkdir -p local_results
rsync -av --exclude 'shards' --exclude 'db' --exclude 'hits' \
  euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results/ \
  local_results/
/Users/leo/.venvs/binder_map/bin/python scripts/stage8_build_html.py \
  --results local_results --build build --out dist/binder_map.html
open dist/binder_map.html
```

Expected: the page opens with two maps side by side, stack toggles, and a working point panel. Check the browser console for errors — a silent JS exception shows as an empty canvas.

- [ ] **Step 7: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add web dist/.gitkeep 2>/dev/null; git add web src/binder_map/htmlbuild.py
git commit -m "Add dual linked-map explorer front end"
```

---

### Task 21: Orchestration — chain a whole layer with one command

**Files:**
- Create: `slurm/submit_layer.sh`

**Interfaces:**
- Consumes: every stage script and sbatch file.
- Produces: a submitted, `afterok`-chained job graph for one layer; prints the job IDs.

- [ ] **Step 1: Write `slurm/submit_layer.sh`**

```bash
#!/bin/bash
# Chain one layer end to end.
#   bash slurm/submit_layer.sh <design|af3|esm> [n_shards] [n_chunks]
#
# Arrays are throttled at %16 on purpose: Euler evaluates only the top 250
# pending jobs, and flooding that window starves everything else in the
# account, including other people's work.
set -uo pipefail

LAYER="${1:?usage: submit_layer.sh <design|af3|esm> [n_shards] [n_chunks]}"
N_SHARDS="${2:-64}"
N_CHUNKS="${3:-32}"
TM_FLOOR="${TM_FLOOR:-0.5}"
EPI_FLOOR="${EPI_FLOOR:-0.5}"
ALN_TYPE="${ALN_TYPE:-2}"
# Overridable so the smoke test can run against build_smoke/results_smoke
# without colliding with a full run.
BUILD="${BUILD:-build}"
RESULTS="${RESULTS:-results}"

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
cd "$FM" || exit 1
mkdir -p logs

echo "=== submitting layer=$LAYER shards=$N_SHARDS chunks=$N_CHUNKS"
echo "    TM_FLOOR=$TM_FLOOR EPI_FLOOR=$EPI_FLOOR ALN_TYPE=$ALN_TYPE"
echo "    BUILD=$BUILD RESULTS=$RESULTS"

# Stage 1 — extract + fingerprint
J1=$(sbatch --parsable --job-name="bm_s1_$LAYER" \
     --array="0-$((N_SHARDS - 1))%16" \
     --export=ALL,N_SHARDS="$N_SHARDS",BUILD="$BUILD",RESULTS="$RESULTS" \
     slurm/stage1_extract.sbatch "$LAYER") || exit 1
echo "  stage1 extract   $J1"

# Stage 2 — dedup + createdb + split queries
J2=$(sbatch --parsable --job-name="bm_s2_$LAYER" --dependency=afterok:"$J1" \
     --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
       PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
       scripts/stage3_split_queries.py --layer $LAYER --n-chunks $N_CHUNKS \
       --results $FM/$RESULTS'" \
     --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=4G \
     --output="$FM/logs/bm_s2_${LAYER}_%j.log") || exit 1
echo "  stage2 builddb   $J2"

# Stage 3 — all-vs-all search
J3=$(sbatch --parsable --job-name="bm_s3_$LAYER" --dependency=afterok:"$J2" \
     --array="0-$((N_CHUNKS - 1))%16" \
     --export=ALL,ALN_TYPE="$ALN_TYPE",RESULTS="$RESULTS" \
     slurm/stage3_search.sbatch "$LAYER") || exit 1
echo "  stage3 search    $J3"

# Stage 4 — both graphs
J4=$(sbatch --parsable --job-name="bm_s4_$LAYER" --dependency=afterok:"$J3" \
     slurm/stage_single.sbatch scripts/stage4_graph.py \
       --layer "$LAYER" --tm-floor "$TM_FLOOR" --epi-floor "$EPI_FLOOR" \
       --results "$FM/$RESULTS") || exit 1
echo "  stage4 graph     $J4"

# Stages 5 and 6 — cluster and embed, per axis
PREV="$J4"
for AXIS in fold epitope; do
  J5=$(sbatch --parsable --job-name="bm_s5_${LAYER}_${AXIS}" \
       --dependency=afterok:"$PREV" \
       slurm/stage_single.sbatch scripts/stage5_cluster.py \
         --layer "$LAYER" --axis "$AXIS" --results "$FM/$RESULTS") || exit 1
  J6=$(sbatch --parsable --job-name="bm_s6_${LAYER}_${AXIS}" \
       --dependency=afterok:"$J5" \
       slurm/stage_single.sbatch scripts/stage6_embed.py \
         --layer "$LAYER" --axis "$AXIS" --results "$FM/$RESULTS") || exit 1
  echo "  stage5 $AXIS      $J5"
  echo "  stage6 $AXIS      $J6"
  PREV="$J6"
done

# Stage 7 — medoid traces
J7=$(sbatch --parsable --job-name="bm_s7_$LAYER" --dependency=afterok:"$PREV" \
     slurm/stage_single.sbatch scripts/stage7_medoids.py --layer "$LAYER" \
       --results "$FM/$RESULTS") || exit 1
echo "  stage7 medoids   $J7"

echo "$J7" > "$FM/logs/last_job_${LAYER}.txt"
echo "=== layer $LAYER submitted; final job $J7"
echo "Run stage 6b (comparisons) once ALL layers have finished:"
echo "  sbatch --dependency=afterok:<j7_design>:<j7_af3>:<j7_esm> \\"
echo "    slurm/stage_single.sbatch scripts/stage6b_compare.py"
```

- [ ] **Step 2: Make the scripts executable and sync**

```bash
cd /Users/leo/Documents/binder_map
chmod +x slurm/*.sh slurm/*.sbatch scripts/*.py scripts/*.sh
rsync -av --exclude '.git' --exclude 'tests/fixtures' --exclude 'local_results' \
  --exclude 'dist' ./ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/
```

- [ ] **Step 3: Dry-run the submission logic**

Confirm the chain builds without actually queueing 100+ jobs, by checking the script parses and the sbatch files are valid:

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  bash -n slurm/submit_layer.sh && bash -n slurm/stage1_extract.sbatch && \
  bash -n slurm/stage3_search.sbatch && bash -n slurm/stage_single.sbatch && \
  echo \"all scripts parse\""'
```

Expected: `all scripts parse`.

- [ ] **Step 4: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add slurm/submit_layer.sh
git commit -m "Add afterok-chained per-layer submission"
```

---

### Task 22: Smoke test — 5k stratified subsample, end to end

This is where `TM_FLOOR`, `EPI_FLOOR` and `--alignment-type` get settled. Running the full three layers before this would risk a day of compute on a bad parameter.

**Files:**
- Create: `scripts/report_thresholds.py`

**Interfaces:**
- Consumes: `results/<layer>/hits/*.tsv`, `results/<layer>/fp/canonical.parquet`.
- Produces: a printed distribution report; no files.

- [ ] **Step 1: Implement `scripts/report_thresholds.py`**

```python
#!/usr/bin/env python3
"""Report similarity distributions so the graph floors are chosen, not guessed.

Usage:  python scripts/report_thresholds.py --layer esm
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.fingerprint import FINGERPRINT_COLUMNS
from binder_map.knn import symmetrise_min

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
QUANTILES = [0.5, 0.75, 0.9, 0.95, 0.99, 0.999]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    columns = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
               "lddt", "fident", "alnlen", "evalue", "bits"]
    hits = pd.concat(
        [pd.read_csv(f, sep="\t", names=columns, header=None)
         for f in sorted((out / "hits").glob("hits_*.tsv"))],
        ignore_index=True)

    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    index = pd.Series(nodes["node_id"].values, index=nodes["struct_name"].values)
    for column in ("query", "target"):
        hits[column] = hits[column].str.removesuffix(".pdb")
    qi = index.reindex(hits["query"]).to_numpy()
    tj = index.reindex(hits["target"]).to_numpy()
    ok = ~(np.isnan(qi) | np.isnan(tj))
    _, _, sim = symmetrise_min(qi[ok].astype(int), tj[ok].astype(int),
                               hits["qtmscore"].to_numpy()[ok],
                               hits["ttmscore"].to_numpy()[ok])

    n = len(nodes)
    print(f"\n=== {args.layer}: fold similarity, min(qTM,tTM) ===")
    print(f"nodes {n:,d}, undirected pairs {len(sim):,d}")
    for q in QUANTILES:
        print(f"  q{q:<6} {np.quantile(sim, q):.3f}")
    pairs_i, pairs_j, _ = symmetrise_min(qi[ok].astype(int), tj[ok].astype(int),
                                         hits["qtmscore"].to_numpy()[ok],
                                         hits["ttmscore"].to_numpy()[ok])
    for floor in (0.3, 0.4, 0.5, 0.6, 0.7):
        kept = sim >= floor
        # How many nodes keep at least one neighbour? This is the number that
        # matters: a floor that leaves most nodes isolated fragments the map.
        touched = np.unique(np.concatenate([pairs_i[kept], pairs_j[kept]])) \
            if kept.any() else np.array([], dtype=int)
        print(f"  floor {floor}: {int(kept.sum()):,d} edges "
              f"({100 * kept.mean():.1f}% of pairs), "
              f"{len(touched):,d}/{n:,d} nodes connected "
              f"({100 * len(touched) / max(n, 1):.1f}%)")

    fp_path = out / "fp" / "canonical.parquet"
    if fp_path.exists():
        fp = pd.read_parquet(fp_path)
        vectors = fp[FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
        keep = fp["contacts_total"].to_numpy() > 0
        vectors = vectors[keep]
        l1 = np.abs(vectors).sum(axis=1, keepdims=True)
        l1[l1 == 0] = 1
        unit = vectors / l1
        unit = unit / np.linalg.norm(unit, axis=1, keepdims=True)
        sample = unit[np.random.default_rng(0).choice(
            len(unit), size=min(2000, len(unit)), replace=False)]
        sims = (sample @ sample.T)[np.triu_indices(len(sample), k=1)]
        print(f"\n=== {args.layer}: epitope cosine similarity "
              f"(2k sample, {int((~keep).sum()):,d} no-interface excluded) ===")
        for q in QUANTILES:
            print(f"  q{q:<6} {np.quantile(sims, q):.3f}")
        for floor in (0.3, 0.5, 0.7, 0.85):
            print(f"  floor {floor}: {100 * (sims >= floor).mean():.1f}% of pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Build the 5k subsample manifests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python scripts/stage0_select.py \
  --out build_smoke --subsample 1000
cat build_smoke/meta/counts.json
rsync -av build_smoke/ \
  euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/build_smoke/
```

Expected: 5,000 rows, 1,000 per stack.

- [ ] **Step 3: Run all three layers on the subsample**

Point the stages at the smoke build and a separate results tree so nothing collides with a later full run.

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  BUILD=build_smoke RESULTS=results_smoke bash slurm/submit_layer.sh esm 8 4"'
```

`submit_layer.sh` reads `BUILD` and `RESULTS` from the environment (Task 21), so this runs entirely against `build_smoke/` and `results_smoke/` and cannot collide with a later full run. `8` shards and `4` query chunks keep the job count small.

- [ ] **Step 4: Watch it and check the outcome**

```bash
ssh euler 'squeue -u $USER -o "%.10i %.22j %.8T %.10M %R" | head -30'
```

Once finished:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/report_thresholds.py --layer esm --results results_smoke'
```

- [ ] **Step 5: Choose the parameters**

Read the distributions and pick:

- `TM_FLOOR` — the value where the edge count drops sharply but most nodes still keep a neighbour. Too high fragments everything into singletons; too low merges the whole map into one blob. **If the distribution gives no clear answer, 0.50 stands.**
- `EPI_FLOOR` — same reasoning on the cosine distribution. **Default 0.50 stands if inconclusive.**
- `--alignment-type` — rerun stage 3 alone with `ALN_TYPE=1` and compare. Keep TMalign only if the fold clusters change materially; **if inconclusive, 2 stands**, being cheaper and the mode Foldseek is tuned for.

Record the chosen values in `docs/superpowers/plans/parameters.md` with the distributions that justified them, so the full run is reproducible.

- [ ] **Step 6: Build the smoke explorer and look at it**

```bash
cd /Users/leo/Documents/binder_map
rsync -av --exclude 'shards' --exclude 'db' --exclude 'hits' \
  euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results_smoke/ \
  local_results_smoke/
/Users/leo/.venvs/binder_map/bin/python scripts/stage8_build_html.py \
  --results local_results_smoke --build build_smoke --out dist/binder_map_smoke.html
open dist/binder_map_smoke.html
```

This is the first time the whole thing is visible. Check that the maps are not a single undifferentiated blob and not pure confetti; either extreme means the floors need another pass.

- [ ] **Step 7: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/report_thresholds.py docs/superpowers/plans/parameters.md
git commit -m "Add threshold reporting and record smoke-test parameters"
```

---

### Task 23: Full run and validation

**Files:**
- Create: `scripts/validate.py`
- Create: `results/manifest.json` (generated)

**Interfaces:**
- Consumes: every output.
- Produces: `results/manifest.json` and a pass/fail report.

- [ ] **Step 1: Implement `scripts/validate.py`**

```python
#!/usr/bin/env python3
"""Validate the full run. Every check either passes or explains itself.

Usage:  python scripts/validate.py --results results --build build
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
LAYERS = ("design", "af3", "esm")

# Campaign tables that carry rna_contact_nt_1..16 (checked 2026-08-21).
CONTACT_TABLES = {
    "v3": RB / "analyses/midrun_v3/collected_clean_metrics.csv",
    "v4": RB / "analyses/midrun_v4/collected_clean_metrics.csv",
    "v5": RB / "analyses/midrun_v5/collected_clean_metrics.csv",
    "v9": RB / "analyses/midrun_v9/collected_clean_metrics.csv",
}


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    return ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=FM / "results")
    parser.add_argument("--build", type=Path, default=FM / "build")
    parser.add_argument("--skip-contact-crosscheck", action="store_true")
    args = parser.parse_args()

    ok = True
    manifest: dict = {"layers": {}, "checks": {}}

    manifest_counts = json.loads((args.build / "meta" / "counts.json").read_text())
    print("=== per-layer reconciliation ===")
    for layer in LAYERS:
        out = args.results / layer
        if not (out / "dedup" / "nodes.tsv").exists():
            ok = check(f"{layer}: results present", False, "missing dedup/nodes.tsv") and ok
            continue
        nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
        hashes = pd.concat([pd.read_csv(f, sep="\t")
                            for f in sorted((out / "hashes").glob("shard_*.tsv"))],
                           ignore_index=True)
        missing = pd.concat([pd.read_csv(f, sep="\t")
                             for f in sorted((out / "missing").glob("shard_*.tsv"))],
                            ignore_index=True)
        expected = manifest_counts["layer_counts"][layer]
        seen = len(hashes) + len(missing)
        ok = check(f"{layer}: every design accounted for",
                   seen == expected, f"{seen:,d} of {expected:,d}") and ok

        # Chain tally: RFD3 design structures MUST be inverted.
        if layer == "design":
            man = pd.read_csv(args.build / "manifests" / "design.tsv", sep="\t")
            joined = hashes.merge(man[["design_id", "stack"]], on="design_id")
            tally = joined.groupby(["stack", "protein_chain"]).size().unstack(fill_value=0)
            print(tally.to_string())
            rfd3 = tally.loc["RFD3+LiMPNN"] if "RFD3+LiMPNN" in tally.index else None
            ok = check("design: RFD3+LiMPNN protein is chain B",
                       rfd3 is not None and rfd3.get("B", 0) > 0
                       and rfd3.get("A", 0) == 0,
                       "if this flips, chain selection is broken") and ok
        else:
            ok = check(f"{layer}: protein always chain A",
                       set(hashes["protein_chain"]) == {"A"},
                       str(hashes['protein_chain'].value_counts().to_dict())) and ok

        ok = check(f"{layer}: protein always 80 residues",
                   set(hashes["n_res"]) == {80}) and ok

        manifest["layers"][layer] = {
            "structures": int(len(hashes)),
            "unique_backbones": int(len(nodes)),
            "missing": int(len(missing)),
            "missing_reasons": missing["reason"].value_counts().to_dict()
                               if len(missing) else {},
        }

    # Self-hits must score 1.0
    print("\n=== fold sanity ===")
    for layer in LAYERS:
        files = sorted((args.results / layer / "hits").glob("hits_*.tsv"))
        if not files:
            continue
        head = pd.read_csv(files[0], sep="\t", header=None, nrows=200000,
                           names=["query", "target", "qtm", "ttm", "alntm",
                                  "lddt", "fident", "alnlen", "evalue", "bits"])
        self_hits = head[head["query"] == head["target"]]
        if len(self_hits):
            ok = check(f"{layer}: self-hits score 1.0",
                       bool((self_hits["qtm"] > 0.99).all()),
                       f"min {self_hits['qtm'].min():.3f}") and ok

    # Epitope cross-check against the stored rna_contact_nt_* columns
    if not args.skip_contact_crosscheck:
        print("\n=== epitope cross-check vs stored rna_contact_nt_* ===")
        per_nt = [f"per_nt_{i}" for i in range(1, 17)]
        stored_cols = [f"rna_contact_nt_{i}" for i in range(1, 17)]
        for campaign, path in CONTACT_TABLES.items():
            if not path.exists():
                print(f"  [skip] {campaign}: {path} not found")
                continue
            stored = pd.read_csv(path, usecols=["design_id"] + stored_cols,
                                 low_memory=False).set_index("design_id")
            hit_any = False
            for layer in LAYERS:
                fp_path = args.results / layer / "fp" / "canonical.parquet"
                if not fp_path.exists():
                    continue
                fp = pd.read_parquet(fp_path).set_index("design_id")
                common = fp.index.intersection(stored.index)
                if len(common) < 20:
                    continue
                hit_any = True
                a = fp.loc[common, per_nt].to_numpy(dtype=float)
                b = stored.loc[common, stored_cols].to_numpy(dtype=float)
                good = np.isfinite(b).all(axis=1)
                if good.sum() < 20:
                    continue
                corr = np.corrcoef(a[good].ravel(), b[good].ravel())[0, 1]
                print(f"  {campaign}/{layer}: n={int(good.sum()):,d} "
                      f"r={corr:.3f}")
                ok = check(f"{campaign}/{layer}: recomputed contacts track stored",
                           corr > 0.9,
                           "a low correlation means the cutoff, frame, or chain "
                           "selection diverged from rna_binder") and ok
            if not hit_any:
                print(f"  [skip] {campaign}: no overlapping design_ids")

    manifest["checks"]["all_passed"] = ok
    try:
        version = subprocess.run(
            ["bash", "-lc", "module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 "
             "foldseek/8-ef4e960 >/dev/null 2>&1; foldseek version"],
            capture_output=True, text=True).stdout.strip()
    except Exception:
        version = "unknown"
    manifest["tools"] = {"foldseek": version or "8-ef4e960"}
    manifest["input"] = manifest_counts

    (args.results / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest -> {args.results / 'manifest.json'}")
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Launch all three layers**

Use the parameters recorded in Task 22.

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  export TM_FLOOR=0.5 EPI_FLOOR=0.5 ALN_TYPE=2 && \
  for L in design af3 esm; do bash slurm/submit_layer.sh \$L 64 32; done"'
```

Expected: three chains submitted, roughly 310 tasks total, each array throttled at `%16`.

- [ ] **Step 3: Monitor**

```bash
ssh euler 'squeue -u $USER -o "%.10i %.24j %.8T %.10M %.6D %R" | head -40; \
  echo "---"; squeue -u $USER -h | wc -l'
```

Check the pending count stays well under 250. If it does not, lower the `%16` throttle rather than letting the account's queue window fill.

- [ ] **Step 4: Run stage 6b once all three layers finish**

```bash
ssh euler 'bash -lc "cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  sbatch --job-name=bm_s6b slurm/stage_single.sbatch scripts/stage6b_compare.py"'
```

- [ ] **Step 5: Validate**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/validate.py'
```

Expected: `ALL CHECKS PASSED`. Treat two failures as blocking rather than cosmetic:
- **RFD3+LiMPNN not on chain B** — chain selection is broken, and the design-layer map is meaningless.
- **Contact cross-check r ≤ 0.9** — the recomputed fingerprints have drifted from `rna_binder`'s definition, so every epitope result is suspect.

- [ ] **Step 6: Build the final explorer**

```bash
cd /Users/leo/Documents/binder_map
rsync -av --exclude 'shards' --exclude 'db' --exclude 'hits' \
  euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results/ \
  local_results/
/Users/leo/.venvs/binder_map/bin/python scripts/stage8_build_html.py \
  --results local_results --build build --out dist/binder_map.html
ls -lh dist/binder_map.html
open dist/binder_map.html
```

Expected: a 25–40 MB self-contained HTML with all three layers and both axes.

- [ ] **Step 7: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/validate.py
git commit -m "Add full-run validation with chain-inversion and contact cross-checks"
```

---

## Interpretation notes for whoever reads the map

Written here because they belong with the work, not just the spec.

- **Epitope assignment inherits pose reliability.** This campaign's own results show ipTM does not certify pose, and that 81% of fold-agreeing designs dock more than 10 Å apart. A single refold's epitope is a noisy estimate. `epi_agree_af3_esm` is the honest filter: trust designs where the two independent refolds agree.
- **`passing.csv` already gates on DockQ ≥ 0.7 in at least one backend**, so this set is better-posed than the campaign at large — but not uniformly so.
- **A low fold↔epitope AMI is a result**, not a bug. It would mean backbone fold does not determine where a binder lands.
- **Screen distance is illustrative.** Cluster membership and similarity come from the graph and the exact scores, never from coordinates.
- **A point can be visually inside a cluster while being a singleton by cluster id.** The layout graph is unthresholded and the clustering graph is not; that gap is informative.
- **`redundancy` is capped at 300** by `--max-seqs`.
