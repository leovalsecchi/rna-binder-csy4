# BG+BIF Library Rescore Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine how many of the 2,400 BG+BIF designs in the delivered 12k library fail the campaign's own gate when scored against the true BoltzGen backbone instead of against BIF's own re-prediction.

**Architecture:** All pure logic lives in one unit-tested module (`rescore_lib.py`); five thin CLIs drive it. The expensive archive read is paid exactly once into a small persistent working set (~174 GB read → ~0.45 GB kept), after which all scoring re-runs touch no archive. Scoring reuses the campaign's own scripts unchanged so results are comparable by construction.

**Tech Stack:** Python 3.12, biotite 1.6.0, zstandard 0.23.0, pytest 8.3.4, DockQ (vendored), SLURM.

**Spec:** `~/docs/2026-09-30-bif-library-rescore-design.md` (mirrored to `$RB/docs/`)

## Global Constraints

- `RB = /cluster/work/wicky/lvalsecchi/rna_binder`. All work happens under `$RB/analyses/bif_library_rescore/` (call it `$W`).
- **There is no git repository on Euler for `rna_binder`** (verified: no `.git` at `$RB` or `$RB/code`). The skill's "commit" step is replaced by an explicit verification checkpoint in each task. Never run `git` against `$RB`.
- Two different pythons, and they are not interchangeable:
  - `PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python` — biotite, numpy, zstandard, pytest. Use for everything except DockQ.
  - `python3` (Euler host, 3.12) with `PYTHONPATH=$RB/code/pylibs_dockq:$RB/code/scratch_dockq` and `PATH=$RB/code/bin:$PATH` — **required** by `dockq_self.py` (boltzgen.sif is py3.11 and DockQ ships no cp311 wheel).
- `lrmsd_rna_aligned.py` needs `PYTHONPATH=$RB/code/src:$RB/code/scripts`.
- SLURM on this cluster **rejects `--mem`**; use `--mem-per-cpu`.
- Never write scratch under `/cluster` from a job; `archive_read.py` handles node-local scratch itself. Set `ARCH_READ_CACHE_MAX=2` if any archive-aware tool is invoked.
- Campaign gate terms, exact: `min(sc_RMSD_af3, sc_RMSD_esm) <= 3.0` **AND** `max(DockQ_af3, DockQ_esm) >= 0.7`.
- Recombination formula, exact (from `dockq_rna_recombine.py`):
  `LRMSD_rna = dockq_LRMSD_stock * (lrmsd_rna / lrmsd_prot)`;
  `DockQ_rna = (fnat + 1/(1+(iRMSD/1.5)**2) + 1/(1+(LRMSD_rna/8.5)**2)) / 3`
- Expected counts, asserted throughout: **2,400** designs · **1,353** uuids · **2,292** `spec` · **107** `limpnn_fallback` · **1** `unresolvable`.
- Identify protein vs RNA chains **by residue composition, never by chain ID**.

### Transferring files to Euler

The `ssh euler '...'` single-quote form **strips single quotes inside the remote
command**, which silently turns `r['pair']` into `r[pair]` and fails at runtime
(`ast.parse` does not catch it — this already cost one failed array job). Always
transfer via base64, in two separate commands:

```bash
base64 < local.py | tr -d '\n' > /tmp/x.b64
B=$(cat /tmp/x.b64); ssh euler "echo $B | base64 -d > /path/on/euler.py"
```

Then verify with `md5sum` on both sides before running anything.

---

## File Structure

| file | responsibility |
|---|---|
| `$W/rescore_lib.py` | **all** pure logic: id parsing, ref-kind classification, sequence extraction, backbone transplant, recombination formula, gate. No CLI, no SLURM. |
| `$W/tests/test_rescore_lib.py` | pytest unit tests for the above |
| `$W/build_manifest.py` | library CSV + backbone/LiMPNN indexes → `manifest.tsv` |
| `$W/extract_backbones.py` | loose `_backbones/` archives → `backbones/` |
| `$W/extract_refolds.py` | one solid midrun archive → `cifs/` (one task per archive) |
| `$W/extract.sbatch` | SLURM array driving `extract_refolds.py` |
| `$W/build_refs.py` | transplant → `refs/`, enforces the sequence gate |
| `$W/score.py` | drives `dockq_self.py` + `lrmsd_rna_aligned.py`, recombines |
| `$W/report.py` | re-gate → `rescore_2400.csv` + `REPORT.md` |

Directories produced: `manifest.tsv`, `backbones/`, `cifs/`, `refs/`, `score_in/`, `score_out/`, `logs/`.

---

## Task 1: `rescore_lib.py` — pure logic, unit-tested

**Files:**
- Create: `$W/rescore_lib.py`
- Test: `$W/tests/test_rescore_lib.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `parse_design_id(design_id: str) -> tuple[str, int, int]` → `(uuid8, design_idx, seq_idx)`
  - `classify_ref_kind(uuid8: str, design_idx: int, spec_index: dict[str, set[int]], limpnn_index: set[tuple[str, int]]) -> str` → `"spec"` | `"limpnn_fallback"` | `"unresolvable"`
  - `dockq_rna(fnat: float, irmsd: float, lrmsd: float) -> float`
  - `lrmsd_rna_corrected(stock: float, rna: float, prot: float) -> float | None`
  - `gate_pass(sc_af3: float, sc_esm: float, dq_af3: float, dq_esm: float) -> bool`
  - `seq_one_letter(atoms) -> str` (biotite AtomArray → 1-letter protein sequence, ordered by `res_id`)
  - `split_prot_rna(atoms) -> tuple[AtomArray, AtomArray]` (by residue composition)
  - `prot_backbone(prot) -> AtomArray` (N, CA, C, O only)
  - `transplant(prot_bb, seq: str) -> AtomArray` (residue names set positionally)

- [ ] **Step 1: Create the test directory and write the failing tests**

Create `$W/tests/test_rescore_lib.py`:

```python
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
import numpy as np
from rescore_lib import (parse_design_id, classify_ref_kind, dockq_rna,
                         lrmsd_rna_corrected, gate_pass, seq_one_letter,
                         split_prot_rna, prot_backbone, transplant)


def test_parse_design_id():
    assert parse_design_id("csy4free_BG_BIF_len80_0e7c80af_design5_seq1") == ("0e7c80af", 5, 1)
    assert parse_design_id("csy4free_BG_BIF_len80_021f712e_design17_seq8") == ("021f712e", 17, 8)


def test_parse_design_id_rejects_other_stacks():
    with pytest.raises(ValueError):
        parse_design_id("csy4free_BG_LiMPNN_len80_0e7c80af_design5_seq1")


def test_classify_ref_kind():
    spec = {"aaaaaaaa": {0, 5, 7}}
    limpnn = {("bbbbbbbb", 3)}
    assert classify_ref_kind("aaaaaaaa", 5, spec, limpnn) == "spec"
    assert classify_ref_kind("bbbbbbbb", 3, spec, limpnn) == "limpnn_fallback"
    assert classify_ref_kind("cccccccc", 9, spec, limpnn) == "unresolvable"
    # present uuid but the needed design index is absent -> not "spec"
    assert classify_ref_kind("aaaaaaaa", 4, spec, limpnn) == "unresolvable"


def test_dockq_rna_reproduces_a_real_campaign_value():
    # pair01 / af3 / reference R1, from the 20-pair probe:
    # fnat 0.7551, iRMSD 1.5250, sc_RMSD 2.9129 -> stored DockQ 0.7139
    assert dockq_rna(0.7551, 1.5250, 2.9129) == pytest.approx(0.7139, abs=5e-5)


def test_dockq_rna_bounds():
    assert dockq_rna(1.0, 0.0, 0.0) == pytest.approx(1.0)
    assert dockq_rna(0.0, 1e6, 1e6) == pytest.approx(0.0, abs=1e-6)


def test_lrmsd_rna_corrected():
    assert lrmsd_rna_corrected(2.0, 3.0, 1.5) == pytest.approx(4.0)
    assert lrmsd_rna_corrected(2.0, 3.0, 0.0) is None      # guard on prot ~ 0


def test_gate_pass_uses_either_backend():
    # sc passes on esm only, dockq passes on af3 only -> still a pass
    assert gate_pass(9.0, 2.5, 0.72, 0.30) is True
    # sc fails on both
    assert gate_pass(9.0, 4.0, 0.99, 0.99) is False
    # dockq fails on both
    assert gate_pass(1.0, 1.0, 0.69, 0.10) is False


def test_gate_pass_boundaries_are_inclusive():
    assert gate_pass(3.0, 3.0, 0.7, 0.7) is True
    assert gate_pass(3.001, 3.001, 0.7, 0.7) is False
    assert gate_pass(3.0, 3.0, 0.699, 0.699) is False


def _toy():
    """2 protein residues (N,CA,C,O + a sidechain atom) + 1 nucleotide."""
    import biotite.structure as struc
    a = struc.AtomArray(11)
    a.coord = np.arange(33, dtype=float).reshape(11, 3)
    a.chain_id = np.array(["A"] * 9 + ["B"] * 2)
    a.res_id = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1])
    a.res_name = np.array(["GLY"] * 5 + ["ALA"] * 4 + ["G", "G"])
    a.atom_name = np.array(["N", "CA", "C", "O", "CB", "N", "CA", "C", "O", "P", "C1'"])
    a.element = np.array(["N", "C", "C", "O", "C", "N", "C", "C", "O", "P", "C"])
    return a


def test_split_prot_rna_by_composition():
    p, r = split_prot_rna(_toy())
    assert len(p) == 9 and len(r) == 2
    assert set(np.unique(r.res_name)) == {"G"}


def test_prot_backbone_drops_sidechains():
    p, _ = split_prot_rna(_toy())
    bb = prot_backbone(p)
    assert len(bb) == 8
    assert "CB" not in set(bb.atom_name)


def test_seq_one_letter():
    p, _ = split_prot_rna(_toy())
    assert seq_one_letter(p) == "GA"


def test_transplant_sets_names_positionally_and_keeps_coords():
    p, _ = split_prot_rna(_toy())
    bb = prot_backbone(p)
    out = transplant(bb, "WK")
    assert seq_one_letter(out) == "WK"
    assert np.array_equal(out.coord, bb.coord)


def test_transplant_rejects_length_mismatch():
    p, _ = split_prot_rna(_toy())
    with pytest.raises(ValueError):
        transplant(prot_backbone(p), "WKD")
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python -m pytest tests/ -q'
```

Expected: collection error — `ModuleNotFoundError: No module named 'rescore_lib'`.

- [ ] **Step 3: Write `rescore_lib.py`**

```python
#!/usr/bin/env python3
"""Pure logic for the BG+BIF library rescore.

Kept free of I/O and SLURM so it can be unit-tested; every CLI in this
directory is a thin wrapper over these functions.
"""
from __future__ import annotations

import re

import numpy as np
import biotite.structure as struc

AA3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
       "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
       "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
       "Y": "TYR", "V": "VAL"}
ONE = {v: k for k, v in AA3.items()}
AA = set(AA3.values())
NT = {"A", "C", "G", "U", "ADE", "CYT", "GUA", "URI"}
PROT_BB = ("N", "CA", "C", "O")

DESIGN_RE = re.compile(r"^csy4free_BG_BIF_len80_([0-9a-f]{8})_design(\d+)_seq(\d+)$")


def parse_design_id(design_id):
    """(uuid8, design_idx, seq_idx) for a BG+BIF design id."""
    m = DESIGN_RE.match(design_id)
    if not m:
        raise ValueError(f"not a BG+BIF design id: {design_id!r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def classify_ref_kind(uuid8, design_idx, spec_index, limpnn_index):
    """Which reference this design can be scored against.

    spec_index:   {uuid8 -> set of spec_<D> indices present in _backbones/}
    limpnn_index: {(uuid8, design_idx)} that have a BG+LiMPNN design
    """
    if design_idx in spec_index.get(uuid8, ()):
        return "spec"
    if (uuid8, design_idx) in limpnn_index:
        return "limpnn_fallback"
    return "unresolvable"


def dockq_rna(fnat, irmsd, lrmsd):
    """DockQ in the binder-campaign frame (dockq_rna_recombine.py)."""
    return (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3


def lrmsd_rna_corrected(stock, rna, prot):
    """Anchor the RNA-aligned LRMSD to DockQ's own output via the ratio, so the
    two implementations' residue-matching differences cancel. None if the
    denominator is degenerate."""
    if prot is None or prot <= 0.01:
        return None
    return stock * (rna / prot)


def gate_pass(sc_af3, sc_esm, dq_af3, dq_esm):
    """Campaign base gate, the two terms a reference change can move.
    Each term passes on EITHER backend; the terms are ANDed."""
    return bool(min(sc_af3, sc_esm) <= 3.0 and max(dq_af3, dq_esm) >= 0.7)


def split_prot_rna(atoms):
    """Protein and RNA sub-arrays, chosen by residue COMPOSITION.
    Chain lettering is not consistent across stacks, so never key on chain_id."""
    return atoms[np.isin(atoms.res_name, list(AA))], atoms[np.isin(atoms.res_name, list(NT))]


def prot_backbone(prot):
    return prot[np.isin(prot.atom_name, list(PROT_BB))]


def _res_order(atoms):
    return sorted({int(r) for r in atoms.res_id})


def seq_one_letter(prot):
    """1-letter protein sequence, residues ordered by res_id."""
    first = {}
    for rid, rname in zip(prot.res_id, prot.res_name):
        first.setdefault(int(rid), str(rname))
    return "".join(ONE[first[r]] for r in _res_order(prot))


def transplant(prot_bb, seq):
    """Set residue NAMES positionally from `seq`, leaving coordinates untouched.

    Legitimate because the backbone is the same chain and inverse folding assigns
    a sequence to those same positions -- the relationship the pipeline already
    assumes when it treats a design file as ground truth. Validated against real
    BG+LiMPNN design files (0.0005 A, 320/320 residue names, 20/20 pairs).
    """
    ids = _res_order(prot_bb)
    if len(ids) != len(seq):
        raise ValueError(f"{len(ids)} residues vs {len(seq)}-aa sequence")
    out = prot_bb.copy()
    m = {rid: AA3[seq[i]] for i, rid in enumerate(ids)}
    out.res_name = np.array([m[int(r)] for r in out.res_id], dtype=out.res_name.dtype)
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python -m pytest tests/ -q'
```

Expected: `13 passed`.

- [ ] **Step 5: Checkpoint**

Record in `$W/PROGRESS.md`: `Task 1 done — rescore_lib.py, 13 unit tests passing.`
The `test_dockq_rna_reproduces_a_real_campaign_value` test is the load-bearing one: it pins the recombination formula to a value the campaign actually stored.

---

## Task 2: `build_manifest.py` — the 2,400 rows

**Files:**
- Create: `$W/build_manifest.py`
- Output: `$W/manifest.tsv`

**Interfaces:**
- Consumes: `rescore_lib.parse_design_id`, `rescore_lib.classify_ref_kind`.
- Produces: `manifest.tsv`, tab-separated, columns in this exact order:
  `design_id, uuid8, design_idx, seq_idx, campaign, ref_kind, sequence, camp_sc_rmsd_af3, camp_sc_rmsd_esm, camp_dockq_af3, camp_dockq_esm`

- [ ] **Step 1: Write `build_manifest.py`**

```python
#!/usr/bin/env python3
"""Build the 2,400-row manifest for the BG+BIF library rescore.

Pure metadata: reads the delivered library CSV, the loose _backbones/ .list
sidecars, and the filter store's design_id column. Touches no archive.
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescore_lib import classify_ref_kind, parse_design_id

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
LIB = f"{RB}/analyses/library_selection/csy4_binder_library_12k_v6_k30.csv"
BBDIR = f"{RB}/campaigns/csy4_crrna_free/_backbones"
STORE = f"{RB}/analyses/cross_campaign/filter_store_v2v9_hdtop"
SPEC_RE = re.compile(r"/spec_(\d+)\.cif$")
LIMPNN_RE = re.compile(r"^csy4free_BG_LiMPNN_len80_([0-9a-f]{8})_design(\d+)_seq(\d+)$")

COLS = ["design_id", "uuid8", "design_idx", "seq_idx", "campaign", "ref_kind",
        "sequence", "camp_sc_rmsd_af3", "camp_sc_rmsd_esm",
        "camp_dockq_af3", "camp_dockq_esm"]


def spec_index(uuids):
    """{uuid8 -> set of spec_<D> indices}, read from the loose .list sidecars."""
    out = {}
    for u in uuids:
        lst = f"{BBDIR}/BG_len80_{u}.tar.zst.list"
        if not os.path.exists(lst):
            continue
        idx = set()
        with open(lst) as f:
            for line in f:
                m = SPEC_RE.search(line.strip())
                if m:
                    idx.add(int(m.group(1)))
        out[u] = idx
    return out


def limpnn_index():
    """{(uuid8, design_idx)} that carry a BG+LiMPNN design."""
    meta = json.load(open(f"{STORE}/meta.json"))
    stack = np.load(f"{STORE}/stack.npy")
    did = np.load(f"{STORE}/design_id.npy")
    out = set()
    for i in np.flatnonzero(stack == meta["vocabs"]["stack"].index("BG+LiMPNN")):
        m = LIMPNN_RE.match(did[i].decode())
        if m:
            out.add((m.group(1), int(m.group(2))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(open(LIB)) if r["stack"] == "BG+BIF"]
    assert len(rows) == 2400, f"expected 2400 BG+BIF library rows, got {len(rows)}"

    parsed = []
    for r in rows:
        u, d, s = parse_design_id(r["design_id"])
        parsed.append((u, d, s, r))
    uuids = sorted({u for u, _, _, _ in parsed})
    assert len(uuids) == 1353, f"expected 1353 uuids, got {len(uuids)}"

    si = spec_index(uuids)
    li = limpnn_index()

    out = []
    for u, d, s, r in parsed:
        out.append({
            "design_id": r["design_id"], "uuid8": u, "design_idx": d, "seq_idx": s,
            "campaign": r["campaign"], "ref_kind": classify_ref_kind(u, d, si, li),
            "sequence": r["sequence"],
            "camp_sc_rmsd_af3": r["sc_rmsd_af3"], "camp_sc_rmsd_esm": r["sc_rmsd_esm"],
            "camp_dockq_af3": r["dockq_af3"], "camp_dockq_esm": r["dockq_esm"]})

    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    kinds = collections.Counter(r["ref_kind"] for r in out)
    print(f"wrote {a.out}: {len(out)} rows")
    print(f"ref_kind: {dict(kinds)}")
    print(f"campaigns: {dict(sorted(collections.Counter(r['campaign'] for r in out).items()))}")
    assert kinds["spec"] == 2292, f"expected 2292 spec, got {kinds['spec']}"
    assert kinds["limpnn_fallback"] == 107, f"expected 107 fallback, got {kinds['limpnn_fallback']}"
    assert kinds["unresolvable"] == 1, f"expected 1 unresolvable, got {kinds['unresolvable']}"
    print("counts match the spec")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python build_manifest.py --out manifest.tsv'
```

Expected, exactly:
```
wrote manifest.tsv: 2400 rows
ref_kind: {'spec': 2292, 'limpnn_fallback': 107, 'unresolvable': 1}
campaigns: {'v2': 235, 'v4': 138, 'v5': 158, 'v6': 133, 'v7': 542, 'v8': 629, 'v9': 565}
counts match the spec
```

If the counts differ, **stop** — the underlying data moved since 2026-09-30 and the spec's facts need re-establishing before anything else runs.

- [ ] **Step 3: Checkpoint**

Append to `$W/PROGRESS.md`: `Task 2 done — manifest.tsv, 2292/107/1, counts asserted.`

---

## Task 3: `extract_backbones.py` — the true backbones

**Files:**
- Create: `$W/extract_backbones.py`
- Output: `$W/backbones/<uuid8>_spec<D>.cif`

**Interfaces:**
- Consumes: `manifest.tsv`.
- Produces: one CIF per `ref_kind == "spec"` row, named `<uuid8>_spec<D>.cif`. `build_refs.py` (Task 5) reads exactly this name.

These archives are loose per-uuid files that the 2026-09-11 repack did not touch, so this is cheap and needs no SLURM.

- [ ] **Step 1: Write `extract_backbones.py`**

```python
#!/usr/bin/env python3
"""Pull spec_<D>.cif for every manifest row whose ref_kind == 'spec'.

Source archives are the loose per-uuid _backbones/BG_len80_<uuid>.tar.zst, which
the 2026-09-11 repack left alone. Index padding is inconsistent across runs
(spec_7 vs spec_07), so members are matched on the PARSED INTEGER.
"""
import argparse
import collections
import csv
import os
import re
import tarfile

import zstandard as zstd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
BBDIR = f"{RB}/campaigns/csy4_crrna_free/_backbones"
SPEC = re.compile(r"^spec_(\d+)\.cif$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--expect", type=int, default=None,
                    help="assert this many backbones end up on disk; omit when "
                         "running on a subset of the manifest")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    want = collections.defaultdict(set)
    for r in csv.DictReader(open(a.manifest), delimiter="\t"):
        if r["ref_kind"] == "spec":
            want[r["uuid8"]].add(int(r["design_idx"]))

    dctx = zstd.ZstdDecompressor()
    n_ok = n_skip = 0
    for u, idxs in sorted(want.items()):
        todo = {d for d in idxs
                if not os.path.exists(os.path.join(a.outdir, f"{u}_spec{d}.cif"))}
        n_skip += len(idxs) - len(todo)
        if not todo:
            continue
        arc = f"{BBDIR}/BG_len80_{u}.tar.zst"
        with open(arc, "rb") as fh, dctx.stream_reader(fh) as rd:
            with tarfile.open(fileobj=rd, mode="r|") as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    mm = SPEC.match(os.path.basename(m.name))
                    if not mm:
                        continue
                    d = int(mm.group(1))
                    if d not in todo:
                        continue
                    with open(os.path.join(a.outdir, f"{u}_spec{d}.cif"), "wb") as o:
                        o.write(tf.extractfile(m).read())
                    todo.discard(d)
                    n_ok += 1
                    if not todo:
                        break
        if todo:
            print(f"  MISSING for {u}: spec_{sorted(todo)}")

    total = sum(len(v) for v in want.values())
    print(f"backbones: {n_ok} written, {n_skip} already present, {total} wanted")
    have = sum(1 for r in csv.DictReader(open(a.manifest), delimiter="\t")
               if r["ref_kind"] == "spec"
               and os.path.exists(os.path.join(a.outdir, f"{r['uuid8']}_spec{r['design_idx']}.cif")))
    print(f"coverage: {have}/{total}")
    if a.expect is not None:
        assert have == a.expect, f"expected {a.expect} backbones on disk, got {have}"


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test on a few uuids before the full run**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && head -4 manifest.tsv > /dev/shm/m3.tsv && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python extract_backbones.py --manifest /dev/shm/m3.tsv --outdir backbones 2>&1 | head'
```

Expected: writes up to 3 files and exits 0, printing `coverage: N/N` for the subset.
No `--expect` is passed, so no assertion runs. Confirm the files appeared in `backbones/`.

- [ ] **Step 3: Full run**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python extract_backbones.py --manifest manifest.tsv --outdir backbones --expect 2292'
```

Expected final lines: `coverage: 2292/2292` and no assertion error.

- [ ] **Step 4: Checkpoint**

Append: `Task 3 done — backbones/ 2292 files.`

---

## Task 4: `extract_refolds.py` + `extract.sbatch` — the one archive pass

**Files:**
- Create: `$W/extract_refolds.py`, `$W/extract.sbatch`
- Output: `$W/cifs/<design_id>_{af3,esm}.cif`, optionally `$W/designs/<design_id>.cif`

**Interfaces:**
- Consumes: `manifest.tsv`.
- Produces: `cifs/<design_id>_af3.cif` and `cifs/<design_id>_esm.cif`. `build_refs.py` and `score.py` read exactly these names.

Member layout inside a midrun archive, and the reason a full stream is unavoidable:
tar order is `designs/ → refold_af3/ → _refold_arrays/ → refold_b2/ → refold_esm/ → run_*/`,
and the archive is a **single** zstd frame.

- [ ] **Step 1: Write `extract_refolds.py`**

```python
#!/usr/bin/env python3
"""Pull AF3 + ESM refold CIFs for one midrun archive, in a single stream.

Post-2026-09-11 the refs are DOUBLE-NESTED:
  <stack>/<midrun>.tar.zst ! <midrun>/refold_<be>/msaoff_<uuid>/designs.tar.zst ! designs/<id>_<BE>....cif
A solid zstd stream has no random access, so: one task = one outer archive,
inner blobs are decompressed IN MEMORY and dropped (never expanded to scratch --
that is what filled /local/scratch in July 2026).

Idempotent: a design whose output files already exist is not requested again, so
a failed task can be re-run alone.
"""
import argparse
import collections
import csv
import io
import json
import os
import sys
import tarfile

import zstandard as zstd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = f"{RB}/campaigns/csy4_crrna_free"
BIF_DIR = "boltzgen__boltzif"
MPNN_DIR = "boltzgen__ligandmpnn"


def plan(rows, midrun, outdir, designs_dir, with_design):
    """{inner blob member -> [(wanted basename, output path)]}"""
    want = collections.defaultdict(list)
    for r in rows:
        u, did = r["uuid8"], r["design_id"]
        a_out = os.path.join(outdir, f"{did}_af3.cif")
        e_out = os.path.join(outdir, f"{did}_esm.cif")
        if not os.path.exists(a_out):
            want[f"{midrun}/refold_af3/msaoff_{u}/designs.tar.zst"].append(
                (f"{did}_AF3_seed1_sample0.cif", a_out))
        if not os.path.exists(e_out):
            want[f"{midrun}/refold_esm/msaoff_{u}/designs.tar.zst"].append(
                (f"{did}_ESM_sample0.cif", e_out))
        if with_design:
            d_out = os.path.join(designs_dir, f"{did}.cif")
            if not os.path.exists(d_out):
                want[f"{midrun}/designs/{u}.tar.zst"].append((f"{did}.cif", d_out))
    return want


def pull(blob, targets, log, blob_name):
    by_base = dict(targets)
    got = 0
    dctx = zstd.ZstdDecompressor()
    with dctx.stream_reader(io.BytesIO(blob)) as rd:
        with tarfile.open(fileobj=rd, mode="r|") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                out = by_base.get(os.path.basename(m.name))
                if out is None:
                    continue
                data = tf.extractfile(m).read()
                os.makedirs(os.path.dirname(out), exist_ok=True)
                with open(out, "wb") as fh:
                    fh.write(data)
                got += 1
                if got == len(by_base):
                    break
    for base, out in targets:
        if not os.path.exists(out):
            log.append(["MISS_MEMBER", blob_name, base])
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--designs-dir", default=None)
    ap.add_argument("--stack", choices=["bif", "limpnn"], default="bif")
    ap.add_argument("--with-design", action="store_true",
                    help="also pull the BIF design file (T2 correctness gate 3 only)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    designs_dir = a.designs_dir or os.path.join(os.path.dirname(a.outdir), "designs")
    rows = [r for r in csv.DictReader(open(a.manifest), delimiter="\t")
            if r["campaign"] == a.campaign and r["ref_kind"] != "unresolvable"]
    if not rows:
        print(f"no rows for campaign {a.campaign}")
        return
    os.makedirs(a.outdir, exist_ok=True)

    midrun = f"midrun_{a.campaign}"
    sub = BIF_DIR if a.stack == "bif" else MPNN_DIR
    outer = f"{CAMP}/{sub}/{midrun}.tar.zst"
    want = plan(rows, midrun, a.outdir, designs_dir, a.with_design)
    nfiles = sum(len(v) for v in want.values())
    print(f"[{a.campaign}/{a.stack}] outer={outer}")
    print(f"[{a.campaign}/{a.stack}] {len(rows)} designs, {len(want)} blobs, {nfiles} files to pull")
    if a.dry_run:
        for b in sorted(want):
            print(f"  {b}  ({len(want[b])} files)")
        return
    if nfiles == 0:
        print("nothing to do (already complete)")
        return

    log, seen = [], 0
    dctx = zstd.ZstdDecompressor()
    with open(outer, "rb") as fh, dctx.stream_reader(fh) as rd:
        with tarfile.open(fileobj=rd, mode="r|") as tf:
            for m in tf:
                if m.name not in want:
                    continue
                blob = tf.extractfile(m).read()
                n = pull(blob, want[m.name], log, m.name)
                del blob
                seen += 1
                print(f"  [{seen}/{len(want)}] {m.name} -> {n}", flush=True)
                if seen == len(want):
                    break   # stops before the run_*/ tail

    got = sum(1 for b in want for _, o in want[b] if os.path.exists(o))
    if seen != len(want):
        log.append(["MISS_BLOB", f"{len(want) - seen} inner blobs never seen", ""])
    print(f"[{a.campaign}/{a.stack}] pulled {got}/{nfiles}; blobs hit {seen}/{len(want)}")
    with open(os.path.join(a.outdir, f"_log_{a.stack}_{a.campaign}.json"), "w") as fh:
        json.dump({"outer": outer, "blobs_wanted": len(want), "blobs_seen": seen,
                   "files_wanted": nfiles, "files_got": got, "log": log}, fh, indent=1)
    if got != nfiles:
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Dry-run to verify the member paths before spending any read**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python extract_refolds.py --manifest manifest.tsv --campaign v6 --outdir cifs --dry-run | head -5'
```

Expected: `133 designs, 266 blobs, 266 files to pull` (blob count may be lower if uuids repeat).

- [ ] **Step 3: Cross-check every planned blob against the archive's `.list` sidecar**

This is the cheap pre-flight that catches a wrong member path without streaming
174 GB. It found 0 missing on the 20-pair probe.

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python extract_refolds.py --manifest manifest.tsv --campaign v6 --outdir cifs --dry-run | grep "^  midrun" | awk "{print \$1}" | sort -u > /dev/shm/w.txt; comm -23 /dev/shm/w.txt <(sort -u /cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/boltzgen__boltzif/midrun_v6.tar.zst.list) | wc -l'
```

Expected: `0`.

- [ ] **Step 4: Write `extract.sbatch`**

```bash
#!/bin/bash
#SBATCH --job-name=bif_rescore_extract
#SBATCH --array=0-7
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore/logs/extract_%A_%a.log
set -euo pipefail

W=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python

# task 7 is the LiMPNN v2 pass that supplies backbones for the 107 fallbacks
CAMPS=(v2 v4 v5 v6 v7 v8 v9 v2)
STACKS=(bif bif bif bif bif bif bif limpnn)
C=${CAMPS[$SLURM_ARRAY_TASK_ID]}
S=${STACKS[$SLURM_ARRAY_TASK_ID]}

OUT=$W/cifs
EXTRA=""
if [ "$S" = "limpnn" ]; then OUT=$W/limpnn_designs; fi
if [ "$C" = "v6" ] && [ "$S" = "bif" ]; then EXTRA="--with-design"; fi

echo "task $SLURM_ARRAY_TASK_ID -> $C/$S  node=$(hostname)  start=$(date)"
"$PY" "$W/extract_refolds.py" --manifest "$W/manifest.tsv" \
    --campaign "$C" --stack "$S" --outdir "$OUT" $EXTRA
echo "done=$(date)"
```

Note: the LiMPNN v2 task needs a different member pattern (it pulls
`designs/<uuid>.tar.zst` members named `csy4free_BG_LiMPNN_...`). Task 10 handles
that; this array entry exists so the archive is streamed once alongside the rest.

- [ ] **Step 5: Checkpoint**

Append: `Task 4 done — extractor written, v6 dry-run + .list pre-flight clean (0 missing).`
Do **not** submit the array yet; Task 9 (T2) runs v6 alone first.

---

## Task 5: `build_refs.py` — transplant, with the sequence gate

**Files:**
- Create: `$W/build_refs.py`
- Output: `$W/refs/<design_id>.pdb`

**Interfaces:**
- Consumes: `rescore_lib.{split_prot_rna, prot_backbone, transplant, seq_one_letter}`, `backbones/<uuid8>_spec<D>.cif`, `cifs/<design_id>_af3.cif`.
- Produces: `refs/<design_id>.pdb` — chain A protein backbone, chain B RNA.

**Correctness gate 1 lives here.** The transplant takes its sequence from the
library CSV; the AF3 refold was folded from the design's own sequence. If those
disagree the design is mislabelled, and BG+LiMPNN sequences were mislabelled as
BG+BIF in **v2** once before — v2 is 235 of these 2,400.

- [ ] **Step 1: Write `build_refs.py`**

```python
#!/usr/bin/env python3
"""Build the corrected reference for each design: the TRUE BoltzGen backbone
carrying the design's own sequence.

ref_kind == 'spec'            -> backbones/<uuid>_spec<D>.cif
ref_kind == 'limpnn_fallback' -> limpnn_designs/<limpnn_design_id>.cif
        (exact, not approximate: a BG+LiMPNN design file IS the spec backbone)

The spec CANNOT be used directly -- it carries BoltzGen's own co-designed
sequence, ~50% identical to the BIF sequence, which trips lrmsd_rna_aligned.py's
MIN_REGISTER=0.9 name-agreement gate and makes DockQ's alignment unreliable.
"""
import argparse
import collections
import csv
import glob
import os
import sys

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescore_lib import prot_backbone, seq_one_letter, split_prot_rna, transplant


def load(path):
    return get_structure(CIFFile.read(path), model=1)


def write_pdb(prot, rna, path):
    prot, rna = prot.copy(), rna.copy()
    prot.chain_id[:] = "A"
    rna.chain_id[:] = "B"
    st = prot + rna
    st.b_factor = np.zeros(len(st), dtype=float)
    st.occupancy = np.ones(len(st), dtype=float)
    f = PDBFile()
    f.set_structure(st)
    f.write(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--campaign", default=None, help="restrict to one campaign")
    a = ap.parse_args()
    R = a.root
    os.makedirs(f"{R}/refs", exist_ok=True)

    rows = [r for r in csv.DictReader(open(a.manifest), delimiter="\t")
            if r["ref_kind"] != "unresolvable"]
    if a.campaign:
        rows = [r for r in rows if r["campaign"] == a.campaign]

    n_ok = 0
    seq_bad, missing = [], []
    for r in rows:
        did, u, d = r["design_id"], r["uuid8"], int(r["design_idx"])
        out = f"{R}/refs/{did}.pdb"
        if r["ref_kind"] == "spec":
            src = f"{R}/backbones/{u}_spec{d}.cif"
        else:
            g = glob.glob(f"{R}/limpnn_designs/csy4free_BG_LiMPNN_len80_{u}_design{d}_seq*.cif")
            src = sorted(g)[0] if g else None
        af3 = f"{R}/cifs/{did}_af3.cif"
        if not src or not os.path.exists(src) or not os.path.exists(af3):
            missing.append(did)
            continue

        sp, sr = split_prot_rna(load(src))
        bb = prot_backbone(sp)

        # --- correctness gate 1: the reference must be the same molecule ---
        ap_, _ = split_prot_rna(load(af3))
        refold_seq = seq_one_letter(ap_)
        if refold_seq != r["sequence"]:
            seq_bad.append((did, r["campaign"], r["sequence"], refold_seq))
            continue

        write_pdb(transplant(bb, r["sequence"]), sr, out)
        n_ok += 1

    print(f"refs written: {n_ok}/{len(rows)}")
    if missing:
        print(f"MISSING inputs: {len(missing)} (first 5: {missing[:5]})")
    if seq_bad:
        print(f"\nSEQUENCE MISMATCH -- manifest vs AF3 refold: {len(seq_bad)}")
        print(f"  by campaign: {dict(collections.Counter(x[1] for x in seq_bad))}")
        for did, camp, a_, b_ in seq_bad[:3]:
            print(f"  {did} ({camp})\n    manifest {a_}\n    refold   {b_}")
        sys.exit(2)
    if missing:
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Regression-test the transplant against a known-good file**

The recipe was validated on the 20-pair probe by reproducing real BG+LiMPNN
design files. Re-assert that here, using the probe's bundle, before trusting it
at scale:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_backbone_probe && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python build_refs.py --bundle bundle --refs /dev/shm/refs_check 2>&1 | tail -3'
```

Expected: `recipe validation: 20/20 pairs reproduce the real BG+LiMPNN design file exactly`.

- [ ] **Step 3: Checkpoint**

Append: `Task 5 done — build_refs.py written; transplant recipe re-validated 20/20.`

---

## Task 6: `score.py` — drive the campaign's scorers and recombine

**Files:**
- Create: `$W/score.py`
- Output: `$W/score_in/*.csv`, `$W/score_out/*.csv`

**Interfaces:**
- Consumes: `refs/`, `cifs/`, `rescore_lib.{dockq_rna, lrmsd_rna_corrected}`.
- Produces: `score_out/scored_<campaign>.csv` with columns
  `design_id, backend, fnat, irmsd, lrmsd_std, lrmsd_rna, lrmsd_prot, sc_rmsd, dockq`.

`score.py` shells out to the two campaign scripts rather than reimplementing
them — that is what makes the numbers comparable to the stored ones.

- [ ] **Step 1: Write `score.py`**

```python
#!/usr/bin/env python3
"""Score the corrected references with the campaign's OWN scorers, then
recombine into the binder-campaign frame.

dockq_self.py       -> fnat, iRMSD, LRMSD(stock, binder-aligned), DockQ(stock)
lrmsd_rna_aligned.py-> lrmsd_rna, lrmsd_prot  (the ratio that re-frames LRMSD)
recombination        -> sc_RMSD (= dockq_LRMSD in the master) and DockQ_rna
"""
import argparse
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescore_lib import dockq_rna, lrmsd_rna_corrected

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
PY_BIO = "/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python"
BE = ("af3", "esm")


def write_inputs(rows, R, tag, ref_dir, ref_ext):
    os.makedirs(f"{R}/score_in", exist_ok=True)
    lr, dq = [], []
    for r in rows:
        did = r["design_id"]
        ref = f"{R}/{ref_dir}/{did}{ref_ext}"
        a, e = f"{R}/cifs/{did}_af3.cif", f"{R}/cifs/{did}_esm.cif"
        if not (os.path.exists(ref) and os.path.exists(a) and os.path.exists(e)):
            continue
        lr.append({"design_id": did, "designed_cif_path": ref,
                   "af3_msaoff_cif_path": a, "esm_cif_path": e, "b2_cif_path": ""})
        for be, p in (("af3", a), ("esm", e)):
            dq.append({"design_id": did, "backend": be, "abs_path": p,
                       "designed_source": ref})
    if not lr:
        sys.exit("no design has all of refs/, cifs/_af3, cifs/_esm -- nothing to score")
    for rows_, name in ((lr, f"lrmsd_in_{tag}.csv"), (dq, f"dockq_in_{tag}.csv")):
        with open(f"{R}/score_in/{name}", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_[0]))
            w.writeheader()
            w.writerows(rows_)
    print(f"score_in: {len(lr)} designs, {len(dq)} design x backend")
    return len(lr)


def run(R, tag, jobs):
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{RB}/code/pylibs_dockq:{RB}/code/scratch_dockq"
    env["PATH"] = f"{RB}/code/bin:" + env.get("PATH", "")
    subprocess.run(["python3", f"{RB}/code/scripts/dockq_self.py",
                    "--csv", f"{R}/score_in/dockq_in_{tag}.csv",
                    "--output", f"{R}/score_out/dockq_{tag}.csv",
                    "--path-col", "abs_path", "--ref-col", "designed_source",
                    "--jobs", str(jobs), "--progress"], check=True, env=env)
    env2 = dict(os.environ)
    env2["PYTHONPATH"] = f"{RB}/code/src:{RB}/code/scripts"
    env2["PATH"] = f"{RB}/code/bin:" + env2.get("PATH", "")
    subprocess.run([PY_BIO, f"{RB}/code/scripts/lrmsd_rna_aligned.py",
                    "--csv", f"{R}/score_in/lrmsd_in_{tag}.csv",
                    "--output", f"{R}/score_out/lrmsd_{tag}.csv",
                    "--jobs", str(jobs)], check=True, env=env2)


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


def recombine(R, tag):
    dq = {}
    for r in csv.DictReader(open(f"{R}/score_out/dockq_{tag}.csv")):
        dq[(r["design_id"], r["backend"])] = r
    lr = {r["design_id"]: r for r in csv.DictReader(open(f"{R}/score_out/lrmsd_{tag}.csv"))}

    out, bad = [], []
    for (did, be), d in sorted(dq.items()):
        l = lr.get(did)
        fnat, irms = fnum(d["dockq_fnat"]), fnum(d["dockq_iRMSD"])
        stock = fnum(d["dockq_LRMSD"])
        rna = fnum(l[f"lrmsd_rna_{be}"]) if l else None
        prot = fnum(l[f"lrmsd_prot_{be}"]) if l else None
        sc = lrmsd_rna_corrected(stock, rna, prot) if None not in (stock, rna, prot) else None
        if sc is None or fnat is None or irms is None:
            bad.append((did, be))          # correctness gate 2
            continue
        out.append({"design_id": did, "backend": be, "fnat": f"{fnat:.4f}",
                    "irmsd": f"{irms:.3f}", "lrmsd_std": f"{stock:.3f}",
                    "lrmsd_rna": f"{rna:.3f}", "lrmsd_prot": f"{prot:.3f}",
                    "sc_rmsd": f"{sc:.4f}", "dockq": f"{dockq_rna(fnat, irms, sc):.4f}"})

    p = f"{R}/score_out/scored_{tag}.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {p}: {len(out)} rows")
    if bad:
        print(f"REGISTER/SCORE FAILURES (gate 2): {len(bad)}")
        for x in bad[:10]:
            print(f"  {x[0]} {x[1]}")
    return len(out), len(bad)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--campaign", default=None)
    ap.add_argument("--tag", default=None,
                    help="names the output score_out/scored_<tag>.csv. REQUIRED when "
                         "several runs share a manifest scope, otherwise they overwrite "
                         "each other (T1/T2b/T3 all cover 'all').")
    ap.add_argument("--ref-dir", default="refs",
                    help="directory of reference structures, relative to --root")
    ap.add_argument("--ref-ext", default=".pdb")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    tag = a.tag or a.campaign
    if not tag:
        sys.exit("--tag is required when --campaign is not given")
    os.makedirs(f"{a.root}/score_out", exist_ok=True)
    rows = [r for r in csv.DictReader(open(a.manifest), delimiter="\t")
            if r["ref_kind"] != "unresolvable"]
    if a.campaign:
        rows = [r for r in rows if r["campaign"] == a.campaign]
    write_inputs(rows, a.root, tag, a.ref_dir, a.ref_ext)
    run(a.root, tag, a.jobs)
    n, bad = recombine(a.root, tag)
    if bad:
        sys.exit(2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Checkpoint**

Append: `Task 6 done — score.py written (not yet run; Task 9 runs it on v6).`

---

## Task 7: `report.py` — re-gate and emit the rate

**Files:**
- Create: `$W/report.py`
- Output: `$W/rescore_2400.csv`, `$W/REPORT.md`

**Interfaces:**
- Consumes: `manifest.tsv`, `score_out/scored_*.csv`, `rescore_lib.gate_pass`.
- Produces: the deliverable.

- [ ] **Step 1: Write `report.py`**

```python
#!/usr/bin/env python3
"""Re-gate the rescored designs and emit the flip rate."""
import argparse
import collections
import csv
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescore_lib import gate_pass

COLS = ["design_id", "campaign", "ref_kind",
        "camp_sc_rmsd_af3", "camp_sc_rmsd_esm", "camp_dockq_af3", "camp_dockq_esm",
        "true_sc_rmsd_af3", "true_sc_rmsd_esm", "true_dockq_af3", "true_dockq_esm",
        "camp_gate", "true_gate", "flip"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--tags", required=True,
                    help="comma-separated score_out tags to read, e.g. 'final'. "
                         "Explicit rather than a glob: T1/T2b/gate-3 runs also live in "
                         "score_out/ and must never be mixed into the headline result.")
    a = ap.parse_args()
    R = a.root

    scored = collections.defaultdict(dict)
    for t in a.tags.split(","):
        p = f"{R}/score_out/scored_{t.strip()}.csv"
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        for r in csv.DictReader(open(p)):
            scored[r["design_id"]][r["backend"]] = r

    man = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    rows, n_unres, n_fail = [], 0, 0
    for m in man:
        did = m["design_id"]
        if m["ref_kind"] == "unresolvable":
            n_unres += 1
            continue
        s = scored.get(did, {})
        if "af3" not in s or "esm" not in s:
            n_fail += 1
            continue
        c = [float(m[f"camp_{k}"]) for k in
             ("sc_rmsd_af3", "sc_rmsd_esm", "dockq_af3", "dockq_esm")]
        t = [float(s["af3"]["sc_rmsd"]), float(s["esm"]["sc_rmsd"]),
             float(s["af3"]["dockq"]), float(s["esm"]["dockq"])]
        cg, tg = gate_pass(*c), gate_pass(*t)
        rows.append(dict(zip(COLS, [did, m["campaign"], m["ref_kind"],
                                    *[f"{x:.4f}" for x in c], *[f"{x:.4f}" for x in t],
                                    int(cg), int(tg), int(cg and not tg)])))

    with open(f"{R}/rescore_2400.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    flips = sum(int(r["flip"]) for r in rows)
    camp_pass = sum(int(r["camp_gate"]) for r in rows)
    sc_c = [float(r["camp_sc_rmsd_af3"]) for r in rows] + [float(r["camp_sc_rmsd_esm"]) for r in rows]
    sc_t = [float(r["true_sc_rmsd_af3"]) for r in rows] + [float(r["true_sc_rmsd_esm"]) for r in rows]
    dq_c = [float(r["camp_dockq_af3"]) for r in rows] + [float(r["camp_dockq_esm"]) for r in rows]
    dq_t = [float(r["true_dockq_af3"]) for r in rows] + [float(r["true_dockq_esm"]) for r in rows]

    # correctness gate 4: the accounting must close
    total = n + n_unres + n_fail
    assert total == len(man), f"accounting: {n}+{n_unres}+{n_fail} != {len(man)}"

    by_c = collections.Counter(r["campaign"] for r in rows if int(r["flip"]))
    tot_c = collections.Counter(r["campaign"] for r in rows)

    lines = [
        "# BG+BIF library rescore — result", "",
        f"Scored **{n}** of {len(man)} BG+BIF designs in the delivered 12k library "
        f"against the TRUE BoltzGen backbone.", "",
        "## Accounting", "",
        f"- scored: **{n}**", f"- unresolvable (no backbone, no LiMPNN twin): **{n_unres}**",
        f"- failed to score: **{n_fail}**", f"- total: {total} = {len(man)}", "",
        "## Result", "",
        f"- passed the gate under the campaign's reference: **{camp_pass}/{n}**",
        f"- still pass against the true backbone: **{camp_pass - flips}/{n}**",
        f"- **flipped PASS -> FAIL: {flips}/{camp_pass} "
        f"({100.0 * flips / camp_pass:.1f}% of those that passed)**", "",
        "| metric | campaign | true backbone |", "|---|---|---|",
        f"| sc_RMSD median | {st.median(sc_c):.3f} A | {st.median(sc_t):.3f} A |",
        f"| DockQ median | {st.median(dq_c):.3f} | {st.median(dq_t):.3f} |", "",
        "## Flips by campaign", "", "| campaign | flipped | scored |", "|---|---|---|",
    ]
    for c in sorted(tot_c):
        lines.append(f"| {c} | {by_c.get(c, 0)} | {tot_c[c]} |")
    lines += [
        "", "## What this does NOT cover", "",
        "- BG+BIF designs **outside** the 12k library. This rate is conditioned on "
        "having been selected under the wrong reference and must not be quoted as "
        "a stack-wide figure.",
        "- The other four stacks. Only BG+BIF stores a re-prediction as its design, "
        "so they are unaffected.", "",
        "Per-design table: `rescore_2400.csv`.",
    ]
    open(f"{R}/REPORT.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Checkpoint**

Append: `Task 7 done — report.py written.`

---

## Task 8: T1 — regression against the 20-pair probe

**Files:** none created. Uses `$W/*` against the probe's known-good numbers.

**Purpose:** prove the rewritten pipeline reproduces the probe before it is
trusted on new data. The probe's R3 numbers are in
`$RB/analyses/bif_backbone_probe/rescore_comparison.csv` (`R3_sc_rmsd`, `R3_dockq`).

- [ ] **Step 1: Build a 20-row manifest from the probe pairs**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && head -1 manifest.tsv > t1_manifest.tsv && for d in $(cut -f5 /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_backbone_probe/pairs_manifest.tsv | tail -n +2); do grep -P "^$d\t" manifest.tsv >> t1_manifest.tsv || true; done && wc -l < t1_manifest.tsv'
```

Expected: `21` (header + 20). If fewer, some probe designs are not in the 2,400 — list them before continuing.

- [ ] **Step 2: Run the full pipeline on those 20**

The probe already extracted their refolds; symlink rather than re-stream:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && mkdir -p cifs && for p in /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_backbone_probe/bundle/pair*/; do for f in $p/02_BIF_af3_*.cif; do b=$(basename $f); d=${b#02_BIF_af3_}; ln -sf "$f" "cifs/${d%.cif}_af3.cif"; done; for f in $p/03_BIF_esm_*.cif; do b=$(basename $f); d=${b#03_BIF_esm_}; ln -sf "$f" "cifs/${d%.cif}_esm.cif"; done; done && ls cifs | wc -l'
```

Expected: `40`.

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python && $PY extract_backbones.py --manifest t1_manifest.tsv --outdir backbones 2>&1 | tail -2; $PY build_refs.py --manifest t1_manifest.tsv --root . && $PY score.py --manifest t1_manifest.tsv --root . --tag t1 --jobs 4'
```

Expected: `refs written: 20/20`, then `wrote .../scored_t1.csv: 40 rows`, no gate-2 failures.

- [ ] **Step 3: Compare to the probe's stored R3 numbers**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python - <<"PYEOF"
import csv
new = {(r["design_id"], r["backend"]): r for r in csv.DictReader(open("score_out/scored_t1.csv"))}
old = list(csv.DictReader(open("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_backbone_probe/rescore_comparison.csv")))
ds = dd = 0.0
n = 0
for r in old:
    k = (r["design_id"], r["backend"])
    if k not in new:
        continue
    ds = max(ds, abs(float(new[k]["sc_rmsd"]) - float(r["R3_sc_rmsd"])))
    dd = max(dd, abs(float(new[k]["dockq"]) - float(r["R3_dockq"])))
    n += 1
print(f"compared {n}; max|d sc_RMSD| = {ds:.6f}; max|d DockQ| = {dd:.6f}")
assert n == 40, n
assert ds < 1e-3 and dd < 1e-3, "T1 REGRESSION FAILED"
print("T1 PASS")
PYEOF'
```

Expected: `compared 40; max|d sc_RMSD| = 0.000000; max|d DockQ| = 0.000000` then `T1 PASS`.
A non-zero difference means the rewrite changed behaviour — **stop and diff against the probe scripts** before going further.

- [ ] **Step 4: Checkpoint**

Append: `Task 8 done — T1 regression PASS, 40/40 identical to the probe.`

---

## Task 9: T2 — v6 end-to-end, and the first real timing

**Purpose:** the first full archive stream. Measures what the probe's 26–37 s
could not (it broke early on low-hex uuids), and exercises correctness gate 3.

- [ ] **Step 1: Extract v6 alone, timed**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && mkdir -p logs && sbatch --array=3 --job-name=rescore_v6 extract.sbatch'
```

(Array index 3 is v6 in the `CAMPS` list; it carries `--with-design`.)

- [ ] **Step 2: Confirm completion and record the elapsed time**

```bash
ssh euler 'sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS -P -n | grep -v extern'
```

Expected: `COMPLETED`. **Record `Elapsed`** — this is the number that sets the T3
estimate. v7 is 4× v6's compressed size, so expect roughly 4× for the long pole.
If v6 takes more than ~30 min, re-plan T3 (raise `--time`, consider splitting).

- [ ] **Step 3: Verify extraction completeness**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && ls cifs/*_af3.cif | wc -l && ls cifs/*_esm.cif | wc -l && ls designs/*.cif 2>/dev/null | wc -l'
```

Expected: 133 v6 designs each for af3/esm (plus the 20 T1 symlinks — subtract them), and 133 design files.

- [ ] **Step 4: Build refs and score v6**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python && $PY extract_backbones.py --manifest manifest.tsv --outdir backbones 2>&1 | tail -2; $PY build_refs.py --manifest manifest.tsv --root . --campaign v6 && $PY score.py --manifest manifest.tsv --root . --campaign v6 --jobs 8'
```

Expected: `refs written: 133/133`, `wrote .../scored_v6.csv: 266 rows`, no gate-2 failures, **no sequence mismatch** (gate 1).

- [ ] **Step 5: Correctness gate 3 — reproduce the campaign numbers on v6**

Score v6 again, but against the **BIF design files** (`--ref-dir designs`) instead
of the corrected references. This is the same code path, so it is one command:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python score.py --manifest manifest.tsv --root . --campaign v6 --tag v6repro --ref-dir designs --ref-ext .cif --jobs 8'
```

Expected: `wrote .../scored_v6repro.csv: 266 rows`.

Then compare those to the stored campaign values:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python - <<"PYEOF"
import csv
man = {r["design_id"]: r for r in csv.DictReader(open("manifest.tsv"), delimiter="\t")}
ds = dd = 0.0
n = 0
for r in csv.DictReader(open("score_out/scored_v6repro.csv")):
    m = man[r["design_id"]]
    be = r["backend"]
    ds = max(ds, abs(float(r["sc_rmsd"]) - float(m[f"camp_sc_rmsd_{be}"])))
    dd = max(dd, abs(float(r["dockq"]) - float(m[f"camp_dockq_{be}"])))
    n += 1
print(f"compared {n}; max|d sc_RMSD| = {ds:.6f}; max|d DockQ| = {dd:.6f}")
assert n == 266, n
assert ds < 0.01 and dd < 0.01, "GATE 3 FAILED -- the library CSV does not reproduce"
print("GATE 3 PASS")
PYEOF'
```

Expected: `compared 266`, both deltas < 0.01 (the probe achieved 0.0005 / 0.0000),
then `GATE 3 PASS`. If it does not reproduce, **stop** — the "trust the library
CSV" assumption behind approach A is invalid and the plan needs revisiting.

- [ ] **Step 6: Checkpoint**

Append: `Task 9 done — T2 v6: extract <ELAPSED>, 133/133 refs, 266 scored, gates 1/2/3 pass.`

---

## Task 10: T2b — the 107 LiMPNN-fallback designs

**Purpose:** the fallback is a separate code path and is tested separately.

- [ ] **Step 1: Extract BIF v2 refolds and LiMPNN v2 designs**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && sbatch --array=0,7 extract.sbatch'
```

Array 0 = BIF v2 (refolds for all 235 v2 designs), array 7 = LiMPNN v2
(design files that serve as backbones for the 107 fallbacks).

- [ ] **Step 2: Adapt the LiMPNN pull**

`extract_refolds.py` as written pulls `refold_*` members. For `--stack limpnn`
the wanted member is the **design** instead. Add to `plan()`:

```python
    if stack == "limpnn":
        want = collections.defaultdict(list)
        for r in rows:
            if r["ref_kind"] != "limpnn_fallback":
                continue
            u, d = r["uuid8"], int(r["design_idx"])
            want[f"{midrun}/designs/{u}.tar.zst"].append(
                (None, os.path.join(outdir, f"{u}_design{d}")))
        return want
```

and in `pull()`, when the wanted basename is `None`, match any member matching
`csy4free_BG_LiMPNN_len80_<uuid>_design<D>_seq*.cif` and write it under its own
name into `limpnn_designs/`. Keep the first seq index encountered — all seqs on
one backbone share coordinates (LigandMPNN moves no atoms), so any of them gives
the same reference.

- [ ] **Step 3: Build refs and score the fallback designs only**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && head -1 manifest.tsv > t2b_manifest.tsv && awk -F"\t" "\$6==\"limpnn_fallback\"" manifest.tsv >> t2b_manifest.tsv && wc -l < t2b_manifest.tsv && PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python && $PY build_refs.py --manifest t2b_manifest.tsv --root . && $PY score.py --manifest t2b_manifest.tsv --root . --tag t2b --jobs 8'
```

Expected: `108` lines (header + 107), `refs written: 107/107`, `wrote .../scored_t2b.csv: 214 rows`.

- [ ] **Step 4: Cross-validate the fallback against the authoritative source**

The fallback asserts that a BG+LiMPNN design file is the `spec_<D>` backbone. That
was validated on the 20 probe pairs; validate it here on real v2 data from this
run, at zero extra I/O — the LiMPNN v2 archive is already streamed for Step 1.

**v2 has 235 library designs, of which 108 lack a spec backbone — so 127 have
BOTH** a `spec_<D>.cif` in `backbones/` and a LiMPNN counterpart. Pull the LiMPNN
designs for those 127 as well (they come from the same stream), build a reference
from each, and compare to the `spec`-derived reference.

Criterion: protein backbone atoms must agree to **< 0.01 Å** and residue names
must match 100%, matching the 0.0005 Å / 320-of-320 seen on the probe. Any design
exceeding that invalidates the fallback for the other 107 — **stop and report**
rather than proceeding.

If the count of designs with both sources turns out to be 0, say so explicitly and
skip; do not silently pass a check that never ran.

- [ ] **Step 5: Checkpoint**

Append: `Task 10 done — T2b: 107/107 fallback refs built and scored.`

---

## Task 11: T3 — the remaining archives and the final report

- [ ] **Step 1: Submit the remaining campaigns**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && sbatch --array=1,2,4,5,6 extract.sbatch'
```

(v4, v5, v7, v8, v9 — v2 and v6 are already done.)

- [ ] **Step 2: Confirm all tasks COMPLETED and count the working set**

```bash
ssh euler 'sacct -j <JOBID> --format=JobID,State,Elapsed -P -n | grep -v extern; cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && ls cifs/*_af3.cif | wc -l && ls cifs/*_esm.cif | wc -l && du -sh cifs'
```

Expected: 2,399 of each (2,400 minus the 1 unresolvable), ~0.45 GB.

- [ ] **Step 3: Build all refs and score everything**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python && $PY extract_backbones.py --manifest manifest.tsv --outdir backbones --expect 2292 && $PY build_refs.py --manifest manifest.tsv --root . && $PY score.py --manifest manifest.tsv --root . --tag final --jobs 8'
```

Expected: `refs written: 2399/2399`, `wrote .../scored_final.csv: 4798 rows`, no gate-1 or gate-2 failures.

- [ ] **Step 4: Produce the report**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore && /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python report.py --manifest manifest.tsv --root . --tags final'
```

Expected: `REPORT.md` printed, with the accounting assertion passing
(`2399 + 1 + 0 == 2400`).

- [ ] **Step 5: Copy the deliverables off Euler**

```bash
rsync -a euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/bif_library_rescore/{REPORT.md,rescore_2400.csv,manifest.tsv} ~/Downloads/bif_library_rescore/
```

- [ ] **Step 6: Final checkpoint**

Append: `Task 11 done — full run complete, REPORT.md produced.`
Update the memory `rna-binder-bgbif-design-is-reprediction.md` with the measured
rate, replacing the 20-pair figure as the headline, and run
`bash ~/sync_claude_memory.sh`.

---

## Self-Review Notes

- **Spec coverage.** §1–§4 → Tasks 1, 5, 6. §3 facts → Task 2's assertions. §5 architecture → Tasks 3, 4. §5 gate → Task 1 `gate_pass` + Task 7. §6 gates: 1 → Task 5 Step 1; 2 → Task 6 `recombine`; 3 → Task 9 Step 5; 4 → Task 7 assertion; 5 → Task 1 `split_prot_rna`. §7 ladder → Tasks 8–11. §9 success criteria → Task 7 output.
- **Known rough edge, flagged deliberately:** Task 10 Step 2 modifies `extract_refolds.py` after Task 4 wrote it. That is real coupling; it is placed there rather than in Task 4 so the LiMPNN path is written and tested against actual v2 data rather than speculatively.
