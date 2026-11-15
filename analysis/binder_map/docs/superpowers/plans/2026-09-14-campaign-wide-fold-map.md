# Campaign-Wide Fold Map Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce three fold clusterings under identical conditions differing only in sample size — a uniform 500,000-node random sample of the whole campaign, the existing 31,763-node filtered map, and the 9,759-node map of the ordered 12k library.

**Architecture:** Reuse `foldseek_map` stages 2–5 unchanged, run once per map. Add a sampler (stage 0) that draws two node sets and rewrites design refs for the 2026-09-11 nested archives; rewrite stage 1 to stream one big archive per task over the UNION of both node sets, so 683 GB is read once; split the shards per map; fix two allocation sites that break above ~100k nodes; add a stage that compares the three finished maps.

**Tech Stack:** Python 3.11 + numpy/pandas/pyarrow/biotite, Foldseek v8-ef4e960, leidenalg, Slurm on Euler, zstd/tar for archive streaming.

## Global Constraints

- **Author locally, test on Euler, commit locally.** The laptop has NO pytest and
  no pyarrow; Euler's venv has both plus biotite and zstd. So the loop is:
  1. Write/edit files in the local repo `~/Documents/binder_map` (normal tools).
  2. `bash sync.sh` — pushes src/, scripts/, tests/, slurm/ to Euler (already exists).
  3. Run the tests over ssh on Euler.
  4. `git add <explicit paths> && git commit` in the local repo.
  Never `git add -A`: the repo carries large untracked result trees.
- **Euler paths**: `FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map`,
  `PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python`.
- **Run tests as** (single command from the laptop):
  `ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/<file> -v"`
- The two trees are identical for every file this plan touches (verified 2026-09-14
  by md5: 36/38 match; the 2 that differ, `stage8_build_html.py` and
  `test_embedding.py`, are out of scope).
- **Clustering parameters are frozen** and must not be changed anywhere:
  `--alignment-type 1`, `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a`,
  `TM = min(qTM, tTM)`, `--tm-floor 0.5`, `--epi-floor 0.6`, `--k 30`, `--mutual`,
  Leiden fixed seed, **operative resolution `cluster_r80`**.
- **Sampling seed is 42**; map-1 size 500,000; universe 2,952,126 coord-hash
  nodes; sampling fraction 0.169372. Map 3 is a census of 9,759 nodes, not a draw.
- **Map 2 is never rebuilt.** `results_sw_k30_f50_mut` is the reference.
- **The dedup key over-splits BG+BIF ~9.8x and HD ~10x** (rigid-body offset and
  0.022 A packer jitter respectively — measured, neither is a fold difference).
  This is KEPT, because map 2 has it too and matching it is what makes the three
  maps comparable. Do not "fix" it.
- **Never write loose per-structure files** under `results/` — the project has a
  hard inode cap. Everything streams into a per-task `.tar.gz` shard.
- **Never expand a nested archive to scratch.** Decompress in memory and drop it.
  Unbounded scratch expansion filled `/local/scratch` to ~800 GB in July 2026 and
  produced blank-but-correct-length outputs that defeated the resume check.
- **Boltz-2 must stay absent.** No b2 column may enter any table.

---

### Task 1: Backbone universe rule

The design→backbone mapping is the foundation of the node set. It is per-stack and
non-obvious, so it gets its own module and its own tests.

**Files:**
- Create: `src/binder_map/universe.py`
- Test: `tests/test_universe.py`

**Interfaces:**
- Consumes: nothing
- Produces: `backbone_key(design_id: str, stack: str) -> str`,
  `ONE_BACKBONE_PER_DESIGN: frozenset[str]`

- [ ] **Step 1: Write the failing test**

Create `tests/test_universe.py`:

```python
from binder_map.universe import ONE_BACKBONE_PER_DESIGN, backbone_key


def test_limpnn_stacks_strip_the_seq_suffix():
    assert backbone_key(
        "csy4free_RFD3_LiMPNN_len80_a546ac7c_design31_seq0", "RFD3+LiMPNN"
    ) == "csy4free_RFD3_LiMPNN_len80_a546ac7c_design31"


def test_seq_variants_of_one_backbone_collapse_to_one_key():
    ids = [f"csy4free_BG_LiMPNN_len80_b2f7f57e_design13_seq{i}" for i in range(10)]
    assert len({backbone_key(i, "BG+LiMPNN") for i in ids}) == 1


def test_bd1_also_collapses():
    ids = [f"csy4free_BD1_LiMPNN_len80_000158e9_design10_seq{i}" for i in range(3)]
    assert len({backbone_key(i, "BD1+LiMPNN") for i in ids}) == 1


def test_bg_bif_is_one_backbone_per_design():
    # BoltzIF is generative all-atom: 1,629 of 6,082 name-prefixes in the v6
    # ground truth split into multiple coordinate hashes.
    ids = [f"csy4free_BG_BIF_len80_0000df81_design0_seq{i}" for i in range(4)]
    assert len({backbone_key(i, "BG+BIF") for i in ids}) == 4


def test_hd_is_one_backbone_per_design():
    # All 802 HD designs in the v6 pool carry distinct coord hashes across
    # 447 (uuid, recycle) groups -- sample0..9 are not one backbone.
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample1_pack1"
    assert backbone_key(a, "HD+LiMPNN") != backbone_key(b, "HD+LiMPNN")


def test_key_is_idempotent():
    key = backbone_key("csy4free_RFD3_LiMPNN_len80_a546ac7c_design31_seq0", "RFD3+LiMPNN")
    assert backbone_key(key, "RFD3+LiMPNN") == key


def test_one_backbone_stacks_are_exactly_these_two():
    assert ONE_BACKBONE_PER_DESIGN == frozenset({"BG+BIF", "HD+LiMPNN"})


def test_unknown_stack_raises():
    import pytest
    with pytest.raises(KeyError):
        backbone_key("whatever_seq0", "NOT+A+STACK")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_universe.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'binder_map.universe'`

- [ ] **Step 3: Write the implementation**

Create `src/binder_map/universe.py`:

```python
"""Design -> backbone mapping, per generation stack.

Foldseek clusters backbones, not designs, and the design->backbone relation is
different for each stack. Verified 2026-09-14 against the v6 coordinate-hash
ground truth (results_v6/design/dedup/members.tsv, 55,706 designs -> 31,763
backbones), counting how many name-prefixes split across multiple coord hashes:

    RFD3+LiMPNN   14,455 prefixes,     0 split
    BG+LiMPNN      5,694 prefixes,     0 split
    BD1+LiMPNN     1,774 prefixes,     0 split
    BG+BIF         6,082 prefixes, 1,629 SPLIT (27%)
    HD+LiMPNN        802 designs -> 802 distinct hashes across 447
                   (uuid, recycle) groups

LigandMPNN is inverse folding: it writes the input backbone with residue names
swapped and moves no atoms, so its _seqN variants are coordinate-identical.
BoltzIF is generative all-atom and moves the backbone, so BG+BIF _seqN variants
are different backbones. HD design ids carry no _seq at all.

This rule only decides WHICH MEMBER TO EXTRACT. The authoritative dedup remains
the stage-2 MD5 coordinate hash over the protein chain's N/CA/C/O atoms, which
also catches cross-uuid collisions the name rule cannot see. Measured residual:
applying this rule to the 55,706-design v6 pool yields 31,767 keys against 31,763
true backbones -- an over-count of exactly 4, the four BG+LiMPNN backbones that
span two name-prefixes.
"""
from __future__ import annotations

import re

#: Stacks where every design is its own backbone.
ONE_BACKBONE_PER_DESIGN = frozenset({"BG+BIF", "HD+LiMPNN"})

#: Stacks where LigandMPNN sequence variants share one backbone.
SEQ_VARIANT_STACKS = frozenset({"RFD3+LiMPNN", "BG+LiMPNN", "BD1+LiMPNN"})

_SEQ_SUFFIX = re.compile(r"_seq\d+$")


def backbone_key(design_id: str, stack: str) -> str:
    """Return the backbone this design sits on.

    Raises KeyError on an unrecognised stack rather than guessing: a silent
    wrong branch here would merge distinct backbones or split identical ones,
    and neither is visible downstream until the cluster sizes are wrong.
    """
    if stack in ONE_BACKBONE_PER_DESIGN:
        return design_id
    if stack in SEQ_VARIANT_STACKS:
        return _SEQ_SUFFIX.sub("", design_id)
    raise KeyError(f"unknown stack: {stack!r}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_universe.py -v`
Expected: PASS, 11 passed

- [ ] **Step 5: Verify the rule against real data**

Run:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import re, collections
from binder_map.universe import backbone_key
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
def stack_of(d):
    for tok, name in (("RFD3_LiMPNN","RFD3+LiMPNN"), ("BG_BIF","BG+BIF"),
                      ("BG_LiMPNN","BG+LiMPNN"), ("BD1_LiMPNN","BD1+LiMPNN"),
                      ("HD_LiMPNN","HD+LiMPNN")):
        if "_" + tok + "_" in d: return name
    raise KeyError(d)
keys, groups = set(), set()
with open(f"{FM}/results_v6/design/dedup/members.tsv") as f:
    f.readline()
    for line in f:
        canon, ids = line.rstrip("\n").split("\t")
        groups.add(canon)
        for d in ids.split(","):
            keys.add(backbone_key(d.strip(), stack_of(d.strip())))
print(f"rule keys {len(keys)}  true backbones {len(groups)}  over-count {len(keys)-len(groups)}")
assert len(groups) == 31763, len(groups)
assert len(keys) == 31767, len(keys)
print("PASS: over-count is exactly the 4 known cross-prefix backbones")
EOF
```

Expected: `rule keys 31767  true backbones 31763  over-count 4` then `PASS`

- [ ] **Step 6: Commit**

```bash
cd ~/Documents/binder_map
git add src/binder_map/universe.py tests/test_universe.py
git commit -m "Add per-stack design->backbone rule

BG+BIF and HD are one backbone per design; the three LiMPNN stacks dedup by
stripping _seq<N>. Verified against the v6 coord-hash ground truth: the rule
over-counts by exactly 4 on the 55,706-design pool."
```

---

### Task 2: Nested archive reader

The 2026-09-11 repack put design structures two levels deep inside solid zstd
streams. This module is the only place that knows that.

**Files:**
- Create: `src/binder_map/nested.py`
- Test: `tests/test_nested.py`

**Interfaces:**
- Consumes: nothing
- Produces: `nested_ref(design_ref: str) -> tuple[str, str, str]` returning
  `(outer_path, inner_path, member_name)`;
  `iter_members(outer: Path, wanted: dict[str, set[str]]) -> Iterator[tuple[str, str, bytes]]`
  yielding `(inner_name, member_name, blob)`

- [ ] **Step 1: Write the failing test**

Create `tests/test_nested.py`:

```python
import io
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest

from binder_map.nested import iter_members, nested_ref

pytestmark = pytest.mark.skipif(shutil.which("zstd") is None, reason="needs zstd")


def test_nested_ref_splits_a_store_design_ref():
    ref = ("campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3/designs/"
           "ab12cd34.tar.zst!csy4free_RFD3_LiMPNN_len80_ab12cd34_design0_seq0.cif")
    outer, inner, member = nested_ref(ref)
    assert outer == "campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3.tar.zst"
    assert inner == "midrun_v3/designs/ab12cd34.tar.zst"
    assert member == "csy4free_RFD3_LiMPNN_len80_ab12cd34_design0_seq0.cif"


def test_nested_ref_handles_a_suffixed_midrun():
    ref = ("campaigns/c/boltzgen__boltzif/midrun_v8_50e6da5c/designs/"
           "ff00.tar.zst!d.cif")
    outer, inner, member = nested_ref(ref)
    assert outer == "campaigns/c/boltzgen__boltzif/midrun_v8_50e6da5c.tar.zst"
    assert inner == "midrun_v8_50e6da5c/designs/ff00.tar.zst"
    assert member == "d.cif"


def _zstd(src: Path, dst: Path) -> None:
    subprocess.run(["zstd", "-q", "-f", str(src), "-o", str(dst)], check=True)


def _build_nested(tmp_path: Path) -> Path:
    """Build midrun_v9.tar.zst ! midrun_v9/designs/<u>.tar.zst ! <m>.cif"""
    work = tmp_path / "work"
    (work / "midrun_v9" / "designs").mkdir(parents=True)

    for uuid, members in (("aaaa", {"d0.cif": b"ZERO\n", "d1.cif": b"ONE\n"}),
                          ("bbbb", {"d2.cif": b"TWO\n"})):
        plain = work / f"{uuid}.tar"
        with tarfile.open(plain, "w") as tar:
            for name, data in members.items():
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
        _zstd(plain, work / "midrun_v9" / "designs" / f"{uuid}.tar.zst")
        plain.unlink()

    outer_plain = tmp_path / "midrun_v9.tar"
    with tarfile.open(outer_plain, "w") as tar:
        tar.add(work / "midrun_v9", arcname="midrun_v9")
    outer = tmp_path / "midrun_v9.tar.zst"
    _zstd(outer_plain, outer)
    return outer


def test_iter_members_pulls_only_what_was_asked_for(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"d1.cif"},
              "midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    got = {(inner, member): blob for inner, member, blob in iter_members(outer, wanted)}
    assert got == {
        ("midrun_v9/designs/aaaa.tar.zst", "d1.cif"): b"ONE\n",
        ("midrun_v9/designs/bbbb.tar.zst", "d2.cif"): b"TWO\n",
    }


def test_iter_members_ignores_inner_archives_not_requested(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/bbbb.tar.zst": {"d2.cif"}}
    inners = {inner for inner, _, _ in iter_members(outer, wanted)}
    assert inners == {"midrun_v9/designs/bbbb.tar.zst"}


def test_iter_members_is_silent_about_absent_members(tmp_path):
    outer = _build_nested(tmp_path)
    wanted = {"midrun_v9/designs/aaaa.tar.zst": {"nope.cif"}}
    assert list(iter_members(outer, wanted)) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'binder_map.nested'`

- [ ] **Step 3: Write the implementation**

Create `src/binder_map/nested.py`:

```python
"""Read design structures out of the post-2026-09-11 nested archives.

On 2026-09-11 every midrun directory was collapsed into one solid
`midrun_vN.tar.zst` (683 GB across 35 archives). A design that the filter store
still records as

    <stack>/<midrun>/designs/<uuid>.tar.zst!<design_id>.cif

now lives two levels down:

    <stack>/<midrun>.tar.zst ! <midrun>/designs/<uuid>.tar.zst ! <design_id>.cif

Verified 2026-09-14: all 35 midruns referenced by the store's 37 design
templates are packed, none are live, and every design template has this one
shape. 40 of 40 sampled design archives are absent at their recorded paths.

A solid zstd stream has no random access, so the only sane access pattern is one
sequential pass per outer archive, pulling everything needed from it. That is
why callers group work by outer archive rather than by design.
"""
from __future__ import annotations

import io
import subprocess
import tarfile
from pathlib import Path
from typing import Iterator

_DESIGNS = "/designs/"


def nested_ref(design_ref: str) -> tuple[str, str, str]:
    """Rewrite a store design ref into (outer, inner, member) paths."""
    archive, member = design_ref.split("!", 1)
    head, uuid_tar = archive.split(_DESIGNS, 1)
    midrun = head.rsplit("/", 1)[-1]
    return f"{head}.tar.zst", f"{midrun}{_DESIGNS}{uuid_tar}", member


def _decompress(blob: bytes) -> bytes:
    return subprocess.run(
        ["zstd", "-dc"], input=blob, stdout=subprocess.PIPE, check=True
    ).stdout


def iter_members(outer: Path, wanted: dict[str, set[str]]) -> Iterator[tuple[str, str, bytes]]:
    """Stream `outer` exactly once, yielding (inner_name, member_name, blob).

    `wanted` maps an inner archive's name to the member names needed from it.

    Inner archives are decompressed in memory and dropped immediately. Nothing
    is written to scratch: unbounded expansion of archives filled /local/scratch
    to ~800 GB in July 2026 and produced blank-but-correct-length outputs, which
    is the failure mode a resume check cannot catch.

    Members that are absent are skipped silently; the caller reconciles counts
    against its manifest and records the misses with a reason.
    """
    proc = subprocess.Popen(
        ["zstd", "-dc", "-T0", str(outer)],
        stdout=subprocess.PIPE,
        bufsize=1 << 20,
    )
    try:
        with tarfile.open(fileobj=proc.stdout, mode="r|") as outer_tar:
            for info in outer_tar:
                if not info.isfile() or info.name not in wanted:
                    continue
                need = wanted[info.name]
                raw = _decompress(outer_tar.extractfile(info).read())
                with tarfile.open(fileobj=io.BytesIO(raw), mode="r|") as inner_tar:
                    for member in inner_tar:
                        if member.isfile() and member.name in need:
                            yield info.name, member.name, inner_tar.extractfile(member).read()
                del raw
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        proc.wait()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v`
Expected: PASS, 5 passed

- [ ] **Step 5: Verify against a real archive**

Run:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
from pathlib import Path
from binder_map.nested import iter_members
RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
outer = RB / "campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v5.tar.zst"
# name one inner archive from the .list sidecar, take whatever members it holds
listing = [l.strip() for l in open(str(outer) + ".list")
           if "/designs/" in l and l.strip().endswith(".tar.zst")]
inner = listing[0]
print("inner:", inner)
seen = [(i, m, len(b)) for i, m, b in iter_members(outer, {inner: {"__none__"}})]
print("members yielded for a name that cannot exist (expect 0):", len(seen))
assert seen == []
print("PASS: outer stream parsed and the inner archive was reached")
EOF
```

Expected: prints the inner path, then `PASS`. This confirms the outer stream
parses and the named inner archive is reachable, without paying for a full
extraction.

- [ ] **Step 6: Commit**

```bash
cd ~/Documents/binder_map
git add src/binder_map/nested.py tests/test_nested.py
git commit -m "Add nested-archive reader for the 2026-09-11 repack

Design refs are now two levels deep inside solid zstd streams with no random
access, so access inverts to one sequential pass per outer archive. Inner
archives are decompressed in memory and never expanded to scratch."
```

---

### Task 3: Stage 0 sampler

Draws the node set and writes the stage-1 manifest.

**Files:**
- Create: `scripts/stage0_sample.py`
- Test: `tests/test_stage0_sample.py`

**Interfaces:**
- Consumes: `backbone_key` (Task 1), `nested_ref` (Task 2)
- Produces: `build_nodes(store_dir, library_csv, passing_csv) -> pandas.DataFrame`
  with columns `bb, design_id, stack, campaign, gen, invfold, uuid8, outer_path,
  inner_path, member_path, n_designs, in_passing, in_12k`;
  `draw_sample(nodes, n, seed) -> numpy.ndarray` (boolean mask over `nodes` rows);
  `write_build(nodes, mask, out_dir)`.
  Writes three builds: `build_full/` (map 1), `build_lib12k/` (map 3) and
  `build_shared/` (their union — what stage 1 extracts), each with
  `manifests/design.tsv` and `meta/nodes.parquet`; plus `meta/sample.json`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_stage0_sample.py`:

```python
import numpy as np
import pandas as pd
import pytest

from stage0_sample import draw_sample


def _nodes(n: int) -> pd.DataFrame:
    return pd.DataFrame({"bb": [f"bb{i:06d}" for i in range(n)]})


def test_draw_is_exact_size_and_without_replacement():
    mask = draw_sample(_nodes(10_000), n=1_000, seed=42)
    assert mask.sum() == 1_000
    assert mask.dtype == bool
    assert len(mask) == 10_000


def test_draw_is_reproducible_for_a_seed():
    a = draw_sample(_nodes(10_000), n=1_000, seed=42)
    b = draw_sample(_nodes(10_000), n=1_000, seed=42)
    assert np.array_equal(a, b)


def test_a_different_seed_gives_a_different_draw():
    a = draw_sample(_nodes(10_000), n=1_000, seed=42)
    b = draw_sample(_nodes(10_000), n=1_000, seed=43)
    assert not np.array_equal(a, b)


def test_draw_larger_than_universe_raises():
    with pytest.raises(ValueError):
        draw_sample(_nodes(10), n=100, seed=42)


def test_draw_of_the_whole_universe_takes_everything():
    mask = draw_sample(_nodes(500), n=500, seed=42)
    assert mask.all()


def test_draw_is_uniform_across_an_arbitrary_label():
    nodes = _nodes(100_000)
    grp = np.arange(len(nodes)) % 4                 # four equal groups
    mask = draw_sample(nodes, n=20_000, seed=42)
    shares = pd.Series(grp[mask]).value_counts() / 20_000
    assert shares.between(0.24, 0.26).all()
```


- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_stage0_sample.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'stage0_sample'`

- [ ] **Step 3: Write the implementation**

Create `scripts/stage0_sample.py`:

```python
#!/usr/bin/env python3
"""Stage 0: draw the node set for the campaign-wide fold map.

Enumerates every coord-hash node in the filter store (2,952,126 across 6,992,515
designs) and writes THREE builds: a uniform random sample (map 1), the census of
backbones behind the shipped 12k library (map 3), and their union, which is what
stage 1 actually extracts so 683 GB is streamed once rather than twice.

Map 3 is NOT injected into map 1. Doing so would give its nodes inclusion
probability 1 against 0.169 for the rest and make map 1 cluster sizes unreadable;
the two are independent runs compared after the fact.

Uniform, not stack-stratified, and deliberately so: cluster size as generated is
what the map is read against, and stratifying would rescale exactly that. The
consequence is that the sample is ~85% BG+BIF and HD, because those two stacks
emit one backbone per design while the other three emit one per ten.

Usage:
    python scripts/stage0_sample.py --n 500000 --seed 42
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.nested import nested_ref
from binder_map.universe import backbone_key

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
FM = RB / "analyses" / "foldseek_map"
STORE = RB / "analyses" / "cross_campaign" / "filter_store_v2v9_hdtop"
LIBRARY = RB / "analyses" / "library_selection" / "csy4_binder_library_12k_v6_k30.csv"
PASSING = FM / "passing_v6_final.csv"

MANIFEST_COLUMNS = [
    "row_id", "design_id", "uuid8", "stack", "campaign", "gen", "invfold",
    "outer_path", "inner_path", "member_path",
]


def read_design_ids(csv: Path) -> set[str]:
    """Read design_id from a binder_filter export.

    The export's first two lines are comments and the header is line 3, so
    pandas needs skiprows=2. Reading it any other way yields a frame whose
    only column is a comment string, which then silently matches nothing.
    """
    return set(pd.read_csv(csv, skiprows=2, usecols=["design_id"])["design_id"])


def build_nodes(store: Path = STORE, library_csv: Path = LIBRARY,
                passing_csv: Path = PASSING) -> pd.DataFrame:
    """One row per unique backbone, with its representative design."""
    meta = json.loads((store / "meta.json").read_text())
    stacks = meta["vocabs"]["stack"]
    campaigns = meta["vocabs"]["campaign"]
    gens = meta["vocabs"]["gen"]
    invfolds = meta["vocabs"]["invfold"]

    load = lambda name: np.load(store / f"{name}.npy")
    design_id = load("design_id").astype(str)
    stack_code = load("stack")
    stack = np.asarray(stacks, dtype=object)[stack_code]

    frame = pd.DataFrame({
        "design_id": design_id,
        "stack": stack,
        "campaign": np.asarray(campaigns, dtype=object)[load("campaign")],
        "gen": np.asarray(gens, dtype=object)[load("gen")],
        "invfold": np.asarray(invfolds, dtype=object)[load("invfold")],
        "uuid8": load("uuid8").astype(str),
        "tpl": load("tpl_idx"),
    })
    frame["bb"] = [backbone_key(d, s) for d, s in zip(design_id, stack)]

    lib = read_design_ids(library_csv)
    passing = read_design_ids(passing_csv)
    frame["in_12k"] = frame["design_id"].isin(lib)
    frame["in_passing"] = frame["design_id"].isin(passing)

    # Sort so the representative is deterministic: the lexicographically
    # smallest design_id on each backbone.
    frame = frame.sort_values("design_id", kind="stable")
    grouped = frame.groupby("bb", sort=False)
    nodes = grouped.first().reset_index()
    nodes["n_designs"] = grouped.size().to_numpy()
    nodes["in_12k"] = grouped["in_12k"].any().to_numpy()
    nodes["in_passing"] = grouped["in_passing"].any().to_numpy()

    templates = meta["templates"]
    refs = [
        templates[t][2].replace("<ID>", d).replace("<U>", u)
        for t, d, u in zip(nodes["tpl"], nodes["design_id"], nodes["uuid8"])
    ]
    triples = [nested_ref(r) for r in refs]
    nodes["outer_path"] = [t[0] for t in triples]
    nodes["inner_path"] = [t[1] for t in triples]
    nodes["member_path"] = [t[2] for t in triples]
    return nodes.drop(columns=["tpl"])


def draw_sample(nodes: pd.DataFrame, n: int, seed: int) -> np.ndarray:
    """Uniform draw of `n` rows, returned as a boolean mask.

    A mask rather than a filtered frame because the same `nodes` table also
    feeds the library census and the union build, and re-deriving row identity
    from a filtered copy is how the two drift apart.
    """
    if n > len(nodes):
        raise ValueError(f"cannot draw {n} from a universe of {len(nodes)}")
    rng = np.random.default_rng(seed)
    mask = np.zeros(len(nodes), dtype=bool)
    mask[rng.choice(len(nodes), size=n, replace=False)] = True
    return mask


def write_build(nodes: pd.DataFrame, mask: np.ndarray, out_dir: Path) -> int:
    """Write one build's manifest and node table. Returns the row count."""
    selected = nodes.loc[mask].reset_index(drop=True)
    manifest = selected.assign(row_id=range(len(selected)))[MANIFEST_COLUMNS]
    manifest = manifest.sort_values(["outer_path", "inner_path"], kind="stable")
    (out_dir / "manifests").mkdir(parents=True, exist_ok=True)
    (out_dir / "meta").mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out_dir / "manifests" / "design.tsv", sep="\t", index=False)
    selected.to_parquet(out_dir / "meta" / "nodes.parquet", index=False)
    return len(selected)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--store", type=Path, default=STORE)
    parser.add_argument("--root", type=Path, default=FM)
    args = parser.parse_args()

    nodes = build_nodes(args.store)
    print(f"universe: {len(nodes):,d} coord-hash nodes from "
          f"{int(nodes['n_designs'].sum()):,d} designs")
    print(nodes.groupby("stack")["n_designs"].agg(["size", "sum"]))

    map1 = draw_sample(nodes, args.n, args.seed)
    map3 = nodes["in_12k"].to_numpy()
    shared = map1 | map3

    n1 = write_build(nodes, map1, args.root / "build_full")
    n3 = write_build(nodes, map3, args.root / "build_lib12k")
    ns = write_build(nodes, shared, args.root / "build_shared")
    print(f"map 1 (sample) {n1:,d} | map 3 (library) {n3:,d} | shared extract {ns:,d}")

    (args.root / "build_shared" / "meta").mkdir(parents=True, exist_ok=True)
    (args.root / "build_shared" / "meta" / "sample.json").write_text(json.dumps({
        "seed": args.seed,
        "universe": int(len(nodes)),
        "map1_n": int(n1),
        "map1_fraction": n1 / len(nodes),
        "map3_n": int(n3),
        "shared_n": int(ns),
        "overlap": int((map1 & map3).sum()),
        "outer_archives": int(nodes.loc[shared, "outer_path"].nunique()),
        "map1_per_stack": nodes.loc[map1].groupby("stack").size().to_dict(),
        "map3_per_stack": nodes.loc[map3].groupby("stack").size().to_dict(),
    }, indent=2))
    print(f"wrote builds under {args.root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_stage0_sample.py -v`
Expected: PASS, 6 passed

- [ ] **Step 5: Run it for real and check the universe reproduces**

Run: `cd $FM && PYTHONPATH=src:scripts $PY scripts/stage0_sample.py --n 500000 --seed 42`

Expected output must include `universe: 2,952,126 backbones from 6,992,515 designs`
and a per-stack table matching:

| stack | backbones | designs |
|---|---|---|
| BD1+LiMPNN | 144,468 | 1,444,614 |
| BG+BIF | 1,480,140 | 1,480,140 |
| BG+LiMPNN | 151,999 | 1,519,958 |
| HD+LiMPNN | 1,023,007 | 1,023,007 |
| RFD3+LiMPNN | 152,512 | 1,524,796 |

Then verify the selection:

```bash
cd $FM && $PY -c "
import json; m=json.load(open('build_shared/meta/sample.json'))
print(json.dumps(m, indent=1))
assert m['universe']==2952126, m['universe']
assert m['map1_n']==500000
assert m['map3_n']==9760, m['map3_n']
assert m['shared_n']==m['map1_n']+m['map3_n']-m['overlap']
assert m['outer_archives']==35, m['outer_archives']
print('PASS')"
```

Expected: `PASS`, with `overlap` ≈ 1,650 and `shared_n` ≈ 508,100.
`map3_n` **must** be exactly 9,760 (stage 2's coord-hash dedup later collapses
these to 9,759 graph nodes — do not compare the two counts) — a smaller number means the library CSV join
failed, almost always because the export header is on line 3 and `skiprows=2`
was dropped.

- [ ] **Step 6: Commit**

```bash
cd ~/Documents/binder_map
git add scripts/stage0_sample.py tests/test_stage0_sample.py
git commit -m "Add stage 0 sampler emitting three builds

Uniform 500k draw over all 2,952,126 coord-hash nodes at seed 42 (map 1), the
9,759-node library census (map 3), and their union so stage 1 streams 683 GB
once. The two maps stay independent runs rather than one map with forced nodes,
which would make map 1 cluster sizes unreadable."
```

---

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

### Task 5: Split the shared shards into two result trees

Stage 1 extracts map 1 and map 3 in one 683 GB pass. This splits its output so
stages 2–5 see two ordinary result trees and need no changes at all.

**Files:**
- Create: `scripts/stage1b_split.py`
- Test: `tests/test_split.py`

**Interfaces:**
- Consumes: `results_shared/design/{shards,hashes,fp,missing}` (Task 4), the
  per-map manifests from Task 3
- Produces: `manifest_ids(build_dir: Path) -> set[str]`;
  `split_tables(frame: pandas.DataFrame, keep: set[str]) -> pandas.DataFrame`.
  Writes `results_full/design/…` and `results_lib12k/design/…` with the same
  four subdirectories.

- [ ] **Step 1: Write the failing test**

Create `tests/test_split.py`:

```python
import io
import tarfile

import pandas as pd

from stage1b_split import split_tables, split_tar


def test_split_tables_keeps_only_requested_designs():
    frame = pd.DataFrame({"design_id": ["a", "b", "c"], "x": [1, 2, 3]})
    out = split_tables(frame, {"a", "c"})
    assert list(out["design_id"]) == ["a", "c"]


def test_split_tables_on_empty_keep_returns_empty_with_columns():
    frame = pd.DataFrame({"design_id": ["a"], "x": [1]})
    out = split_tables(frame, set())
    assert len(out) == 0
    assert list(out.columns) == ["design_id", "x"]


def test_split_tables_ignores_ids_that_are_not_present():
    frame = pd.DataFrame({"design_id": ["a"], "x": [1]})
    assert list(split_tables(frame, {"a", "ghost"})["design_id"]) == ["a"]


def _tar(path, members):
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            tar.addfile(info, io.BytesIO(data))


def test_split_tar_copies_only_matching_members(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"AAA\n", "b.pdb": b"BBB\n", "c.pdb": b"CCC\n"})
    dst = tmp_path / "out.tar.gz"
    n = split_tar(src, dst, {"a", "c"})
    assert n == 2
    with tarfile.open(dst, "r:gz") as tar:
        assert sorted(m.name for m in tar) == ["a.pdb", "c.pdb"]


def test_split_tar_returns_zero_when_nothing_matches(tmp_path):
    src = tmp_path / "shard_00.tar.gz"
    _tar(src, {"a.pdb": b"AAA\n"})
    assert split_tar(src, tmp_path / "out.tar.gz", {"zz"}) == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_split.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'stage1b_split'`

- [ ] **Step 3: Write the implementation**

Create `scripts/stage1b_split.py`:

```python
#!/usr/bin/env python3
"""Stage 1b: split the shared extraction into one result tree per map.

Stage 1 streams 683 GB once and extracts the union of map 1 (the 500k sample)
and map 3 (the 9,759-node library census). Streaming that volume twice for the
sake of 9,759 structures would be indefensible, so the split happens here
instead -- one cheap pass over ~2.5 GB of shards.

The output is two ordinary result trees, so stages 2-5 run against them
unmodified.

Usage:
    python scripts/stage1b_split.py --shared results_shared \
        --map full:build_full:results_full --map lib12k:build_lib12k:results_lib12k
"""
from __future__ import annotations

import argparse
import io
import tarfile
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def manifest_ids(build_dir: Path) -> set[str]:
    frame = pd.read_csv(build_dir / "manifests" / "design.tsv",
                        sep="\t", usecols=["design_id"])
    return set(frame["design_id"])


def split_tables(frame: pd.DataFrame, keep: set[str]) -> pd.DataFrame:
    return frame[frame["design_id"].isin(keep)].reset_index(drop=True)


def split_tar(src: Path, dst: Path, keep: set[str]) -> int:
    """Copy members whose stem is in `keep`. Returns the number copied."""
    written = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(src, "r:gz") as source, tarfile.open(dst, "w:gz") as target:
        for info in source:
            if not info.isfile():
                continue
            stem = info.name[:-4] if info.name.endswith(".pdb") else info.name
            if stem not in keep:
                continue
            blob = source.extractfile(info).read()
            out = tarfile.TarInfo(name=info.name)
            out.size = len(blob)
            out.mtime = 0
            target.addfile(out, io.BytesIO(blob))
            written += 1
    return written


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared", type=Path, default=FM / "results_shared")
    parser.add_argument("--map", action="append", required=True,
                        help="name:build_dir:results_dir")
    parser.add_argument("--root", type=Path, default=FM)
    args = parser.parse_args()

    src = args.shared / "design"
    for spec in args.map:
        name, build, results = spec.split(":")
        keep = manifest_ids(args.root / build)
        out = args.root / results / "design"
        for sub in ("shards", "hashes", "fp", "missing"):
            (out / sub).mkdir(parents=True, exist_ok=True)

        copied = 0
        for shard in sorted((src / "shards").glob("shard_*.tar.gz")):
            copied += split_tar(shard, out / "shards" / shard.name, keep)

        rows = 0
        for table in sorted((src / "hashes").glob("shard_*.tsv")):
            frame = split_tables(pd.read_csv(table, sep="\t"), keep)
            frame.to_csv(out / "hashes" / table.name, sep="\t", index=False)
            rows += len(frame)
        for table in sorted((src / "fp").glob("shard_*.parquet")):
            split_tables(pd.read_parquet(table), keep).to_parquet(
                out / "fp" / table.name, index=False)
        for table in sorted((src / "missing").glob("shard_*.tsv")):
            split_tables(pd.read_csv(table, sep="\t"), keep).to_csv(
                out / "missing" / table.name, sep="\t", index=False)

        print(f"{name}: {rows:,d} of {len(keep):,d} manifest ids parsed, "
              f"{copied:,d} backbone structures")
        # Map 3 is a census, not a sample: a shortfall there is a real loss and
        # must be visible, not averaged away.
        if name == "lib12k" and rows < len(keep):
            print(f"  WARNING: {len(keep) - rows:,d} library designs did not parse")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_split.py -v`
Expected: PASS, 5 passed

- [ ] **Step 5: Commit**

```bash
cd ~/Documents/binder_map
git add scripts/stage1b_split.py tests/test_split.py
git commit -m "Split the shared extraction into one result tree per map

Stage 1 streams 683 GB once for the union of the 500k sample and the 9,759-node
library census; this splits its shards and tables so stages 2-5 run against two
ordinary result trees unmodified."
```

---

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

### Task 7: Stage 9 map comparison

Reads all three finished maps and compares them. No gate is recomputed anywhere —
membership comes from joins against exports that already exist.

**Files:**
- Create: `scripts/stage9_compare.py`
- Test: `tests/test_compare_maps.py`

**Interfaces:**
- Consumes: `clusters/fold/clusters.parquet`, `dedup/nodes.tsv`,
  `dedup/members.tsv` for each of the three maps; `meta/nodes.parquet` per build
- Produces: `expand_members(members) -> DataFrame[struct_name, design_id]`;
  `design_to_cluster(results: Path, resolution: str) -> pandas.Series`;
  `summarise(labels: pandas.DataFrame, name: str) -> dict`;
  `crosswalk_scores(old: pandas.Series, new: pandas.Series) -> dict` with keys
  `n_shared, ari, ami`.
  Writes `compare/map_summary.parquet`, `compare/cluster_table_<map>.parquet`,
  `compare/crosswalk.parquet`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_compare_maps.py`:

```python
import pandas as pd
import pytest

from stage9_compare import crosswalk_scores, expand_members, fold_key, summarise


def _labels():
    return pd.DataFrame({
        "cluster":   [0, 0, 0, 1, 1, 2],
        "stack":     ["RFD3+LiMPNN"] * 3 + ["BG+BIF"] * 2 + ["HD+LiMPNN"],
        "in_12k":    [True, False, True, False, False, False],
        "in_passing":[True, True, True, True, False, False],
        "n_designs": [10, 10, 10, 1, 1, 1],
    })


def test_summarise_counts_nodes_and_clusters():
    got = summarise(_labels(), "m1")
    assert got["map"] == "m1"
    assert got["nodes"] == 6
    assert got["clusters"] == 3


def test_summarise_counts_singletons():
    # cluster 2 has exactly one member
    assert summarise(_labels(), "m1")["singletons"] == 1


def test_summarise_reports_membership_totals():
    got = summarise(_labels(), "m1")
    assert got["n_12k"] == 2
    assert got["n_passing"] == 4


def test_summarise_reports_largest_cluster():
    assert summarise(_labels(), "m1")["largest_cluster"] == 3


def test_fold_key_collapses_bg_bif_seq_variants():
    ids = [f"csy4free_BG_BIF_len80_0000df81_design0_seq{i}" for i in range(4)]
    assert len({fold_key(i) for i in ids}) == 1


def test_fold_key_collapses_hd_sample_variants():
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample8_pack1"
    assert fold_key(a) == fold_key(b)


def test_fold_key_keeps_different_hd_recycles_apart():
    a = "csy4free_HD_LiMPNN_len80_001a4983_recycle10_sample0_pack1"
    b = "csy4free_HD_LiMPNN_len80_001a4983_recycle9_sample0_pack1"
    assert fold_key(a) != fold_key(b)


def test_expand_members_splits_the_design_id_list():
    members = pd.DataFrame({"struct_name": ["s0", "s1"],
                            "design_ids": ["d0,d1, d2", "d3"]})
    out = expand_members(members)
    assert list(out["design_id"]) == ["d0", "d1", "d2", "d3"]
    assert list(out["struct_name"]) == ["s0", "s0", "s0", "s1"]


def test_crosswalk_scores_are_one_for_an_identical_partition():
    a = pd.Series([0, 0, 1, 1], index=["a", "b", "c", "d"])
    b = pd.Series([7, 7, 9, 9], index=["a", "b", "c", "d"])
    got = crosswalk_scores(a, b)
    assert got["n_shared"] == 4
    assert got["ari"] == pytest.approx(1.0)
    assert got["ami"] == pytest.approx(1.0)


def test_crosswalk_scores_use_only_shared_designs():
    a = pd.Series([0, 0, 1], index=["a", "b", "c"])
    b = pd.Series([5, 5, 6], index=["b", "c", "zz"])
    assert crosswalk_scores(a, b)["n_shared"] == 2


def test_crosswalk_scores_on_no_overlap_are_nan():
    import math
    a = pd.Series([0], index=["a"])
    b = pd.Series([1], index=["z"])
    got = crosswalk_scores(a, b)
    assert got["n_shared"] == 0
    assert math.isnan(got["ari"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_compare_maps.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'stage9_compare'`

- [ ] **Step 3: Write the implementation**

Create `scripts/stage9_compare.py`:

```python
#!/usr/bin/env python3
"""Stage 9: compare the three fold maps.

Three clusterings built under identical conditions, differing only in sample
size: the 500k campaign sample, the 31,763-node filtered pool, and the
9,759-node ordered library.

Nothing here recomputes a gate. `in_passing` and `in_12k` were set in stage 0 by
a join against exports that already exist and are already audited.

Cluster IDS ARE NOT COMPARABLE ACROSS MAPS -- three independent Leiden runs on
three node sets. The crosswalk reports ARI/AMI on designs present in more than
one map as a RESULT, not as a correspondence to be relied on.

Usage:
    python scripts/stage9_compare.py --out compare
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")

#: The Leiden resolution the 12k selection consumed (v6hier_q.py:61).
RESOLUTION = "cluster_r80"

#: name -> (results dir, build dir or None for the pre-existing reference map)
MAPS = {
    "sample500k": (FM / "results_full" / "design", FM / "build_full"),
    "passing60k": (FM / "results_sw_k30_f50_mut" / "design", None),
    "library12k": (FM / "results_lib12k" / "design", FM / "build_lib12k"),
}


def expand_members(members: pd.DataFrame) -> pd.DataFrame:
    """One row per design: (struct_name, design_id).

    A graph node is a coordinate hash standing for one or more design_ids, so
    per-design labels must be expanded before they can be joined.
    """
    out = members.assign(design_id=members["design_ids"].str.split(",")).explode("design_id")
    out["design_id"] = out["design_id"].str.strip()
    return out[["struct_name", "design_id"]].reset_index(drop=True)


def design_to_cluster(results: Path, resolution: str = RESOLUTION) -> pd.Series:
    """design_id -> fold cluster id, for one completed map."""
    clusters = pd.read_parquet(results / "clusters" / "fold" / "clusters.parquet")
    nodes = pd.read_csv(results / "dedup" / "nodes.tsv", sep="\t")
    members = expand_members(pd.read_csv(results / "dedup" / "members.tsv", sep="\t"))
    joined = (members.merge(nodes[["struct_name", "node_id"]], on="struct_name")
                     .merge(clusters[["node_id", resolution]], on="node_id"))
    return joined.set_index("design_id")[resolution]


#: Collapse the dedup key's known over-splitting: BG+BIF seq variants are the same
#: fold rigid-body displaced (0.94 A superposed RMSD), and HD sample variants
#: differ by 0.022 A of packer jitter. Both measured 2026-09-14.
_SEQ = re.compile(r"_seq\d+$")
_SAMPLE_PACK = re.compile(r"_sample\d+_pack\d+$")


def fold_key(design_id: str) -> str:
    """The distinct FOLD a node belongs to, as opposed to its coord-hash node."""
    if "_HD_LiMPNN_" in design_id:
        return _SAMPLE_PACK.sub("", design_id)
    return _SEQ.sub("", design_id)


def summarise(labels: pd.DataFrame, name: str) -> dict:
    sizes = labels.groupby("cluster").size()
    return {
        "map": name,
        "nodes": int(len(labels)),
        "designs": int(labels["n_designs"].sum()),
        "clusters": int(labels["cluster"].nunique()),
        "singletons": int((sizes == 1).sum()),
        "largest_cluster": int(sizes.max()),
        "median_cluster": float(sizes.median()),
        "p90_cluster": float(sizes.quantile(0.90)),
        "n_passing": int(labels["in_passing"].sum()),
        "n_12k": int(labels["in_12k"].sum()),
    }


def crosswalk_scores(old: pd.Series, new: pd.Series) -> dict:
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

    shared = old.index.intersection(new.index)
    if len(shared) == 0:
        return {"n_shared": 0, "ari": float("nan"), "ami": float("nan")}
    return {
        "n_shared": int(len(shared)),
        "ari": float(adjusted_rand_score(old.loc[shared], new.loc[shared])),
        "ami": float(adjusted_mutual_info_score(old.loc[shared], new.loc[shared])),
    }


def load_labels(results: Path, build: Path | None) -> pd.DataFrame:
    clusters = pd.read_parquet(results / "clusters" / "fold" / "clusters.parquet")
    nodes = pd.read_csv(results / "dedup" / "nodes.tsv", sep="\t")
    members = expand_members(pd.read_csv(results / "dedup" / "members.tsv", sep="\t"))

    if build is not None:
        meta = pd.read_parquet(build / "meta" / "nodes.parquet")[
            ["design_id", "stack", "in_passing", "in_12k", "n_designs"]]
    else:
        # The reference map predates stage 0, so rebuild its labels the same way:
        # every design in it is by definition in `passing`, and 12k membership
        # comes from the library export.
        lib = set(pd.read_csv(
            FM.parent / "library_selection" / "csy4_binder_library_12k_v6_k30.csv",
            usecols=["design_id"])["design_id"])
        meta = pd.DataFrame({"design_id": members["design_id"].unique()})
        meta["stack"] = "unknown"
        meta["in_passing"] = True
        meta["in_12k"] = meta["design_id"].isin(lib)
        meta["n_designs"] = 1

    joined = members.merge(meta, on="design_id", how="left")
    per_node = joined.groupby("struct_name").agg(
        stack=("stack", lambda x: x.dropna().iloc[0] if x.notna().any() else "unknown"),
        in_passing=("in_passing", "any"),
        in_12k=("in_12k", "any"),
        n_designs=("n_designs", "sum"),
    ).reset_index()

    labels = (nodes.merge(per_node, on="struct_name", how="left")
                   .merge(clusters[["node_id", RESOLUTION]], on="node_id", how="left")
                   .rename(columns={RESOLUTION: "cluster"}))
    labels[["in_passing", "in_12k"]] = labels[["in_passing", "in_12k"]].fillna(False)
    labels["n_designs"] = labels["n_designs"].fillna(0)
    labels["fold_key"] = [fold_key(d) for d in labels["struct_name"]]
    return labels.dropna(subset=["cluster"]).astype({"cluster": int})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=FM / "compare")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    summaries, per_design = [], {}
    for name, (results, build) in MAPS.items():
        if not (results / "clusters" / "fold" / "clusters.parquet").exists():
            print(f"SKIP {name}: not built yet")
            continue
        labels = load_labels(results, build)
        summaries.append(summarise(labels, name))

        # Node count AND distinct-fold count, side by side. The dedup key
        # over-splits BG+BIF ~9.8x and HD ~10x (rigid-body offset and 0.022 A
        # packer jitter), so `size` alone is an inflated denominator for any
        # coverage claim. `folds` is what a fold-level reading needs.
        table = labels.groupby("cluster").agg(
            size=("struct_name", "size"),
            folds=("fold_key", "nunique"),
            n_designs=("n_designs", "sum"),
            n_passing=("in_passing", "sum"),
            n_12k=("in_12k", "sum"),
            top_stack=("stack", lambda x: x.value_counts().idxmax()),
        ).reset_index().sort_values("size", ascending=False)
        table.to_parquet(args.out / f"cluster_table_{name}.parquet", index=False)
        per_design[name] = design_to_cluster(results)

    summary = pd.DataFrame(summaries)
    summary.to_parquet(args.out / "map_summary.parquet", index=False)
    print(summary.to_string(index=False))

    rows = []
    names = list(per_design)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            scores = crosswalk_scores(per_design[a], per_design[b])
            scores.update({"map_a": a, "map_b": b})
            rows.append(scores)
            print(f"{a} vs {b}: n={scores['n_shared']:,d} "
                  f"ARI={scores['ari']:.3f} AMI={scores['ami']:.3f}")
    pd.DataFrame(rows).to_parquet(args.out / "crosswalk.parquet", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_compare_maps.py -v`
Expected: PASS, 11 passed

- [ ] **Step 5: Commit**

```bash
cd ~/Documents/binder_map
git add scripts/stage9_compare.py tests/test_compare_maps.py
git commit -m "Add three-map comparison stage

Summarises each map's cluster-size distribution and membership, and reports
pairwise ARI/AMI on designs present in more than one map. Cluster ids are not
comparable across independent Leiden runs, so the crosswalk is a result rather
than a correspondence."
```

---

### Task 8: Scaling probe, then both runs

The ~2,100 CPU-h stage-3 figure is an N² extrapolation from the v6 run
(31,763 nodes → 8.5 CPU-h), not a measurement. Measure it before committing a
128-task array.

**Files:**
- Create: `probe_scaling.sh`
- Create: `submit_maps.sh`

**Interfaces:**
- Consumes: everything above
- Produces: nothing importable; these are drivers.

- [ ] **Step 1: Write the probe**

Create `probe_scaling.sh`:

```bash
#!/bin/bash
# Measure how Foldseek search actually scales on this node set before
# committing the full run. Builds two sub-DBs from the already-built stage-2
# database and times one query chunk against each.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
RESULTS=${RESULTS:-results_full}
cd "$FM" || exit 1

for N in 50000 150000; do
  sbatch --job-name="probe_$N" --partition=normal.4h --time=04:00:00 \
    --cpus-per-task=16 --mem-per-cpu=4G \
    --output="$FM/logs/probe_${N}_%j.log" \
    --wrap "bash -lc '
      set -e
      module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960
      OUT=$FM/$RESULTS/design
      WORK=\${TMPDIR:-/tmp}/probe_$N
      mkdir -p \$WORK; trap \"rm -rf \$WORK\" EXIT
      head -n $N \$OUT/db/uniq.lookup | cut -f1 > \$WORK/keys.txt
      head -n 2000 \$WORK/keys.txt > \$WORK/qkeys.txt
      for suffix in \"\" _ss _ca; do
        foldseek createsubdb \$WORK/keys.txt  \$OUT/db/uniq\$suffix \$WORK/sub\$suffix
        foldseek createsubdb \$WORK/qkeys.txt \$OUT/db/uniq\$suffix \$WORK/q\$suffix
      done
      echo \"=== N=$N : 2000 queries against a $N-target DB ===\"
      /usr/bin/time -v foldseek search \$WORK/q \$WORK/sub \$WORK/aln \$WORK/tmp \
        -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16
    '"
done
echo "Two probes submitted. When both finish:"
echo "  grep -H 'Elapsed (wall clock)' $FM/logs/probe_*.log"
echo "p = log(t150/t50)/log(3). Full map-1 search ~ t50 * (500000/50000)^p * (500000/2000)."
```

- [ ] **Step 2: Write the submit driver**

Create `submit_maps.sh`:

```bash
#!/bin/bash
# Three-map build: one shared 683 GB extraction, split per map, then stages 2-5
# per map, then the comparison. Map 2 (results_sw_k30_f50_mut) already exists
# and is never rebuilt.
# Parameters are frozen to match results_sw_k30_f50_mut exactly.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
cd "$FM" || exit 1
LAYER=design

# --- stage 1: one task per outer archive, over the UNION manifest -------------
J1=$(sbatch --parsable --job-name=cwfm_s1 --array=0-34%35 \
     --export=ALL,N_SHARDS=35,BUILD=build_shared,RESULTS=results_shared \
     slurm/stage1_extract_nested.sbatch) || exit 1
echo "stage1 extract (shared)  $J1"

# --- stage 1b: split into the two per-map trees ------------------------------
J1B=$(sbatch --parsable --job-name=cwfm_s1b --dependency=afterok:$J1 \
      --partition=normal.4h --time=02:00:00 --cpus-per-task=4 --mem-per-cpu=4G \
      --output="$FM/logs/cwfm_s1b_%j.log" \
      --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage1b_split.py \
        --shared $FM/results_shared \
        --map full:build_full:results_full \
        --map lib12k:build_lib12k:results_lib12k'") || exit 1
echo "stage1b split            $J1B"

# --- stages 2-5, once per map ------------------------------------------------
submit_map () {
  local TAG=$1 RESULTS=$2 CHUNKS=$3 THROTTLE=$4 DEP=$5
  local J2 J3 J4 J5
  J2=$(sbatch --parsable --job-name=${TAG}_s2 --dependency=afterok:$DEP \
       --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
         PYTHONPATH=src $PY scripts/stage3_split_queries.py --layer $LAYER \
         --n-chunks $CHUNKS --results $FM/$RESULTS'" \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s2_%j.log") || exit 1
  J3=$(sbatch --parsable --job-name=${TAG}_s3 --dependency=afterok:$J2 \
       --array=0-$((CHUNKS-1))%${THROTTLE} \
       --export=ALL,ALN_TYPE=1,RESULTS=$RESULTS \
       slurm/stage3_search.sbatch $LAYER) || exit 1
  J4=$(sbatch --parsable --job-name=${TAG}_s4 --dependency=afterok:$J3 \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=16 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s4_%j.log" \
       --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage4_graph.py \
         --layer $LAYER --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual \
         --results $FM/$RESULTS'") || exit 1
  J5=$(sbatch --parsable --job-name=${TAG}_s5 --dependency=afterok:$J4 \
       slurm/stage_single.sbatch scripts/stage5_cluster.py \
         --layer $LAYER --axis fold --results $FM/$RESULTS) || exit 1
  echo "$TAG: s2=$J2 s3=$J3 s4=$J4 s5=$J5" >&2
  echo "$J5"
}

J5_LIB=$(submit_map lib12k results_lib12k 4   4  $J1B)
J5_FULL=$(submit_map full  results_full  128 64 $J1B)

# --- stage 9: compare all three ---------------------------------------------
J9=$(sbatch --parsable --job-name=cwfm_s9 \
     --dependency=afterok:$J5_FULL:$J5_LIB \
     --partition=normal.4h --time=01:00:00 --cpus-per-task=8 --mem-per-cpu=4G \
     --output="$FM/logs/cwfm_s9_%j.log" \
     --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage9_compare.py'") || exit 1
echo "stage9 compare           $J9"
echo "$J9" > $FM/logs/last_job_cwfm.txt
```

- [ ] **Step 3: Run stage 0, then launch**

```bash
cd $FM && PYTHONPATH=src:scripts $PY scripts/stage0_sample.py --n 500000 --seed 42
cd $FM && bash submit_maps.sh
```

Expected: seven job ids printed, `map3_n` in `build_shared/meta/sample.json`
exactly 9,759.

- [ ] **Step 4: Gate on hash identity before paying for the search**

When `cwfm_s1b` completes — **before** `full_s3` starts — verify the repack
rewrite did not change what gets extracted:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import glob
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"

def hashes(results):
    files = sorted(glob.glob(f"{FM}/{results}/design/hashes/shard_*.tsv"))
    return pd.concat([pd.read_csv(f, sep="\t") for f in files]).set_index("design_id")["coord_hash"]

new, old = hashes("results_shared"), hashes("results_v6")
shared = new.index.intersection(old.index)
print(f"design_ids in both runs: {len(shared):,d}")
assert len(shared) > 5000, "too few shared ids to be a real check"
mismatch = int((new.loc[shared] != old.loc[shared]).sum())
print(f"coord_hash mismatches: {mismatch}")
assert mismatch == 0, "extraction changed -- STOP, do not run the search"
print("PASS: extraction is byte-identical to the v6 run")
EOF
```

Expected: `coord_hash mismatches: 0` and `PASS`. If it fails, `scancel` the
`full_s3`/`s4`/`s5` ids printed by `submit_maps.sh` before ~2,000 CPU-h is spent.

- [ ] **Step 4b: Check the shard member names match the manifest ids**

Stage 0 writes the manifests; stage 1 writes the PDB member names; stage 1b
matches one against the other by stripping `.pdb`. If those two ever disagree on
id format by a single character, **every shard produces a well-formed, EMPTY
output tar and nothing raises** — the map would come out valid and empty. Check
it before the search, not after:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import glob, tarfile
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"

man = set(pd.read_csv(f"{FM}/build_shared/manifests/design.tsv",
                      sep="\t", usecols=["design_id"])["design_id"])
shard = sorted(glob.glob(f"{FM}/results_shared/design/shards/shard_*.tar.gz"))[0]
with tarfile.open(shard, "r:gz") as tar:
    stems = [m.name[:-4] for m in tar if m.isfile() and m.name.endswith(".pdb")]
print(f"shard {shard.split('/')[-1]}: {len(stems):,d} members")
assert stems, "shard has no .pdb members at all"
hit = sum(s in man for s in stems)
print(f"members whose stem is a manifest design_id: {hit:,d} / {len(stems):,d}")
print("examples:", stems[:2])
assert hit == len(stems), "ID FORMAT MISMATCH -- stage 1b would emit empty tars"

# and the split itself must not be silently empty
for tree in ("results_full", "results_lib12k"):
    tars = sorted(glob.glob(f"{FM}/{tree}/design/shards/shard_*.tar.gz"))
    assert tars, f"{tree}: no shards written"
    with tarfile.open(tars[0], "r:gz") as tar:
        n = sum(1 for m in tar if m.isfile())
    print(f"{tree}: first shard holds {n:,d} members")
    assert n > 0, f"{tree}: split produced an EMPTY tar"
print("PASS")
EOF
```

Expected: every stem is a manifest `design_id`, and both split trees have
non-empty shards. A mismatch here is the one failure that produces a complete,
plausible, empty result — fix it before spending the search.

- [ ] **Step 5: Check the library map is a complete census**

```bash
cd $FM && $PY -c "
import glob, pandas as pd
FM='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map'
n=pd.read_csv(f'{FM}/results_lib12k/design/dedup/nodes.tsv',sep='\t')
miss=pd.concat([pd.read_csv(f,sep='\t') for f in
                sorted(glob.glob(f'{FM}/results_lib12k/design/missing/shard_*.tsv'))])
print('library nodes:', len(n), 'missing:', len(miss))
if len(miss): print(miss.reason.value_counts())
assert len(n) + 0 >= 9700, 'library census is short -- investigate before comparing'
print('PASS')"
```

Expected: `library nodes: 9759` (or a shortfall fully explained by `missing/`)
and `PASS`. Map 3 is a census, not a sample — a silent shortfall here biases
every comparison that follows.

- [ ] **Step 6: Spot-check TM identity against the v6 hits**

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
COLS = ["query","target","qtmscore","ttmscore","alntmscore",
        "lddt","fident","alnlen","evalue","bits"]
read = lambda f: pd.read_csv(f, sep="\t", names=COLS, header=None,
                             usecols=["query","target","qtmscore","ttmscore"])
new = pd.concat(read(f"{FM}/results_lib12k/design/hits/hits_0{i}.tsv") for i in range(4))
old = pd.concat(read(f"{FM}/results_v6/design/hits/hits_{i:02d}.tsv") for i in range(8))
m = new.merge(old, on=["query","target"], suffixes=("_new","_old"))
print(f"shared query/target pairs: {len(m):,d}")
assert len(m) > 100, "too few shared pairs -- widen the chunks compared"
d = (m["qtmscore_new"] - m["qtmscore_old"]).abs()
print(f"max |dTM| = {d.max():.6f}")
assert d.max() < 1e-3, "Foldseek scores drifted -- STOP"
print("PASS: TM scores reproduce the v6 run")
EOF
```

Expected: `max |dTM|` at or near 0 and `PASS`. The library map is used here
because it finishes first and its node set overlaps map 2 almost entirely.

- [ ] **Step 7: Read the comparison**

When `cwfm_s9` completes:

```bash
cd $FM && $PY -c "
import pandas as pd
s=pd.read_parquet('compare/map_summary.parquet'); print(s.to_string(index=False))
c=pd.read_parquet('compare/crosswalk.parquet');  print(c.to_string(index=False))
assert len(s)==3, 'expected all three maps'
print('PASS')"
```

Expected: three rows — `sample500k`, `passing60k`, `library12k` — with node
counts 500,000 / 31,763 / 9,759, and three pairwise ARI/AMI rows.

- [ ] **Step 8: Commit**

```bash
cd ~/Documents/binder_map
for f in probe_scaling.sh submit_maps.sh; do
  done
git add probe_scaling.sh submit_maps.sh
git commit -m "Add scaling probe and the three-map submit driver

One shared extraction feeds both new maps; stages 2-5 then run once per map with
frozen parameters. The chain gates on a coord-hash identity check against
results_v6 before the 500k search is paid for."
```

---

## Notes for whoever executes this

- **Stages 2, 3 and 5 are not modified.** If a change to them looks necessary,
  stop and re-read the spec — the clustering parameters are the one thing this
  work must not perturb, because comparability with map 2 is the whole point.
- **Map 2 is never rebuilt.** `results_sw_k30_f50_mut` is the reference.
- **`map3_n` must be exactly 9,760**, which stage 2 then dedups to 9,759 graph
  nodes. These are different stages; comparing them is an error the spec made once.
  The library-selection export has its header on **line 1** (no comment preamble),
  unlike the binder_filter exports which have two — count leading `#` lines rather
  than hardcoding `skiprows`.
- **Watch `missing/` after stage 1.** A large `member_missing` count means the
  manifest's nested paths disagree with the archives; check `nested_ref` against
  a real `.list` sidecar before blaming the data.
- **Do not "fix" the dedup key.** It over-splits BG+BIF ~9.8× (rigid-body offset,
  0.94 Å superposed RMSD) and HD ~10× (0.022 Å packer jitter). Map 2 has the same
  defect; matching it is what makes the three maps comparable.
- **Cluster ids are not comparable across maps.** Three independent Leiden runs on
  three node sets. Compare size distributions and composition; the crosswalk's
  ARI/AMI is a result to report, not a correspondence to rely on.
