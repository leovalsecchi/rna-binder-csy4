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
Expected: PASS, 8 passed

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

