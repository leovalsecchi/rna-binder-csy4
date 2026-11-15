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
