#!/usr/bin/env python3
"""Task 9: the hard correctness gate that stands between stage 2 and stage 3.

Stage 3 for the `full` map is a 512-task array costing ~2,100 CPU-h. The plan
put three correctness checks in the window between stage 2 finishing and that
array starting, to be run by hand and aborted with `scancel` if they failed.
That is a race over a VPN that drops, and one of the failures it guards is
silent: if manifest `design_id`s and shard `.pdb` member stems disagree by a
single character, every shard writes a well-formed EMPTY tar, nothing raises,
and a 500k map comes out valid-looking and built on nothing.

This module turns that manual window into a Slurm `afterok` dependency. The
whole design rests on failing closed: a non-zero exit, a crash, or a timeout
all leave stage 3 in `DependencyNeverSatisfied`. So there is deliberately NO
broad exception handling here -- an unreadable table or a corrupt tar must
propagate, not be swallowed into a pass.

The checks are pure functions over paths so they are testable without a
cluster; `main()` is a thin verdict-printing wrapper.

Usage:
    python scripts/gate_checks.py --shared results_shared --ref results_v6 \
        --manifest build_shared/manifests/design.tsv \
        --tree results_full --tree results_lib12k
"""
from __future__ import annotations

import argparse
import tarfile
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")

#: Below this many shared design_ids the hash comparison is not evidence of
#: anything -- an empty or near-empty intersection would otherwise "pass".
#: Measured 2026-09-24 on Euler: build_shared/manifests/design.tsv holds 508,207
#: design_ids, results_v6 holds 55,706, and they intersect in 6,377. Since
#: results_shared's ids are a subset of the manifest, 6,377 is a hard upper bound
#: on n_shared -- a 5000 floor left 27% of margin and would hard-block a good run
#: after ~22% extraction losses. The floor only has to make a vacuous comparison
#: impossible; a single hash mismatch fails the check at any n_shared.
MIN_SHARED_HASHES = 3000

#: check_split_nonempty's sentinel for "this tree has no shards at all",
#: distinct from "it has a shard and the shard is empty" (0).
NO_SHARDS = -1


def _read_hashes(results: Path) -> dict[str, str]:
    """design_id -> coord_hash over every `design/hashes/shard_*.tsv` in a tree.

    A missing directory yields an empty mapping rather than an error: the
    caller's `n_shared` then falls under MIN_SHARED_HASHES and the gate fails,
    which is the correct outcome for a tree that was never built.
    """
    out: dict[str, str] = {}
    for table in sorted((Path(results) / "design" / "hashes").glob("shard_*.tsv")):
        frame = pd.read_csv(table, sep="\t", usecols=["design_id", "coord_hash"])
        out.update(zip(frame["design_id"], frame["coord_hash"]))
    return out


def check_hash_identity(shared_results: Path, ref_results: Path) -> tuple[int, int]:
    """Compare coord_hashes for the design_ids two trees have in common.

    Returns (n_shared, n_mismatch). The reference tree is `results_v6`, built
    before the 2026-09-11 repack; agreement proves the repack did not change
    what gets extracted. Caller fails on n_shared < MIN_SHARED_HASHES or
    n_mismatch > 0.
    """
    have = _read_hashes(Path(shared_results))
    ref = _read_hashes(Path(ref_results))
    shared_ids = have.keys() & ref.keys()
    n_mismatch = sum(1 for design_id in shared_ids if have[design_id] != ref[design_id])
    return len(shared_ids), n_mismatch


def check_id_format(manifest_tsv: Path, shard_targz: Path) -> tuple[int, int]:
    """Check that every `.pdb` member stem in a shard is a manifest design_id.

    Returns (n_members, n_hits). This is the silent failure: stage 1b keeps a
    member only when its stem is in the manifest set, so a one-character
    disagreement makes every split shard empty without raising. Caller fails
    on n_members == 0 or n_hits != n_members.
    """
    ids = set(pd.read_csv(manifest_tsv, sep="\t", usecols=["design_id"])["design_id"])
    stems = []
    with tarfile.open(shard_targz, "r:gz") as tar:
        for info in tar:
            if info.isfile() and info.name.endswith(".pdb"):
                stems.append(info.name[:-4])
    n_hits = sum(1 for stem in stems if stem in ids)
    return len(stems), n_hits


def check_split_nonempty(tree_results: Path) -> int:
    """Members in a tree's first shard, or NO_SHARDS if it has none.

    Stage 2 feeds `foldseek createdb` straight from these tars, so an empty
    first shard is the fingerprint of a split that matched nothing. Caller
    fails on NO_SHARDS or 0.
    """
    shards = sorted((Path(tree_results) / "design" / "shards").glob("shard_*.tar.gz"))
    if not shards:
        return NO_SHARDS
    n_members = 0
    with tarfile.open(shards[0], "r:gz") as tar:
        for info in tar:
            if info.isfile() and info.name.endswith(".pdb"):
                n_members += 1
    return n_members


def first_shard(results: Path) -> Path | None:
    shards = sorted((Path(results) / "design" / "shards").glob("shard_*.tar.gz"))
    return shards[0] if shards else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shared", type=Path, default=FM / "results_shared",
                        help="tree whose hashes and first shard are checked")
    parser.add_argument("--ref", type=Path, default=FM / "results_v6",
                        help="read-only reference tree (pre-repack hashes)")
    parser.add_argument("--manifest", type=Path,
                        default=FM / "build_shared" / "manifests" / "design.tsv")
    parser.add_argument("--tree", type=Path, action="append", default=None,
                        help="result tree to check for non-empty shards (repeatable)")
    args = parser.parse_args(argv)
    trees = args.tree or []

    ok = True

    # --- check 1: the repack did not change what gets extracted --------------
    # Comparing a tree with itself is guaranteed to agree and proves nothing, so
    # it must not be a way to pass. Latent today (the sbatch passes distinct
    # paths) and cheap to keep that way.
    if Path(args.shared).resolve() == Path(args.ref).resolve():
        print(f"[hash identity] --shared and --ref both resolve to "
              f"{Path(args.shared).resolve()}")
        print("  FAIL: a tree compared with itself is not evidence of anything")
        print("GATE FAIL")
        return 1
    n_shared, n_mismatch = check_hash_identity(args.shared, args.ref)
    print(f"[hash identity] {args.shared} vs {args.ref}: "
          f"shared_ids={n_shared:,d} mismatches={n_mismatch:,d} "
          f"(need shared_ids >= {MIN_SHARED_HASHES:,d} and mismatches == 0)")
    if n_shared < MIN_SHARED_HASHES:
        print(f"  FAIL: only {n_shared:,d} design_ids in common -- too few to "
              f"be evidence of anything")
        ok = False
    if n_mismatch > 0:
        print(f"  FAIL: {n_mismatch:,d} coord_hash mismatches against the "
              f"reference tree")
        ok = False
    if n_shared >= MIN_SHARED_HASHES and n_mismatch == 0:
        print("  PASS")

    # --- check 2: manifest ids and shard member stems are the same strings ---
    shard = first_shard(args.shared)
    if shard is None:
        print(f"[id format] {args.shared}: no shard_*.tar.gz found")
        print("  FAIL: nothing to check the manifest against")
        ok = False
    else:
        n_members, n_hits = check_id_format(args.manifest, shard)
        print(f"[id format] {args.manifest} vs {shard}: "
              f"pdb_members={n_members:,d} matching_manifest_ids={n_hits:,d}")
        if n_members == 0:
            print("  FAIL: shard holds no .pdb members")
            ok = False
        elif n_hits != n_members:
            print(f"  FAIL: {n_members - n_hits:,d} member stems are not manifest "
                  f"design_ids -- the split would write empty shards")
            ok = False
        else:
            print("  PASS")

    # --- check 3: the split trees actually hold structures -------------------
    for tree in trees:
        n_first = check_split_nonempty(tree)
        if n_first == NO_SHARDS:
            print(f"[split non-empty] {tree}: no shard_*.tar.gz found")
            print("  FAIL: tree has no shards")
            ok = False
        else:
            print(f"[split non-empty] {tree}: first shard holds "
                  f"{n_first:,d} .pdb members")
            if n_first == 0:
                print("  FAIL: first shard is empty")
                ok = False
            else:
                print("  PASS")
    if not trees:
        print("[split non-empty] no --tree given")
        print("  FAIL: refusing to pass a gate that checked no trees")
        ok = False

    print("GATE PASS" if ok else "GATE FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
