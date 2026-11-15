#!/usr/bin/env python3
"""Stage U0: assemble one union result tree out of two finished ones.

The three maps of this campaign were built separately -- the 500k campaign
sample and the 9,759-node library census out of `results_shared`, the 31,763
passing pool out of `results_v6` -- so they each have their own embedding and
cannot be looked at in one picture. This stage puts all three populations into
a single tree of 526,467 backbones so stages 2-6 can co-embed them, and records
which of the three each node belongs to; a node carrying more than one flag is
conserved across populations, which is the whole point of the exercise.

No archive is re-extracted: the shards are already on disk, so they are
hard-linked, not copied.

Two things have to be reconciled on the way in, and both of them are failures
that have already happened:

1. **More hash rows than structures.** Stage 1 appends a hashes row for every
   design but writes only ONE structure per coord_hash. `results_shared` holds
   508,207 rows against 508,205 structures; `results_v6` holds 55,706 rows
   against 31,767. Stage 2 elects the lexicographically first design_id per
   coord_hash, so an unfiltered table elects ids the Foldseek DB has no key
   for and stage 2 aborts with
       FATAL: N canonical structures are not in the Foldseek DB.
   That killed a run on 2026-09-24. `filter_hashes` drops those rows here,
   against the structures actually present, and they are recorded in missing/.

2. **The same structure in both trees.** 6,376 design_ids have a `.pdb` in
   `results_shared` *and* in `results_v6` -- the library census designs, which
   are a subset of the passing pool. Copying both would put two entries under
   one name in `db/full.lookup`, `stage2_builddb.sh` would resolve each
   canonical name to two keys, and its `resolved N / KEEP` check would fail on
   the overshoot instead of the shortfall. `install_shard` keeps the first
   copy and rewrites the later shard without it. Measured: all 6,376 of those
   ids are elected canonical, so this is not a corner case.

Usage:
    python scripts/stageU0_union.py \
        --source results_shared --source results_v6 --out results_union \
        --layer-nodes sample500k:results_full/design/dedup/nodes.tsv \
        --layer-nodes passing60k:results_sw_k30_f50_mut/design/dedup/nodes.tsv \
        --layer-nodes library12k:results_lib12k/design/dedup/nodes.tsv
"""
from __future__ import annotations

import argparse
import io
import os
import shutil
import tarfile
from pathlib import Path

import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")

#: Reason recorded in missing/ for a design whose hashes row survives stage 1
#: but whose structure does not: another design with the same coord_hash was
#: written first, so this one has no `.pdb` of its own. Same reason string as
#: stage 1b, and for the same event.
DUPLICATE_REASON = "duplicate_structure"

#: missing/ rows written here carry a coord_hash on top of stage 1's schema.
MISSING_COLUMNS = ["design_id", "reason", "coord_hash"]

#: Layer name on the command line -> flag column in meta/layers.parquet. The
#: library column is `in_12k`, not `in_library12k`: it is what the explorer
#: page reads.
LAYER_COLUMNS = {
    "sample500k": "in_sample500k",
    "passing60k": "in_passing60k",
    "library12k": "in_12k",
}

#: Measured on 2026-09-28 against the trees on disk. These are asserted, not
#: assumed: if the assembly does not reproduce them something upstream moved,
#: and the tree must not be handed to stage 2 until that is understood.
EXPECTED_UNION = 526_467
EXPECTED_LAYER = {"sample500k": 499_999, "passing60k": 31_763, "library12k": 9_759}
EXPECTED_CONSERVED = {
    ("sample500k", "passing60k"): 5_295,
    ("sample500k", "library12k"): 1_553,
    ("passing60k", "library12k"): 9_759,
    ("sample500k", "passing60k", "library12k"): 1_553,
}


def source_tag(source: str | Path) -> str:
    """`results_shared` -> `shared`; the tag that keeps two trees apart."""
    name = Path(source).name
    prefix = "results_"
    return name[len(prefix):] if name.startswith(prefix) else name


def shard_name(tag: str, original: str) -> str:
    """`shard_00.tar.gz` + `shared` -> `shard_shared_00.tar.gz`.

    The tag goes *after* `shard_`, not in front of it. Every consumer in this
    repo globs `shard_*` -- stage2_builddb.py on hashes, stage4_graph.py on fp,
    stage7_medoids.py and stage7b_all_traces.py and gate_checks.py and
    validate.py on shards -- so a `shared_shard_00.tar.gz` would be invisible
    to all of them while still looking right in a directory listing.
    """
    head, sep, rest = original.partition("_")
    if not sep:
        raise ValueError(f"not a shard filename: {original}")
    return f"{head}_{tag}_{rest}"


def shard_stems(path: Path) -> set[str]:
    """Stems of the `.pdb` members of one shard."""
    stems: set[str] = set()
    with tarfile.open(path, "r:gz") as tar:
        for info in tar:
            if info.isfile() and info.name.endswith(".pdb"):
                stems.add(info.name[:-4])
    return stems


def present_stems(shard_dir: Path) -> set[str]:
    """Every structure that exists across a tree's shards.

    This is the set the hashes table has to be reconciled against, and it is
    read back off the assembled shards rather than accumulated during the
    assembly so that what is on disk, not what was intended, has the last word.
    """
    stems: set[str] = set()
    for shard in sorted(Path(shard_dir).glob("shard_*.tar.gz")):
        stems |= shard_stems(shard)
    return stems


def link_or_copy(src: Path, dst: Path) -> bool:
    """Hard-link `src` at `dst`; copy only if the link is impossible.

    The shards are 2.5 GB and the project is under a hard inode cap, so the
    union tree borrows them rather than duplicating them. Returns True if a
    hard link was made.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
        return True
    except OSError:
        shutil.copy2(src, dst)
        return False


def filter_shard(src: Path, dst: Path, drop: set[str]) -> int:
    """Rewrite `src` at `dst` without the members whose stem is in `drop`."""
    written = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(src, "r:gz") as source, tarfile.open(dst, "w:gz") as target:
        for info in source:
            if not info.isfile():
                continue
            stem = info.name[:-4] if info.name.endswith(".pdb") else info.name
            if stem in drop:
                continue
            blob = source.extractfile(info).read()
            out = tarfile.TarInfo(name=info.name)
            out.size = len(blob)
            out.mtime = 0
            target.addfile(out, io.BytesIO(blob))
            written += 1
    return written


def install_shard(src: Path, dst: Path, claimed: set[str]) -> tuple[set[str], bool]:
    """Put one source shard into the union tree exactly once.

    Hard-links the shard when none of its structures is already in the tree,
    which is the common case and costs nothing. When a structure is already
    there -- the 6,376 designs `results_shared` and `results_v6` both carry --
    the shard is rewritten without it, because the same `.pdb` name twice in
    the Foldseek DB makes stage 2's canonical-name lookup resolve to two keys.

    Returns (the stems this shard contributed, whether it had to be rewritten).
    """
    stems = shard_stems(src)
    overlap = stems & claimed
    if not overlap:
        link_or_copy(src, dst)
        return stems, False
    filter_shard(src, dst, overlap)
    return stems - overlap, True


def dedup_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse the rows two sources both carry for the same design.

    6,377 design_ids appear in both source hash tables. Left alone each would
    be counted twice in its node's `n_designs` and listed twice in
    `dedup/members.tsv`. A design that carries two *different* coord_hashes is
    a different problem -- two backbones under one name -- and is refused here
    rather than silently resolved by row order.
    """
    pairs = frame[["design_id", "coord_hash"]].drop_duplicates()
    clashing = pairs.loc[pairs["design_id"].duplicated(), "design_id"]
    if len(clashing):
        raise ValueError(
            f"{clashing.nunique()} design_id(s) carry more than one coord_hash, "
            f"e.g. {sorted(set(clashing))[:3]}")
    return frame.drop_duplicates(subset="design_id").reset_index(drop=True)


def filter_hashes(frame: pd.DataFrame,
                  present: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Partition a hashes table into (rows with a structure, rows without).

    Every coord_hash keeps at least one row: stage 1 dropped the structure of a
    design only because another design with the *same* coord_hash already had
    one, and that one is in `present`. So this narrows each node's member list
    without ever removing a node.
    """
    have = frame["design_id"].isin(present)
    kept = frame[have].reset_index(drop=True)
    dropped = frame.loc[~have, ["design_id", "coord_hash"]].copy()
    dropped["reason"] = DUPLICATE_REASON
    return kept, dropped[MISSING_COLUMNS].reset_index(drop=True)


def merge_missing(missing: pd.DataFrame, dropped: pd.DataFrame) -> pd.DataFrame:
    """Append this stage's dropped rows to a source shard's missing table.

    The source rows carry design_id + reason; these carry a coord_hash too, so
    a dropped design still points at the node it collapsed into.
    """
    if dropped.empty:
        return missing.reset_index(drop=True)
    return pd.concat([missing, dropped], ignore_index=True)


def layer_flags(nodes_by_layer: dict[str, set[str]],
                hash_of: dict[str, str]) -> pd.DataFrame:
    """One row per union node, flagged with the populations it belongs to.

    Each map elected its own canonical design per backbone, so the same node
    reaches this function under up to three different names; mapping every
    name through the union hash table is what collapses them onto one row. A
    name with no coord_hash would silently lose its flag and undercount the
    conserved set, so it is refused.
    """
    unknown = set(nodes_by_layer) ^ set(LAYER_COLUMNS)
    if unknown:
        raise ValueError(f"layers must be exactly {sorted(LAYER_COLUMNS)}, "
                         f"mismatched: {sorted(unknown)}")

    hashes_by_layer: dict[str, set[str]] = {}
    for layer, names in nodes_by_layer.items():
        stranded = {n for n in names if n not in hash_of}
        if stranded:
            raise ValueError(
                f"{len(stranded)} {layer} node(s) have no coord_hash in the "
                f"union hash table, e.g. {sorted(stranded)[:3]}")
        hashes_by_layer[layer] = {hash_of[n] for n in names}

    every = sorted(set().union(*hashes_by_layer.values()))
    frame = pd.DataFrame({"coord_hash": every})
    for layer, column in LAYER_COLUMNS.items():
        frame[column] = frame["coord_hash"].isin(hashes_by_layer[layer])
    return frame


def check_counts(union: int, layers: dict[str, int],
                 conserved: dict[tuple[str, ...], int]) -> list[str]:
    """Compare the assembled tree against the measured sizes.

    Returns one message per mismatch; an empty list means the tree is the one
    the plan was written against.
    """
    problems = []
    if union != EXPECTED_UNION:
        problems.append(f"union is {union:,d} nodes, expected {EXPECTED_UNION:,d}")
    for layer, want in EXPECTED_LAYER.items():
        got = layers.get(layer)
        if got != want:
            problems.append(f"{layer} is {got:,d} nodes, expected {want:,d}")
    for combo, want in EXPECTED_CONSERVED.items():
        got = conserved.get(combo)
        if got != want:
            problems.append(f"conserved {' n '.join(combo)} is {got:,d}, "
                            f"expected {want:,d}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description="Assemble the union tree.")
    parser.add_argument("--source", action="append", required=True,
                        help="a finished result tree, e.g. results_shared. "
                             "Repeatable; the order decides which tree's copy "
                             "of a shared structure is the one kept.")
    parser.add_argument("--out", type=Path, default=Path("results_union"))
    parser.add_argument("--layer-nodes", action="append", required=True,
                        metavar="NAME:PATH",
                        help="layer name and that map's dedup/nodes.tsv")
    parser.add_argument("--layer", default="design",
                        help="structure layer (design/af3/esm); default design")
    parser.add_argument("--root", type=Path, default=FM)
    args = parser.parse_args()

    out = args.root / args.out / args.layer
    for sub in ("shards", "hashes", "fp", "missing", "meta"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    claimed: set[str] = set()
    hash_frames: list[pd.DataFrame] = []
    hash_names: list[str] = []
    source_missing: dict[str, pd.DataFrame] = {}
    linked = rewritten = 0

    for source in args.source:
        src = args.root / source / args.layer
        tag = source_tag(source)
        for shard in sorted((src / "shards").glob("shard_*.tar.gz")):
            stems, was_rewritten = install_shard(
                shard, out / "shards" / shard_name(tag, shard.name), claimed)
            claimed |= stems
            rewritten += was_rewritten
            linked += not was_rewritten
        for table in sorted((src / "hashes").glob("shard_*.tsv")):
            frame = pd.read_csv(table, sep="\t")
            frame["_shard"] = shard_name(tag, table.name)
            hash_names.append(shard_name(tag, table.name))
            hash_frames.append(frame)
        for table in sorted((src / "missing").glob("shard_*.tsv")):
            source_missing[shard_name(tag, table.name)] = pd.read_csv(table, sep="\t")
        for table in sorted((src / "fp").glob("shard_*.parquet")):
            link_or_copy(table, out / "fp" / shard_name(tag, table.name))
        print(f"{source}: shards installed as shard_{tag}_*")

    print(f"shards: {linked:,d} hard-linked, {rewritten:,d} rewritten "
          f"(a rewrite means the shard held a structure another source "
          f"already contributed)")

    # Read the structures back off the assembled tree, so what is on disk --
    # not what the install thought it wrote -- decides which hashes rows live.
    present = present_stems(out / "shards")
    if present != claimed:
        print(f"FATAL: the assembled shards hold {len(present):,d} structures "
              f"but the install accounted for {len(claimed):,d}")
        return 1

    hashes = pd.concat(hash_frames, ignore_index=True)
    rows_in = len(hashes)
    hashes = dedup_rows(hashes)
    route = hashes.set_index("design_id")["_shard"]
    kept, dropped = filter_hashes(hashes, present)
    dropped["_shard"] = dropped["design_id"].map(route)

    print(f"hashes: {rows_in:,d} rows in, {rows_in - len(hashes):,d} duplicated "
          f"across sources, {len(kept):,d} kept, {len(dropped):,d} -> missing/ "
          f"({DUPLICATE_REASON})")

    columns = [c for c in kept.columns if c != "_shard"]
    kept_by = dict(tuple(kept.groupby("_shard", sort=False)))
    for name in hash_names:
        frame = kept_by.get(name, kept.iloc[0:0])
        frame[columns].to_csv(out / "hashes" / name, sep="\t", index=False)

    empty = pd.DataFrame(columns=MISSING_COLUMNS)
    dropped_by = {name: group.drop(columns="_shard")
                  for name, group in dropped.groupby("_shard", sort=False)}
    for name, frame in sorted(source_missing.items()):
        merge_missing(frame, dropped_by.pop(name, empty)).to_csv(
            out / "missing" / name, sep="\t", index=False)
    # A hashes shard with no companion missing table is belt and braces, but
    # its dropped rows still have to land somewhere.
    for name, frame in sorted(dropped_by.items()):
        frame.to_csv(out / "missing" / name, sep="\t", index=False)

    hash_of = dict(zip(kept["design_id"], kept["coord_hash"]))
    nodes_by_layer: dict[str, set[str]] = {}
    for spec in args.layer_nodes:
        name, sep, path = spec.partition(":")
        if not sep:
            print(f"FATAL: --layer-nodes wants NAME:PATH, got {spec!r}")
            return 1
        frame = pd.read_csv(args.root / path, sep="\t", usecols=["struct_name"])
        nodes_by_layer[name] = set(frame["struct_name"])

    flags = layer_flags(nodes_by_layer, hash_of)
    flags.to_parquet(out / "meta" / "layers.parquet", index=False)

    union = int(kept["coord_hash"].nunique())
    layers = {layer: int(flags[column].sum())
              for layer, column in LAYER_COLUMNS.items()}
    conserved = {combo: int(flags[[LAYER_COLUMNS[k] for k in combo]]
                            .all(axis=1).sum())
                 for combo in EXPECTED_CONSERVED}

    print(f"union: {union:,d} nodes, {len(present):,d} structures")
    for layer, count in layers.items():
        print(f"  {layer}: {count:,d}")
    for combo, count in conserved.items():
        print(f"  conserved {' n '.join(combo)}: {count:,d}")

    if len(flags) != union:
        print(f"FATAL: meta/layers.parquet has {len(flags):,d} rows for "
              f"{union:,d} union nodes; a node with no layer cannot be drawn")
        return 1

    problems = check_counts(union, layers, conserved)
    for problem in problems:
        print(f"FATAL: {problem}")
    if problems:
        print("The union tree does not match the measured sizes. Do not run "
              "stage 2 on it until that is understood.")
        return 1

    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
