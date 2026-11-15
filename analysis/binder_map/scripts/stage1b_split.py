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

#: Reason recorded in missing/ for a design whose hashes row survives stage 1
#: but whose structure does not: another design with the same coord_hash was
#: written first, so this one has no `.pdb` of its own. The row carries the
#: coord_hash, so the representative structure is still recoverable.
DUPLICATE_REASON = "duplicate_structure"

#: missing/ rows written here carry a coord_hash on top of stage 1's schema.
MISSING_COLUMNS = ["design_id", "reason", "coord_hash"]


def manifest_ids(build_dir: Path) -> set[str]:
    frame = pd.read_csv(build_dir / "manifests" / "design.tsv",
                        sep="\t", usecols=["design_id"])
    return set(frame["design_id"])


def split_tables(frame: pd.DataFrame, keep: set[str]) -> pd.DataFrame:
    return frame[frame["design_id"].isin(keep)].reset_index(drop=True)


def split_tar(src: Path, dst: Path, keep: set[str]) -> tuple[int, set[str]]:
    """Copy members whose stem is in `keep`.

    Returns (number copied, the stems actually written). The caller needs the
    stems and not just the count: stage 1 appends a hashes row for every
    design but writes one structure per coord_hash per task, so the hashes
    table is a superset of the shards and this loop is the only place that
    knows which side of that split each design landed on.
    """
    written = 0
    stems: set[str] = set()
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
            stems.add(stem)
    return written, stems


def split_hashes(frame: pd.DataFrame,
                 present: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Partition a hashes table into (rows with a structure, rows without).

    Stage 1 appends a hashes row for every design_id but skips *writing* the
    structure when that coord_hash was already seen in the task, so the hashes
    table is a superset of the structures on disk. Stage 2 then elects the
    lexicographically first design_id per coord_hash as canonical; when that
    election lands on a design whose structure was deduped away, the Foldseek
    DB has no key for it and stage 2 aborts with a misleading "createdb
    silently skipped them". Drop those rows here -- this pass has just been
    told exactly which structures exist -- and hand them back so the caller
    can record them in missing/ instead of losing them.
    """
    have = frame["design_id"].isin(present)
    kept = frame[have].reset_index(drop=True)
    dropped = frame.loc[~have, ["design_id", "coord_hash"]].copy()
    dropped["reason"] = DUPLICATE_REASON
    return kept, dropped[MISSING_COLUMNS].reset_index(drop=True)


def merge_missing(missing: pd.DataFrame, dropped: pd.DataFrame) -> pd.DataFrame:
    """Append the deduped-away hashes rows to a shard's missing table.

    Stage 1's missing rows carry design_id + reason; these carry a coord_hash
    too, so a dropped design still points at the node it collapsed into. The
    concat leaves coord_hash empty for stage 1's own rows, which never had one.
    """
    if dropped.empty:
        return missing.reset_index(drop=True)
    return pd.concat([missing, dropped], ignore_index=True)


def reconcile(keep: set[str], hashes_rows: int, missing_rows: int) -> tuple[int, int]:
    """Return (explained, unexplained) counts.

    explained = rows present in hashes + rows in missing (accounting for all designs)
    unexplained = designs in the manifest that appear in neither hashes nor missing
    """
    accounted = hashes_rows + missing_rows
    explained = accounted
    unexplained = len(keep) - accounted
    return explained, unexplained


def should_warn_census(is_census: bool, unexplained: int) -> bool:
    """Return True if a census warning should be issued.

    Census maps must have all designs accounted for. An unexplained shortfall
    in a census is a data loss event and must be loudly reported.
    """
    return is_census and unexplained > 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared", type=Path, default=FM / "results_shared")
    parser.add_argument("--map", action="append", required=True,
                        help="name:build_dir:results_dir[:census]")
    parser.add_argument("--root", type=Path, default=FM)
    args = parser.parse_args()

    src = args.shared / "design"
    for spec in args.map:
        parts = spec.split(":")
        if len(parts) < 3 or len(parts) > 4:
            raise ValueError(f"--map must be 'name:build:results' or "
                           f"'name:build:results:census', got: {spec}")
        name, build, results = parts[:3]
        is_census = len(parts) == 4 and parts[3] == "census"
        if len(parts) == 4 and parts[3] != "census":
            raise ValueError(f"--map optional 4th component must be 'census', "
                           f"got: {parts[3]}")

        keep = manifest_ids(args.root / build)
        out = args.root / results / "design"
        for sub in ("shards", "hashes", "fp", "missing"):
            (out / sub).mkdir(parents=True, exist_ok=True)

        copied = 0
        present: set[str] = set()
        for shard in sorted((src / "shards").glob("shard_*.tar.gz")):
            written, stems = split_tar(shard, out / "shards" / shard.name, keep)
            copied += written
            present |= stems

        # Keep a hashes row only when its structure is on disk in THIS tree, so
        # stage 2's election cannot pick a design Foldseek has no key for.
        hashes_rows = 0
        deduped = 0
        dropped: dict[str, pd.DataFrame] = {}
        for table in sorted((src / "hashes").glob("shard_*.tsv")):
            frame = split_tables(pd.read_csv(table, sep="\t"), keep)
            frame, gone = split_hashes(frame, present)
            frame.to_csv(out / "hashes" / table.name, sep="\t", index=False)
            hashes_rows += len(frame)
            if len(gone):
                dropped[table.name] = gone
                deduped += len(gone)

        empty = pd.DataFrame(columns=MISSING_COLUMNS)
        missing_rows = 0
        for table in sorted((src / "missing").glob("shard_*.tsv")):
            frame = split_tables(pd.read_csv(table, sep="\t"), keep)
            frame = merge_missing(frame, dropped.pop(table.name, empty))
            frame.to_csv(out / "missing" / table.name, sep="\t", index=False)
            missing_rows += len(frame)

        # Stage 1 writes hashes/ and missing/ together, so a hashes shard with
        # no companion missing table is belt and braces -- but its dropped rows
        # still have to be written, or reconcile() reports them as a shortfall.
        for table_name, gone in sorted(dropped.items()):
            gone.to_csv(out / "missing" / table_name, sep="\t", index=False)
            missing_rows += len(gone)

        for table in sorted((src / "fp").glob("shard_*.parquet")):
            split_tables(pd.read_parquet(table), keep).to_parquet(
                out / "fp" / table.name, index=False)

        print(f"{name}: {hashes_rows:,d} of {len(keep):,d} manifest ids parsed, "
              f"{copied:,d} backbone structures "
              f"({deduped:,d} deduped -> missing/)")

        explained, unexplained = reconcile(keep, hashes_rows, missing_rows)
        if should_warn_census(is_census, unexplained):
            print(f"  WARNING: {unexplained:,d} library designs unaccounted for "
                  f"({missing_rows:,d} in missing/, {hashes_rows:,d} parsed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
