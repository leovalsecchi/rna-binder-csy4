#!/usr/bin/env python3
"""Stage 1: extract backbones and compute epitope fingerprints.

One pass over the archives does both jobs, because opening 8,931 .tar.zst
files twice would double the most expensive part of the pipeline for no gain.

Everything is streamed through $TMPDIR into a single .tar.gz shard. No loose
per-structure files are ever written under `results/`: the project has a
hard inode cap, and 215k loose files would breach it.

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


def read_archive_members(archive: Path, members: list[str], workdir: Path) -> dict[str, bytes]:
    """Extract the needed members from one .tar.zst into workdir, keyed by name.

    Python's tarfile has no zstd codec here, so tar does the decompression
    once and we sieve out only the requested members - but by NAME, into
    real files, not via `-O` to a single stdout stream.

    `-O` concatenates the matched members in ARCHIVE order, not the order
    they were requested in (verified: `tar --zstd -xOf a.tar.zst b.cif
    a.cif` still prints a.cif's content first). Zipping `members` against
    the resulting blobs positionally therefore silently swaps content
    between design_ids whenever the manifest's row order - inherited from
    passing.csv, not from the tar - differs from the archive's own member
    order. That is not a corner case here: most archives hold more than one
    member (84% of them, per build/manifests/design.tsv), so this would
    have quietly mislabelled the majority of extracted structures. Named
    extraction into files sidesteps the ordering assumption entirely.

    Also handles partial misses correctly: if some requested members are
    absent, `tar` exits non-zero but still writes out the ones it found. We
    don't gate on the exit code - we just check which files landed.
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
    extract_dir = tmp / "members"

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
            blobs = read_archive_members(archive, members, extract_dir)

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

    shutil.rmtree(tmp, ignore_errors=True)

    print(f"done in {time.time() - t0:.0f}s: {len(hashes):,d} parsed, "
          f"{len(seen_local):,d} unique backbones, {len(missing):,d} missing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
