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

from binder_map import RNA_LEN
from binder_map.backbone import coord_hash, write_backbone_pdb
from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.fingerprint import FINGERPRINT_COLUMNS, fingerprint
from binder_map.geometry import approach
from binder_map.nested import iter_members
from binder_map.rnaframe import FrameError

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
FM = RB / "analyses" / "foldseek_map"

# Explicit schema for the fingerprint output, so an all-missing shard (e.g.
# every outer archive for it is absent) still writes a parquet with the right
# columns instead of a schema-less, zero-column one.
FP_COLUMNS = (["design_id"] + FINGERPRINT_COLUMNS
              + [f"per_nt_{i + 1}" for i in range(RNA_LEN)]
              + ["contacts_offframe", "contacts_total",
                 "cx", "cy", "cz", "ax", "ay", "az"])


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

    # Keyed by (outer_path, inner_path, member_path), not just (inner_path,
    # member_path): the same inner-archive/member pair can legitimately recur
    # across two different outer archives (different midruns), and keying on
    # the pair alone silently collapses two distinct manifest rows onto one.
    #
    # Built and validated before any directory (scratch included) is created,
    # so a collision fails loudly with nothing left behind to clean up.
    by_member = {(r.outer_path, r.inner_path, r.member_path): r.design_id
                 for r in mine.itertuples()}
    assert len(by_member) == len(mine), (
        f"by_member key collision: {len(mine)} manifest rows collapsed to "
        f"{len(by_member)} unique (outer_path, inner_path, member_path) keys "
        "- two rows share the same outer/inner/member triple")

    out_dir = args.results / "design"
    for sub in ("shards", "hashes", "fp", "missing"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"bm_s1n_{args.shard}"
    tmp.mkdir(parents=True, exist_ok=True)
    scratch = tmp / "one.cif"

    shard_name = f"shard_{args.shard:02d}"
    tar_path = out_dir / "shards" / f"{shard_name}.tar.gz"
    # Written under this name while open; only renamed to `tar_path` once the
    # accounting assert and the FATAL checks below have passed, so a shard
    # that fails validation never leaves a complete-looking tarball behind.
    partial_tar_path = tar_path.with_name(tar_path.name + ".partial")

    hashes, fps, missing = [], [], []
    seen_local: set[str] = set()
    # Every manifest row this shard has resolved a final outcome for -
    # delivered and parsed, delivered but rejected, or accounted for via an
    # archive_missing/archive_read_error branch - so the reconciliation loop
    # below doesn't re-add it a second time as member_missing.
    accounted: set[tuple[str, str, str]] = set()
    t0 = time.time()

    with tarfile.open(partial_tar_path, "w:gz") as tar:
        for outer_path, group in mine.groupby("outer_path", sort=True):
            outer = RB / outer_path
            if not outer.exists():
                for row in group.itertuples():
                    key = (row.outer_path, row.inner_path, row.member_path)
                    accounted.add(key)
                    missing.append({"design_id": row.design_id, "reason": "archive_missing"})
                continue
            wanted = plan_archive_work(group)
            print(f"  streaming {outer_path} "
                  f"({len(wanted)} inner archives, {len(group):,d} members)", flush=True)

            # The try/except below is scoped to ONLY the generator advance
            # (the `next(it)` call), not the per-member work that follows it.
            # A bug in write_backbone_pdb, tar.addfile or record construction
            # must not be relabelled archive_read_error and silently abandon
            # the rest of this archive - it stays caught by the existing
            # per-member `try` further down, exactly as before.
            it = iter_members(outer, wanted)
            while True:
                try:
                    inner_name, member_name, blob = next(it)
                except StopIteration:
                    break
                except Exception as exc:  # noqa: BLE001 - a corrupt archive must not kill the shard
                    # Raised by iter_members itself (a CalledProcessError from
                    # a corrupt inner archive, a tarfile.ReadError on a
                    # truncated outer stream, ...) while advancing the
                    # generator. Whatever of this archive's rows are already
                    # accounted for stay accounted for; the rest are recorded
                    # as missing so the shard can still finish and a rerun
                    # isn't doomed to fail identically.
                    print(f"  WARNING: {outer_path} raised {type(exc).__name__} while "
                          "streaming; marking its undelivered rows as archive_read_error",
                          flush=True)
                    for row in group.itertuples():
                        key = (row.outer_path, row.inner_path, row.member_path)
                        if key not in accounted:
                            accounted.add(key)
                            missing.append({"design_id": row.design_id,
                                            "reason": f"archive_read_error:{type(exc).__name__}"})
                    break

                key = (outer_path, inner_name, member_name)
                design_id = by_member.get(key)
                if design_id is None:
                    continue
                if key in accounted:
                    # iter_members yielded this exact (outer, inner, member)
                    # a second time - e.g. an inner tar carrying duplicate
                    # member names. Already accounted for; a second delivery
                    # must be a no-op, not a second row in hashes/fp and not
                    # a second hit against the final accounting assert.
                    continue
                accounted.add(key)
                if not blob.strip():
                    # Present in the archive but blank - distinct from
                    # never having been delivered at all (member_missing
                    # below). Blank-but-right-length outputs are this
                    # project's own known failure mode; keep them tellable
                    # apart from a genuinely absent member.
                    missing.append({"design_id": design_id, "reason": "member_empty"})
                    continue
                try:
                    scratch.write_bytes(blob)
                    atoms = load_structure(scratch)
                    protein, rna = split_chains(atoms)
                    chain = str(protein.chain_id[0])
                    n_res = len(set(zip(protein.chain_id, protein.res_id)))
                    if n_res != 80:
                        # Chain selection succeeded; the length gate is
                        # what failed.
                        missing.append({"design_id": design_id,
                                        "reason": f"wrong_length:n_res={n_res}"})
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
        if key not in accounted:
            missing.append({"design_id": design_id, "reason": "member_missing"})

    shutil.rmtree(tmp, ignore_errors=True)

    # Validate BEFORE writing anything. Fail loudly rather than leave a
    # plausible-looking empty shard on disk: a blank output with the right
    # row count is the failure mode that defeated the resume check in July
    # 2026, so assert on CONTENT, not on file existence - and don't write the
    # four output files (or promote the tarball out of `.partial`) until
    # that content has been validated, or a failed shard leaves a
    # complete-looking hashes/fp/missing/shards set for a resume check to
    # wrongly trust.
    assert len(hashes) + len(missing) == len(mine), (
        f"accounting mismatch: {len(hashes)} parsed + {len(missing)} missing "
        f"!= {len(mine)} manifest rows")
    if len(mine) and not hashes:
        raise SystemExit(f"FATAL: shard {args.shard} parsed 0 of {len(mine):,d} structures")
    if len(mine) and partial_tar_path.stat().st_size < 200:
        raise SystemExit(f"FATAL: shard tarball {partial_tar_path} is empty")

    # Write the three small metadata files FIRST, and promote the tarball out
    # of `.partial` LAST. Downstream stages glob on `shard_NN.tar.gz` as the
    # de-facto completion marker, so a crash between these lines must never
    # leave a promoted tarball with no hashes/fp/missing companion - the same
    # complete-looking-output shape the validate-before-writing change above
    # already set out to eliminate, just moved to a different pair of lines.
    pd.DataFrame(hashes, columns=["design_id", "coord_hash", "protein_chain", "n_res"]).to_csv(
        out_dir / "hashes" / f"{shard_name}.tsv", sep="\t", index=False)
    pd.DataFrame(fps, columns=FP_COLUMNS).to_parquet(
        out_dir / "fp" / f"{shard_name}.parquet", index=False)
    pd.DataFrame(missing, columns=["design_id", "reason"]).to_csv(
        out_dir / "missing" / f"{shard_name}.tsv", sep="\t", index=False)
    os.replace(partial_tar_path, tar_path)

    print(f"done in {time.time() - t0:.0f}s: {len(hashes):,d} parsed, "
          f"{len(seen_local):,d} unique backbones, {len(missing):,d} missing", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
