#!/usr/bin/env python3
"""Shared helpers for the gpubench scripts: read designs out of the packed
`designs.tar.zst` the gen pipeline publishes.

WHY THIS EXISTS (discovered 2026-11-11, after the gen jobs finished): a FINISHED
gen run has no `designs/` directory at all. The finalize tail packs it into
`designs.tar.zst` + `.tar.zst.list` and deletes the tree -- that is the
inode-quota discipline the whole campaign runs under, and the refold worker
already knows how to read it (refold_af3_scratch.sh:60-63 unpacks a packed
source_run to node-local scratch via archive_pipeline.sh's unpack_unit). The
bench scripts were written against an unpacked `designs/`, so they all have to
come through here instead.

The tarball doubles as the completeness signal, because only the finalize tail
publishes it. A cancelled or still-running job leaves a partial `designs/` and
NO tarball -- so keying on the tarball excludes exactly the runs that must never
enter a measurement. That is not incidental: this bench lost an HD job to a
DeepSpeed JIT-lock hang, and its half-written run dir must not leak into a pool.

Pack/unpack mirror archive_pipeline.sh:158,179 (`tar -c | zstd` / `zstd -d |
tar -x`) so the files are readable by the production readers. Level 3 not 19:
a pool is transient, and -19 buys nothing but minutes.
"""
from __future__ import annotations

import pathlib
import subprocess

DESIGN_EXT = (".cif", ".pdb", ".cif.gz")
TAR = "designs.tar.zst"


def is_finished(run: pathlib.Path) -> bool:
    """A run is usable iff the finalize tail published its tarball + sidecar."""
    return (run / TAR).is_file() and (run / f"{TAR}.list").is_file()


def finished_runs(parent: pathlib.Path, label: str) -> list[pathlib.Path]:
    return sorted(r for r in parent.glob(f"{label}_*") if is_finished(r))


def unfinished_runs(parent: pathlib.Path, label: str) -> list[pathlib.Path]:
    """Reported, never silently skipped -- a missing run means a lost job."""
    return sorted(r for r in parent.glob(f"{label}_*") if not is_finished(r))


def members(run: pathlib.Path) -> list[str]:
    """Design members inside the tarball ("designs/<name>"), in sidecar order.

    Filtered on extension ALONE, deliberately. Every stack names its designs
    differently -- RFD3/BG/BD1 use `_design<D>_seq<S>`, HD uses
    `_recycle<N>_sample<M>_pack1` -- so any name-shaped filter silently drops a
    whole stack (an earlier `"_design" in ln` test zeroed out HD entirely). The
    sidecars are pure: verified 640/640, 600/600 and 150/150 design files with
    nothing else in them, because the finalize tail canonicalises the tree
    before packing it."""
    return [ln.strip() for ln in (run / f"{TAR}.list").read_text().splitlines()
            if ln.strip().endswith(DESIGN_EXT)]


def extract(run: pathlib.Path, want: list[str], dest: pathlib.Path) -> None:
    """Extract exactly `want` from the run's tarball into dest/ (-> dest/designs/)."""
    if not want:
        return
    dest.mkdir(parents=True, exist_ok=True)
    lst = dest / "_members.txt"
    lst.write_text("\n".join(want) + "\n")
    subprocess.run(
        f"zstd -dq -c {run / TAR} | tar -x -C {dest} --files-from {lst}",
        shell=True, check=True)
    lst.unlink()


def pack(stage: pathlib.Path, out_run: pathlib.Path) -> int:
    """Pack stage/designs -> out_run/designs.tar.zst + .list. Returns n members."""
    d = stage / "designs"
    if not d.is_dir():
        raise ValueError(f"no staged designs at {d} -- nothing was extracted")
    names = sorted(p.name for p in d.iterdir() if p.is_file())
    if not names:
        raise ValueError(f"nothing to pack in {stage}/designs")
    out_run.mkdir(parents=True, exist_ok=True)
    lst = stage / "_pack.txt"
    lst.write_text("".join(f"designs/{n}\n" for n in names))
    subprocess.run(
        f"cd {stage} && tar -c --files-from {lst} | zstd -q -3 -T0 -o {out_run / TAR}",
        shell=True, check=True)
    subprocess.run(["zstd", "-tq", str(out_run / TAR)], check=True)
    (out_run / f"{TAR}.list").write_text(lst.read_text())
    lst.unlink()
    return len(names)
