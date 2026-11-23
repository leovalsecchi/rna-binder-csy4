#!/usr/bin/env python3
"""gpubench_v1 REFOLD -- an identical 600-design batch per stack, all 3 backends.

WHY POOLS: AF3 and ESM load the model ONCE and run over --input_dir
(refold_af3_scratch.sh:262-290). That is where their 10-13x speedup lives, so
per-design cost depends on how many designs share one model load. The campaign's
batch was one uuid (~600 designs). HD, though, makes only 150 designs per uuid
and needs 4 uuids to reach 600 -- refolding per-uuid would hand HD a 150-design
batch and everyone else 600, inflating HD's per-design refold cost for a reason
that has nothing to do with HD.

So each stack gets a pool of exactly N designs, round-robined across its bench
uuids, used as the refold `source_run`. Every stack then gets ONE AF3 + ONE B2 +
ONE ESM job over the same batch size, and the numbers are directly comparable.

THE POOL IS A PACKED SOURCE, NOT A SYMLINK DIR (corrected 2026-11-11): a finished
gen run has no `designs/` -- the finalize tail packs it to `designs.tar.zst` and
deletes the tree (see bench_lib.py). So the pool is built by extracting the
chosen members out of each run's tarball and REPACKING them as one
`pools/<stack>/designs.tar.zst`, which the refold worker already knows how to
read (refold_af3_scratch.sh:60-63 unpacks a packed source_run to node-local
scratch). Inode cost: 2 files per stack instead of 600 symlinks, and nothing is
written to Lustre per design.

BG+BIF'S TRUE BACKBONE IS APPLIED HERE, during staging, so production gen output
is never mutated -- see bif_true_designs.py for what and why.

WHY WE CALL launch.py's BUILDERS INSTEAD OF ITS CLI: launch.py has standalone
branches for config_kind refold_af3 (line 1356) and refold_b2 (1426) but NOT
refold_esm -- an ESM config falls through to the generator path and dies on
KeyError: 'generator'. All three builder functions exist, though, so we import
them and drive all three backends through one uniform path. That reuses the exact
production sbatch text rather than reimplementing it.

Deviations from the stock refold configs, all required:
  * stages/refold_{af3,b2,esm}.sh -> *_scratch.sh. Every campaign driver does
    this rewrite before submit (mv8_refold_bulk.py:227), so the _scratch
    variants are what actually ran. Node-local scratch also keeps the bench off
    the Lustre inode quota.
  * slurm.nodelist cleared. The shipped configs pin eu-g4-[001-032] and
    eu-g4-032 was DECOMMISSIONED -- submitting with it fails "Invalid node name".
    We pin by card model instead (--gpus=rtx_4090:1, which is eu-g6 only).
  * slurm.gres cleared: the card model already fixes VRAM, and gpumem is a
    discrete match key rather than a >= threshold.
  * --requeue / --open-mode stripped. The AF3 worker skips already-done designs,
    so a requeued attempt measures only the remainder -- a corrupt clock.
  * --tmp=30G and --exclude, copied from the campaign driver.

NOT changed, deliberately: B2's --ntasks-per-node=4. The campaign never fixed
that one (only HD's --ntasks), so leaving it reproduces campaign conditions. It
is identical across stacks, so it cannot affect the ranking.

Usage:
  python3 bench_refold.py                    # build pools, print what would go
  python3 bench_refold.py --submit
"""
from __future__ import annotations

import argparse
import csv
import os
import pathlib
import re
import shutil
import subprocess
import sys

import yaml

RB = pathlib.Path("/cluster/work/wicky/lvalsecchi/rna_binder")
CAMP = RB / "campaigns" / "csy4_crrna_free"
BENCH = RB / "bench_gpueff"

sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(RB / "code" / "scripts"))
import bench_lib as BL  # noqa: E402  (path must be set first)
import launch as L  # noqa: E402

PART = "gpuhe.24h"
GPUS = "rtx_4090:1"
WALL = "23:30:00"          # AF3 TIMEOUTs at 4h; headroom, and non-preemptible
TMP = "30G"
EXCLUDE = "eu-g6-064"      # campaign driver's permanent exclude, inside our pool

STACKS = {
    "rfd3_limpnn": "rfd3__ligandmpnn",
    "bg_bif": "boltzgen__boltzif",
    "bg_limpnn": "boltzgen__ligandmpnn",
    "bd1_limpnn": "boltzdesign__ligandmpnn",
    "hd_limpnn": "halludesign__ligandmpnn",
}
BACKENDS = {
    "af3": ("refold__af3", L.build_refold_af3_sbatch),
    "b2": ("refold__b2", L.build_refold_b2_sbatch),
    "esm": ("refold__esm", L.build_refold_esm_sbatch),
}


# The worker extracts each design's sequence with BioPython inside boltzgen.sif
# (refold_af3_scratch.sh:118-120), NOT with biotite and NOT on the host. Its
# MMCIFParser hard-requires _atom_site.B_iso_or_equiv, which biotite omits unless
# the annotation was registered -- and a pool whose designs it cannot read dies in
# stage 1 after burning the queue wait. That happened once (all three BG+BIF
# refolds, 2026-11-12), so every pool is now sampled with the real parser in the
# real container before it is packed.
PROBE_SIF = "/cluster/work/wicky/containers/boltzgen.sif"
PROBE = BENCH / "probe_biopython.py"


def verify_parsable(ddir: pathlib.Path, k: int = 10) -> None:
    files = sorted(p for p in ddir.iterdir() if p.is_file())[:k]
    if not files or not PROBE.is_file():
        print(f"  !! cannot verify {ddir}: no files or missing {PROBE.name}")
        raise SystemExit(3)
    r = subprocess.run(
        ["apptainer", "exec", "--bind", "/cluster/work/wicky",
         "--bind", "/cluster/scratch", PROBE_SIF, "python3", str(PROBE),
         *[str(f) for f in files]],
        capture_output=True, text=True)
    bad = [ln for ln in r.stdout.splitlines() if not ln.startswith("OK")]
    if r.returncode != 0 or bad:
        print(f"  !! {ddir}: designs are NOT readable by the refold worker")
        for ln in (bad or r.stderr.splitlines())[:5]:
            print(f"     {ln}")
        raise SystemExit(3)
    print(f"  parse check: {len(files)}/{len(files)} readable by "
          f"{pathlib.Path(PROBE_SIF).name} BioPython")


def stage_root() -> pathlib.Path:
    """Node-local scratch when we have it, so staging costs no Lustre inodes."""
    t = os.environ.get("TMPDIR")
    if t and pathlib.Path(t).is_dir():
        return pathlib.Path(t) / "gpubench_stage"
    return pathlib.Path(f"/cluster/scratch/{os.environ.get('USER','lvalsecchi')}") \
        / "gpubench_stage"


def build_pool(stack: str, sub: str, n: int, label: str):
    """Round-robin n designs across this stack's FINISHED runs, then repack them
    as one packed source_run the refold worker can read."""
    runs = BL.finished_runs(CAMP / sub, label)
    lost = BL.unfinished_runs(CAMP / sub, label)
    per_run = [(r, BL.members(r)) for r in runs]
    per_run = [(r, m) for r, m in per_run if m]
    if not per_run:
        raise SystemExit(
            f"{stack}: no finished run under {CAMP / sub} for label {label}"
            f"{' (unfinished: ' + ' '.join(r.name for r in lost) + ')' if lost else ''}")

    picked: dict[pathlib.Path, list[str]] = {}
    total = i = 0
    while total < n and any(i < len(m) for _, m in per_run):
        for r, m in per_run:
            if i < len(m):
                picked.setdefault(r, []).append(m[i])
                total += 1
                if total == n:
                    break
        i += 1

    stage = stage_root() / stack
    shutil.rmtree(stage, ignore_errors=True)
    for r, want in picked.items():
        BL.extract(r, want, stage)

    if stack == "bg_bif":
        import bif_true_designs as BT
        BT.rewrite_dir(stage / "designs")
    verify_parsable(stage / "designs")

    pool = BENCH / "pools" / stack
    shutil.rmtree(pool, ignore_errors=True)
    got = BL.pack(stage, pool)
    shutil.rmtree(stage, ignore_errors=True)
    return pool, got, [r.name.split("_")[-1] for r in runs], \
        [r.name.split("_")[-1] for r in lost]


def bench_rewrites(txt: str) -> str:
    txt = re.sub(r"(refold_(?:af3|b2|esm))\.sh", r"\1_scratch.sh", txt)
    txt = re.sub(r"^#SBATCH --requeue\n", "", txt, flags=re.M)
    txt = re.sub(r"^#SBATCH --open-mode=append\n", "", txt, flags=re.M)
    if "--tmp=" not in txt:
        txt = re.sub(r"(^#SBATCH --partition=.*\n)", rf"\1#SBATCH --tmp={TMP}\n",
                     txt, count=1, flags=re.M)
    if "--exclude=" not in txt:
        txt = re.sub(r"(^#SBATCH --partition=.*\n)",
                     rf"\1#SBATCH --exclude={EXCLUDE}\n", txt, count=1, flags=re.M)
    return txt


def make_sbatch(be: str, pool: pathlib.Path, label: str) -> tuple[str, pathlib.Path]:
    sub, builder = BACKENDS[be]
    cfg = yaml.safe_load((RB / "code" / "configs" / sub / "default.yaml").read_text())
    cfg = L.deep_merge(cfg, L.parse_overrides([
        f"slurm.partition={PART}", f"slurm.gpus={GPUS}", f"slurm.time={WALL}",
        "slurm.nodelist=", "slurm.gres=",
    ]))
    cfg["source_run"] = str(pool)
    rdir = pool / f"refold_{be}" / label
    return bench_rewrites(builder(cfg, rdir, pool)), rdir


def submit(stack: str, be: str, pool: pathlib.Path, label: str, do_submit: bool):
    txt, rdir = make_sbatch(be, pool, label)

    # Refuse anything that would corrupt the clock or write to Lustre per design.
    bad = [w for w in ("--requeue", "--open-mode") if w in txt]
    if "_scratch.sh" not in txt:
        bad.append("non-scratch worker")
    if bad:
        print(f"  !! {stack} {be}: rewrite failed ({', '.join(bad)}) -- skipping")
        return None

    if not do_submit:
        print(f"  DRY {stack:<13} {be:<4} -> {rdir}")
        for ln in txt.splitlines():
            if ln.startswith("#SBATCH") or "_scratch.sh" in ln:
                print(f"       {ln.strip()[:150]}")
        return None

    if rdir.exists():
        shutil.rmtree(rdir)
    rdir.mkdir(parents=True)
    sb = rdir / "run.sbatch"
    sb.write_text(txt)
    sb.chmod(0o755)

    p = subprocess.run(["sbatch", "--parsable", str(sb)],
                       capture_output=True, text=True)
    jid = p.stdout.strip()
    if not jid.isdigit():
        print(f"  !! {stack} {be}: sbatch failed: {p.stdout}{p.stderr}")
        return None
    print(f"  {stack:<13} {be:<4} jobid={jid}")
    return jid, str(rdir)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="gpubench_v1")
    ap.add_argument("-n", type=int, default=600,
                    help="designs per stack (identical across stacks)")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--only", default="",
                    help="comma-separated stacks to (re)build and submit, e.g. "
                         "bg_bif -- for resubmitting one stack without "
                         "disturbing pools that already refolded")
    a = ap.parse_args()
    want = [s for s in a.only.split(",") if s] or list(STACKS)
    bad = [s for s in want if s not in STACKS]
    if bad:
        ap.error(f"unknown stack(s): {' '.join(bad)}")

    BENCH.mkdir(parents=True, exist_ok=True)
    print(f"=== pools: {a.n} designs per stack  (stage: {stage_root()}) ===")
    pools, short = {}, []
    for stack in want:
        sub = STACKS[stack]
        pool, got, uuids, lost = build_pool(stack, sub, a.n, a.label)
        print(f"  {stack:<13} {got:>4}/{a.n} from {len(uuids)} finished uuid(s): "
              f"{' '.join(uuids) or '-'}")
        if lost:
            print(f"                 !! unfinished/lost, excluded: {' '.join(lost)}")
        pools[stack] = pool
        if got != a.n:
            short.append((stack, got))

    if short:
        print("\n!! stacks short of n -- gen incomplete or a job was lost:")
        for s, got in short:
            print(f"   {s}: {got}")
        print("   Unequal batch sizes defeat the comparison. Fix gen, re-run.")
        if a.submit:
            sys.exit(2)

    print(f"\n=== refold jobs ({'SUBMIT' if a.submit else 'DRY'}) ===")
    rows = []
    for stack in want:
        for be in BACKENDS:
            got = submit(stack, be, pools[stack], a.label, a.submit)
            if got:
                rows.append({"role": "refold", "stack": stack, "backend": be,
                             "cfg": BACKENDS[be][0], "uuid8": "pool",
                             "jobid": got[0], "run_dir": got[1],
                             "n_designs": a.n})

    if rows:
        man = BENCH / "refold_manifest.tsv"
        new = not man.exists()
        with man.open("a", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
            if new:
                w.writeheader()
            w.writerows(rows)
        print(f"\n{len(rows)} job(s) -> {man}")
    elif not a.submit:
        print("\nDRY RUN -- nothing submitted.")


if __name__ == "__main__":
    main()
