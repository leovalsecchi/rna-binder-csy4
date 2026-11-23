#!/usr/bin/env python3
"""gpubench_v1 collector: GPU-h per design, per stack x role.

Reads the submit manifests, asks sacct for each job's real elapsed time, counts
the designs that job produced or refolded, and emits GPU-h/design.

Correctness rules, each one a scar from the campaign:
  * sacct -D -X. Without -D requeued attempts VANISH and the surviving row
    under-reports. With it, a job that ran twice shows up twice -- which we want
    to SEE, not quietly average.
  * Any job not COMPLETED, or with more than one attempt, is REJECTED rather than
    down-weighted. A preempted or timed-out clock is a wrong clock, not a slow one.
  * GPU-h = ElapsedRaw x gres/gpu parsed from AllocTRES. Bench jobs are all 1 GPU,
    but assert it rather than assume it.
  * The BoltzGen gen job feeds BG+BIF and BG+LiMPNN from ONE allocation. Reported
    as its own `bg_shared` row and split three ways: the measured stage split
    (from output-file mtimes, since the stage banners carry no timestamps), plus
    the 50/50 and full-charge-each readings derived from it.

Usage:  python3 collect_gpuh.py [--label gpubench_v1]
"""
from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import re
import statistics as st
import subprocess
import sys

RB = pathlib.Path("/cluster/work/wicky/lvalsecchi/rna_binder")
CAMP = RB / "campaigns" / "csy4_crrna_free"
BENCH = RB / "bench_gpueff"

sys.path.insert(0, str(BENCH))
import bench_lib as BL  # noqa: E402  (designs live in designs.tar.zst)

# gen stack -> campaign subdirs holding its designs. bg_shared maps to BOTH
# inverse folders on purpose: one allocation, two products.
GEN_DIRS = {
    "rfd3_limpnn": ["rfd3__ligandmpnn"],
    "bg_shared": ["boltzgen__boltzif", "boltzgen__ligandmpnn"],
    "bd1_limpnn": ["boltzdesign__ligandmpnn"],
    "hd_limpnn": ["halludesign__ligandmpnn"],
}
DESIGN_EXT = (".cif", ".pdb", ".cif.gz")
SKIP = re.compile(r"sequences\.fa|\.npz$|\.json$|\.csv$")


def sacct(jobids: list[str]) -> dict[str, list[dict]]:
    if not jobids:
        return {}
    cmd = ["sacct", "-D", "-X", "-n", "-P", "-j", ",".join(sorted(set(jobids))),
           "--format=JobID,JobName%60,Partition,State,ElapsedRaw,Start,End,"
           "NodeList%30,AllocTRES%140"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"sacct failed: {r.stderr.strip()}")
    rows: dict[str, list[dict]] = collections.defaultdict(list)
    for line in r.stdout.splitlines():
        f = line.split("|")
        if len(f) < 9:
            continue
        rows[f[0].split(".")[0]].append({
            "jobid": f[0], "name": f[1], "part": f[2], "state": f[3],
            "elapsed": int(f[4]) if f[4].isdigit() else 0,
            "start": f[5], "end": f[6], "node": f[7], "tres": f[8],
        })
    return rows


def n_gpus(tres: str) -> int:
    m = re.search(r"gres/gpu=(\d+)", tres)
    return int(m.group(1)) if m else 0


def count_designs(stack: str, uuid8: str, label: str) -> int:
    """Designs a gen job actually produced, read from the tarball sidecar.

    A finished run has no `designs/` -- the finalize tail packs it and deletes
    the tree. The sidecar is written by that same tail, so counting it also
    means an unfinished run counts 0 and gets REJECTED rather than averaged in
    with a partial denominator. That matters: this bench lost an HD job to a
    JIT-lock hang, and its run dir still holds one stray file."""
    total = 0
    for sub in GEN_DIRS.get(stack, []):
        for rdir in (CAMP / sub).glob(f"{label}_{uuid8}"):
            if BL.is_finished(rdir):
                total += len(BL.members(rdir))
    return total


# A clean HD job reaches DeepSpeed in 6-10 min (measured: 6, 9, 10). A job that
# loses the race for the SHARED evoformer_attn JIT lock sits there for HOURS
# (measured: 119 min) -- see [[rna-binder-hd-evoformer-jit-fix]]. That wall time
# is a WRONG clock, not a slow one, and sacct cannot see the difference: the job
# still reports COMPLETED on 1 GPU with a full design count. So the startup gap
# is measured from the log and such a job is rejected, the same way a requeued
# or preempted one is.
HD_STARTUP_REJECT_S = 20 * 60


def hd_startup_s(stack: str, uuid8: str, label: str,
                 start_epoch: float | None) -> float | None:
    """Seconds from job start to DeepSpeed init, for HD gen jobs.

    The start comes from SACCT, not from the log's own "=== Job started:" banner:
    that banner is `date`'s default format ("Wed Nov 11 01:26:34 PM CET 2026"),
    which `datetime.fromisoformat` cannot parse -- it returns None, and an
    earlier version of this check therefore measured nothing at all while
    appearing to work. DeepSpeed's own "[2026-11-11 16:10:05,xxx]" stamp IS
    isoformat-compatible.
    """
    if stack != "hd_limpnn" or not start_epoch:
        return None
    for sub in GEN_DIRS.get(stack, []):
        for rdir in (CAMP / sub).glob(f"{label}_{uuid8}"):
            f = rdir / "slurm.log"
            if not f.is_file():
                continue
            with f.open(errors="replace") as fh:
                for line in fh:
                    if line.startswith("[") and "ds_accelerator" in line:
                        t1 = epoch(line[1:].split(",", 1)[0])
                        return (t1 - start_epoch) if t1 else None
    return None


# Stage banners in a timestamped BoltzGen log, e.g.
#   [10:39:40] === Stage 2: BoltzGen BIF (inverse fold + analysis) ===
BG_STAGE = re.compile(r"^\[(\d\d):(\d\d):(\d\d)\] === Stage (\d)")


def bg_stage_seconds(label: str) -> dict[str, float] | None:
    """Seconds per BoltzGen stage, from a run submitted with TS_STAGES=1.

    The stock worker prints its stage banners with NO timestamps, so the split
    cannot be recovered after the fact (and both inverse folders' designs are
    packed at job end with identical mtimes). bench_gen.sh's TS_STAGES=1 wraps
    the worker's stdout in a timestamper; this reads that log.

    Returns seconds for backbone / bif / limpnn. The backbone stage is SHARED;
    the other two are each one inverse folder's own work.
    """
    for rdir in sorted((CAMP / "boltzgen__boltzif").glob(f"{label}_*")):
        for log in sorted(rdir.glob("slurm*.log")):
            marks, rc_t = {}, None
            with log.open(errors="replace") as fh:
                for line in fh:
                    m = BG_STAGE.match(line)
                    if m:
                        h, mi, sec, st = (int(x) for x in m.groups())
                        marks.setdefault(st, h * 3600 + mi * 60 + sec)
                    elif line.startswith("[") and "__RC=" in line:
                        h, mi, sec = (int(x) for x in line[1:9].split(":"))
                        rc_t = h * 3600 + mi * 60 + sec
            if {1, 2, 3} <= set(marks) and rc_t:
                return {"backbone": marks[2] - marks[1],
                        "bif": marks[3] - marks[2],
                        "limpnn": rc_t - marks[3],
                        "total": rc_t - marks[1],
                        "log": str(log)}
    return None


def _max_mtime(paths) -> float | None:
    ms = [p.stat().st_mtime for p in paths if p.exists()]
    return max(ms) if ms else None


def bg_stage_split(uuid8: str, label: str, start_epoch: float,
                   end_epoch: float) -> dict[str, float] | None:
    """Seconds split into the genuinely SHARED part of a BoltzGen gen job and
    the part that serves the two inverse folders.

    boltzgen_multif.sh prints its stage banners with no timestamps, so the log
    cannot date them, and production is left untouched rather than patched. The
    outputs can date one boundary: the `_backbones/BG_len80_<uuid8>/spec_*`
    files are written when backbone design ends, and they survive on Lustre
    unpacked.

    That is the ONLY boundary left. The finer split (boltzif vs ligandmpnn) is
    not recoverable: `bif_native/` is not published, and both inverse folders'
    `designs/` trees are packed into tarballs written at the very end of the
    job, so their mtimes are identical. So this reports the shared backbone
    stage and one combined inverse-folding remainder -- an honest two-way split
    instead of a fabricated three-way one. The caller still prints 50/50 and
    full-charge-each alongside it.
    """
    bb = CAMP / "_backbones" / f"BG_len80_{uuid8}"
    if not bb.is_dir():
        return None
    t1 = _max_mtime(list(bb.glob("spec_*.pdb")) or list(bb.glob("spec_*.cif")))
    if not t1 or not (start_epoch < t1 < end_epoch):
        return None  # mtimes inconsistent with the job window; say nothing
    return {
        "1_backbone_design(SHARED)": t1 - start_epoch,
        "2_inverse_folding(bif+limpnn, unsplittable)": end_epoch - t1,
    }


def epoch(s: str) -> float | None:
    import datetime as dt
    try:
        return dt.datetime.fromisoformat(s).timestamp()
    except ValueError:
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="gpubench_v1")
    ap.add_argument("--out", default=str(BENCH / "gpuh_per_design.tsv"))
    a = ap.parse_args()

    mans = sorted(BENCH.glob("*_manifest.tsv"))
    if not mans:
        sys.exit(f"no *_manifest.tsv in {BENCH} -- nothing submitted yet?")
    entries = []
    for m in mans:
        with m.open() as f:
            entries += [r for r in csv.DictReader(f, delimiter="\t") if r.get("jobid")]
    if not entries:
        sys.exit("manifests are empty")

    att = sacct([e["jobid"] for e in entries])
    recs, rejected = [], []
    for e in entries:
        jid, role, stack = e["jobid"], e["role"], e["stack"]
        be = e.get("backend", "-") or "-"
        tries = att.get(jid, [])
        if not tries:
            rejected.append((jid, role, stack, be, "absent from sacct"))
            continue
        if len(tries) > 1:
            rejected.append((jid, role, stack, be, f"{len(tries)} attempts (requeued)"))
            continue
        t = tries[0]
        if t["state"] != "COMPLETED":
            rejected.append((jid, role, stack, be, t["state"]))
            continue
        g = n_gpus(t["tres"])
        if g != 1:
            rejected.append((jid, role, stack, be, f"gres/gpu={g}"))
            continue
        n = int(e["n_designs"]) if e.get("n_designs") else \
            count_designs(stack, e["uuid8"], a.label)
        if n <= 0:
            rejected.append((jid, role, stack, be, "0 designs packed"))
            continue
        su = hd_startup_s(stack, e["uuid8"], a.label,
                          epoch(t["start"])) if role == "gen" else None
        if su is not None and su > HD_STARTUP_REJECT_S:
            rejected.append((jid, role, stack, be,
                             f"{su/60:.0f} min to DeepSpeed init "
                             f"(evoformer JIT lock; clock unusable)"))
            continue
        gpuh = t["elapsed"] * g / 3600.0
        recs.append({"role": role, "stack": stack, "backend": be,
                     "uuid8": e["uuid8"], "jobid": jid, "node": t["node"],
                     "elapsed_s": t["elapsed"], "gpuh": round(gpuh, 4),
                     "n_designs": n,
                     "gpuh_per_design": round(gpuh / n, 6),
                     "startup_s": round(su) if su is not None else ""})

    if not recs:
        print("no accepted jobs yet.")
        for r in rejected:
            print("   REJECTED", *r)
        return

    outp = pathlib.Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)
    print(f"{len(recs)} accepted job(s) -> {outp}")

    if rejected:
        print(f"\n!! {len(rejected)} REJECTED (deliberately NOT averaged in):")
        for jid, role, stack, be, why in rejected:
            print(f"   {jid:>12}  {role:<6} {stack:<13} {be:<4} {why}")

    by = collections.defaultdict(list)
    for r in recs:
        by[(r["stack"], r["role"], r["backend"])].append(r["gpuh_per_design"])
    print(f"\n{'stack':<13} {'role':<6} {'be':<4} {'n':>2} {'GPU-h/design':>13}"
          f" {'min':>10} {'max':>10}")
    print("-" * 64)
    for (stack, role, be), v in sorted(by.items()):
        print(f"{stack:<13} {role:<6} {be:<4} {len(v):>2} {st.median(v):>13.6f}"
              f" {min(v):>10.6f} {max(v):>10.6f}")

    # per-stack totals: gen + all three refold backends
    gen = {s: st.median(v) for (s, r, _b), v in by.items() if r == "gen"}
    ref = collections.defaultdict(float)
    for (s, r, _b), v in by.items():
        if r == "refold":
            ref[s] += st.median(v)
    print("\n--- per-design totals (gen + af3 + b2 + esm) ---")
    gen_for = {"rfd3_limpnn": "rfd3_limpnn", "bg_bif": "bg_shared",
               "bg_limpnn": "bg_shared", "bd1_limpnn": "bd1_limpnn",
               "hd_limpnn": "hd_limpnn"}
    # BoltzGen: charge each inverse folder the FULL shared backbone stage plus
    # its OWN inverse folding -- what that stack would cost run on its own, which
    # is the question "which stack is most efficient" actually asks. A 50/50 of
    # the whole job is wrong twice over: it credits BG with a backbone saving that
    # only exists if you want both products, and it splits an inverse-folding
    # stage that is 6.8:1 in BIF's favour straight down the middle.
    # gen["bg_shared"] is per-design over all 1200 designs, so the job's own GPU-h
    # is gen*1200 and a stack's own 600 designs give: 2 * fraction * gen.
    st_s = bg_stage_seconds(a.label)
    bg_frac = None
    if st_s:
        bg_frac = {"bg_bif": (st_s["backbone"] + st_s["bif"]) / st_s["total"],
                   "bg_limpnn": (st_s["backbone"] + st_s["limpnn"]) / st_s["total"]}
    rows = []
    for stack, gkey in gen_for.items():
        if gkey not in gen or stack not in ref:
            continue
        g = gen[gkey]
        note = ""
        if gkey == "bg_shared":
            if bg_frac:
                g = 2.0 * bg_frac[stack] * g
                note = f"  (full shared backbone + own inverse fold)"
            else:
                note = "  (gen 50/50 -- no timestamped BG run to split on)"
        rows.append((stack, g, ref[stack], g + ref[stack], note))
    if rows:
        base = min(r[3] for r in rows)
        print(f"{'stack':<13} {'gen':>10} {'refold':>10} {'total':>10} {'xbest':>7}")
        for stack, g, rf, tot, note in sorted(rows, key=lambda x: x[3]):
            print(f"{stack:<13} {g:>10.6f} {rf:>10.6f} {tot:>10.6f}"
                  f" {tot/base:>7.2f}{note}")

    if st_s:
        t = st_s["total"]
        print("\n--- BoltzGen stage split (timestamped run) ---")
        print(f"  1 backbone design  SHARED        {st_s['backbone']:>6.0f} s "
              f"({100*st_s['backbone']/t:4.1f}%)")
        print(f"  2 BIF inverse fold + fold        {st_s['bif']:>6.0f} s "
              f"({100*st_s['bif']/t:4.1f}%)")
        print(f"  3 LigandMPNN over the backbones  {st_s['limpnn']:>6.0f} s "
              f"({100*st_s['limpnn']/t:4.1f}%)")
        print(f"  BIF costs {st_s['bif']/max(st_s['limpnn'],1):.1f}x LigandMPNN; "
              f"source {st_s['log']}")

    for r in recs:
        if st_s:
            break   # superseded by the timestamped split printed just above
        if r["stack"] == "bg_shared" and r["role"] == "gen":
            se = epoch(att[r["jobid"]][0]["start"])
            ee = epoch(att[r["jobid"]][0]["end"])
            sp = bg_stage_split(r["uuid8"], a.label, se, ee) if (se and ee) else None
            tag = f"bg_shared gen {r['uuid8']} ({r['gpuh']:.3f} GPU-h)"
            if sp:
                tot = sum(sp.values()) or 1
                print(f"\n{tag} measured stage split (from output mtimes):")
                for k, v in sp.items():
                    print(f"   {k:<44} {v/3600:7.3f} h  ({100*v/tot:5.1f}%)")
                print(f"   -> 50/50 {r['gpuh_per_design']:.6f} GPU-h/design each,"
                      f" full-charge-each {2*r['gpuh_per_design']:.6f}")
            else:
                print(f"\n{tag}: stage split unavailable -> "
                      f"50/50 {r['gpuh_per_design']:.6f} GPU-h/design each, "
                      f"full-charge-each {2*r['gpuh_per_design']:.6f}")


if __name__ == "__main__":
    main()
