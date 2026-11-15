#!/usr/bin/env python3
# refold_gap_fill.py — fill missing AF3/ESM refolds for ANY campaign stack/run,
# DESIGN-LEVEL (folds only the gaps), array-native, on 24G bulk.
#
# Gaps are derived LIVE from the .list sidecars (current on-disk truth, no stale
# collect CSV): for each uuid, gap = design stems - refold stems, per backend.
# Reuses the REFOLD_ALLOWLIST-patched workers + separate "_gap" leaves (collect.py
# dedups refolds per design across leaves, so gaps + originals merge cleanly).
#
# Placement (both backends on 24G — ESM OOMs on 11G 2080s): gpu.bulk pinned
# eu-g3 (24G Turing) + gpuhe.bulk gpumem:24g (Ampere/4090; 24g excludes 2080).
# Reached via sbatch --hold -> scontrol move partition/qos -> release.
#
#   module load stack/2024-06 python/3.12.8
#   python3 refold_gap_fill.py --stacks rfd3__ligandmpnn --run midrun_v3 --backends esm,af3
#   python3 refold_gap_fill.py --stacks rfd3__ligandmpnn --run midrun_v3 --backends esm --submit
import sys, re, argparse, shlex, subprocess, glob, os, collections
from datetime import datetime
from pathlib import Path

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/scripts")
import yaml
import launch

CAMP = Path(RB) / "campaigns" / "csy4_crrna_free"
POOL_G3 = "eu-g3-[062-077],eu-g3-[083-084]"
# lane = (submit_part, submit_qos, wall, nodelist, gres, bulk_part, bulk_qos)
BACKEND_LANES = {
    "af3": [
        ("gpu.120h",   "es_wicky/gpu/120",   "48:00:00", POOL_G3, "gpu:1",      "gpu.bulk",   "es_wicky/gpu/bulk"),
        ("gpuhe.120h", "es_wicky/gpuhe/120",  "48:00:00", None,    "gpumem:24g", "gpuhe.bulk", "es_wicky/gpuhe/bulk"),
    ],
    "esm": [
        ("gpu.4h",   "es_wicky/gpu/4",   "04:00:00", POOL_G3, "gpu:1",      "gpu.bulk",   "es_wicky/gpu/bulk"),
        ("gpuhe.4h", "es_wicky/gpuhe/4",  "04:00:00", None,    "gpumem:24g", "gpuhe.bulk", "es_wicky/gpuhe/bulk"),
    ],
}
TAG = {"af3": "AF3", "esm": "ESM", "b2": "B2"}

ap = argparse.ArgumentParser()
ap.add_argument("--stacks", required=True, help="comma list of stack dir names")
ap.add_argument("--run", required=True)
ap.add_argument("--backends", default="af3,esm")
ap.add_argument("--submit", action="store_true")
ap.add_argument("--limit", type=int, default=0, help="cap uuids PER BACKEND")
ap.add_argument("--order", choices=["largest", "smallest"], default="largest")
args = ap.parse_args()
STACKS = [s.strip() for s in args.stacks.split(",") if s.strip()]
RUN = args.run
BACKENDS = [b.strip() for b in args.backends.split(",") if b.strip()]


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True)


def design_stems(list_path):
    s = set()
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        if "_ESM_" in b or "_AF3_" in b or "_B2_" in b:
            continue
        for e in (".cif.gz", ".cif", ".pdb"):
            if b.endswith(e):
                s.add(b[:-len(e)]); break
    return s


def refold_stems(list_path, tag):
    s = set()
    if not os.path.exists(list_path):
        return s
    marker = f"_{tag}_"
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        i = b.find(marker)
        if i > 0:
            s.add(b[:i])
    return s


# --- LIVE gap detection from .list sidecars -------------------------------------
gap = {b: collections.defaultdict(list) for b in BACKENDS}
for stack in STACKS:
    pool = CAMP / stack / RUN
    for dl in glob.glob(str(pool / "designs" / "*.tar.zst.list")):
        m = re.search(r"([0-9a-f]{8})\.tar", dl)
        if not m:
            continue
        u = m.group(1)
        ds = design_stems(dl)
        for be in BACKENDS:
            rs = refold_stems(str(pool / f"refold_{be}" / f"msaoff_{u}" / "designs.tar.zst.list"), TAG[be])
            miss = sorted(ds - rs)
            if miss:
                gap[be][(stack, u)] = miss

# --- live jobs ------------------------------------------------------------------
active_sb = set(l.strip() for l in sh("squeue","-u","lvalsecchi","-h","-o","%o").stdout.splitlines() if l.strip())
active_leaves = set()
for o in active_sb:
    if o.endswith(".sbatch"):
        try:
            m = re.search(r"SCRIPTS=\((.*?)\n\)", Path(o).read_text(), re.DOTALL)
        except OSError:
            m = None
        if m:
            active_leaves.update(x.strip() for x in m.group(1).splitlines() if x.strip())


def leaf_live(rdir):
    sbp = str(rdir / "run.sbatch")
    return sbp in active_sb or sbp in active_leaves


def clean_sbatch(sb, part, wall, nodelist, gres):
    sb = sb.replace("#SBATCH --dependency=afterok:0\n", "")
    sb = re.sub(r"#SBATCH --signal=\S+\n", "", sb)
    sb = re.sub(r"#SBATCH --nodelist=.*\n", "", sb)
    sb = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={part}", sb)
    sb = re.sub(r"#SBATCH --time=\S+", f"#SBATCH --time={wall}", sb)
    sb = re.sub(r"#SBATCH --gres=\S+", f"#SBATCH --gres={gres}", sb)
    if nodelist:
        sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", rf"\1#SBATCH --nodelist={nodelist}\n", sb, count=1)
    if "--requeue" not in sb:
        sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", r"\1#SBATCH --requeue\n", sb, count=1)
    return sb


def inject_allowlist(sb, allow_path):
    i = sb.rfind("\nbash ")
    return sb[:i] + f"\nexport REFOLD_ALLOWLIST={shlex.quote(str(allow_path))}" + sb[i:]


# --- build gap leaves -----------------------------------------------------------
buckets = collections.defaultdict(list)
skipped_live = 0
for backend in BACKENDS:
    lanes = BACKEND_LANES[backend]
    needing = {k: v for k, v in gap[backend].items() if v}
    items = sorted(needing.items(), key=lambda kv: len(kv[1]), reverse=(args.order == "largest"))
    if args.limit:
        items = items[: args.limit]
    print(f"[{backend}] uuids needing={len(needing)}  gap designs={sum(len(v) for v in needing.values())}"
          f"  building={len(items)} ({args.order} first)")
    rr = 0
    for (d, uuid), stems in items:
        HD = CAMP / d / RUN
        if not (HD / "designs" / f"{uuid}.tar.zst").exists():
            continue
        jobs = [j for j in launch._build_refold_jobs([HD], RUN, "0", None,
                                                     designs_filter=uuid, leaf_suffix=f"_{uuid}_gap")
                if j[0] == backend]
        if not jobs:
            continue
        _rn, rdir, rcfg, sb = jobs[0]
        if (rdir / "designs.tar.zst").exists():
            continue
        if leaf_live(rdir):
            skipped_live += 1
            continue
        li = rr % len(lanes); rr += 1
        sub_part, sub_qos, wall, nodelist, gres, _bp, _bq = lanes[li]
        sb = clean_sbatch(sb, sub_part, wall, nodelist, gres)
        allow_path = rdir / "allowlist.txt"
        sb = inject_allowlist(sb, allow_path)
        if args.submit:
            rdir.mkdir(parents=True, exist_ok=True)
            allow_path.write_text("".join(f"{s}\n" for s in stems))
            (rdir / "resolved_config.yaml").write_text(yaml.dump(rcfg))
            p = rdir / "run.sbatch"
            p.write_text(sb)
            p.chmod(0o755)
        buckets[(backend, li)].append((rdir, rdir / "run.sbatch"))

print(f"skipped(live)={skipped_live}")
print("plan (submit tier -> scontrol-move to bulk):")
for backend in BACKENDS:
    for li, lane in enumerate(BACKEND_LANES[backend]):
        sub_part, sub_qos, wall, nodelist, gres, bp, bq = lane
        n = len(buckets.get((backend, li), []))
        pin = f" -w {nodelist}" if nodelist else ""
        print(f"  {backend:<4}[{li}] {sub_part:<10}{gres:<11}{wall}{pin}  ->  {bp}: {n} task(s)")

if not args.submit:
    print("\nDRY-RUN — nothing submitted.")
    sys.exit(0)

TS = datetime.now().strftime("%Y%m%d_%H%M%S")
ARRROOT = CAMP / "_gapfill_arrays"
for (backend, li), leaves in buckets.items():
    if not leaves:
        continue
    sub_part, sub_qos, wall, nodelist, gres, bulk_part, bulk_qos = BACKEND_LANES[backend][li]
    adir = ARRROOT / f"{RUN}_{backend}_{re.sub(r'[^A-Za-z0-9]', '_', bulk_part)}_{TS}"
    adir.mkdir(parents=True, exist_ok=True)
    (adir / "logs").mkdir(exist_ok=True)
    (adir / "manifest.txt").write_text("".join(f"{p}\n" for _rd, p in leaves))
    hdr = ["#!/bin/bash",
           f"#SBATCH --job-name=rfgap_{backend}",
           f"#SBATCH --partition={sub_part}", f"#SBATCH --qos={sub_qos}",
           f"#SBATCH --time={wall}", "#SBATCH --cpus-per-task=4", "#SBATCH --mem-per-cpu=8G",
           "#SBATCH --gpus=1", f"#SBATCH --gres={gres}"]
    if nodelist:
        hdr.append(f"#SBATCH --nodelist={nodelist}")
    hdr += ["#SBATCH --requeue", "#SBATCH --open-mode=append",
            f"#SBATCH --array=0-{len(leaves)-1}", f"#SBATCH --chdir={RB}/code",
            f"#SBATCH --output={adir}/logs/%A_%a.log", f"#SBATCH --error={adir}/logs/%A_%a.err",
            "", "set -euo pipefail", "SCRIPTS=("]
    hdr += ["  %s" % shlex.quote(str(p)) for _rd, p in leaves]
    hdr += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
            'echo "[rfgap] task=$SLURM_ARRAY_TASK_ID host=$(hostname) $S"', 'bash "$S"']
    arrsh = adir / "array.sbatch"
    arrsh.write_text("\n".join(hdr) + "\n")
    aid = sh("sbatch", "--hold", "--parsable", str(arrsh)).stdout.strip().split(";")[0]
    if not aid.isdigit():
        print(f"  SUBMIT FAILED {backend} {sub_part}: {aid}")
        continue
    mv = sh("scontrol", "update", f"jobid={aid}", f"partition={bulk_part}", f"qos={bulk_qos}")
    rel = sh("scontrol", "release", f"jobid={aid}")
    ok = ("moved->" + bulk_part) if mv.returncode == 0 else f"MOVE FAILED ({mv.stderr.strip()})"
    if rel.returncode != 0:
        ok += f"; RELEASE FAILED ({rel.stderr.strip()})"
    print(f"  SUBMITTED {backend} array {aid} ({len(leaves)} tasks) {sub_part}->{bulk_part}; {ok}")
