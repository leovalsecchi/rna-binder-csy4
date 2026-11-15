#!/usr/bin/env python3
# v2_af3_complete.py — complete the missing AF3 refolds in the midrun_v2
# BoltzGen/BoltzDesign stacks, array-native, across FOUR AF3 lanes.
#
# WHY: the cross-campaign consensus (AF3 & B2) can only rank designs that have
# BOTH refolds. B2 is ~100% present in v2; AF3 is only ~75% (BG+LiMPNN 51%). This
# fills the AF3 gap so BG/BD become consensus-rankable.
#
# MECHANISM (per-uuid, resume-native):
#   * refold_af3.sh ALREADY resumes: it restores the leaf's packed outputs.tar.zst,
#     skips designs whose model.cif exists, and re-batches only the remainder. So
#     re-running a PARTIAL uuid's AF3 leaf folds only its missing designs — no
#     design-level filter needed for full completion.
#   * The one trap: a leaf with designs.tar.zst present makes refold_af3.sh exit 0
#     ("already complete"), even when the pack is partial (falsely-done leaves —
#     ~49 in BG+LiMPNN). needs_work() re-derives completeness from the collect
#     gap-set, and on --submit we REOPEN a falsely-done leaf (mv designs.tar.zst ->
#     .bak_reopen; refold_af3.sh then resumes from outputs.tar.zst).
#
# PLACEMENT — 4 AF3 lanes (node mapping verified 2026-07-01):
#   public-A  gpu.24h   es_wicky/gpu/24   eu-g3-[062-070]          gres=gpu:1  (gpumem dropped)
#   public-B  gpu.120h  es_wicky/gpu/120  eu-g3-[071-077,083-084]  gres=gpu:1
#   bulk-gpu  gpu.bulk  es_wicky/gpu/bulk eu-g3-[062-077,083-084]  gres=gpu:1  (MUST pin: bulk also has 11G eu-lo)
#   bulk-ghe  gpuhe.bulk es_wicky/gpuhe/bulk  (no pin)             gres=gpumem:24g (routes eu-g4/g6/g3 >=24G)
# Big leaves (many missing designs -> long) are steered to the long-wall lanes.
#
#   module load stack/2024-06 python/3.12.8
#   python3 v2_af3_complete.py                 # DRY-RUN: print the plan, submit nothing
#   python3 v2_af3_complete.py --submit        # build leaves + fire one array per lane
#   python3 v2_af3_complete.py --limit 5       # cap uuids (smoke test)
import sys, re, argparse, shlex, subprocess, csv, collections
from datetime import datetime
from pathlib import Path

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/scripts")
import yaml
import launch

CAMP = Path(RB) / "campaigns" / "csy4_crrna_free"
STACKS = ["boltzgen__boltzif", "boltzgen__ligandmpnn", "boltzdesign__ligandmpnn"]
RUN = "midrun_v2"
GAP_CSV = f"{RB}/analyses/cross_v2_bg_bd/collected_raw.csv"  # authoritative has_af3 per design

# eu-g3 = titan_rtx + quadro_rtx_6000, 24G. Split public tiers by node bin so the two
# QOS (gpu/24 cap50, gpu/120 cap29) don't contend. gpu.bulk also carries 11G eu-lo -> pin.
POOL_A = "eu-g3-[062-070]"
POOL_B = "eu-g3-[071-077],eu-g3-[083-084]"
POOL_ALL = "eu-g3-[062-077],eu-g3-[083-084]"
# (partition, qos, walltime, nodelist_or_None, keep_gpumem_bool)
LANES = [
    ("gpu.24h",   "es_wicky/gpu/24",    "24:00:00", POOL_A,   False),  # public-A
    ("gpu.120h",  "es_wicky/gpu/120",   "48:00:00", POOL_B,   False),  # public-B (long)
    ("gpu.bulk",  "es_wicky/gpu/bulk",  "24:00:00", POOL_ALL, False),  # bulk-gpu
    ("gpuhe.bulk","es_wicky/gpuhe/bulk","24:00:00", None,     True),   # bulk-gpuhe (long-capable)
]
LONG_LANES = [1, 3]   # indices allowed to take big leaves (120h + gpuhe.bulk)
BIG_MISSING = 700     # a leaf with more missing designs than this -> long lane only

ap = argparse.ArgumentParser()
ap.add_argument("--submit", action="store_true")
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--gap-csv", default=GAP_CSV)
args = ap.parse_args()


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True)


def tf(v):
    return str(v).strip().lower() in ("1", "true", "yes")


# --- 1. gap set: per (stack_dirname, uuid) the count of designs missing AF3 -----
# collect's stack labels are display form (BG+BIF); map to campaign dir names.
STACK_DISPLAY_TO_DIR = {
    "BG+BIF": "boltzgen__boltzif",
    "BG+LiMPNN": "boltzgen__ligandmpnn",
    "BD1+LiMPNN": "boltzdesign__ligandmpnn",
}
missing = collections.Counter()   # (dir, uuid) -> n missing af3
totals = collections.Counter()    # (dir, uuid) -> n designs
with open(args.gap_csv) as fh:
    for r in csv.DictReader(fh):
        d = STACK_DISPLAY_TO_DIR.get(r["stack"])
        if not d:
            continue
        key = (d, r["uuid8"])
        totals[key] += 1
        if not tf(r.get("has_af3_msaoff")):
            missing[key] += 1
needing = {k: m for k, m in missing.items() if m > 0}
print(f"[v2_af3_complete] gap-csv={args.gap_csv}")
print(f"  uuids needing AF3: {len(needing)}   missing designs: {sum(needing.values())}")

# --- 2. live jobs -> don't double-submit ---------------------------------------
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


# --- 3. build per-uuid AF3 leaves, route to lanes, bucket ----------------------
todo = sorted(needing.items(), key=lambda kv: -kv[1])   # biggest gaps first
if args.limit:
    todo = todo[: args.limit]

buckets = collections.defaultdict(list)   # lane_idx -> [(rdir, run.sbatch)]
rr = 0
reopened = skipped_live = 0
for (d, uuid), nmiss in todo:
    HD = CAMP / d / RUN
    if not (HD / "designs" / f"{uuid}.tar.zst").exists():
        continue
    # AF3-only build via the project builder (same call the HD loop uses)
    jobs = [j for j in launch._build_refold_jobs([HD], RUN, "0", None,
                                                 designs_filter=uuid, leaf_suffix=f"_{uuid}")
            if j[0] == "af3"]
    if not jobs:
        continue
    _rn, rdir, rcfg, sb = jobs[0]
    if leaf_live(rdir):
        skipped_live += 1
        continue
    # pick lane: big gaps -> long lanes only; else round-robin over all 4
    if nmiss > BIG_MISSING:
        li = LONG_LANES[rr % len(LONG_LANES)]
    else:
        li = rr % len(LANES)
    rr += 1
    part, qos, wall, nodelist, keep_gpumem = LANES[li]

    # clean the sbatch string: drop dep + USR1 signal + stale eu-g4 nodelist; set
    # partition/walltime; drop gpumem -> gpu:1 unless this lane keeps it; pin nodes.
    sb = sb.replace("#SBATCH --dependency=afterok:0\n", "")
    sb = re.sub(r"#SBATCH --signal=\S+\n", "", sb)
    sb = re.sub(r"#SBATCH --nodelist=.*\n", "", sb)
    sb = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={part}", sb)
    sb = re.sub(r"#SBATCH --time=\S+", f"#SBATCH --time={wall}", sb)
    if not keep_gpumem:
        sb = re.sub(r"#SBATCH --gres=gpumem:\S+", "#SBATCH --gres=gpu:1", sb)
    if nodelist:
        sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", rf"\1#SBATCH --nodelist={nodelist}\n", sb, count=1)
    if "--requeue" not in sb:
        sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", r"\1#SBATCH --requeue\n", sb, count=1)

    if not args.submit:
        buckets[li].append((rdir, rdir / "run.sbatch"))
        continue

    # REOPEN a falsely-done leaf so refold_af3.sh doesn't exit-0. outputs.tar.zst
    # (kept) lets its resume skip the already-folded designs; canonical pack is
    # regenerated at the end. mv (not rm) the done-marker so nothing is lost.
    done = rdir / "designs.tar.zst"
    if done.exists():
        done.rename(rdir / "designs.tar.zst.bak_reopen")
        reopened += 1
    rdir.mkdir(parents=True, exist_ok=True)
    (rdir / "resolved_config.yaml").write_text(yaml.dump(rcfg))
    p = rdir / "run.sbatch"
    p.write_text(sb)
    p.chmod(0o755)
    buckets[li].append((rdir, p))

# --- 4. plan / submit one array per lane ---------------------------------------
print(f"  skipped (live): {skipped_live}   reopened falsely-done leaves: {reopened}")
print("  lane plan:")
for li, (part, qos, wall, nodelist, keep) in enumerate(LANES):
    n = len(buckets.get(li, []))
    pin = f" -w {nodelist}" if nodelist else ""
    gm = "gpumem:24g" if keep else "gpu:1"
    print(f"    [{li}] {part:<11} {qos:<20} {wall} {gm:<11}{pin}: {n} task(s)")

if not args.submit:
    print("\nDRY-RUN — nothing submitted. Re-run with --submit to fire the arrays.")
    sys.exit(0)

TS = datetime.now().strftime("%Y%m%d_%H%M%S")
ARRROOT = CAMP / "_v2_af3_complete_arrays"
for li, leaves in buckets.items():
    if not leaves:
        continue
    part, qos, wall, nodelist, keep = LANES[li]
    adir = ARRROOT / f"af3_{re.sub(r'[^A-Za-z0-9]', '_', part)}_{TS}"
    adir.mkdir(parents=True, exist_ok=True)
    (adir / "logs").mkdir(exist_ok=True)
    (adir / "manifest.txt").write_text("".join(f"{p}\n" for _rd, p in leaves))
    gdir = "#SBATCH --gpus=1\n#SBATCH --gres=gpumem:24g" if keep else "#SBATCH --gpus=1\n#SBATCH --gres=gpu:1"
    ndir = f"#SBATCH --nodelist={nodelist}\n" if nodelist else ""
    arrsh = adir / "array.sbatch"
    lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=rfaf3_v2complete",
        f"#SBATCH --partition={part}",
        f"#SBATCH --qos={qos}",
        f"#SBATCH --time={wall}",
        "#SBATCH --cpus-per-task=4",
        "#SBATCH --mem-per-cpu=8G",
        gdir,
        ndir.rstrip("\n") if ndir else "",
        "#SBATCH --requeue",
        "#SBATCH --open-mode=append",
        f"#SBATCH --array=0-{len(leaves)-1}",
        f"#SBATCH --chdir={RB}/code",
        f"#SBATCH --output={adir}/logs/%A_%a.log",
        f"#SBATCH --error={adir}/logs/%A_%a.err",
        "", "set -euo pipefail", "SCRIPTS=(",
    ]
    lines += ["  %s" % shlex.quote(str(p)) for _rd, p in leaves]
    lines += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
              'echo "[v2af3] task=$SLURM_ARRAY_TASK_ID host=$(hostname) $S"', 'bash "$S"']
    arrsh.write_text("\n".join(x for x in lines if x != "") + "\n")
    aid = sh("sbatch", str(arrsh)).stdout
    print(f"  SUBMITTED {part} array over {len(leaves)} tasks: {aid.strip()}")
