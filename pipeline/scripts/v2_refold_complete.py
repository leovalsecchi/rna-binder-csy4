#!/usr/bin/env python3
# v2_refold_complete.py — complete missing AF3 and/or ESM refolds in midrun_v2
# BoltzGen/BoltzDesign, array-native, DESIGN-LEVEL (folds only the gaps).
#
# APPROACH (no data-loss merge): for each uuid with gaps, fold ONLY the missing
# design stems into a SEPARATE "_gap" leaf (refold_<be>/msaoff_<uuid>_gap) via a
# REFOLD_ALLOWLIST file (workers patched 2026-07-01 to honor it). The original leaf
# is left untouched; collect.py dedups AF3/ESM per design across leaves, so gaps +
# originals combine cleanly. Gap stems come from the collect gap-CSV (per-design
# has_af3_msaoff / has_esm).
#
# PLACEMENT — all on BULK via scontrol-move (Euler's submit router remaps *.bulk ->
# time-tier, so we submit to a tier then `scontrol update Partition=<bulk>`):
#   AF3  -> bulk, EXCLUDING 2080s (needs 24G):
#           gpu.bulk   pinned eu-g3-[062-077,083-084] (24G Turing, no 2080)  + gpu:1
#           gpuhe.bulk gpumem:24g (Ampere/4090/eu-g3; 24g req excludes 11G 2080s)
#   ESM  -> bulk, ALL nodes incl 2080 (gpumem:8g fits 11G): gpu.bulk + gpuhe.bulk, no pin
#
#   module load stack/2024-06 python/3.12.8
#   python3 v2_refold_complete.py                       # DRY-RUN both backends
#   python3 v2_refold_complete.py --order smallest --limit 2 --submit   # CANARY
#   python3 v2_refold_complete.py --submit              # full
import sys, re, argparse, shlex, subprocess, csv, collections
from datetime import datetime
from pathlib import Path

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/scripts")
import yaml
import launch

CAMP = Path(RB) / "campaigns" / "csy4_crrna_free"
RUN = "midrun_v2"
GAP_CSV = f"{RB}/analyses/cross_v2_bg_bd/collected_raw.csv"
STACK_DISPLAY_TO_DIR = {
    "BG+BIF": "boltzgen__boltzif",
    "BG+LiMPNN": "boltzgen__ligandmpnn",
    "BD1+LiMPNN": "boltzdesign__ligandmpnn",
}
POOL_G3 = "eu-g3-[062-077],eu-g3-[083-084]"    # the 24G Turing nodes (no 2080)
# lane = (submit_part, submit_qos, wall, nodelist, gres, bulk_part, bulk_qos)
BACKEND_LANES = {
    "af3": [
        ("gpu.120h",   "es_wicky/gpu/120",   "48:00:00", POOL_G3, "gpu:1",      "gpu.bulk",   "es_wicky/gpu/bulk"),
        ("gpuhe.120h", "es_wicky/gpuhe/120",  "48:00:00", None,    "gpumem:24g", "gpuhe.bulk", "es_wicky/gpuhe/bulk"),
    ],
    # ESM needs ~10.5G -> OOMs on 11G 2080s (fine on idle owner 2080s, dies when
    # contended; verified 2026-07-01: 133/133 fail on shared public 2080s). So ESM
    # goes on 24G like AF3 (eu-g3 pin / gpumem:24g excludes 2080). ESM is ~40min/leaf
    # so 24G throughput is fine.
    "esm": [
        ("gpu.4h",   "es_wicky/gpu/4",   "04:00:00", POOL_G3, "gpu:1",      "gpu.bulk",   "es_wicky/gpu/bulk"),
        ("gpuhe.4h", "es_wicky/gpuhe/4",  "04:00:00", None,    "gpumem:24g", "gpuhe.bulk", "es_wicky/gpuhe/bulk"),
    ],
}

ap = argparse.ArgumentParser()
ap.add_argument("--backends", default="af3,esm")
ap.add_argument("--submit", action="store_true")
ap.add_argument("--limit", type=int, default=0, help="cap uuids PER BACKEND")
ap.add_argument("--order", choices=["largest", "smallest"], default="largest")
ap.add_argument("--gap-csv", default=GAP_CSV)
args = ap.parse_args()
BACKENDS = [b.strip() for b in args.backends.split(",") if b.strip()]


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True)


def tf(v):
    return str(v).strip().lower() in ("1", "true", "yes")


# --- gap STEMS per (dir,uuid) per backend ---------------------------------------
gap = {"af3": collections.defaultdict(list), "esm": collections.defaultdict(list)}
with open(args.gap_csv) as fh:
    for r in csv.DictReader(fh):
        d = STACK_DISPLAY_TO_DIR.get(r["stack"])
        if not d:
            continue
        key = (d, r["uuid8"])
        did = r["design_id"]
        if not tf(r.get("has_af3_msaoff")):
            gap["af3"][key].append(did)
        if not tf(r.get("has_esm")):
            gap["esm"][key].append(did)

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
    """Prepend `export REFOLD_ALLOWLIST=...` before the worker `bash` invocation."""
    i = sb.rfind("\nbash ")
    exp = f"\nexport REFOLD_ALLOWLIST={shlex.quote(str(allow_path))}"
    return sb[:i] + exp + sb[i:]


# --- build gap leaves, bucket by (backend, lane) --------------------------------
buckets = collections.defaultdict(list)   # (backend, li) -> [(rdir, run.sbatch)]
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
        if (rdir / "designs.tar.zst").exists():   # this gap leaf already completed
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

# --- plan ------------------------------------------------------------------------
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

# --- submit one array per (backend, lane), then move to bulk --------------------
TS = datetime.now().strftime("%Y%m%d_%H%M%S")
ARRROOT = CAMP / "_v2_refold_complete_arrays"
for (backend, li), leaves in buckets.items():
    if not leaves:
        continue
    sub_part, sub_qos, wall, nodelist, gres, bulk_part, bulk_qos = BACKEND_LANES[backend][li]
    adir = ARRROOT / f"{backend}_{re.sub(r'[^A-Za-z0-9]', '_', bulk_part)}_{TS}"
    adir.mkdir(parents=True, exist_ok=True)
    (adir / "logs").mkdir(exist_ok=True)
    (adir / "manifest.txt").write_text("".join(f"{p}\n" for _rd, p in leaves))
    hdr = ["#!/bin/bash",
           f"#SBATCH --job-name=rf{backend}_v2complete",
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
            'echo "[v2rf] task=$SLURM_ARRAY_TASK_ID host=$(hostname) $S"', 'bash "$S"']
    arrsh = adir / "array.sbatch"
    arrsh.write_text("\n".join(hdr) + "\n")
    # Submit HELD so the job can't start before the partition move applies (else a
    # free GPU starts it instantly on the submit tier and the move fails). Move to
    # bulk while held, then release -> it schedules on bulk.
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
