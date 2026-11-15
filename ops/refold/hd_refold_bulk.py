#!/usr/bin/env python3
# hd_refold_bulk.py — submit refolds for FINISHED HalluDesign generation uuids,
# placing them on the BULK queues of gpu + gpuhe.
#
#   * "Finished" = a packed design slice  <run>/designs/<uuid8>.tar.zst  exists AND
#     the uuid is NOT currently generating (no live  hd_<uuid8>  gen job).
#   * Per finished uuid it builds the 3 refold leaves via launch._build_refold_jobs
#     (AF3 + B2 + ESM, msaoff), submits with NO afterok dependency (gen is done),
#     then places each on the right tier:
#         AF3 -> PUBLIC gpu.24h/gpu.120h, eu-g3 Turing 24G pool, gpumem DROPPED
#                (-> --gres=gpu:1) + --nodelist pin. The gpumem:24g constraint is what
#                forces the router onto the saturated owner gpuhe.*; dropping it +
#                pinning the (physically-24G) eu-g3 nodes lands+RUNS on the public pool
#                immediately. Binned across two tiers/node-groups so the gpu/24 (cap 50
#                gpu) and gpu/120 (cap 29) QOS don't contend. AF3 runs on Turing via xla.
#                (See claude_memory rna-binder-public-gpu-pool-routing; rescue 2026-06-30.)
#         ESM -> gpuhe.bulk   (gpumem:24g, unchanged)
#         B2  -> gpu.bulk     (gpumem:8g, public Turing pool, unchanged)
#   * Strips the stale `--nodelist` the refold configs bake in (eu-g4-[..] is invalid),
#     then re-pins it to the eu-g3 pool for the public-pool (af3) path. Adds
#     `--requeue` (these tiers are PreemptMode=REQUEUE).
#   * Orphan-aware idempotency: a backend leaf is "handled" only if it has
#     designs.tar.zst (done) OR a live job (its run.sbatch is in squeue). A leftover
#     leaf with neither (e.g. a failed submit) is rebuilt and resubmitted.
#
# Run on Euler. Needs the pipeline python:  module load stack/2024-06 python/3.12.8
#   python3 hd_refold_bulk.py [--run midrun_v3] [--submit] [--limit N]
import sys, subprocess, re, argparse, shutil, time
from datetime import datetime
from pathlib import Path

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/scripts")
import yaml
import launch

CAMP = Path(RB) / "campaigns" / "csy4_crrna_free"
STACK = CAMP / "halludesign__ligandmpnn"

# PUBLIC gpu.* Turing 24G pool (eu-g3 = titan_rtx + quadro_rtx_6000), split into two
# node bins so the two QOS tiers (gpu/24 cap=50 gpu, gpu/120 cap=29) don't contend for
# the same nodes. Everything else in gpu.24h/gpu.120h is 11G 2080 Ti -> too small for AF3.
POOL_A = "eu-g3-[062-070]"                    # gpu.24h  bin (9 nodes)
POOL_B = "eu-g3-[071-077],eu-g3-[083-084]"    # gpu.120h bin (9 nodes)

# Per-backend placement(s); submissions round-robin across the list.
# Tuple = (partition, qos, walltime, nodelist).  nodelist set => "public-pool path":
#   drop gpumem (-> --gres=gpu:1) + pin nodes + bake partition/walltime so the router
#   keeps it on the public gpu.* tier (walltime picks .24h/.120h; no gpumem => no gpuhe
#   force).  nodelist None => legacy path (keep gpumem, scontrol-move to the tier).
PLACE = {
    "af3": [("gpu.bulk", "es_wicky/gpu/bulk", "04:00:00", "eu-g3-[062-077,083-084]")],
    "esm": [("gpu.bulk", "es_wicky/gpu/bulk", "04:00:00", "eu-g3-[062-077,083-084]")],
    "b2":  [("gpu.bulk",   "es_wicky/gpu/bulk",   None, None)],
}
rr = {}  # per-backend round-robin counter

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="")
ap.add_argument("--submit", action="store_true")
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--gpuhe", action="store_true",
                help="route ALL backends onto the freed gpuhe 24G pool (gpuhe.bulk, "
                     "es_wicky/gpuhe/bulk, gpumem:24g, no node pin) instead of public gpu.*")
args = ap.parse_args()

# --gpuhe: send everything to the (currently freed) gpuhe bulk tier. nodelist=None keeps
# us off the public-pool rewrite; the elif-gpuhe sbatch branch bakes partition+gpumem:24g
# and the scontrol-move below confirms qos. gpuhe.bulk is preemptible => yields to gen.
if args.gpuhe:
    PLACE = {
        # af3/esm need gpuhe 24G; ~80% bulk + ~20% timed via round-robin so the submitter
        # emits ONE clean array per tier (NO post-hoc per-task scontrol -> no fragmentation).
        # BALANCED across BOTH pools — the QOS caps are PER POOL (gpu/* 96/50/29/240 AND
        # gpuhe/* 96/50/29/240), so spanning both reaches ~830 GPU-slots instead of ~415.
        # ~50/50 gpu.* vs gpuhe.*, ~75% bulk / ~25% timed inside each pool.
        # af3 has NO 4h lane in EITHER pool: p99 is 229 min (3.8 h) over 11,280 COMPLETED
        # rfaf3, so a 4 h cap TIMEOUTs the tail and loses the whole chunk. ESM keeps 4h
        # (p90 55 min). gpu.* lanes keep the eu-g3 pin — the dispatch swaps gpumem->gpu:1
        # for pinned lanes, so unpinned they would land on 11 G 2080 Tis.
        "af3": [("gpu.bulk",    "es_wicky/gpu/bulk",   "06:00:00",   "eu-g3-[062-077,083-084]")] * 6
               + [("gpu.24h",    "es_wicky/gpu/24",    "23:30:00",   "eu-g3-[062-077,083-084]"),
                  ("gpu.120h",   "es_wicky/gpu/120",   "47:30:00",   "eu-g3-[062-077,083-084]")]
               + [("gpuhe.bulk", "es_wicky/gpuhe/bulk", "06:00:00",   None)] * 6
               + [("gpuhe.24h",  "es_wicky/gpuhe/24",  "06:00:00",   None),
                  ("gpuhe.120h", "es_wicky/gpuhe/120", "1-06:00:00", None)],
        "esm": [("gpu.bulk",    "es_wicky/gpu/bulk",   "06:00:00",   "eu-g3-[062-077,083-084]")] * 5
               + [("gpu.4h",     "es_wicky/gpu/4",     "04:00:00",   "eu-g3-[062-077,083-084]"),
                  ("gpu.24h",    "es_wicky/gpu/24",    "23:30:00",   "eu-g3-[062-077,083-084]"),
                  ("gpu.120h",   "es_wicky/gpu/120",   "47:30:00",   "eu-g3-[062-077,083-084]")]
               + [("gpuhe.bulk", "es_wicky/gpuhe/bulk", "06:00:00",   None)] * 5
               + [("gpuhe.4h",   "es_wicky/gpuhe/4",   "04:00:00",   None),
                  ("gpuhe.24h",  "es_wicky/gpuhe/24",  "06:00:00",   None),
                  ("gpuhe.120h", "es_wicky/gpuhe/120", "1-06:00:00", None)],
        # b2 fits the 11G 2080 Ti (gpumem:8g) -> gpu.bulk on the PUBLIC 2080 pool, OFF the
        # precious gpuhe 24G pool that af3/esm need. gpumem:8g auto-routes to public Turing.
        # b2 fits 11 G (gpumem:8g) -> PUBLIC 2080 Ti pool, no pin, and stays OFF the 24 G
        # cards af3/esm/apo need. Spread 4h/24h/bulk per distribute_refolds.sh (b2 IS fine
        # on a 4h tier: median 102 min, p90 109 min).
        "b2": [("gpu.bulk", "es_wicky/gpu/bulk", "06:00:00", None)] * 2
              + [("gpu.4h",  "es_wicky/gpu/4",   "04:00:00", None),
                 ("gpu.24h", "es_wicky/gpu/24",  "23:30:00", None)],
    }


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True)


# resolve run
if args.run:
    RUN = args.run
else:
    cands = []
    for p in STACK.glob("midrun_v*"):
        m = re.search(r"v(\d+)", p.name)
        if m and (p / "designs").is_dir():
            cands.append((int(m.group(1)), p.name))
    RUN = sorted(cands)[-1][1] if cands else "midrun_v3"
HD = STACK / RUN

print(f"[hd_refold_bulk] run={RUN}  HD={HD}  exists={HD.exists()}  submit={args.submit}")
if not (HD / "designs").is_dir():
    print("  no designs/ dir yet -> nothing to do this pass")
    sys.exit(0)

# live HD gen uuids (not finished -> skip) + active refold run.sbatch paths (job exists)
inflight = set()
# --- consolidation: cancel our OWN still-pending arrays from the previous tick so their
# designs re-enter `todo` and get re-bucketed into a few big arrays. Without this, each
# tick leaves a fresh micro-array pending -> hundreds accumulate and flood the 250-deep
# scheduler window (low-prio bulk refolds then starve while GPUs sit idle). Scoped by a
# per-(submitter,run) job-id state file so we never cancel the sibling loop's arrays even
# when they share this pool (HD) or the rf<name>_<run> job-name.
_CONS_STATE = HD / f".refold_consolidate.hd.{RUN}.jobids"
if args.submit and _CONS_STATE.exists():
    _old = [x for x in _CONS_STATE.read_text().split() if x.isdigit()]
    if _old:
        sh("scancel", "--state=PENDING", *_old)
        time.sleep(5)  # let the controller reflect the cancel before reading squeue

active_sb = set()
for line in sh("squeue", "-u", "lvalsecchi", "-h", "-o", "%j").stdout.splitlines():
    m = re.match(r"hd_([0-9a-f]{8})$", line.strip())
    if m:
        inflight.add(m.group(1))
for line in sh("squeue", "-u", "lvalsecchi", "-h", "-o", "%o").stdout.splitlines():
    s = line.strip()
    if s:
        active_sb.add(s)

# Array-aware idempotency: a leaf submitted inside a live job array shows up in
# squeue as the ARRAY's *.sbatch (not the leaf's run.sbatch). So read each live
# array's batch script and collect the leaf run.sbatch paths it covers — from an
# inline SCRIPTS=() block (self-contained arrays, INCLUDING the parallel rec_*
# recovery arrays) or, failing that, a sibling manifest.txt (legacy). Without this
# the next pass resubmits those leaves as duplicates AND rmtree's a leaf mid-refold.
active_leaves = set()
for o in active_sb:
    if not o.endswith(".sbatch"):
        continue
    try:
        txt = Path(o).read_text()
    except OSError:
        continue
    m = re.search(r"SCRIPTS=\((.*?)\n\)", txt, re.DOTALL)
    if m:
        for ln in m.group(1).splitlines():
            if ln.strip():
                active_leaves.add(ln.strip())
    elif o.endswith("/array.sbatch"):
        man = Path(o).parent / "manifest.txt"
        try:
            for ln in man.read_text().splitlines():
                if ln.strip():
                    active_leaves.add(ln.strip())
        except OSError:
            pass

finished = []
for f in sorted((HD / "designs").glob("*.tar.zst")):
    u = f.name.split(".")[0]
    if re.fullmatch(r"[0-9a-f]{8}", u) and u not in inflight:
        finished.append(u)


def leaf_handled(rdir):
    """Done (packed) or has a live job — its run.sbatch is queued/running as an
    individual job (active_sb) or as a task inside a live array (active_leaves)."""
    sbp = str(rdir / "run.sbatch")
    return ((rdir / "designs.tar.zst").exists()
            or sbp in active_sb or sbp in active_leaves)


def needs_work(u):
    return not all(leaf_handled(HD / f"refold_{k}" / f"msaoff_{u}") for k in ("af3", "b2", "esm"))


todo = [u for u in finished if needs_work(u)]
if args.limit:
    todo = todo[: args.limit]
print(f"  finished uuids={len(finished)}  in-flight-skipped={len(inflight)}  needing-work={len(todo)}")

nsub = nmoved = 0
log = []
# Bucket the leaves to submit by placement, then fire ONE job array per bucket
# (queue-depth fix) instead of one sbatch per leaf. Each leaf's run.sbatch is
# written now; the array task runs `bash <leaf>/run.sbatch`. The array header
# reuses the cleaned per-leaf #SBATCH directives (identical within a bucket).
buckets = {}  # (part,qos,wall,nodelist) -> {"rname", "sb", "leaves": [(rdir, run.sbatch)]}
for u in todo:
    for rname, rdir, rcfg, sb in launch._build_refold_jobs(
        [HD], RUN, "0", None, designs_filter=u, leaf_suffix=f"_{u}"
    ):
        if leaf_handled(rdir):
            log.append(f"{u} {rname} OK(done/active)")
            continue
        # pick placement (round-robin within backend)
        places = PLACE[rname]
        idx = rr.get(rname, 0); rr[rname] = idx + 1
        part, qos, wall, nodelist = places[idx % len(places)]

        # clean the strings: no dep, no signal, drop stale nodelist, requeue
        sb = sb.replace("#SBATCH --dependency=afterok:0\n", "")
        sb = re.sub(r"#SBATCH --signal=\S+\n", "", sb)           # USR1 grace = unschedulable on preemptible tiers (2026-06-28)
        sb = re.sub(r"#SBATCH --nodelist=.*\n", "", sb)          # drop the stale eu-g4-[..] the config bakes in
        if nodelist:
            # PUBLIC-POOL path: bake partition, drop gpumem (-> gpu:1), pin eu-g3 nodes.
            sb = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={part}", sb)
            sb = re.sub(r"#SBATCH --gres=gpumem:\S+", "#SBATCH --gres=gpu:1", sb)
            sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", rf"\1#SBATCH --nodelist={nodelist}\n", sb, count=1)
        elif part.startswith("gpuhe"):
            # GPUHE path: bake partition + qos + force gpumem:24g (matches every gpuhe
            # 24G node: 4090/3090/quadro), NO node pin. scontrol-move below confirms qos.
            sb = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={part}", sb)
            if re.search(r"#SBATCH --qos=", sb):
                sb = re.sub(r"#SBATCH --qos=\S+", f"#SBATCH --qos={qos}", sb)
            else:
                sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", rf"\1#SBATCH --qos={qos}\n", sb, count=1)
            if re.search(r"#SBATCH --gres=gpumem:\S+", sb):
                sb = re.sub(r"#SBATCH --gres=gpumem:\S+", "#SBATCH --gres=gpumem:24g", sb)
            else:
                sb = re.sub(r"#SBATCH --gres=\S+", "#SBATCH --gres=gpumem:24g", sb)
        elif part.startswith("gpu."):
            # PUBLIC gpu.* pool (2080 Ti): bake partition + qos, KEEP gpumem as-is
            # (b2 stays gpumem:8g -> routes to the 11G Turing public pool, off gpuhe).
            sb = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={part}", sb)
            if re.search(r"#SBATCH --qos=", sb):
                sb = re.sub(r"#SBATCH --qos=\S+", f"#SBATCH --qos={qos}", sb)
            else:
                sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", rf"\1#SBATCH --qos={qos}\n", sb, count=1)
        if wall:
            sb = re.sub(r"#SBATCH --time=\S+", f"#SBATCH --time={wall}", sb)
        if "--requeue" not in sb:
            sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", r"\1#SBATCH --requeue\n", sb, count=1)

        # node-local scratch worker + node-local disk + black-hole exclude
        sb = re.sub(r"(refold_(?:af3|b2|esm))\.sh", r"\1_scratch.sh", sb)
        if "--tmp=" not in sb:
            sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", r"\1#SBATCH --tmp=30G\n", sb, count=1)
        if "--exclude=" not in sb:
            sb = re.sub(r"(#SBATCH --partition=[^\n]*\n)", r"\1#SBATCH --exclude=eu-lo-g2-028,eu-lo-g3-059,eu-g6-064\n", sb, count=1)

        b = buckets.setdefault((part, qos, wall, nodelist),
                               {"rname": rname, "sb": sb, "leaves": []})
        if not args.submit:
            log.append(f"{u} {rname} would-submit -> {part}" + (f"  -w {nodelist}" if nodelist else ""))
            b["leaves"].append((rdir, rdir / "run.sbatch"))
            continue
        # write per-leaf run.sbatch (rebuild a stale/failed leaf)
        if rdir.exists():
            shutil.rmtree(rdir, ignore_errors=True)
        rdir.mkdir(parents=True, exist_ok=True)
        (rdir / "resolved_config.yaml").write_text(yaml.dump(rcfg))
        p = rdir / "run.sbatch"
        p.write_text(sb)
        p.chmod(0o755)
        b["leaves"].append((rdir, p))

# Submit one job array per bucket.
TS = datetime.now().strftime("%Y%m%d_%H%M%S")
ARRROOT = HD / "_refold_arrays"
# split each routing bucket into <=CHUNK-task arrays so no single array exceeds
# Slurm MaxArraySize (15000); still a handful of BIG arrays, never micro-arrays.
CHUNK = 10000
_jobs = []
for _sig, _b in buckets.items():
    _lv = _b["leaves"]
    for _ci in range(0, len(_lv), CHUNK):
        _jobs.append((_sig, _b, _lv[_ci:_ci + CHUNK], _ci // CHUNK))
_sub_jids = []
for (part, qos, wall, nodelist), b, leaves, _cix in _jobs:
    rname = b["rname"]
    n = len(leaves)
    if not args.submit:
        log.append(f"  [array] {rname} -> {part}{(' -w ' + nodelist) if nodelist else ''}: {n} task(s)")
        continue
    adir = ARRROOT / f"{rname}_{re.sub(r'[^A-Za-z0-9]', '_', part)}_{TS}_{_cix}"
    adir.mkdir(parents=True, exist_ok=True)
    (adir / "array_logs").mkdir(exist_ok=True)
    # manifest.txt is kept for provenance AND for this script's own array-aware
    # idempotency (active_leaves reads live arrays' manifest.txt). The array itself
    # runs from an INLINE SCRIPTS=() list, so deleting manifest.txt can never orphan
    # a live array (Slurm spools the batch script at submit). 2026-06-30.
    manifest = adir / "manifest.txt"
    manifest.write_text("".join(f"{p}\n" for _rd, p in leaves))
    # array header = the cleaned per-leaf #SBATCH directives (identical within a
    # bucket), minus job-name/output/error, plus --array + array-level logs.
    header_part = b["sb"].split("\nset -euo pipefail", 1)[0]
    hdr = [ln for ln in header_part.splitlines()
           if ln.startswith("#SBATCH ")
           and not ln.startswith(("#SBATCH --job-name",
                                  "#SBATCH --output", "#SBATCH --error"))]
    scripts = "\n".join(f"  {p}" for _rd, p in leaves)
    arr = (
        "#!/bin/bash\n"
        f"#SBATCH --job-name=rf{rname}_{RUN}\n"
        + "\n".join(hdr) + "\n"
        f"#SBATCH --array=0-{n - 1}%200\n"
        f"#SBATCH --output={adir}/array_logs/%A_%a.log\n"
        f"#SBATCH --error={adir}/array_logs/%A_%a.err\n\n"
        "set -euo pipefail\n"
        f"SCRIPTS=(\n{scripts}\n)\n"
        'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"\n'
        'RUNDIR=$(dirname "$S")\n'
        'echo "[refold-array] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"\n'
        'bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"\n'
    )
    ap_ = adir / "array.sbatch"
    ap_.write_text(arr)
    ap_.chmod(0o755)
    res = sh("sbatch", "--parsable", str(ap_))
    jid = res.stdout.strip().split(";")[0]
    if not jid.isdigit():
        log.append(f"[array] {rname} {part} SUBMIT_ERR {res.stderr.strip()[:120]}")
        continue
    nsub += n
    _sub_jids.append(jid)
    # confirm partition/qos on the whole array: legacy bulk path NEEDS this move;
    # public-pool path already lands correct (no-op confirm + stray-reroute repair).
    mv = sh("scontrol", "update", f"jobid={jid}", f"partition={part}", f"qos={qos}")
    if mv.returncode == 0:
        nmoved += n
        log.append(f"[array] {rname} jid={jid} ({n} tasks) -> {part}")
    else:
        log.append(f"[array] {rname} jid={jid} ({n} tasks) MOVE_FAIL {mv.stderr.strip()[:100]}")

if args.submit:
    _CONS_STATE.write_text("\n".join(_sub_jids) + ("\n" if _sub_jids else ""))
print(f"  submitted={nsub}  moved_to_bulk={nmoved}  arrays={len(buckets)}")
if args.submit:
    with open(HD / "_hd_refold_bulk.log", "a") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"  log -> {HD/'_hd_refold_bulk.log'}")
else:
    for line in log[:20]:
        print("   ", line)
    print("  DRY-RUN — add --submit to fire.")
