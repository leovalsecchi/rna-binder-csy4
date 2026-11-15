"""ck.cli — python entry behind `campaign probe` (history-harvest).

  campaign probe --model M [--days N] [--jobname REGEX] [--write]

Loads models/M.model.yaml + cluster/euler.yaml (site-resolved), harvests the model's
completed sacct jobs, derives a compat entry, prints a per-pool validation table, and
(with --write) upserts it into compat/euler.compat.yaml.
"""
from __future__ import annotations
import argparse
import datetime
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # put lib/ on path
import yaml                                   # noqa: E402
from ck import site as cksite                 # noqa: E402
from ck import probe as ckprobe               # noqa: E402


def _kit_root() -> str:
    return os.environ.get("KIT") or os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _fmt_h(sec):
    return "-" if not sec else f"{sec/3600:.2f}h"


def cmd_probe(argv) -> int:
    ap = argparse.ArgumentParser(prog="campaign probe")
    ap.add_argument("--model", required=True)
    ap.add_argument("--days", type=int, default=14, help="sacct lookback window (history-harvest)")
    ap.add_argument("--jobname", default=None, help="regex to match sacct JobName (overrides model's probe.jobname_match)")
    ap.add_argument("--write", action="store_true", help="upsert the result into compat/euler.compat.yaml")
    # cold-probe (calibration) options — for a model with no history
    ap.add_argument("--calibrate", action="store_true", help="submit tiny calibration jobs instead of harvesting history")
    ap.add_argument("--units", type=int, default=2, help="calibrate: units per calibration job")
    ap.add_argument("--pools", default=None, help="calibrate: comma list restricting candidate pools")
    ap.add_argument("--tier", default="4h", help="calibrate: walltime tier for calibration jobs")
    ap.add_argument("--timeout", type=int, default=3600, help="calibrate: seconds to wait per job")
    ap.add_argument("--target", default=None, help="calibrate: target pdb for the run (gens)")
    ap.add_argument("--rna", default=None, help="calibrate: rna sequence")
    ap.add_argument("--gres", default=None, help="calibrate: override gres, e.g. gpumem:8g (schedule on more/idle nodes)")
    ap.add_argument("--dry-run", action="store_true", help="calibrate: write scripts, do not submit")
    a = ap.parse_args(argv)

    kit = _kit_root()
    S = cksite.load_site(kit)
    inv = cksite.resolve(_load(os.path.join(kit, "cluster", "euler.yaml")), S)
    mpath = os.path.join(kit, "models", a.model + ".model.yaml")
    if not os.path.exists(mpath):
        print(f"no model spec: {mpath}", file=sys.stderr)
        return 1
    model = cksite.resolve(_load(mpath), S)

    if a.calibrate:
        return _calibrate(kit, S, inv, model, a)

    pattern = a.jobname or (model.get("probe") or {}).get("jobname_match") or a.model
    rx = re.compile(pattern)
    user = S.get("user") or os.environ.get("USER")

    rows = ckprobe.harvest(user, a.days)
    mine = [r for r in rows if r.get("name") and rx.search(r["name"])]
    caps = {p["name"]: p.get("caps", {}) for p in inv["pools"]}

    print(f"probe --model {a.model}  user={user}  window={a.days}d  jobname~/{pattern}/")
    print(f"  matched {len(mine)} jobs "
          f"({sum(1 for r in mine if r.get('state')=='COMPLETED')} COMPLETED)")

    by = ckprobe.summarize_by_pool(mine, inv["pools"])
    if by:
        print(f"  {'pool':<8}{'done':>5}{'fail':>5}{'VRAM_p95':>10}{'med_runtime':>13}")
        for pl, b in sorted(by.items()):
            vram = f"{b['vram_p95']:.1f}G" if b['vram_p95'] else "-"
            print(f"  {pl:<8}{b['done']:>5}{b['fail']:>5}{vram:>10}{_fmt_h(b['med_runtime_s']):>13}")

    entry = ckprobe.derive(model, mine, inv["pools"], caps)
    if not entry:
        print(f"  -> no COMPLETED jobs matched. Widen --days, fix --jobname, or cold-probe:\n"
              f"       campaign probe --model {a.model} --calibrate --write", file=sys.stderr)
        return 2

    print("  derived compat entry:")
    print("    " + json.dumps(entry, indent=2).replace("\n", "\n    "))

    if a.write:
        cp = os.path.join(kit, "compat", "euler.compat.yaml")
        doc = _load(cp) if os.path.exists(cp) else {}
        doc = doc or {}
        doc.setdefault("_meta", {"cluster": inv.get("cluster", "euler")})
        doc["_meta"]["probed_at"] = str(datetime.date.today())
        entry["measured_at"] = str(datetime.date.today())
        entry["source"] = f"history-harvest ({a.days}d)"
        doc[a.model] = entry
        with open(cp, "w") as fh:
            yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)
        print(f"  wrote {a.model} -> {cp}")
    else:
        print("  (dry run — add --write to save into compat/euler.compat.yaml)")
    return 0


def _calibrate(kit, S, inv, model, a) -> int:
    """Cold-probe: submit a tiny calibration job per candidate pool, wait, measure, derive."""
    import json as _json
    import subprocess as sp
    from ck import executor, probe as ckp

    pools = ckp.candidate_pools(model, inv)
    if a.pools:
        want = set(a.pools.split(","))
        pools = [p for p in pools if p["name"] in want]
    if not pools:
        print("no candidate pools (check compat_hints.min_arch or --pools)", file=sys.stderr)
        return 1
    caps = {p["name"]: p.get("caps", {}) for p in inv["pools"]}
    rb = S.get("rb_root") or kit
    logdir = os.path.join(rb, "ops", "logs", "campaign_kit", "_calib")
    outdir = os.path.join(kit, "plans", "_calib")
    os.makedirs(logdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    min_arch = (model.get("compat_hints") or {}).get("min_arch") or "any"
    print(f"cold-probe --model {model['name']}: {a.units}u calibration on "
          f"{[p['name'] for p in pools]} (arch >= {min_arch}), tier={a.tier}")

    jobs = []
    for p in pools:
        name, text = executor.calibration_sbatch(model, p, a.tier, a.units, logdir,
                                                  a.target or "", a.rna or "", gres_override=a.gres)
        path = os.path.join(outdir, name + ".sbatch")
        with open(path, "w") as fh:
            fh.write(text)
        if a.dry_run:
            print(f"  [dry] {p['name']:<7} -> {path}")
            continue
        out = sp.run(["sbatch", "--parsable", path], capture_output=True, text=True)
        jid = out.stdout.strip()
        if jid:
            print(f"  submitted {p['name']:<7} -> job {jid}")
            jobs.append((p, jid))
        else:
            print(f"  FAILED to submit {p['name']}: {out.stderr.strip()[:100]}")

    if a.dry_run:
        print("\n  (dry run — calibration scripts written, not submitted. Drop --dry-run to run the cold-probe.)")
        return 0

    meas = []
    for p, jid in jobs:
        print(f"  waiting on {p['name']} (job {jid}, up to {a.timeout}s) ...")
        r = ckp.poll_until_done(jid, a.timeout)
        meas.append({"pool": p["name"], "state": r.get("state"), "gpumem_g": r.get("gpumem_g"),
                     "gpuutil": r.get("gpuutil"), "elapsed_s": r.get("elapsed_s")})
        print(f"    {p['name']:<7} {r.get('state')}  vram={r.get('gpumem_g')}G  elapsed={r.get('elapsed_s')}s")

    entry = ckp.derive_from_calibration(model, meas, inv, caps, a.units)
    if not entry:
        print("  no pool COMPLETED — model failed everywhere (check run_contract + the _calib logs).", file=sys.stderr)
        return 2
    print("  derived compat entry:")
    print("    " + _json.dumps(entry, indent=2).replace("\n", "\n    "))
    if a.write:
        cp = os.path.join(kit, "compat", "euler.compat.yaml")
        doc = _load(cp) or {}
        doc.setdefault("_meta", {"cluster": inv.get("cluster", "euler")})
        doc["_meta"]["probed_at"] = str(datetime.date.today())
        entry["measured_at"] = str(datetime.date.today())
        doc[model["name"]] = entry
        with open(cp, "w") as fh:
            yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)
        print(f"  wrote {model['name']} -> {cp}")
    return 0


def cmd_plan(argv) -> int:
    ap = argparse.ArgumentParser(prog="campaign plan")
    ap.add_argument("campaign")
    ap.add_argument("--dry-run", action="store_true", help="(default; launch is not wired yet)")
    a = ap.parse_args(argv)
    kit = _kit_root()
    S = cksite.load_site(kit)
    inv = cksite.resolve(_load(os.path.join(kit, "cluster", "euler.yaml")), S)
    cpath = os.path.join(kit, "campaigns", a.campaign + ".campaign.yaml")
    if not os.path.exists(cpath):
        print(f"no campaign: {cpath}", file=sys.stderr)
        return 1
    camp = cksite.resolve(_load(cpath), S)
    compat = _load(os.path.join(kit, "compat", "euler.compat.yaml")) or {}
    from ck import planner

    plan = planner.plan_campaign(camp, compat, inv)
    tgt = (camp.get("target") or {}).get("id")
    print(f"plan: {camp.get('name')}   target={tgt}   compat={camp.get('compat')}")
    cols = [("workload", 28), ("model", 13), ("partition", 11), ("qos", 18),
            ("gres", 12), ("mem/cpu", 8), ("count", 8), ("u/job", 7), ("n_jobs", 8)]
    print("  " + "".join(h.ljust(w) for h, w in cols))
    missing = []
    for w in plan:
        if w.get("error"):
            print(f"  {w['label']:<22}{w['model']:<12}!! {w['error']}")
            missing.append(w["model"])
            continue
        vals = [w["label"], w["model"], w["partition"], w["qos"], w["gres"],
                w["mem_per_cpu"], str(w["count"]), str(w["units_per_job"]), str(w["n_jobs"])]
        line = "  " + "".join(str(v).ljust(w2) for v, (h, w2) in zip(vals, cols))
        if w.get("memory_elastic"):
            line += "  [elastic]"
        print(line)
    if missing:
        print(f"\n  {len(missing)} model(s) unprobed: {', '.join(sorted(set(missing)))}")
        print(f"    -> campaign probe --model <name> --write")
    print("\n  (dry run — routing derived from measured compat.yaml; submit/bucket not wired yet)")
    return 0


def _load_all(kit, camp_name):
    """Shared loader: (site, inventory, campaign, compat, models-by-name)."""
    S = cksite.load_site(kit)
    inv = cksite.resolve(_load(os.path.join(kit, "cluster", "euler.yaml")), S)
    camp = cksite.resolve(_load(os.path.join(kit, "campaigns", camp_name + ".campaign.yaml")), S)
    compat = _load(os.path.join(kit, "compat", "euler.compat.yaml")) or {}
    models = {}
    mdir = os.path.join(kit, "models")
    for fn in os.listdir(mdir):
        if fn.endswith(".model.yaml") and not fn.startswith("_"):
            m = cksite.resolve(_load(os.path.join(mdir, fn)), S)
            models[m["name"]] = m
    return S, inv, camp, compat, models


_STATUS_MARK = {"PASS": "  ok ", "FAIL": "FAIL ", "WARN": "warn ", "INFO": "  .. "}


def cmd_doctor(argv) -> int:
    ap = argparse.ArgumentParser(prog="campaign doctor")
    ap.add_argument("campaign")
    a = ap.parse_args(argv)
    kit = _kit_root()
    S, inv, camp, compat, _ = _load_all(kit, a.campaign)
    from ck import planner, doctor
    plan = planner.plan_campaign(camp, compat, inv)
    checks = doctor.run(camp, plan, inv, S.get("user") or os.environ.get("USER"))
    print(f"doctor: {camp.get('name')}")
    for name, status, detail in checks:
        print(f"  [{_STATUS_MARK.get(status, status)}] {name:<20} {detail}")
    good = doctor.ok(checks)
    print("  => " + ("PASS — safe to launch" if good else "BLOCKED — fix FAIL(s) before launch"))
    return 0 if good else 1


def cmd_launch(argv) -> int:
    ap = argparse.ArgumentParser(prog="campaign launch")
    ap.add_argument("campaign")
    ap.add_argument("--submit", action="store_true", help="actually sbatch (default: write scripts only)")
    ap.add_argument("--role", choices=["all", "gen", "refold"], default="all", help="restrict to gen or refold workloads")
    a = ap.parse_args(argv)
    kit = _kit_root()
    S, inv, camp, compat, models = _load_all(kit, a.campaign)
    from ck import planner, doctor, executor
    plan = planner.plan_campaign(camp, compat, inv)
    if a.role != "all":
        plan = [w for w in plan if w.get("role") == a.role]
    checks = doctor.run(camp, plan, inv, S.get("user") or os.environ.get("USER"))
    if not doctor.ok(checks):
        print("launch BLOCKED by doctor:", file=sys.stderr)
        for n, st, d in checks:
            if st == "FAIL":
                print(f"  FAIL {n}: {d}", file=sys.stderr)
        return 1
    rb_root = S.get("rb_root") or kit
    outdir = os.path.join(kit, "plans", camp["name"])
    logdir = os.path.join(rb_root, "ops", "logs", "campaign_kit", camp["name"])
    res = executor.write_scripts(plan, camp, models, inv, rb_root, outdir, logdir, submit=a.submit)
    print(f"launch: {camp['name']}  ({'SUBMITTED' if a.submit else 'dry-run — scripts written, not submitted'})")
    print(f"  scripts -> {outdir}")
    for path, cmd, jobid in res:
        print(f"  {os.path.basename(path):<26} {'-> ' + jobid if jobid else cmd}")
    if not a.submit:
        print("\n  review the scripts, then: campaign launch " + camp["name"] + " --submit")
    return 0


def main(argv) -> int:
    if not argv:
        print("usage: cli.py <command> ...  (commands: probe, plan, doctor, launch)", file=sys.stderr)
        return 1
    cmd, rest = argv[0], argv[1:]
    dispatch = {"probe": cmd_probe, "plan": cmd_plan, "doctor": cmd_doctor, "launch": cmd_launch}
    fn = dispatch.get(cmd)
    if not fn:
        print(f"unknown command: {cmd}", file=sys.stderr)
        return 1
    return fn(rest)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
