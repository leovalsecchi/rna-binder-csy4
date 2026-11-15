import sys, os, json, datetime, yaml
sys.path.insert(0, "lib")
from ck import probe, site
jid = sys.argv[1]
pool = sys.argv[2] if len(sys.argv) > 2 else "gpu"
S = site.load_site(".")
inv = site.resolve(yaml.safe_load(open("cluster/euler.yaml")), S)
m = site.resolve(yaml.safe_load(open("models/rfd3.model.yaml")), S)
r = probe.read_peak_vram_util(jid)
state = (r.get("state") or "").split()[0] if r.get("state") else None
print("JOB %s -> state=%s vram=%s util=%s elapsed=%s" %
      (jid, state, r.get("gpumem_g"), r.get("gpuutil"), r.get("elapsed_s")))
meas = [{"pool": pool, "state": state, "gpumem_g": r.get("gpumem_g"),
         "gpuutil": r.get("gpuutil"), "elapsed_s": r.get("elapsed_s")}]
caps = {p["name"]: p.get("caps", {}) for p in inv["pools"]}
e = probe.derive_from_calibration(m, meas, inv, caps, 1)
if not e:
    print("no COMPLETED measurement -> compat NOT written")
    sys.exit(0)
e["source"] = "calibration (live, job %s)" % jid
e["measured_at"] = str(datetime.date.today())
cp = "compat/euler.compat.yaml"
doc = yaml.safe_load(open(cp)) or {}
doc["rfd3"] = e
yaml.safe_dump(doc, open(cp, "w"), sort_keys=False, default_flow_style=False)
print("WROTE rfd3 (live calibration) -> %s" % cp)
print(json.dumps(e, indent=2))
