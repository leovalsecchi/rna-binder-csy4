"""Derive bucket_mv9_hd_topup.py from bucket_midrun_v9.py by adding an ALLOWLIST."""
src = "ops/submit/bucket_midrun_v9.py"
dst = "ops/submit/bucket_mv9_hd_topup.py"
s = open(src).read()
orig = s

# 1. isolated output dir so the existing _mv9_arrays is never touched
s = s.replace('OUT = C + "/_mv9_arrays"', 'OUT = C + "/_mv9_hd_topup_arrays"')

# 2. ALLOWLIST from the topup manifest — the point of an isolated topup
anchor = 'SUBMIT = "--submit" in sys.argv'
allow = '''_ALLOW_FILE = os.environ.get("MV9_TOPUP_MANIFEST",
                            RB + "/ops/submit/mv9_hd_topup_seed_manifest.txt")
if not os.path.exists(_ALLOW_FILE):
    sys.exit("FATAL: %s missing - refusing to bucket without an allowlist" % _ALLOW_FILE)
ALLOW = set()
with open(_ALLOW_FILE) as _fh:
    for _l in _fh:
        _p = _l.split()
        if len(_p) >= 2 and _p[1] != "NA":
            ALLOW.add(_p[1])
print("ALLOWLIST: %d topup uuids (%s)" % (len(ALLOW), _ALLOW_FILE))

'''
s = s.replace(anchor, allow + anchor)

# 3. apply the allowlist to the globbed run.sbatch set
old_filter = "paths = [p for p in paths if uuid_of(p) not in EXCLUDE_UUIDS]"
new_filter = (
    "paths = [p for p in paths if uuid_of(p) not in EXCLUDE_UUIDS]\n"
    "    paths = [p for p in paths if uuid_of(p) in ALLOW]\n"
    '    print("  %s: matched %d of %d topup uuids" % (stack, len(paths), len(ALLOW)))'
)
s = s.replace(old_filter, new_filter)

# 4. distinct job name so these are never confused with the main v9 arrays
s = s.replace('"#SBATCH --job-name=mv9_%s" % stack', '"#SBATCH --job-name=mv9_%s_topup" % stack')

assert s != orig, "no substitutions applied"
for probe in ("_mv9_hd_topup_arrays", "ALLOWLIST", "matched %d of %d", "mv9_%s_topup"):
    assert probe in s, "missing: " + probe
open(dst, "w").write(s)
print("wrote", dst)
