"""test_filters.py — plain-assert tests for the filter-tree evaluator.

Run:  ~/miniconda3/bin/python test_filters.py
"""
import numpy as np

import filters as F


class FakeStore:
    """Minimal Store stand-in: just n and col()."""

    def __init__(self, cols):
        self._c = {k: np.asarray(v, dtype=np.float32) for k, v in cols.items()}
        self.n = len(next(iter(self._c.values())))
        self._c["iptm_min"] = np.fmin(self._c["iptm_af3"], self._c["iptm_esm"])
        self._c["iptm_mean"] = (self._c["iptm_af3"] + self._c["iptm_esm"]) / 2.0

    def col(self, name):
        return self._c[name]


NAN = float("nan")
#                     0     1     2     3     4     5     6     7
S = FakeStore({
    "rmsd_af3":  [1.0,  1.0,  5.0,  5.0,  1.0,  NAN,  0.5,  9.0],
    "rmsd_esm":  [1.0,  5.0,  1.0,  5.0,  1.0,  1.0,  0.5,  9.0],
    "dockq_af3": [0.9,  0.9,  0.1,  0.1,  0.6,  0.9,  0.95, 0.0],
    "dockq_esm": [0.9,  0.1,  0.9,  0.1,  0.6,  0.9,  0.95, 0.0],
    "iptm_af3":  [0.9,  0.8,  0.7,  0.6,  0.5,  0.4,  0.95, 0.1],
    "iptm_esm":  [0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.99, 0.05],
    "hb_af3":    [5.0,  0.0,  5.0,  0.0,  3.0,  3.0,  9.0,  0.0],
    "hb_esm":    [5.0,  5.0,  0.0,  0.0,  3.0,  3.0,  9.0,  0.0],
    "seq_kmer":  [3.0,  3.0,  3.0,  3.0,  12.0, 3.0,  3.0,  3.0],
})

checks = 0


def ck(cond, msg):
    global checks
    checks += 1
    assert cond, msg


def idx(mask):
    return sorted(np.flatnonzero(mask).tolist())


def atom(metric, backend, cmp, value, on=True):
    return {"metric": metric, "backend": backend, "cmp": cmp, "value": value, "on": on}


def group(key, op, *atoms):
    return {"key": key, "op": op, "atoms": list(atoms)}


# ---------------------------------------------------------------- atoms & NaN

m = F.atom_mask(S, atom("rmsd", "af3", "<=", 2.0))
ck(idx(m) == [0, 1, 4, 6], f"rmsd_af3<=2 -> {idx(m)}")
ck(not m[5], "row 5 has NaN rmsd_af3 and must FAIL, not pass")

m = F.atom_mask(S, atom("dockq", "esm", ">=", 0.5))
ck(idx(m) == [0, 2, 4, 5, 6], f"dockq_esm>=0.5 -> {idx(m)}")

m = F.atom_mask(S, atom("seq_kmer", None, "<=", 6))
ck(idx(m) == [0, 1, 2, 3, 5, 6, 7], "backend-free metric resolves to its bare column")

for bad, kw in ((ValueError, {"metric": "nope", "backend": None, "cmp": "<=", "value": 1}),
                (ValueError, {"metric": "rmsd", "backend": "af3", "cmp": "~", "value": 1}),
                (ValueError, {"metric": "rmsd", "backend": "b2", "cmp": "<=", "value": 1})):
    try:
        F.atom_mask(S, kw)
        ck(False, f"expected {bad.__name__} for {kw}")
    except bad:
        ck(True, "raised")
ck("b2" not in F.BACKENDS, "Boltz-2 must not be a selectable backend")

# --------------------------------------------------------------------- tiers

# An empty tier is NEUTRAL — represented as None so no all-True array is allocated. What
# matters is that it does not filter, under AND *or* OR.
ck(F.tier_mask(S, {}) is F.NEUTRAL, "absent tier is neutral")
ck(F.tier_mask(S, {"op": "AND", "atoms": []}) is F.NEUTRAL, "empty AND tier neutral")
ck(F.tier_mask(S, {"op": "OR", "atoms": []}) is F.NEUTRAL,
   "empty OR tier must be neutral too, NOT all-False")
ck(idx(F.materialize(F.NEUTRAL, S.n)) == list(range(8)), "neutral materialises to all-True")

# combine_two: neutral drops out under AND and ABSORBS under OR (neutral means "everything",
# so 'anything OR everything' is everything). Getting the OR case wrong would silently
# narrow the result instead of widening it.
some = F.atom_mask(S, atom("rmsd", "af3", "<=", 2.0))
ck(F.combine_two(some, F.NEUTRAL, "AND") is some, "AND drops a neutral operand")
ck(F.combine_two(F.NEUTRAL, some, "AND") is some, "AND drops a neutral operand (left)")
ck(F.combine_two(some, F.NEUTRAL, "OR") is F.NEUTRAL, "OR with neutral absorbs to neutral")
ck(F.combine_two(F.NEUTRAL, F.NEUTRAL, "AND") is F.NEUTRAL, "neutral AND neutral")
m2 = F.atom_mask(S, atom("rmsd", "esm", "<=", 2.0))
ck(idx(F.combine_two(some, m2, "AND")) == [0, 4, 6], "AND of two real masks")
ck(idx(F.combine_two(some, m2, "OR")) == [0, 1, 2, 4, 5, 6], "OR of two real masks")

t = {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 2.0)]}
ck(idx(F.tier_mask(S, t)) == [0, 4, 6], "AND over both refolders")
t["op"] = "OR"
ck(idx(F.tier_mask(S, t)) == [0, 1, 2, 4, 5, 6], "OR over both refolders")

t = {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0),
                            atom("rmsd", "esm", "<=", 2.0, on=False)]}
ck(idx(F.tier_mask(S, t)) == [0, 1, 4, 6], "disabled atom drops out")

# force_backend keeps the user's own per-side thresholds
t = {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 6.0),
                            atom("seq_kmer", None, "<=", 6)]}
ck(idx(F.tier_mask(S, t, force_backend="af3")) == [0, 1, 6],
   "af3 side uses the af3 threshold and keeps backend-free atoms")
ck(idx(F.tier_mask(S, t, force_backend="esm")) == [0, 1, 2, 3, 5, 6],
   "esm side uses the esm threshold")

# ------------------------------------------------------- groups (per-metric AND/OR)

# A tier carrying only `atoms` must evaluate exactly like one group under the tier's op —
# otherwise old specs (and tier 3) would change meaning.
flat = {"op": "OR", "atoms": [atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 2.0)]}
grpd = {"op": "AND", "groups": [group("rmsd", "OR", *flat["atoms"])]}
ck(idx(F.tier_mask(S, flat)) == idx(F.tier_mask(S, grpd)), "atoms-only == one group")
ck(len(F.tier_groups(flat)) == 1 and F.tier_groups(flat)[0]["op"] == "OR",
   "atoms-only normalises to a single group carrying the tier op")
ck(len(F.tier_atoms(grpd)) == 2, "tier_atoms reaches inside groups")

# The whole point: two metrics, two different ops. rmsd OR -> {0,1,2,4,5,6},
# dockq AND -> {0,4,5,6}; tier AND -> {0,4,5,6}.
mixed = {"op": "AND", "groups": [
    group("rmsd", "OR", atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 2.0)),
    group("dockq", "AND", atom("dockq", "af3", ">=", 0.5), atom("dockq", "esm", ">=", 0.5))]}
ck(idx(F.tier_mask(S, mixed)) == [0, 4, 5, 6], f"rmsd OR + dockq AND -> {idx(F.tier_mask(S, mixed))}")
# ...and that result is reachable NO other way: one op for all four atoms gives {0,4,6} or
# {0,1,2,4,5,6}, never {0,4,5,6}.
allfour = F.tier_atoms(mixed)
ck(idx(F.tier_mask(S, {"op": "AND", "atoms": allfour})) == [0, 4, 6], "tier-wide AND")
ck(idx(F.tier_mask(S, {"op": "OR", "atoms": allfour})) == [0, 1, 2, 4, 5, 6], "tier-wide OR")

mixed["op"] = "OR"
ck(idx(F.tier_mask(S, mixed)) == [0, 1, 2, 4, 5, 6], "groups combine with the tier op")

# op_for_metric drives the plot shading — each plot must read its OWN group's op
mixed["op"] = "AND"
ck(F.op_for_metric(mixed, "rmsd") == "OR" and F.op_for_metric(mixed, "dockq") == "AND",
   "each metric reports its own op, not the tier's")
ck(F.op_for_metric(mixed, "hb") == "AND", "a metric with no group falls back to the tier op")
ck(F.op_for_metric(flat, "rmsd") == "OR", "atoms-only tier reports the tier op")

# force_backend reaches inside groups
fb = {"op": "AND", "groups": [
    group("rmsd", "AND", atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 6.0)),
    group("dockq", "AND", atom("dockq", "af3", ">=", 0.5))]}
ck(idx(F.tier_mask(S, fb, force_backend="af3")) == [0, 1, 4, 6],
   "af3 side keeps only af3 atoms, group by group")

# An empty group ABSORBS under OR — same rule as an empty tier, one level down
empty_or = {"op": "OR", "groups": [
    group("rmsd", "AND", atom("rmsd", "af3", "<=", 2.0), atom("rmsd", "esm", "<=", 2.0)),
    group("dockq", "AND", atom("dockq", "af3", ">=", 0.5, on=False))]}
ck(F.tier_mask(S, empty_or) is F.NEUTRAL, "empty group OR'd into a tier neutralises it")
empty_or["op"] = "AND"
ck(idx(F.tier_mask(S, empty_or)) == [0, 4, 6], "under AND the empty group just drops out")
ck(F.n_enabled(empty_or) == 2, "n_enabled counts only switched-on atoms across groups")

# ------------------------------------------------------------------ ipTM gate

base = np.ones(8, dtype=bool)
m, i = F.iptm_mask(S, {"on": False}, base)
ck(i["applied"] is False and idx(m) == list(range(8)), "gate off is a no-op")

m, i = F.iptm_mask(S, {"on": True, "mode": "threshold", "value": "min", "threshold": 0.7},
                   base)
ck(idx(m) == [0, 1, 6], f"min(af3,esm)>=0.7 -> {idx(m)}")

m, _ = F.iptm_mask(S, {"on": True, "mode": "threshold", "value": "af3", "threshold": 0.7},
                   base)
ck(idx(m) == [0, 1, 2, 6], "af3-only ipTM threshold")

# rank, scope=filtered: exactly k out of the base set
sub = np.zeros(8, dtype=bool)
sub[[1, 2, 3, 7]] = True
m, i = F.iptm_mask(S, {"on": True, "mode": "rank", "value": "af3", "rank_kind": "n",
                       "rank_n": 2, "scope": "filtered"}, sub)
ck(idx(m) == [1, 2], f"top-2 within the filtered set -> {idx(m)}")
ck(i["pool"] == 4 and i["k"] == 2, f"pool/k reported from the filtered set: {i}")
ck(int(m.sum()) == 2, "scope=filtered yields exactly k")

# rank, scope=global: rank everything, then AND
m, i = F.iptm_mask(S, {"on": True, "mode": "rank", "value": "af3", "rank_kind": "n",
                       "rank_n": 2, "scope": "global"}, sub)
ck(idx(m) == [], "global top-2 by af3 ipTM is rows 6,0 - neither is in the base set")
ck(i["pool"] == 8, "global pool is the whole store")

m, _ = F.iptm_mask(S, {"on": True, "mode": "rank", "value": "af3", "rank_kind": "n",
                       "rank_n": 3, "scope": "global"}, base)
ck(idx(m) == [0, 1, 6], "global top-3 by af3 ipTM")

m, i = F.iptm_mask(S, {"on": True, "mode": "rank", "value": "af3", "rank_kind": "pct",
                       "rank_pct": 25.0, "scope": "global"}, base)
ck(i["k"] == 2 and int(m.sum()) == 2, f"top 25% of 8 = 2 -> {i}")

m, _ = F.iptm_mask(S, {"on": True, "mode": "rank", "value": "af3", "rank_kind": "n",
                       "rank_n": 999, "scope": "global"}, base)
ck(int(m.sum()) == 8, "k larger than the pool keeps everything, does not error")

# ------------------------------------------------------------------- evaluate

spec = {
    "tiers": {"t1": {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0),
                                            atom("rmsd", "esm", "<=", 2.0)]},
              "t3": {"op": "AND", "atoms": [atom("seq_kmer", None, "<=", 6)]}},
    "combine": "AND",
    "iptm": {"on": False},
}
mask, info = F.evaluate(S, spec)
ck(idx(mask) == [0, 6], f"tier1 AND tier3 -> {idx(mask)}")
ck(info["n_tier1"] == 3 and info["n_tier3"] == 7 and info["n_pass"] == 2, f"info {info}")

spec["combine"] = "OR"
mask, info = F.evaluate(S, spec)
ck(idx(mask) == [0, 1, 2, 3, 4, 5, 6, 7], "tier1 OR tier3")
ck(not info["unfiltered_or"], "both tiers populated -> not the surprising OR case")

# tier1 set, tier3 empty, top-level OR => everything passes. Correct, but the user needs
# telling, or it looks like the filter broke.
spec_or = {"tiers": {"t1": {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0)]},
                     "t3": {"op": "AND", "atoms": []}},
           "combine": "OR", "iptm": {"on": False}}
mask, info = F.evaluate(S, spec_or)
ck(int(mask.sum()) == S.n, "OR with an empty tier passes everything")
ck(info["unfiltered_or"], "and that case is flagged")
spec_or["combine"] = "AND"
mask, info = F.evaluate(S, spec_or)
ck(idx(mask) == [0, 1, 4, 6] and not info["unfiltered_or"], "AND with an empty tier filters")
ck(info["neutralized"] == [], "nothing neutralised when the filter still constrains")

# The same trap one level down: tier 1's own OR with an empty metric group. The tier stops
# constraining even though four boxes are ticked in the other group — reported, not silent.
spec_g = {"tiers": {"t1": {"op": "OR", "groups": [
              group("rmsd", "AND", atom("rmsd", "af3", "<=", 2.0),
                    atom("rmsd", "esm", "<=", 2.0)),
              group("dockq", "AND", atom("dockq", "af3", ">=", 0.5, on=False))]},
                    "t3": {"op": "AND", "atoms": []}},
          "combine": "AND", "iptm": {"on": False}}
mask, info = F.evaluate(S, spec_g)
ck(int(mask.sum()) == S.n, "an empty group inside an OR passes everything")
ck(info["neutralized"] == ["t1"], f"and tier 1 is named as the cause: {info['neutralized']}")
ck(info["unfiltered_or"], "the top-level flag fires too — atoms are on but nothing constrains")
spec_g["tiers"]["t1"]["op"] = "AND"
mask, info = F.evaluate(S, spec_g)
ck(idx(mask) == [0, 4, 6] and info["neutralized"] == [] and not info["unfiltered_or"],
   "switching that group op to AND restores the filter")

spec["combine"] = "AND"
spec["iptm"] = {"on": True, "mode": "threshold", "value": "min", "threshold": 0.9}
mask, info = F.evaluate(S, spec)
ck(idx(mask) == [6], f"ipTM gate narrows to row 6 -> {idx(mask)}")
ck(info["n_base"] == 2 and info["n_pass"] == 1, f"base vs pass {info}")

# ---------------------------------------------------------------- agreement

spec = {"tiers": {"t1": {"op": "AND", "atoms": [atom("rmsd", "af3", "<=", 2.0),
                                                atom("rmsd", "esm", "<=", 2.0)]}},
        "combine": "AND", "iptm": {"on": False}}
sp = F.agreement_split(S, spec)
# af3 side passes {0,1,4,6}; esm side passes {0,2,4,5,6}
ck(sp["af3_total"] == 4 and sp["esm_total"] == 5, f"marginals {sp}")
ck(sp["both"] == 3, f"both = rows 0,4,6 -> {sp}")
ck(sp["af3_only"] == 1 and sp["esm_only"] == 2, f"disagreement {sp}")
ck(sp["both"] + sp["af3_only"] + sp["esm_only"] + sp["neither"] == S.n, "split partitions")

# ----------------------------------------------------------------- quadrants

q = F.quadrants(S, "rmsd", 2.0, 2.0, "AND")
ck(q["both"]["n"] == 3 and q["both"]["included"], f"AND keeps only 'both': {q['both']}")
ck(q["af3_only"]["n"] == 1 and not q["af3_only"]["included"], "AND excludes af3_only")
ck(q["_nan"]["n"] == 1, "row 5 has a NaN axis and is counted separately")
tot = sum(q[k]["n"] for k in ("both", "af3_only", "esm_only", "neither")) + q["_nan"]["n"]
ck(tot == S.n, f"quadrants + NaN partition the store: {tot} != {S.n}")

q = F.quadrants(S, "rmsd", 2.0, 2.0, "OR")
ck(all(q[k]["included"] for k in ("both", "af3_only", "esm_only")), "OR keeps three regions")
ck(not q["neither"]["included"], "OR still excludes 'neither'")

q = F.quadrants(S, "dockq", 0.5, 0.5, "AND")
# af3>=0.5 -> {0,1,4,5,6}; esm>=0.5 -> {0,2,4,5,6}; both -> {0,4,5,6}
ck(q["both"]["n"] == 4, f"higher-is-better metric flips the comparison: {q['both']}")
ck(q["af3_only"]["n"] == 1 and q["esm_only"]["n"] == 1, f"dockq disagreement: {q}")

pm = np.zeros(8, dtype=bool)
pm[[0, 6]] = True
q = F.quadrants(S, "rmsd", 2.0, 2.0, "AND", pass_mask=pm)
ck(q["both"]["n_passing"] == 2, "n_passing counts the current filter inside the region")

ck(abs(q["both"]["pct"] - 100.0 * 3 / 8) < 1e-9, "pct is over the full population")


# ---- variant metrics (footn): one filter row, a family of columns -------------------------
# The variant picks WHICH distance the atoms are counted within. An unknown variant must RAISE
# rather than silently fall back, or the filter would gate on a different distance than the UI
# displays — the silent-no-op failure mode this module exists to prevent.
fm = F.METRICS["footnb"]
ck(fm.column("af3") == "footnb_010_af3", "default variant resolves to 1.0 A column")
ck(fm.column("esm", "030") == "footnb_030_esm", "explicit variant resolves")
ck(fm.column("af3", "100") == "footnb_100_af3", "10.0 A tag is not confused with 1.0 A")
for _bad in ("999", "10", "", "1.0"):
    try:
        fm.column("af3", _bad)
        ck(False, f"variant {_bad!r} should have raised")
    except ValueError:
        ck(True, "unknown variant raises")
try:
    F.METRICS["dsasa"].column("af3", "010")
    ck(False, "non-variant metric accepted a variant")
except ValueError:
    ck(True, "non-variant metric rejects a variant")
_fd = fm.as_dict()
ck(len(_fd["variants"]) == 11, "11 distance variants exposed to the UI")
ck(_fd["variant_default"] == "010", "default variant is 1.0 A")
ck([v["key"] for v in _fd["variants"]][:3] == ["005", "010", "015"], "variant order")
ck(fm.tier == 3 and fm.cmp == "<=", "footnb is a tier-3 upper-bound gate")
ck("footn" not in F.METRICS, "foot-side mirror metric is not exposed")

print(f"all {checks} checks passed")
