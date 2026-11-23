"""Re-derive the scaffold design space and its fold checks from primary files.

Writes three tables into `inputs/`, so the figure rebuilds without the NAS
construct folder or Euler:

  supp_scaffold_design_space.csv   one row per designed construct (12), from
                                   the SnapGene maps in LV_Cloning
  supp_scaffold_design_folds.csv   one row per PyFurnace/ROAD run, from
                                   `verify_fold.py`'s engine on Euler
  supp_scaffold_design_layouts.csv one row per nucleotide of the five drawn
                                   scaffolds, laid out by ViennaRNA's naview
  supp_scaffold_design_ordered.csv one row per PANEL (d) CELL: the ORDERED
                                   D-series construct at that copy number,
                                   re-folded against its own design target
  supp_scaffold_design_layouts_ordered.csv  one row per nucleotide of the four
                                   ORDERED 1x fragments, laid out by naview

Nothing here is transcribed from a written record. The site counts come from
counting the Csy4 extended repeat in the construct sequence; the fold numbers
come from re-folding each run's `best_seq.txt` with ViennaRNA and taking the
base-pair distance to that run's own target.

Run:  ~/.venvs/binder_map/bin/python3 code/supp/scaffold_design/make_inputs.py
"""
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from common import paths, seqfiles

NAS = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi")
ASSEMBLED = NAS / "LV_Cloning/Projects/F30_mod/Assembled_Final"
PLV009 = NAS / "LV_Cloning/pLV_NNN/pLV009/pLV009_POS_A_example.fa"
DBN = paths.INPUTS / "csy4_dotbracket"

# The Euler side of `pyfurnace_road`. verify_fold.py lives there and uses the
# same venv; this walks every run rather than the five it hard-codes.
EULER_DIR = "/cluster/work/wicky/lvalsecchi/rna_binder/pyfurnace_road"
REMOTE = r'''
import glob, os, RNA
B = "%s"
print("run,best_nt,target_nt,target_bp,bp_distance,pct_pairs,n_csy4,n_neg_loop,mfe")
CS, NEG = "GUUCACUGCCGUAUAGGCAGCUAAGAAA", "AUUUCGAUCAAUGGAUUAU"
for d in sorted(glob.glob(B + "/runs/*/")):
    n = os.path.basename(d.rstrip("/"))
    bs, tg = os.path.join(d, "best_seq.txt"), os.path.join(d, "target.txt")
    if not os.path.exists(tg):
        continue
    T = [l.strip() for l in open(tg) if l.strip()]
    targ = T[1] if len(T) > 1 else T[0]
    npair = targ.count("(")
    if not os.path.exists(bs):
        print("%%s,,%%d,%%d,,,,," %% (n, len(targ), npair)); continue
    s = [l.strip() for l in open(bs) if l.strip()][0].replace("T", "U")
    if len(s) != len(targ):
        print("%%s,%%d,%%d,%%d,,,,," %% (n, len(s), len(targ), npair)); continue
    st, mfe = RNA.fold(s)
    d2 = RNA.bp_distance(st, targ)
    pct = "" if not npair else "%%.1f" %% (100 * (1 - d2 / (2 * npair)))
    print("%%s,%%d,%%d,%%d,%%d,%%s,%%d,%%d,%%.1f"
          %% (n, len(s), len(targ), npair, d2, pct, s.count(CS), s.count(NEG), mfe))
''' % EULER_DIR


def design_space():
    """One row per SnapGene map in Assembled_Final, nothing transcribed.

    Two fields are taken from the SEQUENCE rather than from the file name,
    because the names are labels and the names are wrong in two places:

    * `copies` is `n_csy4 / 8`. pART254's name reads `(mod4)x2_(mod3)`, which
      is its composition, not its copy number: it is the 24-site, three-unit
      construct of the mod4 row.
    * `architecture` for `pART249_notcloned.dna`, whose name carries no mod,
      is read off its site-spacing fingerprint, which is mod3's exactly.
    """
    rows = []
    for f in sorted(ASSEMBLED.glob("pART*.dna")):
        name = re.match(r"pART\d+", f.name).group(0)
        mods = re.findall(r"F30_(mod\d)", f.name[len(name):])
        s = seqfiles.read(f)
        sites = seqfiles.find_all(s, seqfiles.CSY4)
        spacings = [b - a for a, b in zip(sites, sites[1:])]
        rows.append(dict(
            construct=name,
            architecture="+".join(dict.fromkeys(mods)),
            copies=len(sites) // 8,
            insert_nt=len(s),
            n_csy4=len(sites),
            site_spacings=" ".join(str(x) for x in spacings),
            cloned="no" if "notcloned" in f.name else "yes",
            source_file=f.name,
        ))

    # Fingerprint the one unnamed map against the first unit (7 gaps between
    # its 8 sites) of every map that names a single architecture.
    def unit_fp(r):
        return " ".join(r["site_spacings"].split()[:7])

    one = {r["architecture"]: unit_fp(r) for r in rows
           if r["architecture"] and "+" not in r["architecture"]}
    for r in rows:
        if r["architecture"]:
            continue
        hit = [a for a, fp in one.items() if fp == unit_fp(r)]
        assert len(hit) == 1, f"{r['construct']} matches {hit}, expected one"
        r["architecture"] = hit[0]
        r["source_file"] += "  [architecture from site-spacing fingerprint]"
    return rows


def targets():
    """The four dot-bracket design targets, parsed, not described."""
    rows = []
    for f in sorted(DBN.glob("F30_mod*_pos.dbn")):
        lines = [l.strip() for l in f.read_text().splitlines() if l.strip()]
        seq = next(l for l in lines if set(l) <= set("ACGUN"))
        ss = next(l for l in lines if set(l) <= set("().") and len(l) == len(seq))
        unit = len(ss) // 3
        rows.append(dict(
            target=f.stem,
            total_nt=len(ss),
            unit_nt=unit,
            target_bp=ss.count("("),
            n_csy4=len(seqfiles.find_all(seq, seqfiles.CSY4.replace("T", "U"))),
            # top-level helices in ONE repeat: 4 for the arrays, 1 otherwise.
            top_level_helices=_top_level(ss[:unit]),
            source_file=f.name,
        ))
    return rows


def _top_level(ss):
    depth = n = 0
    for c in ss:
        if c == "(":
            if depth == 0:
                n += 1
            depth += 1
        elif c == ")":
            depth -= 1
    return n


def folds():
    out = subprocess.run(["ssh", "-o", "ConnectTimeout=25", "euler",
                          f"cd {EULER_DIR} && ./venv/bin/python -"],
                         input=REMOTE, capture_output=True, text=True, timeout=900)
    if out.returncode or "run,best_nt" not in out.stdout:
        raise SystemExit("Euler fold check failed; the existing CSV is kept.\n"
                         + out.stderr[-2000:])
    return list(csv.DictReader(out.stdout.splitlines()))


# --- the ORDERED constructs, folded at their real copy number --------------
# Panel (d)'s twelve cells. Every one is a real designed sequence that was
# actually ordered: the D-series GenBank files of 22-JUL-2026, which are a
# straight relabel of the mod series (D1=mod1, D2=mod2, D3=mod3, D4=mod4 -
# confirmed here by the per-unit Csy4 site-spacing fingerprint, which matches
# the pART maps exactly).
#
# A construct at k copies is the first k of its three fragments, assembled
# head to tail: 8xCsy4 = f1, 16x = f1+f2, 24x = f1+f2+f3. The fragments are
# DISTINCT sequences, not one unit tripled - `design_fragments_seq.py` designed
# each in the locked context of the previous ones - so the construct is built
# by concatenation and never by repeating a unit.
#
# ⚠️ This REPLACED a thirds-cut of `runs/F30_mod*_pos/best_seq.txt`, which was
# wrong twice over: those are the FIRST-GENERATION targets of 20-JUL-2026,
# superseded two days later and present in no ordered map, and cutting a 3x
# design into thirds reconstructs a repeat that does not exist.
ORDERED = NAS / "LV_Cloning/Projects/F30_mod/WIP/F30_final_assembly"
DTARG = paths.INPUTS / "dseries_targets"
DMAP = {"mod1": "D1", "mod2": "D2", "mod3": "D3", "mod4": "D4"}
ARCHS = ("mod1", "mod2", "mod3", "mod4")


def _gb(path):
    """Sequence of a GenBank file (the repo's reader handles .gb)."""
    return seqfiles.read(path)


def assemble():
    """{mod: (assembled RNA, [cut after 1, 2, 3 fragments], [fragment nt])}.

    Asserts that the 8x / 16x / 24x file sets really are the same fragments, so
    a re-export that changed one of them cannot pass silently.
    """
    out = {}
    for mod, d in DMAP.items():
        frag = [_gb(ORDERED / f"F30_{d}-(24xCsy4)-f{k}.gb") for k in (1, 2, 3)]
        for cn, k in ((8, 1), (16, 2), (24, 3)):
            for i in range(k):
                got = _gb(ORDERED / f"F30_{d}-({cn}xCsy4)-f{i + 1}.gb")
                assert got == frag[i], f"{d}: {cn}x f{i + 1} differs from the 24x set"
        cuts = [sum(len(x) for x in frag[:k]) for k in (1, 2, 3)]
        asm = "".join(frag)
        for k in (1, 2, 3):
            n = len(seqfiles.find_all(asm[:cuts[k - 1]], seqfiles.CSY4))
            assert n == 8 * k, f"{mod} at {k} copies: {n} Csy4 sites, expected {8 * k}"
        out[mod] = (asm.replace("T", "U"), cuts, [len(x) for x in frag])
    return out


def _target(mod, cuts, fragnt):
    """The 3x design target, built from the vendored dot-bracket.

    Two kinds, because the two pairs of architectures were designed differently.
    `flanked_unit`: mod1 and mod3 never converged as one molecule, so each
    ordered fragment is 5 nt flank + the 1x design + 5 nt flank (the GenBank
    annotates both flanks) and the 3x target is that unit three times.
    `whole_3x`: mod2 and mod4 did converge whole, so the 3x design IS the
    target and the ordered fragments are its three synthesis cuts.
    """
    lines = [l.rstrip("\n") for l in (DTARG / f"{mod}.dbn").read_text().splitlines()]
    kind = next(l.split(":", 1)[1].strip() for l in lines if l.startswith("# kind:"))
    body = [l for l in lines if l and not l.startswith("#")]
    st = next(l for l in body if set(l) <= set("()."))
    mseq = next((l for l in body if set(l) <= set("ACGUN")), None)
    if kind == "flanked_unit":
        pad = (fragnt[0] - len(st)) // 2
        assert 2 * pad + len(st) == fragnt[0], (mod, "flanks are not symmetric")
        unit, useq = "." * pad + st + "." * pad, None
        if mseq is not None:
            useq = "N" * pad + mseq + "N" * pad
        assert len(set(fragnt)) == 1, (mod, "flanked fragments differ in length")
        return unit * 3, (useq * 3 if useq else None), kind
    assert kind == "whole_3x", (mod, kind)
    return st, None, kind


def ordered_folds():
    """Fold each ordered construct against its own target, on Euler."""
    asm = assemble()
    data = {}
    for mod in ARCHS:
        a, cuts, fragnt = asm[mod]
        t3, tseq, kind = _target(mod, cuts, fragnt)
        assert len(t3) == cuts[-1], (mod, len(t3), cuts[-1], "target != assembly")
        if tseq is not None:
            assert (seqfiles.find_all(tseq, seqfiles.CSY4.replace("T", "U"))
                    == seqfiles.find_all(a, seqfiles.CSY4.replace("T", "U"))), \
                (mod, "Csy4 sites out of register with the target")
        data[mod] = dict(asm=a, cuts=cuts, targ=t3, kind=kind, frag=fragnt)
    script = "DATA = " + repr(json.dumps(data)) + "\n" + ORDERED_REMOTE
    out = subprocess.run(["ssh", "-o", "ConnectTimeout=25", "euler",
                          f"{EULER_DIR}/venv/bin/python -u -"],
                         input=script, capture_output=True, text=True, timeout=1800)
    if out.returncode or "arch,copies" not in out.stdout:
        raise SystemExit("Euler fold of the ordered constructs failed; the "
                         "existing CSV is kept.\n" + out.stderr[-2000:])
    print(out.stderr.strip())
    rows = list(csv.DictReader(out.stdout.splitlines()))
    assert len(rows) == 12, f"expected 12 cells, got {len(rows)}"
    return rows


ORDERED_REMOTE = r"""
import json, sys, RNA
CS = "GUUCACUGCCGUAUAGGCAGCUAAGAAA"
D = json.loads(DATA)
def pset(x):
    o, stk = set(), []
    for i, c in enumerate(x):
        if c == "(": stk.append(i)
        elif c == ")": o.add((stk.pop(), i))
    return o
def sites(s):
    return [k for k in range(len(s)) if s.startswith(CS, k)]
print("arch,copies,nt,target_bp,bp_distance,pct_pairs,n_csy4,hairpins_ok,mfe,kind,engine")
for mod in ("mod1", "mod2", "mod3", "mod4"):
    d = D[mod]
    for k in (1, 2, 3):
        c = d["cuts"][k - 1]
        t, s = d["targ"][:c], d["asm"][:c]
        assert t.count("(") == t.count(")"), (mod, k, "cut splits a base pair")
        assert len(sites(s)) == 8 * k, (mod, k)
        npair = t.count("(")
        st, mfe = RNA.fold(s)
        bd = RNA.bp_distance(st, t)
        assert bd <= 2 * npair, (mod, k, bd, 2 * npair, "distance exceeds maximum")
        P = pset(st)
        hp = sum(1 for p in sites(s)
                 if all((p + 5 + j, p + 19 - j) in P for j in range(5)))
        print("%s,%d,%d,%d,%d,%.1f,%d,%d,%.1f,%s,ViennaRNA %s"
              % (mod, k, c, npair, bd, 100 * (1 - bd / (2.0 * npair)),
                 8 * k, hp, mfe, d["kind"], RNA.__version__))
"""


# --- the five drawn scaffolds, laid out by ViennaRNA ------------------------
# naview is RNAplot's own layout: the conventional RNA secondary-structure
# drawing, stems as ladders and loops opened out, which is the whole point of
# drawing a fold rather than an arc diagram. Of the three engines the 2.7.2
# Python module exposes, `simple` overlaps its own arms on these targets and
# `turtle` writes one polyline per element rather than one per chain, so naview
# is both the nicest and the least code.
#
# The layout is computed on Euler because ViennaRNA is not importable in
# ~/.venvs/binder_map; the coordinates are cached here so the figure redraws
# without it. Data travels INSIDE the script on stdin - the login banner breaks
# a separate scp or a second stdin stream (memory: euler-ssh-transfer-quirks).
LAYOUT_REMOTE = r'''
import json, sys, RNA
CS = "GUUCACUGCCGUAUAGGCAGCUAAGAAA"
d = json.loads(DATA)
arr = d.pop("_f30_array")[0]
st, _ = RNA.fold(arr)
pt = list(RNA.ptable(st))
sites = [i for i in range(len(arr)) if arr.startswith(CS, i)]
# One F30 unit, as the array itself folds it: the OUTERMOST closed domain
# holding exactly two Csy4 sites, taken from the repeat nearest the middle of
# the array so neither end's dangling sequence is in it.
cand = [(i - 1, pt[i] - 1) for i in range(1, len(pt)) if pt[i] > i
        and sum(1 for x in sites if i - 1 <= x < pt[i]) == 2]
outer = [c for c in cand if not any(a < c[0] and c[1] < b for a, b in cand)]
a, b = min(outer, key=lambda c: abs((c[0] + c[1]) / 2 - len(arr) / 2))
d["f30"] = [arr[a:b + 1], st[a:b + 1]]
sys.stderr.write("F30 unit = array %d-%d, %d nt\n" % (a, b, b - a + 1))
print("scaffold,i,base,x,y,pair,csy4")
for name, (seq, ss) in d.items():
    p = list(RNA.ptable(ss))
    xy = RNA.naview_xy_coordinates(ss)
    mark = [0] * len(seq)
    for i in range(len(seq)):
        if seq.startswith(CS, i):
            for k in range(i, i + len(CS)):
                mark[k] = 1
    for i in range(len(seq)):
        print("%s,%d,%s,%.3f,%.3f,%d,%d"
              % (name, i, seq[i], xy[i].X, xy[i].Y, p[i + 1] - 1, mark[i]))
'''


def _f30_array():
    """The eight-unit F30 tandem array of pLV009, as RNA.

    The repeat period is MEASURED, not assumed: the sixteen Csy4 sites sit at
    alternating 48 / 68 nt spacings, so one unit is 116 nt and the array is
    eight of them. (The README's "98 nt F30 repeat unit" is the one number here the
    sequence contradicts: 6,427 - 6,185 = 242 = 116 + 14 x 9, one whole unit
    plus the fourteen 9 nt loop swaps, not 98 + 126.)
    """
    s = seqfiles.read(PLV009)
    h = seqfiles.find_all(s, seqfiles.CSY4)
    gaps = [b - a for a, b in zip(h, h[1:])]
    assert len(h) == 16 and set(gaps) == {48, 68}, f"unexpected F30 array: {gaps}"
    period = gaps[0] + gaps[1]
    arr = s[h[0]: h[0] + period * (len(h) // 2)]
    assert len(seqfiles.find_all(arr, seqfiles.CSY4)) == 16
    return arr.replace("T", "U")


def layouts():
    """naview coordinates for the F30 unit and for one repeat of each target."""
    data = {}
    for f in sorted(DBN.glob("F30_mod*_pos.dbn")):
        lines = [l.strip() for l in f.read_text().splitlines() if l.strip()]
        seq = next(l for l in lines if set(l) <= set("ACGUN"))
        ss = next(l for l in lines if set(l) <= set("().") and len(l) == len(seq))
        n = len(ss) // 3                      # the files hold three repeats
        assert ss[:n].count("(") == ss[:n].count(")"), f"{f.name} repeat unbalanced"
        data[re.search(r"mod\d", f.stem).group(0)] = [seq[:n], ss[:n]]
    data["_f30_array"] = [_f30_array(), None]

    script = "DATA = " + repr(json.dumps(data)) + "\n" + LAYOUT_REMOTE
    out = subprocess.run(["ssh", "-o", "ConnectTimeout=25", "euler",
                          f"{EULER_DIR}/venv/bin/python -"],
                         input=script, capture_output=True, text=True, timeout=900)
    if out.returncode or "scaffold,i,base" not in out.stdout:
        raise SystemExit("Euler layout failed; the existing CSV is kept.\n"
                         + out.stderr[-2000:])
    print(out.stderr.strip())
    return list(csv.DictReader(out.stdout.splitlines()))


# --- the ORDERED 1x fragments, laid out for panel (c) ----------------------
# Panel (c) draws design TARGETS, so these are the target dot-brackets of the
# four ordered 1x fragments - the same objects panel (d)'s 1x column scores -
# laid out by the same naview the superseded panel used. Written to its OWN
# table so `supp_scaffold_design_layouts.csv`, which panel (a) and the existing
# panel (c) both read, is untouched.
#
# The FLANKS ARE DRAWN. mod1 and mod3 are `flanked_unit`, so their ordered
# fragment is 5 nt + design + 5 nt; drawing the fragment as ordered keeps all
# four panels on one definition of "1x" and makes the length printed under each
# fold (396 / 409 / 442 / 402) the same number panel (d) scores. The flanks show
# as two short unpaired tails, which is what they are.
LAYOUT_ORDERED_REMOTE = r"""
import json, sys, RNA
CS = "GUUCACUGCCGUAUAGGCAGCUAAGAAA"
d = json.loads(DATA)
print("scaffold,i,base,x,y,pair,csy4")
for mod, (seq, ss) in sorted(d.items()):
    # BOTH drawings of the same molecule, from ONE fold call, so panel (c) and
    # panel (d) cannot disagree about what the sequence folds into: the design
    # TARGET it was optimised toward, and the MFE RNAfold actually predicts.
    mfe, e = RNA.fold(seq)
    sys.stderr.write("%s bp_distance %d  mfe %.1f\n"
                     % (mod, RNA.bp_distance(mfe, ss), e))
    mark = [0] * len(seq)
    for i in range(len(seq)):
        if seq.startswith(CS, i):
            for k in range(i, i + len(CS)):
                mark[k] = 1
    for tag, st in (("target", ss), ("mfe", mfe)):
        p = list(RNA.ptable(st))
        xy = RNA.naview_xy_coordinates(st)
        for i in range(len(seq)):
            print("%s_%s,%d,%s,%.3f,%.3f,%d,%d"
                  % (mod, tag, i, seq[i], xy[i].X, xy[i].Y, p[i + 1] - 1, mark[i]))
"""


def ordered_units():
    """{mod: (1x ordered sequence, its design target)} - both real, neither cut
    from a repeat. For `whole_3x` architectures the 1x target is the first
    synthesis cut of the converged 3x design, not an independently converged
    design of its own; `flanked_unit` ones carry their 5 nt flanks."""
    asm = assemble()
    out = {}
    for mod in ARCHS:
        a, cuts, fragnt = asm[mod]
        t3, _, kind = _target(mod, cuts, fragnt)
        c = cuts[0]
        seq, ss = a[:c], t3[:c]
        assert len(seq) == len(ss) == fragnt[0], (mod, len(seq), len(ss))
        assert ss.count("(") == ss.count(")"), (mod, "1x target unbalanced")
        assert len(seqfiles.find_all(seq, seqfiles.CSY4.replace("T", "U"))) == 8, mod
        out[mod] = (seq, ss, kind)
    return out


def layouts_ordered():
    """naview coordinates for the four ORDERED 1x fragments, TWICE each.

    `<mod>_target` is the design target, `<mod>_mfe` the structure ViennaRNA
    predicts for that same ordered sequence - the two drawings panel (c) puts
    above and below each other. The MFE comes from the SAME `RNA.fold` call
    that reports the distance, and the distance is asserted against
    `supp_scaffold_design_ordered.csv`, so the two panels cannot disagree
    about what the sequence folds into.
    """
    u = ordered_units()
    data = {m: [s, ss] for m, (s, ss, _) in u.items()}
    script = "DATA = " + repr(json.dumps(data)) + "\n" + LAYOUT_ORDERED_REMOTE
    out = subprocess.run(["ssh", "-o", "ConnectTimeout=25", "euler",
                          f"{EULER_DIR}/venv/bin/python -"],
                         input=script, capture_output=True, text=True, timeout=900)
    if out.returncode or "scaffold,i,base" not in out.stdout:
        raise SystemExit("Euler layout of the ordered fragments failed; the "
                         "existing CSV is kept.\n" + out.stderr[-2000:])
    rows = list(csv.DictReader(out.stdout.splitlines()))
    for m, (s, _, _) in u.items():
        for tag in ("target", "mfe"):
            n = sum(1 for r in rows if r["scaffold"] == f"{m}_{tag}")
            assert n == len(s), f"{m}_{tag}: {n} layout rows for {len(s)} nt"
    # the distances this fold implies must be the ones panel (d) already ships
    got = {}
    for ln in out.stderr.strip().splitlines():
        f = ln.split()
        if len(f) >= 3 and f[1] == "bp_distance":
            got[f[0]] = int(f[2])
    prev = {r["arch"]: int(r["bp_distance"]) for r in
            csv.DictReader((paths.INPUTS / "supp_scaffold_design_ordered.csv")
                           .read_text().splitlines()) if r["copies"] == "1"}
    assert got == prev, f"panel (c) MFE disagrees with panel (d): {got} vs {prev}"
    print("panel (c) MFE agrees with panel (d) at 1x:", got)
    return rows


def write(rows, name):
    p = paths.INPUTS / name
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", p, f"({len(rows)} rows)")


if __name__ == "__main__":
    space = design_space()
    assert len(space) == 12, f"expected 12 designed constructs, got {len(space)}"
    write(space, "supp_scaffold_design_space.csv")
    write(targets(), "supp_scaffold_design_targets.csv")
    if "--skip-folds" in sys.argv:          # the Euler re-fold takes ~3 min
        print("skipped the Euler fold check and layout; existing CSVs kept")
    else:
        write(folds(), "supp_scaffold_design_folds.csv")
        write(ordered_folds(), "supp_scaffold_design_ordered.csv")
        write(layouts(), "supp_scaffold_design_layouts.csv")
        write(layouts_ordered(), "supp_scaffold_design_layouts_ordered.csv")
