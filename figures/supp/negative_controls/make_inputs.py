"""Re-derive the negative controls from the construct files and the Sanger reads.

Three tables into `inputs/`:

  supp_negative_controls_design.csv   what the negative puts where the Csy4
                                      repeat was, per construct pair
  supp_negative_controls_fold.csv     ViennaRNA MFE fold of the positive and
                                      negative scaffold inserts (run on Euler)
  supp_negative_controls_foldxy.csv   naview layout coordinates for the two
                                      ELEMENTS, so panel (b) can draw the folds
                                      instead of asserting them as numbers
  supp_negative_controls_colonies.csv one row per screened colony, with the
                                      repeat count read off its own Sanger file

Nothing is transcribed. The replacement sequence is found by walking the
flanks of each Csy4 site in the positive into the negative; the colony repeat
counts come from counting that replacement in each `.fasta` read.

Run:  ~/.venvs/binder_map/bin/python3 code/supp/negative_controls/make_inputs.py
      [--skip-fold]
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from common import paths, seqfiles

# The share mounts under EITHER root, depending on how it was connected, so the
# one that is actually there wins. A hardcoded /Volumes/wicky fails outright
# under the other one and takes every table in this file with it.
NAS = next(p for p in (Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi"),
                       Path("/Volumes/bsse/wicky/users/lvalsecchi/Leo_Valsecchi"))
           if p.is_dir())
PLV = NAS / "LV_Cloning/pLV_NNN"
# Positive, negative, the name of the pair as the screen used it, and the
# architecture. The architecture is carried as DATA from palette v0.8 on
# (2026-11-20): all four scaffold architectures are one colour now, so every
# figure that draws them has to NAME them, and a name a panel prints should not
# be a string the panel invents.
#
# A dict rather than a tuple since 2026-11-20, because the third pair needs two
# things the first two do not: its own annotation label for the F30 unit (the
# as-built map was annotated by pLannotate, which names units
# `F30 unit 1 of 4 (as-built)` and not `F30_three_way_junction_scaffold`), and a
# separate DESIGN map, because the negative it draws is an AS-BUILT sequence
# that departs from the map it was ordered as.
#
# `neg_f` is always the sequence that was actually in the well. For the first
# two pairs that IS the designed map - pART257 is itself the seven-unit clone
# the screen adopted, so its map is its as-built record. For the third it is
# not: pLV16 was ordered at eight F30 units and came back at four, and the
# whole point of drawing it is that difference.
PAIRS = [
    dict(pos="pART252", pos_f="pLV002/pLV002_part252.gb",
         neg="pART256", neg_f="pLV003/pLV003_part256.gb",
         pair="B: mod4 ring, 8 positions", arch="mod4"),
    dict(pos="pART158", pos_f="pLV001/pLV001_part158.gb",
         neg="pART257", neg_f="pLV004/pLV004_part257.gb",
         pair="A: F30 tandem array", arch="F30"),
    # Leo, 2026-11-20: "for d add the negative with less copies in the 1 by 1
    # resulting in pLV15", and on being shown that pLV15 appears zero times in
    # the transfection records: "ah yeah pLV16". So this is pLV16, the negative
    # the ARRAYED (one-by-one) screen actually used, against posA/pLV13, the
    # positive it was paired with there. See the README for what pLV15 is and
    # why it is not this bar.
    dict(pos="pLV13", pos_f="pLV013/pLV013.gb",
         neg="pLV16", neg_f="pLV016/pLV16_asbuilt_LV2-26-004c_annotated.gb",
         pair="posA / NegA, arrayed screen", arch="F30",
         neg_unit_label="F30 unit ",
         neg_design_f="pLV016/pLV16.gb"),
]

# The two colony-screen plates. Plate 49761 is journal LV1.24 (colonies 1-24),
# plate 49780 is journal LV1.28 (colonies 25-56); the reads for the second one
# were filed under LV1.30.
PLATES = [("49761", "LV1.24", NAS / "LV_Sequencing/LV1.24_NightSeq"),
          ("49780", "LV1.28", NAS / "LV_Sequencing/LV1.30_NightSeq/Sequencing_Results")]

FLANK = 18


F30MOD = NAS / "LV_Cloning/Projects/F30_mod/Assembled_Final"
SANGER = NAS / "LV_Sequencing/LV1.16_SangerSeq"
ARCH_FROM_NAME = re.compile(r"F30_(mod\d)")


def scaffolds():
    """The twelve designed F30_mod combinations: length, and which were obtained.

    Leo: "would be good to see two failed scaffolds and how much shorter they
    were". The designed SnapGene map exists for all twelve, the three never
    obtained included, so the length question is answerable from sequence where
    no plasmid ever existed. It answers in the OTHER DIRECTION - see the README.

    Nothing here is read off a label that sequence can settle:

      copies        csy4_sites // 8, each unit presenting eight sites, and
                    cross-checked against the `)xN` in the map name where the
                    name carries one
      architecture  from the map name, EXCEPT pART249, whose file is the bare
                    `pART249_notcloned.dna`; its mod3 identity is proved by its
                    Csy4 array being contained in pART250's
      outcome       `notcloned` in the file name -> the fragment was never
                    delivered; otherwise present among the LV1.16 Sanger maps
                    -> assembled; otherwise the assembly was not obtained
    """
    got = {m.group(0) for p in SANGER.glob("*.dna") if not p.name.startswith("._")
           for m in [re.search(r"pART\d+", p.name)] if m}
    arrays, rows = {}, []
    for p in sorted(F30MOD.glob("pART*.dna")):
        if p.name.startswith("._"):
            continue
        part = re.match(r"pART\d+", p.name).group(0)
        s = seqfiles.read(p)
        sites = seqfiles.find_all(s, seqfiles.CSY4)
        a, b = sites[0], sites[-1] + len(seqfiles.CSY4)
        arrays[part] = s[a:b]
        copies = len(sites) // 8
        assert len(sites) == copies * 8, f"{part}: {len(sites)} Csy4 sites"
        xn = re.findall(r"\)x(\d)", p.name)
        if xn:                               # )x2 / )x3, where the name has one
            assert copies == sum(int(k) for k in xn) + len(re.findall(
                r"\(F30_mod\d", p.name)) - len(xn), \
                f"{part}: {copies} copies from sequence, name says x{xn}"
        m = ARCH_FROM_NAME.search(p.name)
        rows.append(dict(part=part, arch=m.group(1) if m else "", copies=copies,
                         plasmid_nt=len(s), csy4_sites=len(sites),
                         array_nt=b - a,
                         outcome=("fragment not delivered" if "notcloned" in p.name
                                  else "assembled" if part in got
                                  else "assembly not obtained"),
                         source_file=p.name))
    # pART249's map name carries no architecture at all. Its one-copy mod3 array
    # has to be the first unit of pART250's two-copy array, which settles it
    # from sequence rather than from the 243-245/246-248/249-251/252-254
    # numbering convention.
    bare = [r for r in rows if not r["arch"]]
    assert [r["part"] for r in bare] == ["pART249"], f"no architecture: {bare}"
    assert arrays["pART249"] in arrays["pART250"], \
        "pART249's Csy4 array is not contained in pART250's, so it is not mod3"
    bare[0]["arch"] = "mod3"
    n = {}
    for r in rows:
        n[r["outcome"]] = n.get(r["outcome"], 0) + 1
    assert len(rows) == 12 and n == {"assembled": 9, "fragment not delivered": 1,
                                     "assembly not obtained": 2}, n
    return rows



def _span(pos, motif):
    """First hairpin position to the end of the last one."""
    return pos[-1] + len(motif) - pos[0]


def _period(pos):
    """The repeat period in nt, or 0 where the positions are not periodic.

    An F30 unit carries TWO hairpin positions, so the period is the sum of a
    short and a long gap and is the same for every unit: 48 + 68 = 116 nt in
    pART158, 39 + 59 = 98 in pART257. The mod4 ring's eight positions sit at
    irregular gaps (43 to 47), so the sum differs unit to unit and there is no
    repeat unit to count - which is the whole reason only the F30 pair can lose
    one.
    """
    g = [b - a for a, b in zip(pos, pos[1:])]
    per = {g[i] + g[i + 1] for i in range(0, len(g) - 1, 2)}
    return per.pop() if len(g) > 1 and len(per) == 1 else 0


def _labels(path, label):
    """How many features in the .gb carry a `/label=` STARTING with `label`.

    The annotation count, which is an INDEPENDENT record of the sequence counts
    above - the maps were annotated by hand in SnapGene. Asserted against the
    motif counts below, so "seven F30 units" is two sources agreeing and not
    one source read twice.

    A PREFIX since 2026-11-20, for pLV16's as-built map: pLannotate numbers each
    feature (`NegativeControlLoop 3 of 8 (as-built)`, `F30 unit 3 of 4
    (as-built)`) where SnapGene repeats one label, so an exact match counts zero
    on a file that is in fact fully annotated. The four prefixes used here -
    `Csy4 Site (Extended)`, `NegativeControlLoop`, `F30_three_way_junction_
    scaffold` and `F30 unit ` - are each unambiguous: none is a prefix of
    another label in any of these files (`F30_Stem1` and `NegCtrlLoop array` do
    not start with any of them).
    """
    return Path(path).read_text(errors="replace").count(f'/label="{label}')


def design():
    rows = []
    for spec in PAIRS:
        pos_name, pos_f = spec["pos"], spec["pos_f"]
        neg_name, neg_f = spec["neg"], spec["neg_f"]
        pair, arch = spec["pair"], spec["arch"]
        neg_unit_label = spec.get("neg_unit_label",
                                  "F30_three_way_junction_scaffold")
        P, N = seqfiles.read(PLV / pos_f), seqfiles.read(PLV / neg_f)
        sites = seqfiles.find_all(P, seqfiles.CSY4)
        reps = seqfiles.find_all(N, seqfiles.NEG_LOOP)
        # Confirm the replacement by walking one site's own flanks across.
        i = sites[0]
        a = N.find(P[i - FLANK:i])
        got = N[a + FLANK:a + FLANK + len(seqfiles.NEG_LOOP)]
        assert got == seqfiles.NEG_LOOP, f"{neg_name}: flank walk gave {got}"
        # The hairpin array, and the F30 repeat unit where there is one.
        per_p, per_n = _period(sites), _period(reps)
        units_p = len(sites) // 2 if per_p else 0
        units_n = len(reps) // 2 if per_n else 0
        # Both counts, from the annotations, against both counts from sequence.
        assert _labels(PLV / pos_f, "Csy4 Site (Extended)") == len(sites)
        assert _labels(PLV / neg_f, "NegativeControlLoop") == len(reps)
        assert _labels(PLV / pos_f, "F30_three_way_junction_scaffold") == units_p
        assert _labels(PLV / neg_f, neg_unit_label) == units_n, \
            f"{neg_name}: {_labels(PLV / neg_f, neg_unit_label)} units annotated "\
            f"as {neg_unit_label!r}, {units_n} from the sequence"
        # The shortfall, closed: whole units the negative never got, plus 9 nt
        # at every hairpin position it did get. Asserted rather than asserted in
        # prose - this is the arithmetic panel (d) draws.
        step = len(seqfiles.CSY4) - len(seqfiles.NEG_LOOP)
        assert len(P) - len(N) == (units_p - units_n) * per_p + step * len(reps)
        assert len(P) - len(N) == _span(sites, seqfiles.CSY4) \
            - _span(reps, seqfiles.NEG_LOOP), \
            f"{neg_name}: the length difference is not all inside the array"
        # What the negative was ORDERED as, where that is a different file from
        # what came back. pART256 and pART257 are their own as-built record, so
        # these three columns repeat the as-built numbers and panel (d) prints
        # nothing extra; pLV16 is not, and the difference is the bar's point.
        if "neg_design_f" in spec:
            D = seqfiles.read(PLV / spec["neg_design_f"])
            d_reps = seqfiles.find_all(D, seqfiles.NEG_LOOP)
            d_per = _period(d_reps)
            d_units = len(d_reps) // 2 if d_per else 0
            assert _labels(PLV / spec["neg_design_f"], "NegativeControlLoop") \
                == len(d_reps)
            # The design was VALENCY-MATCHED to its positive: same unit count,
            # same position count, shorter only by the 9 nt at each position.
            assert d_units == units_p and len(d_reps) == len(sites), \
                f"{neg_name}: its design is not valency-matched to {pos_name}"
            assert len(P) - len(D) == step * len(d_reps), \
                f"{neg_name}: design is not {step} nt per position shorter"
            # And the as-built lost whole units of the NEGATIVE's own period.
            assert len(D) - len(N) == (d_units - units_n) * per_n, \
                f"{neg_name}: the as-built shortfall does not close on {per_n} nt"
            design_nt, design_units, design_loops = len(D), d_units, len(d_reps)
        else:
            design_nt, design_units, design_loops = len(N), units_n, len(reps)
        rows.append(dict(
            pair=pair, arch=arch, positive=pos_name, negative=neg_name,
            positive_nt=len(P), negative_nt=len(N),
            delta_nt=len(P) - len(N),
            csy4_in_positive=len(sites), csy4_in_negative=len(seqfiles.find_all(N, seqfiles.CSY4)),
            loops_in_negative=len(reps),
            positive_array_nt=_span(sites, seqfiles.CSY4),
            negative_array_nt=_span(reps, seqfiles.NEG_LOOP),
            positive_units=units_p, negative_units=units_n,
            positive_unit_nt=per_p, negative_unit_nt=per_n,
            negative_design_nt=design_nt, negative_design_units=design_units,
            negative_design_loops=design_loops,
            nt_per_position=step,
            csy4_motif=seqfiles.CSY4, csy4_nt=len(seqfiles.CSY4),
            csy4_gc=round(seqfiles.gc(seqfiles.CSY4), 1),
            replacement=seqfiles.NEG_LOOP, replacement_nt=len(seqfiles.NEG_LOOP),
            replacement_gc=round(seqfiles.gc(seqfiles.NEG_LOOP), 1),
            source_positive=pos_f, source_negative=neg_f,
            source_negative_design=spec.get("neg_design_f", neg_f)))
    # Which pairs came back a whole repeat short, and by how many. That is the
    # comparison panel (d) exists to make, so it is checked, not assumed: the
    # mod4 ring lost none (it has no repeat unit to lose), pART257 lost one of
    # eight, pLV16 four of eight.
    lost = {r["negative"]: r["positive_units"] - r["negative_units"]
            for r in rows if r["positive_units"] > r["negative_units"]}
    assert lost == {"pART257": 1, "pLV16": 4}, lost
    return rows


def insert_span(seq, motif, pad=25):
    """First to last occurrence of `motif`, plus a pad either side."""
    p = seqfiles.find_all(seq, motif)
    return seq[max(0, p[0] - pad): p[-1] + len(motif) + pad]


def fold():
    """Fold both pairs on Euler, at two scales. ViennaRNA is not in this venv.

    Scale one is the whole scaffold insert, which is what a reader means by
    "is the negative structured". Scale two is the replaced element on its own
    (the 28 nt Csy4 repeat against the 19 nt loop), which is the only level at
    which the negative is unstructured. The two answers differ, so both are
    recorded rather than the convenient one.
    """
    items = [
        ("pART252 insert", insert_span(seqfiles.read(PLV / "pLV002/pLV002_part252.gb"),
                                       seqfiles.CSY4)),
        ("pART256 insert", insert_span(seqfiles.read(PLV / "pLV003/pLV003_part256.gb"),
                                       seqfiles.NEG_LOOP)),
        ("pART158 insert", insert_span(seqfiles.read(PLV / "pLV001/pLV001_part158.gb"),
                                       seqfiles.CSY4)),
        ("pART257 insert", insert_span(seqfiles.read(PLV / "pLV004/pLV004_part257.gb"),
                                       seqfiles.NEG_LOOP)),
        ("Csy4 repeat alone", seqfiles.CSY4),
        ("replacement alone", seqfiles.NEG_LOOP),
    ]
    # Two tables out of one ssh call. Plain rows are the fold table; rows
    # prefixed `XY` are the per-nucleotide naview layout, which only the two
    # ELEMENTS get (<= 40 nt) - a naview layout of a 938 nt insert is not a
    # drawable object, and panel (b) draws the elements only.
    script = (
        "import RNA, sys\n"
        "for name, s in %r:\n" % (items,)) + (
        "    r = s.replace('T','U')\n"
        "    st, mfe = RNA.fold(r)\n"
        "    print('%s\\t%d\\t%d\\t%.1f\\t%.4f\\t%s' % "
        "(name, len(r), st.count('('), mfe, mfe/len(r), st))\n"
        "    if len(r) > 40: continue\n"
        "    for i, c in enumerate(RNA.naview_xy_coordinates(st)[:len(r)]):\n"
        "        print('XY\\t%s\\t%d\\t%s\\t%s\\t%.3f\\t%.3f' % "
        "(name, i + 1, r[i], st[i], c.X, c.Y))\n")
    out = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=25", "euler",
         "/cluster/work/wicky/lvalsecchi/rna_binder/pyfurnace_road/venv/bin/python -"],
        input=script, capture_output=True, text=True, timeout=600)
    if out.returncode or not out.stdout.strip():
        raise SystemExit("Euler fold failed; existing CSV kept.\n" + out.stderr[-2000:])
    rows, xy = [], []
    for ln in out.stdout.strip().splitlines():
        f = ln.split("\t")
        if f[0] == "XY":
            xy.append(dict(element=f[1], i=int(f[2]), base=f[3], dot_bracket=f[4],
                           x=float(f[5]), y=float(f[6])))
            continue
        name, n, pairs, mfe, per_nt, struct = f
        rows.append(dict(insert=name, nt=int(n), predicted_pairs=int(pairs),
                         mfe_kcal_mol=float(mfe), mfe_per_nt=float(per_nt),
                         structure=struct))
    assert xy, "no naview layout came back; panel (b) cannot draw the folds"
    return rows, xy


def _loops_1mm(s, motif=seqfiles.NEG_LOOP):
    """Count the replacement allowing ONE miscalled base, non-overlapping.

    A Sanger read of a 7-unit clone is ~700 nt of near-identical repeat, so a
    single miscall inside one copy turns an even loop count odd. Counting at
    <= 1 mismatch recovers those copies; advancing past each hit stops one
    copy being counted twice.
    """
    n, i, L = 0, 0, len(motif)
    while i <= len(s) - L:
        mm = 0
        for a, b in zip(s[i:i + L], motif):
            if a != b:
                mm += 1
                if mm > 1:
                    break
        else:
            n += 1
            i += L
            continue
        i += 1
    return n


def colonies():
    """Repeat count per colony, from that colony's own Sanger read.

    One F30 unit carries TWO hairpin positions, so a clone's read holds an even
    number of replacement loops and `units = loops / 2`. Two counts are kept:

      loops       exact matches of the 19 nt replacement
      loops_1mm   matches allowing one miscalled base, non-overlapping

    `call` is `clean` when the exact count is already even, `rescued` when it
    is odd and the 1-mismatch count is even (one miscall inside one copy), and
    `ambiguous` otherwise. `units` is left empty for an ambiguous colony rather
    than rounded.
    """
    rows = []
    for plate, entry, root in PLATES:
        for f in sorted(_fastas(root)):
            m = re.match(r"LV1\.\d+\.(\d+)", f.name)
            if not m:
                continue
            s = "".join(l.strip() for l in f.read_text(errors="replace").splitlines()
                        if not l.startswith(">")).upper()
            n = len(seqfiles.find_all(s, seqfiles.NEG_LOOP))
            n1 = _loops_1mm(s)
            if n % 2 == 0:
                call, units = "clean", n // 2
            elif n1 % 2 == 0 and n1 >= n:
                call, units = "rescued", n1 // 2
            else:
                call, units = "ambiguous", ""
            rows.append(dict(plate=plate, journal_entry=entry,
                             colony=int(m.group(1)), read_nt=len(s), loops=n,
                             loops_1mm=n1, units=units, call=call,
                             source_file=f.name))
    rows.sort(key=lambda r: (r["plate"], r["colony"]))
    return rows


def _fastas(root):
    if any(root.glob("*.fasta")):
        return [p for p in root.glob("*.fasta") if not p.name.startswith("._")]
    out = []
    for d in root.iterdir():                 # LV1.24 splits over two report folders
        if d.is_dir() and d.name.startswith("Sequence Report"):
            out += [p for p in d.glob("*.fasta") if not p.name.startswith("._")]
    return out


def write(rows, name):
    p = paths.INPUTS / name
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", p, f"({len(rows)} rows)")


if __name__ == "__main__":
    write(design(), "supp_negative_controls_design.csv")
    write(scaffolds(), "supp_negative_controls_scaffolds.csv")
    col = colonies()
    assert len(col) == 56, f"expected 56 screened colonies, got {len(col)}"
    write(col, "supp_negative_controls_colonies.csv")
    if "--skip-fold" in sys.argv:
        print("skipped the Euler fold")
    else:
        f, xy = fold()
        write(f, "supp_negative_controls_fold.csv")
        write(xy, "supp_negative_controls_foldxy.csv")
