"""Per-nucleotide contact profile for the four backbone generators.

    python code/Fig2d_binder_density/contacts_table.py [nt_contacts.csv]

Writes inputs/binder_density/nt_contacts_by_group.csv - one row per group per
nucleotide, which is all make_figure.py reads.

WHERE THE NUMBERS COME FROM.  `nt_contacts.py`, run on Euler over the SAME
design structures the renders are built from: protein->RNA heavy-atom pairs
within 4.5 A, counted per nucleotide of the 16-mer core.  4.5 A and
pair-counting to match the epitope panels elsewhere in the thesis (Fig 2d's
funnel, Fig 4g's hit footprint), not the campaign's 4.0 A `n_heavy_contacts_4A`.

Every structure counted is FULL-ATOM.  Three of the stacks write backbone-only
design files, so those designs are read from the LigandMPNN side-chain
retro-pack (`analyses/library_selection/packed_designs/`) - sequence frozen,
backbone preserved - and the rest from the design file itself.  Without that the
comparison is meaningless: on raw design files the gated means run 55 contacts
for RFD3 against 277 for HalluDesign, a 5x split that falls exactly on the
backbone-only / full-atom line.  See the README.

Coverage is not equal across groups: the retro-pack was built for a slightly
earlier gate, so some sampled designs have no packed structure.  `n` per group
is carried through to the figure and the caption for that reason.
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(os.environ.get("LV_REPO") or Path(__file__).resolve().parents[2])
GROUPS = ("RFD3", "BoltzGen", "BoltzDesign", "HalluDesign")
NT = 16


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        REPO / "inputs" / "binder_density" / "nt_contacts.csv"
    rows = list(csv.DictReader(src.open()))
    if not rows:
        raise SystemExit(f"no rows in {src}")

    tot = defaultdict(lambda: [0.0] * NT)
    n = defaultdict(int)
    src_mix = defaultdict(lambda: defaultdict(int))
    for r in rows:
        g = r["group"]
        n[g] += 1
        src_mix[g][r["source"]] += 1
        for i in range(NT):
            tot[g][i] += float(r["nt%02d" % (i + 1)])

    out = [("group", "n", "n_packed", "n_design", "nt", "mean_contacts")]
    for g in GROUPS:
        if not n[g]:
            raise SystemExit(f"no designs for {g}")
        means = [t / n[g] for t in tot[g]]
        for i, m in enumerate(means, 1):
            out.append((g, n[g], src_mix[g]["packed"], src_mix[g]["design"],
                        i, "%.4f" % m))
        print("%-12s n=%4d (packed %4d, design %4d)  total %6.1f  peak nt%d  %s"
              % (g, n[g], src_mix[g]["packed"], src_mix[g]["design"], sum(means),
                 means.index(max(means)) + 1,
                 " ".join("%.0f" % x for x in means)))

    p = REPO / "inputs" / "binder_density" / "nt_contacts_by_group.csv"
    with p.open("w", newline="") as fh:
        csv.writer(fh).writerows(out)
    print(f"wrote {p} ({len(out) - 1} rows)")


if __name__ == "__main__":
    main()
