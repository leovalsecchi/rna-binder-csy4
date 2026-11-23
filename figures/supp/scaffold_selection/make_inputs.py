"""Collect both ONT splicing readouts from the analysis outputs, not the journal.

Source of truth is the counter's own CSV in each run folder on the NAS:

  LV1.15_ONT/RawReadsAnalysis/LV1_15_AllResults.csv            raw counts
  LV1.15_ONT/RawReadsAnalysis/LV1_15_AllResults_CorrectedpART.csv
                                                 the corrected sample->construct
                                                 assignment written after LV1.17
  LV1.27.1_ONT/RawReadsAnalysis/LV1_27_1_AllResults.csv        raw counts

Writes `inputs/supp_scaffold_selection_ont.csv`, one row per sample, with the
ratio recomputed here as removed / (intact + removed) rather than copied from
the file's own `removed_fraction`, so the denominator the figure labels is the
denominator the figure divides by. That IS the counter's convention: its own
`total` column excludes `ambiguous_splicing` (LV1.15.001 is 426 + 299 = 725
with one ambiguous read set aside), and every row is asserted against it
below. Ambiguous reads are kept as their own column, not folded in.

Run:  ~/.venvs/binder_map/bin/python3 code/supp/scaffold_selection/make_inputs.py
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd

from common import paths

NAS = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_Sequencing")
RUN1 = NAS / "LV1.15_ONT/RawReadsAnalysis"
RUN2 = NAS / "LV1.27.1_ONT/RawReadsAnalysis"

# The repeat run's five samples are not labelled with a construct anywhere in
# the counter's output: every sample is aligned to one merged reference
# (MLS12_18_1_AllSamples), so the CSV cannot say which plasmid a sample was.
# This mapping is the transfection list in journal LV1.17, in its own order.
# VERIFY against the platemap before it reaches the thesis.
RUN2_MAP = {
    "LV1.27.B.001": ("pART158", "F30, 8 units"),
    "LV1.27.B.002": ("pART252", "mod4, 1x"),
    "LV1.27.B.003": ("pART243", "mod1, 1x"),
    "LV1.27.B.004": ("pART215", "negative control, 4 hairpin positions"),
    "LV1.27.B.005": ("pART242", "negative control, 12 hairpin positions"),
}


def rows():
    raw1 = pd.read_csv(RUN1 / "LV1_15_AllResults.csv")
    cor1 = pd.read_csv(RUN1 / "LV1_15_AllResults_CorrectedpART.csv")
    raw2 = pd.read_csv(RUN2 / "LV1_27_1_AllResults.csv")

    # The correction is the whole point of this figure, so it is checked, not
    # assumed: the corrected file must carry the same counts as the raw one,
    # only a different construct label per sample.
    m = raw1.merge(cor1, on="sample_name", suffixes=("", "_c"))
    assert len(m) == 10
    assert (m.intact_count == m.intact_count_c).all()
    assert (m.removed_count == m.removed_count_c).all()

    out = []
    for _, r in m.iterrows():
        out.append(_row("LV1.15", r.sample_name, r.pART, "", r))
    for _, r in raw2.iterrows():
        c, note = RUN2_MAP[r.sample_name]
        out.append(_row("LV1.27.1", r.sample_name, c, note, r))
    return out


def _row(run, sample, construct, note, r):
    total = int(r.intact_count) + int(r.removed_count)
    assert total == int(r.total), (
        f"{sample}: intact + removed = {total} but the counter wrote {r.total}")
    return dict(run=run, sample=sample, construct=construct, note=note,
                raw_reads=int(r.raw_reads), intact=int(r.intact_count),
                removed=int(r.removed_count),
                ambiguous_splicing=int(r.ambiguous_splicing),
                ambiguous_transcript=int(r.ambiguous_transcript),
                evaluable=total,
                splicing_ratio=round(int(r.removed_count) / total, 6))


if __name__ == "__main__":
    data = rows()
    p = paths.INPUTS / "supp_scaffold_selection_ont.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)
    print("wrote", p, f"({len(data)} rows)")
    for r in data:
        print(f"  {r['run']:9s} {r['sample']:13s} {r['construct']:8s} "
              f"{r['removed']:5d}/{r['evaluable']:5d} = {r['splicing_ratio']:.3f}")
