"""Cache what arrayed_flow_dotplots draws, and check its gate against both
existing analyses of the same wells.

Reads the raw FCS of the LV2.33 / LV2.34.1 arrayed re-test (140 wells x 100,000
events, ~900 MB over the lab NAS) and writes into `working/arrayed_flow_dotplots/`:

  well_stats.tsv   140 rows: the gate hierarchy's counts, this figure's
                   mScarlet-positive percentage, and the two statistics needed
                   to check it against FlowJo.
  events.npz       mTagBFP2 and mScarlet for every single-cell event of the 24
                   wells drawn, plus the full channel set of the one well whose
                   gating hierarchy panel (a) draws.

Run: ~/.venvs/binder_map/bin/python3 code/supp/arrayed_flow_dotplots/make_inputs.py
`--force` re-reads wells already cached. The NAS flaps, so the pass is
RESUMABLE: each well is appended as it is read.

## No second FCS parser, and no second gate

`LV_FlowCytometry/analysis/lv2_33_screen.py` already parses these files and
already gates them, so it is IMPORTED here - `read_fcs`, `gate`, the channel
names, the spillover constant and both thresholds all come from it. A parser or
a gate that disagreed with the one the lab's own analysis used would be worse
than no figure. What this module adds is the per-event cache and the two
cross-checks; it decides nothing about the gating.

    Cells          FSC-A 3e4-2.5e5, SSC-A 1e4-2.2e5
    Single Cells   FSC-A / FSC-H in 1.0-1.9
    transfected    405(E)-A > 1000
    mScarlet+      561(D)-A > 500, after 0.00059 x 488(C)-A spillover

⚠️ That module's own `FCS_DIR` is DEAD: `.../expdata/flow/lvalsecchi/Matteo/
LV2_X_1_RNAbindersCandidates_85ng_48h` does not exist any more, so the script
cannot be re-run as written and its `out/wells.csv` is frozen at the 46-well
partial export it was written against (plates 5 and 6 plus two stragglers; its
docstring's "V1 to V12 drop out" is that export, not this experiment). Importing
it is unaffected - nothing touches `FCS_DIR` at import.

The COMPLETE 140-well export is `FCS_DIR` below, beside the FlowJo workspace
that `code/supp/arrayed_flow` already cites. Checked, not assumed: for every
well present in more than one copy the parsed frames are equal, so this is one
acquisition stored three times, and all 32 designs are drawable.

## The two checks this script prints

1. **Against `lv2_33_screen.py`'s own `out/wells.csv`**, well by well, on the
   96 wells it covers. Same code, same gate, so this must be an identity; it
   fails if the export this points at is not the one that table was built from.
2. **Against the ON/OFF of `code/supp/arrayed_flow`**, which is a DIFFERENT
   readout of the same wells - FlowJo's normalized mScarlet, `(561+ mean x
   561+ freq) / (405+ mean x 405+ freq)`, with all three gates taken on the
   single cells as siblings rather than nested. The thresholds behind it are
   not in any file, so they are RECOVERED: per well, the quantile of the single
   cells that reproduces FlowJo's `Freq. of Parent` exactly. They come back
   tight (see what the run prints), which is also the check on the two gates
   above them - a wrong single-cell population would make the implied threshold
   drift from well to well.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import paths  # noqa: E402

ANALYSIS = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_FlowCytometry"
                "/analysis")
sys.path.insert(0, str(ANALYSIS))
import lv2_33_screen as S  # noqa: E402  the lab's own reader and gate

FIG = Path(__file__).resolve().parent.name
WORK = paths.ROOT / "working" / FIG
FCS_DIR = Path("/Volumes/wicky/expdata/flow/lvalsecchi/LV2_34_1/Matteo"
               "/LV2_34_1_RNAbindersCandidates_85ng_48h")
FCS_STEM = FCS_DIR.name
FLOWJO = FCS_DIR.parents[1] / "LV2_34_1_Analysis.txt"

CHANNELS = {"bfp": S.BFP, "rfp": S.RFP, "gfp": S.GFP}
ROW_TARGET = S.TARGETS                     # A posA, B posB, C NegA, D Neg_noScaf
# arrayed_flow's column names for the same four targets, so the two figures'
# tables can be joined without anyone re-deriving the mapping.
AF_TARGET = {"posA": "PosA", "posB": "PosB", "NegA": "NegA",
             "Neg_noScaf": "FullNeg"}
GATE_WELL = "Plate6_A5"       # Csy4-H29A 10 ng on posA: the clearest ON population


def wells():
    """Every well in the export, as (well, plate, row, col, path)."""
    out = []
    for p in sorted(FCS_DIR.glob(f"{FCS_STEM}_Plate?_[A-D]?.fcs")):
        plate, rc = p.stem.rsplit("_", 2)[-2:]
        out.append((f"{plate}_{rc}", int(plate[-1]), rc[0], int(rc[1:]), p))
    return out


def flowjo():
    """FlowJo's per-well gate frequencies and means, indexed by well.

    The export header names the hierarchy it measured them on:
    `Cells/Single Cells/<channel>+ | Freq. of Parent` and `| Mean`.
    """
    d = pd.read_csv(FLOWJO, sep="\t")
    d.columns = ["sample", "bfp_freq", "bfp_mean", "rfp_mean", "rfp_freq",
                 "gfp_freq", "gfp_mean"]
    d["well"] = d["sample"].str.extract(r"_(Plate\d_[A-D]\d)\.fcs$")
    return d.dropna(subset=["well"]).set_index("well").drop(columns="sample")


def scan(force=False):
    """Read every well once; append as we go so a dropped mount can resume."""
    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / "well_stats.tsv"
    done = (set() if force or not out.exists()
            else set(pd.read_csv(out, sep="\t").well))
    fj = flowjo()
    todo = [w for w in wells() if w[0] not in done]
    for i, (well, plate, row, col, path) in enumerate(todo, 1):
        g = S.gate(S.read_fcs(path))
        tx = g[g[S.BFP] > 1000]
        rec = dict(well=well, plate=plate, row=row, col=col,
                   binder=f"P{plate}_C{col}", target=ROW_TARGET[row],
                   af_target=AF_TARGET[ROW_TARGET[row]],
                   n_gated=len(g), n_tx=len(tx),
                   pct_tx=100 * len(tx) / max(len(g), 1),
                   # THE figure's number: lv2_33_screen.py's own definition.
                   pct_mscarlet=100 * (tx[S.RFP] > S.RFP_POS).mean(),
                   # Background in the readout channel, from the same well.
                   rfp_p9999=float(np.percentile(g[S.RFP].values, 99.99)))
        for tag, ch in CHANNELS.items():
            v = np.sort(g[ch].values)
            want = fj.loc[well, f"{tag}_freq"] / 100.0
            # The threshold that reproduces FlowJo's Freq. of Parent exactly.
            rec[f"{tag}_implied_T"] = float(np.quantile(v, 1 - want))
            rec[f"{tag}_fj_freq"] = fj.loc[well, f"{tag}_freq"]
            rec[f"{tag}_fj_mean"] = fj.loc[well, f"{tag}_mean"]
        pd.DataFrame([rec]).to_csv(out, sep="\t", index=False, mode="a",
                                   header=not out.exists())
        print(f"  [{i}/{len(todo)}] {well}  {len(g):,} single cells", flush=True)
    d = pd.read_csv(out, sep="\t").drop_duplicates("well", keep="last")
    d = d.sort_values("well")
    d.to_csv(out, sep="\t", index=False)
    return d.set_index("well")


# ------------------------------------------------------------------- checks
def check_wells_csv(d):
    """Identity against lv2_33_screen.py's frozen out/wells.csv, where it overlaps."""
    w = pd.read_csv(ANALYSIS / "out" / "wells.csv")
    w["well"] = "Plate" + w.plate.astype(str) + "_" + w.well
    w = w.set_index("well")
    both = d.index.intersection(w.index)
    for col in ("n_gated", "n_tx", "pct_mscarlet"):
        delta = (d.loc[both, col] - w.loc[both, col]).abs().max()
        print(f"  vs out/wells.csv  {col:13s} max |diff| {delta:.3g} over "
              f"n={len(both)}")
        assert delta < 1e-6, (
            f"{col} differs from the lab's own table by {delta:g}: this is the "
            f"same code on the same gate, so the export must differ")
    return len(both)


def check_arrayed_flow(d, onoff):
    """This figure's percentage against `code/supp/arrayed_flow`'s ON/OFF.

    Two things are asked of the same wells and they are NOT the same statistic,
    so the test is whether they rank and scale together, not whether they are
    equal. Both are printed, per design, by main().
    """
    T = {tag: float(np.median(d[f"{tag}_implied_T"].dropna()))
         for tag in CHANNELS}
    for tag in CHANNELS:
        v = d[f"{tag}_implied_T"].dropna()
        print(f"  FlowJo {tag}+ threshold recovered at {T[tag]:6.1f}  "
              f"(per-well {v.min():.0f}-{v.max():.0f}, "
              f"IQR {v.quantile(.25):.0f}-{v.quantile(.75):.0f}, n={len(v)})")
    print(f"  this figure's mScarlet+ gate is {S.RFP_POS:.0f}; untransfected "
          f"P6_C3 reaches {d.loc[d.binder == 'P6_C3', 'rfp_p9999'].max():.0f} "
          f"at its 99.99th percentile")
    # FlowJo's own statistic, re-derived from its own exported columns, so the
    # join below compares like with like.
    fj_norm = ((d.rfp_fj_mean * d.rfp_fj_freq) / (d.bfp_fj_mean * d.bfp_fj_freq))
    piv = d.assign(fj_norm=fj_norm).pivot_table(
        index="binder", columns="af_target", values=["pct_mscarlet", "fj_norm"])
    cols = ["PosA", "PosB", "NegA", "FullNeg"]
    onoff_fj = (piv["fj_norm"][["PosA", "PosB"]].mean(1)
                / piv["fj_norm"][["NegA", "FullNeg"]].mean(1))
    onoff_mine = (piv["pct_mscarlet"][["PosA", "PosB"]].mean(1)
                  / piv["pct_mscarlet"][["NegA", "FullNeg"]].mean(1))
    pub = onoff.ONOFF_avg
    both = onoff_fj.index.intersection(pub.index)
    r = (onoff_fj[both] / pub[both])
    print(f"  FlowJo's columns re-derived -> ON/OFF matches the published "
          f"ONOFF_avg to {np.abs(r - 1).max() * 100:.3f}% over n={len(both)}")
    assert np.abs(r - 1).max() < 0.01, "the published CSV is not this export"
    lo = np.log2(pd.concat([onoff_mine[both], pub[both]], axis=1,
                           keys=["mine", "published"])).dropna()
    rho = float(np.corrcoef(lo["mine"].rank(), lo["published"].rank())[0, 1])
    pear = float(np.corrcoef(lo["mine"], lo["published"])[0, 1])
    print(f"  THIS figure's gate -> ON/OFF vs the published one, log2: "
          f"Pearson {pear:+.3f}, Spearman {rho:+.3f}, n={len(lo)}")
    return onoff_mine, onoff_fj, piv["pct_mscarlet"][cols]


# --------------------------------------------------------------- well choice
def choose(onoff):
    """The 6 binder columns drawn, picked by ON/OFF, never by eye.

    The three controls the experiment has, plus three designs spanning what it
    measured, each by a one-line rule on the ranked table.
    """
    c = onoff[onoff.design_id.notna()].sort_values("ONOFF_avg")
    four = c[["PosA", "PosB", "NegA", "FullNeg"]]
    # Designs whose four wells all carry signal. This EXCLUDES the two failed
    # wells arrayed_flow's README names (binder 4 posB 0.0061, binder 6 posA
    # 0.0029), so a dead well is never shown as an inactive design.
    live = c[(four.min(axis=1) / four.max(axis=1)) > 0.2]
    return [
        ("P6_C3", "untransfected", "untransfected control"),
        ("P6_C4", "Csy4-H29A 100 ng", "positive control, higher dose"),
        ("P6_C5", "Csy4-H29A 10 ng", "positive control, lower dose"),
        (c.index[-1], "proposed binder", "highest ON/OFF of the 32"),
        (live.index[len(live) // 2], "inactive design", "median ON/OFF of the 32"),
        (live.index[0], "inactive design", "lowest ON/OFF with no failed well"),
    ]


def cache_events(onoff):
    """mTagBFP2 + mScarlet for the single cells of the 24 wells drawn."""
    sel = choose(onoff)
    store = {}
    for binder, _, _ in sel:
        plate, col = int(binder[1]), int(binder.split("_C")[1])
        for row in ROW_TARGET:
            well = f"Plate{plate}_{row}{col}"
            g = S.gate(S.read_fcs(FCS_DIR / f"{FCS_STEM}_{well}.fcs"))
            store[f"{well}/bfp"] = g[S.BFP].values.astype(np.float32)
            store[f"{well}/rfp"] = g[S.RFP].values.astype(np.float32)
            print(f"  events {well} {len(g):,}", flush=True)
    # Panel (a) draws the hierarchy, so it needs the events the gates REJECT:
    # everything, before `gate` is applied.
    raw = S.read_fcs(FCS_DIR / f"{FCS_STEM}_{GATE_WELL}.fcs")
    for ch in ("FSC-A", "FSC-H", "SSC-A", S.BFP, S.RFP):
        store[f"gating/{ch}"] = raw[ch].values.astype(np.float32)
    np.savez_compressed(WORK / "events.npz", **store)
    print(f"  events.npz {len(store)} arrays, "
          f"{(WORK / 'events.npz').stat().st_size / 1e6:.1f} MB")
    return sel


def main():
    onoff = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv", index_col=0)
    print(f"scanning {len(wells())} wells from {FCS_DIR} (resumable)")
    d = scan("--force" in sys.argv)
    assert len(d) == 140, f"{len(d)} wells read, the export has 140"
    assert set(d.binder) == set(onoff.index) | {"P6_C3"}, "plate map disagrees"
    print("check 1 - the lab's own table, same code and same gate:")
    n = check_wells_csv(d)
    print(f"  (identical on all {n}; out/wells.csv is the 46-well partial export)")
    print("check 2 - the figure this one must not contradict:")
    mine, fj, pct = check_arrayed_flow(d, onoff)
    print("caching the events the figure draws")
    sel = cache_events(onoff)
    print(f"{'binder':8s} {'kind':17s} {'published':>9s} {'this gate':>9s}  why")
    for binder, kind, why in sel:
        v = onoff.ONOFF_avg.get(binder, np.nan)
        print(f"{binder:8s} {kind:17s} "
              f"{'n/a' if pd.isna(v) else f'{v:9.3f}'} "
              f"{mine.get(binder, np.nan):9.3f}  {why}")


if __name__ == "__main__":
    main()
