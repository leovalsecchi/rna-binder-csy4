"""T5 - the 32 candidates: id, design route, selection metrics, splicing ratio,
fold change, arrayed on/off.

NOT BLOCKED. Every column except the two recorded NEED:s comes off a primary file.
The pooled screen was never missing: its count tables are on the NAS at
LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment/counts__*.csv and the published gate
(reads >= 75 in POS and matched NEG, splicing >= 0.07, fold change >= 10) reproduces
exactly 7 + 25 = 32 off them. That reproduction is asserted here.

TWO CONVENTIONS ARE KEPT APART, never merged into one column:

1. *Selection metrics.* ligand_RMSD and sc_DockQ for BoltzGen+BoltzIF were scored against
   BoltzIF's own Boltz re-prediction, not against the generated backbone. The STORED
   columns are what the campaign gate actually saw and are filled for all 32; the
   *_true columns are the corrected true-backbone rescore and are filled for the
   BoltzGen+BoltzIF rows only, blank elsewhere because no other route has the defect.
2. *Fold change.* fold_change_uncapped is the raw ratio; fold_change_capped50 is the
   screen report's display cap at 50, which is what gated_Q1_75Reads_OrderedHits.xlsx
   exports. Nine of the 32 are at the cap, so the two columns are not interchangeable.

scipy is unavailable in ~/.venvs/binder_map (its _propack extension fails to dlopen).
Nothing statistical is computed here, so none is needed.

Run from the repo root:
    ~/.venvs/binder_map/bin/python3 code/supp/tables/make_T5_candidates.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd

from common import paths

HERE = Path(__file__).resolve().parent
SCREEN = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_BinderDesign"
              "/LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment")

READ_FLOOR, MIN_SPLICING, MIN_FC, FC_CAP = 75, 0.07, 10.0, 50.0
PROPOSED_BINDER_FC = 1.2       # >= 1.2 proposed binder, below inactive design

ROUTE = {("RFD3", "LiMPNN"): "RFdiffusion3 + LigandMPNN",
         ("BG", "LiMPNN"): "BoltzGen + LigandMPNN",
         ("BG", "BIF"): "BoltzGen + BoltzIF",
         ("BD1", "LiMPNN"): "BoltzDesign1 + LigandMPNN",
         ("HD", "LiMPNN"): "HalluDesign + LigandMPNN"}
ROUTE_ORDER = ["RFdiffusion3 + LigandMPNN", "BoltzGen + LigandMPNN",
               "BoltzGen + BoltzIF", "BoltzDesign1 + LigandMPNN",
               "HalluDesign + LigandMPNN"]


# --- the pooled screen, from the count tables ------------------------------
def pooled(sample):
    parts = [pd.read_csv(SCREEN / f"counts__{sample}_rep{r}.csv",
                         usecols=["transcript", "intact_count", "removed_count", "total"])
             .set_index("transcript") for r in (1, 2)]
    m = parts[0] + parts[1]
    assert len(m) == 12000, f"{sample}: {len(m)} members, expected 12,000"
    m["splicing"] = np.where(m.total > 0, m.removed_count / m.total, np.nan)
    return m


def screen():
    frames = []
    for ctx in ("A", "B"):
        j = pooled(f"POS_{ctx}").join(pooled(f"NEG_{ctx}"), rsuffix="_neg")
        # The report's own pseudocount: NEG splicing = 1/NEG reads when nothing spliced.
        neg_spl = np.where(j.removed_count_neg > 0, j.splicing_neg,
                           np.where(j.total_neg > 0, 1.0 / j.total_neg, np.nan))
        j["neg_splicing"] = neg_spl
        j["fold_change"] = j.splicing / neg_spl
        j["screen_context"] = ctx
        frames.append(j)
    pairs = pd.concat(frames)
    hits = pairs[(pairs.total >= READ_FLOOR) & (pairs.total_neg >= READ_FLOOR)
                 & (pairs.splicing >= MIN_SPLICING) & (pairs.fold_change >= MIN_FC)]
    n = hits.screen_context.value_counts()
    assert (len(hits), n.get("A"), n.get("B")) == (32, 7, 25), (
        f"gate gives {len(hits)} (A={n.get('A')}, B={n.get('B')}), expected 32 = 7 + 25")
    # The count tables key on the alignment reference id, not the design id. The
    # reference manifest is the only place the two are tied together, and the tie
    # is per context (each context has its own 12,000 references, all distinguishable
    # within it: context_equivalence_size == 1 for every row).
    man = pd.read_csv(SCREEN / "reference_manifest.csv",
                      usecols=["reference_id", "binder_design", "context_id",
                               "context_equivalence_size"])
    man = man[man.context_id.str.startswith("POS_")]
    assert (man.context_equivalence_size == 1).all(), (
        "a reference is indistinguishable from another inside its own context; "
        "the reference -> design mapping is not one-to-one")
    man["screen_context"] = man.context_id.str[-1]
    hits = (hits.reset_index()
            .merge(man, left_on=["transcript", "screen_context"],
                   right_on=["reference_id", "screen_context"], how="left"))
    assert hits.binder_design.notna().all(), "a gated reference has no design id"
    return hits.rename(columns={"binder_design": "design_id"}).set_index("design_id")


def build():
    inp = paths.INPUTS
    cand = pd.read_csv(inp / "candidates32_mapped.csv").set_index("design_id")
    lib = pd.read_csv(inp / "csy4_binder_library_12k_v6_k30_CORRECTED.csv").set_index("design_id")
    gatem = pd.read_csv(inp / "large" / "gated_pool_v6_gate_metrics.csv.gz").set_index("design_id")
    flow = pd.read_csv(inp / "lv2_34_1_flow_on_off.csv").set_index("binder")
    bif = pd.read_csv(inp / "large" / "bif_true_per_backend.csv.gz").set_index("design_id")
    sc = screen()

    assert len(cand) == 32 and set(cand.index) == set(sc.index), (
        "candidates32_mapped.csv and the screen gate disagree on which 32 designs")
    lb, gm = lib.loc[cand.index], gatem.loc[cand.index]

    t = pd.DataFrame(index=cand.index)
    t["id"] = [f"Q1_{i:02d}" for i in cand.twist_id]
    t["design_route"] = [ROUTE[(g, f)] for g, f in zip(cand.generator, cand.inverse_folding)]
    t["binder_length_aa"] = lb.length.values
    t["fold_cluster"] = lb.fold_cluster.values
    t["epitope_cluster"] = lb.epitope_cluster.values

    # Selection metrics, in the sense each gate term actually used:
    # rmsd and sc_DockQ were OR over the two gate backends, footnb was AND over both.
    t["iptm_max_af3_esm"] = gm.iptm_max.values
    t["sc_rmsd_min_af3_esm_A"] = gm[["rmsd_af3", "rmsd_esm"]].min(axis=1).values
    t["dockq_max_af3_esm"] = gm[["dockq_af3", "dockq_esm"]].max(axis=1).values
    t["seq_run_max_aa"] = gm.seq_run.values
    t["footnb_010_max_af3_esm"] = gm[["footnb_010_af3", "footnb_010_esm"]].max(axis=1).values

    is_bif = t.design_route == "BoltzGen + BoltzIF"
    bifr = bif.reindex(t.index)
    assert bifr.loc[is_bif].notna().all().all(), "a BoltzGen+BoltzIF row has no true-backbone rescore"
    t["sc_rmsd_min_af3_esm_true_A"] = np.where(
        is_bif, bifr[["rmsd_af3_true", "rmsd_esm_true"]].min(axis=1), np.nan)
    t["dockq_max_af3_esm_true"] = np.where(
        is_bif, bifr[["dockq_af3_true", "dockq_esm_true"]].max(axis=1), np.nan)

    s = sc.loc[t.index]
    t["screen_context"] = s.screen_context.values
    t["pos_reads"] = s.total.values.astype(int)
    t["neg_reads"] = s.total_neg.values.astype(int)
    t["splicing_ratio_pos"] = s.splicing.values
    t["splicing_ratio_neg"] = s.neg_splicing.values
    t["fold_change_uncapped"] = s.fold_change.values
    t["fold_change_capped50"] = np.minimum(s.fold_change.values, FC_CAP)

    f = flow.loc[cand.flow_well]
    t["flow_well"] = cand.flow_well.values
    t["flow_on_mScarlet_per_BFP"] = f.ON_avg.values
    t["flow_off_mScarlet_per_BFP"] = f.OFF_avg.values
    t["flow_on_off"] = f.ONOFF_avg.values
    t["flow_call"] = np.where(t.flow_on_off >= PROPOSED_BINDER_FC,
                              "proposed binder", "inactive design")
    assert np.allclose(t.flow_on_off, cand.flow_FC, atol=5e-4), (
        "flow on/off from lv2_34_1_flow_on_off.csv disagrees with candidates32_mapped.flow_FC")

    t["design_route"] = pd.Categorical(t.design_route, ROUTE_ORDER, ordered=True)
    return t.sort_values(["design_route", "id"]).reset_index()[
        ["id", "design_id", "design_route", "binder_length_aa", "fold_cluster",
         "epitope_cluster", "iptm_max_af3_esm", "sc_rmsd_min_af3_esm_A",
         "dockq_max_af3_esm", "seq_run_max_aa", "footnb_010_max_af3_esm",
         "sc_rmsd_min_af3_esm_true_A", "dockq_max_af3_esm_true", "screen_context",
         "pos_reads", "neg_reads", "splicing_ratio_pos", "splicing_ratio_neg",
         "fold_change_uncapped", "fold_change_capped50", "flow_well",
         "flow_on_mScarlet_per_BFP", "flow_off_mScarlet_per_BFP", "flow_on_off",
         "flow_call"]]


HEADERS = {
    "id": "ID", "design_id": "Design", "design_route": "Design route",
    "binder_length_aa": "Len (aa)", "fold_cluster": "Fold cl.",
    "epitope_cluster": "Epitope cl.", "iptm_max_af3_esm": "ipTM (max)",
    "sc_rmsd_min_af3_esm_A": "ligand_RMSD (Å, min, stored)",
    "dockq_max_af3_esm": "sc_DockQ (max, stored)", "seq_run_max_aa": "seq_run",
    "footnb_010_max_af3_esm": "footnb_010 (max)",
    "sc_rmsd_min_af3_esm_true_A": "ligand_RMSD (Å, min, corrected)",
    "dockq_max_af3_esm_true": "sc_DockQ (max, corrected)",
    "screen_context": "Context", "pos_reads": "POS reads", "neg_reads": "NEG reads",
    "splicing_ratio_pos": "Splicing (POS)", "splicing_ratio_neg": "Splicing (NEG)",
    "fold_change_uncapped": "FC (uncapped)", "fold_change_capped50": "FC (cap 50)",
    "flow_well": "Well", "flow_on_mScarlet_per_BFP": "Flow ON",
    "flow_off_mScarlet_per_BFP": "Flow OFF", "flow_on_off": "Flow ON/OFF",
    "flow_call": "Flow call",
}
ROUND = {"iptm_max_af3_esm": 3, "sc_rmsd_min_af3_esm_A": 2, "dockq_max_af3_esm": 3,
         "footnb_010_max_af3_esm": 0, "sc_rmsd_min_af3_esm_true_A": 2,
         "dockq_max_af3_esm_true": 3, "splicing_ratio_pos": 4,
         "splicing_ratio_neg": 4, "fold_change_uncapped": 1,
         "fold_change_capped50": 1, "flow_on_mScarlet_per_BFP": 4,
         "flow_off_mScarlet_per_BFP": 4, "flow_on_off": 3}


def write(t):
    t.to_csv(HERE / "T5_candidates.csv", index=False)

    d = t.copy()
    for c, n in ROUND.items():
        d[c] = d[c].round(n)
    d = d.astype(object).where(d.notna(), "-")
    lines = [
        "# T5 - the 32 candidates",
        "",
        "Selected on the pooled screen at a 75-read floor (POS and matched NEG), "
        "splicing ratio >= 0.07 and fold change >= 10; 7 in context A, 25 in context B. "
        "**A selection on ranked splicing ratio and fold change, not a hit rate.**",
        "",
        "`ligand_RMSD (stored)` / `sc_DockQ (stored)` are what the campaign gate saw. For "
        "BoltzGen+BoltzIF those were scored against a Boltz re-prediction rather than "
        "the generated backbone, so the `corrected` columns carry the true-backbone "
        "rescore for those rows and are blank for every other route. **The two "
        "conventions are never mixed in one column.** `FC (uncapped)` is the raw ratio; "
        "`FC (cap 50)` is the screen report's export cap.",
        "",
        "Flow on/off is normalised mScarlet (mScarlet mean x %mScarlet+ / BFP mean x "
        "%BFP+) averaged over the two positive targets, divided by the same averaged "
        "over the negative target and the full negative; one well per condition, no "
        "replicates. It is a reporter proxy for receptor dimerisation and for avidity "
        "on a multivalent array, **not** a measure of affinity - no Kd was measured. "
        f"`proposed binder` is ON/OFF >= {PROPOSED_BINDER_FC}, below it `inactive design`.",
        "",
        "| " + " | ".join(HEADERS[c] for c in t.columns) + " |",
        "|" + "|".join("---" for _ in t.columns) + "|",
    ]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in t.columns) + " |")
    (HERE / "T5_candidates.md").write_text("\n".join(lines) + "\n")


def report(t):
    print(f"T5: {len(t)} rows")
    print(t.design_route.value_counts().reindex(ROUTE_ORDER).to_string())
    print("context:", t.screen_context.value_counts().to_dict())
    print(f"at the FC cap of {FC_CAP:.0f}: "
          f"{(t.fold_change_uncapped >= FC_CAP).sum()}/{len(t)}")
    print(f"flow on/off >= {PROPOSED_BINDER_FC} (proposed binders): "
          f"{(t.flow_call == 'proposed binder').sum()}/{len(t)} "
          f"({100 * (t.flow_call == 'proposed binder').mean():.1f}%)")
    print("flow on/off range: %.3f to %.3f" % (t.flow_on_off.min(), t.flow_on_off.max()))
    print("POS reads range: %d to %d; splicing %.3f to %.3f" % (
        t.pos_reads.min(), t.pos_reads.max(),
        t.splicing_ratio_pos.min(), t.splicing_ratio_pos.max()))
    bif = t[t.design_route == "BoltzGen + BoltzIF"]
    print(f"BoltzGen+BoltzIF rows {len(bif)}: median stored sc_DockQ "
          f"{bif.dockq_max_af3_esm.median():.3f} -> corrected "
          f"{bif.dockq_max_af3_esm_true.median():.3f}; stored ligand_RMSD "
          f"{bif.sc_rmsd_min_af3_esm_A.median():.2f} -> "
          f"{bif.sc_rmsd_min_af3_esm_true_A.median():.2f} A")
    print(f"  of those, still passing the gate on corrected values: "
          f"{int(((bif.dockq_max_af3_esm_true >= 0.7) & (bif.sc_rmsd_min_af3_esm_true_A <= 3)).sum())}"
          f"/{len(bif)}")


if __name__ == "__main__":
    t = build()
    write(t)
    report(t)
