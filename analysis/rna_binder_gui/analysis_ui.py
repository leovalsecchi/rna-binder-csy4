"""Analysis view: load an aggregated-metrics CSV, filter designs live, rank, plot.

The expensive step (aggregation) runs offline/over SSH and writes a CSV under
analyses/. Here we just read that CSV via the mount and filter it in pandas —
so dragging thresholds is instant, no re-aggregation. Logic lives in analysis.py.
"""
from __future__ import annotations

import os

import streamlit as st

import analysis as A
import opendesign

# Filter dimensions -> per-design column suffix (prefixed af3_/b2_ at use).
_RANK_COLS = ["iptm", "rmsd", "dockq", "lrmsd", "bsa", "contacts", "hbonds", "composite"]


@st.cache_data(show_spinner=False)
def _load_and_reshape(path: str, mtime: float):
    """Load long CSV -> per-design table. Cached on (path, mtime) so filter
    changes don't reload; a fresh aggregation (new mtime) invalidates it."""
    return A.to_per_design(A.load_metrics(path))


def render(mount: str):
    st.subheader("Analysis — filter & rank designs")

    csvs = A.find_metric_csvs(mount)
    if not csvs:
        st.info("No structural-metrics CSVs found under `analyses/`. "
                "Run the aggregator (`aggregate_metrics.py --structural`) first.")
        return

    base = mount.rstrip("/") + "/"
    choice = st.selectbox("Metrics CSV", csvs,
                          format_func=lambda p: p.replace(base, ""))
    try:
        mtime = os.path.getmtime(choice)
    except OSError:
        st.error("Can't stat the CSV — is the eulermount up?")
        return

    with st.spinner("Loading + reshaping…"):
        try:
            designs = _load_and_reshape(choice, mtime)
        except Exception as e:
            st.error(f"Failed to load/reshape CSV: {e}")
            return

    st.caption(f"{len(designs)} designs across {designs['stack'].nunique()} stacks "
               f"· {os.path.basename(choice)}")

    prefix = st.radio("Refolder metrics to filter on", ["af3", "b2", "esm"],
                      horizontal=True, format_func=str.upper)
    if f"{prefix}_iptm" not in designs:
        st.warning(f"This CSV has no {prefix.upper()} metrics. Pick the other refolder.")
        return

    c1, c2, c3, c4, c5, c6 = st.columns(6)
    thr = dict(
        iptm_min=c1.number_input("iPTM >", 0.0, 1.0, 0.5, 0.05),
        rmsd_max=c2.number_input("Cα RMSD <", 0.0, 30.0, 2.0, 0.5),
        dockq_min=c3.number_input("DockQ ≥", 0.0, 1.0, 0.23, 0.01,
                                  help="DockQ score (0-1); 0.23 = CAPRI acceptable. "
                                       "0 disables. (This is the score, not LRMSD.)"),
        bsa_min=c4.number_input("BSA >", 0.0, 5000.0, 800.0, 50.0),
        contacts_min=c5.number_input("contacts >", 0, 300, 10, 1),
        hbonds_min=c6.number_input("H-bonds ≥", 0, 50, 0, 1,
                                   help="0 disables the H-bond filter"),
    )

    # Per-stack pass rates (updates live with the thresholds).
    rates = A.pass_rate_table(designs, prefix, **thr)
    st.markdown("**Pass rate per stack**")
    st.dataframe(rates, hide_index=True, use_container_width=True)

    # Filtered, ranked design table.
    scored = A.add_composite(designs, prefix)
    passing = scored[A.filter_mask(designs, prefix, **thr)].sort_values(
        f"{prefix}_composite", ascending=False)
    st.markdown(f"**{len(passing)} passing designs** (ranked by composite score)")
    cols = ["stack", "design_key"] + [f"{prefix}_{c}" for c in _RANK_COLS]
    cols = [c for c in cols if c in passing]
    st.dataframe(passing[cols].round(3), hide_index=True, use_container_width=True)
    st.download_button("Download passing designs (CSV)",
                       passing[cols].to_csv(index=False),
                       file_name="passing_designs.csv", mime="text/csv")

    # Open a design's 4 variants (design + 3 refolds) in ChimeraX, superposed on
    # the binder. The GUI runs locally, so it can launch ChimeraX directly.
    openable = (passing[passing["designed_source"].notna()]
                if "designed_source" in passing else passing.iloc[:0])
    if openable.empty:
        st.caption("ℹ️ Open-in-ChimeraX needs designed-file paths — use a `--structural` CSV.")
    else:
        key2src = dict(zip(openable["design_key"], openable["designed_source"]))
        oc1, oc2 = st.columns([4, 1])
        sel = oc1.selectbox(
            "Open in ChimeraX (design + 3 refolds, superposed on binder)",
            list(key2src), key="open_sel")
        if oc2.button("🧬 Open", use_container_width=True):
            src = key2src[sel]
            if not str(src).startswith("/"):
                src = os.path.join(mount, str(src))
            try:
                opendesign.open_design(src)
                st.success(f"Opening **{sel}** — first open ~30 s, then leave ChimeraX "
                           "open and later opens are instant.")
            except BaseException as e:
                st.error(f"Open failed: {e}")

    # Plots — none of this exists in the offline pipeline, so it's the value-add.
    st.markdown("**Plots**")
    p1, p2 = st.columns(2)
    with p1:
        st.caption(f"{prefix.upper()} iPTM vs Cα-RMSD (all designs, colored by stack)")
        sc = scored.dropna(subset=[f"{prefix}_iptm", f"{prefix}_rmsd"])
        if not sc.empty:
            st.scatter_chart(sc, x=f"{prefix}_rmsd", y=f"{prefix}_iptm", color="stack")
    with p2:
        st.caption("Pass % per stack")
        bar = rates[rates["stack"] != "ALL"].set_index("stack")["pass %"]
        if not bar.empty:
            st.bar_chart(bar)
