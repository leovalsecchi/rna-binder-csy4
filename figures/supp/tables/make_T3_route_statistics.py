#!/usr/bin/env python3
"""T3 - per-route campaign statistics. Writes T3_route_statistics.{csv,md}.

Run from the repo root:
    ~/.venvs/binder_map/bin/python3 code/supp/tables/make_T3_route_statistics.py

Every number is read from a file; nothing is transcribed from prose. Sources and
their mtimes are listed in code/supp/tables/README.md.

The one trap this script exists to prevent: GPU-hours per route exist in TWO
quantities that differ by ~3.7x and must never be mixed. Both are emitted, each
in its own labelled column. See README.md section "T3".
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]          # repo root
SRC = Path(__file__).resolve().parent / "sources"
OUT = Path(__file__).resolve().parent

# Fixed route order, used everywhere. Keys are the short labels the primary
# files use; ROUTES maps each to the prose name for the rendered table.
ROUTES = [
    ("RFD3+LiMPNN", "RFdiffusion3 + LigandMPNN"),
    ("BG+LiMPNN", "BoltzGen + LigandMPNN"),
    ("BG+BIF", "BoltzGen + BoltzIF"),
    ("BD1+LiMPNN", "BoltzDesign1 + LigandMPNN"),
    ("HD+LiMPNN", "HalluDesign + LigandMPNN"),
]
# inputs/Fig3e_gpu_efficiency.csv spells the routes out; map to the short label.
BENCH_LABEL = {
    "RFD3 + LigandMPNN": "RFD3+LiMPNN",
    "BoltzGen + LigandMPNN": "BG+LiMPNN",
    "BoltzGen + BoltzIF": "BG+BIF",
    "BoltzDesign + LigandMPNN": "BD1+LiMPNN",
    "HalluDesign + LigandMPNN": "HD+LiMPNN",
}
# sacct ledger spells them differently again.
LEDGER_KEY = {
    "RFD3+LiMPNN": "rfd3__ligandmpnn",
    "BG+LiMPNN": "boltzgen__ligandmpnn",
    "BG+BIF": "boltzgen__boltzif",
    "BD1+LiMPNN": "boltzdesign__ligandmpnn",
    "HD+LiMPNN": "halludesign__ligandmpnn",
}
BG_SHARED = "boltzgen (BIF+LiMPNN shared)"
REFOLD_ROLES = ("refold_AF3", "refold_B2", "refold_ESM")

# The gate cascade has two variants; only old_bif is campaign history and only
# old_bif sums to exactly 55,706. Asserted below, not assumed.
CASCADE_VARIANT = "old_bif"
GATE_STEP = "g5_iptm"


def read_ledger(path):
    """Parse the `== stack x role` block of gpu_hours/result.txt -> {(key, role): gpu_h}."""
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("== attribution"):
            break
        m = re.match(r"  (\S.*?)\s{2,}(\S+)\s+([\d.]+)\s+(\d+)$", line)
        if m:
            out[(m.group(1).strip(), m.group(2))] = float(m.group(3))
    return out


def read_cascade(path):
    """-> {(variant, route, step): remaining}"""
    with path.open() as fh:
        return {(r["variant"], r["stack"], r["step"]): int(r["remaining"])
                for r in csv.DictReader(fh)}


def read_bench(path):
    """inputs/Fig3e_gpu_efficiency.csv -> {route: row}. Skips the '#' header block."""
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    return {BENCH_LABEL[r["stack"]]: r for r in csv.DictReader(lines)}


def count_library_routes(path):
    """Panel share, counted from the delivered 12k library itself."""
    with path.open() as fh:
        counts = {}
        for r in csv.DictReader(fh):
            counts[r["stack"]] = counts.get(r["stack"], 0) + 1
    return counts


def rate(n, d):
    """Every rate is written n/N (%). Percent to 3 s.f. below 1, else 1 dp."""
    pct = 100.0 * n / d
    return f"{n:,}/{d:,} ({pct:.4f}%)" if pct < 1 else f"{n:,}/{d:,} ({pct:.1f}%)"


def main():
    ledger = read_ledger(SRC / "result.txt")
    cascade = read_cascade(SRC / "cascade_by_stack.csv")
    bench = read_bench(ROOT / "inputs" / "Fig3e_gpu_efficiency.csv")
    panel = count_library_routes(
        ROOT / "inputs" / "csy4_binder_library_12k_v6_k30_CORRECTED.csv")

    bg_shared = ledger[(BG_SHARED, "generation")]
    rows = []
    for short, pretty in ROUTES:
        k = LEDGER_KEY[short]
        gen = ledger.get((k, "generation"), 0.0)
        refold = sum(ledger.get((k, r), 0.0) for r in REFOLD_ROLES)
        # Ledger convention for the two BoltzGen routes: each is charged the FULL
        # shared backbone job, so the two per-route figures deliberately
        # double-count it. This is NOT the 16.7/72.6/10.7 benchmark split.
        shared = bg_shared if short.startswith("BG+") else 0.0
        ledger_gpuh = gen + refold + shared

        generated = int(bench[short]["generated"])
        survivors = cascade[(CASCADE_VARIANT, short, GATE_STEP)]
        bench_dpg = 1.0 / (float(bench[short]["gpuh_gen"])
                           + float(bench[short]["gpuh_refold"]))
        rows.append(dict(
            route=pretty,
            designs_generated=generated,
            gpu_h_ledger_sacct=round(ledger_gpuh, 1),
            designs_per_gpu_h_ledger_sacct=round(generated / ledger_gpuh, 1),
            designs_per_gpu_h_bench_rtx4090=round(bench_dpg),
            gate_survivors=survivors,
            gate_pass_rate=rate(survivors, generated),
            panel_designs=panel[short],
            panel_share=rate(panel[short], sum(panel.values())),
        ))

    # --- checks: anything that fails here means a source changed. ---
    gen_total = sum(r["designs_generated"] for r in rows)
    assert gen_total == 6_992_515, gen_total
    surv_total = sum(r["gate_survivors"] for r in rows)
    assert surv_total == 55_706, surv_total
    assert surv_total == cascade[(CASCADE_VARIANT, "ALL", GATE_STEP)], "variant mismatch"
    assert sum(panel.values()) == 12_000, sum(panel.values())
    # Routes counted once = 184,053.2; the per-route column sums HIGHER by
    # exactly the shared BoltzGen job, because both BG routes carry it.
    routes_once = sum(v for (s, _), v in ledger.items() if s != "-")
    assert round(routes_once, 1) == 184_053.2, routes_once
    col_sum = sum(r["gpu_h_ledger_sacct"] for r in rows)
    assert abs(col_sum - (routes_once + bg_shared)) < 0.2, (col_sum, routes_once)
    grand = round(sum(ledger.values()), 1)
    assert grand == 197_796.5, grand
    # The benchmark reads faster than the ledger for every route, but NOT by a
    # constant factor: the ratio spans 2.4x (BoltzDesign1+LigandMPNN) to 3.8x
    # (BoltzGen+LigandMPNN). A single "~3.7x" would misdescribe the table.
    ratios = [r["designs_per_gpu_h_bench_rtx4090"]
              / r["designs_per_gpu_h_ledger_sacct"] for r in rows]
    assert all(x > 1 for x in ratios), ratios
    assert 2.3 < min(ratios) < 2.5 and 3.7 < max(ratios) < 3.9, ratios

    cols = list(rows[0])
    with (OUT / "T3_route_statistics.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    hdr = {
        "route": "Design route",
        "designs_generated": "Designs generated",
        "gpu_h_ledger_sacct": "GPU-h (campaign ledger, sacct)",
        "designs_per_gpu_h_ledger_sacct": "Designs/GPU-h (campaign ledger, sacct)",
        "designs_per_gpu_h_bench_rtx4090": "Designs/GPU-h (pinned-RTX-4090 benchmark)",
        "gate_survivors": "Gate survivors",
        "gate_pass_rate": "Gate pass rate",
        "panel_designs": "Library designs",
        "panel_share": "Library share",
    }
    fmt = lambda c, v: f"{v:,}" if isinstance(v, int) else (
        f"{v:,.1f}" if isinstance(v, float) else str(v))
    md = ["# T3 - per-route campaign statistics", "",
          f"Generated by `code/supp/tables/make_T3_route_statistics.py`. "
          f"Gate cascade variant: `{CASCADE_VARIANT}`.", "",
          "| " + " | ".join(hdr[c] for c in cols) + " |",
          "|" + "|".join("---" for _ in cols) + "|"]
    md += ["| " + " | ".join(fmt(c, r[c]) for c in cols) + " |" for r in rows]
    md += ["| **Total** | "
           + " | ".join([f"**{gen_total:,}**",
                         f"{col_sum:,.1f}[^bg]",
                         "-", "-",
                         f"**{surv_total:,}**",
                         rate(surv_total, gen_total),
                         "**12,000**", "12,000/12,000 (100.0%)"]) + " |",
           "",
           f"[^bg]: The GPU-h column sums to {col_sum:,.1f}, which exceeds the "
           f"five-route total of {routes_once:,.1f} by exactly the shared "
           f"BoltzGen generation job ({bg_shared:,.1f} GPU-h): the ledger "
           "convention charges that whole job to **each** BoltzGen route. "
           f"Total project GPU time is {grand:,.1f} GPU-h. See `README.md`, "
           "\"The shared-BoltzGen convention\" and \"Totals\".",
           "",
           "Other notes this table depends on - the two GPU-hour quantities and "
           "why they must not be mixed, the `old_bif` cascade variant, and "
           "HalluDesign's panel share being a ceiling rather than a choice - are "
           "in `README.md`."]
    (OUT / "T3_route_statistics.md").write_text("\n".join(md) + "\n")
    print(f"wrote T3_route_statistics.csv / .md  ({len(rows)} routes)")
    print(f"  ledger: routes counted once {routes_once:,.1f} GPU-h; "
          f"column sum {col_sum:,.1f} (= +{bg_shared:,.1f} shared BoltzGen); "
          f"grand total {grand:,.1f}")
    print("  bench/ledger designs-per-GPU-h ratio by route: "
          + ", ".join(f"{s} {x:.2f}x" for (s, _), x in zip(ROUTES, ratios)))


if __name__ == "__main__":
    main()
