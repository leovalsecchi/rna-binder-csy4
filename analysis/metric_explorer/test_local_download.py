#!/usr/bin/env python3
"""Unit tests for local_download — the bulk copy engine (no real ssh).

Run:  python3 test_local_download.py      (plain asserts, no pytest needed)

We point DOWNLOADS_ROOT at a temp dir and stub nas_collections._extract_member so "CIF
extraction" just drops a placeholder file — so the batch mechanics (name validation, job
lifecycle, idempotent skip, metrics.csv, missing-refold handling) are exercised without a cluster.
"""
import sys
import time
import tempfile
from pathlib import Path

import pandas as pd

import local_download as ld
import nas_collections as nc

_FAILS = []


def check(name, cond):
    print(("  ok  " if cond else "FAIL  ") + name)
    if not cond:
        _FAILS.append(name)


def expect_raises(name, fn):
    try:
        fn()
        check(name + " (should raise)", False)
    except ValueError:
        check(name, True)
    except Exception as e:                               # noqa: BLE001
        check(name + f" (wrong exc {type(e).__name__})", False)


def wait_done(job_id, timeout=5.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        st = ld.job_status(job_id)
        if st and st["status"] in ("done", "error"):
            return st
        time.sleep(0.02)
    return ld.job_status(job_id)


def mk(i, drop=None):
    did = f"csy4_design_{i}"
    refs = {"af3": f"af3_{i}.tar.zst!{did}_af3.cif",
            "b2":  f"b2_{i}.tar.zst!{did}_b2.cif",
            "esm": f"esm_{i}.tar.zst!{did}_esm.cif"}
    if drop:
        refs[drop] = ""
    row = {"design_id": did, "stack": "BG+LiMPNN", "consensus_iptm": 0.5 + 0.01 * i}
    return {"design_id": did, "refs": refs, "row": row}


def main():
    # ---- name validation (reuses nas_collections rules) ------------------------------
    check("accepts normal name", ld.sanitize_name("My Set 1") == "My Set 1")
    expect_raises("rejects empty name", lambda: ld.sanitize_name(""))
    expect_raises("rejects '..'", lambda: ld.sanitize_name(".."))
    expect_raises("rejects path separator", lambda: ld.sanitize_name("a/b"))
    expect_raises("rejects leading dot", lambda: ld.sanitize_name(".hidden"))
    check("CAP is 500", ld.CAP == 500)

    with tempfile.TemporaryDirectory() as td:
        ld.DOWNLOADS_ROOT = Path(td) / "downloads"
        calls = {"n": 0}

        def fake_extract(ref, dest):                     # stub the ssh tar.zst extraction
            calls["n"] += 1
            Path(dest).parent.mkdir(parents=True, exist_ok=True)
            Path(dest).write_text("CIF " + Path(dest).name)
            return True
        nc._extract_member = fake_extract

        # ---- happy path: 3 designs x 3 refolds + metrics.csv ------------------------
        designs = [mk(0), mk(1), mk(2)]
        jid = ld.start("Set A", designs)
        st = wait_done(jid)
        check("job completes", st["status"] == "done")
        check("done == total == 3", st["done"] == 3 and st["total"] == 3)
        check("9 cifs copied", st["cifs"] == 9)
        check("no missing", st["missing"] == [])
        cifdir = Path(td) / "downloads" / "Set A" / "cif"
        check("af3 file on disk", (cifdir / "csy4_design_0_af3.cif").exists())
        check("all 9 cif files present", len(list(cifdir.glob("*.cif"))) == 9)
        mcsv = Path(td) / "downloads" / "Set A" / "metrics.csv"
        check("metrics.csv written", mcsv.exists())
        df = pd.read_csv(mcsv)
        check("metrics.csv has 3 rows", len(df) == 3)
        check("metrics.csv has expected cols",
              {"design_id", "stack", "consensus_iptm"}.issubset(df.columns))

        # ---- idempotent re-run: present cifs are NOT re-fetched ----------------------
        n_before = calls["n"]
        jid2 = ld.start("Set A", designs)
        st2 = wait_done(jid2)
        check("re-run completes", st2["status"] == "done")
        check("re-run skips present cifs (no new extract calls)", calls["n"] == n_before)
        check("re-run still reports 9 cifs", st2["cifs"] == 9)

        # ---- a missing refold is recorded, batch still completes ---------------------
        jid3 = ld.start("Set B", [mk(5, drop="b2")])
        st3 = wait_done(jid3)
        check("missing-refold job completes", st3["status"] == "done")
        check("one missing recorded", len(st3["missing"]) == 1)
        check("2 cifs for that design", st3["cifs"] == 2)

    print("\n" + ("ALL PASSED" if not _FAILS else f"{len(_FAILS)} FAILED: " + ", ".join(_FAILS)))
    return 1 if _FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
