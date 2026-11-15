#!/usr/bin/env python3
"""Unit tests for nas_collections — pure logic only (no real NAS / no ssh euler).

Run:  python3 test_nas_collections.py      (plain asserts, no pytest needed)

We point LV_ROOT at a temp dir, force nas_available() True, and stub _extract_member so
"CIF extraction" just drops a placeholder file — so the folder mechanics (sanitize, upsert /
dedup, column-union, remove-by-prefix, path-escape guards) are exercised without a cluster.
"""
import sys
import tempfile
from pathlib import Path

import pandas as pd

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


def ref(design_id, backend):
    """A ref whose member basename starts with '<design_id>_' (mirrors real names)."""
    return f"campaigns/x/refold_{backend}/msa/designs.tar.zst!designs/{design_id}_{backend.upper()}_s0.cif"


def refs(design_id):
    return {"af3": ref(design_id, "af3"), "b2": ref(design_id, "b2"), "esm": ref(design_id, "esm")}


def fake_extract(r, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(b"CIFDATA\n")
    return True


def main():
    # ---- sanitize_folder --------------------------------------------------------------
    print("sanitize_folder:")
    check("plain", nc.sanitize_folder("promising") == "promising")
    check("spaces collapsed+trimmed", nc.sanitize_folder("  a   b  ") == "a b")
    check("dots/dash/underscore ok", nc.sanitize_folder("v2.best-set_1") == "v2.best-set_1")
    check("leading underscore ok", nc.sanitize_folder("_drafts") == "_drafts")
    for bad in ["", "   ", ".", "..", "a/b", "a\\b", ".hidden", "/etc/passwd", "../x", "a\0b"]:
        expect_raises(f"reject {bad!r}", lambda b=bad: nc.sanitize_folder(b))

    # ---- _parse_ref -------------------------------------------------------------------
    print("_parse_ref:")
    a, m, b = nc._parse_ref("dir/x.tar.zst!designs/foo_AF3_s0.cif")
    check("archive", a == "dir/x.tar.zst")
    check("member", m == "designs/foo_AF3_s0.cif")
    check("basename", b == "foo_AF3_s0.cif")
    a2, m2, b2 = nc._parse_ref("loose/path/bar.cif")
    check("loose basename", b2 == "bar.cif" and m2 == "")

    # ---- validate_ref (RB confinement / option-injection guard) ------------------------
    print("validate_ref:")
    okarch, okmem = nc.validate_ref("campaigns/x/designs.tar.zst!designs/d_AF3_s0.cif")
    check("valid ref passes", okarch == "campaigns/x/designs.tar.zst" and okmem.endswith("d_AF3_s0.cif"))
    check("valid loose ref passes", nc.validate_ref("campaigns/x/y.cif") == ("campaigns/x/y.cif", ""))
    for bad in ["/etc/passwd", "../../etc/passwd", "a/../../etc/x", "x.tar.zst!/abs/m.cif",
                "x.tar.zst!-oremote.cif", "-oProxyCommand.cif", "~/secret"]:
        expect_raises(f"reject ref {bad!r}", lambda b=bad: nc.validate_ref(b))

    # ---- folder mechanics (temp LV_ROOT, stubbed extract) -----------------------------
    with tempfile.TemporaryDirectory() as td:
        nc.LV_ROOT = Path(td) / "LV_Coding"
        nc.nas_available = lambda: True
        nc._extract_member = fake_extract

        print("add_design / upsert:")
        r1 = nc.add_design(["setA"], "d_seq5", refs("d_seq5"),
                           {"design_id": "d_seq5", "stack": "RFD3", "sc_value": 0.6},
                           note="nice", dataset="100k", driven_metric="sc")
        check("add ok", r1[0]["ok"] and r1[0]["cifs"] == 3)
        csv = nc.LV_ROOT / "setA" / "metrics.csv"
        df = pd.read_csv(csv)
        check("1 row after first add", len(df) == 1)
        check("extra cols present", {"note", "added_utc", "dataset", "driven_metric"} <= set(df.columns))
        check("note stored", df.iloc[0]["note"] == "nice")
        cifdir = nc.LV_ROOT / "setA" / "cif"
        check("3 cif files", len(list(cifdir.glob("*.cif"))) == 3)

        # re-add same design_id -> still 1 row (idempotent), note updated
        nc.add_design(["setA"], "d_seq5", refs("d_seq5"),
                      {"design_id": "d_seq5", "stack": "RFD3", "sc_value": 0.7}, note="better")
        df = pd.read_csv(csv)
        check("still 1 row after re-add", len(df) == 1)
        check("note updated", df.iloc[0]["note"] == "better")
        check("value updated", abs(float(df.iloc[0]["sc_value"]) - 0.7) < 1e-9)

        # second design with an extra column -> 2 rows, column union
        nc.add_design(["setA"], "d_seq50", refs("d_seq50"),
                      {"design_id": "d_seq50", "stack": "BG", "sc_value": 0.5, "newcol": 42})
        df = pd.read_csv(csv)
        check("2 rows", len(df) == 2)
        check("column union added newcol", "newcol" in df.columns)

        print("folder_members / list_folders:")
        check("members", nc.folder_members("setA") == {"d_seq5", "d_seq50"})
        (nc.LV_ROOT / "not_a_collection").mkdir(parents=True, exist_ok=True)   # no metrics.csv
        folders = {f["name"]: f["count"] for f in nc.list_folders()}
        check("listed with count", folders.get("setA") == 2)
        check("non-collection dir hidden", "not_a_collection" not in folders)

        print("remove_design (prefix-safe):")
        rem = nc.remove_design("setA", "d_seq5")
        check("remove ok", rem["ok"])
        df = pd.read_csv(csv)
        check("1 row after remove", len(df) == 1 and df.iloc[0]["design_id"] == "d_seq50")
        names = sorted(p.name for p in cifdir.glob("*.cif"))
        check("d_seq5 cifs gone", not any(n.startswith("d_seq5_") for n in names))
        check("d_seq50 cifs KEPT (prefix not clobbered)",
              sum(n.startswith("d_seq50_") for n in names) == 3)

        print("browse_folder / serve_cif_path:")
        recs = nc.browse_folder("setA")
        check("browse returns the surviving design", len(recs) == 1 and recs[0]["design_id"] == "d_seq50")
        good = sorted(cifdir.glob("*.cif"))[0].name
        p = nc.serve_cif_path("setA", good)
        check("serve valid cif", p.is_file())
        expect_raises("serve rejects ../escape", lambda: nc.serve_cif_path("setA", "../../metrics.csv"))
        expect_raises("serve rejects absolute", lambda: nc.serve_cif_path("setA", "/etc/passwd"))

        # removing the last design empties the folder -> flagged + the dir is cleaned up
        rem2 = nc.remove_design("setA", "d_seq50")
        check("emptied flag set on last removal", rem2.get("emptied") is True)
        check("metrics.csv removed when empty", not csv.exists())
        check("empty folder dir cleaned up", not (nc.LV_ROOT / "setA").exists())

        print("corrupt metrics.csv does NOT silently wipe the folder:")
        nc.add_design(["setB"], "keepA", refs("keepA"), {"design_id": "keepA", "sc_value": 1})
        nc.add_design(["setB"], "keepB", refs("keepB"), {"design_id": "keepB", "sc_value": 2})
        bcsv = nc.LV_ROOT / "setB" / "metrics.csv"
        bcsv.write_bytes(b"design_id,sc_value\n\xff\xfe not valid utf-8 \x80\x81\n")  # -> UnicodeDecodeError
        corrupt_bytes = bcsv.read_bytes()
        res = nc.add_design(["setB"], "keepC", refs("keepC"), {"design_id": "keepC", "sc_value": 3})
        check("add reports failure (not silent success)", res[0]["ok"] is False)
        check("corrupt metrics.csv NOT overwritten", bcsv.read_bytes() == corrupt_bytes)

    print()
    if _FAILS:
        print(f"FAILED ({len(_FAILS)}): " + ", ".join(_FAILS))
        sys.exit(1)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
