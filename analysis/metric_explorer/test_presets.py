#!/usr/bin/env python3
"""Unit tests for presets — named filter presets stored in presets.json.

Run:  python3 test_presets.py      (plain asserts, no pytest needed)

PRESETS_FILE is pointed at a temp path, so the real presets.json is never touched.
"""
import sys
import json
import tempfile
from pathlib import Path

import presets as pr

_FAILS = []


def check(name, cond):
    print(("  ok  " if cond else "FAIL  ") + name)
    if not cond:
        _FAILS.append(name)


def expect_raises(name, fn, exc=ValueError):
    try:
        fn()
        check(name + " (should raise)", False)
    except exc:
        check(name, True)
    except Exception as e:                               # noqa: BLE001
        check(name + f" (wrong exc {type(e).__name__})", False)


def main():
    # ---- name validation -------------------------------------------------------------
    check("trims name", pr.sanitize_preset_name("  my set  ") == "my set")
    check("collapses inner whitespace", pr.sanitize_preset_name("a   b") == "a b")
    expect_raises("rejects empty", lambda: pr.sanitize_preset_name(""))
    expect_raises("rejects whitespace-only", lambda: pr.sanitize_preset_name("   "))
    expect_raises("rejects over-long", lambda: pr.sanitize_preset_name("x" * 200))

    with tempfile.TemporaryDirectory() as td:
        pr.PRESETS_FILE = Path(td) / "presets.json"

        # ---- empty state ------------------------------------------------------------
        check("empty when no file", pr.list_presets() == [])

        # ---- save -> list round-trip, raw filter text preserved ---------------------
        out = pr.save_preset("self-consistent", {"apo_rmsd": "<=2", "hbonds": ">=3"},
                             "apo_rmsd", 2.0, 1.0)
        check("save returns list of 1", len(out) == 1)
        got = pr.list_presets()
        check("listed after save", len(got) == 1 and got[0]["name"] == "self-consistent")
        p = got[0]
        check("filter text preserved verbatim (<=2 not 2)", p["filters"]["apo_rmsd"] == "<=2")
        check("second filter preserved", p["filters"]["hbonds"] == ">=3")
        check("driven persisted", p["driven"] == "apo_rmsd")
        check("thresh persisted", p["thresh"] == 2.0)
        check("delta persisted", p["delta"] == 1.0)
        check("saved_utc stamped", bool(p.get("saved_utc")))

        # ---- upsert by name: replaces, does not duplicate ---------------------------
        pr.save_preset("self-consistent", {"apo_rmsd": "<=1.5"}, "iptm", 0.6, 0.1)
        got = pr.list_presets()
        check("upsert keeps 1 record", len(got) == 1)
        check("upsert updated filters", got[0]["filters"] == {"apo_rmsd": "<=1.5"})
        check("upsert updated driven", got[0]["driven"] == "iptm")

        # ---- a second preset, sorted by name ---------------------------------------
        pr.save_preset("another", {"sc": ">0.5"}, "sc", 0.5, 0.05)
        got = pr.list_presets()
        check("two presets", len(got) == 2)
        check("sorted by name", [x["name"] for x in got] == ["another", "self-consistent"])

        # ---- zero-filter preset is legal (a 'clear all' preset) --------------------
        pr.save_preset("clear all", {}, None, None, None)
        got = {x["name"]: x for x in pr.list_presets()}
        check("empty-filters preset saved", got["clear all"]["filters"] == {})
        check("driven None persisted", got["clear all"]["driven"] is None)

        # ---- delete ----------------------------------------------------------------
        left = pr.delete_preset("another")
        names = [x["name"] for x in left]
        check("delete removed only that one", "another" not in names and len(names) == 2)
        check("delete is idempotent-safe", isinstance(pr.delete_preset("nope"), list))

        # ---- corrupt file: list is lenient, save ABORTS (no data loss) -------------
        pr.PRESETS_FILE.write_text("{not json at all")
        check("corrupt -> list_presets() == []", pr.list_presets() == [])
        expect_raises("corrupt -> save_preset raises (no wipe)",
                      lambda: pr.save_preset("x", {"a": "1"}, None, None, None), Exception)
        check("corrupt file left intact", pr.PRESETS_FILE.read_text() == "{not json at all")

        # ---- filters values coerced to str -----------------------------------------
        pr.PRESETS_FILE.unlink()
        pr.save_preset("coerce", {"iptm": 0.5}, None, None, None)
        check("numeric filter coerced to str", pr.list_presets()[0]["filters"]["iptm"] == "0.5")
        check("file is valid json list", isinstance(json.loads(pr.PRESETS_FILE.read_text()), list))

    print("\n" + ("ALL PASSED" if not _FAILS else f"{len(_FAILS)} FAILED: " + ", ".join(_FAILS)))
    return 1 if _FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
