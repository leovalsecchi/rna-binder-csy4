"""store.py — load the columnar campaign store and reconstruct CIF paths.

The store is one .npy per column, written by build_filter_store.py on Euler. Numeric columns
load into RAM (~370 MB for 6.5M designs) because every threshold change re-scans them; the
two string columns stay memory-mapped since they are only touched when a design is clicked or
exported.

CIF paths are not stored. Each row carries a template index; the path is rebuilt by literal
substitution of design_id and uuid8 into the template. The builder asserted a byte-exact
round-trip for every row, so this is reconstruction, not inference.

BOLTZ-2 IS ABSENT BY CONSTRUCTION — the store has no b2 column, so no filter can reach one.
"""
from __future__ import annotations

import json
import os

import numpy as np

# BF_STORE lets a deployment point at a store it does not own — on Euler the 775 MB store
# already exists and is group-readable, so it is used in place rather than copied per user.
STORE_DIR = os.environ.get(
    "BF_STORE", os.path.join(os.path.dirname(os.path.abspath(__file__)), "store"))
# design_id is 342 MB and only touched on a click or an export, so it stays mapped.
# uuid8 is 85 MB and is read on EVERY filter request (run-diversity), where the access is
# scattered — on Lustre that meant thousands of separate page faults and a 12 s first request.
# Small enough to just hold.
MMAP_COLS = ("design_id",)
# BF_MMAP=1 memory-maps the numeric columns too. Slower on first touch, but resident size
# drops from ~370 MB to near zero and the pages are shared by every user on the node via the
# page cache — which is what you want when several people run this on shared hardware.
MMAP_ALL = os.environ.get("BF_MMAP", "").strip() not in ("", "0", "false", "no")


class Store:
    def __init__(self, path: str = STORE_DIR):
        self.path = path
        with open(os.path.join(path, "meta.json")) as f:
            self.meta = json.load(f)
        self.n = int(self.meta["n"])
        self.templates = [tuple(t) for t in self.meta["templates"]]
        self.path_cols = list(self.meta["path_cols"])
        self.vocabs = {k: list(v) for k, v in self.meta["vocabs"].items()}

        self._cols: dict[str, np.ndarray] = {}
        for fn in sorted(os.listdir(path)):
            if not fn.endswith(".npy"):
                continue
            name = fn[:-4]
            mm = "r" if (name in MMAP_COLS or MMAP_ALL) else None
            self._cols[name] = np.load(os.path.join(path, fn), mmap_mode=mm)
        for c in ("design_id", "uuid8", "tpl_idx"):
            if c not in self._cols:
                raise RuntimeError(f"store {path} is missing required column {c}")

        # 2-backend consensus, computed here so it can never be read from a B2-tainted column
        self._cols["iptm_min"] = np.fmin(self._cols["iptm_af3"], self._cols["iptm_esm"])
        self._cols["iptm_mean"] = (self._cols["iptm_af3"] + self._cols["iptm_esm"]) / 2.0

    # ---------------------------------------------------------------- columns

    def __contains__(self, name: str) -> bool:
        return name in self._cols

    def col(self, name: str) -> np.ndarray:
        try:
            return self._cols[name]
        except KeyError:
            raise KeyError(f"no column {name!r} in store; have {sorted(self._cols)}")

    def finite(self, name: str) -> np.ndarray:
        """Cached np.isfinite for a column. Every atom and every quadrant needs it, and the
        answer cannot change while the store is loaded, so computing it per request was
        costing more than the comparisons themselves."""
        key = "~fin~" + name
        if key not in self._cols:
            self._cols[key] = np.isfinite(self._cols[name])
        return self._cols[key]

    def has_data(self, name: str) -> bool:
        """False for placeholder columns that are entirely NaN (e.g. globularity before its
        backfill). The UI greys these out instead of silently failing every design."""
        return name in self._cols and bool(np.isfinite(self._cols[name]).any())

    def numeric_cols(self) -> list[str]:
        return sorted(k for k, v in self._cols.items()
                      if v.dtype.kind == "f" and not k.startswith("~"))

    # ------------------------------------------------------------------ rows

    def design_ids(self, idx) -> list[str]:
        return [b.decode() for b in np.asarray(self._cols["design_id"])[idx]]

    def paths(self, i: int) -> dict[str, str]:
        """RB-relative '<archive>.tar.zst!<member>' refs for one row, keyed by the master CSV
        column name (af3_msaoff_cif_path / esm_cif_path / designed_cif_path)."""
        did = self._cols["design_id"][i].decode()
        u = self._cols["uuid8"][i].decode()
        tpl = self.templates[int(self._cols["tpl_idx"][i])]
        return {c: t.replace("<ID>", did).replace("<U>", u)
                for c, t in zip(self.path_cols, tpl)}

    def label(self, vocab: str, code: int) -> str:
        v = self.vocabs.get(vocab, [])
        return v[code] if 0 <= code < len(v) else ""

    # --------------------------------------------------------------- display

    def subsample(self, size: int, seed: int = 12345, stratify: str | None = None):
        """Row indices for the scatter plots. Drawn once with a fixed seed so the point cloud
        stays put while thresholds move — dragging a threshold recolours, never redraws.
        Uniform by default (honest density); stratify='stack' gives rare stacks visible mass
        at the cost of distorting that density."""
        size = min(int(size), self.n)
        rng = np.random.default_rng(seed)
        if stratify and stratify in self._cols:
            codes = self._cols[stratify]
            groups = np.unique(codes)
            per = max(1, size // len(groups))
            out = []
            for g in groups:
                where = np.flatnonzero(codes == g)
                take = min(per, where.size)
                out.append(rng.choice(where, take, replace=False) if take < where.size
                           else where)
            return np.sort(np.concatenate(out))
        return np.sort(rng.choice(self.n, size, replace=False))
