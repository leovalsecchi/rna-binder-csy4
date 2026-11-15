#!/usr/bin/env python3
"""open_cofolds.py — open ALL scaffold co-folds in ChimeraX with every binder aligned on the
shared scaffold RNA (chain B), so you see at a glance which binders dock a hairpin vs drift onto
the spacer. One reference scaffold (grey) + every binder around it, colored by whether it binds a
hairpin (green) or drifts (red) — read from the scaffold_cofold_metric sidecar.

    python3 ~/Documents/rna_binder_gui/open_cofolds.py \
        [--outputs <euler AF3 outputs dir>] [--metric <sidecar.csv on euler>] \
        [--limit N] [--only binders|drifters]

Defaults to the pilot (analyses/scaffold_fit/pilot/outputs + its sidecar). Needs a working
`ssh euler` and ChimeraX in /Applications; CIFs stay on Euler, only the opened ones are streamed.
"""
from __future__ import annotations
import argparse
import csv
import os
import shlex
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import opendesign as OD  # reuse its ChimeraX launch / REST helpers

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SSH = ["ssh", "euler"]
RNA_CHAIN = "B"                         # scaffold RNA is chain B; binder is chain A (AF3 input ids)
PILOT_OUT = f"{RB}/analyses/scaffold_fit/pilot/outputs"
PILOT_METRIC = f"{RB}/analyses/scaffold_fit/pilot/scaffold_cofold_pilot.csv"
# distinct fallback palette when there's no metric to color by
PALETTE = ["cornflower blue", "orange", "medium purple", "teal", "goldenrod", "hot pink",
           "sea green", "chocolate", "steel blue", "dark khaki"]


def _ssh_out(cmd):
    return subprocess.run(SSH + [cmd], capture_output=True, text=True)


def list_models(outdir):
    """Top-ranked AF3 model per design: <outdir>/<name>/<name>_model.cif (not the per-seed ones)."""
    r = _ssh_out(f"find {shlex.quote(outdir)} -name '*_model.cif' ! -path '*seed-*' 2>/dev/null")
    return sorted(l for l in r.stdout.splitlines() if l.strip())


def load_metric(path):
    """{design_id_lower: binds_bool} from the sidecar (empty if absent)."""
    r = _ssh_out(f"cat {shlex.quote(path)} 2>/dev/null")
    if r.returncode != 0 or not r.stdout.strip():
        return {}
    out = {}
    for row in csv.DictReader(r.stdout.splitlines()):
        did = (row.get("design_id") or "").strip().lower()
        b = (row.get("scaffold_binds_hairpin") or "").strip()
        if did and b != "":
            out[did] = (b == "1")
    return out


def _design_of(remote_cif):
    return os.path.basename(remote_cif)[:-len("_model.cif")].lower()


def _short(name):
    return name.replace("csy4free_", "").replace("_len80", "")[:26]


def extract(remote_cif):
    tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False, prefix="cofold_")
    tmp.close()
    with open(tmp.name, "wb") as fh:
        r = subprocess.run(SSH + [f"cat {shlex.quote(remote_cif)}"], stdout=fh, stderr=subprocess.PIPE)
    if r.returncode == 0 and os.path.getsize(tmp.name) > 0:
        return tmp.name
    try:
        os.unlink(tmp.name)
    except OSError:
        pass
    return None


def build_cxc(models, binds):
    """models = [(design_name, local_cif)] -> ChimeraX script aligning all binders on the scaffold."""
    L = ["close session", "set bgColor white"]
    for i, (name, path) in enumerate(models, 1):
        L.append(f'open "{path}"')
        L.append(f"rename #{i} {_short(name)}")
    for i in range(2, len(models) + 1):
        L.append(f"matchmaker #{i}/{RNA_CHAIN} to #1/{RNA_CHAIN}")   # superpose on the shared scaffold
    L += ["hide atoms", "show cartoon"]
    # keep ONE reference scaffold (model #1's RNA); hide every other copy so it's not N-fold clutter
    for i in range(2, len(models) + 1):
        L.append(f"hide #{i}/{RNA_CHAIN} cartoon")
    L.append(f"color #1/{RNA_CHAIN} gray")
    # color each binder (chain A): green = docks a hairpin, red = drifts, palette when no metric
    for i, (name, _p) in enumerate(models, 1):
        b = binds.get(name)
        if b is True:
            col = "forest green"
        elif b is False:
            col = "firebrick"
        else:
            col = PALETTE[(i - 1) % len(PALETTE)]
        L.append(f"color #{i}/A {col}")
    L += ["view", "set bgColor white"]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", default=PILOT_OUT)
    ap.add_argument("--metric", default=PILOT_METRIC)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", choices=["binders", "drifters"], default=None,
                    help="show only hairpin-binders or only drifters (needs the metric sidecar)")
    a = ap.parse_args()

    cifs = list_models(a.outputs)
    if not cifs:
        print(f"no *_model.cif found under {a.outputs}", file=sys.stderr); return
    binds = load_metric(a.metric) if a.metric else {}

    picked = []
    for c in cifs:
        name = _design_of(c)
        b = binds.get(name)
        if a.only == "binders" and b is not True:
            continue
        if a.only == "drifters" and b is not False:
            continue
        picked.append((name, c))
    if a.limit:
        picked = picked[:a.limit]
    if not picked:
        print("nothing to show after filtering", file=sys.stderr); return

    nb = sum(1 for n, _ in picked if binds.get(n) is True)
    nd = sum(1 for n, _ in picked if binds.get(n) is False)
    print(f"opening {len(picked)} co-folds (scaffold-aligned)"
          + (f" — {nb} bind a hairpin (green), {nd} drift (red)" if binds else ""))
    models = []
    for name, c in picked:
        p = extract(c)
        if p:
            models.append((name, p))
        else:
            print(f"  extract failed: {name}", file=sys.stderr)
    if not models:
        print("nothing extracted (check ssh euler)", file=sys.stderr); return

    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".cxc", delete=False, prefix="opencofolds_")
    tf.write(build_cxc(models, binds))
    tf.close()
    if OD._rest_run(f"open {tf.name}"):
        print(f"opened {len(models)} co-folds in ChimeraX (grey scaffold + binders on it)")
    else:
        subprocess.Popen([OD._chimerax_bin(), "--cmd",
                          f"remotecontrol rest start port {OD.REST_PORT}", tf.name])
        print(f"launched ChimeraX with {len(models)} co-folds")


if __name__ == "__main__":
    main()
