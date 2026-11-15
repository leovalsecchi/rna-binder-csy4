#!/usr/bin/env python3
"""Open a design + its 3 refold variants superposed in ChimeraX.

The four variants of one design:
    design       <run>/designs/<DESIGN>.{cif,pdb,cif.gz}        (original prediction)
    af3_msaoff   <run>/refold_af3/msaoff/designs/<DESIGN>_AF3_top.cif
    af3_msaon    <run>/refold_af3/msaon/designs/<DESIGN>_AF3_top.cif
    b2           <run>/refold_b2/msaoff/designs/<DESIGN>_B2_sample0.cif
    esm          <run>/refold_esm/msaoff/designs/<DESIGN>_ESM_sample0.cif

All refolds are superposed onto the design's binder (ChimeraX `matchmaker`, which
anchors on the protein chain), and each variant's RNA is coloured differently so
the pose disagreement is obvious. Binder shown grey (it overlaps after the fit).

Usage:
    opendesign.py RUN DESIGN [--print] [--b2-sample N]
      RUN     a design-run dir (absolute, or relative to the mount base)
      DESIGN  a pattern matching the design in <run>/designs/ (e.g. design84,
              design0_seq0, or a full stem)
    --print   resolve + emit the ChimeraX script, do NOT launch (for testing)

Can also be imported: open_design(designed_file, launch=True) takes a designed
file path directly (used by the GUI's "open in ChimeraX" button).
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import re
import shlex
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

MOUNT_BASE = "/Volumes/eulermount/lvalsecchi/rna_binder"
CAMPAIGN = "csy4_crrna_free"        # default campaign for --short
REST_PORT = 45000                    # reuse a running ChimeraX (instant repeat opens)
STACK_ALIASES = {
    "rfd3": "rfd3__ligandmpnn",
    "bg": "boltzgen__ligandmpnn", "boltzgen": "boltzgen__ligandmpnn",
    "bif": "boltzgen__boltzif",
    "hd": "halludesign__ligandmpnn", "halludesign": "halludesign__ligandmpnn",
    "bd1": "boltzdesign__ligandmpnn", "boltzdesign": "boltzdesign__ligandmpnn",
}
_DESIGN_EXTS = (".cif", ".pdb", ".cif.gz")


def _strip_ext(name: str) -> str:
    for e in (".cif.gz", ".cif", ".pdb"):
        if name.endswith(e):
            return name[: -len(e)]
    return name


def _slice_uuid(stem: str) -> str | None:
    """The run uuid8 embedded in a canonical design stem, if any."""
    m = re.search(r"_([0-9a-f]{8})(?:_|$)", stem)
    return m.group(1) if m else None


def _extract_member(archive: str, member: str) -> str | None:
    """Extract one member from a .tar.zst to a temp .cif (local zstd). None on miss."""
    tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False,
                                      prefix=_strip_ext(os.path.basename(member)) + "_")
    tmp.close()
    cmd = (f"zstd -dq -c {shlex.quote(archive)} | "
           f"tar -xOf - {shlex.quote(member)} > {shlex.quote(tmp.name)}")
    r = subprocess.run(["bash", "-c", cmd], stderr=subprocess.PIPE)
    if r.returncode == 0 and os.path.getsize(tmp.name) > 0:
        return tmp.name
    os.unlink(tmp.name)
    if b"command not found" in r.stderr and b"zstd" in r.stderr:
        sys.exit("opendesign needs `zstd` to read packed runs — install it: brew install zstd")
    return None


def _find_design_stem(run: Path, pattern: str) -> str:
    """Unique design stem matching `pattern`: loose files, per-uuid slice .lists,
    or a whole-dir designs.tar.zst.list."""
    designs = run / "designs"
    stems: set[str] = set()
    if designs.is_dir():
        for ext in _DESIGN_EXTS:                                   # loose
            for h in glob.glob(str(designs / f"*{pattern}*{ext}")):
                stems.add(_strip_ext(os.path.basename(h)))
        for lst in designs.glob("*.tar.zst.list"):                 # per-uuid slices
            for line in lst.read_text().splitlines():
                b = os.path.basename(line.strip())
                if pattern in b and b.endswith(_DESIGN_EXTS):
                    stems.add(_strip_ext(b))
    wd = run / "designs.tar.zst.list"                              # whole-dir
    if wd.exists():
        for line in wd.read_text().splitlines():
            b = os.path.basename(line.strip())
            if pattern in b and b.endswith(_DESIGN_EXTS):
                stems.add(_strip_ext(b))
    if not stems:
        sys.exit(f"no design matching '{pattern}' in {run}/designs")
    if len(stems) > 1:
        sys.exit("pattern matches multiple designs — be more specific:\n  " +
                 "\n  ".join(sorted(stems)))
    return next(iter(stems))


def _resolve_design(run: Path, stem: str) -> str | None:
    """Real local .cif for source design <stem>: loose, per-uuid slice, or whole-dir."""
    designs = run / "designs"
    for ext in _DESIGN_EXTS:                                       # loose
        p = designs / f"{stem}{ext}"
        if p.exists():
            return _prep_path(str(p))
    uu = _slice_uuid(stem)
    if uu and designs.is_dir():                                    # slice (flat members)
        arch = designs / f"{uu}.tar.zst"
        if arch.exists():
            for ext in _DESIGN_EXTS:
                got = _extract_member(str(arch), f"{stem}{ext}")
                if got:
                    return got
    wd = run / "designs.tar.zst"                                   # whole-dir (designs/ member)
    if wd.exists():
        for ext in _DESIGN_EXTS:
            got = _extract_member(str(wd), f"designs/{stem}{ext}")
            if got:
                return got
    return None


def _af3_from_outputs(base: Path, stem: str) -> str | None:
    """Fallback for UN-canonicalized AF3: pull the top model straight from the raw
    outputs.tar.zst. Member is `outputs/<STEM>/<stem.lower()>_<ts>/<stem.lower()>_model.cif`
    (the inner path is lowercased, unlike the design stem)."""
    arch = base / "outputs.tar.zst"
    lst = base / "outputs.tar.zst.list"
    if not (arch.exists() and lst.exists()):
        return None
    target = f"{stem.lower()}_model.cif"          # top model has the stem prefix; seed samples are bare model.cif
    for line in lst.read_text().splitlines():
        m = line.strip()
        if m.endswith(target):
            return _extract_member(str(arch), m)
    return None


def _resolve_refold(run: Path, stem: str, refolder: str, msa: str,
                    b2_sample: int) -> str | None:
    """Real local .cif for a refold variant; handles multi-run msaX_<uuid> leaves + packing.
    For AF3, falls back to the raw outputs.tar.zst when designs/ isn't canonicalized."""
    if refolder == "af3":
        member = f"{stem}_AF3_top.cif"
    elif refolder == "esm":
        member = f"{stem}_ESM_sample0.cif"
    else:
        member = f"{stem}_B2_sample{b2_sample}.cif"
    uu = _slice_uuid(stem)
    leaves = [f"{msa}_{uu}", msa] if uu else [msa]                 # multi-run leaf first
    for leaf in leaves:
        base = run / f"refold_{refolder}" / leaf
        loose = base / "designs" / member
        if loose.exists():
            return str(loose)
        arch = base / "designs.tar.zst"                            # whole-dir refold output
        if arch.exists():
            got = _extract_member(str(arch), f"designs/{member}")
            if got:
                return got
        if refolder == "af3":                                      # fallback: raw, un-canonicalized outputs
            got = _af3_from_outputs(base, stem)
            if got:
                return got
    return None


# (label, colour) for the design + each refold variant, in display order. Source
# of truth for make_html_plots' click-to-open path list (OD._VARIANTS).
_VARIANTS = [
    ("design", "forest"),
    ("af3_msaoff", "cornflower blue"),
    ("af3_msaon", "orange"),
    ("b2", "magenta"),
    ("esm", "gold"),
]


def resolve_variants(run: Path, stem: str, b2_sample: int = 0) -> list[tuple[str, str, str]]:
    """Return [(label, colour, local_path)] for the design + the refold variants present."""
    out = []
    d = _resolve_design(run, stem)
    if d:
        out.append(("design", "forest", d))
    else:
        print(f"  (missing design: {stem})", file=sys.stderr)
    for label, colour, refolder, msa in (
        ("af3_msaoff", "cornflower blue", "af3", "msaoff"),
        ("af3_msaon", "orange", "af3", "msaon"),
        ("b2", "magenta", "b2", "msaoff"),
        ("esm", "gold", "esm", "msaoff"),
    ):
        p = _resolve_refold(run, stem, refolder, msa, b2_sample)
        if p:
            out.append((label, colour, p))
        else:
            print(f"  (missing {label})", file=sys.stderr)
    return out


def _prep_path(path: str) -> str:
    """ChimeraX-openable path; decompress .cif.gz to a temp .cif."""
    if path.endswith(".cif.gz"):
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False,
                                          prefix=_strip_ext(os.path.basename(path)) + "_")
        with gzip.open(path, "rt") as f:
            tmp.write(f.read())
        tmp.close()
        return tmp.name
    return path


def build_cxc(variants: list[tuple[str, str, str]]) -> str:
    """ChimeraX command script: open, name, superpose on binder, colour RNA."""
    lines = ["set bgColor white"]
    for i, (label, _colour, path) in enumerate(variants, start=1):
        lines.append(f'open "{_prep_path(path)}"')
        lines.append(f"rename #{i} {label}")
    n = len(variants)
    if n > 1:
        # Superpose models 2..n onto #1 (matchmaker anchors on the protein binder).
        lines.append(f"matchmaker #2-{n} to #1")
    lines.append("hide atoms")
    lines.append("show cartoon")
    lines.append(f"color #1-{n} & protein gray")
    for i, (_label, colour, _path) in enumerate(variants, start=1):
        lines.append(f"color #{i} & nucleic {colour}")
    lines.append("view")
    return "\n".join(lines) + "\n"


def _chimerax_bin() -> str:
    apps = sorted(glob.glob("/Applications/ChimeraX*.app"))
    if not apps:
        sys.exit("ChimeraX.app not found in /Applications")
    return f"{apps[-1]}/Contents/MacOS/ChimeraX"


def _rest_run(command: str, port: int = REST_PORT, timeout: float = 2.0) -> bool:
    """Send a command to an already-running ChimeraX REST server; True if it answered."""
    url = f"http://127.0.0.1:{port}/run?command=" + urllib.parse.quote(command)
    try:
        urllib.request.urlopen(url, timeout=timeout).read()
        return True
    except Exception:
        return False


def open_design(designed_file: str, *, b2_sample: int = 0, launch: bool = True,
                print_only: bool = False) -> str:
    # designed_file is a logical <run>/designs/<stem>.cif path (it may not exist on
    # disk when the run is packed); we only use it to recover run + stem, then
    # resolve every variant archive-aware.
    dp = Path(designed_file)
    variants = resolve_variants(dp.parent.parent, _strip_ext(dp.name), b2_sample=b2_sample)
    if not variants:
        sys.exit("no variants found to open")
    cxc = "close session\n" + build_cxc(variants)   # replace prior models on reuse
    print("Opening variants:")
    for label, _c, p in variants:
        print(f"  {label:12s} {p}")
    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".cxc", delete=False, prefix="opendesign_")
    tf.write(cxc)
    tf.close()
    if print_only:
        print("\n--- ChimeraX script ---")
        print(cxc)
        print(f"(script: {tf.name})")
    elif launch:
        if _rest_run(f"open {tf.name}"):             # a ChimeraX is already open -> instant
            print(f"\nSent to running ChimeraX ({len(variants)} models) — instant.")
        else:                                         # first open: start ChimeraX + its REST server
            subprocess.Popen([_chimerax_bin(),
                              "--cmd", f"remotecontrol rest start port {REST_PORT}",
                              tf.name])
            print(f"\nLaunched ChimeraX ({len(variants)} models). "
                  "Leave the window open — later opens reuse it instantly.")
    return tf.name


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--short", action="store_true",
                    help="shorthand: STACK RUN UUID NUM (e.g. rfd3 midrun_v2 f2653389 84); _seq0 if no seq given")
    ap.add_argument("--campaign", default=CAMPAIGN, help=f"campaign for --short (default {CAMPAIGN})")
    ap.add_argument("--b2-sample", type=int, default=0, help="B2 sample index (default 0)")
    ap.add_argument("--print", dest="print_only", action="store_true",
                    help="resolve + print the ChimeraX script; do not launch")
    ap.add_argument("pos", nargs="+", help="full: RUN DESIGN   |   --short: STACK RUN UUID NUM")
    args = ap.parse_args(argv)

    if args.short:
        if len(args.pos) == 4:                       # parallel-run pool: scope to a uuid
            stack, run_name, uuid, num = args.pos
        elif len(args.pos) == 3:                     # back-compat: single-run dir
            stack, run_name, num = args.pos; uuid = ""
        else:
            sys.exit("--short needs: STACK RUN UUID NUM  (e.g. rfd3 midrun_v2 f2653389 84)\n"
                     "  (or STACK RUN NUM for an old single-run dir)")
        stack_dir = STACK_ALIASES.get(stack, stack)
        run = Path(MOUNT_BASE) / "campaigns" / args.campaign / stack_dir / run_name
        d = num if str(num).startswith("design") else f"design{num}"
        if "_seq" not in d:                          # default to seq0 when seq not given
            d = f"{d}_seq0"
        design = f"{uuid}_{d}" if uuid else d
    else:
        if len(args.pos) != 2:
            sys.exit("need: RUN DESIGN   (or use --short STACK RUN NUM)")
        run_name, design = args.pos
        run = Path(run_name)
        if not run.is_absolute() and not run.exists():
            run = Path(MOUNT_BASE) / run_name
    try:
        run_ok = run.is_dir()
    except OSError:                                  # stale/dropped sshfs mount
        sys.exit(f"can't reach {run}\n  (eulermount may be down — check: ls /Volumes/eulermount)")
    if not run_ok:
        sys.exit(f"run dir not found: {run}")
    stem = _find_design_stem(run, design)
    open_design(str(run / "designs" / f"{stem}.cif"),
                b2_sample=args.b2_sample, print_only=args.print_only)


if __name__ == "__main__":
    sys.exit(main())
