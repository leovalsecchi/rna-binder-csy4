#!/usr/bin/env python3
"""open_refolds_server.py — local bridge so clicking a dot in the HTML plots opens that
design's THREE refolds (AF3 / Boltz-2 / ESMFold2) in ChimeraX, superposed on the shared
16-nt RNA (so you see how each refolder poses the binder).

Run once on your Mac:

    python3 ~/Documents/rna_binder_gui/open_refolds_server.py

It (1) makes sure ChimeraX is running with its REST server on :45000, and (2) serves a
tiny HTTP endpoint on :45001 that the plot pages call when you click a point. On each
click it **ssh-extracts** the 3 refold CIFs from Euler (mount-INDEPENDENT — works even
when /Volumes/eulermount is down), opens them, and superposes on the RNA (chain B).

Nothing to install beyond a working `ssh euler` and ChimeraX in /Applications. The design
CIFs stay on Euler; only the clicked design's 3 files are streamed to a local temp dir.
Swap EXTRACT_MODE below to 'mount' to read from the sshfs mount instead.
"""
from __future__ import annotations
import gzip
import http.server
import os
import shlex
import shutil
import socketserver
import subprocess
import sys
import tempfile
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import opendesign as OD  # reuse its ChimeraX launch / REST helpers

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"     # Euler repo root (refs are relative to this)
SSH = ["ssh", "euler"]
BRIDGE_PORT = 45001
RNA_CHAIN = "B"                                       # fallback; real chain DETECTED per file
EXTRACT_MODE = "ssh"                                  # 'ssh' (on-demand, robust) or 'mount'
MOUNT_BASE = OD.MOUNT_BASE                            # used only when EXTRACT_MODE == 'mount'
LABELS = ("AF3", "B2", "ESM")
PROT_COLOR = {"AF3": "cornflower blue", "B2": "magenta", "ESM": "gold"}

_NT = {"A", "U", "G", "C"}


def _rna_chain_of(path: str) -> str | None:
    """Chain id holding the RNA in an mmCIF, or None if undetermined.

    Chain conventions are NOT uniform (verified on midrun_v9 2026-08-14):
      RFD3+LiMPNN                               -> A=RNA,     B=protein
      BG+BIF, BG+LiMPNN, BD1+LiMPNN, HD+LiMPNN  -> A=protein, B=RNA
      every refolder (AF3/B2/ESM)               -> A=protein, B=RNA
    A hardcoded "B" therefore anchors the refold RNA on the DESIGN'S BINDER for RFD3.
    """
    try:
        opener = gzip.open if path.endswith(".gz") else open
        hdr: list[str] = []
        chains: dict[str, set] = {}
        ci = ri = None
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("_atom_site."):
                    hdr.append(line.strip().split(".", 1)[1])
                    continue
                if hdr and line.startswith(("ATOM", "HETATM")):
                    if ci is None:
                        ci = (hdr.index("auth_asym_id") if "auth_asym_id" in hdr
                              else hdr.index("label_asym_id"))
                        ri = hdr.index("label_comp_id")
                    f = line.split()
                    if len(f) > max(ci, ri):
                        chains.setdefault(f[ci], set()).add(f[ri])
                elif chains:
                    break
        for cid, comps in chains.items():
            if comps and comps <= _NT:
                return cid
    except Exception:
        pass
    return None


def _extract_ref(ref: str) -> str | None:
    """ref = RB-relative '<archive.tar.zst>!<member>' (or a loose path). Stream the CIF to
    a local temp file. Returns the local path, or None on failure."""
    tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False, prefix="refold_")
    tmp.close()
    if EXTRACT_MODE == "mount":
        # read straight from the sshfs mount (fast when up)
        src = ref if ref.startswith("/") else os.path.join(MOUNT_BASE, ref)
        if "!" in src:
            arch, member = src.split("!", 1)
            cmd = f"zstd -dq -c {shlex.quote(arch)} | tar -xOf - {shlex.quote(member)} > {shlex.quote(tmp.name)}"
        else:
            cmd = f"cp {shlex.quote(src)} {shlex.quote(tmp.name)}"
        r = subprocess.run(["bash", "-c", cmd], stderr=subprocess.PIPE)
    else:
        # ssh-extract on Euler, stream the member back (mount-independent)
        if "!" in ref:
            arch, member = ref.split("!", 1)
            remote = f"cd {shlex.quote(RB)} && tar --zstd -xOf {shlex.quote(arch)} {shlex.quote(member)}"
        else:
            remote = f"cd {shlex.quote(RB)} && cat {shlex.quote(ref)}"
        with open(tmp.name, "wb") as fh:
            r = subprocess.run(SSH + [remote], stdout=fh, stderr=subprocess.PIPE)
    if r.returncode == 0 and os.path.getsize(tmp.name) > 0:
        return tmp.name
    try:
        os.unlink(tmp.name)
    except OSError:
        pass
    sys.stderr.write(f"  extract FAILED: {ref}\n  {r.stderr.decode()[:200]}\n")
    return None


def _build_cxc(models: list[tuple[str, str]], did: str = "") -> str:
    """models = [(label, local_cif)] -> ChimeraX script: open, superpose on the RNA, color.
    `did` (the design's original name/id) is shown as an on-screen 2D label."""
    L = ["close session", "set bgColor white"]
    for i, (label, path) in enumerate(models, 1):
        L.append(f'open "{path}"')
        L.append(f"rename #{i} {label}")
    # Align on the shared RNA. Reference (#1) and mobile chain letters are resolved
    # separately because the design and the refolds do not share a convention.
    # cutoffDistance none = no iterative pruning, matching dockq_self / dockq_LRMSD.
    rna = [_rna_chain_of(p) or RNA_CHAIN for _lab, p in models]
    for i in range(2, len(models) + 1):
        L.append(f"matchmaker #{i}/{rna[i - 1]} to #1/{rna[0]} cutoffDistance none")
    L += ["hide atoms", "show cartoon"]
    for i, (label, _p) in enumerate(models, 1):
        L.append(f"color #{i} & nucleic gray")                       # RNA = common frame
        L.append(f"color #{i} & protein {PROT_COLOR.get(label, 'slate gray')}")
    # show the design's original name in the corner (design_id ~= the on-disk file name)
    did = (did or "").strip().strip('"')
    if did and did != "?":
        L.append(f'2dlabels create designname text "{did}" xpos 0.02 ypos 0.965 size 14 color black')
    L.append("view")
    return "\n".join(L) + "\n"


def open_refolds(refs: dict[str, str], did: str = "") -> str:
    models = []
    for label in LABELS:
        ref = refs.get(label)
        if not ref:
            continue
        p = _extract_ref(ref)
        if p:
            models.append((label, p))
    if not models:
        return "no refolds resolved (check ssh euler + the paths)"
    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".cxc", delete=False, prefix="openrefolds_")
    tf.write(_build_cxc(models, did))
    tf.close()
    if OD._rest_run(f"open {tf.name}"):
        return f"opened {len(models)} refolds (RNA-aligned) in ChimeraX"
    subprocess.Popen([OD._chimerax_bin(),
                      "--cmd", f"remotecontrol rest start port {OD.REST_PORT}", tf.name])
    return f"launched ChimeraX with {len(models)} refolds"


DESIGN_LABELS = ("DESIGN", "AF3", "ESM")            # no B2 — binder_filter excludes it
CYCLE_ROOT = os.path.expanduser("~/Documents/binder_filter/cycle")
DESIGN_COLOR = {"DESIGN": "dim gray", "AF3": "cornflower blue", "ESM": "gold"}


def open_set(designs: list[dict], label: str = "") -> str:
    """Pre-extract a whole SET of designs into one folder, write the manifest DesignCycler
    reads, and open the first one. Extracting up front means stepping through the set costs
    no further ssh round-trips.

    designs = [{id, design, af3, esm}] with RB-relative refs. Model order puts the DESIGN
    first so each refold is superposed onto the structure that was asked for."""
    import concurrent.futures
    import json

    token = tempfile.mkdtemp(prefix="set_", dir=_ensure_cycle_root())
    jobs = []
    for d in designs:
        for lab in DESIGN_LABELS:
            ref = d.get(lab.lower())
            if ref:
                jobs.append((d["id"], lab, ref))

    got: dict[tuple[str, str], str] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(_extract_ref, ref): (did, lab) for did, lab, ref in jobs}
        for fut in concurrent.futures.as_completed(futs):
            did, lab = futs[fut]
            try:
                p = fut.result()
            except Exception:                                # noqa: BLE001
                p = None
            if not p:
                continue
            dest = os.path.join(token, f"{did}__{lab}.cif")
            try:
                os.replace(p, dest)                          # same volume: cheap rename
            except OSError:
                try:
                    shutil.move(p, dest)                     # different volume: copy + unlink
                except OSError:
                    continue
            got[(did, lab)] = dest

    man = {"label": label, "designs": []}
    for d in designs:
        files = [[lab, got[(d["id"], lab)]] for lab in DESIGN_LABELS if (d["id"], lab) in got]
        if files:
            man["designs"].append({"id": d["id"], "files": files})
    if not man["designs"]:
        return "no structures resolved (check ssh euler + the refs)"
    with open(os.path.join(token, "_designset.json"), "w") as f:
        json.dump(man, f)

    n_files = sum(len(x["files"]) for x in man["designs"])
    # two separate REST calls rather than one ';'-joined string: the preset must have finished
    # registering `designset` before it is invoked
    if OD._rest_run("preset DesignCycler") and OD._rest_run(f'designset "{token}"'):
        return (f"opened {len(man['designs'])} designs ({n_files} structures) — "
                f"step with nextd / backd")
    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".cxc", delete=False, prefix="openset_")
    tf.write(f'preset DesignCycler\ndesignset "{token}"\n')
    tf.close()
    subprocess.Popen([OD._chimerax_bin(),
                      "--cmd", f"remotecontrol rest start port {OD.REST_PORT}", tf.name])
    return f"launched ChimeraX with {len(man['designs'])} designs"


def _ensure_cycle_root() -> str:
    os.makedirs(CYCLE_ROOT, exist_ok=True)
    return CYCLE_ROOT


class Handler(http.server.BaseHTTPRequestHandler):
    def _headers(self, code=200, ctype="text/plain"):
        self.send_response(code)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Type", ctype)
        self.end_headers()

    def do_OPTIONS(self):
        # a page served from elsewhere (e.g. binder_filter tunnelled from Euler) may preflight
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Max-Age", "86400")
        self.end_headers()

    def do_POST(self):
        import json
        if urllib.parse.urlparse(self.path).path != "/openset":
            self._headers(404); self.wfile.write(b"not found"); return
        try:
            n = int(self.headers.get("Content-Length") or 0)
            body = json.loads(self.rfile.read(n) or b"{}")
            designs = body.get("designs") or []
            print(f"[openset] {len(designs)} designs  label={body.get('label','')!r}")
            msg = open_set(designs, body.get("label", ""))
        except Exception as e:                               # noqa: BLE001
            msg = f"error: {type(e).__name__}: {e}"
        self._headers(); self.wfile.write(msg.encode())

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        if u.path != "/open":
            self._headers(404); self.wfile.write(b"not found"); return
        q = urllib.parse.parse_qs(u.query)
        refs = {lab: q[lab.lower()][0] for lab in LABELS if lab.lower() in q}
        did = q.get("id", ["?"])[0]
        print(f"[open] {did}: {list(refs)}")
        try:
            msg = open_refolds(refs, did)
        except Exception as e:  # noqa: BLE001
            msg = f"error: {type(e).__name__}: {e}"
        self._headers(); self.wfile.write(msg.encode())

    def log_message(self, *a):  # quiet
        pass


def main():
    if not OD._rest_run("version"):
        print(f"Starting ChimeraX with REST on :{OD.REST_PORT} ...")
        subprocess.Popen([OD._chimerax_bin(), "--cmd",
                          f"remotecontrol rest start port {OD.REST_PORT}"])
    print(f"open_refolds bridge -> http://127.0.0.1:{BRIDGE_PORT}  "
          f"(mode={EXTRACT_MODE}; Ctrl-C to stop)")
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    with socketserver.ThreadingTCPServer(("127.0.0.1", BRIDGE_PORT), Handler) as httpd:
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nbye")


if __name__ == "__main__":
    main()
