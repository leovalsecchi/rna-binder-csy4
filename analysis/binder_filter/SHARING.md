# Sharing binder_filter

The tool runs **on Euler** and you view it in your own browser through an ssh tunnel. Nothing
is copied, nothing leaves ETH, and everyone sees the same live store.

Anyone in `bsse-wicky-euler-storage` can run it — the store and the app are group-readable.

---

## For the person you're sharing with — the whole thing

**1. Log in to Euler and start it:**

```bash
ssh euler
/cluster/work/wicky/lvalsecchi/rna_binder/tools/binder_filter/bf-serve
```

It prints a command. It looks like this, with *your* username and a port that changes each time:

```
ssh -o StrictHostKeyChecking=accept-new -J euler -L 8010:127.0.0.1:44285 <you>@eu-a2p-352
```

**2. Paste that command into a terminal on your OWN machine.** Leave it open — that terminal
*is* the tunnel.

**3. Open <http://127.0.0.1:8010>.**

That's it. Nothing to install.

When you're done: `bf-serve stop` on Euler (it also stops by itself after 4 hours).
`bf-serve status` shows where it's running.

---

## What you get

The whole campaign — **6,821,957 designs, v2–v9** — filtered live, ~150 ms per change.
Two plots (sc_RMSD and sc_DockQ, AF3 × ESMFold2), a grouped AND/OR filter tree, an ipTM
threshold-or-rank gate, and a breakdown of what survives versus the campaign baseline.
CSV export writes to your own Euler home.

Read `README.md` for what the metrics mean and where the data gaps are — in particular that
**sc_RMSD is pose RMSD, a DockQ input**, so the two plots are not independent evidence.

---

## Optional: clicking a design open in ChimeraX

Works, but needs two things **on your own machine**: ChimeraX installed, and `ssh euler`
working without a password (the structures are streamed from Euler on demand).

```bash
python3 ~/Documents/rna_binder_gui/open_refolds_server.py     # or: openrefolds
```

Then shift-click a point to open one design, or shift-click empty space to open 20 random
designs from that region and step through them with `nextd` / `backd`.

The app can't do this for you from Euler — the compute node has no idea your ChimeraX exists.
So in remote mode the page talks to the bridge on your own machine directly. Without the
bridge running, everything else still works; clicking just tells you it isn't there.

---

## Notes

- **Each person runs their own instance.** No shared port, no dependency on someone else's job
  staying alive, and the app only ever listens on the compute node's loopback.
- **The port changes every launch** — it's picked free at start, so two people on the same node
  can't collide. Re-read the printed command each time.
- **CIF download doesn't work in this mode.** It needs the local ssh-extraction helper from
  `metric_explorer`, which isn't on Euler. CSV export is unaffected.
- The store is a snapshot: v2-v8 built 2026-08-07, v9 appended 2026-08-17.
  It does **not** include v9 (~1.19M designs). Rebuilding is
  `sbatch ops/submit/build_filter_store.sbatch` with a different `--csv`.
