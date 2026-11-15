# Click-to-open: dots in the HTML plots → 3 refolds in ChimeraX (RNA-aligned)

Clicking a point in the interactive plots opens that design's **three refolds**
(AF3 / Boltz-2 / ESMFold2) in ChimeraX, superposed on the shared 16-nt RNA — so you
see how each refolder poses the binder against a common RNA frame.

## One-time setup on the Mac
Nothing to install beyond a working `ssh euler` and ChimeraX in `/Applications`.
(No `/Volumes/eulermount` needed — CIFs are pulled on demand over `ssh euler`.)

## Use it
1. Start the bridge (also starts ChimeraX with its REST server if needed):
   ```
   python3 ~/Documents/rna_binder_gui/open_refolds_server.py
   ```
   Leave it running. It prints: `open_refolds bridge -> http://127.0.0.1:45001 (mode=ssh)`.
2. Open a **clickable** plot in your browser, e.g.
   `~/Downloads/rna_binder_plots_v2v3_10k/html/plot_04_scatter_iptm_dockq.html`.
3. **Click a dot.** The bridge ssh-extracts that design's 3 refold CIFs from Euler,
   opens them, and runs `matchmaker #i/B to #1/B` (aligns on RNA chain B). Proteins are
   colored per refolder (AF3 blue / B2 magenta / ESM gold), RNA gray. Repeat clicks reuse
   the same ChimeraX window instantly.

## Which plots are clickable
- `plot_04_scatter_iptm_dockq` (ipTM × DockQ) — one clickable point per design×refolder.
  (Violins / funnels / heatmaps / tables have no per-design dots, so they're not clickable.)
  More scatters can be made clickable the same way: set `CLICK_JS = H.CLICK_JS` on the
  module and put `[design_id, af3_ref, b2_ref, esm_ref]` in each point's `customdata`.

## How it works / knobs
- The plot HTML embeds each point's 3 refold CIF refs (`<archive>.tar.zst!member`,
  relative to the Euler repo) in `customdata`; a small injected script `fetch`es them to
  the bridge on click (`no-cors`, so it works from `file://`).
- The bridge (`open_refolds_server.py`) streams each member back with
  `ssh euler "tar --zstd -xOf <archive> <member>"` → local temp `.cif` → ChimeraX REST :45000.
- To read from the sshfs mount instead of ssh (faster when the mount is up), set
  `EXTRACT_MODE = "mount"` at the top of `open_refolds_server.py`.
- RNA chain assumed to be `B` (`RNA_CHAIN`); alignment/colors are in `_build_cxc`.
