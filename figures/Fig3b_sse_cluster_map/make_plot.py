#!/usr/bin/env python3
"""Helix content vs strand content, one mark per fold cluster per generator.

Reads cluster_sse.tsv (produced on Euler by agg_sse.py over all 551,516 backbone
representatives) and writes a self-contained interactive page.

Design notes:
  - Form: scatter. Two quantitative measures, identity carried by colour.
  - Colour: the thesis stack palette (config/palette.json), validated for CVD
    separation with the dataviz validator (all checks PASS, light surface).
  - Generator is encoded TWICE, colour + marker shape, so the panel survives
    greyscale and colour-vision deficiency (STYLE_NOTES.md section 4).
  - Area is proportional to backbone count, with a size key.
  - No chart title (STYLE_NOTES.md section 1); the finding goes in the caption.
"""
import csv, json, pathlib

HERE = pathlib.Path(__file__).resolve().parent
# Adopted into LV_ThesisFigures 2026-11-09: inputs come from inputs/, outputs go
# to output/<fig>/, like every other figure here. Falls back to the script's own
# directory so the scripts still run standalone if copied out again.
ROOT = HERE.parents[1]
IN = ROOT / "inputs" if (ROOT / "inputs" / "cluster_sse.tsv").exists() else HERE
OUT = ROOT / "output" / HERE.name
OUT.mkdir(parents=True, exist_ok=True)
COLOR = {"RFD3": "#a83902", "BG": "#a376e9", "BD1": "#bd8703", "HD": "#d35182"}
# One shape for everything: the four stack colours carry identity on their own
# (CVD-validated, worst adjacent pair dE 12.8 deutan / 20.2 normal vision).
SHAPE = {"RFD3": "circle", "BG": "circle", "BD1": "circle", "HD": "circle"}
LABEL = {"RFD3": "RFD3 + LigandMPNN", "BG": "BoltzGen", "BD1": "BoltzDesign + LigandMPNN",
         "HD": "HalluDesign + LigandMPNN"}
ORDER = ["RFD3", "BG", "BD1", "HD"]          # canonical stack order, never sorted by value

rows = list(csv.DictReader(open(IN / "cluster_sse.tsv"), delimiter="\t"))
for r in rows:
    for k in ("n_bb", "n_designs", "lib_bb", "lib_designs", "cluster"):
        r[k] = int(r[k])
    for k in ("helix", "strand", "sd_helix", "sd_strand"):
        r[k] = float(r[k])
rows.sort(key=lambda r: -r["n_bb"])          # big marks drawn first, small on top

meta = json.loads((IN / "cluster_sse_meta.json").read_text()) if (IN / "cluster_sse_meta.json").exists() else {}
payload = {"rows": rows, "color": COLOR, "shape": SHAPE, "label": LABEL,
           "order": ORDER, "meta": meta}

HEAD = """<title>Fold Cluster SSE Map</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>
  :root {
    --surface: #fcfcfb; --panel: #ffffff;
    --ink: #171b1e; --ink2: #4e5356; --ink3: #6d7276;
    --grid: #e3e9ec; --axis: #b3b8bc; --limit: #cfd5d9;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --surface: #14181b; --panel: #1b2024;
      --ink: #f0f6fa; --ink2: #cfd5d9; --ink3: #919699;
      --grid: #2c3237; --axis: #4e5356; --limit: #3a4147;
    }
  }
  :root[data-theme="dark"] {
    --surface: #14181b; --panel: #1b2024;
    --ink: #f0f6fa; --ink2: #cfd5d9; --ink3: #919699;
    --grid: #2c3237; --axis: #4e5356; --limit: #3a4147;
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--surface); color: var(--ink);
         font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
  .wrap { max-width: 1060px; margin: 0 auto; padding-block: 28px 56px;
          padding-left: max(16px, env(safe-area-inset-left, 0px));
          padding-right: max(16px, env(safe-area-inset-right, 0px)); }
  h1 { font-size: 19px; font-weight: 600; margin: 0 0 4px; letter-spacing: -0.01em; }
  .sub { color: var(--ink2); font-size: 13px; margin: 0 0 20px; max-width: 70ch; }
  .controls { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 14px; }
  button { font: inherit; font-size: 12px; padding: 5px 11px; border-radius: 999px;
           border: 1px solid var(--grid); background: var(--panel); color: var(--ink2); cursor: pointer; }
  button[aria-pressed="true"] { border-color: currentColor; color: var(--ink); font-weight: 600; }
  button .dot { display: inline-block; width: 9px; height: 9px; border-radius: 2px; margin-right: 6px;
                vertical-align: -1px; }
  .panel { background: var(--panel); border: 1px solid var(--grid); border-radius: 10px; padding: 8px; }
  svg { display: block; width: 100%; height: auto; overflow: visible; }
  .tick { fill: var(--ink3); font-size: 10px; }
  .axlab { fill: var(--ink2); font-size: 11px; font-weight: 500; }
  .mark { cursor: pointer; }
  .mark:hover { stroke: var(--ink); stroke-width: 2; }
  #tip { position: fixed; pointer-events: none; opacity: 0; transition: opacity .08s;
         background: var(--panel); color: var(--ink); border: 1px solid var(--grid);
         border-radius: 8px; padding: 9px 11px; font-size: 12px; line-height: 1.45;
         box-shadow: 0 6px 20px rgba(0,0,0,.14); max-width: 280px; z-index: 9; }
  #tip b { font-weight: 600; }
  #tip .k { color: var(--ink3); }
  table { border-collapse: collapse; width: 100%; font-size: 12px; margin-top: 18px; }
  th, td { text-align: right; padding: 5px 9px; border-bottom: 1px solid var(--grid); }
  th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) { text-align: left; }
  th { color: var(--ink3); font-weight: 600; }
  details { margin-top: 18px; }
  summary { cursor: pointer; color: var(--ink2); font-size: 13px; }
  .note { color: var(--ink3); font-size: 11.5px; margin-top: 14px; max-width: 80ch; }
  @media (max-width: 620px) { .wrap { padding-block: 18px 40px; } }
</style>
"""

BODY = """<div class="wrap">
  <h1>Fold clusters in secondary-structure space</h1>
  <p class="sub">One mark per fold cluster per backbone generator. Area is proportional to the
  number of backbones. Position is the mean P-SEA helix and strand content of that generator's
  own backbones in the cluster, so a mark describes what the generator built, not the cluster medoid.</p>

  <div class="controls" id="legend"></div>
  <div class="panel"><svg id="chart" viewBox="0 0 900 600" role="img"
       aria-label="Scatter of helix content against strand content for fold clusters"></svg></div>

  <details>
    <summary>Data table</summary>
    <table id="tbl"></table>
  </details>

  <p class="note" id="foot"></p>
</div>
<div id="tip" role="tooltip"></div>
<script>
const P = __PAYLOAD__;
const svg = document.getElementById('chart'), tip = document.getElementById('tip');
const W = 900, H = 600, M = {t: 18, r: 18, b: 52, l: 58};
const iw = W - M.l - M.r, ih = H - M.t - M.b;
const x = v => M.l + (v / 100) * iw;
const y = v => M.t + ih - (v / 100) * ih;
const maxBB = Math.max(...P.rows.map(r => r.n_bb));
const rad = n => 4 + 26 * Math.sqrt(n / maxBB);
const on = new Set(P.order);

const NS = 'http://www.w3.org/2000/svg';
const el = (t, a) => { const e = document.createElementNS(NS, t);
  for (const k in a) e.setAttribute(k, a[k]); return e; };

function shapePath(kind, cx, cy, r) {
  if (kind === 'circle') return el('circle', {cx, cy, r});
  if (kind === 'square') { const s = r * 1.77;
    return el('rect', {x: cx - s/2, y: cy - s/2, width: s, height: s, rx: Math.min(3, s/5)}); }
  if (kind === 'triangle') { const s = r * 2.1;
    return el('polygon', {points: `${cx},${cy-s*0.62} ${cx+s*0.56},${cy+s*0.42} ${cx-s*0.56},${cy+s*0.42}`}); }
  const s = r * 1.3;
  return el('polygon', {points: `${cx},${cy-s} ${cx+s},${cy} ${cx},${cy+s} ${cx-s},${cy}`});
}

function draw() {
  svg.textContent = '';
  const gGrid = el('g', {}); svg.appendChild(gGrid);
  for (let v = 0; v <= 100; v += 20) {
    gGrid.appendChild(el('line', {x1: M.l, x2: M.l + iw, y1: y(v), y2: y(v),
      stroke: 'var(--grid)', 'stroke-width': 1}));
    gGrid.appendChild(el('line', {x1: x(v), x2: x(v), y1: M.t, y2: M.t + ih,
      stroke: 'var(--grid)', 'stroke-width': 1}));
    const tx = el('text', {x: x(v), y: M.t + ih + 16, 'text-anchor': 'middle', class: 'tick'});
    tx.textContent = v; gGrid.appendChild(tx);
    const ty = el('text', {x: M.l - 8, y: y(v) + 3.5, 'text-anchor': 'end', class: 'tick'});
    ty.textContent = v; gGrid.appendChild(ty);
  }
  const xl = el('text', {x: M.l + iw / 2, y: H - 12, 'text-anchor': 'middle', class: 'axlab'});
  xl.textContent = 'Helix content (%)'; svg.appendChild(xl);
  const yl = el('text', {x: 14, y: M.t + ih / 2, 'text-anchor': 'middle', class: 'axlab',
    transform: `rotate(-90 14 ${M.t + ih / 2})`});
  yl.textContent = 'Strand content (%)'; svg.appendChild(yl);

  for (const r of P.rows) {
    if (!on.has(r.gen)) continue;
    const g = el('g', {class: 'mark'});
    const s = shapePath(P.shape[r.gen], x(r.helix), y(r.strand), rad(r.n_bb));
    s.setAttribute('fill', P.color[r.gen]);
    s.setAttribute('fill-opacity', 0.62);
    s.setAttribute('stroke', 'var(--panel)');     // 2px surface ring keeps overlaps readable
    s.setAttribute('stroke-width', 2);
    g.appendChild(s);
    g.addEventListener('pointerenter', e => show(e, r));
    g.addEventListener('pointermove', e => place(e));
    g.addEventListener('pointerleave', () => { tip.style.opacity = 0; });
    svg.appendChild(g);
  }
  sizeKey();
}

function sizeKey() {
  const g = el('g', {}); svg.appendChild(g);
  const vals = [Math.round(maxBB / 10), Math.round(maxBB / 2), maxBB];
  let cx = M.l + 18;
  const cy = M.t + 30;
  const t0 = el('text', {x: M.l + 6, y: M.t + 62, class: 'tick'});
  t0.textContent = 'backbones'; g.appendChild(t0);
  for (const v of vals) {
    const r = rad(v);
    const c = el('circle', {cx: cx + r, cy, r, fill: 'none', stroke: 'var(--axis)', 'stroke-width': 1});
    g.appendChild(c);
    const t = el('text', {x: cx + r, y: cy + r + 12, class: 'tick', 'text-anchor': 'middle'});
    t.textContent = v >= 1000 ? (v / 1000).toFixed(0) + 'k' : v;
    g.appendChild(t);
    cx += 2 * r + 16;
  }
}

function show(e, r) {
  const pct = r.lib_bb ? (100 * r.lib_bb / r.n_bb).toFixed(1) : '0';
  tip.innerHTML = `<b>${P.label[r.gen]}</b><br>cluster r05 #${r.cluster}<br>`
    + `<span class="k">helix</span> ${r.helix.toFixed(1)}% <span class="k">± ${r.sd_helix.toFixed(1)}</span><br>`
    + `<span class="k">strand</span> ${r.strand.toFixed(1)}% <span class="k">± ${r.sd_strand.toFixed(1)}</span><br>`
    + `<span class="k">coil</span> ${(100 - r.helix - r.strand).toFixed(1)}%<br>`
    + `<span class="k">backbones</span> ${r.n_bb.toLocaleString()}<br>`
    + `<span class="k">designs</span> ${r.n_designs.toLocaleString()}<br>`
    + `<span class="k">in 12k library</span> ${r.lib_bb.toLocaleString()} bb (${pct}%), `
    + `${r.lib_designs.toLocaleString()} designs`;
  tip.style.opacity = 1; place(e);
}
function place(e) {
  const b = tip.getBoundingClientRect();
  let l = e.clientX + 14, t = e.clientY + 14;
  if (l + b.width > innerWidth - 8) l = e.clientX - b.width - 14;
  if (t + b.height > innerHeight - 8) t = e.clientY - b.height - 14;
  tip.style.left = l + 'px'; tip.style.top = t + 'px';
}

const lg = document.getElementById('legend');
for (const g of P.order) {
  const n = P.rows.filter(r => r.gen === g).length;
  const b = document.createElement('button');
  b.setAttribute('aria-pressed', 'true');
  b.innerHTML = `<span class="dot" style="background:${P.color[g]}"></span>${P.label[g]} <span style="opacity:.6">(${n})</span>`;
  b.onclick = () => { on.has(g) ? on.delete(g) : on.add(g);
    b.setAttribute('aria-pressed', on.has(g)); draw(); };
  lg.appendChild(b);
}

const tbl = document.getElementById('tbl');
tbl.innerHTML = '<thead><tr><th>Generator</th><th>Cluster</th><th>Backbones</th><th>Designs</th>'
  + '<th>Helix %</th><th>Strand %</th><th>Coil %</th><th>Library bb</th></tr></thead><tbody>'
  + P.rows.map(r => `<tr><td>${P.label[r.gen]}</td><td>r05 #${r.cluster}</td>`
      + `<td>${r.n_bb.toLocaleString()}</td><td>${r.n_designs.toLocaleString()}</td>`
      + `<td>${r.helix.toFixed(1)}</td><td>${r.strand.toFixed(1)}</td>`
      + `<td>${(100 - r.helix - r.strand).toFixed(1)}</td>`
      + `<td>${r.lib_bb.toLocaleString()}</td></tr>`).join('') + '</tbody>';

document.getElementById('foot').textContent =
  `Secondary structure is P-SEA (biotite annotate_sse) on the C-alpha trace of every backbone `
  + `representative, the same assignment the campaign's sse_* columns use. Element counts depend `
  + `on the P-SEA thresholds, so absolute percentages shift a little under DSSP. `
  + `Clusters shown have at least ${P.meta.min_bb || 200} backbones within a generator; `
  + `BoltzGen covers BG+LigandMPNN and BG+BoltzIF, which share backbones. `
  + `Cluster ids are campaign r05 and are not comparable to the passing or library maps.`;

draw();
</script>
"""

head = HEAD
body = BODY.replace("__PAYLOAD__", json.dumps(payload))

# Standalone file for file:// viewing. The doctype matters: without it the page
# renders in quirks mode and the SVG comes out blank.
standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
              '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
              + head + "</head>\n<body>\n" + body + "</body>\n</html>\n")
(OUT / "sse_cluster_map.html").write_text(standalone)

# Artifact build: the publish skeleton supplies doctype/head/body, so this one
# carries only the title, meta, style and content.
(OUT / "sse_cluster_map.artifact.html").write_text(head + body)

for f in ("sse_cluster_map.html", "sse_cluster_map.artifact.html"):
    print(f"wrote {OUT/f} ({(OUT/f).stat().st_size/1e3:.0f} kB, {len(rows)} marks)")
