#!/usr/bin/env python3
"""build_library_map.py — self-contained fold-map explorer for the 12k library selection.

One dot per passing design, laid out on the Foldseek fine-run UMAP (design layer, so the
coordinates are the GENERATOR's backbone, not a refold). Clicking a dot posts that design's
whole fold cluster to the open_refolds bridge on :45001, clicked design first, so ChimeraX's
DesignCycler opens design+AF3+ESM and cycling walks the rest of the cluster.

    python3 build_library_map.py            # writes exports/library_map.html

Needs the bridge running:  python3 ~/Documents/rna_binder_gui/open_refolds_server.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

BF   = Path("/Users/leo/Documents/binder_filter")
EXP  = BF / "exports"
BM   = Path("/Users/leo/Documents/binder_map")
EMBED = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-leo/093c41c9-0439-47bd-8cd9-d7a2131eebcc/scratchpad/embed_fine.parquet")
RES  = "cluster_r80"
OUT  = EXP / "library_map.html"
JIT  = 0.055          # UMAP units; ~0.17% of the range. Makes same-backbone designs clickable.

def templatise(paths, uuids, ids):
    """Replace uuid8 and design_id with placeholders -> a handful of templates."""
    tpl, idx = {}, np.empty(len(paths), dtype=np.int16)
    for i, (p, u, d) in enumerate(zip(paths, uuids, ids)):
        t = str(p).replace(str(d), "<ID>").replace(str(u), "<U>")
        idx[i] = tpl.setdefault(t, len(tpl))
    return list(tpl), idx

def main() -> int:
    df = pd.read_csv(EXP / "passing_v1.csv", skiprows=2)
    mem   = pd.read_csv(BM / "local_results/design/dedup/members.tsv", sep="\t")
    nodes = pd.read_csv(BM / "local_results/design/dedup/nodes.tsv", sep="\t")
    rep, ids = [], []
    for sn, x in zip(mem["struct_name"], mem["design_ids"]):
        for i in str(x).split(","): rep.append(sn); ids.append(i)
    nid = pd.Series(nodes["node_id"].to_numpy(), index=nodes["struct_name"])
    df["node"] = pd.Series(nid.reindex(pd.Index(rep)).to_numpy(),
                           index=ids).reindex(df["design_id"]).to_numpy()
    cl = pd.read_parquet(BM / "local_results_fine/design/clusters/fold/clusters.parquet"
                         ).set_index("node_id")[RES]
    df["cluster"] = cl.reindex(df["node"]).to_numpy()
    dropped = int(df["cluster"].isna().sum())
    df = df.dropna(subset=["cluster"]).reset_index(drop=True)
    df["cluster"] = df["cluster"].astype(int); df["node"] = df["node"].astype(int)

    emb = pd.read_parquet(EMBED).set_index("node_id")
    df["x"] = emb["x"].reindex(df["node"]).to_numpy()
    df["y"] = emb["y"].reindex(df["node"]).to_numpy()
    rng = np.random.default_rng(0)                      # deterministic, seeded
    df["x"] += rng.normal(0, JIT, len(df)); df["y"] += rng.normal(0, JIT, len(df))

    for tag, fn in [("inA", "library_12k_base_iptm.csv"), ("inB", "library_12k_cluster_rr_r80.csv")]:
        p = EXP / fn
        df[tag] = df["design_id"].isin(set(pd.read_csv(p)["design_id"])).astype(np.uint8) \
                  if p.exists() else np.zeros(len(df), np.uint8)

    stacks = sorted(df["stack"].unique())
    scode  = {s: i for i, s in enumerate(stacks)}
    tpls, tidx = {}, {}
    for c in ("design_cif", "af3_cif", "esm_cif"):
        tpls[c], tidx[c] = templatise(df[c].to_numpy(), df["uuid8"].to_numpy(),
                                      df["design_id"].to_numpy())
    r3 = lambda s: [round(float(v), 3) for v in s]
    data = {
        "res": RES, "dropped": dropped, "n": len(df),
        "stacks": stacks,
        "ids": "\n".join(df["design_id"]), "uuids": "\n".join(df["uuid8"].astype(str)),
        "x": r3(df["x"]), "y": r3(df["y"]),
        "cluster": [int(v) for v in df["cluster"]],
        "stack": [scode[s] for s in df["stack"]],
        "iptm": r3(df["iptm_mean"]), "rog": r3(df["rog_af3"]),
        "dsasa": [round(float(v), 1) for v in df["dsasa_af3"]], "fbase": r3(df["fbase_af3"]),
        "inA": [int(v) for v in df["inA"]], "inB": [int(v) for v in df["inB"]],
        "tpl": {k: tpls[k] for k in tpls},
        "ti":  {k: [int(v) for v in tidx[k]] for k in tidx},
    }
    html = TEMPLATE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    OUT.write_text(html)
    print(f"{len(df):,} designs  {df['cluster'].nunique():,} clusters  "
          f"({dropped} dropped: no cluster)")
    print(f"templates: " + ", ".join(f"{k} {len(v)}" for k, v in tpls.items()))
    print(f"wrote {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)")
    return 0

TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>Library fold map</title>
<style>
 :root{--bg:#11131a;--fg:#e8eaf0;--dim:#8b93a7;--line:#262a36;--card:#171a23}
 *{box-sizing:border-box} html,body{margin:0;height:100%}
 body{background:var(--bg);color:var(--fg);font:13px/1.5 ui-sans-serif,-apple-system,Segoe UI,Roboto,sans-serif;display:flex}
 #main{flex:1;position:relative;min-width:0}
 canvas{display:block;width:100%;height:100%;cursor:crosshair}
 #side{width:330px;flex:none;background:var(--card);border-left:1px solid var(--line);
       padding:14px;overflow:auto}
 h1{font-size:15px;margin:0 0 2px} .sub{color:var(--dim);font-size:11px;margin-bottom:14px}
 fieldset{border:1px solid var(--line);border-radius:6px;padding:8px 10px;margin:0 0 12px}
 legend{color:var(--dim);font-size:10px;text-transform:uppercase;letter-spacing:.08em;padding:0 4px}
 label{display:block;margin:3px 0;cursor:pointer} input,select{vertical-align:middle}
 select{background:#0e1017;color:var(--fg);border:1px solid var(--line);border-radius:4px;
        padding:3px 6px;width:100%}
 .kv{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;font-size:12px}
 .kv b{color:var(--dim);font-weight:400}
 .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;word-break:break-all}
 #legend{display:flex;flex-wrap:wrap;gap:4px 10px;font-size:11px}
 .sw{width:9px;height:9px;border-radius:2px;display:inline-block;margin-right:4px}
 #status{font-size:11px;color:var(--dim);min-height:2.6em;margin-top:6px}
 #tip{position:absolute;pointer-events:none;background:#000c;border:1px solid var(--line);
      border-radius:4px;padding:4px 7px;font-size:11px;display:none;white-space:nowrap}
 button{background:#222736;color:var(--fg);border:1px solid var(--line);border-radius:5px;
        padding:5px 9px;cursor:pointer;font:inherit} button:hover{background:#2c3242}
 .hint{color:var(--dim);font-size:10.5px;margin-top:8px;line-height:1.45}
</style></head><body>
<div id="main"><canvas id="cv"></canvas><div id="tip"></div></div>
<div id="side">
 <h1>Library fold map</h1>
 <div class="sub" id="hdr"></div>
 <fieldset><legend>colour</legend><select id="colour">
   <option value="stack">generator stack</option>
   <option value="lib">library membership</option>
   <option value="iptm">ipTM</option>
   <option value="cluster">fold cluster</option>
 </select><div id="legend" style="margin-top:8px"></div></fieldset>
 <fieldset><legend>show</legend><select id="filter">
   <option value="all">all designs</option>
   <option value="A">base 12k only</option>
   <option value="B">cluster 12k only</option>
   <option value="both">in both</option>
   <option value="onlyB">in cluster-12k, not base</option>
   <option value="onlyA">in base-12k, not cluster</option>
 </select>
 <label style="margin-top:8px">point size <input type="range" id="psz" min="1" max="6" value="2"></label>
 </fieldset>
 <fieldset><legend>on click — open in ChimeraX</legend>
  <label>cluster-mates to load <select id="ncyc">
    <option>1</option><option>5</option><option selected>10</option><option>20</option>
  </select></label>
  <div class="hint">Posts the cluster to the bridge on :45001, clicked design first.
   ChimeraX opens DESIGN + AF3 + ESM; cycling steps through the rest of the cluster.
   Each design is 3 files streamed from Euler, so keep this modest.</div>
 </fieldset>
 <fieldset><legend>selected</legend><div class="kv" id="info">
   <b>—</b><span>click a dot</span></div><div id="status"></div></fieldset>
 <div class="hint">Coordinates are the Foldseek <b>design-layer</b> UMAP: the generator's own
  backbone, not a refold. Designs sharing a backbone sit at the same spot (tiny jitter added
  so each is clickable).</div>
</div>
<script>
const D = __DATA__;
const IDS = D.ids.split("\n"), UU = D.uuids.split("\n");
const N = D.n;
const PAL = ["#4a7ebb","#4f9d69","#e07b39","#8264a5","#b5514f","#46a5b3","#c2a53a"];
const cv=document.getElementById("cv"), ctx=cv.getContext("2d"), tip=document.getElementById("tip");
let W=0,H=0,DPR=Math.min(2,devicePixelRatio||1);
const xs=D.x, ys=D.y;
const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
let view={s:1,ox:0,oy:0};
function resize(){ W=cv.clientWidth; H=cv.clientHeight; cv.width=W*DPR; cv.height=H*DPR;
  ctx.setTransform(DPR,0,0,DPR,0,0); draw(); }
addEventListener("resize",resize);
const PADX=26,PADY=26;
function sx(x){ return PADX+(x-xmin)/(xmax-xmin)*(W-2*PADX)*view.s+view.ox; }
function sy(y){ return H-PADY-(y-ymin)/(ymax-ymin)*(H-2*PADY)*view.s+view.oy; }
function hue(c){ let h=(c*2654435761)%360; return "hsl("+h+",62%,58%)"; }
function colourOf(i){
  const m=document.getElementById("colour").value;
  if(m==="stack") return PAL[D.stack[i]%PAL.length];
  if(m==="cluster") return hue(D.cluster[i]);
  if(m==="iptm"){ const t=Math.max(0,Math.min(1,(D.iptm[i]-0.55)/0.4));
    return "hsl("+(220-200*t)+",70%,"+(38+26*t)+"%)"; }
  const a=D.inA[i], b=D.inB[i];
  return a&&b ? "#4f9d69" : b ? "#e07b39" : a ? "#4a7ebb" : "#3a3f4d";
}
function visible(i){
  const f=document.getElementById("filter").value, a=D.inA[i], b=D.inB[i];
  return f==="all" || (f==="A"&&a) || (f==="B"&&b) || (f==="both"&&a&&b)
      || (f==="onlyB"&&b&&!a) || (f==="onlyA"&&a&&!b);
}
function draw(){
  ctx.clearRect(0,0,W,H);
  const r=+document.getElementById("psz").value;
  ctx.globalAlpha = r<=2 ? 0.72 : 0.85;
  for(let i=0;i<N;i++){ if(!visible(i)) continue;
    const X=sx(xs[i]), Y=sy(ys[i]); if(X<-5||X>W+5||Y<-5||Y>H+5) continue;
    ctx.fillStyle=colourOf(i); ctx.beginPath(); ctx.arc(X,Y,r,0,6.2832); ctx.fill(); }
  ctx.globalAlpha=1;
  if(sel>=0){ const X=sx(xs[sel]),Y=sy(ys[sel]);
    ctx.strokeStyle="#fff"; ctx.lineWidth=2; ctx.beginPath(); ctx.arc(X,Y,7,0,6.2832); ctx.stroke(); }
}
function nearest(px,py){
  let best=-1,bd=1e9;
  for(let i=0;i<N;i++){ if(!visible(i)) continue;
    const dx=sx(xs[i])-px, dy=sy(ys[i])-py, d=dx*dx+dy*dy;
    if(d<bd){bd=d;best=i;} }
  return bd<=400?best:-1;
}
let sel=-1;
function path(kind,i){ return D.tpl[kind][D.ti[kind][i]].replace("<U>",UU[i]).replace(/<ID>/g,IDS[i]); }
function info(i){
  const st=D.stacks[D.stack[i]], c=D.cluster[i];
  const mates=[]; for(let j=0;j<N;j++) if(D.cluster[j]===c) mates.push(j);
  document.getElementById("info").innerHTML =
    "<b>design</b><span class='mono'>"+IDS[i]+"</span>"+
    "<b>stack</b><span>"+st+"</span>"+
    "<b>cluster</b><span>"+c+" &middot; "+mates.length+" designs</span>"+
    "<b>ipTM</b><span>"+D.iptm[i].toFixed(3)+"</span>"+
    "<b>Rg</b><span>"+D.rog[i].toFixed(2)+" &Aring;</span>"+
    "<b>&Delta;SASA</b><span>"+D.dsasa[i].toFixed(0)+" &Aring;&sup2;</span>"+
    "<b>base frac</b><span>"+D.fbase[i].toFixed(3)+"</span>"+
    "<b>library</b><span>"+((D.inA[i]?"base":"")+(D.inA[i]&&D.inB[i]?" + ":"")+
      (D.inB[i]?"cluster":"")||"neither")+"</span>";
  return mates;
}
function send(i,mates){
  const k=+document.getElementById("ncyc").value;
  mates.sort((a,b)=>D.iptm[b]-D.iptm[a]);
  const order=[i].concat(mates.filter(j=>j!==i)).slice(0,k);
  const body={label:"cluster "+D.cluster[i]+" ("+mates.length+" designs)",
    designs:order.map(j=>({id:IDS[j],design:path("design_cif",j),
                           af3:path("af3_cif",j),esm:path("esm_cif",j)}))};
  const s=document.getElementById("status");
  s.textContent="sending "+order.length+" design(s) to ChimeraX ...";
  fetch("http://127.0.0.1:45001/openset",{method:"POST",
        headers:{"Content-Type":"text/plain"},body:JSON.stringify(body)})
   .then(r=>r.text()).then(t=>{s.textContent=t;})
   .catch(()=>{s.innerHTML="bridge not reachable &mdash; run<br><span class='mono'>python3 ~/Documents/rna_binder_gui/open_refolds_server.py</span>";});
}
cv.addEventListener("click",e=>{ const r=cv.getBoundingClientRect();
  const i=nearest(e.clientX-r.left,e.clientY-r.top); if(i<0) return;
  sel=i; draw(); send(i,info(i)); });
cv.addEventListener("mousemove",e=>{ const r=cv.getBoundingClientRect();
  const i=nearest(e.clientX-r.left,e.clientY-r.top);
  if(i<0){tip.style.display="none";return;}
  tip.style.display="block"; tip.style.left=(e.clientX-r.left+12)+"px";
  tip.style.top=(e.clientY-r.top+12)+"px";
  tip.textContent=D.stacks[D.stack[i]]+" | cl "+D.cluster[i]+" | ipTM "+D.iptm[i].toFixed(3); });
cv.addEventListener("wheel",e=>{ e.preventDefault();
  const f=e.deltaY<0?1.12:1/1.12, r=cv.getBoundingClientRect();
  const mx=e.clientX-r.left, my=e.clientY-r.top;
  view.ox=mx-(mx-view.ox)*f; view.oy=my-(my-view.oy)*f; view.s*=f; draw(); },{passive:false});
let drag=null;
cv.addEventListener("mousedown",e=>drag={x:e.clientX,y:e.clientY,ox:view.ox,oy:view.oy});
addEventListener("mouseup",()=>drag=null);
addEventListener("mousemove",e=>{ if(!drag)return;
  view.ox=drag.ox+e.clientX-drag.x; view.oy=drag.oy+e.clientY-drag.y; draw(); });
function legend(){
  const m=document.getElementById("colour").value, L=document.getElementById("legend");
  const row=(c,t)=>"<span><i class='sw' style='background:"+c+"'></i>"+t+"</span>";
  if(m==="stack") L.innerHTML=D.stacks.map((s,i)=>row(PAL[i%PAL.length],s)).join("");
  else if(m==="lib") L.innerHTML=row("#4f9d69","both")+row("#4a7ebb","base only")+
       row("#e07b39","cluster only")+row("#3a3f4d","neither");
  else if(m==="iptm") L.innerHTML=row("hsl(220,70%,38%)","0.55")+row("hsl(120,70%,51%)","0.75")+
       row("hsl(20,70%,64%)","0.95+");
  else L.innerHTML="<span style='color:var(--dim)'>"+D.res+" &mdash; one hue per cluster</span>";
}
for(const id of ["colour","filter","psz"])
  document.getElementById(id).addEventListener("input",()=>{legend();draw();});
document.getElementById("hdr").textContent =
  N.toLocaleString()+" designs · "+(new Set(D.cluster)).size+" clusters ("+D.res+
  ")"+(D.dropped?" · "+D.dropped+" unmapped dropped":"");
legend(); resize();
</script></body></html>
"""

if __name__ == "__main__":
    raise SystemExit(main())
