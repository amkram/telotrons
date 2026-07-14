#!/usr/bin/env python3
"""Build a self-contained interactive HTML reviewer for telotron-ortholog locus_text files.

Walks   locus_text/<version>/<category>/*.txt
Emits   a single HTML file that renders each locus with IUPAC nt/aa coloring +
        conservation shading (the same scheme as the VSCode extension), lets you
        confirm / reject each locus with the keyboard, persists decisions in the
        browser, and exports the confirmed set to good_orthologs/<version>/<category>/
        either as a copy-script (portable) or directly via the File System Access API.

Usage:
    python scripts/build_ortholog_review.py \
        [--root work/results/telotron_orthologs/locus_text] \
        [--out  work/results/telotron_orthologs/review.html]
"""
import argparse, hashlib, json, os, sys
from pathlib import Path

CATEGORY_ORDER = ["telotron_present", "intron_present_nontelo", "intron_absent", "uncertain"]


def strip_unaligned(text: str):
    """Drop every '### UNALIGNED — … ###' block (header line through to the next
    '### …' section header), keeping the ALIGNED sections and the locus header.
    Returns (stripped_text, n_blocks_dropped)."""
    out, skip, n_blocks = [], False, 0
    for ln in text.split("\n"):
        if ln.startswith("### "):
            was_skipping = skip
            skip = "UNALIGNED" in ln
            if skip and not was_skipping:
                n_blocks += 1
        if not skip:
            out.append(ln)
    return "\n".join(out), n_blocks


def collect(root: Path):
    recs = []
    for version_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        version = version_dir.name
        for cat_dir in sorted(p for p in version_dir.iterdir() if p.is_dir()):
            category = cat_dir.name
            for f in sorted(cat_dir.glob("*.txt")):
                stripped, n_blocks = strip_unaligned(f.read_text(encoding="utf-8", errors="replace"))
                recs.append({
                    "id": f"{version}/{category}/{f.name}",
                    "version": version,
                    "category": category,
                    "file": f.name,
                    "src": str(f.resolve()),
                    "text": stripped,
                    "n_unaligned_dropped": n_blocks,
                })
    return recs


HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Telotron ortholog review</title>
<style>
  :root{
    --bg:#1e2127; --bg2:#23262d; --bg3:#2b2f38; --fg:#c9d1d9; --mut:#8b949e;
    --line:#363b45; --accent:#4f86c6; --ok:#3fb950; --no:#e0625b; --warn:#d4a72c;
  }
  *{box-sizing:border-box}
  html,body{margin:0;height:100%;background:var(--bg);color:var(--fg);
    font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:13px}
  #app{display:grid;grid-template-columns:300px 1fr;grid-template-rows:auto 1fr;height:100vh}
  header{grid-column:1/3;display:flex;align-items:center;gap:14px;flex-wrap:wrap;
    padding:8px 14px;background:var(--bg2);border-bottom:1px solid var(--line)}
  header h1{font-size:14px;margin:0;font-weight:600}
  .stat{color:var(--mut)} .stat b{color:var(--fg)}
  .pill{padding:2px 8px;border-radius:10px;background:var(--bg3);font-size:12px}
  .pill.ok{color:var(--ok)} .pill.no{color:var(--no)} .pill.un{color:var(--mut)}
  select,button,input{background:var(--bg3);color:var(--fg);border:1px solid var(--line);
    border-radius:6px;padding:4px 8px;font-size:12px;font-family:inherit}
  button{cursor:pointer} button:hover{border-color:var(--accent)}
  button.primary{background:var(--accent);border-color:var(--accent);color:#fff}
  button.ok{background:var(--ok);border-color:var(--ok);color:#06210f}
  button.no{background:var(--no);border-color:var(--no);color:#2a0a08}
  .spacer{flex:1}
  .bar{grid-column:1/3;height:3px;background:var(--bg3)}
  .bar>i{display:block;height:100%;background:var(--accent);width:0}
  aside{background:var(--bg2);border-right:1px solid var(--line);overflow:auto}
  .grp{padding:6px 10px 2px;color:var(--mut);font-size:11px;text-transform:uppercase;
    letter-spacing:.04em;position:sticky;top:0;background:var(--bg2)}
  .row{display:flex;align-items:center;gap:8px;padding:4px 10px;cursor:pointer;
    border-left:3px solid transparent;white-space:nowrap;overflow:hidden}
  .row:hover{background:var(--bg3)}
  .row.cur{background:var(--bg3);border-left-color:var(--accent)}
  .dot{width:8px;height:8px;border-radius:50%;background:#4b525d;flex:none}
  .dot.ok{background:var(--ok)} .dot.no{background:var(--no)}
  .row .nm{overflow:hidden;text-overflow:ellipsis;font-size:11px;color:var(--fg)}
  main{overflow:auto;padding:0}
  .meta{padding:8px 16px;border-bottom:1px solid var(--line);background:var(--bg2);
    position:sticky;top:0;z-index:2;display:flex;align-items:center;gap:12px;flex-wrap:wrap}
  .meta .fn{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px;color:var(--fg)}
  .badge{padding:2px 8px;border-radius:6px;font-size:11px;background:var(--bg3);color:var(--mut)}
  .badge.telotron_present{color:#d2b3ff;background:#352c47}
  .badge.intron_present_nontelo{color:#9fd0ff;background:#23344a}
  .badge.intron_absent{color:#ffb3a8;background:#46282a}
  .badge.uncertain{color:#e8d18a;background:#403821}
  pre{margin:0;padding:14px 16px 80px;font-family:ui-monospace,Menlo,Consolas,monospace;
    font-size:12.5px;line-height:1.42;white-space:pre;tab-size:4;color:var(--fg)}
  .tag{opacity:.32}
  footer{position:fixed;left:300px;right:0;bottom:0;display:flex;align-items:center;gap:10px;
    padding:8px 16px;background:var(--bg2);border-top:1px solid var(--line)}
  kbd{background:var(--bg3);border:1px solid var(--line);border-bottom-width:2px;border-radius:4px;
    padding:0 5px;font-family:ui-monospace,monospace;font-size:11px;color:var(--mut)}
  .hint{color:var(--mut);font-size:11px}
  .empty{padding:40px;text-align:center;color:var(--mut)}
</style>
</head>
<body>
<div id="app">
  <header>
    <h1>Telotron ortholog review</h1>
    <span class="stat">total <b id="s-total">0</b></span>
    <span class="pill ok">✓ <b id="s-ok">0</b></span>
    <span class="pill no">✗ <b id="s-no">0</b></span>
    <span class="pill un">• <b id="s-un">0</b></span>
    <span class="spacer"></span>
    <label class="hint">version
      <select id="f-version"><option value="">all</option></select></label>
    <label class="hint">category
      <select id="f-category"><option value="">all</option></select></label>
    <label class="hint">show
      <select id="f-status">
        <option value="">all</option><option value="un">undecided</option>
        <option value="ok">confirmed</option><option value="no">rejected</option>
      </select></label>
    <button id="b-export">Export revise-script</button>
    <button id="b-write" class="primary" title="Chrome/Edge only">Write good set…</button>
    <button id="b-json" title="save/restore decisions">JSON ▾</button>
  </header>
  <div class="bar"><i id="bar"></i></div>
  <aside id="list"></aside>
  <main>
    <div class="meta" id="meta"></div>
    <pre id="view"></pre>
  </main>
</div>
<footer>
  <button class="no" id="b-no">Reject (r)</button>
  <button class="ok" id="b-ok">Confirm (c)</button>
  <button id="b-un">Undecide (u)</button>
  <span class="hint">
    <kbd>j</kbd>/<kbd>k</kbd> next/prev · <kbd>c</kbd> confirm · <kbd>r</kbd> reject ·
    <kbd>u</kbd> undecide · <kbd>n</kbd> next undecided</span>
  <span class="spacer"></span>
  <span class="hint" id="pos"></span>
</footer>

<script>
window.REVIEW_DATA = __DATA__;
window.REVIEW_META = __META__;
</script>
<script>
"use strict";
const DATA = window.REVIEW_DATA;
const META = window.REVIEW_META;

// surface any runtime error on-page (so failures are visible without the console)
window.addEventListener("error", function(e){
  let b=document.getElementById("err");
  if(!b){b=document.createElement("div");b.id="err";
    b.style.cssText="position:fixed;left:0;right:0;top:0;z-index:99999;background:#6a1b1b;color:#fff;padding:8px 14px;font:12px ui-monospace,monospace;white-space:pre-wrap";
    (document.body||document.documentElement).appendChild(b);}
  b.textContent="JS error: "+(e.message||(e.error&&e.error.message)||e)+"  @"+(e.lineno||"")+":"+(e.colno||"");
});
if(!DATA||!DATA.length){ throw new Error("REVIEW_DATA missing/empty — likely a cached old page; hard-reload (Ctrl/Cmd+Shift+R)."); }
const LSKEY = "telotron_ortho_review::" + META.sig;

// ── palettes / detection (ported from the VSCode extension) ──────────────────
const GREY=[144,152,160];
const NT={A:[78,154,107],C:[79,134,198],G:[199,154,62],T:[203,106,99],U:[203,106,99],
  R:GREY,Y:GREY,S:GREY,W:GREY,K:GREY,M:GREY,B:GREY,D:GREY,H:GREY,V:GREY,N:[130,138,146]};
const AA={A:[102,145,196],V:[102,145,196],L:[102,145,196],I:[102,145,196],M:[102,145,196],
  F:[102,145,196],W:[102,145,196],K:[203,106,99],R:[203,106,99],D:[176,124,198],E:[176,124,198],
  N:[90,160,111],Q:[90,160,111],S:[90,160,111],T:[90,160,111],C:[217,143,176],G:[210,154,82],
  P:[182,168,67],H:[79,163,163],Y:[79,163,163]};
const NT_A=new Set("ACGTUNRYSWKMBDHV"),AA_A=new Set("ACDEFGHIKLMNPQRSTVWYBXZ"),
  PROT=new Set("EFILPQZJO"),GAP=new Set("-.*~ ");
const NT_MINLEN=12,NT_FRAC=.9,AA_MINLEN=24,AA_FRAC=.9,AA_MIND=8;
const CONS_MIN_ROWS=3,CONS_RGB=[150,165,190],CONS_BG_MAX=.24;
const RUN_RE=/[A-Za-z.*~-]+/g, TAG_RE=/<\/?color>/g;

function classify(run,forced){
  let n=0,nt=0,aa=0,p=0,lo=0,up=0;const seen=new Set();
  for(const ch of run){ if(GAP.has(ch))continue; const u=ch.toUpperCase();
    if(u<"A"||u>"Z")continue; n++;seen.add(u);
    if(ch>="a"&&ch<="z")lo++;else up++;
    if(NT_A.has(u))nt++; if(AA_A.has(u))aa++; if(PROT.has(u))p++; }
  if(!n)return null;
  if(forced)return p>0?"aa":"nt";
  const uni=lo===0||up===0;
  if(nt/n>=NT_FRAC&&n>=NT_MINLEN&&seen.size>=2&&uni)return"nt";
  if(aa/n>=AA_FRAC&&n>=AA_MINLEN&&seen.size>=AA_MIND&&uni)return"aa";
  return null;
}
const esc=s=>s.replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));

// per-line forced-ranges (inside <color>..</color>) and tag-char ranges
function lineMarks(line){
  const tags=[]; let m; TAG_RE.lastIndex=0;
  while((m=TAG_RE.exec(line))) tags.push([m.index,m.index+m[0].length,m[0][1]!=="/"]);
  const forced=[]; let open=null;
  for(const[s,e,isOpen]of tags){ if(isOpen)open=e; else if(open!==null){forced.push([open,s]);open=null;} }
  return {tags,forced};
}
const inRanges=(c,rs)=>rs.some(([s,e])=>c>=s&&c<e);

// render one locus text -> colored HTML (keeps tags, dims them; conservation = bg glow)
function renderLocus(text){
  const lines=text.split("\n");
  const marks=lines.map(lineMarks);
  // pass 1: collect typed runs for conservation
  const runs=[]; // {li,col,seq,type}
  lines.forEach((line,li)=>{
    RUN_RE.lastIndex=0; let mm;
    while((mm=RUN_RE.exec(line))){
      const col=mm.index, seq=mm[0];
      const forced=inRanges(col,marks[li].forced);
      const t=classify(seq,forced);
      if(t) runs.push({li,col,seq,type:t});
    }
  });
  // pass 2: conservation -> per-run bg alpha array
  const bg=new Map(); // key li:col -> Float array
  const byKey=new Map();
  for(const r of runs){ const k=r.col+"|"+r.seq.length+"|"+r.type;
    (byKey.get(k)||byKey.set(k,[]).get(k)).push(r); }
  for(const g of byKey.values()){
    g.sort((a,b)=>a.li-b.li);
    for(let i=0;i<g.length;){
      let j=i; while(j+1<g.length&&g[j+1].li===g[j].li+1)j++;
      const blk=g.slice(i,j+1); i=j+1;
      if(blk.length<CONS_MIN_ROWS)continue;
      const L=blk[0].seq.length;
      const alpha=new Float32Array(L);
      for(let p=0;p<L;p++){
        const cnt={};let best=0,nn=0;
        for(const r of blk){ const ch=r.seq[p],u=ch.toUpperCase();
          if(GAP.has(ch)||u<"A"||u>"Z")continue; nn++; cnt[u]=(cnt[u]||0)+1; if(cnt[u]>best)best=cnt[u]; }
        if(nn<2){alpha[p]=0;continue;}
        const cons=best/blk.length;
        alpha[p]=Math.max(0,(cons-.5)/.5)*CONS_BG_MAX;
      }
      for(const r of blk) bg.set(r.li+":"+r.col,alpha);
    }
  }
  // pass 3: build HTML
  let out="";
  lines.forEach((line,li)=>{
    const {tags,forced}=marks[li];
    const typed=new Map(); // col -> {seq,type}
    RUN_RE.lastIndex=0; let mm;
    while((mm=RUN_RE.exec(line))){ const col=mm.index,seq=mm[0];
      const t=classify(seq,inRanges(col,forced)); if(t)typed.set(col,{seq,type:t}); }
    let html="",c=0;
    while(c<line.length){
      // tag?
      const tg=tags.find(t=>t[0]===c);
      if(tg){ html+='<span class="tag">'+esc(line.slice(tg[0],tg[1]))+"</span>"; c=tg[1]; continue; }
      const run=typed.get(c);
      if(run){
        const pal=run.type==="nt"?NT:AA, al=bg.get(li+":"+c);
        for(let p=0;p<run.seq.length;p++){
          const ch=run.seq[p],u=ch.toUpperCase(),rgb=pal[u];
          const a=al?al[p]:0;
          let st="";
          if(rgb)st+="color:rgb("+rgb[0]+","+rgb[1]+","+rgb[2]+");";
          if(a>0.001)st+="background:rgba("+CONS_RGB[0]+","+CONS_RGB[1]+","+CONS_RGB[2]+","+a.toFixed(3)+");";
          html+= st? "<span style=\""+st+"\">"+esc(ch)+"</span>" : esc(ch);
        }
        c+=run.seq.length; continue;
      }
      // plain char run until next tag/typed-run start
      let nxt=line.length;
      for(const t of tags) if(t[0]>c&&t[0]<nxt)nxt=t[0];
      for(const k of typed.keys()) if(k>c&&k<nxt)nxt=k;
      html+=esc(line.slice(c,nxt)); c=nxt;
    }
    out+=html+"\n";
  });
  return out;
}

// ── state ────────────────────────────────────────────────────────────────────
let decisions={}; try{decisions=JSON.parse(localStorage.getItem(LSKEY))||{};}catch(e){}
let cur=0, view=[]; // view = filtered indices into DATA

const $=id=>document.getElementById(id);
function save(){ try{localStorage.setItem(LSKEY,JSON.stringify(decisions));}catch(e){} }

function applyFilter(){
  const v=$("f-version").value,c=$("f-category").value,s=$("f-status").value;
  view=[];
  DATA.forEach((d,i)=>{
    if(v&&d.version!==v)return; if(c&&d.category!==c)return;
    const st=decisions[d.id]||"un";
    if(s&&st!==s)return;
    view.push(i);
  });
  if(view.length===0){cur=0;}
  else if(cur>=view.length)cur=view.length-1;
  renderList(); renderCur(); updateStats();
}

function renderList(){
  const byGrp={};
  view.forEach(i=>{const d=DATA[i];const g=d.version+" / "+d.category;(byGrp[g]||(byGrp[g]=[])).push(i);});
  let html="";
  for(const g of Object.keys(byGrp)){
    html+='<div class="grp">'+esc(g)+" ("+byGrp[g].length+")</div>";
    for(const i of byGrp[g]){
      const d=DATA[i],st=decisions[d.id]||"un";
      const isCur = view[cur]===i;
      html+='<div class="row'+(isCur?" cur":"")+'" data-i="'+i+'">'+
            '<span class="dot '+(st==="un"?"":st)+'"></span>'+
            '<span class="nm">'+esc(d.file.replace(/\.txt$/,""))+"</span></div>";
    }
  }
  $("list").innerHTML=html||'<div class="empty">no loci match the filter</div>';
  const c=$("list").querySelector(".row.cur"); if(c&&c.scrollIntoView)c.scrollIntoView({block:"nearest"});
}

function renderCur(){
  if(view.length===0){ $("meta").innerHTML=""; $("view").innerHTML='<div class="empty">nothing to show</div>'; $("pos").textContent=""; return; }
  const d=DATA[view[cur]], st=decisions[d.id]||"un";
  const ndrop = d.n_unaligned_dropped || 0;
  $("meta").innerHTML=
    '<span class="badge '+d.category+'">'+esc(d.category)+"</span>"+
    '<span class="badge">'+esc(d.version)+"</span>"+
    '<span class="fn">'+esc(d.file)+"</span>"+
    '<span class="pill '+(st==="un"?"un":st)+'">'+(st==="ok"?"✓ confirmed":st==="no"?"✗ rejected":"• undecided")+"</span>"+
    (ndrop>0?'<span class="badge" title="UNALIGNED section(s) stripped from view; ALIGNED blocks carry the same content. Open the source .txt to see the raw unaligned material.">⚠ '+ndrop+' UNALIGNED block'+(ndrop===1?'':'s')+' hidden</span>':'');
  $("view").innerHTML=renderLocus(d.text);
  $("view").parentElement.scrollTop=0;
  $("pos").textContent=(cur+1)+" / "+view.length;
}

function updateStats(){
  let ok=0,no=0; for(const k in decisions){ if(decisions[k]==="ok")ok++; else if(decisions[k]==="no")no++; }
  $("s-total").textContent=DATA.length; $("s-ok").textContent=ok; $("s-no").textContent=no;
  $("s-un").textContent=DATA.length-ok-no;
  $("bar").style.width=(100*(ok+no)/DATA.length)+"%";
}

function decide(st){ if(view.length===0)return; const d=DATA[view[cur]];
  if(st==="un")delete decisions[d.id]; else decisions[d.id]=st;
  save(); updateStats();
  // update just the dot + meta, then advance
  const dot=$("list").querySelector('.row[data-i="'+view[cur]+'"] .dot');
  if(dot)dot.className="dot "+(st==="un"?"":st);
  if($("f-status").value){ applyFilter(); } else { renderCur(); if(st!=="un")move(1); }
}
function move(dir){ if(view.length===0)return; cur=Math.max(0,Math.min(view.length-1,cur+dir)); renderList(); renderCur(); }
function nextUndecided(){ for(let k=1;k<=view.length;k++){ const idx=(cur+k)%view.length;
  const d=DATA[view[idx]]; if(!decisions[d.id]){cur=idx;renderList();renderCur();return;} } }

// ── export ───────────────────────────────────────────────────────────────────
function confirmedRecs(){ return DATA.filter(d=>decisions[d.id]==="ok"); }

function exportScript(){
  const recs=confirmedRecs();
  const lines=["#!/usr/bin/env bash",
    "# Revise the original v2 locus_text GROUPS in place (no copying):",
    "#   * drop any locus whose focal-frame aa-flank alignment has <7 columns identical",
    "#     across all orthologs (left+right flank combined)",
    "#   * keep the loci confirmed in this reviewer (heredoc below)",
    "# Dropped loci are MOVED to work/results/telotron_orthologs_v2/locus_text_dropped/<reason>/",
    "# Reversible:  python scripts/curate_locus_text.py --restore",
    "set -euo pipefail",
    'cd "'+META.repo+'"',
    "python scripts/curate_locus_text.py --apply --min 7 --confirmed-file <(cat <<'IDS'"];
  for(const d of recs) lines.push(d.id);
  lines.push("IDS", ")", "");
  dl("apply_good_orthologs.sh", lines.join("\n"), "text/x-shellscript");
}

async function writeFolder(){
  const recs=confirmedRecs();
  if(!recs.length){alert("No confirmed loci yet.");return;}
  if(!window.showDirectoryPicker){alert("This browser lacks the File System Access API.\nUse “Export copy-script” instead (Chrome/Edge support direct write).");return;}
  let base; try{ base=await window.showDirectoryPicker({mode:"readwrite"}); }catch(e){return;}
  if(!confirm("Write "+recs.length+" confirmed loci into good_orthologs/<version>/<category>/ inside the chosen folder?"))return;
  const root=await base.getDirectoryHandle("good_orthologs",{create:true});
  let n=0;
  for(const d of recs){
    const vd=await root.getDirectoryHandle(d.version,{create:true});
    const cd=await vd.getDirectoryHandle(d.category,{create:true});
    const fh=await cd.getFileHandle(d.file,{create:true});
    const w=await fh.createWritable(); await w.write(d.text); await w.close(); n++;
  }
  alert("Wrote "+n+" files into good_orthologs/ under the chosen folder.");
}

function jsonMenu(){
  const choice=prompt("Type 'save' to download decisions, or paste a decisions JSON to import:","save");
  if(choice===null)return;
  if(choice.trim()==="save"){ dl("decisions.json",JSON.stringify(decisions,null,0),"application/json"); return; }
  try{ const obj=JSON.parse(choice); if(obj&&typeof obj==="object"){decisions=obj;save();applyFilter();alert("Imported "+Object.keys(obj).length+" decisions.");} }
  catch(e){ alert("Not valid JSON."); }
}

function dl(name,text,type){ const b=new Blob([text],{type});const u=URL.createObjectURL(b);
  const a=document.createElement("a");a.href=u;a.download=name;a.click();URL.revokeObjectURL(u); }

// ── wire up ──────────────────────────────────────────────────────────────────
(function init(){
  const vers=[...new Set(DATA.map(d=>d.version))].sort();
  const cats=[...new Set(DATA.map(d=>d.category))];
  const catOrder=__CATORDER__;
  cats.sort((a,b)=>catOrder.indexOf(a)-catOrder.indexOf(b));
  for(const v of vers)$("f-version").add(new Option(v,v));
  for(const c of cats)$("f-category").add(new Option(c,c));
  $("f-version").onchange=$("f-category").onchange=$("f-status").onchange=()=>{cur=0;applyFilter();};
  $("list").onclick=e=>{const r=e.target.closest(".row");if(!r)return;
    const i=+r.dataset.i; cur=view.indexOf(i); if(cur<0)cur=0; renderList(); renderCur();};
  $("b-ok").onclick=()=>decide("ok"); $("b-no").onclick=()=>decide("no"); $("b-un").onclick=()=>decide("un");
  $("b-export").onclick=exportScript; $("b-write").onclick=writeFolder; $("b-json").onclick=jsonMenu;
  document.addEventListener("keydown",e=>{
    if(/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName))return;
    const k=e.key.toLowerCase();
    if(k==="j"||e.key==="ArrowDown"){e.preventDefault();move(1);}
    else if(k==="k"||e.key==="ArrowUp"){e.preventDefault();move(-1);}
    else if(k==="c"||k==="y"){decide("ok");}
    else if(k==="r"){decide("no");}
    else if(k==="u"){decide("un");}
    else if(k==="n"){nextUndecided();}
  });
  applyFilter();
})();
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parents[1]
    ap.add_argument("--root", default=here / "work/results/telotron_orthologs_v2/locus_text", type=Path)
    ap.add_argument("--out", default=here / "work/results/telotron_orthologs_v2/review.html", type=Path)
    ap.add_argument("--dest", default=None,
                    help="good_orthologs destination base shown in the copy-script "
                         "(default: alongside locus_text)")
    args = ap.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.exit(f"locus_text root not found: {root}")
    recs = collect(root)
    if not recs:
        sys.exit(f"no .txt files under {root}")

    dest = args.dest or str(root.parent / "good_orthologs")
    # Hash the immutable per-record ID set so adding/dropping a single locus
    # does not invalidate the entire reviewer's localStorage state.
    id_digest = hashlib.sha1("\n".join(sorted(r["id"] for r in recs)).encode("utf-8")).hexdigest()[:16]
    meta = {"sig": id_digest, "dest": dest, "root": str(root), "repo": str(here)}

    data_json = json.dumps(recs, ensure_ascii=False).replace("</", "<\\/")
    meta_json = json.dumps(meta, ensure_ascii=False).replace("</", "<\\/")
    html = (HTML
            .replace("__DATA__", data_json)
            .replace("__META__", meta_json)
            .replace("__CATORDER__", json.dumps(CATEGORY_ORDER)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    n_by = {}
    for r in recs:
        n_by[f"{r['version']}/{r['category']}"] = n_by.get(f"{r['version']}/{r['category']}", 0) + 1
    print(f"wrote {args.out}  ({len(recs)} loci, {args.out.stat().st_size/1048576:.1f} MB)")
    for k in sorted(n_by):
        print(f"   {k:42s} {n_by[k]}")
    print("\nopen the HTML in a browser; 'Export revise-script' downloads apply_good_orthologs.sh")
    print("which prunes the original groups in place (keeps confirmed, drops <7-conserved).")


if __name__ == "__main__":
    main()
