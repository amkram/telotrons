#!/usr/bin/env python3
"""Telotron formation, revised: two deposition modes (dominant single-arm vs central-seed bidirectional),
both filtered by the splice code. ancestral intron -> telotron."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch
EX="#555555"; INTRON="#d9c9a3"; CR="#b2182b"; GR="#2166ac"; LINK="#9e9e9e"; SS="#1a7a1a"; BRK="#d62728"; SEED="#6a51a3"
fig,ax=plt.subplots(figsize=(13.5,16.5)); ax.set_xlim(0,14); ax.set_ylim(0,40); ax.axis("off")
def locus(y,x0,parts,h=0.92):
    x=x0
    for w,c,lab,tc,fs in parts:
        ax.add_patch(Rectangle((x,y-h/2),w,h,facecolor=c,edgecolor="white",lw=.6))
        if lab: ax.text(x+w/2,y,lab,ha="center",va="center",color=tc,fontsize=fs,weight="bold")
        x+=w
    return x
def gt(x,y,t="GT",c=SS): ax.text(x,y+0.72,t,ha="center",fontsize=7.5,color=c,weight="bold")
def ag(x,y,t="AG",c=SS): ax.text(x,y+0.72,t,ha="center",fontsize=7.5,color=c,weight="bold")
def varrow(x,y0,y1,txt="",side="right",fs=8):
    ax.add_patch(FancyArrowPatch((x,y0),(x,y1),arrowstyle="-|>",mutation_scale=16,color="#444",lw=1.6))
    if txt: ax.text(x+(0.3 if side=="right" else -0.3),(y0+y1)/2,txt,ha="left" if side=="right" else "right",
                    va="center",fontsize=fs,style="italic",color="#444")

# ---- Stage 0: ancestral intron (centered) ----
def conn(y0,y1,txt,fs=7.6):
    ax.add_patch(FancyArrowPatch((7,y0),(7,y1),arrowstyle="-|>",mutation_scale=16,color="#444",lw=1.6))
    ax.text(7,y1-0.42,txt,ha="center",va="center",fontsize=fs,style="italic",color="#444")
y=39.0
ax.text(7,y,"0   ordinary GT–AG intron in a conserved gene (retained by telotron-free sisters)",ha="center",fontsize=9.3,weight="bold",color="#222")
y=38.0
locus(y,3.6,[(0.9,EX,"5′ exon","white",7.5),(6.2,INTRON,"ancestral intron","#5a4a1a",8),(0.9,EX,"3′ exon","white",7.5)])
gt(4.5,y); ag(9.7,y)
conn(37.0,36.4,"break at a fragile (TAD-boundary) domain → resection → outward-facing 3′-OH ends")

# ---- Stage 1: de-novo telomerase ----
y=35.3
ax.text(7,y,"1   de-novo telomerase writes TTTAGGG at the break (fixed template; no copy, no microhomology)",ha="center",fontsize=9.3,weight="bold",color="#222")
y=34.2
locus(y,3.6,[(0.9,EX,"5′ exon","white",7.5),(2.5,INTRON,"","",7)])
ax.text(7.0,y,"⚡",ha="center",va="center",fontsize=16,color=BRK)
locus(y,7.6,[(2.5,INTRON,"","",7),(0.9,EX,"3′ exon","white",7.5)])
conn(33.0,32.4,"array orientation depends on WHICH 3′ end is extended  →  TWO modes")

# ===================== SPLIT =====================
LX=3.4; RX=10.6  # column centers
ax.text(LX,30.7,"DOMINANT  (~72%)   single arm",ha="center",fontsize=9.2,weight="bold",color=GR)
ax.text(RX,30.7,"CONVERGENT  (~21–26%)   central seed",ha="center",fontsize=9.2,weight="bold",color=SEED)
ax.plot([7,7],[31.0,9.0],color="#dddddd",lw=1,ls=":")

# ---- LEFT: single arm ----
ax.text(LX,29.7,"telomerase extends ONE 3′ end → a single G-rich tract",ha="center",fontsize=7.6,color="#333")
y=28.4
locus(y,0.4,[(0.8,EX,"5′","white",7),(1.7,LINK,"linker / intron","#222",6.3),(2.3,GR,"G-rich →","white",6.8),(0.8,EX,"3′","white",7)])
gt(1.2,y,"GT"); ag(5.2,y,"AG")
ax.text(1.25,y-0.95,"ancestral donor",ha="center",fontsize=6,color="#777")
ax.text(4.3,y-0.95,"array supplies AG acceptor",ha="center",fontsize=6,color="#777")
ax.text(LX,26.9,"orientation is ~50/50 at deposition\n→ SPLICE-SELECTED: C-rich-only arms have no\nAG acceptor and are purged; G-rich arms survive",
        ha="center",fontsize=7,color="#333")
ax.text(LX,25.3,"(array sits in the intron body; in ~53% it reaches\nNEITHER edge and the ancestral splice sites are used)",
        ha="center",fontsize=6.5,style="italic",color="#888")

# ---- RIGHT: central seed geometry inset ----
bx=FancyBboxPatch((7.7,26.6),5.7,3.1,boxstyle="round,pad=0.1",facecolor="#f4f1f9",edgecolor=SEED,lw=1)
ax.add_patch(bx)
ax.text(10.55,29.35,"telomerase extends BOTH 3′ ends of the central fragment, OUTWARD",ha="center",fontsize=6.9,weight="bold",color=SEED)
# duplex seed with two outward 3' ends
sy=28.1
ax.plot([9.9,11.2],[sy+0.18,sy+0.18],color="#333",lw=2)   # top strand
ax.plot([9.9,11.2],[sy-0.18,sy-0.18],color="#333",lw=2)   # bottom strand
ax.text(10.55,sy+0.62,"captured seed (duplex)",ha="center",fontsize=6,color="#555")
ax.text(9.78,sy-0.18,"3′",ha="right",va="center",fontsize=6.5,color=BRK,weight="bold")   # bottom-left 3'
ax.text(11.32,sy+0.18,"3′",ha="left",va="center",fontsize=6.5,color=BRK,weight="bold")    # top-right 3'
ax.add_patch(FancyArrowPatch((9.85,sy-0.18),(8.55,sy-0.18),arrowstyle="-|>",mutation_scale=10,color=CR,lw=1.6))
ax.add_patch(FancyArrowPatch((11.25,sy+0.18),(12.55,sy+0.18),arrowstyle="-|>",mutation_scale=10,color=GR,lw=1.6))
ax.text(8.5,sy-0.62,"G-rich on bottom\n= C-rich on coding",ha="center",fontsize=5.6,color=CR)
ax.text(12.6,sy+0.62,"G-rich on coding",ha="center",fontsize=5.6,color=GR)
ax.text(10.55,26.85,"polarity FORCED by geometry (same for + and − genes)",ha="center",fontsize=6,style="italic",color=SEED)
varrow(RX,26.4,25.6,"",fs=7)
# right resulting locus
y=25.0
locus(y,7.5,[(0.7,EX,"5′","white",6.5),(2.1,CR,"C-rich","white",6.5),(0.9,LINK,"seed","#222",6),(2.1,GR,"G-rich","white",6.5),(0.7,EX,"3′","white",6.5)])
gt(8.2,y,"GT"); ag(12.6,y,"AG")
ax.text(RX,23.6,"= GT–Crich–linker–Grich–AG  (the convergent architecture)\nmade in ONE event; geometry & splice-need are confluent",
        ha="center",fontsize=7,color="#333")
ax.text(RX,22.3,"(the ~58 bp central A⁻|A⁺ linker = this seed;\nthe ~11 bp A⁺|A⁺ seams = successive same-arm elongation)",
        ha="center",fontsize=6.5,style="italic",color="#888")

# ===================== MERGE =====================
ax.add_patch(FancyArrowPatch((LX,24.6),(6.7,20.4),arrowstyle="-|>",mutation_scale=15,color="#444",lw=1.5))
ax.add_patch(FancyArrowPatch((RX,21.6),(7.3,20.4),arrowstyle="-|>",mutation_scale=15,color="#444",lw=1.5))
ax.text(7,20.0,"spliced out at normal efficiency (≈0.96)  →  protein intact  →  ~neutral  →  FIXES",
        ha="center",fontsize=8.4,weight="bold",style="italic",color="#444")
y=18.7
locus(y,3.6,[(0.9,EX,"5′ exon","white",7.5),(6.2,"#eeeeee","mature mRNA: exons joined, telotron removed","#444",8),(0.9,EX,"3′ exon","white",7.5)])

# WHY box
bx=FancyBboxPatch((0.5,14.3),13,3.0,boxstyle="round,pad=0.12",facecolor="#f3f6fb",edgecolor="#2166ac",lw=1.2); ax.add_patch(bx)
ax.text(7,16.85,"WHY G-rich is always 3′ (the asymmetry behind both modes)",ha="center",fontsize=9.3,weight="bold",color="#2166ac")
ax.text(7,16.0,"The spliceosome needs the intron to read GT…AG.  Only the G-rich strand (TTTAGGG) carries the obligatory 3′ AG acceptor —",ha="center",fontsize=7.8)
ax.text(7,15.35,"C-rich (CCCTAAA) has neither GT nor AG.  So G-rich must end up at the 3′/acceptor end; C-rich can only flank the 5′ donor (GT|AAACCCT).",ha="center",fontsize=7.8)
ax.text(7,14.65,"Single-arm: enforced by SELECTION (wrong orientations purged).   Central-seed: enforced by GEOMETRY (outward extension) — and consistent with it.",ha="center",fontsize=7.8,weight="bold")

# OUTCOME
ax.text(7,13.2,"OUTCOME",ha="center",fontsize=9.3,weight="bold",color="#444")
ax.text(7,12.5,"Recurs at fragile (TAD-boundary) domains, genome-wide, across the Eimeria radiation → thousands of telotrons.",ha="center",fontsize=7.8)
ax.text(7,11.85,"≈90% fill a pre-existing intron;  ≈6% where the array bootstraps BOTH its own splice sites = de-novo intron birth.",ha="center",fontsize=7.8)
ax.text(7,11.2,"Telotrons = the intronic, spliced-out, tolerated subset of genome-wide de-novo telomerase array deposition.",ha="center",fontsize=7.8,style="italic",color="#333")

# legend
ly=9.4
for i,(c,t) in enumerate([(EX,"exon"),(INTRON,"ancestral intron"),(GR,"G-rich array (TTTAGGG)"),(CR,"C-rich array (CCCTAAA)"),(LINK,"linker / seed / stub")]):
    ax.add_patch(Rectangle((1.2+i*2.4,ly),0.35,0.5,facecolor=c,edgecolor="white"))
    ax.text(1.65+i*2.4,ly+0.25,t,va="center",fontsize=7)
ax.text(7,8.5,"open: physical basis of the fragile-domain breaks · why Eimeria · insertion-vs-replacement of the original intron body · central-seed prediction (single contiguous linker, de-novo junctions both sides) untested",
        ha="center",fontsize=6.5,style="italic",color="#999")
ax.set_title("How a telotron forms:  one break, two deposition modes (single-arm vs central-seed), both filtered by the splice code",
             fontsize=12.5,weight="bold",y=0.995)
fig.savefig("analysis/telotron_mechanism_diagram.png",dpi=160,bbox_inches="tight")
fig.savefig("analysis/telotron_mechanism_diagram.pdf",bbox_inches="tight")
print("wrote analysis/telotron_mechanism_diagram.png")
