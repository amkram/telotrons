#!/usr/bin/env python3
"""How the old intron body is removed and replaced: DSB -> resection (removes ~half the body) ->
de-novo telomerase fill -> ligation. Quantified from linker_segmentation (n=343 focal telotrons)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch
EX="#555555"; BODY="#d9c9a3"; CR="#b2182b"; GR="#2166ac"; STUB="#6a51a3"; SS="#1a7a1a"; BRK="#d62728"; REM="#cccccc"
fig,ax=plt.subplots(figsize=(13,15)); ax.set_xlim(0,13); ax.set_ylim(0,32); ax.axis("off")
XL,XR=2.0,11.0; E=1.0; Ia,Ib=XL+E,XR-E   # exon ends; intron span Ia..Ib
def box(x0,w,y,c,lab="",h=0.85,tc="white",fs=8):
    ax.add_patch(Rectangle((x0,y-h/2),w,h,facecolor=c,edgecolor="white",lw=.6))
    if lab: ax.text(x0+w/2,y,lab,ha="center",va="center",color=tc,fontsize=fs,weight="bold")
def strand(x0,x1,y,color="#444",lw=2.4,ls="-"):
    ax.plot([x0,x1],[y,y],color=color,lw=lw,ls=ls,solid_capstyle="butt")
def three_prime(x,y,direction,color="#444"):  # 3' arrowhead + label
    dx=0.5*direction
    ax.add_patch(FancyArrowPatch((x,y),(x+dx,y),arrowstyle="-|>",mutation_scale=12,color=color,lw=2.2))
    ax.text(x+dx+0.12*direction,y+0.0,"3′",ha="left" if direction>0 else "right",va="center",fontsize=6.5,color=color)
def conn(y0,y1,txt,fs=8.2):
    ax.add_patch(FancyArrowPatch((6.5,y0),(6.5,y1),arrowstyle="-|>",mutation_scale=16,color="#333",lw=1.7))
    ax.text(6.5,y1-0.45,txt,ha="center",va="center",fontsize=fs,style="italic",color="#333")
def stage(y,n,t): ax.text(0.15,y,n,ha="left",va="center",fontsize=12,weight="bold",color="#222"); ax.text(11.35,y,t,ha="left",va="center",fontsize=8,color="#444")

# A: ancestral intron (duplex)
y=30.2; stage(y,"0","ancestral intron\n(ordinary, ~200 bp)")
box(XL,E,y,EX,"5′ exon"); box(Ia,Ib-Ia,y,BODY,"ancestral intron body  ~200 bp",tc="#5a4a1a",fs=8); box(Ib,E,y,EX,"3′ exon")
ax.text(Ia+0.25,y+0.7,"GT",color=SS,fontsize=7.5,weight="bold"); ax.text(Ib-0.25,y+0.7,"AG",color=SS,fontsize=7.5,weight="bold")
conn(29.5,28.6,"① double-strand break inside the intron")

# B: DSB
y=28.0; stage(y,"1","DSB → two ends")
box(XL,E,y,EX,"5′ exon"); box(Ia,2.3,y,BODY); box(Ib-2.3,2.3,y,BODY); box(Ib,E,y,EX,"3′ exon")
ax.text((Ia+Ib)/2,y,"⚡",ha="center",va="center",fontsize=17,color=BRK)
conn(27.3,26.0,"② 5′→3′ END RESECTION (MRN–CtIP → EXO1, DNA2/BLM) degrades DNA outward from the break")

# C: resection — strand level. central stub survives, two outward 3' ends; flanking body removed
y=25.0; stage(y,"2","resection removes\n~half the body")
box(XL,E,y,EX,"5′ exon"); box(Ib,E,y,EX,"3′ exon")
# removed (resected) regions: faded hatched boxes
for x0,w in [(Ia,2.4),(Ib-2.4,2.4)]:
    ax.add_patch(Rectangle((x0,y-0.42),w,0.85,facecolor=REM,edgecolor="#aaa",hatch="xx",lw=0.5,alpha=0.6))
ax.text(Ia+1.2,y-0.95,"~50 bp degraded",ha="center",fontsize=6,color="#999")
ax.text(Ib-1.2,y-0.95,"~50 bp degraded",ha="center",fontsize=6,color="#999")
# central stub (short duplex) with two outward 3' overhangs
cx0,cx1=6.0,7.0
strand(cx0,cx1,y+0.13,STUB,3); strand(cx0,cx1,y-0.13,STUB,3)
ax.text((cx0+cx1)/2,y+0.72,"central stub",ha="center",fontsize=6.5,color=STUB,weight="bold")
three_prime(cx0,y-0.13,-1,STUB)   # bottom strand 3' points LEFT (outward)
three_prime(cx1,y+0.13,+1,STUB)   # top strand 3' points RIGHT (outward)
conn(24.0,22.6,"③ de-novo TELOMERASE extends each 3′ end OUTWARD, writing TTTAGGG (no template copy, no microhomology)")

# D: telomerase fills outward from the stub's two 3' ends
y=21.6; stage(y,"3","telomerase fills\nthe resected gaps")
box(XL,E,y,EX,"5′ exon"); box(Ib,E,y,EX,"3′ exon")
strand(cx0,cx1,y+0.13,STUB,3); strand(cx0,cx1,y-0.13,STUB,3)   # stub
# arrays growing outward (G-rich strand drawn as colored bar on the extending strand)
ax.add_patch(Rectangle((Ia+0.05,y-0.4),cx0-(Ia+0.05),0.8,facecolor=CR,edgecolor="white",alpha=.92)); ax.text((Ia+cx0)/2,y,"C-rich on coding\n(G-rich made on bottom)",ha="center",va="center",color="white",fontsize=5.8,weight="bold")
ax.add_patch(Rectangle((cx1,y-0.4),(Ib-0.05)-cx1,0.8,facecolor=GR,edgecolor="white",alpha=.92)); ax.text((cx1+Ib)/2,y,"G-rich on coding",ha="center",va="center",color="white",fontsize=6,weight="bold")
three_prime(Ia+0.05,y-0.13,-1,"#7a1020"); three_prime(Ib-0.05,y+0.13,+1,"#10325a")
ax.text((cx0+cx1)/2,y+0.72,"stub",ha="center",fontsize=6,color=STUB,weight="bold")
conn(20.6,19.2,"④ second strand filled in + ends ligated to the exons (NHEJ); ancestral splice sites reused")

# E: product telotron
y=18.2; stage(y,"4","telotron (~212 bp)")
box(XL,E,y,EX,"5′ exon"); box(Ib,E,y,EX,"3′ exon")
w=Ib-Ia; box(Ia,(cx0-Ia),y,CR,"C-rich array",fs=6.8); box(cx0,cx1-cx0,y,STUB,"stub",tc="#fff",fs=6); box(cx1,Ib-cx1,y,GR,"G-rich array",fs=6.8)
ax.text(Ia+0.25,y+0.7,"GT",color=SS,fontsize=7.5,weight="bold"); ax.text(Ib-0.25,y+0.7,"AG",color=SS,fontsize=7.5,weight="bold")

# balance sheet box
bx=FancyBboxPatch((0.6,12.6),11.8,4.3,boxstyle="round,pad=0.15",facecolor="#f6f4fb",edgecolor=STUB,lw=1.2); ax.add_patch(bx)
ax.text(6.5,16.45,"THE MASS BALANCE  (median, n=343 focal telotrons)",ha="center",fontsize=10,weight="bold",color=STUB)
# bars
def barrow(y,segs,label):
    x=2.6
    for w,c,t in segs:
        ax.add_patch(Rectangle((x,y-0.32),w/40,0.64,facecolor=c,edgecolor="white"));
        if t: ax.text(x+(w/40)/2,y,t,ha="center",va="center",color="white",fontsize=6.5,weight="bold")
        x+=w/40
    ax.text(2.4,y,label,ha="right",va="center",fontsize=8,weight="bold")
    return x
ax.text(2.4,15.6,"ancestral intron",ha="right",va="center",fontsize=8,weight="bold")
ax.add_patch(Rectangle((2.6,15.28),200/40,0.64,facecolor=BODY,edgecolor="white")); ax.text(2.6+200/80,15.6,"ordinary body  ~200 bp",ha="center",va="center",fontsize=6.8,color="#5a4a1a",weight="bold")
ax.text(2.4,14.3,"telotron",ha="right",va="center",fontsize=8,weight="bold")
x=2.6
ax.add_patch(Rectangle((x,13.98),94/40,0.64,facecolor=STUB,edgecolor="white")); ax.text(x+94/80,14.3,"retained 94 bp",ha="center",va="center",color="white",fontsize=6.3,weight="bold"); x+=94/40
ax.add_patch(Rectangle((x,13.98),112/40,0.64,facecolor=GR,edgecolor="white")); ax.text(x+112/80,14.3,"new array 112 bp",ha="center",va="center",color="white",fontsize=6.3,weight="bold"); x+=112/40
ax.text(6.5,13.2,"≈106 bp of the old body RESECTED away  →  replaced ~1:1 by ~112 bp telomeric array;  ~94 bp of the body survives as the linker/stub.",
        ha="center",fontsize=7.6,color="#333")
ax.text(6.5,12.85,"net length barely changes — a balanced half-swap, not a big insertion or a clean deletion.",ha="center",fontsize=7.3,style="italic",color="#666")

# caveats
ax.text(6.5,11.6,"Mechanism (resection) is inferred from standard DSB biology + the length deficit — not directly observed; per-locus extent unknown.",ha="center",fontsize=7,style="italic",color="#999")
ax.text(6.5,11.05,"Ancestral length is a within-Eimeria proxy (varies a lot); if the true ancestor was the longer outgroup intron (~400 bp), MORE was removed.",ha="center",fontsize=7,style="italic",color="#999")
ax.text(6.5,10.5,"Single-arm telotrons (~72% genome-wide) are insertion-dominant (little removal); this balanced removal is the convergent/central-seed case.",ha="center",fontsize=7,style="italic",color="#999")

# legend
ly=9.2
for i,(c,t) in enumerate([(EX,"exon"),(BODY,"ancestral intron body"),(GR,"G-rich array"),(CR,"C-rich array"),(STUB,"stub / linker (retained)"),(REM,"resected (removed)")]):
    ax.add_patch(Rectangle((1.0+i*1.95,ly),0.32,0.45,facecolor=c,edgecolor="#999",hatch="xx" if c==REM else None))
    ax.text(1.4+i*1.95,ly+0.22,t,va="center",fontsize=6.6)
ax.set_title("How the old intron body is removed:  break → end-resection (≈half the body) → de-novo telomerase fills the gap",
             fontsize=12.5,weight="bold",y=0.995)
fig.savefig("analysis/telotron_removal_fill_diagram.png",dpi=160,bbox_inches="tight")
fig.savefig("analysis/telotron_removal_fill_diagram.pdf",bbox_inches="tight")
print("wrote analysis/telotron_removal_fill_diagram.png")
