#!/usr/bin/env python3
"""Dissect the telotron gene-bias at the INTRON level (the correct unit: a gene with N introns has N
chances). Decisive question: is the intron-rich/large host-gene bias just OPPORTUNITY (uniform per-intron
rate -> intron-rich genes win by having more trials) or genuine PREFERENCE (intron-rich/large genes have a
HIGHER per-intron telotron rate)? Build every Eimeria intron from the GFFs, mark telotron by coordinate
overlap, then: (1) per-intron rate vs gene intron-count, (2) vs gene length, (3) logistic partialling
log(N)+log(len)+position, (4) per-intron rate vs position within gene. NOTE: intron GC/length are NOT used
as predictors — they are altered by telomerization (circular)."""
import csv, re, glob, bisect
from collections import defaultdict
import numpy as np
from scipy.stats import chi2_contingency
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
GFF={"GCF_000499545.2":"data/raw/refseq/GCF_000499545.2/GCF_000499545.2_ETH001_genomic.gff",
     "GCF_000499385.1":"data/raw/refseq/GCF_000499385.1/GCF_000499385.1_ENH001_genomic.gff",
     "GCF_000499425.1":"data/raw/refseq/GCF_000499425.1/GCF_000499425.1_EAH001_genomic.gff",
     "GCF_000499605.1":"data/raw/refseq/GCF_000499605.1/GCF_000499605.1_EMW001_genomic.gff",
     "GCF_000499745.2":"data/raw/refseq/GCF_000499745.2/GCF_000499745.2_EMH001_genomic.gff"}
# telotron intervals per (genome,seqid)
telo=defaultdict(list)
for r in csv.DictReader(open("work/results/final_telotron_set_architecture.tsv"),delimiter="\t"):
    if r["genome_id"] in GFF:
        telo[(r["genome_id"],r["seqid"])].append((int(r["start"]),int(r["end"])))
for k in telo: telo[k].sort()
def is_telo(genome,seqid,s,e):
    ivs=telo.get((genome,seqid))
    if not ivs: return False
    import bisect as bs
    i=bs.bisect_right([x[0] for x in ivs],e)
    for ts,te in ivs[max(0,i-30):i+1]:
        ov=min(e,te)-max(s,ts)
        if ov>0 and ov/min(e-s,te-ts+1)>0.5: return True
    return False
rows=[]  # genome, gene, N, genelen, pos_idx, ilen, is_telo
for gid,path in GFF.items():
    mex=defaultdict(list); mg={}; glen={}; gseq={}
    for l in open(path):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2]=="gene":
            m=re.search(r'ID=([^;]+)',f[8])
            if m: glen[m.group(1)]=int(f[4])-int(f[3])
        elif f[2]=="mRNA":
            mi=re.search(r'ID=([^;]+)',f[8]); pa=re.search(r'Parent=([^;]+)',f[8])
            if mi and pa: mg[mi.group(1)]=pa.group(1)
        elif f[2]=="exon":
            pa=re.search(r'Parent=([^;]+)',f[8])
            if pa: mex[pa.group(1)].append((int(f[3]),int(f[4]),f[0],f[6]))
    # unique introns per gene
    gintr=defaultdict(set); gstrand={}
    for mr,exs in mex.items():
        g=mg.get(mr)
        if not g or len(exs)<2: continue
        exs=sorted(exs); seqid=exs[0][2]; gstrand[g]=exs[0][3]
        for i in range(len(exs)-1):
            istart=exs[i][1]; iend=exs[i+1][0]
            if iend-istart>=20: gintr[g].add((istart,iend,seqid))
    for g,introns in gintr.items():
        ins=sorted(introns)
        # order 5'->3' by strand for position
        if gstrand.get(g)=="-": ins=ins[::-1]
        N=len(ins)
        for idx,(s,e,seqid) in enumerate(ins):
            rows.append((gid,g,N,glen.get(g,e-s),idx/(N-1) if N>1 else 0.5,e-s,is_telo(gid,seqid,s,e)))
print(f"total Eimeria introns: {len(rows)}; telotron introns matched: {sum(r[6] for r in rows)}")
N=np.array([r[2] for r in rows]); GL=np.array([r[3] for r in rows],float); POS=np.array([r[4] for r in rows]); Y=np.array([1 if r[6] else 0 for r in rows])
# (1) per-intron rate vs gene intron count
print("\n=== per-intron telotron rate vs GENE INTRON COUNT (opportunity=flat; preference=rising) ===")
bins=[(1,1),(2,3),(4,6),(7,10),(11,20),(21,200)]
rate_n=[]
for lo,hi in bins:
    m=(N>=lo)&(N<=hi); tot=m.sum(); pos=Y[m].sum()
    rate_n.append((f"{lo}-{hi}",tot,pos,1e4*pos/tot if tot else 0))
    print(f"  N={lo:>2}-{hi:<3}: {tot:6d} introns, {pos:4d} telotrons, rate {1e4*pos/tot if tot else 0:6.1f} /10k")
# (2) per-intron rate vs gene length
print("\n=== per-intron telotron rate vs GENE LENGTH ===")
qs=np.quantile(GL,[0,.2,.4,.6,.8,1.0]); rate_l=[]
for i in range(5):
    m=(GL>=qs[i])&(GL<=qs[i+1]); tot=m.sum(); pos=Y[m].sum()
    rate_l.append((f"{int(qs[i])}-{int(qs[i+1])}",1e4*pos/tot if tot else 0))
    print(f"  len {int(qs[i]):>5}-{int(qs[i+1]):<6}: {tot:6d} introns, {pos:4d} telotrons, rate {1e4*pos/tot if tot else 0:6.1f} /10k")
# (3) logistic
print("\n=== logistic: is_telotron ~ log(N) + log(genelen) + position (which is independent?) ===")
X=np.column_stack([np.log(N), np.log(GL+1), POS])
Xz=(X-X.mean(0))/X.std(0)
try:
    import statsmodels.api as sm
    mod=sm.Logit(Y, sm.add_constant(Xz)).fit(disp=0)
    for nm,c,p in zip(["const","log(N introns)","log(gene length)","position(5'->3')"],mod.params,mod.pvalues):
        print(f"  {nm:20s} coef {c:+.3f}  p={p:.1e}")
except Exception as ex:
    from sklearn.linear_model import LogisticRegression
    lr=LogisticRegression().fit(Xz,Y)
    for nm,c in zip(["log(N introns)","log(gene length)","position"],lr.coef_[0]): print(f"  {nm:20s} coef {c:+.3f}")
# (4) position within gene
print("\n=== per-intron rate vs POSITION within gene (0=5',1=3') ===")
rate_p=[]
for i in range(5):
    m=(POS>=i/5)&(POS<(i+1)/5); tot=m.sum(); pos=Y[m].sum()
    rate_p.append((f"{i/5:.1f}-{(i+1)/5:.1f}",1e4*pos/tot if tot else 0))
    print(f"  pos {i/5:.1f}-{(i+1)/5:.1f}: {tot:6d} introns, {pos:4d} telotrons, rate {1e4*pos/tot if tot else 0:6.1f} /10k")
# ---- figure ----
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
ax[0].plot([r[3] for r in rate_n],"o-",color="#b2182b"); ax[0].set_xticks(range(len(rate_n))); ax[0].set_xticklabels([r[0] for r in rate_n],fontsize=8)
ax[0].set_xlabel("introns per host gene"); ax[0].set_ylabel("telotron rate /10k introns")
ax[0].set_title("A  per-intron rate is INVERTED-U\npeaks at 4-10 introns, falls for N≥21 (rate 12)",fontsize=9.5,weight="bold")
[ax[0].spines[s].set_visible(False) for s in ("top","right")]
ax[1].plot([r[1] for r in rate_l],"o-",color="#2166ac"); ax[1].set_xticks(range(len(rate_l))); ax[1].set_xticklabels([r[0] for r in rate_l],fontsize=7,rotation=30)
ax[1].set_xlabel("gene length quintile (bp)"); ax[1].set_ylabel("telotron rate /10k introns")
ax[1].set_title("B  rate FALLS with gene length\nbig genes host telotrons by trial-count, NOT per-intron preference",fontsize=9.5,weight="bold")
[ax[1].spines[s].set_visible(False) for s in ("top","right")]
ax[2].plot([r[1] for r in rate_p],"o-",color="#1a7a1a"); ax[2].set_xticks(range(len(rate_p))); ax[2].set_xticklabels([r[0] for r in rate_p],fontsize=8)
ax[2].set_xlabel("position within gene (5'->3')"); ax[2].set_ylabel("telotron rate /10k introns")
ax[2].set_title("C  genuine 5' bias (logistic p=3e-4)\nrate 40-46 in 5' half vs 21 at 3' end",fontsize=9.5,weight="bold")
[ax[2].spines[s].set_visible(False) for s in ("top","right")]
fig.suptitle("Per-intron dissection (Eimeria): the host-gene bias is mostly OPPORTUNITY (intron trial-count) + a genuine 5' positional preference — the LARGEST/longest genes are DISFAVORED per intron",fontsize=10.5,weight="bold",y=1.02)
fig.tight_layout(); fig.savefig("work/results/figures/telotron_per_intron.png",dpi=150,bbox_inches="tight"); print("\nwrote work/results/figures/telotron_per_intron.png")
