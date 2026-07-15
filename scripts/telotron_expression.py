#!/usr/bin/env python3
"""Higher-power EXPRESSION vs TELOTRON-PRESENCE test, now with E. necatrix RNA-seq (174 telotron genes,
~8x E. tenella) + E. tenella (21). Per species and POOLED (expression normalised to within-species median
so depths are comparable). Gene-level host-vs-non-host and per-intron telotron rate vs host-gene expression quintile.
Expression = samtools bedcov depth-sum / gene length.
HEADLINE = the size-controlled OLS coefficient (log10 expr ~ is_host + log N introns), NOT the raw
host-vs-non-host MWU: host genes are intron-rich/long and large genes are generically lower-expression,
so the raw gap is a gene-architecture confound. The figure title states whichever the OLS supports."""
import csv, re, bisect as bs, math
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
np.random.seed(0)
import os
# Persistent gene-coverage TSVs (produced by the rnaseq_gene_coverage Snakemake rule); override dir
# with TELOTRON_RNASEQ_DIR. Previously read /tmp/eten_gene_cov.tsv (ephemeral, non-reproducible).
RNASEQ_DIR=os.environ.get("TELOTRON_RNASEQ_DIR","work/results/rnaseq")
SP={  # species: (gene_cov, gff, genome_id, label)
 "necatrix":(f"{RNASEQ_DIR}/necatrix_gene_cov.tsv","data/raw/refseq/GCF_000499385.1/GCF_000499385.1_ENH001_genomic.gff","GCF_000499385.1","E. necatrix"),
 "tenella":(f"{RNASEQ_DIR}/eten_gene_cov.tsv","data/raw/refseq/GCF_000499545.2/GCF_000499545.2_ETH001_genomic.gff","GCF_000499545.2","E. tenella"),
}
def load_telo(genome_id):
    telo=defaultdict(list)
    for r in csv.DictReader(open("work/results/final_telotron_set_architecture.tsv"),delimiter="\t"):
        if r["genome_id"]==genome_id: telo[r["seqid"]].append((int(r["start"]),int(r["end"])))
    for k in telo: telo[k].sort()
    return telo
def make_is_telo(telo):
    def f(seqid,s,e):
        ivs=telo.get(seqid)
        if not ivs: return False
        i=bs.bisect_right([x[0] for x in ivs],e)
        for ts,te in ivs[max(0,i-30):i+1]:
            ov=min(e,te)-max(s,ts)
            if ov>0 and ov/min(e-s,te-ts+1)>0.5: return True
        return False
    return f
def species_data(cov,gff,genome_id):
    expr={}
    for l in open(cov):
        f=l.rstrip().split("\t"); ln=int(f[2])-int(f[1])
        if ln>0: expr[f[3]]=int(f[6])/ln
    it=make_is_telo(load_telo(genome_id))
    mex=defaultdict(list); mg={}
    for l in open(gff):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2]=="mRNA":
            mi=re.search(r'ID=([^;]+)',f[8]); pa=re.search(r'Parent=([^;]+)',f[8])
            if mi and pa: mg[mi.group(1)]=pa.group(1)
        elif f[2]=="exon":
            pa=re.search(r'Parent=([^;]+)',f[8])
            if pa: mex[pa.group(1)].append((int(f[3]),int(f[4]),f[0]))
    gintr=defaultdict(set)
    for mr,exs in mex.items():
        g=mg.get(mr)
        if not g or len(exs)<2: continue
        exs=sorted(exs)
        for i in range(len(exs)-1):
            if exs[i+1][0]-exs[i][1]>=20: gintr[g].add((exs[i][1],exs[i+1][0],exs[i][2]))
    genes=[]; introns=[]
    for g,ins in gintr.items():
        if g not in expr: continue
        h=any(it(seq,s,e) for s,e,seq in ins); N=len(ins)
        genes.append((expr[g],N,h))
        for s,e,seq in ins: introns.append((expr[g],it(seq,s,e)))
    return genes,introns
allg={}; alli={}
for sp,(cov,gff,gid,lab) in SP.items():
    import os
    if not os.path.exists(cov): print(f"[skip {sp}: {cov} missing]"); continue
    g,i=species_data(cov,gff,gid); allg[sp]=g; alli[sp]=i
    he=[x[0] for x in g if x[2]]; ne=[x[0] for x in g if not x[2]]
    print(f"\n=== {lab} ===  genes {len(g)}  telotron-host {len(he)}")
    print(f"  host median expr {np.median(he):.2f}  non-host {np.median(ne):.2f}  ratio {np.median(he)/np.median(ne):.2f}  MWU p={mannwhitneyu(he,ne).pvalue:.2e}")
    # size-controlled OLS
    E=np.array([x[0] for x in g]); N=np.array([x[1] for x in g]); H=np.array([1 if x[2] else 0 for x in g])
    try:
        import statsmodels.api as sm
        m=sm.OLS(np.log10(E+0.1), sm.add_constant(np.column_stack([H,np.log(N+1)]))).fit()
        print(f"  OLS log10(expr)~host+log(N): host coef {m.params[1]:+.3f} p={m.pvalues[1]:.3f} (size-controlled)")
    except Exception as ex: print("  (no statsmodels)",ex)
    # per-intron rate vs expression quintile
    ie=np.array([x[0] for x in i]); iy=np.array([1 if x[1] else 0 for x in i]); qs=np.quantile(ie,[0,.2,.4,.6,.8,1.0])
    print("  per-intron telotron rate by host-gene expression quintile (/10k):", end=" ")
    for k in range(5):
        m=(ie>=qs[k])&(ie<=qs[k+1]); print(f"{1e4*iy[m].sum()/m.sum():.1f}",end=" ")
    print()
# pooled (normalise expr to within-species median)
pg=[]; pi=[]
for sp in allg:
    med=np.median([x[0] for x in allg[sp]])
    for e,N,h in allg[sp]: pg.append((e/med,N,h))
    for e,h in alli[sp]: pi.append((e/med,h))
if pg:
    he=[x[0] for x in pg if x[2]]; ne=[x[0] for x in pg if not x[2]]
    print(f"\n=== POOLED (within-species-normalised) ===  telotron-host {len(he)}")
    print(f"  host median {np.median(he):.2f}  non-host {np.median(ne):.2f}  ratio {np.median(he)/np.median(ne):.2f}  MWU p={mannwhitneyu(he,ne).pvalue:.2e}")
    E=np.array([x[0] for x in pg]); N=np.array([x[1] for x in pg]); H=np.array([1 if x[2] else 0 for x in pg])
    # pooled size-controlled OLS = the headline (raw host-vs-non-host gap is gene-architecture-confounded)
    ols_coef=ols_p=None
    try:
        import statsmodels.api as sm
        m=sm.OLS(np.log10(E+0.01), sm.add_constant(np.column_stack([H,np.log(N+1)]))).fit()
        ols_coef=m.params[1]; ols_p=m.pvalues[1]
        print(f"  OLS host coef {ols_coef:+.3f} p={ols_p:.2e} (size-controlled, pooled HEADLINE)")
    except Exception: pass
    if ols_coef is not None and ols_p is not None and ols_p<0.05 and ols_coef<0:
        HEADLINE="telotrons sit in LOWER-expression genes even after size control (pooled OLS p=%.1e)"%ols_p
    elif ols_coef is not None and ols_p is not None:
        HEADLINE="host-expression gap NOT significant after size control (pooled OLS coef %+.3f p=%.2g) — gene-architecture confound"%(ols_coef,ols_p)
    else:
        HEADLINE="size-controlled OLS unavailable — raw gap is gene-architecture-confounded"
    # ---- figure ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.3))
    # A necatrix host vs non-host
    if "necatrix" in allg:
        g=allg["necatrix"]; he2=[x[0] for x in g if x[2]]; ne2=[x[0] for x in g if not x[2]]
        ax[0].boxplot([np.log10(np.array(ne2)+.1),np.log10(np.array(he2)+.1)],labels=[f"non-host\nn={len(ne2)}",f"telotron host\nn={len(he2)}"],showfliers=False)
        ax[0].set_title(f"A  E. necatrix (174 host genes)\nmed {np.median(he2):.0f} vs {np.median(ne2):.0f}  p={mannwhitneyu(he2,ne2).pvalue:.0e}",fontsize=9.5,weight="bold")
    ax[0].set_ylabel("log10 gene expression (depth/bp)"); [ax[0].spines[s].set_visible(False) for s in ("top","right")]
    # B pooled
    ax[1].boxplot([np.log10(np.array(ne)+.01),np.log10(np.array(he)+.01)],labels=[f"non-host\nn={len(ne)}",f"telotron host\nn={len(he)}"],showfliers=False)
    ax[1].set_title(f"B  POOLED necatrix+tenella (norm.)\nratio {np.median(he)/np.median(ne):.2f}  p={mannwhitneyu(he,ne).pvalue:.0e}",fontsize=9.5,weight="bold")
    ax[1].set_ylabel("log10 normalised expression"); [ax[1].spines[s].set_visible(False) for s in ("top","right")]
    # C pooled per-intron rate vs expr quintile
    ie=np.array([x[0] for x in pi]); iy=np.array([1 if x[1] else 0 for x in pi]); qs=np.quantile(ie,[0,.2,.4,.6,.8,1.0])
    rate=[1e4*iy[(ie>=qs[k])&(ie<=qs[k+1])].sum()/((ie>=qs[k])&(ie<=qs[k+1])).sum() for k in range(5)]
    ax[2].plot(rate,"o-",color="#b2182b"); ax[2].set_xticks(range(5)); ax[2].set_xticklabels(["Q1\nlow","Q2","Q3","Q4","Q5\nhigh"],fontsize=8)
    ax[2].set_xlabel("host-gene expression quintile"); ax[2].set_ylabel("telotron rate /10k introns")
    ax[2].set_title("C  pooled per-intron rate vs expression\n(necatrix+tenella; ~195 telotrons)",fontsize=9.5,weight="bold")
    [ax[2].spines[s].set_visible(False) for s in ("top","right")]
    fig.suptitle("Expression vs telotron presence (E. necatrix RNA-seq added, 174 host genes).  "+HEADLINE,fontsize=9.5,weight="bold",y=1.04)
    fig.tight_layout(); fig.savefig("work/results/figures/telotron_expression.png",dpi=150,bbox_inches="tight"); print("\nwrote work/results/figures/telotron_expression.png")
