#!/usr/bin/env python3
"""Test EXPRESSION vs TELOTRON PRESENCE in E. tenella (the only species with RNA-seq; N capped at 21
telotron genes -> low power, stated). Two independent measures:
 (1) GENE-LEVEL: bedcov depth/bp per gene, telotron-host vs non-host. Host genes are intron-rich, and
     intron count weakly tracks expression, so we also control for size (OLS log10(expr) ~ is_host + log(N)).
 (2) LOCUS-LEVEL (independent): the splicing-analysis loci (per_locus_counts.tsv) pair each telotron with
     matched non-telotron control introns; `total` junction-spanning reads = local expression proxy. Compare
     TELO vs CTRL total + fraction expressed (>0). Bootstrap CIs throughout given small N."""
import csv, re, bisect as bs, random
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
random.seed(0); np.random.seed(0)
import os
GID="GCF_000499545.2"; GFFp=f"data/raw/refseq/{GID}/{GID}_ETH001_genomic.gff"
# Persistent gene-coverage TSV (rnaseq_gene_coverage rule); was /tmp/eten_gene_cov.tsv (ephemeral).
COV=os.path.join(os.environ.get("TELOTRON_RNASEQ_DIR","work/results/rnaseq"),"eten_gene_cov.tsv")
def boot_med_diff(a,b,n=10000):
    a=np.array(a); b=np.array(b); d=[]
    for _ in range(n): d.append(np.median(np.random.choice(a,len(a)))-np.median(np.random.choice(b,len(b))))
    return np.percentile(d,[2.5,97.5])
# ---- (1) gene-level ----
expr={}
for l in open(COV):
    f=l.rstrip().split("\t"); ln=int(f[2])-int(f[1])
    if ln>0: expr[f[3]]=int(f[6])/ln
telo=defaultdict(list)
for r in csv.DictReader(open("work/results/final_telotron_set_architecture.tsv"),delimiter="\t"):
    if r["genome_id"]==GID: telo[r["seqid"]].append((int(r["start"]),int(r["end"])))
for k in telo: telo[k].sort()
def is_telo(seqid,s,e):
    ivs=telo.get(seqid)
    if not ivs: return False
    i=bs.bisect_right([x[0] for x in ivs],e)
    for ts,te in ivs[max(0,i-30):i+1]:
        ov=min(e,te)-max(s,ts)
        if ov>0 and ov/min(e-s,te-ts+1)>0.5: return True
    return False
mex=defaultdict(list); mg={}
for l in open(GFFp):
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
host_e=[]; non_e=[]; Hn=[]; Nn=[]; En=[]; isH=[]
for g,introns in gintr.items():
    if g not in expr: continue
    h=any(is_telo(seq,s,e) for s,e,seq in introns); N=len(introns)
    Nn.append(N); En.append(expr[g]); isH.append(1 if h else 0)
    (host_e if h else non_e).append(expr[g])
print(f"(1) GENE-LEVEL  host n={len(host_e)} non-host n={len(non_e)}")
print(f"    host median expr {np.median(host_e):.2f}  non-host {np.median(non_e):.2f}  ratio {np.median(host_e)/np.median(non_e):.2f}")
print(f"    MWU p={mannwhitneyu(host_e,non_e).pvalue:.3f}  bootstrap 95%CI median diff {boot_med_diff(host_e,non_e)}")
# size-matched: match each host gene to 20 non-host genes with same intron-count, compare
Nn=np.array(Nn); En=np.array(En); isH=np.array(isH)
bynN=defaultdict(list)
for N,e,h in zip(Nn,En,isH):
    if not h: bynN[N].append(e)
matched=[]
for g,introns in gintr.items():
    if g not in expr: continue
    if any(is_telo(seq,s,e) for s,e,seq in introns):
        pool=bynN.get(len(introns),[])
        if pool: matched.append(np.median(pool))
print(f"    SIZE-MATCHED (same intron-count non-host median): host {np.median(host_e):.2f} vs matched-non-host {np.median(matched):.2f}")
# OLS control for size — this (not the raw MWU) is the headline: it removes the gene-architecture
# confound (host genes are intron-rich, and intron count tracks expression).
ols_coef=ols_p=None
try:
    import statsmodels.api as sm
    y=np.log10(En+0.1); X=sm.add_constant(np.column_stack([isH,np.log(Nn+1)]))
    m=sm.OLS(y,X).fit()
    ols_coef=m.params[1]; ols_p=m.pvalues[1]
    print(f"    OLS log10(expr) ~ is_host + log(N): is_host coef {ols_coef:+.3f} p={ols_p:.3f} (size-controlled HEADLINE)")
except Exception as ex: print("    (statsmodels unavailable)",ex)
# size-controlled verdict drives the headline (raw host-vs-non-host gap is architecture-confounded)
if ols_coef is not None and ols_p is not None and ols_p<0.05 and ols_coef<0:
    HEADLINE="telotrons sit in LOWER-expression genes even after controlling for gene size (OLS p=%.3f)"%ols_p
elif ols_coef is not None and ols_p is not None:
    HEADLINE="raw host-expression gap is NOT significant once gene size is controlled (OLS coef %+.3f p=%.3f) — likely a gene-architecture confound"%(ols_coef,ols_p)
else:
    HEADLINE="size-controlled test unavailable (statsmodels missing) — raw gap is gene-architecture-confounded, interpret with caution"
# ---- (2) locus-level (optional: needs the manual splice-counts build; panel skipped if absent) ----
LOCUS_TSV=os.path.join(os.environ.get("TELOTRON_SPLICE_DIR","data/raw/rnaseq_splice_2026"),"per_locus_counts.tsv")
tt=[]; ct=[]
if os.path.exists(LOCUS_TSV):
    rows=list(csv.DictReader(open(LOCUS_TSV),delimiter="\t"))
    tt=[int(r["total"]) for r in rows if r["class"]=="TELO"]; ct=[int(r["total"]) for r in rows if r["class"]=="CTRL"]
    print(f"\n(2) LOCUS-LEVEL  TELO n={len(tt)} CTRL n={len(ct)}  (`total` junction reads = expression proxy)")
    print(f"    TELO median total {np.median(tt):.0f}  CTRL {np.median(ct):.0f}  MWU p={mannwhitneyu(tt,ct).pvalue:.3f}")
    print(f"    bootstrap 95%CI median diff (TELO-CTRL) {boot_med_diff(tt,ct)}")
    print(f"    fraction expressed (total>0): TELO {np.mean([x>0 for x in tt]):.2f}  CTRL {np.mean([x>0 for x in ct]):.2f}")
else:
    print(f"\n(2) LOCUS-LEVEL skipped: {LOCUS_TSV} not found (manual splice-counts build)")
# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(10.5,4.3))
ax[0].boxplot([np.log10(np.array(non_e)+.1),np.log10(np.array(host_e)+.1)],labels=[f"non-host\nn={len(non_e)}",f"telotron host\nn={len(host_e)}"],showfliers=False)
ax[0].set_ylabel("log10 gene expression (depth/bp)")
ax[0].set_title(f"(1) GENE-LEVEL: host vs non-host\nmed {np.median(host_e):.1f} vs {np.median(non_e):.1f} (p={mannwhitneyu(host_e,non_e).pvalue:.2f}); host LOWER",fontsize=9.5,weight="bold")
[ax[0].spines[s].set_visible(False) for s in ("top","right")]
if tt and ct:
    ax[1].boxplot([np.log10(np.array(ct)+1),np.log10(np.array(tt)+1)],labels=[f"CTRL introns\nn={len(ct)}",f"telotrons\nn={len(tt)}"],showfliers=False)
    ax[1].set_ylabel("log10 junction reads (`total`+1)")
    ax[1].set_title(f"(2) LOCUS-LEVEL: telotron vs control intron\nmed {np.median(tt):.0f} vs {np.median(ct):.0f} (p={mannwhitneyu(tt,ct).pvalue:.2f})",fontsize=9.5,weight="bold")
else:
    ax[1].text(0.5,0.5,"(2) LOCUS-LEVEL skipped\n(splice-counts file not built)",ha="center",va="center",fontsize=9,transform=ax[1].transAxes); ax[1].axis("off")
[ax[1].spines[s].set_visible(False) for s in ("top","right")]
fig.suptitle("Expression vs telotron presence (E. tenella; low power, N=21).  "+HEADLINE,fontsize=9.5,weight="bold",y=1.04)
fig.tight_layout(); fig.savefig("analysis/telotron_expr_presence.png",dpi=150,bbox_inches="tight"); print("\nwrote analysis/telotron_expr_presence.png")
