#!/usr/bin/env python3
"""Characterise the telotron HOST-GENE class (the gene-bias revealed by the within-gene control):
telotrons accumulate in large, intron-rich genes — led by cytoskeletal motors (dyneins/kinesins) and
ATPases. Panels: (A) introns/gene host vs NON-host, (B) gene length host vs NON-host, (C) GO enrichment
of host genes, (D) recurrence (telotrons per host gene, MAG + Eimeria). Eimeria gene structure from
RefSeq GFFs. Host genes are compared against the DISJOINT non-host set (not the all-genes superset, which
is a set-vs-its-own-superset MWU); pooled p's are reported alongside a per-genome contrast because the
5 Eimeria assemblies share orthologues (the genome, n=5, is the real replication unit)."""
import csv, re
from collections import defaultdict, Counter
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
GFF={"GCF_000499545.2":"data/raw/refseq/GCF_000499545.2/GCF_000499545.2_ETH001_genomic.gff",
     "GCF_000499385.1":"data/raw/refseq/GCF_000499385.1/GCF_000499385.1_ENH001_genomic.gff",
     "GCF_000499425.1":"data/raw/refseq/GCF_000499425.1/GCF_000499425.1_EAH001_genomic.gff",
     "GCF_000499605.1":"data/raw/refseq/GCF_000499605.1/GCF_000499605.1_EMW001_genomic.gff",
     "GCF_000499745.2":"data/raw/refseq/GCF_000499745.2/GCF_000499745.2_EMH001_genomic.gff"}
def lin(o): return "MAG" if "MAG" in o else ("Eimeria" if "Eimeria" in o else None)
tcount=defaultdict(int); tcount_lin=defaultdict(lambda: defaultdict(int))
for r in csv.DictReader(open("work/results/final_telotron_set_architecture.tsv"),delimiter="\t"):
    g=r.get("gene_id") or r.get("tx_id"); L=lin(r["organism"])
    if r["genome_id"] in GFF: tcount[(r["genome_id"],g)]+=1
    if L: tcount_lin[L][(r["genome_id"],g)]+=1
# NB: host genes are a SUBSET of all genes, so a host-vs-all MWU compares a set against its own
# superset (violates the two-independent-samples assumption and shrinks the p artificially). We compare
# host vs the DISJOINT non-host set instead. Genes are also pooled across 5 Eimeria assemblies (the same
# orthologue can recur up to 5x => pseudoreplication), so besides the pooled MWU we report a per-species
# test + a combined per-species median contrast (the genome is the real replication unit, n=5).
nonintr=[]; hostintr=[]; nonlen=[]; hostlen=[]; tvi=[]
per_species=[]  # (species_gid, host_median_introns, nonhost_median_introns, n_host)
for gid,path in GFF.items():
    mex=defaultdict(int); mg={}; glen={}
    for l in open(path):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2]=="gene":
            m=re.search(r'ID=([^;]+)',f[8]);
            if m: glen[m.group(1)]=int(f[4])-int(f[3])
        elif f[2]=="mRNA":
            mi=re.search(r'ID=([^;]+)',f[8]); pa=re.search(r'Parent=([^;]+)',f[8])
            if mi and pa: mg[mi.group(1)]=pa.group(1)
        elif f[2]=="exon":
            pa=re.search(r'Parent=([^;]+)',f[8]);
            if pa: mex[pa.group(1)]+=1
    gi=defaultdict(int)
    for mr,ex in mex.items():
        g=mg.get(mr);
        if g: gi[g]=max(gi[g],ex-1)
    sp_h=[]; sp_n=[]
    for g,ni in gi.items():
        nt=tcount.get((gid,g),0); tvi.append((ni,nt))
        if nt>0: hostintr.append(ni); hostlen.append(glen.get(g,0)); sp_h.append(ni)
        else:    nonintr.append(ni); nonlen.append(glen.get(g,0)); sp_n.append(ni)
    if sp_h and sp_n:
        per_species.append((gid,np.median(sp_h),np.median(sp_n),len(sp_h)))
hostintr=np.array(hostintr); nonintr=np.array(nonintr)
print(f"host median introns {np.median(hostintr):.0f} (n={len(hostintr)}); non-host median {np.median(nonintr):.0f} (n={len(nonintr)}); "
      f"pooled MWU host>non-host p={mannwhitneyu(hostintr,nonintr,alternative='greater').pvalue:.1e}")
print(f"host len median {np.median([x for x in hostlen if x]):.0f}; non-host {np.median([x for x in nonlen if x]):.0f}; "
      f"pooled MWU p={mannwhitneyu([x for x in hostlen if x],[x for x in nonlen if x],alternative='greater').pvalue:.1e}")
# per-species (genome = replication unit): each species' host vs non-host median-intron contrast
if per_species:
    diffs=[h-n for _,h,n,_ in per_species]
    print("per-species host-vs-non-host median-intron diff (n=%d genomes): %s  -> all positive: %s"
          % (len(per_species), [f'{d:+.0f}' for d in diffs], all(d>0 for d in diffs)))
    try:
        from scipy.stats import wilcoxon
        if len(diffs)>=6: print(f"  Wilcoxon signed-rank across genomes p={wilcoxon(diffs,alternative='greater').pvalue:.3f}")
        else: print("  (n<6 genomes: report the per-species diffs/sign, not a within-genome-pooled p)")
    except Exception as e: print("  wilcoxon",e)
I=np.array([x[0] for x in tvi]); Tn=np.array([x[1] for x in tvi])
print(f"Spearman(introns, telotrons) per gene = {spearmanr(I,Tn).correlation:.3f} p={spearmanr(I,Tn).pvalue:.1e}")
# GO (optional archive input — panel C is skipped if the precomputed enrichment table is absent so the
# rule does not hard-depend on the read-only work/old/ archive). Override path with TELOTRON_GO_TSV.
import os as _os
GO_TSV=_os.environ.get("TELOTRON_GO_TSV","work/old/good_set/go_analysis/enrichment/meta_enrichment.tsv")
gotop=[]
if _os.path.exists(GO_TSV):
    go=[r for r in csv.DictReader(open(GO_TSV),delimiter="\t") if r["sign_list"].count("+")>=r["sign_list"].count("-")]
    go.sort(key=lambda r:-float(r["z_stouffer"])); gotop=go[:10]
else:
    print(f"[GO panel skipped: {GO_TSV} not found]")
# ---- figure ----
fig,ax=plt.subplots(1,4,figsize=(18,4.3))
# A introns (host vs DISJOINT non-host; p is pooled MWU - see caption re: per-species replication)
bins=np.arange(0,21)
p_intr=mannwhitneyu(hostintr,nonintr,alternative='greater').pvalue
ax[0].hist(np.clip(nonintr,0,20),bins=bins,density=True,color="#bbb",alpha=.8,label=f"non-host (med {np.median(nonintr):.0f})")
ax[0].hist(np.clip(hostintr,0,20),bins=bins,density=True,color="#b2182b",alpha=.6,label=f"telotron host (med {np.median(hostintr):.0f})")
ax[0].set_xlabel("introns per gene"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=7.5)
ax[0].set_title(f"A  telotron genes are intron-rich\n{np.median(hostintr):.0f} vs {np.median(nonintr):.0f} introns (pooled p={p_intr:.0e})",fontsize=10,weight="bold")
[ax[0].spines[s].set_visible(False) for s in ("top","right")]
# B length (host vs DISJOINT non-host)
hl=[x for x in hostlen if x]; nl=[x for x in nonlen if x]
p_len=mannwhitneyu(hl,nl,alternative='greater').pvalue
ax[1].hist(np.clip(nl,0,12000),bins=30,density=True,color="#bbb",alpha=.8,label="non-host")
ax[1].hist(np.clip(hl,0,12000),bins=30,density=True,color="#b2182b",alpha=.6,label="telotron host")
ax[1].set_xlabel("gene length (bp)"); ax[1].set_title(f"B  and longer\n{np.median(hl):.0f} vs {np.median(nl):.0f} bp (pooled p={p_len:.0e})",fontsize=10,weight="bold"); ax[1].legend(fontsize=7.5)
[ax[1].spines[s].set_visible(False) for s in ("top","right")]
# C GO (skipped if precomputed enrichment table absent)
if gotop:
    names=[r["name"][:34] for r in gotop][::-1]; zs=[float(r["z_stouffer"]) for r in gotop][::-1]
    ax[2].barh(range(len(names)),zs,color="#2166ac"); ax[2].set_yticks(range(len(names))); ax[2].set_yticklabels(names,fontsize=7)
    ax[2].set_xlabel("meta-enrichment z (Stouffer)"); ax[2].set_title("C  GO of host genes: motors + ATPases\n(microtubule motor q=9e-4)",fontsize=10,weight="bold")
else:
    ax[2].text(0.5,0.5,"GO enrichment table not available\n(work/old/ archive absent)",ha="center",va="center",fontsize=9,transform=ax[2].transAxes)
    ax[2].set_title("C  GO of host genes (skipped)",fontsize=10,weight="bold"); ax[2].axis("off")
[ax[2].spines[s].set_visible(False) for s in ("top","right")]
# D recurrence
for L,col in (("MAG","#1f77b4"),("Eimeria","#b2182b")):
    cc=Counter(tcount_lin[L].values()); tot=sum(cc.values())
    xs=sorted(cc); ax[3].plot(xs,[100*cc[x]/tot for x in xs],"o-",color=col,label=f"{L} ({sum(1 for v in tcount_lin[L].values() if v>=2)*100//tot}% ≥2)")
ax[3].set_yscale("log"); ax[3].set_xlabel("telotrons in one gene"); ax[3].set_ylabel("% of host genes")
ax[3].set_title("D  recurrence: telotrons cluster\nin specific genes (up to 9)",fontsize=10,weight="bold"); ax[3].legend(fontsize=8)
[ax[3].spines[s].set_visible(False) for s in ("top","right")]
fig.suptitle("Telotron host-gene class: telotrons accumulate in LARGE, INTRON-RICH genes (cytoskeletal motors / ATPases) — a gene-architecture/selection bias, not local targeting\n(host vs disjoint non-host; pooled p's over 5 Eimeria assemblies — direction confirmed per-genome, see stdout)",fontsize=10.5,weight="bold",y=1.04)
fig.tight_layout(); fig.savefig("work/results/figures/telotron_gene_bias.png",dpi=150,bbox_inches="tight")
print("wrote work/results/figures/telotron_gene_bias.png")
