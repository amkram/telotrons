#!/usr/bin/env python3
"""Decisive test: is the telotron insertion-site feature signal a HOST-GENE-CLASS bias (selection/
detection) or LOCAL insertion targeting? For each verified-real feature (GC, CpG O/E, 10-bp WW
periodicity) we compute it on the exon flanks of three intron sets, per lineage:
  (T) telotron introns
  (S) SAME-GENE sibling introns (other, non-telomeric introns of the telotron's own gene) -> controls
      gene identity AND expression
  (R) RANDOM introns from non-telotron genes
Logic: telotron-vs-R is the original effect. If T ~= S (paired, within gene) while S != R, the signal
is a whole-gene property (gene-class). If T != S, the telotron site is special beyond its gene = local."""
import csv, os, re, math, random, glob, statistics as st
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon
random.seed(0); W=100
def load_g(p):
    s={};h=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if h:s[h]="".join(b)
            h=l[1:].split()[0];b=[]
        else:b.append(l.strip().upper())
    if h:s[h]="".join(b); return s
def resolve(g,rs="data/raw/refseq",ta="data/raw/tara"):
    if g.startswith(("GCF_","GCA_")):
        h=glob.glob(f"{rs}/{g}/{g}_*_genomic.fna"); return h[0] if h else None
    p=f"{ta}/Contigs/{g}.fa"; return p if os.path.exists(p) else None
def f_gc(s):
    n=sum(s.count(b) for b in "ACGT") or 1; return (s.count("G")+s.count("C"))/n
def f_cpg(s):
    n=len(s); c=s.count("C"); g=s.count("G"); return (s.count("CG"))/((c*g)/n) if c and g and n>1 else np.nan
def f_ww(s):
    if len(s)<60: return np.nan
    w=np.array([1.0 if (s[i] in "AT" and s[i+1] in "AT") else 0.0 for i in range(len(s)-1)]); w=w-w.mean()
    if w.std()==0: return np.nan
    ac=lambda L:np.dot(w[:-L],w[L:])/len(w[:-L]); return ac(10)-np.mean([ac(L) for L in (6,7,8,13,14)])
FEATS=[("GC",f_gc),("CpG O/E",f_cpg),("10-bp WW periodicity",f_ww)]
GENOMES=["GCF_000499545.2","GCF_000499385.1","GCF_000499745.2","GCF_000499425.1","GCF_000499605.1","TARA_PSW_86_MAG_00284"]
G={g:load_g(resolve(g)) for g in GENOMES if resolve(g)}
def lineage(o): return "MAG" if "MAG" in o else ("Eimeria" if "Eimeria" in o else None)
# telotron intron keys (genome, seqid, start, end) + host gene
telo=set(); telo_gene=defaultdict(list); arch_org={}
for r in csv.DictReader(open("work/results/final_telotron_set_architecture.tsv"),delimiter="\t"):
    if r["genome_id"] not in G: continue
    key=(r["genome_id"],r["seqid"],int(r["start"]),int(r["end"])); telo.add(key)
    gid=r.get("gene_id") or r.get("tx_id")
    telo_gene[(r["genome_id"],gid)].append((key,r["organism"]))
# all introns grouped by gene; classify telotron / sibling(non-telo) ; collect random pool
by_gene=defaultdict(list); telo_genes=set(telo_gene.keys())
def flankfeat(genome,seqid,s,e,fn):
    chrom=G[genome].get(seqid)
    if not chrom or s-W<0 or e+W>len(chrom): return np.nan
    up=chrom[s-W:s]; dn=chrom[e:e+W]
    if "NNNNN" in up or "NNNNN" in dn: return np.nan
    return np.nanmean([fn(up),fn(dn)])
rand_pool=defaultdict(list)
for r in csv.DictReader(open("work/results/non_telotron_controls.tsv"),delimiter="\t"):  # confident non-telo introns, has gene_id
    if r["genome_id"] not in G: continue
    lin=lineage(r["organism"]);
    if lin is None: continue
    try: tf=float(r.get("telomeric_frac",0))
    except: tf=0
    s,e=int(r["start"]),int(r["end"])
    if e-s<30 or e-s>5000: continue
    gid=r.get("gene_id") or r.get("tx_id"); gk=(r["genome_id"],gid)
    rec=(r["genome_id"],r["seqid"],s,e,lin,tf)
    by_gene[gk].append(rec)
    if tf<0.1 and gk not in telo_genes: rand_pool[lin].append(rec)
# build the three sets per feature/lineage
print(f"telotron genes: {len(telo_genes)}  loaded genomes: {list(G)}")
res={}  # (feat,lin)-> dict of arrays
for fname,fn in FEATS:
    cache={}
    def feat(rec):
        k=rec[:4]
        if k not in cache: cache[k]=flankfeat(rec[0],rec[1],rec[2],rec[3],fn)
        return cache[k]
    pairs=defaultdict(list)  # lineage -> (telo_val, sibling_mean) paired within gene
    Tvals=defaultdict(list); Svals=defaultdict(list)
    for gk,(keys) in telo_gene.items():
        sibs=[rec for rec in by_gene.get(gk,[]) if rec[5]<0.1 and (rec[0],rec[1],rec[2],rec[3]) not in telo]
        sv=[feat(r) for r in sibs]; sv=[x for x in sv if not np.isnan(x)]
        for key,org in keys:
            lin=lineage(org); tv=feat((key[0],key[1],key[2],key[3]))
            if np.isnan(tv): continue
            Tvals[lin].append(tv)
            if sv:
                Svals[lin].append(np.mean(sv)); pairs[lin].append((tv,np.mean(sv)))
    Rvals={lin:[feat(r) for r in random.sample(rand_pool[lin],min(800,len(rand_pool[lin])))] for lin in ("MAG","Eimeria")}
    Rvals={lin:[x for x in v if not np.isnan(x)] for lin,v in Rvals.items()}
    print(f"\n=== {fname} ===")
    for lin in ("MAG","Eimeria"):
        T=np.array(Tvals[lin]); R=np.array(Rvals[lin]); pr=pairs[lin]
        if len(T)<5 or len(R)<5: continue
        tvR=mannwhitneyu(T,R).pvalue
        msg=f"  {lin:8s} T med {np.median(T):+.4f} (n={len(T)}) | sibling med {np.median([s for _,s in pr]) if pr else float('nan'):+.4f} | random med {np.median(R):+.4f} (n={len(R)})"
        print(msg)
        print(f"           T vs RANDOM: telo-rand {np.median(T)-np.median(R):+.4f} p={tvR:.1e}")
        if len(pr)>=10:
            d=[t-s for t,s in pr]; wp=wilcoxon(d).pvalue if any(d) else 1.0
            S=np.array([s for _,s in pr])
            svR=mannwhitneyu(S,R).pvalue
            print(f"           SIBLING vs RANDOM (is the whole gene special?): sib-rand {np.median(S)-np.median(R):+.4f} p={svR:.1e}")
            print(f"           T vs SIBLING (within-gene; local?): median paired diff {np.median(d):+.4f} p={wp:.1e}  [~0 => gene-class; !=0 => local]")
