#!/usr/bin/env python3
"""Matched non-telotron-intron CONTROL for the telotron NuPoP analysis. Samples confident
non-telotron introns (telomeric_frac<0.1) from the SAME genomes, builds the identical adaptive-flank
with/without NuPoP inputs, so we can ask whether the telotron occupancy pattern is telotron-specific
or generic to introns. Control count per species is matched to the telotron count (clamped [50,400])."""
import csv, os, re, random, argparse, glob
from collections import Counter, defaultdict
def load(p):
    s={};h=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if h:s[h]="".join(b)
            h=l[1:].split()[0];b=[]
        else:b.append(l.strip().upper())
    if h:s[h]="".join(b); return s
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def resolve(gid, refseq_dir, tara_dir):
    if gid.startswith("GCF_") or gid.startswith("GCA_"):
        hit=glob.glob(f"{refseq_dir}/{gid}/{gid}_*_genomic.fna"); return hit[0] if hit else None
    p=f"{tara_dir}/Contigs/{gid}.fa"; return p if os.path.exists(p) else None
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--controls", default="work/results/non_telotron_controls.tsv")
    ap.add_argument("--telo-manifest", default="work/results/nucleosome/manifest.tsv")
    ap.add_argument("--refseq-dir", default="data/raw/refseq"); ap.add_argument("--tara-dir", default="data/raw/tara")
    ap.add_argument("--out", default="work/results/nucleosome")
    ap.add_argument("--fmin", type=int, default=1000); ap.add_argument("--fmax", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    a=ap.parse_args(); random.seed(a.seed); FMAX,FMIN=a.fmax,a.fmin
    def gapfree(chrom,s,e):
        up=chrom[:s]; m=list(re.finditer("N{5,}",up)); left=s-(m[-1].end() if m else 0)
        dn=chrom[e:]; m=re.search("N{5,}",dn); right=(m.start() if m else len(dn))
        return min(left,right,FMAX)
    telo=Counter(r["organism"] for r in csv.DictReader(open(a.telo_manifest),delimiter="\t"))
    target={sp:max(50,min(400,n)) for sp,n in telo.items()}
    rows=list(csv.DictReader(open(a.controls),delimiter="\t"))
    gids=set(r["genome_id"] for r in rows if r["organism"] in target)
    G={}
    for g in gids:
        p=resolve(g,a.refseq_dir,a.tara_dir)
        if p: G[g]=load(p)
    bysp=defaultdict(list)
    for r in rows:
        if r["genome_id"] in G and r["organism"] in target: bysp[r["organism"]].append(r)
    for d in ("control_seqs/with","control_seqs/without"): os.makedirs(f"{a.out}/{d}",exist_ok=True)
    man=open(f"{a.out}/control_manifest.tsv","w"); man.write("locus_id\tgenome_id\torganism\tarchitecture\tstrand\tilen\tflank\n")
    kept=Counter()
    for sp,want in target.items():
        cand=bysp.get(sp,[]); random.shuffle(cand)
        for r in cand:
            if kept[sp]>=want: break
            try:
                if float(r.get("telomeric_frac",0))>=0.1: continue
            except: pass
            g=r["genome_id"]
            if r["seqid"] not in G[g]: continue
            chrom=G[g][r["seqid"]]; s,e=int(r["start"])-1,int(r["end"]); ilen=e-s
            if ilen<10 or ilen>5000: continue
            flank=gapfree(chrom,s,e)
            if flank<FMIN: continue
            up=chrom[s-flank:s]; arr=chrom[s:e]; dn=chrom[e:e+flank]
            if "NNNNN" in up or "NNNNN" in dn: continue
            strand=r["strand"]
            if strand=="-": up,arr,dn=rc(dn),rc(arr),rc(up)
            lid=f"CTRL__{g}__{r['seqid']}__{r['start']}_{r['end']}__{strand}"
            open(f"{a.out}/control_seqs/with/{lid}.fa","w").write(f">{lid}\n{up+arr+dn}\n")
            open(f"{a.out}/control_seqs/without/{lid}.fa","w").write(f">{lid}\n{up+dn}\n")
            man.write(f"{lid}\t{g}\t{r['organism']}\tcontrol\t{strand}\t{ilen}\t{flank}\n")
            kept[sp]+=1
    man.close(); print("control loci kept per species:",dict(kept),"total",sum(kept.values()))
if __name__=="__main__": main()
