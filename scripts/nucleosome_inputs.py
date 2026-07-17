#!/usr/bin/env python3
"""Inputs for the telotron NuPoP nucleosome-occupancy analysis (replicating Gozashti et al. 2022,
PNAS, doi:10.1073/pnas.2209766119). ADAPTIVE flank for best accuracy with the data available: each
locus gets its maximum symmetric GAP-FREE flank up to --fmax (NuPoP's accuracy plateau); loci with
flank < --fmin are dropped (--fmin past NuPoP's edge-effect ramp). This includes the Tara MAG (short
contigs) that a fixed 5 kb flank would exclude. Per locus, emits two single-sequence FASTAs (NuPoP
needs one sequence per file):
  with/<lid>.fa    = 5'flank + telotron + 3'flank   (telotron at [flank, flank+ilen])
  without/<lid>.fa = 5'flank + 3'flank              (element removed; insertion site at position=flank)
oriented to the transcript strand (5' splice site at position=flank). Manifest carries species +
architecture (for faceting) + the per-locus flank."""
import csv, os, re, argparse, glob, sys, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import find_genome_fasta  # noqa: E402
def load(p):
    """Load a FASTA (plain or .gz) into a dict of seqid -> uppercase sequence."""
    opener = gzip.open if p.endswith(".gz") else open
    s={};h=None;b=[]
    for l in opener(p, "rt"):
        if l[0]=='>':
            if h:s[h]="".join(b)
            h=l[1:].split()[0];b=[]
        else:b.append(l.strip().upper())
    if h:s[h]="".join(b); return s
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def resolve(genome_id, refseq_dir, tara_dir):
    """Resolve GCF_/GCA_ + Tara through _common (handles .fna.gz + data/raw/genbank/)."""
    return find_genome_fasta(genome_id, refseq_dir, tara_dir, required=False)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--table", default="work/results/final_telotron_set_architecture.tsv")
    ap.add_argument("--refseq-dir", default="data/raw/refseq")
    ap.add_argument("--tara-dir", default="data/raw/tara")
    ap.add_argument("--out", default="work/results/nucleosome")
    ap.add_argument("--fmin", type=int, default=1000)
    ap.add_argument("--fmax", type=int, default=5000)
    # Genomes to analyse. Default = EVERY confident bearer species, read from
    # confident_species.tsv — not a hardcoded panel.
    #
    # This used to default to 6 literal IDs (5 Eimeria + the Tara MAG) and the
    # Snakefile passed no --genomes, so the default always won. That silently
    # falsified two claims: CLAUDE.md's "Every downstream analysis keys off
    # this file — new bearer species flow through automatically without
    # touching Python", and nucleosome_features' own "cluster count comes from
    # the actual manifest ... so new bearer species scale correctly". The
    # manifest could only ever contain those 6, so N_GENOMES_T was permanently
    # {Eimeria: 5, MAG: 1}, MIN_CLUSTERS=3 permanently suppressed the MAG
    # p-value, and admitting a new bearer (e.g. Monocercomonoides exilis)
    # changed the figure by exactly nothing.
    ap.add_argument("--confident-species", default="work/results/confident_species.tsv",
                    help="TSV from the confident_species rule; its genome_id column is "
                         "the default analysis set.")
    ap.add_argument("--genomes", nargs="*", default=None,
                    help="explicit genome_ids, overriding --confident-species "
                         "(use only for debugging a specific panel)")
    a=ap.parse_args()
    if not a.genomes:
        if not os.path.exists(a.confident_species):
            sys.exit(f"--confident-species {a.confident_species!r} not found and no "
                     f"--genomes given: refusing to fall back to a hardcoded panel. "
                     f"Run the confident_species rule first.")
        import csv as _csv
        with open(a.confident_species) as _fh:
            a.genomes = [r["genome_id"] for r in _csv.DictReader(_fh, delimiter="\t")]
        if not a.genomes:
            sys.exit(f"{a.confident_species}: no confident bearer species — nothing to analyse.")
        print(f"nucleosome inputs: {len(a.genomes)} confident bearer species "
              f"from {a.confident_species}", file=sys.stderr)
    FMAX,FMIN=a.fmax,a.fmin
    def gapfree(chrom,s,e):
        up=chrom[:s]; m=list(re.finditer("N{5,}",up)); left=s-(m[-1].end() if m else 0)
        dn=chrom[e:]; m=re.search("N{5,}",dn); right=(m.start() if m else len(dn))
        return min(left,right,FMAX)
    G={}
    for g in a.genomes:
        p=resolve(g,a.refseq_dir,a.tara_dir)
        if p: G[g]=load(p)
        else: print(f"WARN genome not found: {g}")
    for d in ("seqs/with","seqs/without"): os.makedirs(f"{a.out}/{d}",exist_ok=True)
    man=open(f"{a.out}/manifest.tsv","w"); man.write("locus_id\tgenome_id\torganism\tarchitecture\tstrand\tilen\tflank\n")
    from collections import Counter
    kept=0; spc=Counter(); arc=Counter()
    for r in csv.DictReader(open(a.table),delimiter="\t"):
        g=r["genome_id"]
        if g not in G or r["seqid"] not in G[g]: continue
        chrom=G[g][r["seqid"]]; s,e=int(r["start"])-1,int(r["end"]); ilen=e-s
        if ilen<10 or ilen>5000: continue
        flank=gapfree(chrom,s,e)
        if flank<FMIN: continue
        up=chrom[s-flank:s]; arr=chrom[s:e]; dn=chrom[e:e+flank]
        if "NNNNN" in up or "NNNNN" in dn: continue
        strand=r["strand"]
        if strand=="-": up,arr,dn=rc(dn),rc(arr),rc(up)
        lid=f"{g}__{r['seqid']}__{r['start']}_{r['end']}__{strand}"
        open(f"{a.out}/seqs/with/{lid}.fa","w").write(f">{lid}\n{up+arr+dn}\n")
        open(f"{a.out}/seqs/without/{lid}.fa","w").write(f">{lid}\n{up+dn}\n")
        man.write(f"{lid}\t{g}\t{r['organism']}\t{r['architecture']}\t{strand}\t{ilen}\t{flank}\n")
        kept+=1; spc[r['organism']]+=1; arc[r['architecture']]+=1
    man.close()
    print(f"FMIN={FMIN} FMAX={FMAX} ; kept {kept} telotron loci")
    print("by species:",dict(spc)); print("by architecture:",dict(arc))
if __name__=="__main__": main()
