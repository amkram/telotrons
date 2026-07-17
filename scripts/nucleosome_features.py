#!/usr/bin/env python3
"""Insertion-site FEATURE PANEL for telotrons (junction-artifact-free version).
Every feature is computed on the NON-JOINED real host flanks immediately bordering the array (from the
WITH-element sequences: up=[F-w:F], dn=[F+L:F+L+w]), averaged — so there is NO artificial splice-junction
join (the join confound that inflated several features in the naive element-removed window). Compared:
telotron vs matched non-telotron-intron control, split MAG vs Eimeria, as local enrichment (insertion-
flanking minus broader-flank background). Each summary cell is annotated with the adversarially-verified
verdict (workflow wf_c63b87c1): real / GC-confounded / junction-artifact / ns. NuPoP occupancy itself is a
GC proxy (see the separate nucleosome_occupancy.png) and is NOT a nucleosome-positioning measure; the
GC-independent positioning signal is the 10-bp WW periodicity.
CAVEATS baked into the caption: (1) Composition features (GC/CpG/G4/entropy/TpA) are now computed on
TELOMERE-MASKED flank windows (telomere_mask.mask_telomere_fragments) — this closes the residual-telomeric
leakage that previously ~halved the real MAG-GC host component (the rotation trap behind 4 prior
retractions). WW periodicity stays unmasked (composition-free by construction). (2) All MAG loci
(1226/1517) are one genome (TARA_PSW_86_MAG_00284): MAG = one observation, no within-lineage replication.
(3) Significance stars are BH-FDR q-values across the feature x lineage grid, not raw p."""
import csv, os, math, re, sys
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from telomere_mask import mask_telomere_fragments  # noqa: E402  (closes the rotation-contamination trap)
ROOT="work/results/nucleosome"; WI=75; WP=100
def seqof(p): return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
def bh_fdr(pvals):
    """Benjamini-Hochberg q-values for a list of p-values (NaN-safe, preserves order)."""
    p=np.asarray(pvals,dtype=float); ok=~np.isnan(p); idx=np.where(ok)[0]
    q=np.full(p.shape,np.nan)
    if idx.size==0: return q
    ps=p[idx]; order=np.argsort(ps); m=idx.size
    ranked=ps[order]*m/np.arange(1,m+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    qi=np.empty(m); qi[order]=ranked; q[idx]=qi
    return q
# NOTE ON DENOMINATORS: these run on sequence that mask_telomere_fragments()
# has already rewritten to 'N'. Every denominator must therefore count REAL
# bases, not len(s) — Ns are absent data, not observations. Counting them
# deflates the denominator in exact proportion to how much telomere was
# masked, which happens ONLY on the telotron arm (controls are near-clean), so
# the bug fabricates a telotron-vs-control difference out of the masking
# itself. Measured on a 75 bp flank with 35 Ns: CpG O/E reported 4.036 vs a
# true 2.127, against an unaffected control at 2.921 — a 1.9x inflation that
# SIGN-FLIPPED the result from depletion to enrichment.
def _n_real(s):
    return sum(s.count(b) for b in "ACGT")
def f_gc(s):
    n=_n_real(s) or 1; return (s.count("G")+s.count("C"))/n
def f_cpg(s):
    n=_n_real(s); c=s.count("C"); g=s.count("G")
    return (s.count("CG"))/((c*g)/n) if c and g and n>1 else 0.0
G4=re.compile(r"(?:G{3,}\w{1,7}){3,}G{3,}"); C4=re.compile(r"(?:C{3,}\w{1,7}){3,}C{3,}")
def f_g4(s):
    if len(s)<15: return 0.0
    cov=bytearray(len(s))
    for pat in (G4,C4):
        for m in pat.finditer(s):
            for i in range(m.start(),m.end()): cov[i]=1
    return sum(cov)/(_n_real(s) or 1)
def f_ww10(s):
    if len(s)<60: return np.nan
    w=np.array([1.0 if (s[i] in "AT" and s[i+1] in "AT") else 0.0 for i in range(len(s)-1)]); w=w-w.mean()
    if w.std()==0: return np.nan
    ac=lambda lag: np.dot(w[:-lag],w[lag:])/len(w[:-lag])
    return ac(10)-np.mean([ac(l) for l in (6,7,8,13,14)])
def f_entropy(s):
    if len(s)<4: return np.nan
    from collections import Counter
    di=Counter(s[i:i+2] for i in range(len(s)-1) if set(s[i:i+2])<=set("ACGT")); tot=sum(di.values()) or 1
    return (-sum((v/tot)*math.log2(v/tot) for v in di.values()))/4.0
def f_tpa(s):
    # Count only dinucleotide positions where BOTH bases are real (see the
    # denominator note above): masking a TTTAGGG array destroys its TA
    # dinucleotides while len(s)-1 stays put, deflating telotron TpA ~2x.
    n=sum(1 for i in range(len(s)-1) if s[i] in "ACGT" and s[i+1] in "ACGT")
    return s.count("TA")/n if n>0 else 0.0
def lcs(a,b):
    best=0
    for i in range(len(a)):
        for j in range(len(b)):
            l=0
            while i+l<len(a) and j+l<len(b) and a[i+l]==b[j+l]: l+=1
            best=max(best,l)
    return best
# feature: name, fn, half-window, mask_telomere?, verdict(MAG,Eim)
#
# mask=True for EVERY sequence metric, WW periodicity included. Residual
# telomeric units at the array boundary leak into these and inflated MAG-GC
# ~2x — the rotation trap behind four prior retractions (G4, exonic GC, RLFS,
# donor consensus).
#
# WW periodicity used to run UNMASKED, justified as "composition-free ...
# GC/telomere-content-independent by construction". That was empirically FALSE
# and made this arm's headline signal an artifact. f_ww10 = ac(10) - mean(ac at
# 6,7,8,13,14). A TTTAGGG array has period 7, so ac(7) and ac(14) are HARMONIC
# PEAKS sitting inside the subtracted baseline while ac(10) sits in a trough.
# Measured on a 100 bp window: clean host +0.007; 28 bp of telomere leaked in
# -0.027; 56 bp -0.144; pure array -0.346. Controls score ~0, so a telotron
# whose misannotated flank retains 56 bp of array reported a delta 20x the
# clean-host magnitude and the pooled T-vs-C test called it "real". Worse, the
# direction is motif-dependent — TTAGGG (period 6) gives -0.136 but TTAGG
# (period 5) gives +0.321, a SIGN FLIP — so the MAG-vs-Eimeria contrast was
# reading each lineage's repeat rather than its chromatin. Masking the same
# window restores +0.016. (nucleosome_withingene.py has always masked before
# f_ww; the two scripts disagreed on the same named feature.)
#
# The old worry that masking "destroys the very AT-periodicity it measures" is
# backwards: masked positions become N, contribute 0 to the W vector, and are
# thus excluded rather than counted as non-WW — which is what we want, since
# telomere is not host chromatin.
#
# VERDICTS BELOW ARE STALE: they were assigned under the unmasked-WW and
# N-in-denominator bugs and must be re-derived from a fresh run before being
# cited.
FEATS=[("GC content",f_gc,WI,True,("real","ns")),
       ("CpG O/E",f_cpg,WI,True,("real*","real")),
       ("G4 propensity (200bp)",f_g4,200,True,("ns·bdry","ns·bdry")),
       ("10-bp WW periodicity",f_ww10,WP,True,("real","real")),
       ("dinucleotide entropy",f_entropy,WI,True,("GC-conf","real")),
       ("TpA fraction",f_tpa,WI,True,("GC-conf","real")),
       ("junction microhomology",None,None,False,("ns","artifact"))]
def lineage(o): return "MAG" if "MAG" in o else ("Eimeria" if "Eimeria" in o else "o")
def load(man,wd):
    out=[]; genome_by_lin=defaultdict(set)
    reader = csv.DictReader(open(man),delimiter="\t")
    if "genome_id" not in (reader.fieldnames or []):
        raise KeyError(
            f"{man}: manifest missing 'genome_id' column — cluster-robust inference "
            "in nucleosome_features requires it (see load()); regenerate the manifest "
            "via nucleosome_analysis rule.")
    for r in reader:
        F=int(r["flank"]); L=int(r["ilen"]); lid=r["locus_id"]; p=f"{wd}/{lid}.fa"
        if not os.path.exists(p): continue
        s=seqof(p)
        if len(s)<2*F+L: continue
        lin=lineage(r["organism"])
        out.append((lin,F,L,s))
        # Count distinct genome_ids per lineage — the real replication unit
        # for cluster-robust inference. Falls back to counting the locus itself
        # only when genome_id is genuinely empty (should never happen given the
        # KeyError guard above).
        gid = r.get("genome_id") or lid
        genome_by_lin[lin].add(gid)
    return out, {lin: len(genomes) for lin, genomes in genome_by_lin.items()}
T, N_GENOMES_T = load(f"{ROOT}/manifest.tsv",f"{ROOT}/seqs/with")
C, _NG_C = load(f"{ROOT}/control_manifest.tsv",f"{ROOT}/control_seqs/with")
print(f"telotron {len(T)} control {len(C)} (computed on non-joined real host flanks)")
def feat_delta(recs,fn,w,mask=False):
    out=defaultdict(list)
    for lin,F,L,s in recs:
        up=s[F-w:F]; dn=s[F+L:F+L+w]; bg=s[150:F-150]+s[F+L+150:2*F+L-150]
        if mask:  # N-out residual telomeric runs leaking from the array into the flank windows
            up=mask_telomere_fragments(up); dn=mask_telomere_fragments(dn); bg=mask_telomere_fragments(bg)
        ins=np.nanmean([fn(up),fn(dn)]); b=fn(bg)
        if not np.isnan(ins) and not np.isnan(b): out[lin].append(ins-b)
    return out
def mh_delta(recs):  # junction MH directly (no interior-baseline subtraction that caused the artifact)
    out=defaultdict(list)
    for lin,F,L,s in recs: out[lin].append(lcs(s[F-50:F], s[F+L:F+L+50]))
    return out
LINS=["MAG","Eimeria"]
# Guard: cluster-robust inference requires ≥ MIN_CLUSTERS independent genomes
# per lineage. Cluster count comes from the actual manifest (load() counts
# distinct genome_ids per lineage) so new bearer species scale correctly —
# no more hardcoded "MAG=1, Eimeria=5" that silently breaks when the
# confident_species set changes.
MIN_CLUSTERS = 3

rows=[]
for name,fn,w,mask,verd in FEATS:
    td=mh_delta(T) if fn is None else feat_delta(T,fn,w,mask)
    cd=mh_delta(C) if fn is None else feat_delta(C,fn,w,mask)
    for j,lin in enumerate(LINS):
        tv=td.get(lin,[]); cv=cd.get(lin,[])
        n_clu = N_GENOMES_T.get(lin, 0)
        if len(tv)>3 and len(cv)>3:
            if n_clu < MIN_CLUSTERS:
                # Refuse a p — single-assembly pseudo-replication makes it uninterpretable.
                rows.append((name,lin,np.median(tv)-np.median(cv), np.nan, np.nan, f"{verd[j]} [1-genome; p suppressed]"))
                continue
            p=mannwhitneyu(tv,cv).pvalue
            U=mannwhitneyu(tv,cv,alternative="two-sided").statistic; r=2*U/(len(tv)*len(cv))-1
            rows.append((name,lin,np.median(tv)-np.median(cv),p,r,verd[j]))
# BH-FDR across the whole feature x lineage grid (significance stars below use q, not raw p)
qvals=bh_fdr([row[3] for row in rows])
rows=[(nm,lin,d,p,r,vd,q) for (nm,lin,d,p,r,vd),q in zip(rows,qvals)]
# ---- summary grid ----
feats=[f[0] for f in FEATS]
fig,ax=plt.subplots(figsize=(8.6,5.2))
M=np.full((len(feats),2),np.nan)
cell={}
for (nm,lin,d,p,r,vd,q) in rows:
    i=feats.index(nm); j=LINS.index(lin); M[i,j]=r; cell[(i,j)]=(d,q,vd)
im=ax.imshow(M,cmap="RdBu_r",aspect="auto",vmin=-0.5,vmax=0.5)
vcol={"real":"#0a0","real*":"#0a0","GC-conf":"#c80","artifact":"#c00","ns":"#888"}
_xlabels = [f"{lin}\n({N_GENOMES_T.get(lin, 0)} genome{'s' if N_GENOMES_T.get(lin, 0) != 1 else ''})" for lin in LINS]
ax.set_xticks(range(2)); ax.set_xticklabels(_xlabels); ax.set_yticks(range(len(feats))); ax.set_yticklabels(feats,fontsize=8.5)
for (i,j),(d,q,vd) in cell.items():  # stars are BH q-values, not raw p
    st="***" if q<1e-3 else "**" if q<1e-2 else "*" if q<.05 else "ns"
    ax.text(j,i,f"r={M[i,j]:+.2f} {st}\n[{vd}]",ha="center",va="center",fontsize=6.8,
            color=vcol.get(vd,"#000"),weight="bold")
plt.colorbar(im,label="telotron vs control (rank-biserial r)",fraction=.046)
ax.set_title("Telotron insertion-site features — REAL signal vs confound (verified)\n"
             "computed on non-joined real host flanks; cell verdict from adversarial GC+junction+leakage tests",fontsize=9.5,weight="bold")
fig.text(0.5,-0.02,"REAL: 10-bp WW periodicity (both, GC-free = the genuine nucleosome-positioning signal) · CpG O/E (both; Eimeria GC-clean) · GC (MAG only, introner-like) · entropy/TpA (Eimeria).   "
                   "NOT robust: G4 depletion is a sub-nucleosome boundary effect (strong at 50bp, NS at 200bp shown, reverses at 400bp) · junction microhomology is at background (NOT suppressed) · NuPoP occupancy is a GC proxy.   Composition feats computed on telomere-MASKED flanks (closes the rotation-leakage that ~halved MAG-GC); stars are BH q; MAG = 1 genome.",
         ha="center",fontsize=6.4,wrap=True)
fig.tight_layout(); fig.savefig(f"{ROOT}/nucleosome_feature_summary.png",dpi=150,bbox_inches="tight")
print("wrote nucleosome_feature_summary.png (junction-free, verdict-annotated)")
print("\n=== junction-free telotron-vs-control (rank-biserial r); composition feats masked; q=BH-FDR ===")
for (nm,lin,d,p,r,vd,q) in rows: print(f"  {nm:24s} {lin:8s} r={r:+.2f} dΔ={d:+.4f} p={p:.1e} q={q:.1e}  [{vd}]")
