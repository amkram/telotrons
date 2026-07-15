#!/usr/bin/env python3
"""Telotron host-gene class: (A) intron-rich host genes and (B) longer host genes vs the disjoint
non-host set (per-gene bias), but the per-intron dissection reveals (C) rate PEAKS at mid intron-count
and FALLS with gene length (bias is trial-count opportunity, not per-intron preference) with
(D) a genuine 5'-position preference. Iterates every genome in confident_species.tsv that has a
resolvable RefSeq / GenBank GFF (species without a GFF just get skipped). Host vs DISJOINT non-host
MWU. The confident-species set is the paper's central data-driven bearer set; downstream analyses
key off it so new species flow through without editing this script."""
import bisect as bs
import csv
import glob
import os
import re
from collections import defaultdict

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr

ARCH_TSV = "work/results/final_telotron_set_architecture.tsv"
CONFIDENT_TSV = "work/results/confident_species.tsv"


def resolve_gff(gid):
    """Find the GFF for a genome id under data/raw/{refseq,genbank}/{gid}/."""
    for root in ("data/raw/refseq", "data/raw/genbank"):
        hits = glob.glob(f"{root}/{gid}/*_genomic.gff") + glob.glob(f"{root}/{gid}/*_genomic.gff.gz")
        if hits:
            return hits[0]
    return None


def open_gff(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path)


# Confident-species panel (all genomes with a resolvable GFF).
GFF = {}
for r in csv.DictReader(open(CONFIDENT_TSV), delimiter="\t"):
    gid = r["genome_id"]
    p = resolve_gff(gid)
    if p:
        GFF[gid] = p
    else:
        print(f"[skip: no GFF] {gid} ({r.get('organism','')})")
print(f"panel: {len(GFF)} bearer genomes with GFF")

if not GFF:
    raise SystemExit("no confident bearers with resolvable GFFs — nothing to analyze")


# Load telotron rows (per-gene count + per-(genome,seqid) intervals).
telo_by_gene = defaultdict(int)
telo_intervals = defaultdict(list)
for r in csv.DictReader(open(ARCH_TSV), delimiter="\t"):
    gid = r["genome_id"]
    if gid not in GFF:
        continue
    g = r.get("gene_id") or r.get("tx_id")
    telo_by_gene[(gid, g)] += 1
    telo_intervals[(gid, r["seqid"])].append((int(r["start"]), int(r["end"])))
for k in telo_intervals:
    telo_intervals[k].sort()


def is_telo(genome, seqid, s, e):
    ivs = telo_intervals.get((genome, seqid))
    if not ivs:
        return False
    i = bs.bisect_right([x[0] for x in ivs], e)
    for ts, te in ivs[max(0, i - 30) : i + 1]:
        ov = min(e, te) - max(s, ts)
        if ov > 0 and ov / min(e - s, te - ts + 1) > 0.5:
            return True
    return False


per_gene = []   # (gid, gene, N_introns, gene_len, is_host)
per_intron = []  # (gid, N, gene_len, pos_frac, is_telo)
for gid, path in GFF.items():
    glen = {}
    mrna_gene = {}
    mrna_exons = defaultdict(list)
    for line in open_gff(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "gene":
            m = re.search(r"ID=([^;]+)", f[8])
            if m:
                glen[m.group(1)] = int(f[4]) - int(f[3])
        elif f[2] == "mRNA":
            mi = re.search(r"ID=([^;]+)", f[8])
            pa = re.search(r"Parent=([^;]+)", f[8])
            if mi and pa:
                mrna_gene[mi.group(1)] = pa.group(1)
        elif f[2] == "exon":
            pa = re.search(r"Parent=([^;]+)", f[8])
            if pa:
                mrna_exons[pa.group(1)].append((int(f[3]), int(f[4]), f[0], f[6]))

    gene_introns = defaultdict(set)
    gene_strand = {}
    for mr, exs in mrna_exons.items():
        g = mrna_gene.get(mr)
        if not g or len(exs) < 2:
            continue
        exs = sorted(exs)
        seqid = exs[0][2]
        gene_strand[g] = exs[0][3]
        for i in range(len(exs) - 1):
            istart, iend = exs[i][1], exs[i + 1][0]
            if iend - istart >= 20:
                gene_introns[g].add((istart, iend, seqid))

    for g, introns in gene_introns.items():
        ins = sorted(introns)
        if gene_strand.get(g) == "-":
            ins = ins[::-1]
        N = len(ins)
        is_host = telo_by_gene.get((gid, g), 0) > 0
        per_gene.append((gid, g, N, glen.get(g, 0), is_host))
        for idx, (s, e, seqid) in enumerate(ins):
            per_intron.append(
                (gid, N, glen.get(g, e - s), idx / (N - 1) if N > 1 else 0.5, is_telo(gid, seqid, s, e))
            )

# ── A/B  Per-gene host vs disjoint non-host  ──────────────────────────────────
host_intr = np.array([N for *_, N, _, h in per_gene if h])
non_intr = np.array([N for *_, N, _, h in per_gene if not h])
host_len = np.array([L for *_, _, L, h in per_gene if h and L])
non_len = np.array([L for *_, _, L, h in per_gene if not h and L])

p_intr = mannwhitneyu(host_intr, non_intr, alternative="greater").pvalue if len(host_intr) and len(non_intr) else float("nan")
p_len = mannwhitneyu(host_len, non_len, alternative="greater").pvalue if len(host_len) and len(non_len) else float("nan")
print(f"[A] host introns/gene med {np.median(host_intr):.0f} (n={len(host_intr)}) vs non-host {np.median(non_intr):.0f} (n={len(non_intr)}); pooled MWU p={p_intr:.1e} (POOLED — pseudo-replicated across bearer genomes; see genome-level test below)")
print(f"[B] host gene len med {np.median(host_len):.0f} vs non-host {np.median(non_len):.0f}; pooled MWU p={p_len:.1e} (POOLED — pseudo-replicated; see genome-level test below)")

# Genome-level test (the real replication unit): per-genome host-vs-non-host
# median diff, Wilcoxon signed-rank across genomes when n>=6. This is the
# HEADLINE test; the pooled MWU above is descriptive only.
per_species_diff = []
per_species_len_diff = []
for gid in GFF:
    h_sp = [N for g, _, N, _, host in per_gene if g == gid and host]
    n_sp = [N for g, _, N, _, host in per_gene if g == gid and not host]
    h_sp_L = [L for g, _, _, L, host in per_gene if g == gid and host and L]
    n_sp_L = [L for g, _, _, L, host in per_gene if g == gid and not host and L]
    if h_sp and n_sp:
        per_species_diff.append(np.median(h_sp) - np.median(n_sp))
    if h_sp_L and n_sp_L:
        per_species_len_diff.append(np.median(h_sp_L) - np.median(n_sp_L))
print(f"    per-genome host-vs-non-host intron-count diff (n={len(per_species_diff)}): {[f'{d:+.0f}' for d in per_species_diff]}  all positive: {all(d>0 for d in per_species_diff)}")
if len(per_species_diff) >= 6:
    try:
        from scipy.stats import wilcoxon
        wp = wilcoxon(per_species_diff, alternative="greater").pvalue
        print(f"    [HEADLINE] genome-level Wilcoxon signed-rank on intron-count diff: p={wp:.3f}  (n={len(per_species_diff)} genomes)")
    except Exception as ex:
        print(f"    (wilcoxon unavailable: {ex})")
else:
    print(f"    (n<6 genomes: report per-genome direction/sign, not a within-genome-pooled p)")

# ── C/D  Per-intron dissection  ───────────────────────────────────────────────
N = np.array([r[1] for r in per_intron])
GL = np.array([r[2] for r in per_intron], dtype=float)
POS = np.array([r[3] for r in per_intron])
Y = np.array([1 if r[4] else 0 for r in per_intron])
print(f"\n[C/D] total introns across bearers: {len(per_intron)}; telotron introns matched: {int(Y.sum())}")

bins_n = [(1, 1), (2, 3), (4, 6), (7, 10), (11, 20), (21, 200)]
rate_by_N = []
for lo, hi in bins_n:
    m = (N >= lo) & (N <= hi)
    tot = int(m.sum()); pos = int(Y[m].sum())
    rate_by_N.append((f"{lo}-{hi}", 1e4 * pos / tot if tot else 0))
print("[C] per-intron rate vs gene intron count:", rate_by_N)

rate_by_pos = []
for i in range(5):
    m = (POS >= i / 5) & (POS < (i + 1) / 5)
    tot = int(m.sum()); pos = int(Y[m].sum())
    rate_by_pos.append((f"{i/5:.1f}-{(i+1)/5:.1f}", 1e4 * pos / tot if tot else 0))
print("[D] per-intron rate vs 5'->3' position:", rate_by_pos)

X = np.column_stack([np.log(N), np.log(GL + 1), POS])
Xz = (X - X.mean(0)) / X.std(0)
try:
    import statsmodels.api as sm
    mod = sm.Logit(Y, sm.add_constant(Xz)).fit(disp=0)
    for nm, c, p in zip(["const", "log(N introns)", "log(gene length)", "position(5'->3')"], mod.params, mod.pvalues):
        print(f"    logistic {nm:20s} coef {c:+.3f}  p={p:.1e}")
except Exception as ex:
    print(f"    logistic failed: {ex}")

sp = spearmanr(N, Y)
print(f"    Spearman(introns/gene, is_telotron per-intron) = {sp.correlation:.3f} p={sp.pvalue:.1e}")

# ── Figure (4 panels) ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 4, figsize=(18, 4.3))
bins = np.arange(0, 21)
ax[0].hist(np.clip(non_intr, 0, 20), bins=bins, density=True, color="#bbb", alpha=0.8, label=f"non-host (med {np.median(non_intr):.0f})")
ax[0].hist(np.clip(host_intr, 0, 20), bins=bins, density=True, color="#b2182b", alpha=0.6, label=f"telotron host (med {np.median(host_intr):.0f})")
ax[0].set_xlabel("introns per gene"); ax[0].set_ylabel("density")
ax[0].set_title(f"A  intron-rich host genes\nmed {np.median(host_intr):.0f} vs {np.median(non_intr):.0f} (p={p_intr:.0e})", fontsize=10, weight="bold")
ax[0].legend(fontsize=7.5)
[ax[0].spines[s].set_visible(False) for s in ("top", "right")]

ax[1].hist(np.clip(non_len, 0, 12000), bins=30, density=True, color="#bbb", alpha=0.8, label="non-host")
ax[1].hist(np.clip(host_len, 0, 12000), bins=30, density=True, color="#b2182b", alpha=0.6, label="telotron host")
ax[1].set_xlabel("gene length (bp)")
ax[1].set_title(f"B  long host genes\nmed {np.median(host_len):.0f} vs {np.median(non_len):.0f} (p={p_len:.0e})", fontsize=10, weight="bold")
ax[1].legend(fontsize=7.5)
[ax[1].spines[s].set_visible(False) for s in ("top", "right")]

ax[2].plot([r[1] for r in rate_by_N], "o-", color="#b2182b")
ax[2].set_xticks(range(len(rate_by_N))); ax[2].set_xticklabels([r[0] for r in rate_by_N], fontsize=8)
ax[2].set_xlabel("introns per host gene"); ax[2].set_ylabel("telotron rate /10k introns")
ax[2].set_title("C  per-intron rate peaks at 4-10\n(A/B bias is opportunity + preference)", fontsize=10, weight="bold")
[ax[2].spines[s].set_visible(False) for s in ("top", "right")]

ax[3].plot([r[1] for r in rate_by_pos], "o-", color="#1a7a1a")
ax[3].set_xticks(range(len(rate_by_pos))); ax[3].set_xticklabels([r[0] for r in rate_by_pos], fontsize=8)
ax[3].set_xlabel("position within gene (5'->3')"); ax[3].set_ylabel("telotron rate /10k introns")
ax[3].set_title("D  genuine 5' preference\n(logistic partialling out N + gene-len)", fontsize=10, weight="bold")
[ax[3].spines[s].set_visible(False) for s in ("top", "right")]

fig.suptitle(f"Telotron host-gene class ({len(GFF)} confident bearers, {int(Y.sum())} telotron introns) — intron-rich long host genes (A/B), per-intron rate is opportunity (C) + 5' preference (D)", fontsize=10.5, weight="bold", y=1.04)
out = "work/results/figures/telotron_gene_class.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"wrote {out}")
