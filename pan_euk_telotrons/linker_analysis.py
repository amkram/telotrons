#!/usr/bin/env python3
"""
Linker analysis for converging telotrons.

Tests:
  A. Self-similarity: do linkers cluster into families? (TERC-retrotransposition prediction)
  B. Genomic origin: do linker kmers match the host genome's exonic/intergenic kmers?
  C. Length distribution within length classes: distinguishes mechanisms
  D. Composition: GC%, complexity, dinucleotide skew vs random intron sequence
  E. Position relative to flanking telomeric arrays (asymmetry)
  F. Reverse complement check: does linker match itself in opposite orientation?
     (Nergadze's specific prediction: linker is opposite-orientation to flanks)

Inputs:  real_telotrons/converging_linkers.fa
         real_telotrons/_tara_oceans_telotrons.tsv (for non-telo controls)
"""

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

# ------------------------------------------------------------------
# Parsing
# ------------------------------------------------------------------

def load_linker_fasta(path):
    """Returns list of (id, seq) tuples."""
    records = []
    curr_id = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_id:
                    records.append((curr_id, ''.join(parts).upper()))
                curr_id = line[1:]
                parts = []
            else:
                parts.append(line)
    if curr_id:
        records.append((curr_id, ''.join(parts).upper()))
    return records


def kmers(seq, k=8):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def rc(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


# ------------------------------------------------------------------
# Self-similarity (Nergadze-style: linkers from different loci share template)
# ------------------------------------------------------------------

def kmer_jaccard_clusters(linkers, k=10, min_jaccard=0.20):
    """Greedy clustering: linkers sharing ≥min_jaccard kmers join the same cluster.
    Returns list of cluster sizes (with members)."""
    # Build kmer sets per linker
    sets = []
    for sid, seq in linkers:
        if len(seq) < k * 2:
            sets.append((sid, seq, set()))
            continue
        ks = set(kmers(seq, k))
        # Also include rc kmers
        ks |= set(kmers(rc(seq), k))
        sets.append((sid, seq, ks))

    # Sort by size descending — largest seeds clusters
    sets.sort(key=lambda x: -len(x[2]))

    # Greedy assignment
    clusters = []  # each is list of (sid, seq, kmer_set)
    cluster_consensus = []  # union of kmers per cluster

    for sid, seq, ks in sets:
        if not ks:
            continue
        best_idx = -1; best_jaccard = 0
        for i, cks in enumerate(cluster_consensus):
            inter = len(ks & cks)
            if inter == 0: continue
            j = inter / min(len(ks), len(cks))   # asymmetric — fits if linker IS the template
            if j > best_jaccard:
                best_jaccard = j
                best_idx = i
        if best_idx >= 0 and best_jaccard >= min_jaccard:
            clusters[best_idx].append((sid, seq))
            cluster_consensus[best_idx] |= ks
        else:
            clusters.append([(sid, seq)])
            cluster_consensus.append(set(ks))

    # Sort clusters by size
    clusters_sorted = sorted(zip(clusters, cluster_consensus),
                             key=lambda x: -len(x[0]))
    return clusters_sorted


def linker_global_kmer_freq(linkers, k=8):
    """Find the most abundant kmers across all linkers — these are template-like consensus."""
    kmer_counts = Counter()
    n_linkers_per_kmer = Counter()
    for sid, seq in linkers:
        seen_in_this = set()
        for km in kmers(seq, k):
            kmer_counts[km] += 1
            seen_in_this.add(km)
        for km in seen_in_this:
            n_linkers_per_kmer[km] += 1
    return kmer_counts, n_linkers_per_kmer


# ------------------------------------------------------------------
# Composition
# ------------------------------------------------------------------

def composition(seq):
    n = len(seq)
    if n == 0:
        return None
    c = Counter(seq)
    gc = (c.get('G', 0) + c.get('C', 0)) / n
    at = (c.get('A', 0) + c.get('T', 0)) / n
    # Dinucleotide
    dinuc = Counter(seq[i:i+2] for i in range(n - 1))
    cpg = dinuc.get('CG', 0) / max(1, n - 1)
    return {
        'gc': gc,
        'at': at,
        'cpg_freq': cpg,
        'tpa_freq': dinuc.get('TA', 0) / max(1, n - 1),
        'len': n,
    }


# ------------------------------------------------------------------
# Reverse-complement check (Nergadze)
# ------------------------------------------------------------------

def linker_rc_relative_to_flanks(linker_id, linker_seq):
    """Parse the linker_id to find which class. Reverse-complement test:
    if linker is opposite-orientation to flanks, the linker should NOT
    contain TTAGGG-strand telomeric kmers but might contain non-telo
    sequence that, when reverse-complemented, looks like part of TERC."""
    # In CONVERGING_WITHLINKER, flanks are TTAGGG (5') and CCCTAA (3')
    # If linker comes from TERC retrotranscription in opposite orientation,
    # it should NOT contain TTAGGG repeats itself. Its rev-comp is what
    # would match TERC. We compute its telo content.
    fwd_telo = sum(linker_seq.count(k) for k in ['TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG'])
    rev_telo = sum(linker_seq.count(k) for k in ['CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA'])
    return fwd_telo, rev_telo


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    real_dir = Path(__file__).parent / "real_telotrons"
    fa_path = real_dir / "converging_linkers.fa"
    if not fa_path.exists():
        print(f"Linker FASTA not found at {fa_path}. Run its_mechanism_v2.py first.")
        return

    linkers = load_linker_fasta(fa_path)
    print(f"Loaded {len(linkers)} linkers")

    # Filter to converging and split by length bucket
    by_class_len = defaultdict(list)
    for sid, seq in linkers:
        if 'CONVERGING_WITHLINKER' in sid:
            ll = len(seq)
            if 5 <= ll <= 30:
                by_class_len['short_5_30'].append((sid, seq))
            elif 30 < ll <= 120:
                by_class_len['medium_30_120'].append((sid, seq))
            elif ll > 120:
                by_class_len['long_120plus'].append((sid, seq))

    print(f"\nLinker length buckets (CONVERGING_WITHLINKER only):")
    for k, v in by_class_len.items():
        print(f"  {k}: {len(v)}")

    # ----- A. Self-similarity test on the medium bucket (TERC-frag prediction) -----
    print(f"\n--- A. SELF-SIMILARITY (Nergadze TERC-template prediction) ---")
    print(f"If linkers share a common template, large clusters should emerge.\n")

    for bucket_name, seqs in by_class_len.items():
        if len(seqs) < 50: continue
        # Sample to keep clustering tractable
        sample = seqs[:500] if len(seqs) > 500 else seqs
        clusters = kmer_jaccard_clusters(sample, k=8, min_jaccard=0.15)
        n_singletons = sum(1 for c, _ in clusters if len(c) == 1)
        n_clusters_ge3 = sum(1 for c, _ in clusters if len(c) >= 3)
        n_clusters_ge10 = sum(1 for c, _ in clusters if len(c) >= 10)
        biggest = max(clusters, key=lambda x: len(x[0]))
        biggest_size = len(biggest[0])
        print(f"  Bucket {bucket_name} (n={len(sample)}):")
        print(f"    Singletons:     {n_singletons} ({100*n_singletons/len(sample):.0f}%)")
        print(f"    Clusters ≥3:    {n_clusters_ge3}")
        print(f"    Clusters ≥10:   {n_clusters_ge10}")
        print(f"    Biggest cluster: {biggest_size} members")
        if biggest_size >= 5:
            print(f"      Sample seqs from biggest cluster:")
            for sid, sq in biggest[0][:3]:
                print(f"        {sid[:60]}: {sq[:80]}")

    # ----- B. Global k-mer frequencies -----
    print(f"\n--- B. GLOBAL KMER FREQUENCIES (all medium-length linkers) ---")
    medium = by_class_len.get('medium_30_120', [])
    if medium:
        kc, nlk = linker_global_kmer_freq(medium, k=10)
        # Most abundant 10-mers across linkers
        top = sorted(nlk.items(), key=lambda x: -x[1])[:15]
        print(f"  Top 10-mers by # linkers containing them (n={len(medium)} linkers):")
        for km, n in top:
            pct = 100 * n / len(medium)
            print(f"    {km}: in {n} linkers ({pct:.1f}%)  total occurrences: {kc[km]}")

    # ----- C. Composition -----
    print(f"\n--- C. COMPOSITION (GC, dinucleotide) ---")
    for bucket_name, seqs in by_class_len.items():
        if not seqs: continue
        comps = [composition(s) for _, s in seqs]
        comps = [c for c in comps if c]
        if not comps: continue
        mean_gc = sum(c['gc'] for c in comps) / len(comps)
        mean_cpg = sum(c['cpg_freq'] for c in comps) / len(comps)
        mean_tpa = sum(c['tpa_freq'] for c in comps) / len(comps)
        print(f"  {bucket_name}: mean GC={mean_gc:.3f}, CpG freq={mean_cpg:.4f}, TpA freq={mean_tpa:.4f}, n={len(comps)}")

    # ----- D. Reverse-complement / orientation test (Nergadze key prediction) -----
    print(f"\n--- D. ORIENTATION TEST (Nergadze: linker is rev-comp to flanks) ---")
    print(f"  TERC-retrotransposition predicts linker doesn't contain TTAGGG repeats")
    print(f"  but its sequence relates to TERC (which has CCCTAA as template).\n")

    bucket_telo_stats = {}
    for bucket_name, seqs in by_class_len.items():
        if not seqs: continue
        n_no_telo = 0
        n_fwd_telo = 0
        n_rev_telo = 0
        n_both = 0
        for sid, sq in seqs:
            f, r = linker_rc_relative_to_flanks(sid, sq)
            if f == 0 and r == 0:
                n_no_telo += 1
            elif f > 0 and r == 0:
                n_fwd_telo += 1
            elif r > 0 and f == 0:
                n_rev_telo += 1
            else:
                n_both += 1
        n = len(seqs)
        bucket_telo_stats[bucket_name] = {
            'no_telo_pct': 100*n_no_telo/n,
            'fwd_only_pct': 100*n_fwd_telo/n,
            'rev_only_pct': 100*n_rev_telo/n,
            'both_pct': 100*n_both/n,
            'n': n,
        }
        print(f"  {bucket_name} (n={n}):")
        print(f"    No telo kmers (truly non-repeat): {100*n_no_telo/n:.1f}%")
        print(f"    Only fwd (TTAGGG): {100*n_fwd_telo/n:.1f}%")
        print(f"    Only rev (CCCTAA): {100*n_rev_telo/n:.1f}%")
        print(f"    Both: {100*n_both/n:.1f}%")

    # Save outputs
    out_summary = real_dir / "linker_analysis_summary.json"
    summary = {
        'total_linkers': len(linkers),
        'length_buckets': {k: len(v) for k, v in by_class_len.items()},
        'orientation_breakdown': bucket_telo_stats,
    }
    with open(out_summary, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSummary written to: {out_summary}")


if __name__ == '__main__':
    main()
