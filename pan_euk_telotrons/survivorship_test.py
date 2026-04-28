#!/usr/bin/env python3
"""
Survivorship test: Are telomeric tracts found in INTRONS but absent from EXONS?

If insertions of telomeric repeats happen genome-wide but only those that
land in introns survive (because exonic insertions disrupt protein and get
purged), we should find:
  - Telomeric tracts COMMON in introns (telotrons)
  - Telomeric tracts RARE/ABSENT in CDS (purged by selection)
  - Telomeric tracts at NEUTRAL frequency in INTERGENIC (no selection)

Test for the primary haptophyte MAG:
  1. Extract all introns, all CDS regions, and intergenic windows
  2. For each, count telomeric content (≥10% TTAGGG/CCCTAA)
  3. Compute density per Mb in each compartment
  4. Test: intron telo >> CDS telo (survivorship)
     Test: intergenic telo ≈ intron telo (neutral)
     Test: intergenic telo >> CDS telo (selection acts only on CDS)
"""

import csv
import json
import re
import gzip
from collections import defaultdict, Counter
from pathlib import Path


TARA_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags')
GENOME_DIR = TARA_DIR / 'smags' / 'contigs_individual'
GFF_DIR = TARA_DIR / 'smags' / 'gff_individual' / 'GFF'

# Use top telotron MAGs
TARGET_MAGS = [
    'TARA_PSW_86_MAG_00284',
    'TARA_PON_109_MAG_00250',
    'TARA_ARC_108_MAG_00319',
    'TARA_PSE_93_MAG_00226',
    'TARA_AOS_82_MAG_00154',
]

# Telomeric kmers (TTAGGG family, both strands)
TELO_FWD = {f'TTAGGG'[i:] + 'TTAGGG'[:i] for i in range(6)}
TELO_REV = {f'CCCTAA'[i:] + 'CCCTAA'[:i] for i in range(6)}
TELO_KMERS = TELO_FWD | TELO_REV


def telo_density(seq, window=200, min_kmers=2):
    """Count windows of `window` bp with ≥min_kmers telomeric kmer hits.
    Returns (n_windows_with_telo, total_windows, n_telo_kmer_total)."""
    n = len(seq)
    if n < window:
        return 0, 0, 0
    n_telo_total = 0
    for k in TELO_KMERS:
        i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1:
                break
            n_telo_total += 1
            i = idx + 1
    # Sliding window count
    n_windows = n - window + 1
    n_pos_windows = 0
    # Build cumulative kmer-hit array for fast windowing
    hits = bytearray(n)
    for k in TELO_KMERS:
        i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1:
                break
            hits[idx] = 1
            i = idx + 1
    # Sum each window
    win_count = sum(hits[:window])
    if win_count >= min_kmers:
        n_pos_windows += 1
    for i in range(1, n_windows):
        win_count += hits[i + window - 1] - hits[i - 1]
        if win_count >= min_kmers:
            n_pos_windows += 1
    return n_pos_windows, n_windows, n_telo_total


def load_fasta_streaming(path):
    """Yield (contig_id, seq) tuples from FASTA."""
    curr_id = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if curr_id:
                    yield curr_id, ''.join(parts).upper()
                curr_id = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if curr_id:
        yield curr_id, ''.join(parts).upper()


def parse_gff_features(gff_path):
    """Parse GFF, return per-contig list of (type, start_0, end, strand) intervals.
    type is 'CDS', 'intron', 'gene', or 'mRNA'."""
    contigs = defaultdict(lambda: {'cds': [], 'intron': [], 'gene': [], 'mrna': []})
    mrna_features = defaultdict(list)  # mRNA_id -> [(start, end)]
    mrna_contig = {}
    mrna_strand = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            ftype = parts[2]
            contig = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            attrs = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))

            if ftype == 'gene':
                contigs[contig]['gene'].append((start, end, strand))
            elif ftype in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                mrna_contig[mid] = contig
                mrna_strand[mid] = strand
                contigs[contig]['mrna'].append((start, end, strand))
            elif ftype == 'CDS' or ftype == 'exon':
                parent = attrs.get('Parent', '')
                # Use exons or CDSes to compute introns
                mrna_features[parent].append((start, end, ftype))
                if ftype == 'CDS':
                    contigs[contig]['cds'].append((start, end, strand))

    # Compute introns from exon gaps per mRNA
    for mid, feats in mrna_features.items():
        if mid not in mrna_contig: continue
        contig = mrna_contig[mid]
        strand = mrna_strand[mid]
        # Use exons preferentially; fall back to CDS
        exons = sorted(set((s, e) for s, e, t in feats if t == 'exon'))
        if not exons:
            exons = sorted(set((s, e) for s, e, t in feats if t == 'CDS'))
        if len(exons) < 2: continue
        for i in range(len(exons) - 1):
            i_start = exons[i][1]
            i_end = exons[i+1][0]
            if i_end > i_start:
                contigs[contig]['intron'].append((i_start, i_end, strand))

    return dict(contigs)


def merge_intervals(intervals):
    """Merge overlapping (start, end[, strand]) intervals. Return [(start, end)]."""
    if not intervals: return []
    sorted_iv = sorted((s, e) for s, e, *_ in intervals)
    merged = [sorted_iv[0]]
    for s, e in sorted_iv[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def intergenic_intervals(contig_len, gene_intervals, intron_intervals):
    """Compute regions that are not in any gene or intron."""
    occupied = merge_intervals(gene_intervals + intron_intervals)
    intergenic = []
    last_end = 0
    for s, e in occupied:
        if s > last_end:
            intergenic.append((last_end, s))
        last_end = max(last_end, e)
    if last_end < contig_len:
        intergenic.append((last_end, contig_len))
    return intergenic


def main():
    results = {}

    for mag in TARGET_MAGS:
        gff = GFF_DIR / f"{mag}.gmove.gff"
        fa = GENOME_DIR / f"{mag}.fa"
        if not gff.exists() or not fa.exists():
            print(f"  Skip {mag}: gff or fa missing")
            continue

        print(f"\n=== {mag} ===")
        contig_features = parse_gff_features(gff)

        compartment_stats = {
            'cds': {'bp': 0, 'telo_kmers': 0, 'telo_windows': 0, 'total_windows': 0},
            'intron': {'bp': 0, 'telo_kmers': 0, 'telo_windows': 0, 'total_windows': 0},
            'intergenic': {'bp': 0, 'telo_kmers': 0, 'telo_windows': 0, 'total_windows': 0},
        }

        n_contigs = 0
        for contig_id, seq in load_fasta_streaming(fa):
            n_contigs += 1
            n = len(seq)
            feats = contig_features.get(contig_id, None)
            if not feats:
                continue

            cds = merge_intervals(feats.get('cds', []))
            introns = merge_intervals(feats.get('intron', []))
            genes = merge_intervals(feats.get('gene', []))
            interg = intergenic_intervals(n, genes, introns)

            # Subtract introns from cds in case of overlap
            # (CDS should NOT overlap introns, but some annotations have errors)

            for compartment, intervals in [('cds', cds), ('intron', introns), ('intergenic', interg)]:
                for s, e in intervals:
                    if s < 0: s = 0
                    if e > n: e = n
                    if e <= s: continue
                    sub = seq[s:e]
                    pos_w, tot_w, telo_count = telo_density(sub, window=200, min_kmers=2)
                    compartment_stats[compartment]['bp'] += (e - s)
                    compartment_stats[compartment]['telo_kmers'] += telo_count
                    compartment_stats[compartment]['telo_windows'] += pos_w
                    compartment_stats[compartment]['total_windows'] += tot_w

        # Print results
        print(f"  Contigs: {n_contigs}")
        print(f"  {'Compartment':<12s} {'bp (Mb)':>10s} {'kmers':>10s} {'kmers/Mb':>12s} {'+windows':>10s} {'window%':>10s}")
        for comp, st in compartment_stats.items():
            mb = st['bp'] / 1e6
            kmer_per_mb = st['telo_kmers'] / mb if mb > 0 else 0
            window_pct = 100 * st['telo_windows'] / st['total_windows'] if st['total_windows'] > 0 else 0
            print(f"  {comp:<12s} {mb:>10.2f} {st['telo_kmers']:>10,} {kmer_per_mb:>12.1f} "
                  f"{st['telo_windows']:>10,} {window_pct:>9.3f}%")

        # Compute survivorship ratios
        if compartment_stats['cds']['bp'] > 0 and compartment_stats['intron']['bp'] > 0:
            intron_density = compartment_stats['intron']['telo_kmers'] / (compartment_stats['intron']['bp'] / 1e6)
            cds_density = compartment_stats['cds']['telo_kmers'] / (compartment_stats['cds']['bp'] / 1e6)
            interg_density = compartment_stats['intergenic']['telo_kmers'] / (compartment_stats['intergenic']['bp'] / 1e6) \
                             if compartment_stats['intergenic']['bp'] > 0 else 0
            ratio_intron_cds = intron_density / cds_density if cds_density > 0 else float('inf')
            ratio_interg_cds = interg_density / cds_density if cds_density > 0 else float('inf')
            ratio_intron_interg = intron_density / interg_density if interg_density > 0 else float('inf')
            print(f"\n  Intron telo / CDS telo:        {ratio_intron_cds:>8.1f}× (survivorship signal)")
            print(f"  Intergenic telo / CDS telo:    {ratio_interg_cds:>8.1f}× (selection on CDS)")
            print(f"  Intron telo / Intergenic telo: {ratio_intron_interg:>8.2f}× (intron-specific?)")

        results[mag] = compartment_stats

    # Save
    out_json = Path(__file__).parent / "real_telotrons" / "survivorship_test.json"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {out_json}")


if __name__ == '__main__':
    main()
