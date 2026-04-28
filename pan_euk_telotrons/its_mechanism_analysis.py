#!/usr/bin/env python3
"""
ITS-mechanism diagnostic battery for telotron candidates.

Implements the diagnostics from the synthesized ITS literature review:
  1. Orientation classification (4 classes):
     - SINGLE_G : single-strand TTAGGG-only array
     - SINGLE_C : single-strand CCCTAA-only array (rev-comp orientation)
     - CONVERGING_NOLINKER : head-to-head TTAGGG/CCCTAA, ≤4bp gap
     - CONVERGING_WITHLINKER : converging arrays with non-repetitive middle linker
  2. Variant-repeat composition (TCAGGG/TGAGGG/GTAGGG/etc.) — TTI vs canonical signal
  3. Microhomology at junctions (5' and 3' of intron, vs canonical TTAGGG)
  4. TSDs (target site duplications) at the flanks
  5. NR2C/F binding motif scan (GGGTCA-like DR0/6/7 spacings) in flanks
  6. Linker extraction (for CONVERGING_WITHLINKER class) for downstream BLAST
  7. Age proxy (purity/canonicality)

Inputs:  real_telotrons/*.tsv (validated telotrons with intron_seq + flanks)
Outputs: real_telotrons/its_mechanism_classification.tsv
         real_telotrons/converging_linkers.fa  (for BLAST)
         real_telotrons/its_mechanism_summary.json
"""

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------

def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

# Canonical telomeric repeats (as G-strand)
CANONICAL_FAMILIES = {
    'TTAGGG':   {'fwd': _rotations('TTAGGG'),
                 'rev': _rotations('CCCTAA')},
    'TTTAGGG':  {'fwd': _rotations('TTTAGGG'),
                 'rev': _rotations('CCCTAAA')},
    'TTAGG':    {'fwd': _rotations('TTAGG'),
                 'rev': _rotations('CCTAA')},
    'TTTGGG':   {'fwd': _rotations('TTTGGG'),
                 'rev': _rotations('CCCAAA')},
    'TTGGGG':   {'fwd': _rotations('TTGGGG'),
                 'rev': _rotations('CCCCAA')},
}

# Variant repeats associated with TTI / NR2C/F mechanism (Marzec 2015):
# These are interspersed within ALT telomere arrays
TTI_VARIANTS_FROM_TTAGGG = ['TCAGGG', 'TGAGGG', 'GTAGGG', 'TTGGGG', 'TTCGGG']
TTI_VARIANTS_REV = [s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
                    for s in TTI_VARIANTS_FROM_TTAGGG]

# NR2C/F binding motif: GGGTCA half-sites (or AGGTCA also reported)
# DR0 = direct repeat, 0bp spacing; DR6/7 = 6 or 7bp spacing
NR2CF_HALFSITE = 'GGGTCA'
NR2CF_HALFSITE_RC = 'TGACCC'

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def rc(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def find_runs(seq, kmers, min_run_bp=12):
    """Find runs of any kmer in `kmers`. Return list of (start, end) tuples."""
    n = len(seq)
    if n < min_run_bp:
        return []
    covered = bytearray(n)
    for k in kmers:
        i = 0
        klen = len(k)
        while True:
            idx = seq.find(k, i)
            if idx == -1:
                break
            for j in range(idx, min(idx + klen, n)):
                covered[j] = 1
            i = idx + 1
    runs = []
    in_run = False
    rs = 0
    for i in range(n):
        if covered[i] and not in_run:
            in_run = True
            rs = i
        elif not covered[i] and in_run:
            in_run = False
            if i - rs >= min_run_bp:
                runs.append((rs, i))
    if in_run and n - rs >= min_run_bp:
        runs.append((rs, n))
    return runs


def detect_canonical_family(seq):
    """Return the canonical family with the most coverage in seq.
    Returns (family_name, fwd_bp, rev_bp, fwd_run_bp, rev_run_bp)."""
    best = (None, 0, 0, 0, 0)
    for fam, sets in CANONICAL_FAMILIES.items():
        fwd_runs = find_runs(seq, sets['fwd'], min_run_bp=10)
        rev_runs = find_runs(seq, sets['rev'], min_run_bp=10)
        fwd_bp = sum(e - s for s, e in fwd_runs)
        rev_bp = sum(e - s for s, e in rev_runs)
        fwd_run = max((e - s for s, e in fwd_runs), default=0)
        rev_run = max((e - s for s, e in rev_runs), default=0)
        score = fwd_bp + rev_bp
        if score > best[1] + best[2]:
            best = (fam, fwd_bp, rev_bp, fwd_run, rev_run)
    return best


def classify_orientation(seq, family):
    """
    4-class orientation:
      SINGLE_G                : only fwd repeats (≥10% coverage)
      SINGLE_C                : only rev repeats
      CONVERGING_NOLINKER     : fwd block then rev block with ≤4bp gap
      CONVERGING_WITHLINKER   : fwd block then rev block with >4bp gap
                                (or rev then fwd)
    Also returns: fwd_runs, rev_runs, linker (str or None), linker_pos (start, end) or None
    """
    sets = CANONICAL_FAMILIES[family]
    fwd_runs = find_runs(seq, sets['fwd'], min_run_bp=10)
    rev_runs = find_runs(seq, sets['rev'], min_run_bp=10)

    fwd_bp = sum(e - s for s, e in fwd_runs)
    rev_bp = sum(e - s for s, e in rev_runs)
    n = len(seq)

    has_fwd = fwd_bp >= 10 and (fwd_bp / n) >= 0.05
    has_rev = rev_bp >= 10 and (rev_bp / n) >= 0.05

    if has_fwd and not has_rev:
        return 'SINGLE_G', fwd_runs, rev_runs, None, None
    if has_rev and not has_fwd:
        return 'SINGLE_C', fwd_runs, rev_runs, None, None
    if not has_fwd and not has_rev:
        return 'NONE', fwd_runs, rev_runs, None, None

    # Both present. Decide linker.
    # Find the largest fwd run and the largest rev run, and the gap between them
    fwd_max = max(fwd_runs, key=lambda r: r[1] - r[0])
    rev_max = max(rev_runs, key=lambda r: r[1] - r[0])

    # Order them on the sequence
    if fwd_max[0] < rev_max[0]:
        first, second = fwd_max, rev_max
        order = 'G_then_C'  # 5'-(TTAGGG)n-(CCCTAA)m-3' = converging head-to-head
    else:
        first, second = rev_max, fwd_max
        order = 'C_then_G'  # 5'-(CCCTAA)n-(TTAGGG)m-3' = diverging tail-to-tail

    gap_start = first[1]
    gap_end = second[0]
    linker_len = gap_end - gap_start
    linker = seq[gap_start:gap_end] if linker_len > 0 else ''

    if order == 'G_then_C':
        if linker_len <= 4:
            return 'CONVERGING_NOLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
        else:
            return 'CONVERGING_WITHLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
    else:
        # diverging — uncommon, treat as separate class
        if linker_len <= 4:
            return 'DIVERGING_NOLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
        else:
            return 'DIVERGING_WITHLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)


def variant_repeat_fraction(seq, family='TTAGGG'):
    """Fraction of hexamer windows that are 1-mismatch variants of canonical TTAGGG.
    Returns dict of variant_kmer -> count, plus canonical count and total."""
    if family != 'TTAGGG':
        # only well-characterized for TTAGGG
        return None

    canonical_kmers = CANONICAL_FAMILIES['TTAGGG']['fwd'] | CANONICAL_FAMILIES['TTAGGG']['rev']
    canonical_count = 0
    variant_counts = Counter()
    total_hex_windows = 0
    n = len(seq)

    # Slide 6bp window
    for i in range(0, n - 5):
        kmer = seq[i:i+6]
        total_hex_windows += 1
        if kmer in canonical_kmers:
            canonical_count += 1
            continue
        # Check if 1-mismatch from any canonical
        for ck in canonical_kmers:
            mismatches = sum(a != b for a, b in zip(kmer, ck))
            if mismatches == 1:
                variant_counts[kmer] += 1
                break

    return {
        'canonical': canonical_count,
        'variants': dict(variant_counts),
        'total_windows': total_hex_windows,
        'canonical_frac': canonical_count / total_hex_windows if total_hex_windows else 0,
        'variant_frac': sum(variant_counts.values()) / total_hex_windows if total_hex_windows else 0,
    }


def microhomology_at_junction(flank_3p_end, intron_5p_start, family_kmers,
                               min_match=2, max_match=8):
    """
    Microhomology between exon 3' end and the canonical TTAGGG-pattern at intron 5'.
    flank_3p_end: last N bp of upstream exon (3' end of flank)
    intron_5p_start: first N bp of intron
    Returns: longest match length where the exon's 3' end matches a canonical
             telomeric kmer at the start of the intron.
    """
    best = 0
    matches = []
    for k in family_kmers:
        for shift in range(len(k)):
            # Try aligning the kmer with shift to the boundary
            check_seq = k[shift:] + k[:shift]
            for L in range(max_match, min_match - 1, -1):
                if L > len(flank_3p_end):
                    continue
                # The intron should "start" with check_seq[:L] worth of repeat
                # And the exon's last L bp should match check_seq[-L:]
                tail_of_exon = flank_3p_end[-L:]
                head_of_intron = intron_5p_start[:L]
                # Look for case where exon ends in something that could be telo
                # i.e. last L bp of exon == first L bp of intron's repeat-pattern
                if tail_of_exon == head_of_intron and 'TT' in tail_of_exon + head_of_intron:
                    # additionally check the tail looks like a telo subsequence
                    if any(tail_of_exon in fk for fk in family_kmers) or \
                       tail_of_exon[:min(L,4)] in {k_[:4] for k_ in family_kmers}:
                        if L > best:
                            best = L
                            matches.append((L, tail_of_exon))
                        break
    return best


def find_tsd(left_flank, right_flank, min_len=3, max_len=25):
    """Exact TSD: last N of left flank matches first N of right flank."""
    best = (0, '')
    for L in range(max_len, min_len - 1, -1):
        if L > len(left_flank) or L > len(right_flank):
            continue
        if left_flank[-L:] == right_flank[:L]:
            return (L, left_flank[-L:])
    return best


def scan_nr2cf(seq, max_spacing=10):
    """Scan for NR2C/F-type binding motifs: GGGTCA half-sites in DR0/6/7 spacings."""
    hits = []
    half = NR2CF_HALFSITE
    half_rc = NR2CF_HALFSITE_RC
    L = len(half)

    fwd_pos = [i for i in range(len(seq) - L + 1) if seq[i:i+L] == half]
    rc_pos = [i for i in range(len(seq) - L + 1) if seq[i:i+L] == half_rc]

    n_fwd = len(fwd_pos)
    n_rc = len(rc_pos)

    direct_repeats = []
    for i in range(len(fwd_pos)):
        for j in range(i+1, len(fwd_pos)):
            spacing = fwd_pos[j] - fwd_pos[i] - L
            if 0 <= spacing <= max_spacing:
                direct_repeats.append((fwd_pos[i], fwd_pos[j], spacing))
    inverted_repeats = []
    for i in fwd_pos:
        for j in rc_pos:
            if i < j:
                spacing = j - i - L
            else:
                spacing = i - j - L
            if 0 <= spacing <= max_spacing:
                inverted_repeats.append((min(i, j), max(i, j), spacing))

    return {
        'half_sites_fwd': n_fwd,
        'half_sites_rc': n_rc,
        'direct_repeats': len(direct_repeats),
        'inverted_repeats': len(inverted_repeats),
        'direct_spacings': sorted(set(s for _, _, s in direct_repeats))[:10],
    }


# ------------------------------------------------------------------
# Per-telotron analysis
# ------------------------------------------------------------------

def analyze_telotron(rec):
    """Run the full diagnostic battery on one telotron record."""
    seq = rec.get('intron_seq', '').upper()
    left = rec.get('left_flank_100bp', '').upper()
    right = rec.get('right_flank_100bp', '').upper()

    if not seq:
        return None

    # 1. Detect canonical family
    fam, fwd_bp, rev_bp, fwd_run, rev_run = detect_canonical_family(seq)
    if fam is None:
        return {
            'classification': 'NO_TELO', 'family': None,
            'fwd_bp': 0, 'rev_bp': 0, 'fwd_run_max': 0, 'rev_run_max': 0,
            'n_fwd_runs': 0, 'n_rev_runs': 0,
            'linker_info': None, 'variant_info': None,
            'microhomology_5p': 0, 'microhomology_3p': 0,
            'tsd_length': 0, 'tsd_seq': '', 'nr2cf': None,
        }

    # 2. Classify orientation
    orient_class, fwd_runs, rev_runs, linker, linker_pos = classify_orientation(seq, fam)

    # 3. Variant repeat composition (only for TTAGGG family)
    variant_info = variant_repeat_fraction(seq, fam) if fam == 'TTAGGG' else None

    # 4. Microhomology at boundaries
    fam_kmers = CANONICAL_FAMILIES[fam]['fwd'] | CANONICAL_FAMILIES[fam]['rev']
    mh_5p = microhomology_at_junction(left, seq, fam_kmers) if left else 0
    mh_3p = microhomology_at_junction(seq, right, fam_kmers) if right else 0

    # 5. TSDs
    tsd_len, tsd_seq = find_tsd(left, right)

    # 6. NR2C/F motifs (in flanks combined)
    flank_combined = left + right
    nr2cf = scan_nr2cf(flank_combined) if flank_combined else None

    # 7. Linker properties (if any)
    linker_info = None
    if linker:
        linker_info = {
            'length': len(linker),
            'gc_content': (linker.count('G') + linker.count('C')) / len(linker) if linker else 0,
            'sequence': linker,
        }

    return {
        'family': fam,
        'fwd_bp': fwd_bp,
        'rev_bp': rev_bp,
        'fwd_run_max': fwd_run,
        'rev_run_max': rev_run,
        'classification': orient_class,
        'n_fwd_runs': len(fwd_runs),
        'n_rev_runs': len(rev_runs),
        'linker_info': linker_info,
        'variant_info': variant_info,
        'microhomology_5p': mh_5p,
        'microhomology_3p': mh_3p,
        'tsd_length': tsd_len,
        'tsd_seq': tsd_seq,
        'nr2cf': nr2cf,
    }


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    real_dir = Path(__file__).parent / "real_telotrons"
    out_tsv = real_dir / "its_mechanism_classification.tsv"
    out_linker_fa = real_dir / "converging_linkers.fa"
    out_summary = real_dir / "its_mechanism_summary.json"

    all_records = []
    by_class = defaultdict(int)
    by_class_species = defaultdict(set)

    for tsv in sorted(real_dir.glob("*_real_telotrons.tsv")):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if not row.get('intron_seq'):
                    continue
                analysis = analyze_telotron(row)
                if analysis is None:
                    continue
                rec = {
                    'species': row['species'],
                    'contig': row['contig'],
                    'start': row['start'],
                    'end': row['end'],
                    'length': row['length'],
                    'gene_name': row.get('gene_name', ''),
                    **analysis,
                }
                all_records.append(rec)
                by_class[analysis['classification']] += 1
                by_class_species[analysis['classification']].add(row['species'])

    # Write TSV
    cols = [
        'species', 'contig', 'start', 'end', 'length', 'gene_name',
        'family', 'classification',
        'fwd_bp', 'rev_bp', 'fwd_run_max', 'rev_run_max',
        'n_fwd_runs', 'n_rev_runs',
        'linker_length', 'linker_gc',
        'canonical_frac', 'variant_frac',
        'microhomology_5p', 'microhomology_3p',
        'tsd_length', 'tsd_seq',
        'nr2cf_half_sites', 'nr2cf_direct_repeats', 'nr2cf_inverted_repeats',
        'linker_seq',
    ]
    with open(out_tsv, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in all_records:
            li = r.get('linker_info') or {}
            vi = r.get('variant_info') or {}
            ni = r.get('nr2cf') or {}
            row = [
                r['species'], r['contig'], r['start'], r['end'], r['length'],
                r['gene_name'], r['family'], r['classification'],
                r['fwd_bp'], r['rev_bp'], r['fwd_run_max'], r['rev_run_max'],
                r['n_fwd_runs'], r['n_rev_runs'],
                li.get('length', ''), f"{li.get('gc_content', 0):.3f}" if li else '',
                f"{vi.get('canonical_frac', 0):.3f}" if vi else '',
                f"{vi.get('variant_frac', 0):.3f}" if vi else '',
                r['microhomology_5p'], r['microhomology_3p'],
                r['tsd_length'], r['tsd_seq'],
                (ni.get('half_sites_fwd', 0) + ni.get('half_sites_rc', 0)) if ni else '',
                ni.get('direct_repeats', '') if ni else '',
                ni.get('inverted_repeats', '') if ni else '',
                li.get('sequence', '') if li else '',
            ]
            f.write('\t'.join(str(x) for x in row) + '\n')

    # Write linker FASTA for BLAST
    with open(out_linker_fa, 'w') as f:
        for r in all_records:
            li = r.get('linker_info')
            if li and li.get('length', 0) >= 15:
                seqid = (f"{r['species']}__{r['contig']}__{r['start']}_"
                         f"{r['end']}__{r['classification']}")
                f.write(f">{seqid}\n{li['sequence']}\n")

    # Summary stats
    summary = {
        'total': len(all_records),
        'by_classification': {k: v for k, v in by_class.items()},
        'species_per_class': {k: len(v) for k, v in by_class_species.items()},
        'family_distribution': dict(Counter(r['family'] for r in all_records)),
    }

    # Mean variant_frac per class (TTAGGG family only) — TTI signal
    by_class_variant = defaultdict(list)
    for r in all_records:
        if r['family'] == 'TTAGGG' and r.get('variant_info'):
            by_class_variant[r['classification']].append(r['variant_info']['variant_frac'])
    summary['mean_variant_frac_per_class'] = {
        k: sum(v) / len(v) if v else 0 for k, v in by_class_variant.items()
    }

    # TSD presence per class
    by_class_tsd = defaultdict(list)
    for r in all_records:
        by_class_tsd[r['classification']].append(r['tsd_length'])
    summary['mean_tsd_per_class'] = {
        k: sum(v) / len(v) if v else 0 for k, v in by_class_tsd.items()
    }
    summary['pct_with_tsd_per_class'] = {
        k: 100 * sum(1 for x in v if x >= 3) / len(v) if v else 0
        for k, v in by_class_tsd.items()
    }

    # NR2C/F enrichment per class
    by_class_nr2cf = defaultdict(list)
    for r in all_records:
        ni = r.get('nr2cf')
        if ni:
            total_hs = ni['half_sites_fwd'] + ni['half_sites_rc']
            by_class_nr2cf[r['classification']].append(total_hs)
    summary['mean_nr2cf_halfsites_per_class'] = {
        k: sum(v) / len(v) if v else 0 for k, v in by_class_nr2cf.items()
    }

    # Linker length distribution
    linker_lens = [r['linker_info']['length'] for r in all_records
                   if r.get('linker_info') and r['linker_info']['length'] >= 5]
    summary['linker_length_distribution'] = {
        'n': len(linker_lens),
        'min': min(linker_lens) if linker_lens else 0,
        'max': max(linker_lens) if linker_lens else 0,
        'mean': sum(linker_lens) / len(linker_lens) if linker_lens else 0,
        'median': sorted(linker_lens)[len(linker_lens)//2] if linker_lens else 0,
        # Distribution buckets
        'short_5_30': sum(1 for x in linker_lens if 5 <= x <= 30),
        'medium_30_120': sum(1 for x in linker_lens if 30 < x <= 120),
        'long_120plus': sum(1 for x in linker_lens if x > 120),
    }

    with open(out_summary, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    # Print to stdout
    print(f"\n{'='*70}")
    print(f"  ITS-MECHANISM CLASSIFICATION SUMMARY")
    print(f"{'='*70}\n")
    print(f"Total telotrons analyzed: {summary['total']}\n")
    print(f"Orientation class:")
    for k, v in sorted(summary['by_classification'].items(), key=lambda x: -x[1]):
        n_sp = summary['species_per_class'][k]
        print(f"  {k:<28s} {v:>6,}  ({n_sp} species)")
    print(f"\nRepeat family:")
    for k, v in sorted(summary['family_distribution'].items(), key=lambda x: -x[1]):
        kstr = str(k) if k else 'NONE'
        print(f"  {kstr:<28s} {v:>6,}")
    print(f"\nMean variant fraction per class (TTI signal, TTAGGG family):")
    for k, v in sorted(summary['mean_variant_frac_per_class'].items()):
        print(f"  {k:<28s} {v:>.4f}")
    print(f"\nTSD presence per class:")
    for k, v in sorted(summary['pct_with_tsd_per_class'].items()):
        print(f"  {k:<28s} {v:>5.1f}% have TSD ≥3bp")
    print(f"\nMean NR2C/F half-sites in flanks (200bp combined):")
    for k, v in sorted(summary['mean_nr2cf_halfsites_per_class'].items()):
        print(f"  {k:<28s} {v:>.3f}")
    if linker_lens:
        print(f"\nLinker length distribution (n={len(linker_lens)} CONVERGING_WITHLINKER):")
        print(f"  Range: {summary['linker_length_distribution']['min']}-"
              f"{summary['linker_length_distribution']['max']} bp")
        print(f"  Mean:  {summary['linker_length_distribution']['mean']:.1f} bp")
        print(f"  Median: {summary['linker_length_distribution']['median']} bp")
        print(f"  Short (5-30bp, NHEJ-scar-like): {summary['linker_length_distribution']['short_5_30']}")
        print(f"  Medium (30-120bp, TERC-fragment-like): {summary['linker_length_distribution']['medium_30_120']}")
        print(f"  Long (120+bp, captured-fragment-like): {summary['linker_length_distribution']['long_120plus']}")

    print(f"\nOutputs written to:")
    print(f"  {out_tsv}")
    print(f"  {out_linker_fa}")
    print(f"  {out_summary}")


if __name__ == '__main__':
    main()
