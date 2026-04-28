#!/usr/bin/env python3
"""
Extract TRUE non-repeat linkers (linkers between converging arrays that
contain NO telomeric kmers).
These are the legitimate candidates for:
  - NHEJ scar (random non-templated nt)
  - TERC-retrotransposition (Nergadze 2007)
  - Captured DNA fragment from elsewhere
"""

import json
import re
from collections import Counter, defaultdict
from pathlib import Path


FWD_TELO_KMERS = {'TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG'}
REV_TELO_KMERS = {'CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA'}
ALL_TELO_KMERS = FWD_TELO_KMERS | REV_TELO_KMERS


def has_telo(seq):
    return any(k in seq for k in ALL_TELO_KMERS)


def has_telo_count(seq):
    return sum(seq.count(k) for k in ALL_TELO_KMERS)


def composition(seq):
    if not seq:
        return None
    n = len(seq)
    c = Counter(seq)
    return {
        'gc': (c.get('G', 0) + c.get('C', 0)) / n,
        'a': c.get('A', 0) / n,
        't': c.get('T', 0) / n,
        'g': c.get('G', 0) / n,
        'c': c.get('C', 0) / n,
    }


def main():
    real_dir = Path(__file__).parent / "real_telotrons"
    fa_path = real_dir / "converging_linkers.fa"

    # Load all linkers
    linkers = []
    curr_id = None; parts = []
    with open(fa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_id:
                    linkers.append((curr_id, ''.join(parts).upper()))
                curr_id = line[1:]
                parts = []
            else:
                parts.append(line)
    if curr_id:
        linkers.append((curr_id, ''.join(parts).upper()))

    # Filter to TRUE non-repeat linkers (no telomeric kmers at all)
    pure_nonrepeat = [(sid, seq) for sid, seq in linkers if not has_telo(seq)]
    contaminated = [(sid, seq) for sid, seq in linkers if has_telo(seq)]

    # Among contaminated, those with low telo content (e.g. ≤2 hits)
    low_telo = [(sid, seq) for sid, seq in contaminated if has_telo_count(seq) <= 2]

    print(f"Total converging linkers (≥15bp): {len(linkers)}")
    print(f"  Pure non-repeat (no telo kmers): {len(pure_nonrepeat)}  ({100*len(pure_nonrepeat)/len(linkers):.1f}%)")
    print(f"  Contaminated (≥1 telo kmer):     {len(contaminated)}")
    print(f"  Low-telo (≤2 telo kmers):        {len(low_telo)}")

    # Length distribution of pure non-repeat
    lens = [len(s) for _, s in pure_nonrepeat]
    if not lens:
        print("No pure non-repeat linkers found.")
        return
    lens_sorted = sorted(lens)
    print(f"\nPure non-repeat linker lengths:")
    print(f"  Range: {min(lens)}-{max(lens)} bp")
    print(f"  Median: {lens_sorted[len(lens)//2]}")
    print(f"  Mean: {sum(lens)/len(lens):.1f}")
    print(f"  Buckets:")
    for lo, hi, label in [(15, 30, 'short (NHEJ scar)'),
                           (30, 120, 'medium (TERC-frag zone)'),
                           (120, 1e9, 'long (captured frag)')]:
        n = sum(1 for x in lens if lo <= x < hi)
        print(f"    {label:30s} ({lo}-{hi if hi < 1e9 else 'inf'}): {n}")

    # Examine top medium-bucket pure non-repeats — these are best Nergadze candidates
    medium_pure = [(sid, seq) for sid, seq in pure_nonrepeat if 30 <= len(seq) <= 120]
    print(f"\nMedium pure non-repeat linkers (30-120bp): {len(medium_pure)}")

    if medium_pure:
        # Cluster these for shared template
        # Use simple kmer-set Jaccard
        sample = medium_pure[:300]
        # Build kmer sets
        sets = []
        for sid, seq in sample:
            ks = set()
            for i in range(len(seq) - 7):
                ks.add(seq[i:i+8])
            # Also rc
            rc = seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))
            for i in range(len(rc) - 7):
                ks.add(rc[i:i+8])
            sets.append((sid, seq, ks))

        # Find most "central" linker — one whose kmers overlap most others
        if len(sets) > 1:
            scores = []
            for i, (sid_i, seq_i, ks_i) in enumerate(sets):
                if not ks_i: continue
                overlaps = 0
                for j, (sid_j, seq_j, ks_j) in enumerate(sets):
                    if i == j: continue
                    inter = len(ks_i & ks_j)
                    if inter / min(len(ks_i), len(ks_j) or 1) >= 0.20:
                        overlaps += 1
                scores.append((overlaps, sid_i, seq_i))
            scores.sort(reverse=True)
            print(f"\nTop 5 'central' linkers (most kmer-similar to others) — TERC consensus candidates:")
            for ovl, sid, seq in scores[:5]:
                print(f"  Overlaps {ovl}/{len(sample)-1}:")
                print(f"    ID: {sid[:80]}")
                print(f"    Seq ({len(seq)}bp): {seq}")
                comp = composition(seq)
                print(f"    GC={comp['gc']:.2f}, A={comp['a']:.2f}, T={comp['t']:.2f}, G={comp['g']:.2f}, C={comp['c']:.2f}")
                print()

        # Most common 8-mers across pure non-repeat medium linkers
        kmer_in_n = Counter()
        for sid, seq in medium_pure:
            seen = set()
            for i in range(len(seq) - 7):
                seen.add(seq[i:i+8])
            for km in seen:
                kmer_in_n[km] += 1

        print(f"\nTop 15 most-shared 8-mers across {len(medium_pure)} pure non-repeat medium linkers:")
        for km, n in kmer_in_n.most_common(15):
            pct = 100 * n / len(medium_pure)
            print(f"  {km}: in {n} ({pct:.1f}%)")

    # Also write pure non-repeat linker FASTA for external BLAST
    out_fa = real_dir / "pure_nonrepeat_linkers.fa"
    with open(out_fa, 'w') as f:
        for sid, seq in pure_nonrepeat:
            f.write(f">{sid}\n{seq}\n")
    print(f"\nPure non-repeat linker FASTA: {out_fa} ({len(pure_nonrepeat)} sequences)")

    # Analyze composition vs length
    print(f"\nComposition vs length (pure non-repeat):")
    for lo, hi, label in [(15, 30, 'short'),
                           (30, 60, 'medium-short'),
                           (60, 120, 'medium-long'),
                           (120, 1e9, 'long')]:
        bucket = [(sid, seq) for sid, seq in pure_nonrepeat
                  if lo <= len(seq) < hi]
        if not bucket: continue
        gc = [composition(seq)['gc'] for _, seq in bucket]
        print(f"  {label:15s} ({lo}-{int(hi) if hi < 1e9 else 'inf'}bp): n={len(bucket)}, mean GC={sum(gc)/len(gc):.3f}")

    # Save summary
    out_json = real_dir / "true_nonrepeat_linkers_summary.json"
    summary = {
        'total_linkers': len(linkers),
        'pure_nonrepeat': len(pure_nonrepeat),
        'contaminated': len(contaminated),
        'low_telo': len(low_telo),
        'medium_pure_nonrepeat': len(medium_pure) if medium_pure else 0,
        'pure_nonrepeat_length_distribution': {
            'short_15_30': sum(1 for x in lens if 15 <= x < 30),
            'medium_short_30_60': sum(1 for x in lens if 30 <= x < 60),
            'medium_long_60_120': sum(1 for x in lens if 60 <= x < 120),
            'long_120plus': sum(1 for x in lens if x >= 120),
        }
    }
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary: {out_json}")


if __name__ == '__main__':
    main()
