#!/usr/bin/env python3
"""
Deep analysis of Eimeria linkers.

Eimeria differs dramatically from haptophyte telotrons:
  - Mostly DIVERGING (5'-CCCTAAA...TTTAGGG-3' on genome strand)
  - Long linkers (median 116bp vs 44 in haptophytes)
  - 27% truly non-repeat (vs 8.6%)

Tests:
  1. BLAST Eimeria pure non-repeat linkers against:
     - Eimeria genome itself (any genomic origin?)
     - Each other (shared template like TERC?)
  2. Composition vs random Eimeria intron
  3. Inverted-repeat structure: do flanks form palindrome?
  4. Compare CONVERGING vs DIVERGING linker properties
"""

import csv
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

FAMILIES = {
    'TTTAGGG': {'fwd': _rotations('TTTAGGG'), 'rev': _rotations('CCCTAAA')},
    'TTAGGG':  {'fwd': _rotations('TTAGGG'),  'rev': _rotations('CCCTAA')},
}
ALL_TELO = (FAMILIES['TTTAGGG']['fwd'] | FAMILIES['TTTAGGG']['rev']
            | FAMILIES['TTAGGG']['fwd'] | FAMILIES['TTAGGG']['rev'])


def find_runs(seq, kmers, min_run=10):
    n = len(seq)
    cov = bytearray(n)
    for k in kmers:
        klen = len(k); i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1: break
            for j in range(idx, min(idx+klen, n)): cov[j] = 1
            i = idx + 1
    runs = []; in_r = False; rs = 0
    for i in range(n):
        if cov[i] and not in_r: in_r=True; rs=i
        elif not cov[i] and in_r:
            in_r=False
            if i-rs >= min_run: runs.append((rs, i))
    if in_r and n-rs >= min_run: runs.append((rs, n))
    return runs


def composition(seq):
    n = len(seq)
    c = Counter(seq)
    return {
        'gc': (c.get('G', 0) + c.get('C', 0)) / n if n else 0,
        'len': n,
    }


def has_telo_kmer(seq):
    return any(k in seq for k in ALL_TELO)


def main():
    files = list(Path('telotrons').glob('*Eimeria*.tsv')) + list(Path('ultra_results').glob('*Eimeria*.tsv'))
    print(f"Reading {len(files)} Eimeria telotron files\n")

    converging = []  # 5'-fwd...rev-3' on genome strand
    diverging = []   # 5'-rev...fwd-3'
    seen_loci = set()  # dedupe across v1/ultra

    for f in files:
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                seq = row.get('intron_seq', '').upper()
                if not seq: continue
                key = (row['contig'], row['start'], row['end'])
                if key in seen_loci: continue
                seen_loci.add(key)

                # Use TTTAGGG family
                fwd_runs = find_runs(seq, FAMILIES['TTTAGGG']['fwd'])
                rev_runs = find_runs(seq, FAMILIES['TTTAGGG']['rev'])
                if not (fwd_runs and rev_runs): continue

                fwd_max = max(fwd_runs, key=lambda r: r[1]-r[0])
                rev_max = max(rev_runs, key=lambda r: r[1]-r[0])

                if fwd_max[0] < rev_max[0]:
                    cls = 'CONVERGING'
                    linker_start = fwd_max[1]
                    linker_end = rev_max[0]
                else:
                    cls = 'DIVERGING'
                    linker_start = rev_max[1]
                    linker_end = fwd_max[0]

                gap = linker_end - linker_start
                if gap <= 4:
                    cls += '_NOLINKER'
                else:
                    cls += '_WITHLINKER'

                rec = {
                    'contig': row['contig'], 'start': row['start'], 'end': row['end'],
                    'length': int(row['length']),
                    'strand': row.get('strand', ''),
                    'donor': row.get('donor', ''),
                    'acceptor': row.get('acceptor', ''),
                    'gene': row.get('mrna_id', '') or row.get('gene_name', ''),
                    'cls': cls,
                    'fwd_run': fwd_max,
                    'rev_run': rev_max,
                    'gap': gap,
                    'linker_start': linker_start,
                    'linker_end': linker_end,
                    'linker': seq[linker_start:linker_end] if linker_end > linker_start else '',
                    'seq': seq,
                }
                if cls.startswith('CONVERGING'):
                    converging.append(rec)
                elif cls.startswith('DIVERGING'):
                    diverging.append(rec)

    print(f"CONVERGING: {len(converging)}")
    print(f"DIVERGING:  {len(diverging)}")
    print(f"Ratio:      1:{len(diverging)/max(1,len(converging)):.1f}")

    # Restrict to with-linker
    cwl = [r for r in converging if r['cls'] == 'CONVERGING_WITHLINKER']
    dwl = [r for r in diverging if r['cls'] == 'DIVERGING_WITHLINKER']

    print(f"\nWith-linker subsets:")
    print(f"  CONVERGING_WITHLINKER: {len(cwl)}")
    print(f"  DIVERGING_WITHLINKER:  {len(dwl)}")

    # Linker composition
    print(f"\n--- LINKER PROPERTIES ---")
    for label, recs in [('CONVERGING_WITHLINKER', cwl), ('DIVERGING_WITHLINKER', dwl)]:
        if not recs: continue
        lens = [len(r['linker']) for r in recs]
        gcs = [composition(r['linker'])['gc'] for r in recs if r['linker']]
        pure = sum(1 for r in recs if r['linker'] and not has_telo_kmer(r['linker']))
        print(f"\n  {label} (n={len(recs)}):")
        print(f"    Length range: {min(lens)}-{max(lens)}, median {sorted(lens)[len(lens)//2]}")
        print(f"    Mean GC: {sum(gcs)/len(gcs):.3f}")
        print(f"    Pure non-repeat (no TTAGGG/TTTAGGG/CCCTAA/CCCTAAA kmers): {pure}/{len(recs)} ({100*pure/len(recs):.1f}%)")
        print(f"    Length buckets:")
        for lo, hi, name in [(5, 30, 'short'), (30, 120, 'medium'), (120, 1e9, 'long')]:
            n = sum(1 for x in lens if lo <= x < hi)
            print(f"      {name} ({lo}-{int(hi) if hi<1e9 else 'inf'}): {n} ({100*n/len(recs):.1f}%)")

    # Save pure non-repeat linkers per class for BLAST
    out_dir = Path('real_telotrons')
    for label, recs in [('cwl', cwl), ('dwl', dwl)]:
        out_fa = out_dir / f"eimeria_{label}_pure_linkers.fa"
        n = 0
        with open(out_fa, 'w') as f:
            for r in recs:
                if r['linker'] and len(r['linker']) >= 20 and not has_telo_kmer(r['linker']):
                    f.write(f">eimeria__{r['contig']}_{r['start']}__{r['cls']}__len{len(r['linker'])}\n")
                    f.write(f"{r['linker']}\n")
                    n += 1
        print(f"\n  Wrote {n} pure non-repeat linkers (≥20bp) to {out_fa}")

    # Inverted-repeat / palindrome test on DIVERGING flanks
    print(f"\n--- INVERTED-REPEAT TEST (Diverging hypothesis: captured fragment in inverted orientation) ---")
    if dwl:
        # In DIVERGING (5'-rev_max-linker-fwd_max-3'), check if the rev block and fwd block are
        # inverted repeats (i.e., one is rev-comp of the other, allowing for typical telomere variation)
        # They SHOULD be related (both are TTTAGGG-family) but we want to check exact correspondence
        # of the flanks of the linker.
        # Simpler: check if linker contains an inverted repeat structure spanning a palindrome center

        inverted_count = 0
        for r in dwl[:200]:
            linker = r['linker']
            if len(linker) < 30: continue
            # Check first half vs reverse-complement of second half
            half = len(linker) // 2
            first = linker[:half]
            second_rc = linker[-half:][::-1].translate(str.maketrans('ACGT', 'TGCA'))
            matches = sum(a == b for a, b in zip(first, second_rc))
            if matches / half >= 0.7:
                inverted_count += 1
        print(f"  Linkers with ≥70% palindromic identity (first half ≈ rev-comp of last half):")
        print(f"    {inverted_count} / {min(200, len(dwl))} ({100*inverted_count/min(200,len(dwl)):.1f}%)")

    # Compare CONVERGING vs DIVERGING linker GC%
    if cwl and dwl:
        cwl_gcs = [composition(r['linker'])['gc'] for r in cwl if r['linker']]
        dwl_gcs = [composition(r['linker'])['gc'] for r in dwl if r['linker']]
        print(f"\n  Linker GC%:")
        print(f"    CONVERGING: mean {sum(cwl_gcs)/len(cwl_gcs):.3f}, median {sorted(cwl_gcs)[len(cwl_gcs)//2]:.3f}, n={len(cwl_gcs)}")
        print(f"    DIVERGING:  mean {sum(dwl_gcs)/len(dwl_gcs):.3f}, median {sorted(dwl_gcs)[len(dwl_gcs)//2]:.3f}, n={len(dwl_gcs)}")

    # Save summary
    summary = {
        'n_converging': len(converging),
        'n_diverging': len(diverging),
        'n_cwl': len(cwl),
        'n_dwl': len(dwl),
        'cwl_pure_pct': 100 * sum(1 for r in cwl if r['linker'] and not has_telo_kmer(r['linker'])) / len(cwl) if cwl else 0,
        'dwl_pure_pct': 100 * sum(1 for r in dwl if r['linker'] and not has_telo_kmer(r['linker'])) / len(dwl) if dwl else 0,
    }
    out_json = out_dir / 'eimeria_linker_deep.json'
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary: {out_json}")


if __name__ == '__main__':
    main()
