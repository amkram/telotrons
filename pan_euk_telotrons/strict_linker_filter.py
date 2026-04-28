#!/usr/bin/env python3
"""
Stricter pure-non-repeat linker definition.

A sequence is "telomeric-degenerate" if ≥THRESH of its bases are covered by
EXACT or 1-mismatch (or 2-mismatch over 8+bp window) of any canonical
telomeric kmer. Such sequences are rejected as linkers.

Also rejects:
  - Sequences with > MAX_N_FRAC ambiguous Ns
  - Sequences shorter than MIN_LEN
"""
import csv, json, re
from pathlib import Path
from collections import defaultdict


def _rotations(b):
    return {b[i:] + b[:i] for i in range(len(b))}

FAMILIES = {
    'TTAGGG':  {'fwd': _rotations('TTAGGG'),  'rev': _rotations('CCCTAA')},
    'TTTAGGG': {'fwd': _rotations('TTTAGGG'), 'rev': _rotations('CCCTAAA')},
}

MAX_N_FRAC = 0.05
MIN_LEN = 25
TELO_FRAC_REJECT = 0.10
MAX_EXACT_TELO_KMERS = 1  # also reject if more than 1 exact telo kmer match


def telo_coverage_with_mismatches(seq, family, max_mismatches=1):
    """Fraction of bases covered by any telomeric kmer (exact or up to
    `max_mismatches` mismatches against any rotation/strand)."""
    if not seq: return 0.0
    fwd = FAMILIES[family]['fwd']
    rev = FAMILIES[family]['rev']
    all_kmers = list(fwd | rev)
    klen = len(family)
    n = len(seq)
    cov = bytearray(n)
    # Slide window
    for i in range(n - klen + 1):
        window = seq[i:i+klen]
        # Check exact + 1-mm
        for k in all_kmers:
            mm = sum(a != b for a, b in zip(window, k))
            if mm <= max_mismatches:
                for j in range(i, i + klen):
                    cov[j] = 1
                break
    return sum(cov) / n


def is_pure_nonrepeat(seq, family, telo_frac_reject=TELO_FRAC_REJECT):
    """True if the sequence is mostly non-telomeric:
       - low N content
       - <= MAX_EXACT_TELO_KMERS exact telo kmer matches
       - 1-mismatch coverage < telo_frac_reject
    """
    if not seq: return False
    if len(seq) < MIN_LEN: return False
    n_count = seq.count('N')
    if n_count / len(seq) > MAX_N_FRAC: return False

    # Exact match count
    fwd = FAMILIES[family]['fwd']
    rev = FAMILIES[family]['rev']
    n_exact = 0
    for k in fwd | rev:
        n_exact += seq.count(k)
    if n_exact > MAX_EXACT_TELO_KMERS: return False

    cov = telo_coverage_with_mismatches(seq, family, max_mismatches=1)
    if cov >= telo_frac_reject: return False

    return True


def find_runs(seq, kmers, min_run=10):
    n = len(seq); cov = bytearray(n)
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
            if i-rs >= min_run: runs.append((rs,i))
    if in_r and n-rs >= min_run: runs.append((rs,n))
    return runs


def extract_pure_linker(intron_seq, family):
    """Find the converging linker (fwd_block ... linker ... rev_block) and
    return (linker_seq, linker_start, linker_end) only if it passes
    is_pure_nonrepeat. Otherwise None."""
    sets = FAMILIES[family]
    fwd = find_runs(intron_seq, sets['fwd'])
    rev = find_runs(intron_seq, sets['rev'])
    if not (fwd and rev): return None

    fmax = max(fwd, key=lambda r: r[1] - r[0])
    rmax = max(rev, key=lambda r: r[1] - r[0])
    # Converging on genome strand: fwd before rev
    if fmax[0] < rmax[0]:
        gap_s, gap_e = fmax[1], rmax[0]
    else:
        gap_s, gap_e = rmax[1], fmax[0]
    if gap_e - gap_s < MIN_LEN: return None
    link = intron_seq[gap_s:gap_e]
    if not is_pure_nonrepeat(link, family): return None
    return link, gap_s, gap_e


def reextract_eimeria():
    """Re-extract Eimeria linkers under strict filter."""
    files = list(Path('ultra_results').glob('*Eimeria*.tsv'))
    eim = []
    seen = set()
    for f in files:
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                seq = row.get('intron_seq', '').upper()
                if not seq: continue
                # Filter intron-level: drop if mostly N
                n_frac = seq.count('N') / len(seq)
                if n_frac > MAX_N_FRAC: continue
                key = (acc, row['contig'], int(row['start']), int(row['end']))
                if key in seen: continue
                seen.add(key)
                result = extract_pure_linker(seq, 'TTTAGGG')
                if not result: continue
                link, ls, le = result
                eim.append({
                    'acc': acc,
                    'contig': row['contig'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'linker_seq': link,
                    'linker_offset': ls,
                    'linker_end_offset': le,
                    'linker_genomic_start': int(row['start']) + ls,
                    'linker_genomic_end': int(row['start']) + le,
                    'family': 'TTTAGGG',
                    'intron_seq': seq,
                })
    return eim


def reextract_tara():
    """Re-extract Tara linkers under strict filter."""
    tara = []
    seen = set()
    with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            seq = row.get('intron_seq', '').upper()
            if not seq: continue
            n_frac = seq.count('N') / len(seq)
            if n_frac > MAX_N_FRAC: continue
            try:
                ts, te = int(row['start']), int(row['end'])
            except: continue
            key = (row['mag'], row['contig'], ts, te)
            if key in seen: continue
            seen.add(key)
            result = extract_pure_linker(seq, 'TTAGGG')
            if not result: continue
            link, ls, le = result
            tara.append({
                'mag': row['mag'],
                'contig': row['contig'],
                'start': ts, 'end': te,
                'linker_seq': link,
                'linker_offset': ls,
                'linker_end_offset': le,
                'linker_genomic_start': ts + ls,
                'linker_genomic_end': ts + le,
                'family': 'TTAGGG',
                'intron_seq': seq,
            })
    return tara


if __name__ == '__main__':
    eim = reextract_eimeria()
    print(f"Eimeria pure non-repeat linkers (strict): {len(eim)}")
    by_acc = defaultdict(int)
    for r in eim: by_acc[r['acc']] += 1
    for a, n in sorted(by_acc.items()): print(f"  {a}: {n}")

    tara = reextract_tara()
    print(f"\nTara pure non-repeat linkers (strict): {len(tara)}")
    by_mag = defaultdict(int)
    for r in tara: by_mag[r['mag']] += 1
    for m, n in sorted(by_mag.items(), key=lambda x: -x[1]): print(f"  {m}: {n}")

    with open('real_telotrons/strict_linkers_eimeria.json', 'w') as f:
        json.dump(eim, f, indent=1)
    with open('real_telotrons/strict_linkers_tara.json', 'w') as f:
        json.dump(tara, f, indent=1)
    print("\nWrote strict_linkers_eimeria.json and strict_linkers_tara.json")
