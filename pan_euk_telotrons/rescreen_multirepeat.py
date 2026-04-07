#!/usr/bin/env python3
"""
Re-screen all telotron candidates with species-appropriate telomeric repeats.

For each intron sequence in telotrons/*.tsv, test all known repeat types
and assign the best-matching one. Recalculate purity using the correct
repeat unit. Output updated TSVs.
"""
import csv, re, json, os
from pathlib import Path
from collections import Counter, defaultdict

def _rotations(base):
    """All rotations of base + its reverse complement."""
    s = set()
    for i in range(len(base)):
        s.add(base[i:] + base[:i])
    rc = base[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    for i in range(len(rc)):
        s.add(rc[i:] + rc[:i])
    return s

def _fwd_rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

def _rev_rotations(base):
    rc = base[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    return {rc[i:] + rc[:i] for i in range(len(rc))}

REPEAT_TYPES = {
    'TTAGGG':   {'all': _rotations('TTAGGG'),   'fwd': _fwd_rotations('TTAGGG'),   'rev': _rev_rotations('TTAGGG')},
    'TTTAGGG':  {'all': _rotations('TTTAGGG'),  'fwd': _fwd_rotations('TTTAGGG'),  'rev': _rev_rotations('TTTAGGG')},
    'TTTGGG':   {'all': _rotations('TTTGGG'),   'fwd': _fwd_rotations('TTTGGG'),   'rev': _rev_rotations('TTTGGG')},
    'TTGGGG':   {'all': _rotations('TTGGGG'),   'fwd': _fwd_rotations('TTGGGG'),   'rev': _rev_rotations('TTGGGG')},
    'TTAGG':    {'all': _rotations('TTAGG'),    'fwd': _fwd_rotations('TTAGG'),    'rev': _rev_rotations('TTAGG')},
    'TTTTAGGG': {'all': _rotations('TTTTAGGG'), 'fwd': _fwd_rotations('TTTTAGGG'), 'rev': _rev_rotations('TTTTAGGG')},
    'TTAGGC':   {'all': _rotations('TTAGGC'),   'fwd': _fwd_rotations('TTAGGC'),   'rev': _rev_rotations('TTAGGC')},
}

def coverage(seq, hexset):
    n = len(seq)
    if n < 4: return 0.0
    covered = bytearray(n)
    for h in hexset:
        start = 0
        while True:
            idx = seq.find(h, start)
            if idx == -1: break
            for j in range(idx, min(idx + len(h), n)):
                covered[j] = 1
            start = idx + 1
    return sum(covered) / n

def strand_ratio(seq, fwd_set, rev_set):
    fwd = sum(len(re.findall(re.escape(h), seq)) for h in fwd_set)
    rev = sum(len(re.findall(re.escape(h), seq)) for h in rev_set)
    total = fwd + rev
    if total == 0: return 0.5, 'unknown'
    frac = fwd / total
    if frac > 0.7: return frac, 'coding-strand'
    elif frac < 0.3: return frac, 'template-strand'
    else: return frac, 'converging'

def rc(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))

# Process all telotron TSVs
telo_dir = Path('telotrons')
out_dir = Path('telotrons_multirepeat')
out_dir.mkdir(exist_ok=True)

total_in = 0
total_out = 0
species_summary = []

for tsv in sorted(telo_dir.glob('*_telotrons.tsv')):
    rows_in = []
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        rows_in = list(reader)

    if not rows_in:
        continue

    total_in += len(rows_in)
    rows_out = []

    for row in rows_in:
        seq = row.get('intron_seq', '').upper()
        if not seq: continue

        # Find best repeat type
        best_type = None
        best_cov = 0
        for rtype, rsets in REPEAT_TYPES.items():
            cov = coverage(seq, rsets['all'])
            if cov > best_cov:
                best_cov = cov
                best_type = rtype

        if best_cov < 0.85:
            continue

        # Recalculate orientation with correct repeat
        rsets = REPEAT_TYPES[best_type]
        frac, orient = strand_ratio(seq, rsets['fwd'], rsets['rev'])

        # Normalize splice sites
        d2 = seq[:2]; a2 = seq[-2:]
        norm_seq = seq
        if d2 == 'CT' and a2 == 'AC':
            norm_seq = rc(seq)

        donor_8 = norm_seq[:8] if len(norm_seq) >= 8 else norm_seq
        acceptor_8 = norm_seq[-8:] if len(norm_seq) >= 8 else norm_seq

        row_out = dict(row)
        row_out['best_repeat'] = best_type
        row_out['telo_coverage'] = f"{best_cov:.4f}"
        row_out['fwd_frac'] = f"{frac:.4f}"
        row_out['orientation'] = orient
        row_out['donor_8'] = donor_8
        row_out['acceptor_8'] = acceptor_8
        rows_out.append(row_out)

    if rows_out:
        # Write updated TSV
        out_path = out_dir / tsv.name
        out_headers = [h for h in headers if h != 'best_repeat']
        # Insert best_repeat after telo_coverage
        idx = out_headers.index('telo_coverage') + 1 if 'telo_coverage' in out_headers else len(out_headers)
        out_headers.insert(idx, 'best_repeat')

        with open(out_path, 'w') as f:
            f.write('\t'.join(out_headers) + '\n')
            for r in rows_out:
                f.write('\t'.join(str(r.get(h, '')) for h in out_headers) + '\n')

        total_out += len(rows_out)

        # Summary
        repeat_counts = Counter(r['best_repeat'] for r in rows_out)
        orient_counts = Counter(r['orientation'] for r in rows_out)
        species_name = rows_out[0].get('genome_name', tsv.stem)
        species_summary.append({
            'species': species_name,
            'n_in': len(rows_in),
            'n_out': len(rows_out),
            'repeats': dict(repeat_counts),
            'orientations': dict(orient_counts),
        })

# Print summary
print(f"{'='*70}")
print(f"MULTI-REPEAT RE-SCREENING RESULTS")
print(f"{'='*70}")
print(f"Input:  {total_in} candidates across {len(list(telo_dir.glob('*_telotrons.tsv')))} species")
print(f"Output: {total_out} at ≥85% with correct repeat ({out_dir}/)")
print()

# Overall repeat distribution
all_repeats = Counter()
all_orients = Counter()
for s in species_summary:
    for r, n in s['repeats'].items(): all_repeats[r] += n
    for o, n in s['orientations'].items(): all_orients[o] += n

print("Repeat type distribution:")
for r, n in all_repeats.most_common():
    print(f"  {r:>10s}: {n:>5d}")

print(f"\nOrientation distribution:")
for o, n in all_orients.most_common():
    print(f"  {o:>15s}: {n:>5d}")

print(f"\nPer-species breakdown:")
print(f"{'Species':<45s} {'N':>4s} {'Repeat':>10s} {'T':>3s} {'C':>3s} {'X':>3s}")
print('-' * 75)
for s in sorted(species_summary, key=lambda x: -x['n_out']):
    dominant = max(s['repeats'], key=s['repeats'].get)
    o = s['orientations']
    print(f"  {s['species']:<43s} {s['n_out']:>4d} {dominant:>10s} "
          f"{o.get('template-strand',0):>3d} {o.get('coding-strand',0):>3d} "
          f"{o.get('converging',0):>3d}")
