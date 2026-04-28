#!/usr/bin/env python3
"""
ITS-mechanism diagnostic battery — v2 with Tara Oceans data.

Input:
  - real_telotrons/*_real_telotrons.tsv (pan-euk validated, 186)
  - real_telotrons/_tara_oceans_telotrons.tsv (Tara, 55822 with ≥30% telo)

Output:
  - real_telotrons/its_mechanism_v2.tsv     (per-telotron classification)
  - real_telotrons/its_mechanism_v2.json    (summary stats)
  - real_telotrons/converging_linkers.fa    (for BLAST)
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

CANONICAL_FAMILIES = {
    'TTAGGG':   {'fwd': _rotations('TTAGGG'),   'rev': _rotations('CCCTAA')},
    'TTTAGGG':  {'fwd': _rotations('TTTAGGG'),  'rev': _rotations('CCCTAAA')},
    'TTAGG':    {'fwd': _rotations('TTAGG'),    'rev': _rotations('CCTAA')},
    'TTTGGG':   {'fwd': _rotations('TTTGGG'),   'rev': _rotations('CCCAAA')},
    'TTGGGG':   {'fwd': _rotations('TTGGGG'),   'rev': _rotations('CCCCAA')},
}

# TTI variants (Marzec 2015): variants that mark NR2C/F-bound telomeres
TTI_VARIANTS_FWD = {'TCAGGG', 'TGAGGG', 'GTAGGG', 'TTGGGG', 'TTCGGG', 'CTAGGG'}
TTI_VARIANTS_REV = {''.join({'A':'T','T':'A','G':'C','C':'G'}[b] for b in s[::-1])
                    for s in TTI_VARIANTS_FWD}

# NR2C/F binding half-site
NR2CF_HALFSITE = 'GGGTCA'
NR2CF_HALFSITE_RC = 'TGACCC'

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def rc(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def find_runs(seq, kmers, min_run_bp=10):
    """Return list of (start, end) where any kmer in kmers covers ≥min_run_bp contiguous bp."""
    n = len(seq)
    if n < min_run_bp:
        return []
    cov = bytearray(n)
    for k in kmers:
        klen = len(k)
        i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1:
                break
            for j in range(idx, min(idx + klen, n)):
                cov[j] = 1
            i = idx + 1
    runs = []
    in_run = False
    rs = 0
    for i in range(n):
        if cov[i] and not in_run:
            in_run = True; rs = i
        elif not cov[i] and in_run:
            in_run = False
            if i - rs >= min_run_bp:
                runs.append((rs, i))
    if in_run and n - rs >= min_run_bp:
        runs.append((rs, n))
    return runs


def detect_canonical_family(seq):
    best = (None, 0, 0, 0, 0)
    for fam, sets in CANONICAL_FAMILIES.items():
        fwd_runs = find_runs(seq, sets['fwd'], 10)
        rev_runs = find_runs(seq, sets['rev'], 10)
        fwd_bp = sum(e - s for s, e in fwd_runs)
        rev_bp = sum(e - s for s, e in rev_runs)
        fwd_run = max((e - s for s, e in fwd_runs), default=0)
        rev_run = max((e - s for s, e in rev_runs), default=0)
        score = fwd_bp + rev_bp
        if score > best[1] + best[2]:
            best = (fam, fwd_bp, rev_bp, fwd_run, rev_run)
    return best


def classify_orientation(seq, family):
    sets = CANONICAL_FAMILIES[family]
    fwd_runs = find_runs(seq, sets['fwd'], 10)
    rev_runs = find_runs(seq, sets['rev'], 10)
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

    fwd_max = max(fwd_runs, key=lambda r: r[1] - r[0])
    rev_max = max(rev_runs, key=lambda r: r[1] - r[0])

    if fwd_max[0] < rev_max[0]:
        first, second = fwd_max, rev_max
        order = 'G_then_C'   # head-to-head 5'-(TTAGGG)n-(CCCTAA)m-3'
    else:
        first, second = rev_max, fwd_max
        order = 'C_then_G'   # tail-to-tail diverging

    gap_start = first[1]
    gap_end = second[0]
    linker_len = max(0, gap_end - gap_start)
    linker = seq[gap_start:gap_end] if linker_len > 0 else ''

    if order == 'G_then_C':
        if linker_len <= 4:
            return 'CONVERGING_NOLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
        return 'CONVERGING_WITHLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
    else:
        if linker_len <= 4:
            return 'DIVERGING_NOLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)
        return 'DIVERGING_WITHLINKER', fwd_runs, rev_runs, linker, (gap_start, gap_end)


def variant_repeat_signal(seq, family):
    """For TTAGGG family: count canonical and 1-mismatch variant hexamers."""
    if family != 'TTAGGG':
        return None
    canonical = CANONICAL_FAMILIES['TTAGGG']['fwd'] | CANONICAL_FAMILIES['TTAGGG']['rev']
    n_can = 0; n_var = 0; n_total = 0
    tti_var_count = 0
    n = len(seq)
    for i in range(n - 5):
        kmer = seq[i:i+6]
        n_total += 1
        if kmer in canonical:
            n_can += 1
            continue
        # 1-mismatch from any canonical TTAGGG-family kmer?
        is_var = False
        for ck in canonical:
            mm = sum(a != b for a, b in zip(kmer, ck))
            if mm == 1:
                n_var += 1
                is_var = True
                break
        if is_var and (kmer in TTI_VARIANTS_FWD or kmer in TTI_VARIANTS_REV):
            tti_var_count += 1
    return {
        'canonical': n_can,
        'variants_1mm': n_var,
        'tti_specific_variants': tti_var_count,
        'total_windows': n_total,
        'canonical_frac': n_can / n_total if n_total else 0,
        'variant_frac': n_var / n_total if n_total else 0,
        'tti_frac': tti_var_count / n_total if n_total else 0,
    }


def microhomology_at_5p_boundary(left_flank, intron_seq, family_kmers, max_match=8):
    """How many bp of the upstream exon's 3' end match a canonical telo kmer?"""
    best = 0
    if not left_flank: return 0
    # Slice last N bp of exon
    for L in range(max_match, 1, -1):
        if L > len(left_flank): continue
        tail = left_flank[-L:]
        # Does this tail appear at the START of any rotation of canonical?
        for k in family_kmers:
            for shift in range(len(k)):
                rotated = k[shift:] + k[:shift]
                if rotated.startswith(tail):
                    return L
    return best


def microhomology_at_3p_boundary(intron_seq, right_flank, family_kmers, max_match=8):
    best = 0
    if not right_flank: return 0
    for L in range(max_match, 1, -1):
        if L > len(right_flank): continue
        head = right_flank[:L]
        for k in family_kmers:
            for shift in range(len(k)):
                rotated = k[shift:] + k[:shift]
                if rotated.endswith(head):
                    return L
    return best


def find_tsd(left, right, min_len=3, max_len=25):
    if not left or not right: return (0, '')
    for L in range(min(max_len, len(left), len(right)), min_len - 1, -1):
        if left[-L:] == right[:L]:
            return (L, left[-L:])
    return (0, '')


def scan_nr2cf(seq, max_spacing=10):
    if not seq: return {'half_sites_fwd': 0, 'half_sites_rc': 0,
                       'direct_repeats': 0, 'inverted_repeats': 0}
    L = len(NR2CF_HALFSITE)
    fwd_pos = [i for i in range(len(seq) - L + 1) if seq[i:i+L] == NR2CF_HALFSITE]
    rc_pos = [i for i in range(len(seq) - L + 1) if seq[i:i+L] == NR2CF_HALFSITE_RC]
    direct_repeats = 0
    inverted_repeats = 0
    for i in range(len(fwd_pos)):
        for j in range(i+1, len(fwd_pos)):
            spacing = fwd_pos[j] - fwd_pos[i] - L
            if 0 <= spacing <= max_spacing:
                direct_repeats += 1
    for i in fwd_pos:
        for j in rc_pos:
            spacing = abs(j - i) - L
            if 0 <= spacing <= max_spacing:
                inverted_repeats += 1
    return {
        'half_sites_fwd': len(fwd_pos),
        'half_sites_rc': len(rc_pos),
        'direct_repeats': direct_repeats,
        'inverted_repeats': inverted_repeats,
    }


# ------------------------------------------------------------------
# Per-record analysis
# ------------------------------------------------------------------

def analyze(record):
    seq = record.get('intron_seq', '').upper()
    left = record.get('left_flank_100bp', '').upper()
    right = record.get('right_flank_100bp', '').upper()

    if not seq:
        return None

    fam, fwd_bp, rev_bp, fwd_run, rev_run = detect_canonical_family(seq)
    if fam is None:
        return None

    cls, fwd_runs, rev_runs, linker, linker_pos = classify_orientation(seq, fam)
    var_info = variant_repeat_signal(seq, fam)

    fam_kmers = CANONICAL_FAMILIES[fam]['fwd'] | CANONICAL_FAMILIES[fam]['rev']
    mh5 = microhomology_at_5p_boundary(left, seq, fam_kmers)
    mh3 = microhomology_at_3p_boundary(seq, right, fam_kmers)

    tsd_len, tsd_seq = find_tsd(left, right)
    nr2cf = scan_nr2cf(left + right)

    return {
        'family': fam,
        'classification': cls,
        'fwd_bp': fwd_bp,
        'rev_bp': rev_bp,
        'fwd_run_max': fwd_run,
        'rev_run_max': rev_run,
        'n_fwd_runs': len(fwd_runs),
        'n_rev_runs': len(rev_runs),
        'linker': linker,
        'linker_len': len(linker) if linker else 0,
        'linker_gc': (linker.count('G') + linker.count('C')) / len(linker) if linker else 0,
        'variant_info': var_info,
        'mh_5p': mh5,
        'mh_3p': mh3,
        'tsd_len': tsd_len,
        'tsd_seq': tsd_seq,
        'nr2cf': nr2cf,
    }


# ------------------------------------------------------------------
# Loaders
# ------------------------------------------------------------------

def load_panieuk():
    real_dir = Path(__file__).parent / "real_telotrons"
    records = []
    for tsv in sorted(real_dir.glob("*_real_telotrons.tsv")):
        if tsv.name.startswith('_'):  # skip _tara_oceans_telotrons.tsv
            continue
        with open(tsv) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            for row in rdr:
                if row.get('intron_seq'):
                    row['_dataset'] = 'pan_euk'
                    row['_id'] = f"{row['species']}__{row['contig']}_{row['start']}"
                    records.append(row)
    return records


def load_tara():
    real_dir = Path(__file__).parent / "real_telotrons"
    tara_tsv = real_dir / "_tara_oceans_telotrons.tsv"
    records = []
    if not tara_tsv.exists():
        return records
    with open(tara_tsv) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            if row.get('intron_seq'):
                row['_dataset'] = 'tara'
                row['_id'] = f"{row['mag']}__{row['contig']}_{row['start']}"
                records.append(row)
    return records


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    real_dir = Path(__file__).parent / "real_telotrons"

    pan = load_panieuk()
    tara = load_tara()
    all_records = pan + tara
    print(f"Loaded {len(pan)} pan-euk + {len(tara)} Tara = {len(all_records)} total telotron candidates")

    # Run analysis
    analyzed = []
    for i, rec in enumerate(all_records):
        if i % 5000 == 0:
            print(f"  Analyzed {i}/{len(all_records)}", flush=True)
        a = analyze(rec)
        if a is None: continue
        a['dataset'] = rec['_dataset']
        a['id'] = rec['_id']
        a['length'] = int(rec.get('length', len(rec.get('intron_seq', ''))))
        a['telo_coverage'] = float(rec.get('telo_coverage', 0))
        analyzed.append(a)

    print(f"  Analyzed {len(analyzed)} records\n")

    # Class distributions per dataset
    by_ds_class = defaultdict(lambda: defaultdict(int))
    by_ds_class_purity = defaultdict(lambda: defaultdict(list))
    for r in analyzed:
        by_ds_class[r['dataset']][r['classification']] += 1
        by_ds_class_purity[r['dataset']][r['classification']].append(r['telo_coverage'])

    # Linker length distribution per dataset/class
    linker_lens = defaultdict(list)
    for r in analyzed:
        if r['linker_len'] >= 5:
            linker_lens[r['classification']].append(r['linker_len'])

    # TSD rates
    tsd_rates = defaultdict(list)
    for r in analyzed:
        tsd_rates[r['classification']].append(r['tsd_len'])

    # NR2C/F per class
    nr2cf_rates = defaultdict(list)
    for r in analyzed:
        nr2cf_rates[r['classification']].append(
            r['nr2cf']['half_sites_fwd'] + r['nr2cf']['half_sites_rc']
        )

    # Variant fraction (TTI signal) per class — only TTAGGG family
    var_frac = defaultdict(list)
    tti_frac = defaultdict(list)
    for r in analyzed:
        if r['family'] == 'TTAGGG' and r['variant_info']:
            var_frac[r['classification']].append(r['variant_info']['variant_frac'])
            tti_frac[r['classification']].append(r['variant_info']['tti_frac'])

    # Microhomology
    mh = defaultdict(lambda: {'mh5': [], 'mh3': []})
    for r in analyzed:
        mh[r['classification']]['mh5'].append(r['mh_5p'])
        mh[r['classification']]['mh3'].append(r['mh_3p'])

    # Print summary
    print(f"{'='*72}")
    print(f"  ITS-MECHANISM CLASSIFICATION (v2: pan-euk + Tara Oceans)")
    print(f"{'='*72}\n")
    print(f"Class distribution by dataset:")
    print(f"  {'Class':<28s} {'pan_euk':>10s} {'tara':>10s} {'total':>10s}")
    print(f"  {'-'*28} {'-'*10} {'-'*10} {'-'*10}")
    classes = sorted(set(c for d in by_ds_class.values() for c in d))
    for c in classes:
        n_pan = by_ds_class['pan_euk'].get(c, 0)
        n_tara = by_ds_class['tara'].get(c, 0)
        print(f"  {c:<28s} {n_pan:>10,} {n_tara:>10,} {n_pan+n_tara:>10,}")

    print(f"\nMean telomeric coverage per class (Tara only, ≥30% threshold):")
    for c in sorted(by_ds_class_purity['tara'].keys()):
        v = by_ds_class_purity['tara'][c]
        if v:
            print(f"  {c:<28s} mean={sum(v)/len(v):.3f}  median={sorted(v)[len(v)//2]:.3f}  n={len(v)}")

    print(f"\nLinker length distribution (CONVERGING_WITHLINKER class):")
    for c in ['CONVERGING_WITHLINKER', 'DIVERGING_WITHLINKER']:
        ll = linker_lens.get(c, [])
        if ll:
            ll_sorted = sorted(ll)
            print(f"  {c}:")
            print(f"    n={len(ll)}, range {min(ll)}-{max(ll)}, median={ll_sorted[len(ll)//2]}")
            print(f"    short  (5-30bp,    NHEJ-scar):  {sum(1 for x in ll if 5 <= x <= 30):>5}  ({100*sum(1 for x in ll if 5 <= x <= 30)/len(ll):.1f}%)")
            print(f"    medium (30-120bp,  TERC-frag):  {sum(1 for x in ll if 30 < x <= 120):>5}  ({100*sum(1 for x in ll if 30 < x <= 120)/len(ll):.1f}%)")
            print(f"    long   (120+bp,    captured):   {sum(1 for x in ll if x > 120):>5}  ({100*sum(1 for x in ll if x > 120)/len(ll):.1f}%)")

    print(f"\nTSD presence per class:")
    for c in sorted(tsd_rates.keys()):
        v = tsd_rates[c]
        if v:
            n_with = sum(1 for x in v if x >= 3)
            mean_when_present = sum(x for x in v if x >= 3) / max(1, n_with)
            print(f"  {c:<28s} {100*n_with/len(v):>5.1f}%  (n={len(v)}, mean TSD len when present: {mean_when_present:.1f}bp)")

    print(f"\nMicrohomology to canonical TTAGGG at 5'/3' boundaries:")
    for c in sorted(mh.keys()):
        mh5 = mh[c]['mh5']
        mh3 = mh[c]['mh3']
        if mh5 and mh3:
            print(f"  {c:<28s} mean5'={sum(mh5)/len(mh5):.2f}bp, mean3'={sum(mh3)/len(mh3):.2f}bp")

    print(f"\nVariant repeat fraction (1-mismatch from TTAGGG, TTAGGG-family only) — TTI marker:")
    print(f"  Higher = more variant-rich = TTI/ALT-like signature")
    for c in sorted(var_frac.keys()):
        v = var_frac[c]
        t = tti_frac[c]
        if v:
            print(f"  {c:<28s} variant_frac mean={sum(v)/len(v):.4f}, TTI-specific mean={sum(t)/len(t):.5f}, n={len(v)}")

    print(f"\nMean NR2C/F GGGTCA half-sites in 200bp flanks:")
    print(f"  (Higher = more potential TTI binding)")
    for c in sorted(nr2cf_rates.keys()):
        v = nr2cf_rates[c]
        if v:
            n_with_any = sum(1 for x in v if x > 0)
            print(f"  {c:<28s} mean={sum(v)/len(v):.3f} ({100*n_with_any/len(v):.1f}% have ≥1)")

    # Write outputs
    out_tsv = real_dir / "its_mechanism_v2.tsv"
    cols = ['dataset', 'id', 'length', 'telo_coverage', 'family', 'classification',
            'fwd_bp', 'rev_bp', 'fwd_run_max', 'rev_run_max',
            'n_fwd_runs', 'n_rev_runs',
            'linker_len', 'linker_gc',
            'canonical_frac', 'variant_frac', 'tti_frac',
            'mh_5p', 'mh_3p', 'tsd_len', 'tsd_seq',
            'nr2cf_halfsites', 'nr2cf_direct_repeats', 'nr2cf_inverted_repeats',
            'linker']
    with open(out_tsv, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in analyzed:
            vi = r.get('variant_info') or {}
            ni = r.get('nr2cf') or {}
            row = [
                r['dataset'], r['id'], r['length'], f"{r['telo_coverage']:.4f}",
                r['family'], r['classification'],
                r['fwd_bp'], r['rev_bp'], r['fwd_run_max'], r['rev_run_max'],
                r['n_fwd_runs'], r['n_rev_runs'],
                r['linker_len'], f"{r['linker_gc']:.3f}",
                f"{vi.get('canonical_frac', 0):.4f}",
                f"{vi.get('variant_frac', 0):.4f}",
                f"{vi.get('tti_frac', 0):.5f}",
                r['mh_5p'], r['mh_3p'], r['tsd_len'], r['tsd_seq'],
                ni.get('half_sites_fwd', 0) + ni.get('half_sites_rc', 0),
                ni.get('direct_repeats', 0), ni.get('inverted_repeats', 0),
                r.get('linker', '') or '',
            ]
            f.write('\t'.join(str(x) for x in row) + '\n')

    # Linker FASTA for BLAST
    out_fa = real_dir / "converging_linkers.fa"
    n_linker = 0
    with open(out_fa, 'w') as f:
        for r in analyzed:
            if r['linker_len'] >= 15 and 'CONVERGING' in r['classification']:
                f.write(f">{r['dataset']}__{r['id']}__{r['classification']}__len{r['linker_len']}\n")
                f.write(f"{r['linker']}\n")
                n_linker += 1
    print(f"\nLinker FASTA: {out_fa} ({n_linker} linkers ≥15bp written)")

    # Summary JSON
    summary = {
        'n_total': len(analyzed),
        'n_pan_euk': sum(1 for r in analyzed if r['dataset']=='pan_euk'),
        'n_tara': sum(1 for r in analyzed if r['dataset']=='tara'),
        'class_distribution_by_dataset': {
            ds: dict(d) for ds, d in by_ds_class.items()
        },
        'linker_length_buckets_per_class': {
            c: {
                'n': len(ll),
                'short_5_30': sum(1 for x in ll if 5 <= x <= 30),
                'medium_30_120': sum(1 for x in ll if 30 < x <= 120),
                'long_120plus': sum(1 for x in ll if x > 120),
            }
            for c, ll in linker_lens.items()
        },
        'tsd_presence_per_class': {
            c: {'n': len(v),
                'pct_with_tsd_3plus': 100 * sum(1 for x in v if x >= 3) / len(v) if v else 0,
                'mean_tsd_when_present':
                    sum(x for x in v if x >= 3) / max(1, sum(1 for x in v if x >= 3)),
                'mean_tsd_all': sum(v)/len(v) if v else 0,
            }
            for c, v in tsd_rates.items()
        },
        'variant_frac_per_class': {
            c: {'mean': sum(v)/len(v) if v else 0, 'n': len(v)}
            for c, v in var_frac.items()
        },
        'tti_frac_per_class': {
            c: {'mean': sum(v)/len(v) if v else 0, 'n': len(v)}
            for c, v in tti_frac.items()
        },
        'nr2cf_halfsites_per_class': {
            c: {'mean': sum(v)/len(v) if v else 0,
                'pct_with_any': 100 * sum(1 for x in v if x > 0) / len(v) if v else 0,
                'n': len(v)}
            for c, v in nr2cf_rates.items()
        },
        'microhomology_per_class': {
            c: {'mean_5p': sum(d['mh5'])/len(d['mh5']) if d['mh5'] else 0,
                'mean_3p': sum(d['mh3'])/len(d['mh3']) if d['mh3'] else 0,
                'n': len(d['mh5'])}
            for c, d in mh.items()
        },
    }
    out_json = real_dir / "its_mechanism_v2.json"
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nOutputs:")
    print(f"  {out_tsv}")
    print(f"  {out_json}")
    print(f"  {out_fa}")


if __name__ == '__main__':
    main()
