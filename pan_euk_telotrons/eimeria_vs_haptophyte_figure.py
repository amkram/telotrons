#!/usr/bin/env python3
"""Comparison figure: haptophyte vs Eimeria telotron architecture."""
import csv, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict


def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}


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


def classify(seq, fwd_kmers, rev_kmers):
    fwd_runs = find_runs(seq, fwd_kmers)
    rev_runs = find_runs(seq, rev_kmers)
    if not (fwd_runs and rev_runs): return None, 0, ''
    fwd_max = max(fwd_runs, key=lambda r: r[1]-r[0])
    rev_max = max(rev_runs, key=lambda r: r[1]-r[0])
    if fwd_max[0] < rev_max[0]:
        gap = rev_max[0] - fwd_max[1]
        cls = 'CONVERGING'
        link = seq[fwd_max[1]:rev_max[0]] if gap > 4 else ''
    else:
        gap = fwd_max[0] - rev_max[1]
        cls = 'DIVERGING'
        link = seq[rev_max[1]:fwd_max[0]] if gap > 4 else ''
    if gap > 4: cls += '_WITHLINKER'
    else: cls += '_NOLINKER'
    return cls, gap, link


# Load Tara
TARA_FWD = _rotations('TTAGGG'); TARA_REV = _rotations('CCCTAA')
EIM_FWD = _rotations('TTTAGGG'); EIM_REV = _rotations('CCCTAAA')

tara = {'CONVERGING_WITHLINKER': [], 'DIVERGING_WITHLINKER': [],
        'CONVERGING_NOLINKER': [], 'DIVERGING_NOLINKER': []}
print("Loading Tara...")
with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
    rdr = csv.DictReader(f, delimiter='\t')
    for row in rdr:
        seq = row['intron_seq'].upper()
        cls, gap, link = classify(seq, TARA_FWD, TARA_REV)
        if cls and cls in tara:
            tara[cls].append((gap, link, len(link) if link else 0))

eim = {'CONVERGING_WITHLINKER': [], 'DIVERGING_WITHLINKER': [],
       'CONVERGING_NOLINKER': [], 'DIVERGING_NOLINKER': []}
print("Loading Eimeria...")
seen = set()
for f in Path('ultra_results').glob('*Eimeria*.tsv'):
    with open(f) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            key = (row['contig'], row['start'], row['end'])
            if key in seen: continue
            seen.add(key)
            seq = row['intron_seq'].upper()
            cls, gap, link = classify(seq, EIM_FWD, EIM_REV)
            if cls and cls in eim:
                eim[cls].append((gap, link, len(link) if link else 0))

print(f"Tara: {sum(len(v) for v in tara.values())}")
print(f"Eim:  {sum(len(v) for v in eim.values())}")
print(f"Tara CONV/DIV: {len(tara['CONVERGING_WITHLINKER'])+len(tara['CONVERGING_NOLINKER'])}/{len(tara['DIVERGING_WITHLINKER'])+len(tara['DIVERGING_NOLINKER'])}")
print(f"Eim CONV/DIV:  {len(eim['CONVERGING_WITHLINKER'])+len(eim['CONVERGING_NOLINKER'])}/{len(eim['DIVERGING_WITHLINKER'])+len(eim['DIVERGING_NOLINKER'])}")


def has_telo(seq, fwd, rev):
    return any(k in seq for k in fwd | rev)


def pure_pct(records, fwd, rev):
    pure = [r for r in records if r[1] and not has_telo(r[1], fwd, rev)]
    return 100 * len(pure) / len(records) if records else 0


# === Figure ===
fig, axes = plt.subplots(2, 3, figsize=(17, 10))

# Panel A: Class proportions
ax = axes[0, 0]
classes = ['CONVERGING_NOLINKER', 'CONVERGING_WITHLINKER',
           'DIVERGING_WITHLINKER', 'DIVERGING_NOLINKER']
class_short = ['Conv\nno linker', 'Conv\nwith linker', 'Div\nwith linker', 'Div\nno linker']
tara_counts = [len(tara[c]) for c in classes]
eim_counts = [len(eim[c]) for c in classes]
x = np.arange(len(classes))
w = 0.4
ax.bar(x - w/2, tara_counts, w, label='Haptophytes (Tara)', color='steelblue')
ax.bar(x + w/2, eim_counts, w, label='Eimeria', color='coral')
ax.set_xticks(x); ax.set_xticklabels(class_short, fontsize=8)
ax.set_ylabel('count')
ax.set_yscale('symlog')
ax.set_title('A. Architecture class distribution', loc='left', fontweight='bold')
ax.legend()
for i, (t, e) in enumerate(zip(tara_counts, eim_counts)):
    ax.text(i - w/2, t+0.1, f"{t:,}", ha='center', va='bottom', fontsize=7, rotation=0)
    ax.text(i + w/2, e+0.1, f"{e:,}", ha='center', va='bottom', fontsize=7, rotation=0)

# Panel B: CONV/DIV ratio
ax = axes[0, 1]
ratios = []
labels = []
for ds, data in [('Haptophytes\n(Tara)', tara), ('Eimeria', eim)]:
    conv = len(data['CONVERGING_WITHLINKER']) + len(data['CONVERGING_NOLINKER'])
    div = len(data['DIVERGING_WITHLINKER']) + len(data['DIVERGING_NOLINKER'])
    if div > 0:
        ratio = conv / div
    else:
        ratio = float('inf')
    ratios.append(ratio)
    labels.append(ds)

bars = ax.bar(labels, [r if r < 100 else 100 for r in ratios],
              color=['steelblue', 'coral'])
ax.set_ylabel('CONVERGING : DIVERGING ratio')
ax.set_yscale('log')
ax.set_title('B. CONV vs DIV ratio (haptophyte CONV-dominated, Eimeria DIV-dominated)',
             loc='left', fontweight='bold', fontsize=9)
ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Equal')
for b, r in zip(bars, ratios):
    if r > 1:
        s = f"{r:.1f}:1"
    else:
        s = f"1:{1/r:.1f}"
    ax.text(b.get_x() + b.get_width()/2, b.get_height(), s, ha='center', va='bottom', fontsize=10, fontweight='bold')

# Panel C: Linker length distribution comparison
ax = axes[0, 2]
tara_lens = [g for g, l, ll in tara['CONVERGING_WITHLINKER'] + tara['DIVERGING_WITHLINKER'] if g >= 5]
eim_lens = [g for g, l, ll in eim['CONVERGING_WITHLINKER'] + eim['DIVERGING_WITHLINKER'] if g >= 5]
bins = np.linspace(5, 300, 40)
ax.hist(tara_lens, bins=bins, alpha=0.5, label=f'Haptophytes (n={len(tara_lens):,})', color='steelblue', density=True)
ax.hist(eim_lens, bins=bins, alpha=0.5, label=f'Eimeria (n={len(eim_lens):,})', color='coral', density=True)
ax.axvspan(5, 30, alpha=0.1, color='red')
ax.axvspan(30, 120, alpha=0.1, color='gold')
ax.axvspan(120, 300, alpha=0.1, color='green')
ax.set_xlabel('Linker length (bp)')
ax.set_ylabel('density')
ax.set_title('C. Linker length distribution', loc='left', fontweight='bold')
ax.legend()

# Panel D: Pure non-repeat % per dataset
ax = axes[1, 0]
groups = ['CONV-w/linker', 'DIV-w/linker']
tara_pure = [
    pure_pct(tara['CONVERGING_WITHLINKER'], TARA_FWD, TARA_REV),
    pure_pct(tara['DIVERGING_WITHLINKER'], TARA_FWD, TARA_REV),
]
eim_pure = [
    pure_pct(eim['CONVERGING_WITHLINKER'], EIM_FWD, EIM_REV),
    pure_pct(eim['DIVERGING_WITHLINKER'], EIM_FWD, EIM_REV),
]
x = np.arange(len(groups))
w = 0.35
ax.bar(x - w/2, tara_pure, w, label='Haptophytes (Tara)', color='steelblue')
ax.bar(x + w/2, eim_pure, w, label='Eimeria', color='coral')
ax.set_xticks(x); ax.set_xticklabels(groups)
ax.set_ylabel('% truly non-repeat linkers')
ax.set_title('D. Pure non-repeat linker fraction (Eimeria higher)', loc='left', fontweight='bold')
for i, (t, e) in enumerate(zip(tara_pure, eim_pure)):
    ax.text(i - w/2, t+0.5, f"{t:.1f}%", ha='center', va='bottom', fontsize=9)
    ax.text(i + w/2, e+0.5, f"{e:.1f}%", ha='center', va='bottom', fontsize=9)
ax.legend()

# Panel E: length bucket comparison (CONVERGING_WITHLINKER + DIVERGING_WITHLINKER)
ax = axes[1, 1]
buckets = ['short\n(5-30)', 'medium\n(30-120)', 'long\n(120+)']
def bucketize(lens):
    return [
        sum(1 for x in lens if 5 <= x < 30),
        sum(1 for x in lens if 30 <= x < 120),
        sum(1 for x in lens if x >= 120),
    ]
tara_b = bucketize(tara_lens)
eim_b = bucketize(eim_lens)
total_t = sum(tara_b); total_e = sum(eim_b)
tara_pct = [100 * x / total_t if total_t else 0 for x in tara_b]
eim_pct = [100 * x / total_e if total_e else 0 for x in eim_b]
x = np.arange(len(buckets))
w = 0.35
ax.bar(x - w/2, tara_pct, w, label='Haptophytes', color='steelblue')
ax.bar(x + w/2, eim_pct, w, label='Eimeria', color='coral')
ax.set_xticks(x); ax.set_xticklabels(buckets, fontsize=8)
ax.set_ylabel('% of linkers')
ax.set_title('E. Linker length buckets (Eimeria longer)', loc='left', fontweight='bold')
for i, (t, e) in enumerate(zip(tara_pct, eim_pct)):
    ax.text(i - w/2, t+1, f"{t:.0f}%", ha='center', va='bottom', fontsize=9)
    ax.text(i + w/2, e+1, f"{e:.0f}%", ha='center', va='bottom', fontsize=9)
ax.legend()

# Panel F: Schematic illustration
ax = axes[1, 2]
ax.axis('off')
schem_text = """
ARCHITECTURE COMPARISON

HAPTOPHYTE (CONVERGING):
  5'─exon─[GT]─(TTAGGG)n→...←(CCCTAA)n─[AG]─exon─3'
  G-rich strands point INWARD (toward linker)
  ↳ Classical 2-ended dnTA at DSB
    (Nergadze 2004 model)

EIMERIA (DIVERGING):
  5'─exon─[GT]─(CCCTAAA)n←...→(TTTAGGG)n─[AG]─exon─3'
  G-rich strands point OUTWARD (away from linker)
  ↳ NOT classical dnTA. Likely:
    - Subtelomeric fragment capture, or
    - Inversion after dnTA, or
    - Different DSB repair geometry

KEY DIFFERENCE:
  Haptophyte linker = degenerate telomeric DNA
  (8.6% truly non-repeat, median 44 bp)

  Eimeria linker = often genuine non-repeat
  (25.9% truly non-repeat, median 123 bp)
  → fragment capture from elsewhere?
"""
ax.text(0, 1, schem_text, fontsize=8.5, family='monospace', va='top', ha='left',
        transform=ax.transAxes)
ax.set_title('F. Mechanism interpretation', loc='left', fontweight='bold')

plt.suptitle('Haptophyte vs Eimeria telotron architecture: '
             'two distinct mechanisms?', fontsize=14, fontweight='bold')
plt.tight_layout()

out = Path('real_telotrons/eimeria_vs_haptophyte.png')
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")
