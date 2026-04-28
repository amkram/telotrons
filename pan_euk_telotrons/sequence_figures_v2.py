#!/usr/bin/env python3
"""
Sequence-level figures v2 — improved layouts, larger fonts, all panels working.

Outputs (under real_telotrons/):
  seqfig_classes.png      — Per-class examples with full sequence detail
  seqfig_shared.png       — 5 shared-linker tandem pairs (working)
  seqfig_compare.png      — Haptophyte vs Eimeria side-by-side
  seqfig_linkers.png      — Pure non-repeat vs degenerate linker examples
  seqfig_tandem.png       — Tandem cluster examples (multiple telotrons close together)
"""

import csv
import json
import re
import random
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.font_manager import FontProperties

from sequence_viewer import classify_bases, render_row, COLOR_BG

random.seed(42)


def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

# Repeat sets per family
FAMILIES = {
    'TTAGGG':  {'fwd': _rotations('TTAGGG'),  'rev': _rotations('CCCTAA')},
    'TTTAGGG': {'fwd': _rotations('TTTAGGG'), 'rev': _rotations('CCCTAAA')},
}


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


def classify_orientation(seq, family):
    sets = FAMILIES.get(family)
    if not sets: return 'NONE', None
    fwd = find_runs(seq, sets['fwd'])
    rev = find_runs(seq, sets['rev'])
    if not fwd and not rev: return 'NONE', None
    n = len(seq)
    fbp = sum(e-s for s,e in fwd); rbp = sum(e-s for s,e in rev)
    has_f = fbp >= 10 and fbp/n >= 0.05
    has_r = rbp >= 10 and rbp/n >= 0.05
    if has_f and not has_r: return 'SINGLE_G', None
    if has_r and not has_f: return 'SINGLE_C', None
    if not has_f and not has_r: return 'NONE', None
    fmax = max(fwd, key=lambda r: r[1]-r[0])
    rmax = max(rev, key=lambda r: r[1]-r[0])
    if fmax[0] < rmax[0]:
        gap = rmax[0] - fmax[1]
        return ('CONVERGING_NOLINKER' if gap <= 4 else 'CONVERGING_WITHLINKER',
                (fmax[1], rmax[0]) if gap > 4 else None)
    else:
        gap = fmax[0] - rmax[1]
        return ('DIVERGING_NOLINKER' if gap <= 4 else 'DIVERGING_WITHLINKER',
                (rmax[1], fmax[0]) if gap > 4 else None)


def has_telo_kmer(seq, family):
    kmers = FAMILIES[family]['fwd'] | FAMILIES[family]['rev']
    return any(k in seq for k in kmers)


def find_tsd(left, right, min_len=3, max_len=20):
    if not left or not right: return 0
    L = min(max_len, len(left), len(right))
    for tl in range(L, min_len-1, -1):
        if left[-tl:] == right[:tl]: return tl
    return 0


# ------------------------------------------------------------------
# Loaders
# ------------------------------------------------------------------

def load_tara():
    out = []
    p = Path('real_telotrons/_tara_oceans_telotrons.tsv')
    if not p.exists(): return out
    with open(p) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            if r.get('intron_seq'):
                r['_dataset'] = 'tara'
                r['_family'] = 'TTAGGG'
                out.append(r)
    return out


def load_eimeria_with_flanks():
    """Load Eimeria telotrons and extract flanks from genomes on the fly."""
    # Cache contig sequences per acc
    contig_cache = {}

    eim = []
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                if row.get('intron_seq'):
                    row['acc'] = acc
                    row['_dataset'] = 'eimeria'
                    row['_family'] = 'TTTAGGG'
                    eim.append(row)
    print(f"Loaded {len(eim)} Eimeria telotrons; extracting flanks from genomes...")

    for acc in sorted(set(r['acc'] for r in eim)):
        contig_cache[acc] = {}
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen:
            print(f"  Genome not found for {acc}")
            continue
        cur = None; parts = []
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contig_cache[acc][cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]
                    parts = []
                else:
                    parts.append(line.strip().upper())
            if cur: contig_cache[acc][cur] = ''.join(parts)
        print(f"  {acc}: {len(contig_cache[acc])} contigs loaded")

    n_with_flanks = 0
    for r in eim:
        contig_seq = contig_cache.get(r['acc'], {}).get(r['contig'])
        if not contig_seq: continue
        try:
            s = int(r['start']); e = int(r['end'])
        except: continue
        r['left_flank_100bp'] = contig_seq[max(0, s-100):s]
        r['right_flank_100bp'] = contig_seq[e:e+100]
        n_with_flanks += 1
    print(f"  {n_with_flanks} Eimeria telotrons have flanks")

    return eim


# ------------------------------------------------------------------
# Common rendering helpers
# ------------------------------------------------------------------

def draw_intron_with_flanks(ax, y, intron, left_flank, right_flank,
                             family, char_w=1.0, flank_bp=40, fontsize=3.5,
                             tsd_len=0, highlight=None,
                             x_offset=0):
    """Draw a row: [left_flank] [splice5'] [intron] [splice3'] [right_flank]
    Returns the right edge x-coordinate.
    highlight: dict {position: 'shared'|'tsd'} for per-base override
    """
    lf = left_flank[-flank_bp:] if left_flank else ''
    rf = right_flank[:flank_bp] if right_flank else ''
    x = x_offset

    # Left flank — with TSD highlight
    if lf:
        lf_classes = ['flank'] * len(lf)
        tsd_overrides = []
        if tsd_len > 0:
            tsd_overrides = [(len(lf) - tsd_len, len(lf), 'tsd')]
        render_row(ax, y, lf, lf_classes, x_start=x, char_width=char_w,
                   highlight_ranges=tsd_overrides, fontsize=fontsize)
        x += len(lf) * char_w

    # 5' splice site
    ax.add_patch(Rectangle((x, y), char_w * 0.6, 1.0, facecolor='black', edgecolor='none'))
    x += char_w * 0.6

    # Intron
    cls_per_base = classify_bases(intron, family)
    intron_highlights = []
    if highlight:
        # group by type
        by_type = defaultdict(list)
        for pos, t in highlight.items():
            by_type[t].append(pos)
        for t, positions in by_type.items():
            if not positions: continue
            positions.sort()
            # Make ranges
            start = positions[0]
            prev = positions[0]
            for p in positions[1:] + [None]:
                if p is None or p != prev + 1:
                    intron_highlights.append((start, prev + 1, t))
                    if p is not None:
                        start = p
                if p is not None:
                    prev = p
    render_row(ax, y, intron, cls_per_base, x_start=x, char_width=char_w,
               highlight_ranges=intron_highlights, fontsize=fontsize)
    x += len(intron) * char_w

    # 3' splice site
    ax.add_patch(Rectangle((x, y), char_w * 0.6, 1.0, facecolor='black', edgecolor='none'))
    x += char_w * 0.6

    # Right flank
    if rf:
        rf_classes = ['flank'] * len(rf)
        tsd_overrides2 = []
        if tsd_len > 0:
            tsd_overrides2 = [(0, tsd_len, 'tsd')]
        render_row(ax, y, rf, rf_classes, x_start=x, char_width=char_w,
                   highlight_ranges=tsd_overrides2, fontsize=fontsize)
        x += len(rf) * char_w

    return x


def add_legend(ax, x, y, fontsize=8, vertical=False):
    items = [
        ('fwd_canonical', 'TTAGGG/TTTAGGG (G-strand)'),
        ('rev_canonical', 'CCCTAA/CCCTAAA (C-strand)'),
        ('fwd_variant', '1-mm fwd variant'),
        ('rev_variant', '1-mm rev variant'),
        ('flank', 'Flanking exon'),
        ('non_telo', 'Non-repeat'),
        ('tsd', 'TSD'),
        ('shared', 'Shared with paralog'),
    ]
    for i, (cls, lbl) in enumerate(items):
        if vertical:
            ax.add_patch(Rectangle((x, y - i * 1.5), 3, 1, facecolor=COLOR_BG[cls], edgecolor='black', linewidth=0.3))
            ax.text(x + 4, y - i * 1.5 + 0.5, lbl, fontsize=fontsize, va='center')
        else:
            ax.add_patch(Rectangle((x + i * 30, y), 3, 1, facecolor=COLOR_BG[cls], edgecolor='black', linewidth=0.3))
            ax.text(x + i * 30 + 4, y + 0.5, lbl, fontsize=fontsize, va='center')


# ------------------------------------------------------------------
# FIGURE 1: All classes (improved)
# ------------------------------------------------------------------

def figure_classes(out='real_telotrons/seqfig_classes.png', n_per=4):
    tara = load_tara()
    eim = load_eimeria_with_flanks()
    print(f"\n=== Figure 1: Class examples ===")

    # Bin by class
    by_class = defaultdict(list)
    for r in tara + eim:
        seq = r.get('intron_seq', '').upper()
        if not seq or len(seq) < 40 or len(seq) > 250: continue
        family = r['_family']
        cls, link = classify_orientation(seq, family)
        if cls == 'NONE': continue
        left = r.get('left_flank_100bp', '').upper()
        right = r.get('right_flank_100bp', '').upper()
        if not left or not right: continue
        by_class[cls].append({'rec': r, 'seq': seq, 'left': left, 'right': right,
                              'family': family, 'link': link})

    classes_show = ['SINGLE_G', 'SINGLE_C',
                     'CONVERGING_NOLINKER', 'CONVERGING_WITHLINKER',
                     'DIVERGING_WITHLINKER']

    selected = {}
    for c in classes_show:
        cands = by_class.get(c, [])
        # Sort by intron length ascending for visibility
        cands.sort(key=lambda x: len(x['seq']))
        # Filter linker classes for visible linkers
        if 'WITHLINKER' in c:
            cands = [x for x in cands if x['link'] and 15 <= (x['link'][1]-x['link'][0]) <= 60]
        # Take a few from start (short) and middle (medium)
        if len(cands) > n_per:
            step = len(cands) // n_per
            sel = [cands[i*step] for i in range(n_per)]
        else:
            sel = cands
        selected[c] = sel
        print(f"  {c}: {len(by_class.get(c, []))} avail, {len(sel)} selected")

    # Layout
    char_w = 1.0
    flank_bp = 40
    max_intron = max((len(x['seq']) for cls in selected.values() for x in cls), default=200)
    row_w = flank_bp + 1 + max_intron + 1 + flank_bp + 5
    n_rows = sum(len(v) + 2 for v in selected.values())  # +2 for title and gap

    fig_w = max(20, row_w * 0.06)
    fig_h = max(12, n_rows * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-30, row_w + 5)
    ax.set_ylim(n_rows + 3, -2)
    ax.axis('off')

    cur_y = 0
    class_color = {'SINGLE_G': '#1955a8', 'SINGLE_C': '#a82020',
                   'CONVERGING_NOLINKER': '#1f6e1f', 'CONVERGING_WITHLINKER': '#1f6e1f',
                   'DIVERGING_WITHLINKER': '#5818a8'}
    class_desc = {
        'SINGLE_G': "single-strand TTAGGG (G-rich on gene strand)",
        'SINGLE_C': "single-strand CCCTAA (C-rich on gene strand)",
        'CONVERGING_NOLINKER': "5'-(TTAGGG)→…←(CCCTAA)-3'  arrays meet directly",
        'CONVERGING_WITHLINKER': "5'-(TTAGGG)→ LINKER ←(CCCTAA)-3'  with intervening sequence",
        'DIVERGING_WITHLINKER': "5'-(CCCTAA)← LINKER →(TTAGGG)-3'  arrays point outward",
    }

    for c in classes_show:
        # Title row
        ax.text(-30, cur_y + 0.5, c, ha='left', va='center',
                fontsize=12, fontweight='bold', color=class_color[c])
        ax.text(-30 + 65, cur_y + 0.5, class_desc[c],
                ha='left', va='center', fontsize=9, fontstyle='italic', color='gray')
        n_total = len(by_class.get(c, []))
        ax.text(row_w, cur_y + 0.5, f"(n={n_total:,} total)",
                ha='right', va='center', fontsize=9, color='gray')
        cur_y += 1

        for ex in selected[c]:
            seq = ex['seq']; left = ex['left']; right = ex['right']
            family = ex['family']
            tsd = find_tsd(left[-flank_bp:], right[:flank_bp])
            draw_intron_with_flanks(ax, cur_y, seq, left, right, family,
                                     char_w=char_w, flank_bp=flank_bp,
                                     fontsize=3.0, tsd_len=tsd)
            r = ex['rec']
            sp = r.get('species', r.get('mag', '?'))
            sp = sp.replace('TARA_', '').split('_MAG_')[0] if 'TARA' in sp else sp
            sp = sp.replace('Eimeria_', 'E. ')
            lbl = f"{sp[:15]:<15} {len(seq)}bp"
            ax.text(-30, cur_y + 0.5, lbl, ha='left', va='center', fontsize=7)
            cur_y += 1
        cur_y += 1  # gap

    # Header — column labels
    ax.text(-30 + flank_bp/2, -1, "Upstream exon (-40bp)", ha='center', va='center',
            fontsize=9, color='gray')
    ax.text(-30 + flank_bp + 1 + max_intron / 2, -1, "Telotron intron (with array(s))",
            ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(-30 + flank_bp + 1 + max_intron + 1 + flank_bp / 2, -1, "Downstream exon (+40bp)",
            ha='center', va='center', fontsize=9, color='gray')

    # Legend at bottom
    add_legend(ax, 0, n_rows + 1, fontsize=8, vertical=False)

    plt.suptitle('Sequence-level examples of telotron orientation classes',
                  fontsize=15, fontweight='bold')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")


# ------------------------------------------------------------------
# FIGURE 2: Shared linker pairs (working version)
# ------------------------------------------------------------------

def figure_shared(out='real_telotrons/seqfig_shared.png'):
    print(f"\n=== Figure 2: Shared linker pairs ===")
    pairs_path = Path('real_telotrons/shared_linker_pairs.json')
    if not pairs_path.exists():
        print("  shared_linker_pairs.json not found; skip")
        return
    with open(pairs_path) as f:
        pairs = json.load(f)
    print(f"  Loaded {len(pairs)} pairs")

    # Need flanks for these — load from genome
    contig_cache = defaultdict(dict)
    accs_needed = set(p['acc'] for p in pairs)
    for acc in accs_needed:
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen: continue
        cur = None; parts = []
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contig_cache[acc][cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]
                    parts = []
                else:
                    parts.append(line.strip().upper())
            if cur: contig_cache[acc][cur] = ''.join(parts)

    # Render: 1 figure with 5 sub-panels (one per pair)
    fig, axes = plt.subplots(len(pairs), 1, figsize=(22, len(pairs) * 4))
    if len(pairs) == 1: axes = [axes]

    for idx, p in enumerate(pairs):
        ax = axes[idx]
        ax.axis('off')
        acc = p['acc']
        t1 = p['src_telotron']
        t2 = p['paralog_telotron']

        # Get the linker subsequence from t1
        linker_seq = t1['intron_seq'][t1['linker_offset']:t1['linker_end_offset']]

        # Find best match in t2 (allow 1-2 mismatches)
        seq2 = t2['intron_seq']
        best_score = -1; best_pos = -1
        L = len(linker_seq)
        for i in range(len(seq2) - L + 1):
            window = seq2[i:i+L]
            matches = sum(a == b for a, b in zip(window, linker_seq))
            if matches > best_score:
                best_score = matches
                best_pos = i

        # Get flanks
        contig_seqs = contig_cache.get(acc, {})
        cs1 = contig_seqs.get(t1['contig'], '')
        left1 = cs1[max(0, t1['start']-30):t1['start']] if cs1 else ''
        right1 = cs1[t1['end']:t1['end']+30] if cs1 else ''
        cs2 = contig_seqs.get(t2['contig'], '')
        left2 = cs2[max(0, t2['start']-30):t2['start']] if cs2 else ''
        right2 = cs2[t2['end']:t2['end']+30] if cs2 else ''

        seq1 = t1['intron_seq']

        # Render
        char_w = 1.0
        flank_bp = 30
        max_intron = max(len(seq1), len(seq2)) + 100

        ax.set_xlim(-25, max_intron + 5)
        ax.set_ylim(8, -3)

        same_contig = t1['contig'] == t2['contig']
        gap_info = ''
        if same_contig:
            gap = abs(t2['start'] - t1['end']) if t2['start'] > t1['end'] else abs(t1['start'] - t2['end'])
            gap_info = f"  [SAME CONTIG, gap={gap:,}bp]"

        ax.text(0, -2,
                f"Pair {idx+1}: {acc}{gap_info}",
                fontsize=10, fontweight='bold')
        ax.text(0, -1,
                f"Shared linker: {linker_seq}  (length {L}bp, BLAST identity = {best_score}/{L} = {100*best_score/L:.0f}%)",
                fontsize=8, family='monospace', color='darkgreen')

        # Telotron 1 with linker highlighted
        highlight1 = {pos: 'shared' for pos in range(t1['linker_offset'], t1['linker_end_offset'])}
        end_x1 = draw_intron_with_flanks(ax, 1, seq1, left1, right1, 'TTTAGGG',
                                          char_w=char_w, flank_bp=flank_bp,
                                          fontsize=3.5,
                                          highlight=highlight1)
        ax.text(-25, 1.5,
                f"{t1['contig']}\n:{t1['start']}-{t1['end']}\n({t1['length']}bp)",
                ha='left', va='center', fontsize=6.5)

        # Telotron 2 with paralog highlighted
        if best_pos >= 0:
            highlight2 = {pos: 'shared' for pos in range(best_pos, best_pos + L)}
        else:
            highlight2 = {}
        draw_intron_with_flanks(ax, 4, seq2, left2, right2, 'TTTAGGG',
                                 char_w=char_w, flank_bp=flank_bp,
                                 fontsize=3.5,
                                 highlight=highlight2)
        ax.text(-25, 4.5,
                f"{t2['contig']}\n:{t2['start']}-{t2['end']}\n({t2['length']}bp)",
                ha='left', va='center', fontsize=6.5)

        # Connect with arrow
        # Compute approximate x position of linker in t1
        t1_x = flank_bp + 0.6 + (t1['linker_offset'] + L/2)
        t2_x = flank_bp + 0.6 + (best_pos + L/2) if best_pos >= 0 else flank_bp + 0.6
        arr = FancyArrowPatch((t1_x, 2), (t2_x, 4),
                               connectionstyle="arc3,rad=0.2",
                               arrowstyle='<->', color='green',
                               linewidth=1.5, alpha=0.7)
        ax.add_patch(arr)

    add_legend(axes[0], 50, -3, fontsize=7)

    plt.suptitle('Shared-linker tandem cluster pairs in Eimeria — direct evidence of common origin',
                  fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")


# ------------------------------------------------------------------
# FIGURE 3: Hapt vs Eimeria (improved)
# ------------------------------------------------------------------

def figure_compare(out='real_telotrons/seqfig_compare.png', n_each=8):
    print(f"\n=== Figure 3: Hapt vs Eimeria ===")
    tara = load_tara()
    eim = load_eimeria_with_flanks()

    hapt = []
    for r in tara:
        seq = r.get('intron_seq', '').upper()
        if not seq or len(seq) < 80 or len(seq) > 200: continue
        if not r.get('left_flank_100bp') or not r.get('right_flank_100bp'): continue
        cls, link = classify_orientation(seq, 'TTAGGG')
        if cls == 'CONVERGING_WITHLINKER' and link and 15 <= link[1]-link[0] <= 60:
            hapt.append({'rec': r, 'seq': seq, 'family': 'TTAGGG'})
    random.shuffle(hapt)
    hapt = hapt[:n_each]

    eim_div = []
    for r in eim:
        seq = r.get('intron_seq', '').upper()
        if not seq or len(seq) < 80 or len(seq) > 200: continue
        if not r.get('left_flank_100bp') or not r.get('right_flank_100bp'): continue
        cls, link = classify_orientation(seq, 'TTTAGGG')
        if cls == 'DIVERGING_WITHLINKER' and link and 15 <= link[1]-link[0] <= 60:
            eim_div.append({'rec': r, 'seq': seq, 'family': 'TTTAGGG'})
    random.shuffle(eim_div)
    eim_div = eim_div[:n_each]
    print(f"  Hapt: {len(hapt)}, Eim: {len(eim_div)}")

    char_w = 1.0
    flank_bp = 35
    col_max = max((len(x['seq']) for x in hapt + eim_div), default=200)
    col_w = flank_bp + 1 + col_max + 1 + flank_bp
    gap = 40
    n_rows = max(len(hapt), len(eim_div))

    fig, ax = plt.subplots(figsize=(28, max(9, n_rows * 0.7) + 2))
    ax.set_xlim(-25, 2 * col_w + gap + 5)
    ax.set_ylim(n_rows + 3, -3)
    ax.axis('off')

    # Headers
    ax.text(col_w / 2, -2.5, "HAPTOPHYTES (Tara) — CONVERGING",
            ha='center', va='center', fontsize=14, fontweight='bold', color='steelblue')
    ax.text(col_w / 2, -1.5, "5'-(TTAGGG)n →...← (CCCTAA)m-3'",
            ha='center', va='center', fontsize=10, family='monospace', color='steelblue')

    eim_x = col_w + gap
    ax.text(eim_x + col_w / 2, -2.5, "EIMERIA — DIVERGING",
            ha='center', va='center', fontsize=14, fontweight='bold', color='coral')
    ax.text(eim_x + col_w / 2, -1.5, "5'-(CCCTAAA)n ←...→ (TTTAGGG)m-3'",
            ha='center', va='center', fontsize=10, family='monospace', color='coral')

    # Render rows
    for i, ex in enumerate(hapt):
        r = ex['rec']
        left = r['left_flank_100bp']; right = r['right_flank_100bp']
        tsd = find_tsd(left[-flank_bp:], right[:flank_bp])
        draw_intron_with_flanks(ax, i, ex['seq'], left, right, 'TTAGGG',
                                 char_w=char_w, flank_bp=flank_bp, fontsize=3.0,
                                 tsd_len=tsd)
        mag = r.get('mag', 'TARA').replace('TARA_', '').split('_MAG_')[0]
        ax.text(-25, i + 0.5, f"{mag} {len(ex['seq'])}bp",
                ha='left', va='center', fontsize=7)

    for i, ex in enumerate(eim_div):
        r = ex['rec']
        left = r['left_flank_100bp']; right = r['right_flank_100bp']
        tsd = find_tsd(left[-flank_bp:], right[:flank_bp])
        draw_intron_with_flanks(ax, i, ex['seq'], left, right, 'TTTAGGG',
                                 char_w=char_w, flank_bp=flank_bp, fontsize=3.0,
                                 tsd_len=tsd, x_offset=eim_x)
        sp = r.get('species', '').replace('Eimeria_', 'E.')
        ax.text(eim_x - 25, i + 0.5, f"{sp[:14]} {len(ex['seq'])}bp",
                ha='left', va='center', fontsize=7)

    add_legend(ax, 0, n_rows + 1, fontsize=8)

    plt.suptitle('Sequence-level: Haptophyte CONVERGING vs Eimeria DIVERGING — distinct architectures',
                  fontsize=14, fontweight='bold')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")


# ------------------------------------------------------------------
# FIGURE 4: Linker sequences (pure non-repeat vs degenerate)
# ------------------------------------------------------------------

def figure_linkers(out='real_telotrons/seqfig_linkers.png', n_per=8):
    print(f"\n=== Figure 4: Linker examples ===")
    tara = load_tara()
    eim = load_eimeria_with_flanks()

    pure_eim, mosaic_eim, pure_tara, mosaic_tara = [], [], [], []
    for r in eim:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, link = classify_orientation(seq, 'TTTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if not (20 <= len(linker) <= 80): continue
        if has_telo_kmer(linker, 'TTTAGGG'):
            mosaic_eim.append({'rec': r, 'seq': seq, 'family': 'TTTAGGG', 'link': link})
        else:
            pure_eim.append({'rec': r, 'seq': seq, 'family': 'TTTAGGG', 'link': link})

    for r in tara:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, link = classify_orientation(seq, 'TTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if not (20 <= len(linker) <= 80): continue
        if has_telo_kmer(linker, 'TTAGGG'):
            mosaic_tara.append({'rec': r, 'seq': seq, 'family': 'TTAGGG', 'link': link})
        else:
            pure_tara.append({'rec': r, 'seq': seq, 'family': 'TTAGGG', 'link': link})

    random.shuffle(pure_eim); random.shuffle(mosaic_eim)
    random.shuffle(pure_tara); random.shuffle(mosaic_tara)
    pure_eim, mosaic_eim = pure_eim[:n_per], mosaic_eim[:n_per]
    pure_tara, mosaic_tara = pure_tara[:n_per], mosaic_tara[:n_per]

    sections = [
        ('Tara haptophytes — DEGENERATE telomeric "linker" (mosaic, 91% of cases)',
            mosaic_tara, 'TTAGGG'),
        ('Tara haptophytes — TRULY non-repeat linker (8.6% of cases)',
            pure_tara, 'TTAGGG'),
        ('Eimeria — DEGENERATE telomeric "linker" (mosaic)',
            mosaic_eim, 'TTTAGGG'),
        ('Eimeria — TRULY non-repeat linker (25.9% of CONVERGING_WITHLINKER)',
            pure_eim, 'TTTAGGG'),
    ]
    total = sum(len(s[1]) + 1 for s in sections)

    char_w = 1.0
    max_intron = max((len(x['seq']) for s in sections for x in s[1]), default=200)

    fig, ax = plt.subplots(figsize=(20, max(12, total * 0.45)))
    ax.set_xlim(-30, max_intron + 50)
    ax.set_ylim(total + 2, -1)
    ax.axis('off')

    cur_y = 0
    for title, examples, family in sections:
        ax.text(-30, cur_y + 0.5, title, ha='left', va='center',
                fontsize=11, fontweight='bold')
        cur_y += 1
        for ex in examples:
            seq = ex['seq']
            ls, le = ex['link']
            linker = seq[ls:le]
            cls = classify_bases(seq, family)
            x_start = 0
            render_row(ax, cur_y, seq, cls, x_start=x_start, char_width=char_w, fontsize=3.5)
            # Mark linker boundaries
            ax.plot([ls, ls], [cur_y, cur_y + 1], color='green', linewidth=1.5, linestyle='--')
            ax.plot([le, le], [cur_y, cur_y + 1], color='green', linewidth=1.5, linestyle='--')
            ax.text(-30, cur_y + 0.5,
                    f"{len(seq)}bp / link={len(linker)}bp",
                    ha='left', va='center', fontsize=7)
            # Show the linker sequence as text on the right
            ax.text(max_intron + 5, cur_y + 0.5, f"linker: {linker[:50]}",
                    ha='left', va='center', fontsize=6.5, family='monospace', color='darkgreen')
            cur_y += 1
        cur_y += 0.5

    plt.suptitle('Linker sequences: degenerate telomeric vs truly non-repeat',
                  fontsize=14, fontweight='bold')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")


# ------------------------------------------------------------------
# FIGURE 5: Tandem clusters — multiple telotrons close together
# ------------------------------------------------------------------

def figure_tandem_cluster(out='real_telotrons/seqfig_tandem.png'):
    """Show specific contigs with multiple telotrons clustered close together
    (genomic context view, with the actual sequence between them)."""
    print(f"\n=== Figure 5: Tandem clusters ===")
    eim = load_eimeria_with_flanks()

    # Group by (acc, contig)
    by_contig = defaultdict(list)
    for r in eim:
        seq = r.get('intron_seq', '')
        if not seq: continue
        try:
            by_contig[(r['acc'], r['contig'])].append({
                'start': int(r['start']), 'end': int(r['end']),
                'rec': r, 'seq': seq.upper(),
            })
        except: continue

    # Find contigs with 3-5 closely-spaced telotrons
    interesting = []
    for key, telos in by_contig.items():
        if len(telos) < 3 or len(telos) > 5: continue
        telos.sort(key=lambda x: x['start'])
        # Compute gaps
        gaps = [telos[i+1]['start'] - telos[i]['end'] for i in range(len(telos)-1)]
        if all(0 <= g <= 1000 for g in gaps):
            interesting.append((key, telos, sum(gaps)))

    # Sort by total compactness
    interesting.sort(key=lambda x: x[2])
    print(f"  Found {len(interesting)} contigs with 3-5 tandemly-clustered telotrons")
    if not interesting:
        print("  Skipping figure")
        return
    selected = interesting[:5]

    # Render each as a row: full contig stretch from first telotron - 30bp to last + 30bp
    fig, axes = plt.subplots(len(selected), 1, figsize=(22, len(selected) * 4))
    if len(selected) == 1: axes = [axes]

    for idx, ((acc, contig), telos, total_gap) in enumerate(selected):
        ax = axes[idx]
        ax.axis('off')

        # Get genome contig
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        contig_seq = ''
        if gen:
            cur = None; parts = []; found = False
            with open(gen) as f:
                for line in f:
                    if line.startswith('>'):
                        if found: break
                        cur_c = line[1:].strip().split()[0]
                        if cur_c == contig:
                            found = True
                    elif found:
                        parts.append(line.strip().upper())
            contig_seq = ''.join(parts)

        if not contig_seq:
            ax.text(0, 0, f"Could not load {acc} {contig}", fontsize=10)
            continue

        # Define visualization range
        vis_start = max(0, telos[0]['start'] - 30)
        vis_end = min(len(contig_seq), telos[-1]['end'] + 30)
        vis_len = vis_end - vis_start

        ax.set_xlim(-30, vis_len + 5)
        ax.set_ylim(3, -3)

        # Title
        ax.text(0, -2.5,
                f"Cluster {idx+1}: {acc} {contig} — {len(telos)} tandem telotrons in {vis_len}bp",
                fontsize=11, fontweight='bold')

        # Render the full stretch
        full_seq = contig_seq[vis_start:vis_end]

        # Build base classes: default 'flank' (intergenic/exon), then mark each telotron interior with classify_bases
        full_classes = ['flank'] * len(full_seq)
        for t in telos:
            t_start_in_vis = t['start'] - vis_start
            t_end_in_vis = t['end'] - vis_start
            intron_seq = full_seq[t_start_in_vis:t_end_in_vis]
            intron_classes = classify_bases(intron_seq, 'TTTAGGG')
            for i, c in enumerate(intron_classes):
                if t_start_in_vis + i < len(full_classes):
                    full_classes[t_start_in_vis + i] = c

        # Render
        render_row(ax, 0, full_seq, full_classes, x_start=0, char_width=1.0, fontsize=3.5)

        # Mark splice sites for each telotron with black bars
        for t in telos:
            ts = t['start'] - vis_start
            te = t['end'] - vis_start
            ax.add_patch(Rectangle((ts - 0.3, 0), 0.6, 1.0, facecolor='black', edgecolor='none'))
            ax.add_patch(Rectangle((te - 0.3, 0), 0.6, 1.0, facecolor='black', edgecolor='none'))
            # Telotron label
            ax.text((ts + te) / 2, 1.5, f"T{telos.index(t)+1}\n{t['end']-t['start']}bp",
                    ha='center', va='center', fontsize=6.5, fontweight='bold')
            # Position
            ax.text(ts, 2.5, f"{t['start']}", ha='center', va='center', fontsize=5.5, color='gray')

    add_legend(axes[0], 100, -2, fontsize=7)

    plt.suptitle('Eimeria tandem-cluster examples: multiple telotrons within 1kb on same contig',
                  fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")


if __name__ == '__main__':
    figure_classes()
    figure_shared()
    figure_compare()
    figure_linkers()
    figure_tandem_cluster()
    print("\nAll figures generated.")
