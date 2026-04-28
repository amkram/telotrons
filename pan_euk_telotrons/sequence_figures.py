#!/usr/bin/env python3
"""
Generate sequence-level figures for ITS mechanism analyses.

Figures:
  1. all_classes_examples.png:
     5 rows × 5 examples = 25 representative sequences across all orientation
     classes (SINGLE_G, SINGLE_C, CONVERGING_NOLINKER, CONVERGING_WITHLINKER,
     DIVERGING_WITHLINKER) with their flanks.

  2. shared_linker_pairs.png:
     The 5 specific shared-linker tandem-cluster cases. Shows both telotrons
     of each pair with the shared sequence highlighted in green.

  3. haptophyte_vs_eimeria.png:
     Side-by-side: 8 haptophyte CONVERGING examples vs 8 Eimeria DIVERGING
     examples, each with flanks, to highlight the architecture difference.

  4. linker_origin_examples.png:
     Examples of pure non-repeat linkers (showing the actual linker DNA
     sequence) compared to mosaic-degenerate linkers.
"""

import csv
import json
import re
import sys
import random
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

from sequence_viewer import (
    classify_bases, render_row, render_full_intron_with_flanks,
    COLOR_BG, COLOR_TXT, add_legend,
)

random.seed(42)


def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}


# ------------------------------------------------------------------
# Loaders
# ------------------------------------------------------------------

def load_tara_telotrons():
    """Yield Tara telotron records with intron_seq + flanks."""
    path = Path('real_telotrons/_tara_oceans_telotrons.tsv')
    if not path.exists(): return []
    out = []
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            if row.get('intron_seq'):
                row['_repeat_family'] = 'TTAGGG'
                row['_dataset'] = 'tara'
                out.append(row)
    return out


def load_eimeria_telotrons():
    """Eimeria telotrons from ULTRA results."""
    out = []
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_(.+)_telotrons\.tsv', f.name)
        if not m: continue
        acc = m.group(1)
        species = m.group(2)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                if row.get('intron_seq'):
                    row['acc'] = acc
                    row['species'] = species
                    row['_repeat_family'] = 'TTTAGGG'
                    row['_dataset'] = 'eimeria'
                    out.append(row)
    return out


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
    """Returns (cls_name, fwd_runs, rev_runs, linker_range_or_None)."""
    if family == 'TTAGGG':
        FWD, REV = _rotations('TTAGGG'), _rotations('CCCTAA')
    elif family == 'TTTAGGG':
        FWD, REV = _rotations('TTTAGGG'), _rotations('CCCTAAA')
    else:
        return 'NONE', [], [], None
    fwd_runs = find_runs(seq, FWD)
    rev_runs = find_runs(seq, REV)
    n = len(seq)
    if not fwd_runs and not rev_runs: return 'NONE', [], [], None
    fbp = sum(e-s for s,e in fwd_runs); rbp = sum(e-s for s,e in rev_runs)
    has_f = fbp >= 10 and fbp/n >= 0.05
    has_r = rbp >= 10 and rbp/n >= 0.05
    if has_f and not has_r: return 'SINGLE_G', fwd_runs, rev_runs, None
    if has_r and not has_f: return 'SINGLE_C', fwd_runs, rev_runs, None
    if not has_f and not has_r: return 'NONE', fwd_runs, rev_runs, None

    fwd_max = max(fwd_runs, key=lambda r: r[1]-r[0])
    rev_max = max(rev_runs, key=lambda r: r[1]-r[0])
    if fwd_max[0] < rev_max[0]:
        gap = rev_max[0] - fwd_max[1]
        if gap <= 4:
            return 'CONVERGING_NOLINKER', fwd_runs, rev_runs, None
        return 'CONVERGING_WITHLINKER', fwd_runs, rev_runs, (fwd_max[1], rev_max[0])
    else:
        gap = fwd_max[0] - rev_max[1]
        if gap <= 4:
            return 'DIVERGING_NOLINKER', fwd_runs, rev_runs, None
        return 'DIVERGING_WITHLINKER', fwd_runs, rev_runs, (rev_max[1], fwd_max[0])


def has_telo(seq, family):
    if family == 'TTAGGG':
        kmers = _rotations('TTAGGG') | _rotations('CCCTAA')
    else:
        kmers = _rotations('TTTAGGG') | _rotations('CCCTAAA')
    return any(k in seq for k in kmers)


def find_tsd(left, right, min_len=3, max_len=20):
    if not left or not right: return 0
    L = min(max_len, len(left), len(right))
    for tl in range(L, min_len - 1, -1):
        if left[-tl:] == right[:tl]:
            return tl
    return 0


# ------------------------------------------------------------------
# FIGURE 1: All orientation classes — sequence examples
# ------------------------------------------------------------------

def figure_all_classes(out_path='real_telotrons/seqfig_all_classes.png',
                        n_per_class=4, max_intron_display=180,
                        flank_bp=40):
    """5 classes × N examples each, side by side, showing actual sequences."""

    tara = load_tara_telotrons()
    eim = load_eimeria_telotrons()
    print(f"Loaded {len(tara)} Tara + {len(eim)} Eimeria")

    # Classify each
    by_class = defaultdict(list)
    for r in tara + eim:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        n = len(seq)
        if n < 30 or n > max_intron_display * 2: continue
        family = r['_repeat_family']
        cls, fwd_runs, rev_runs, linker_range = classify_orientation(seq, family)
        if cls == 'NONE': continue
        # Need flanks
        left = r.get('left_flank_100bp', '').upper()
        right = r.get('right_flank_100bp', '').upper()
        if not left or not right:
            # Tara has flanks; Eimeria ULTRA results don't
            continue
        by_class[cls].append({
            'rec': r, 'seq': seq, 'left': left, 'right': right,
            'family': family, 'cls': cls,
            'fwd_runs': fwd_runs, 'rev_runs': rev_runs,
            'linker_range': linker_range,
        })

    # Print available counts
    print(f"\nAvailable per class:")
    for c, v in by_class.items():
        print(f"  {c}: {len(v)}")

    classes_show = ['SINGLE_G', 'SINGLE_C',
                     'CONVERGING_NOLINKER', 'CONVERGING_WITHLINKER',
                     'DIVERGING_WITHLINKER']

    # Pick examples — prefer short and representative
    selected = {}
    for c in classes_show:
        candidates = by_class.get(c, [])
        if not candidates:
            selected[c] = []
            continue
        # Sort: short introns first, but with clear arrays
        candidates.sort(key=lambda x: (len(x['seq']), -len(x['fwd_runs']) - len(x['rev_runs'])))
        # For converging/diverging, need linker examples
        if 'WITHLINKER' in c:
            # Filter to linker length 20-80bp for visibility
            candidates = [x for x in candidates
                          if x['linker_range'] and
                          20 <= (x['linker_range'][1] - x['linker_range'][0]) <= 80]
        selected[c] = candidates[:n_per_class]
        print(f"  Selected {len(selected[c])} for {c}")

    # Build figure
    n_rows_per_class = n_per_class
    n_classes = len(classes_show)
    total_rows = n_classes * (n_rows_per_class + 1)  # +1 for class title
    max_intron = max((len(x['seq']) for cls in selected.values() for x in cls), default=100)
    total_w = flank_bp + 1 + max_intron + 1 + flank_bp + 5  # +5 for label area

    char_w = 1.0
    fig_h = total_rows * 0.32
    fig_w = max(15, total_w * 0.05)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-15, total_w + 5)
    ax.set_ylim(total_rows + 2, -2)
    ax.axis('off')

    cur_y = 0
    for cls in classes_show:
        # Class title
        ax.text(-15, cur_y + 0.5, cls, ha='left', va='center',
                fontsize=11, fontweight='bold',
                color={'SINGLE_G': '#2050cc',
                        'SINGLE_C': '#cc2020',
                        'CONVERGING_NOLINKER': '#206020',
                        'CONVERGING_WITHLINKER': '#206020',
                        'DIVERGING_WITHLINKER': '#502080'}.get(cls, 'black'))
        n_in_class = len([x for x in by_class.get(cls, [])])
        ax.text(total_w - 5, cur_y + 0.5, f"(n={n_in_class:,} total)",
                ha='right', va='center', fontsize=9, fontstyle='italic', color='gray')
        cur_y += 1

        for ex in selected[cls]:
            seq = ex['seq']
            left = ex['left'][-flank_bp:]
            right = ex['right'][:flank_bp]
            family = ex['family']

            # Build label
            r = ex['rec']
            sp = r.get('species', r.get('mag', '?'))
            sp_short = sp.split('_')[0] + '. ' + sp.split('_')[-1] if '_' in sp else sp
            lbl = f"{sp_short[:18]} {r.get('contig','')[:18]}:{r.get('start','')} ({len(seq)}bp)"

            # Highlights: linker range, TSDs
            highlights = []
            if ex.get('linker_range'):
                ls, le = ex['linker_range']
                # We could highlight the linker itself but it'd just be 'non_telo' anyway
                # Instead we'll let the colors show what's repeat vs not
                pass

            # TSDs
            tsd_len = find_tsd(left, right)

            # Render flanks with TSDs
            x_start = 0
            # Left flank
            lf_classes = ['flank'] * len(left)
            tsd_left_overrides = []
            tsd_right_overrides = []
            if tsd_len > 0:
                tsd_left_overrides = [(len(left) - tsd_len, len(left), 'tsd')]
                tsd_right_overrides = [(0, tsd_len, 'tsd')]
            render_row(ax, cur_y, left, lf_classes, x_start=x_start,
                       char_width=char_w,
                       highlight_ranges=tsd_left_overrides,
                       show_text=True, fontsize=2.0)
            x_start += len(left) * char_w

            # 5' splice site
            ax.add_patch(Rectangle((x_start, cur_y), char_w * 0.5, 1.0,
                                    facecolor='black', edgecolor='none'))
            x_start += char_w * 0.5

            # Intron
            intron_classes = classify_bases(seq, family)
            render_row(ax, cur_y, seq, intron_classes, x_start=x_start,
                       char_width=char_w, show_text=True, fontsize=2.0)
            x_start += len(seq) * char_w

            # 3' splice site
            ax.add_patch(Rectangle((x_start, cur_y), char_w * 0.5, 1.0,
                                    facecolor='black', edgecolor='none'))
            x_start += char_w * 0.5

            # Right flank
            rf_classes = ['flank'] * len(right)
            render_row(ax, cur_y, right, rf_classes, x_start=x_start,
                       char_width=char_w,
                       highlight_ranges=tsd_right_overrides,
                       show_text=True, fontsize=2.0)

            # Label
            ax.text(-15, cur_y + 0.5, lbl, ha='left', va='center',
                    fontsize=6.5, family='sans-serif')

            cur_y += 1

        cur_y += 0.5  # gap between classes

    # Legend at top
    legend_y = -1
    legend_items = [
        ('fwd_canonical', 'TTAGGG/TTTAGGG (G-strand)'),
        ('rev_canonical', 'CCCTAA/CCCTAAA (C-strand)'),
        ('fwd_variant', '1-mm variant fwd'),
        ('rev_variant', '1-mm variant rev'),
        ('flank', 'Flanking exon'),
        ('non_telo', 'Non-repeat (linker)'),
        ('tsd', 'TSD (target site dup)'),
    ]
    legend_x = 0
    for cls_, lbl in legend_items:
        ax.add_patch(Rectangle((legend_x, legend_y), 5, 0.8,
                                facecolor=COLOR_BG[cls_], edgecolor='black', linewidth=0.3))
        ax.text(legend_x + 6, legend_y + 0.4, lbl, fontsize=7, va='center')
        legend_x += 30

    # Black bar legend
    ax.add_patch(Rectangle((legend_x, legend_y), 0.5, 0.8,
                            facecolor='black', edgecolor='none'))
    ax.text(legend_x + 1.5, legend_y + 0.4, 'Splice site', fontsize=7, va='center')

    plt.suptitle('Sequence-level examples of telotron orientation classes',
                  fontsize=14, fontweight='bold')
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out_path}")


# ------------------------------------------------------------------
# FIGURE 2: Shared-linker tandem cluster pairs
# ------------------------------------------------------------------

def figure_shared_linkers(out_path='real_telotrons/seqfig_shared_linkers.png'):
    """Show the 5 specific shared-linker pairs."""

    # The 5 shared cases — extract from BLAST results
    eim = load_eimeria_telotrons()
    eim_by_locus = {}
    for r in eim:
        try:
            key = (r.get('acc', ''), r['contig'], int(r['start']), int(r['end']))
            eim_by_locus[key] = r
        except: pass

    # Shared-linker pairs from previous analysis
    pairs = [
        # (acc, src_contig, src_start, src_end, hit_contig, hit_start, hit_end, len, pident)
        ('GCF_000499385.1', 'NW_013654861.1', 79248, 79326,
         'NW_013655283.1', 269, 346, 78, 100.0),
        ('GCF_000499385.1', 'NW_013651890.1', 1453, 1534,
         'NW_013651890.1', 902, 975, 74, 100.0),
        ('GCF_000499605.1', 'NW_013564060.1', 17301, 17370,
         'NW_013564060.1', 17580, 17641, 62, 98.4),
        ('GCF_000499425.1', 'NW_013549356.1', 7752, 7806,
         'NW_013549356.1', 7479, 7524, 46, 100.0),
        ('GCF_000499545.2', 'NW_013545172.1', 1021, 1050,
         'NW_013542093.1', 1229, 1256, 28, 92.9),
    ]

    # Find each telotron's full sequence
    fig, axes = plt.subplots(len(pairs), 1, figsize=(20, len(pairs) * 4))
    if len(pairs) == 1: axes = [axes]

    for idx, (acc, c1, s1, e1, c2, s2, e2, hit_len, pident) in enumerate(pairs):
        ax = axes[idx]
        ax.axis('off')

        # Get telotron 1 (source)
        r1 = eim_by_locus.get((acc, c1, s1, e1))
        if not r1:
            # try with strip
            for k, v in eim_by_locus.items():
                if k[0] == acc and k[1] == c1 and abs(k[2] - s1) < 5 and abs(k[3] - e1) < 5:
                    r1 = v
                    break

        # Get telotron 2 (paralog) — find by matching genomic position to a known telotron
        r2 = None
        for k, v in eim_by_locus.items():
            if k[0] == acc and k[1] == c2:
                if k[2] <= s2 + 50 and k[3] >= e2 - 50:
                    r2 = v
                    break

        if not r1 or not r2:
            ax.text(0, 0, f"Pair {idx+1}: telotron records not found ({acc} {c1}:{s1}-{e1} <-> {c2}:{s2}-{e2})",
                    fontsize=10, ha='left', va='top')
            continue

        seq1 = r1['intron_seq'].upper()
        seq2 = r2['intron_seq'].upper()

        # Find the shared region in each
        # We know from BLAST: in seq1, the linker is at some position; in seq2, the paralog hits at offset
        # Shared region in seq1: linker positions
        # Shared region in seq2: positions corresponding to hit start-end relative to telotron start

        # For seq2, the hit's genomic position s2-e2 maps to seq2 position (s2 - r2_start) to (e2 - r2_start)
        s2_offset = s2 - int(r2['start'])
        e2_offset = e2 - int(r2['start'])
        s2_offset = max(0, min(len(seq2)-1, s2_offset))
        e2_offset = max(0, min(len(seq2), e2_offset))

        # For seq1, the source linker = need to find where the linker is in the seq.
        # The whole intron seq1 is what we BLASTed — but only the linker region produced the hit.
        # We need to identify the linker position in seq1.
        cls1, fr, rr, linker_range = classify_orientation(seq1, 'TTTAGGG')
        if linker_range:
            ls1, le1 = linker_range
        else:
            # Fall back: assume center of intron
            ls1, le1 = len(seq1) // 3, 2 * len(seq1) // 3

        # Render
        char_w = 1.0
        max_w = max(len(seq1), len(seq2)) + 100  # space for flanks
        ax.set_xlim(-15, max_w + 10)
        ax.set_ylim(8, -2)

        # Title
        ax.text(0, -1.5, f"Pair {idx+1}: {acc} | shared 60-80bp linker | "
                          f"{hit_len}bp BLAST hit at {pident:.1f}% identity",
                fontsize=9, fontweight='bold')

        # Telotron 1
        flank_bp = 30
        left1 = r1.get('left_flank_100bp', '').upper()[-flank_bp:]
        right1 = r1.get('right_flank_100bp', '').upper()[:flank_bp]

        # If flanks are empty for Eimeria ULTRA records, that's fine — render without
        x = 0
        if left1:
            render_row(ax, 1, left1, ['flank']*len(left1), x_start=x, char_width=char_w, fontsize=3)
            x += len(left1) * char_w
        ax.add_patch(Rectangle((x, 1), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5

        # Render seq1 with linker highlighted in green (shared)
        cls = classify_bases(seq1, 'TTTAGGG')
        # Highlight the linker region as 'shared'
        highlights = []
        if ls1 < le1:
            highlights = [(ls1, le1, 'shared')]
        render_row(ax, 1, seq1, cls, x_start=x, char_width=char_w,
                   highlight_ranges=highlights, fontsize=3)
        x += len(seq1) * char_w
        ax.add_patch(Rectangle((x, 1), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        if right1:
            render_row(ax, 1, right1, ['flank']*len(right1), x_start=x, char_width=char_w, fontsize=3)

        ax.text(-15, 1.5, f"{c1}:{s1}-{e1}\n({len(seq1)}bp)",
                ha='left', va='center', fontsize=7)

        # Telotron 2
        x = 0
        left2 = r2.get('left_flank_100bp', '').upper()[-flank_bp:]
        right2 = r2.get('right_flank_100bp', '').upper()[:flank_bp]
        if left2:
            render_row(ax, 4, left2, ['flank']*len(left2), x_start=x, char_width=char_w, fontsize=3)
            x += len(left2) * char_w
        ax.add_patch(Rectangle((x, 4), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5

        cls = classify_bases(seq2, 'TTTAGGG')
        highlights2 = [(s2_offset, e2_offset, 'shared')]
        render_row(ax, 4, seq2, cls, x_start=x, char_width=char_w,
                   highlight_ranges=highlights2, fontsize=3)
        x += len(seq2) * char_w
        ax.add_patch(Rectangle((x, 4), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        if right2:
            render_row(ax, 4, right2, ['flank']*len(right2), x_start=x, char_width=char_w, fontsize=3)

        ax.text(-15, 4.5, f"{c2}:{s2}-{e2}\n({len(seq2)}bp)",
                ha='left', va='center', fontsize=7)

        # Draw an arrow / line between the two
        # Show linker sequence text on the side
        if ls1 < le1:
            shared_seq1 = seq1[ls1:le1]
            shared_seq2 = seq2[s2_offset:e2_offset]
            ax.text(max_w + 2, 2.5, f"Shared:\n5' {shared_seq1[:40]}\n   {shared_seq2[:40]}",
                    fontsize=7, family='monospace', va='center')

    plt.suptitle('Tandem cluster pairs: linker shared between adjacent Eimeria telotrons',
                  fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out_path}")


# ------------------------------------------------------------------
# FIGURE 3: Haptophyte CONVERGING vs Eimeria DIVERGING
# ------------------------------------------------------------------

def figure_haptophyte_vs_eimeria(out_path='real_telotrons/seqfig_hapto_vs_eimeria.png',
                                   n_each=8):
    """Side-by-side examples showing the architecture difference."""
    tara = load_tara_telotrons()
    eim = load_eimeria_telotrons()

    # Find Tara CONVERGING_WITHLINKER examples
    hapt_examples = []
    for r in tara:
        seq = r.get('intron_seq', '').upper()
        left = r.get('left_flank_100bp', '').upper()
        right = r.get('right_flank_100bp', '').upper()
        if not seq or not left or not right: continue
        if len(seq) < 80 or len(seq) > 200: continue
        cls, fr, rr, link = classify_orientation(seq, 'TTAGGG')
        if cls == 'CONVERGING_WITHLINKER' and link and 15 <= (link[1]-link[0]) <= 60:
            hapt_examples.append({'seq': seq, 'left': left, 'right': right, 'family': 'TTAGGG',
                                  'rec': r, 'cls': cls, 'link': link})
    random.shuffle(hapt_examples)
    hapt_examples = hapt_examples[:n_each]
    print(f"Selected {len(hapt_examples)} Tara CONVERGING examples")

    # Find Eimeria DIVERGING_WITHLINKER examples (need flanks)
    # Eimeria ULTRA results don't have flanks — extract from genome on the fly
    contig_lens = {}
    print("Loading Eimeria genomes for flank extraction...")
    for f in Path('genomes/GCF_000499385.1').rglob('*.fna'):
        cur, parts = None, []
        with open(f) as fh:
            for line in fh:
                if line.startswith('>'):
                    if cur: contig_lens[cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]
                    parts = []
                else:
                    parts.append(line.strip().upper())
            if cur: contig_lens[cur] = ''.join(parts)

    eim_examples = []
    for r in eim:
        if r.get('acc') != 'GCF_000499385.1': continue
        seq = r.get('intron_seq', '').upper()
        if not seq or len(seq) < 80 or len(seq) > 200: continue
        if r['contig'] not in contig_lens: continue
        contig_seq = contig_lens[r['contig']]
        try:
            s = int(r['start']); e = int(r['end'])
        except: continue
        left = contig_seq[max(0, s-100):s]
        right = contig_seq[e:min(len(contig_seq), e+100)]
        if not left or not right: continue
        cls, fr, rr, link = classify_orientation(seq, 'TTTAGGG')
        if cls == 'DIVERGING_WITHLINKER' and link and 15 <= (link[1]-link[0]) <= 60:
            eim_examples.append({'seq': seq, 'left': left, 'right': right, 'family': 'TTTAGGG',
                                  'rec': r, 'cls': cls, 'link': link})
    random.shuffle(eim_examples)
    eim_examples = eim_examples[:n_each]
    print(f"Selected {len(eim_examples)} Eimeria DIVERGING examples")

    # Render two columns
    char_w = 1.0
    max_intron = max((len(x['seq']) for x in hapt_examples + eim_examples), default=200)
    flank_bp = 40
    col_w = flank_bp + 1 + max_intron + 1 + flank_bp
    n_rows = max(len(hapt_examples), len(eim_examples))
    total_h = n_rows + 4

    fig, ax = plt.subplots(figsize=(28, max(8, n_rows * 0.6)))
    ax.set_xlim(-25, 2 * col_w + 30)
    ax.set_ylim(total_h + 1, -3)
    ax.axis('off')

    # Headers
    ax.text(col_w / 2, -2, "HAPTOPHYTES (Tara) — CONVERGING_WITHLINKER",
            ha='center', va='center', fontsize=13, fontweight='bold', color='steelblue')
    ax.text(col_w + 30 + col_w / 2, -2, "EIMERIA — DIVERGING_WITHLINKER",
            ha='center', va='center', fontsize=13, fontweight='bold', color='coral')

    ax.text(col_w / 2, -1, "5'-(TTAGGG)→...←(CCCTAA)-3'  G-strand points INWARD",
            ha='center', va='center', fontsize=9, color='steelblue', fontstyle='italic')
    ax.text(col_w + 30 + col_w / 2, -1, "5'-(CCCTAAA)←...→(TTTAGGG)-3'  G-strand points OUTWARD",
            ha='center', va='center', fontsize=9, color='coral', fontstyle='italic')

    # Render haptophyte column
    for i, ex in enumerate(hapt_examples):
        seq = ex['seq']; left = ex['left'][-flank_bp:]; right = ex['right'][:flank_bp]
        x = 0
        # Left flank
        render_row(ax, i, left, ['flank']*len(left), x_start=x, char_width=char_w, fontsize=3)
        x += len(left) * char_w
        ax.add_patch(Rectangle((x, i), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        # Intron
        cls_ = classify_bases(seq, 'TTAGGG')
        render_row(ax, i, seq, cls_, x_start=x, char_width=char_w, fontsize=3)
        x += len(seq) * char_w
        ax.add_patch(Rectangle((x, i), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        render_row(ax, i, right, ['flank']*len(right), x_start=x, char_width=char_w, fontsize=3)
        # Label
        sp = ex['rec'].get('mag', 'TARA').split('_')[1] if '_' in ex['rec'].get('mag', '') else 'TARA'
        ax.text(-25, i + 0.5, f"{sp} | {len(seq)}bp", ha='left', va='center', fontsize=6)

    # Render eimeria column
    eim_x_offset = col_w + 30
    for i, ex in enumerate(eim_examples):
        seq = ex['seq']; left = ex['left'][-flank_bp:]; right = ex['right'][:flank_bp]
        x = eim_x_offset
        render_row(ax, i, left, ['flank']*len(left), x_start=x, char_width=char_w, fontsize=3)
        x += len(left) * char_w
        ax.add_patch(Rectangle((x, i), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        cls_ = classify_bases(seq, 'TTTAGGG')
        render_row(ax, i, seq, cls_, x_start=x, char_width=char_w, fontsize=3)
        x += len(seq) * char_w
        ax.add_patch(Rectangle((x, i), char_w * 0.5, 1.0, facecolor='black'))
        x += char_w * 0.5
        render_row(ax, i, right, ['flank']*len(right), x_start=x, char_width=char_w, fontsize=3)
        ax.text(eim_x_offset - 5, i + 0.5,
                f"{ex['rec'].get('contig','')[:14]}:{ex['rec'].get('start','')}",
                ha='left', va='center', fontsize=6)

    plt.suptitle('Sequence-level: haptophyte CONVERGING vs Eimeria DIVERGING telotrons',
                  fontsize=14, fontweight='bold')
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out_path}")


# ------------------------------------------------------------------
# FIGURE 4: Linker examples - pure non-repeat vs degenerate
# ------------------------------------------------------------------

def figure_linker_examples(out_path='real_telotrons/seqfig_linkers.png',
                            n_per_class=8):
    """Show actual linker sequences side-by-side from different categories."""

    eim = load_eimeria_telotrons()

    pure_eim = []   # Eimeria pure non-repeat linkers
    mosaic_eim = []  # Eimeria linkers with telomeric content
    pure_tara = []   # Tara pure non-repeat linkers
    mosaic_tara = []

    for r in eim:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, fr, rr, link = classify_orientation(seq, 'TTTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if len(linker) < 20 or len(linker) > 80: continue
        if has_telo(linker, 'TTTAGGG'):
            mosaic_eim.append({'seq': seq, 'linker': linker, 'family': 'TTTAGGG',
                                'cls': cls, 'rec': r, 'linker_range': link})
        else:
            pure_eim.append({'seq': seq, 'linker': linker, 'family': 'TTTAGGG',
                              'cls': cls, 'rec': r, 'linker_range': link})

    tara = load_tara_telotrons()
    for r in tara:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, fr, rr, link = classify_orientation(seq, 'TTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if len(linker) < 20 or len(linker) > 80: continue
        if has_telo(linker, 'TTAGGG'):
            mosaic_tara.append({'seq': seq, 'linker': linker, 'family': 'TTAGGG',
                                  'cls': cls, 'rec': r, 'linker_range': link})
        else:
            pure_tara.append({'seq': seq, 'linker': linker, 'family': 'TTAGGG',
                                'cls': cls, 'rec': r, 'linker_range': link})

    random.shuffle(pure_eim); random.shuffle(mosaic_eim)
    random.shuffle(pure_tara); random.shuffle(mosaic_tara)

    pure_eim = pure_eim[:n_per_class]
    mosaic_eim = mosaic_eim[:n_per_class]
    pure_tara = pure_tara[:n_per_class]
    mosaic_tara = mosaic_tara[:n_per_class]

    print(f"Linker examples: pure_tara={len(pure_tara)}, mosaic_tara={len(mosaic_tara)}, "
          f"pure_eim={len(pure_eim)}, mosaic_eim={len(mosaic_eim)}")

    # Layout: 4 sections × n_per_class rows
    sections = [
        ('Tara haptophytes — DEGENERATE telomeric linker (mosaic, 91% of cases)', mosaic_tara, 'TTAGGG'),
        ('Tara haptophytes — TRULY non-repeat linker (8.6% of cases)', pure_tara, 'TTAGGG'),
        ('Eimeria — DEGENERATE telomeric linker (mosaic)', mosaic_eim, 'TTTAGGG'),
        ('Eimeria — TRULY non-repeat linker (25.9% of CONVERGING_WITHLINKER)', pure_eim, 'TTTAGGG'),
    ]
    total_rows = sum(len(s[1]) + 1 for s in sections)

    char_w = 1.0
    max_intron = 200
    fig_h = max(10, total_rows * 0.4)

    fig, ax = plt.subplots(figsize=(22, fig_h))
    ax.set_xlim(-30, max_intron + 20)
    ax.set_ylim(total_rows + 1, -1)
    ax.axis('off')

    cur_y = 0
    for title, examples, family in sections:
        ax.text(-30, cur_y + 0.4, title,
                ha='left', va='center', fontsize=10, fontweight='bold')
        cur_y += 1
        for ex in examples:
            seq = ex['seq']
            ls, le = ex['linker_range']

            # Render full intron with linker highlighted
            cls_ = classify_bases(seq, family)
            # No need to highlight differently since classification already shows
            x = 0
            render_row(ax, cur_y, seq, cls_, x_start=x, char_width=char_w, fontsize=3)
            # Mark linker boundaries with small ticks
            ax.plot([ls, ls], [cur_y, cur_y + 1], color='black', linewidth=0.5, linestyle=':')
            ax.plot([le, le], [cur_y, cur_y + 1], color='black', linewidth=0.5, linestyle=':')
            ax.text(-30, cur_y + 0.5,
                    f"{len(seq)}bp intron, {len(ex['linker'])}bp linker",
                    ha='left', va='center', fontsize=6)
            cur_y += 1
        cur_y += 0.5

    plt.suptitle('Linker sequences: degenerate telomeric vs truly non-repeat',
                  fontsize=14, fontweight='bold')
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == '__main__':
    print("Generating sequence figures...\n")
    figure_all_classes()
    print()
    figure_shared_linkers()
    print()
    figure_haptophyte_vs_eimeria()
    print()
    figure_linker_examples()
    print("\nAll figures generated.")
