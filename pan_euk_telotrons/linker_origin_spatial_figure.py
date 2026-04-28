#!/usr/bin/env python3
"""
Spatial summary: for each linker→origin pair, render the actual DNA window
around the origin with full base-resolution characters. Each row = one origin
with 50bp left + hit + 50bp right of genomic flanking sequence.
"""
import json, re
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

from sequence_viewer import classify_bases, render_row, COLOR_BG


def main():
    with open('real_telotrons/linker_origin_pairs.json') as f:
        pairs = json.load(f)

    # Sort by class for visual grouping
    def sort_key(p):
        same_contig = p['source_contig'] == p['origin_contig']
        in_telo = bool(p.get('origin_in_telotron'))
        if in_telo: return (0, -p['pident'])
        if same_contig: return (1, -p['pident'])
        return (2, -p['pident'])

    pairs.sort(key=sort_key)
    for i, p in enumerate(pairs):
        p['display_idx'] = i + 1

    # Compute per-pair distance string
    for p in pairs:
        same_contig = p['source_contig'] == p['origin_contig']
        if same_contig:
            src_mid = (p['source_linker_start'] + p['source_linker_end']) // 2
            org_mid = (p['origin_start'] + p['origin_end']) // 2
            p['_dist_label'] = f"{abs(org_mid - src_mid):,}bp on same contig"
        else:
            p['_dist_label'] = f"different contig"

        if p.get('origin_in_telotron'):
            p['_class_label'] = 'SHARED-LINKER (in another telotron)'
            p['_class_color'] = 'darkgreen'
        elif same_contig:
            p['_class_label'] = 'SAME-CONTIG-INTERGENIC'
            p['_class_color'] = 'darkblue'
        else:
            p['_class_label'] = 'CROSS-CONTIG'
            p['_class_color'] = 'darkred'

    n = len(pairs)
    char_w = 1.0
    flank_bp = 50
    max_aln = max(p['aln_length'] for p in pairs)
    row_w = flank_bp + max_aln + flank_bp + 50  # +50 for label area

    # Each row gets 2 lines of text + sequence
    PAIR_HEIGHT = 3
    fig_h = max(15, n * PAIR_HEIGHT * 0.5)
    fig_w = max(20, row_w * 0.05)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-50, row_w + 5)
    ax.set_ylim(n * PAIR_HEIGHT + 3, -3)
    ax.axis('off')

    plt.suptitle('Linker → genomic origin: actual DNA at each hit position\n'
                 'All 24 Eimeria paralog hits with full sequence context',
                 fontsize=14, fontweight='bold')

    for idx, p in enumerate(pairs):
        y_base = idx * PAIR_HEIGHT
        origin_full = p['origin_with_flanks']

        # Hit position within origin_full
        hit_start_in_full = p.get('origin_flank_left_bp', 100)
        hit_end_in_full = hit_start_in_full + p['aln_length']

        # Show 50bp flank + hit + 50bp flank
        flank_l_actual = min(flank_bp, hit_start_in_full)
        flank_r_actual = min(flank_bp, len(origin_full) - hit_end_in_full)

        display_start = hit_start_in_full - flank_l_actual
        display_end = hit_end_in_full + flank_r_actual
        display_seq = origin_full[display_start:display_end]

        # Hit positions within display_seq
        hit_start_in_disp = flank_l_actual
        hit_end_in_disp = flank_l_actual + p['aln_length']

        # Header line
        cls_label = p['_class_label']
        cls_color = p['_class_color']
        ax.text(-50, y_base + 0.5,
                f"#{p['display_idx']} {cls_label}",
                ha='left', va='center', fontsize=9, fontweight='bold', color=cls_color)
        ax.text(50, y_base + 0.5,
                f"{p['acc']} | {p['origin_contig']}:{p['origin_start']}-{p['origin_end']} | "
                f"src→origin: {p['_dist_label']} | "
                f"end-distance: {p['origin_dist_to_end']:,}bp | "
                f"BLAST: {p['pident']:.0f}% over {p['aln_length']}bp",
                ha='left', va='center', fontsize=7, color='gray')

        # Render the actual DNA sequence with hit region highlighted
        # Classify based on whether origin is in another telotron
        if p.get('origin_in_telotron'):
            classes = classify_bases(display_seq, 'TTTAGGG')
        else:
            # Mostly intergenic — flag flanks as 'flank', and the hit region as 'shared'
            classes = ['flank'] * len(display_seq)

        # Highlight: hit region as 'shared'
        highlights = [(hit_start_in_disp, hit_end_in_disp, 'shared')]

        render_row(ax, y_base + 1, display_seq, classes,
                   x_start=0, char_width=char_w,
                   highlight_ranges=highlights,
                   show_text=True, fontsize=5.5)

        # Position labels
        ax.text(0, y_base + 2.2, f"-{flank_l_actual}bp",
                ha='left', va='center', fontsize=6, color='gray')
        ax.text(hit_start_in_disp + (hit_end_in_disp - hit_start_in_disp) / 2, y_base + 2.2,
                f"BLAST hit ({p['aln_length']}bp)",
                ha='center', va='center', fontsize=7, fontweight='bold', color='darkgreen')
        ax.text(hit_end_in_disp + flank_r_actual, y_base + 2.2,
                f"+{flank_r_actual}bp",
                ha='right', va='center', fontsize=6, color='gray')

        # Mark contig-end distance with a special bar if very close
        if p['origin_start'] - flank_l_actual < 100:
            # Very close to start of contig
            ax.plot([0 - 0.5, 0 - 0.5], [y_base + 1, y_base + 2],
                    color='red', linewidth=2, linestyle=':')
            ax.text(-3, y_base + 1.5, "5'\nedge",
                    ha='right', va='center', fontsize=7, fontweight='bold', color='red')

    # Legend
    legend_y = n * PAIR_HEIGHT + 1
    items = [
        ('shared',          'BLAST hit region'),
        ('fwd_canonical',   'TTTAGGG (when origin in telotron)'),
        ('rev_canonical',   'CCCTAAA (when origin in telotron)'),
        ('fwd_variant',     '1-mm fwd variant'),
        ('rev_variant',     '1-mm rev variant'),
        ('flank',           'Surrounding genomic context'),
    ]
    x = 0
    for cls, lbl in items:
        ax.add_patch(Rectangle((x, legend_y), 4, 1, facecolor=COLOR_BG[cls], edgecolor='black', linewidth=0.3))
        ax.text(x + 5, legend_y + 0.5, lbl, fontsize=8, va='center')
        x += 50

    out = 'real_telotrons/seqfig_linker_origin_spatial.png'
    plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
