#!/usr/bin/env python3
"""
Render: for each linker→origin pair, show:
  Row 1 (top):    Source telotron with linker highlighted in green
  Row 2 (middle): Linker sequence aligned to origin (mismatches in red)
  Row 3 (bottom): Origin location with genomic flanks
                  + classification (subtelomeric / in another telotron / other)

One large figure with all 24 pairs, each pair is a 3-row mini-panel.
"""

import json
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.font_manager import FontProperties
from pathlib import Path

from sequence_viewer import classify_bases, render_row, COLOR_BG, COLOR_TXT


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def smith_waterman_simple(s1, s2, match=2, mismatch=-1, gap=-2):
    """Trivial pairwise local alignment for short sequences (< 200bp).
    Returns (aligned_s1, aligned_s2, mismatches_array)."""
    # For our case the sequences should be very similar (BLAST said >85% identity)
    # so a global alignment with no gaps is fine for visualization.
    # Just match positions directly.
    L = min(len(s1), len(s2))
    aligned_s1 = s1[:L]
    aligned_s2 = s2[:L]
    mismatches = [a != b for a, b in zip(aligned_s1, aligned_s2)]
    return aligned_s1, aligned_s2, mismatches


def render_alignment_row(ax, y, s1, s2, x_start=0, char_w=1.0, fontsize=5.5):
    """Render a base-paired alignment row with mismatches highlighted in red."""
    L = min(len(s1), len(s2))
    mismatches = [a != b for a, b in zip(s1[:L], s2[:L])]
    n_match = sum(1 for m in mismatches if not m)

    # Top row: s1
    for j in range(L):
        x = x_start + j * char_w
        bg = (0.85, 0.92, 1.0, 1.0) if not mismatches[j] else (1.0, 0.7, 0.7, 1.0)
        ax.add_patch(Rectangle((x, y), char_w, 1.0, facecolor=bg, edgecolor='none'))
    mono = FontProperties(family='monospace', size=fontsize)
    for j in range(L):
        x = x_start + j * char_w + char_w / 2
        ax.text(x, y + 0.5, s1[j], ha='center', va='center',
                fontproperties=mono, color='black' if not mismatches[j] else (0.5, 0, 0))

    # Connector row: middle line of | for matches, space for mismatches
    for j in range(L):
        x = x_start + j * char_w + char_w / 2
        if not mismatches[j]:
            ax.text(x, y - 0.4, '|', ha='center', va='center',
                    fontproperties=mono, color='gray')
        else:
            ax.text(x, y - 0.4, '·', ha='center', va='center',
                    fontproperties=mono, color='red')

    # Bottom row: s2
    for j in range(L):
        x = x_start + j * char_w
        bg = (0.92, 1.0, 0.85, 1.0) if not mismatches[j] else (1.0, 0.7, 0.7, 1.0)
        ax.add_patch(Rectangle((x, y - 1.0), char_w, 1.0, facecolor=bg, edgecolor='none'))
    for j in range(L):
        x = x_start + j * char_w + char_w / 2
        ax.text(x, y - 0.5, s2[j], ha='center', va='center',
                fontproperties=mono, color='black' if not mismatches[j] else (0.5, 0, 0))

    return n_match, L


def compute_distance_info(p, contig_telotrons):
    """Returns dict with class, distance_str, and other distance info."""
    same_contig = p['source_contig'] == p['origin_contig']
    src_mid = (p['source_linker_start'] + p['source_linker_end']) // 2
    org_mid = (p['origin_start'] + p['origin_end']) // 2

    src_clen = p.get('src_contig_len', 0)
    org_clen = p.get('origin_contig_len', 0)
    src_dist_end = min(p['source_linker_start'], src_clen - p['source_linker_end']) if src_clen else None
    org_dist_end = min(p['origin_start'], org_clen - p['origin_end']) if org_clen else None

    in_diff_telotron = bool(p.get('origin_in_telotron'))

    if in_diff_telotron:
        cls = 'shared-linker'
    elif same_contig:
        cls = 'same-contig-intergenic'
    else:
        cls = 'cross-contig'

    if same_contig:
        dist = abs(org_mid - src_mid)
        dist_str = f"{dist:,}bp apart on same contig"
    else:
        dist_str = (f"src {src_dist_end:,}bp from end | "
                    f"origin {org_dist_end:,}bp from end")

    return cls, dist_str, src_dist_end, org_dist_end


def render_pairs(pairs_subset, out_path, title):
    """Render a batch of pairs to one figure."""
    char_w = 1.0
    fontsize = 4.5

    PAIR_HEIGHT = 7
    n_pairs = len(pairs_subset)

    max_seq_len = max(
        max(len(p['linker_seq']),
            len(p.get('src_telotron_intron_seq', '')),
            len(p['origin_with_flanks']))
        for p in pairs_subset
    )
    # Cap at reasonable size
    max_seq_len = min(max_seq_len, 800)

    fig_w = min(28, max(15, max_seq_len * 0.025 + 5))
    fig_h = min(40, max(8, n_pairs * PAIR_HEIGHT * 0.35))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-25, max_seq_len + 50)
    ax.set_ylim(n_pairs * PAIR_HEIGHT + 2, -3)
    ax.axis('off')

    plt.suptitle(title, fontsize=13, fontweight='bold')

    for idx, p in enumerate(pairs_subset):
        y_base = idx * PAIR_HEIGHT
        linker = p['linker_seq']

        # Header — compute class and distance
        cls, dist_str, src_d, org_d = compute_distance_info(p, None)
        cls_color = {
            'shared-linker': 'darkgreen',
            'same-contig-intergenic': 'darkblue',
            'cross-contig': 'darkred',
        }.get(cls, 'black')

        in_telo = "INSIDE another telotron" if p.get('origin_in_telotron') else "intergenic/exon"
        subtelo = f"<5kb from contig end ({p['origin_dist_to_end']:,}bp)" if p['origin_dist_to_end'] < 5000 else f"internal ({p['origin_dist_to_end']:,}bp from end)"

        # Strand info
        strand_info = "fwd" if p['origin_strand'] == '+' else "rev-comp"

        ax.text(-25, y_base + 0.5,
                f"#{p.get('global_idx', idx+1)}", ha='left', va='center', fontsize=11, fontweight='bold')
        ax.text(0, y_base + 0.5,
                f"{p['acc']}: {p['source_contig']}:{p['source_linker_start']}-{p['source_linker_end']} → "
                f"{p['origin_contig']}:{p['origin_start']}-{p['origin_end']} ({strand_info})",
                ha='left', va='center', fontsize=8, fontweight='bold')

        # Distance + class on its own line, color-coded
        ax.text(0, y_base + 1.2,
                f"DISTANCE: {dist_str}  |  CLASS: ",
                ha='left', va='center', fontsize=8, fontweight='bold')
        ax.text(180, y_base + 1.2,
                cls.upper(),
                ha='left', va='center', fontsize=8, fontweight='bold', color=cls_color)

        ax.text(0, y_base + 1.9,
                f"BLAST: {p['pident']:.1f}% identity over {p['aln_length']}bp | linker length: {len(linker)}bp | "
                f"origin: {in_telo}, {subtelo}",
                ha='left', va='center', fontsize=7.5, color='gray', fontstyle='italic')

        # Adjust y offsets to account for the extra header line we added
        # Row 1 (y = y_base + 2.5): source telotron with linker highlighted
        src_seq = p.get('src_telotron_intron_seq', '')
        src_start_telo = p.get('src_telotron_start')
        if src_seq and src_start_telo:
            # Linker offset within source telotron
            link_off_start = p['source_linker_start'] - src_start_telo
            link_off_end = p['source_linker_end'] - src_start_telo

            # Determine repeat family — Eimeria uses TTTAGGG
            cls_per_base = classify_bases(src_seq, 'TTTAGGG')
            # Highlight linker as 'shared' green
            highlights = []
            if 0 <= link_off_start < link_off_end <= len(src_seq):
                highlights = [(link_off_start, link_off_end, 'shared')]

            ax.text(-25, y_base + 2.5, "Source telotron:",
                    ha='left', va='center', fontsize=7)
            # Add splice site bars
            ax.add_patch(Rectangle((-0.5, y_base + 2), 0.5, 1.0,
                                    facecolor='black', edgecolor='none'))
            render_row(ax, y_base + 2, src_seq, cls_per_base,
                       x_start=0, char_width=char_w,
                       highlight_ranges=highlights, fontsize=fontsize)
            ax.add_patch(Rectangle((len(src_seq) * char_w, y_base + 2), 0.5, 1.0,
                                    facecolor='black', edgecolor='none'))
        else:
            ax.text(0, y_base + 2.5, "[source telotron sequence unavailable]",
                    fontsize=8, color='red')

        # Row 2 (y = y_base + 4 with offset): linker vs origin alignment
        # The actual aligned region in origin (just the BLAST hit core)
        origin_full = p['origin_with_flanks']
        origin_core = p['origin_core']
        if not origin_core:
            # Reconstruct from origin_with_flanks
            left_flank_n = p.get('origin_flank_left_bp', 100)
            origin_core_len = p['origin_end'] - p['origin_start'] + 1
            origin_core = origin_full[left_flank_n:left_flank_n + origin_core_len]

        # If origin is on - strand, reverse-complement for alignment view
        if p['origin_strand'] == '-':
            display_origin = revcomp(origin_core)
        else:
            display_origin = origin_core

        ax.text(-25, y_base + 4.5,
                f"Linker → Origin\n({p['pident']:.0f}% id)",
                ha='left', va='center', fontsize=7)

        # Show alignment from qstart to qend in linker
        # Use the qstart/qend positions if available
        qstart = max(0, p.get('qstart', 1) - 1)
        qend = p.get('qend', len(linker))
        linker_aligned = linker[qstart:qend]

        # Trim origin to match alignment length
        L = min(len(linker_aligned), len(display_origin))
        s1 = linker_aligned[:L]
        s2 = display_origin[:L]

        n_match, total = render_alignment_row(ax, y_base + 4.5, s1, s2,
                                                x_start=0, char_w=char_w,
                                                fontsize=fontsize)
        ax.text(L * char_w + 2, y_base + 4.5,
                f"{n_match}/{total} matches ({100*n_match/total:.0f}%)",
                ha='left', va='center', fontsize=7, color='gray')

        # Row 3 (y = y_base + 6): Origin in genomic context with flanks
        ax.text(-25, y_base + 6.5,
                f"Genomic context\n({p['origin_contig']})",
                ha='left', va='center', fontsize=7)

        # Render the full origin with flanks, but mark the BLAST hit region
        # For Eimeria, the contig is mostly non-telomeric; classify origin as flank then mark hit
        origin_classes = ['flank'] * len(origin_full)
        # Position of the hit in origin_with_flanks
        hit_start_in_full = p.get('origin_flank_left_bp', 100)
        hit_end_in_full = hit_start_in_full + (p['origin_end'] - p['origin_start'] + 1)

        # If the origin is INSIDE another telotron, color the hit region as TTTAGGG repeats
        # Otherwise the hit region is the linker source — color it as shared
        if p.get('origin_in_telotron'):
            # Re-classify the full sequence as TTTAGGG telotron
            origin_classes = classify_bases(origin_full, 'TTTAGGG')
        # Highlight the BLAST hit region as 'shared'
        highlights2 = [(hit_start_in_full, hit_end_in_full, 'shared')]

        render_row(ax, y_base + 6, origin_full, origin_classes,
                   x_start=0, char_width=char_w,
                   highlight_ranges=highlights2, fontsize=fontsize)

        # Mark contig boundaries if close
        if p['origin_dist_to_end'] < 100:
            # The origin is very close to a contig end — show with a vertical dashed line
            if p['origin_start'] < 100:
                # Hit is near start of contig
                edge_x = -p['origin_start']
                ax.plot([edge_x, edge_x], [y_base + 6, y_base + 7],
                        color='black', linewidth=2, linestyle=':')
                ax.text(edge_x, y_base + 7.5, "← contig start",
                        fontsize=7, color='black')

        # Spacer
        ax.axhline(y=y_base + PAIR_HEIGHT - 0.3, color='gray',
                    linewidth=0.3, alpha=0.4)

    plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def main():
    with open('real_telotrons/linker_origin_pairs.json') as f:
        pairs = json.load(f)
    print(f"Loaded {len(pairs)} linker→origin pairs")

    # Need contig lengths for distance computation
    from collections import defaultdict
    contig_lens = {}
    for acc in set(p['acc'] for p in pairs):
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen: continue
        cur, n = None, 0
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contig_lens[(acc, cur)] = n
                    cur = line[1:].strip().split()[0]
                    n = 0
                else:
                    n += len(line.strip())
            if cur: contig_lens[(acc, cur)] = n

    for p in pairs:
        p['src_contig_len'] = contig_lens.get((p['acc'], p['source_contig']), 0)
        # origin_contig_len already in p

    # Sort: prioritize cases where origin is in another telotron, then by pident
    pairs_sorted = sorted(pairs, key=lambda p: (
        not p.get('origin_in_telotron'),
        -p['pident']
    ))
    # Add global index to track them across batches
    for i, p in enumerate(pairs_sorted):
        p['global_idx'] = i + 1

    # Save in batches of 6
    batch_size = 6
    for i, batch_start in enumerate(range(0, len(pairs_sorted), batch_size)):
        batch = pairs_sorted[batch_start:batch_start + batch_size]
        out = f'real_telotrons/seqfig_linker_origin_batch{i+1}.png'
        title = f"Linker → origin alignments (batch {i+1}/{(len(pairs_sorted) + batch_size - 1)//batch_size}, pairs {batch_start+1}-{batch_start+len(batch)} of {len(pairs_sorted)})"
        print(f"  Rendering {out}...")
        render_pairs(batch, out, title)
        print(f"  Saved: {out}")


if __name__ == '__main__':
    main()
