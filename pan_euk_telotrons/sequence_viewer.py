#!/usr/bin/env python3
"""
Sequence viewer module: render DNA sequences as colored rows with annotations.

Used by mechanism-figure scripts to show actual sequence-level examples for:
  - Each orientation class
  - Shared-linker tandem clusters
  - Haptophyte vs Eimeria comparison
  - Linker origin examples
"""

import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties


# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------

def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

# Forward/reverse repeat sets per family
TTAGGG_FWD = _rotations('TTAGGG')
TTAGGG_REV = _rotations('CCCTAA')
TTTAGGG_FWD = _rotations('TTTAGGG')
TTTAGGG_REV = _rotations('CCCTAAA')
TTAGG_FWD = _rotations('TTAGG')
TTAGG_REV = _rotations('CCTAA')
TTTGGG_FWD = _rotations('TTTGGG')
TTTGGG_REV = _rotations('CCCAAA')
TTGGGG_FWD = _rotations('TTGGGG')
TTGGGG_REV = _rotations('CCCCAA')

ALL_FWD = TTAGGG_FWD | TTTAGGG_FWD | TTAGG_FWD | TTTGGG_FWD | TTGGGG_FWD
ALL_REV = TTAGGG_REV | TTTAGGG_REV | TTAGG_REV | TTTGGG_REV | TTGGGG_REV

# Background colors per base class
COLOR_BG = {
    'non_telo':     np.array([0.96, 0.96, 0.96, 1.0]),  # near-white
    'fwd_canonical': np.array([0.30, 0.50, 0.85, 1.0]),  # blue
    'rev_canonical': np.array([0.85, 0.30, 0.30, 1.0]),  # red
    'fwd_variant':  np.array([0.55, 0.72, 1.00, 1.0]),  # light blue
    'rev_variant':  np.array([1.00, 0.62, 0.62, 1.0]),  # light red
    'flank':        np.array([0.92, 0.92, 0.85, 1.0]),  # cream-gray
    'tsd':          np.array([1.00, 0.88, 0.0, 1.0]),    # yellow highlight
    'shared':       np.array([0.65, 1.00, 0.65, 1.0]),   # green highlight
    'splice':       np.array([0.0, 0.0, 0.0, 1.0]),      # black
}

# Text color on each background
COLOR_TXT = {
    'non_telo':      (0.10, 0.10, 0.10),
    'fwd_canonical': (1.00, 1.00, 1.00),
    'rev_canonical': (1.00, 1.00, 1.00),
    'fwd_variant':   (0.0, 0.0, 0.30),
    'rev_variant':   (0.30, 0.0, 0.0),
    'flank':         (0.30, 0.30, 0.30),
    'tsd':           (0.0, 0.0, 0.0),
    'shared':        (0.0, 0.30, 0.0),
}


def classify_bases(seq, repeat_family='TTAGGG'):
    """Per-base classification: returns array of class names."""
    n = len(seq)
    classes = ['non_telo'] * n

    # Family-specific kmers
    if repeat_family == 'TTAGGG':
        canon = TTAGGG_FWD | TTAGGG_REV
        canon_fwd = TTAGGG_FWD
    elif repeat_family == 'TTTAGGG':
        canon = TTTAGGG_FWD | TTTAGGG_REV
        canon_fwd = TTTAGGG_FWD
    elif repeat_family == 'TTAGG':
        canon = TTAGG_FWD | TTAGG_REV
        canon_fwd = TTAGG_FWD
    elif repeat_family == 'TTTGGG':
        canon = TTTGGG_FWD | TTTGGG_REV
        canon_fwd = TTTGGG_FWD
    elif repeat_family == 'TTGGGG':
        canon = TTGGGG_FWD | TTGGGG_REV
        canon_fwd = TTGGGG_FWD
    else:
        canon = TTAGGG_FWD | TTAGGG_REV
        canon_fwd = TTAGGG_FWD

    canon_rev = canon - canon_fwd

    # Mark exact canonical first (highest priority)
    for k in canon_fwd:
        klen = len(k); i = 0
        while i < n - klen + 1:
            if seq[i:i+klen] == k:
                for j in range(i, i + klen):
                    classes[j] = 'fwd_canonical'
                i += klen  # don't double-mark, but allow next match
            else:
                i += 1
    for k in canon_rev:
        klen = len(k); i = 0
        while i < n - klen + 1:
            if seq[i:i+klen] == k:
                for j in range(i, i + klen):
                    if classes[j] == 'non_telo':
                        classes[j] = 'rev_canonical'
                i += klen
            else:
                i += 1

    # Mark 1-mismatch variants (only on remaining non_telo positions)
    rep_len = len(next(iter(canon_fwd)))
    for i in range(n - rep_len + 1):
        if classes[i] != 'non_telo':
            continue
        kmer = seq[i:i+rep_len]
        # Check 1-mismatch from any canonical fwd
        for c in canon_fwd:
            mm = sum(a != b for a, b in zip(kmer, c))
            if mm == 1:
                for j in range(i, i + rep_len):
                    if classes[j] == 'non_telo':
                        classes[j] = 'fwd_variant'
                break
        else:
            for c in canon_rev:
                mm = sum(a != b for a, b in zip(kmer, c))
                if mm == 1:
                    for j in range(i, i + rep_len):
                        if classes[j] == 'non_telo':
                            classes[j] = 'rev_variant'
                    break
    return classes


def render_row(ax, y, seq, classes, x_start=0, char_width=1.0,
               highlight_ranges=None, bg_override=None,
               text_color_override=None, show_text=True, fontsize=2.5):
    """
    Render a single sequence row.
      ax: matplotlib axis
      y: row index (0-based from top); each row is 1 unit tall
      seq: DNA sequence string
      classes: list of class names (per base)
      x_start: x-coordinate offset
      char_width: width of each character
      highlight_ranges: list of (start, end, type) where type is 'tsd', 'shared', etc.
      bg_override: list of (start, end, color) overriding the class-based bg
      text_color_override: dict {position: rgb}
    """
    n = len(seq)
    # Build per-position background
    bgs = [COLOR_BG[c] for c in classes]
    txts = [COLOR_TXT[c] for c in classes]

    # Apply overrides
    if highlight_ranges:
        for start, end, typ in highlight_ranges:
            color = COLOR_BG.get(typ, COLOR_BG['non_telo'])
            for j in range(max(0, start), min(n, end)):
                bgs[j] = color
                if typ in COLOR_TXT:
                    txts[j] = COLOR_TXT[typ]

    if bg_override:
        for start, end, color in bg_override:
            for j in range(max(0, start), min(n, end)):
                bgs[j] = np.array(color)

    if text_color_override:
        for pos, c in text_color_override.items():
            if 0 <= pos < n:
                txts[pos] = c

    # Draw background rectangles (per base)
    for j in range(n):
        x = x_start + j * char_width
        ax.add_patch(Rectangle((x, y), char_width, 1.0,
                                facecolor=bgs[j], edgecolor='none'))

    # Draw text
    if show_text:
        mono = FontProperties(family='monospace', size=fontsize)
        for j in range(n):
            x = x_start + j * char_width + char_width / 2
            ax.text(x, y + 0.5, seq[j], ha='center', va='center',
                    fontproperties=mono, color=txts[j], clip_on=True)


def add_splice_markers(ax, y, x_5p, x_3p, height=1.0):
    """Add black vertical bars at 5' and 3' splice sites."""
    bw = 0.4
    ax.add_patch(Rectangle((x_5p - bw, y), bw, height,
                            facecolor='black', edgecolor='none'))
    ax.add_patch(Rectangle((x_3p, y), bw, height,
                            facecolor='black', edgecolor='none'))


def add_label(ax, y, text, x=0, fontsize=8, color='black', ha='right'):
    """Add a row label to the left of the row."""
    ax.text(x, y + 0.5, text, ha=ha, va='center',
            fontsize=fontsize, color=color)


def render_full_intron_with_flanks(ax, y, intron_seq, left_flank, right_flank,
                                    repeat_family='TTAGGG',
                                    flank_visible_bp=60, char_width=1.0,
                                    highlight_ranges=None,
                                    label_text=None, label_x=-2,
                                    show_text=True, fontsize=2.5,
                                    x_start=0):
    """
    Render: [left_flank] | [splice_5p] [intron_seq] [splice_3p] | [right_flank]
    Returns total width.
    """
    # Truncate flanks
    lf = left_flank[-flank_visible_bp:] if left_flank else ''
    rf = right_flank[:flank_visible_bp] if right_flank else ''

    # Render left flank with 'flank' class for all
    if lf:
        lf_classes = ['flank'] * len(lf)
        render_row(ax, y, lf, lf_classes, x_start=x_start,
                   char_width=char_width, show_text=show_text, fontsize=fontsize)
        x_intron_start = x_start + len(lf) * char_width
    else:
        x_intron_start = x_start

    # Splice 5p marker
    splice_5p = x_intron_start
    boundary_w = char_width * 0.5

    # Intron
    intron_classes = classify_bases(intron_seq, repeat_family)
    # Adjust highlights to intron-relative coordinates
    intron_highlights = []
    if highlight_ranges:
        for hs, he, typ in highlight_ranges:
            intron_highlights.append((hs, he, typ))
    render_row(ax, y, intron_seq, intron_classes,
               x_start=x_intron_start + boundary_w,
               char_width=char_width,
               highlight_ranges=intron_highlights,
               show_text=show_text, fontsize=fontsize)
    intron_w = len(intron_seq) * char_width
    x_intron_end = x_intron_start + boundary_w + intron_w

    # Splice 3p
    splice_3p = x_intron_end

    # Splice site bars
    ax.add_patch(Rectangle((splice_5p, y), boundary_w, 1.0,
                            facecolor='black', edgecolor='none'))
    ax.add_patch(Rectangle((splice_3p, y), boundary_w, 1.0,
                            facecolor='black', edgecolor='none'))

    # Right flank
    x_rf_start = x_intron_end + boundary_w
    if rf:
        rf_classes = ['flank'] * len(rf)
        render_row(ax, y, rf, rf_classes, x_start=x_rf_start,
                   char_width=char_width, show_text=show_text, fontsize=fontsize)
        total_w = x_rf_start + len(rf) * char_width
    else:
        total_w = x_rf_start

    # Label
    if label_text:
        ax.text(label_x, y + 0.5, label_text, ha='right', va='center',
                fontsize=8)

    return total_w


def add_legend(ax, x, y, fontsize=8):
    """Add a small legend showing what each color means."""
    legend_items = [
        ('non_telo',       'Non-repeat'),
        ('fwd_canonical',  'TTAGGG/TTTAGGG (G-strand)'),
        ('rev_canonical',  'CCCTAA/CCCTAAA (C-strand)'),
        ('fwd_variant',    '1-mm variant fwd'),
        ('rev_variant',    '1-mm variant rev'),
        ('flank',          'Exon flank'),
        ('tsd',            'TSD'),
        ('shared',         'Shared with adjacent telotron'),
        ('splice',         'Splice site (GT/AG)'),
    ]
    for i, (cls, label) in enumerate(legend_items):
        ax.add_patch(Rectangle((x, y - i * 1.2), 1.5, 1.0,
                                facecolor=COLOR_BG[cls], edgecolor='black', linewidth=0.3))
        ax.text(x + 2, y - i * 1.2 + 0.5, label, ha='left', va='center', fontsize=fontsize)
