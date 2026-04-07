#!/usr/bin/env python3
"""
Lightweight analysis of telotron terminal sequences.
"""

import pandas as pd
import numpy as np
import json
import re
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

BASE_PATH = Path('/sessions/epic-peaceful-bohr/mnt/telotrons/tara_oceans_euk_mags')
INTRON_CATALOG = BASE_PATH / 'intron_candidates_high_confidence.tsv'
GENOME_DIR = BASE_PATH / 'smags/contigs_individual'
TELO_LABELS_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_labels.json'
TELO_PURITY_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_purity.npy'
OUTPUT_FIG = BASE_PATH / 'ED_Fig14_terminal_sequences.png'

TOP_8_MAGS = {
    'TARA_PSW_86_MAG_00284', 'TARA_PON_109_MAG_00250', 'TARA_ARC_108_MAG_00319',
    'TARA_PSE_93_MAG_00226', 'TARA_AOS_82_MAG_00154', 'TARA_AON_82_MAG_00318',
    'TARA_MED_95_MAG_00464', 'TARA_AOS_82_MAG_00183'
}

TELOMERIC_HEXAMERS = {
    'TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG',
    'CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA'
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def find_telomeric_array(seq):
    """Find telomeric hexamer runs."""
    runs = []
    i = 0
    while i < len(seq) - 5:
        if seq[i:i+6] in TELOMERIC_HEXAMERS:
            start = i
            j = i
            num_hex = 0
            while j <= len(seq) - 6:
                if seq[j:j+6] in TELOMERIC_HEXAMERS:
                    num_hex += 1
                    j += 6
                else:
                    break
            purity = num_hex / (num_hex + 1) if num_hex > 0 else 0
            runs.append((start, j, num_hex, purity))
            i = j
        else:
            i += 1
    return runs

def get_nuc_comp(seq):
    """Get nucleotide frequencies."""
    if not seq:
        return {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    total = len(seq)
    c = Counter(seq.upper())
    return {nt: c.get(nt, 0) / total for nt in 'ATGC'}

def extract_terminals(intron_seq):
    """Extract leader and trailer sequences."""
    arrays = find_telomeric_array(intron_seq)
    if not arrays:
        return None, None, None, None, None

    main_array = max(arrays, key=lambda x: x[2])
    telo_start, telo_end, num_hex, purity = main_array

    leader = intron_seq[:telo_start] if telo_start > 0 else None
    trailer = intron_seq[telo_end:] if telo_end < len(intron_seq) else None

    return leader, trailer, telo_start, telo_end, purity

def map_splice_sites_in_ttaggg():
    """Map GT and AG positions in TTAGGG."""
    seq = 'TTAGGG'
    sites = {'GT': [], 'AG': []}
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc == 'GT':
            sites['GT'].append(i)
        elif dinuc == 'AG':
            sites['AG'].append(i)
    return sites

def read_fasta_slice(fasta_file, contig_id, start, end):
    """Read a slice from a FASTA file efficiently."""
    try:
        with open(fasta_file, 'r') as f:
            in_target = False
            seq_start = 0
            seq_pos = 0
            current_seq = []

            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    header = line[1:].split()[0]
                    if header == contig_id:
                        in_target = True
                        seq_start = 0
                        seq_pos = 0
                    elif in_target:
                        break
                elif in_target:
                    for char in line:
                        if seq_pos >= start and seq_pos < end:
                            current_seq.append(char)
                        seq_pos += 1
                        if seq_pos >= end:
                            return ''.join(current_seq)

            return ''.join(current_seq) if current_seq else None
    except:
        return None

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

def main():
    print("=" * 80)
    print("TELOTRON TERMINAL SEQUENCE ANALYSIS")
    print("=" * 80)

    # Load telomeric metadata
    print("\n[1/5] Loading telomeric intron metadata...")
    with open(TELO_LABELS_PATH, 'r') as f:
        labels = json.load(f)
    purity_array = np.load(TELO_PURITY_PATH)

    telo_mags = labels['mag_ids']
    print(f"  Total telomeric introns: {len(telo_mags)}")

    # Filter to top 8 MAGs
    telo_top8_indices = [i for i, mag in enumerate(telo_mags) if mag in TOP_8_MAGS]
    print(f"  Telomeric introns in top 8 MAGs: {len(telo_top8_indices)}")

    if len(telo_top8_indices) == 0:
        print("  WARNING: No telomeric introns found in top 8 MAGs!")
        return

    # Load intron catalog
    print("\n[2/5] Loading intron catalog...")
    df = pd.read_csv(INTRON_CATALOG, sep='\t', low_memory=False)
    df_top8 = df[df['mag_id'].isin(TOP_8_MAGS)].copy().reset_index(drop=True)
    print(f"  Total introns in top 8 MAGs: {len(df_top8)}")

    # Extract terminal sequences from matched introns
    print("\n[3/5] Extracting terminal sequences...")

    analyses = []
    genome_cache = {}

    # Build a mapping of which MAGs have telomeric introns
    telo_mag_count = Counter(telo_mags[i] for i in telo_top8_indices)

    for mag_id in telo_mag_count.keys():
        mag_introns = df_top8[df_top8['mag_id'] == mag_id]
        if len(mag_introns) == 0:
            continue

        genome_file = GENOME_DIR / f"{mag_id}.fa"
        if not genome_file.exists():
            print(f"  Genome file not found: {mag_id}")
            continue

        print(f"  Processing {mag_id} ({len(mag_introns)} introns)...")

        # Process introns - limit to first 100 per MAG to avoid memory issues
        for idx, (_, row) in enumerate(mag_introns.iterrows()):
            if idx >= 100:  # Limit
                break

            contig_id = row['contig']
            start = int(row['bed_start'])
            end = int(row['bed_end'])

            # Read sequence from FASTA
            intron_seq = read_fasta_slice(str(genome_file), contig_id, start, end)
            if not intron_seq or len(intron_seq) < 12:
                continue

            # Extract terminals
            leader, trailer, telo_start, telo_end, purity = extract_terminals(intron_seq)

            if leader is None:  # No telomeric array found
                continue

            # Store results
            result = {
                'mag_id': mag_id,
                'contig': contig_id,
                'intron_len': len(intron_seq),
                'leader': leader,
                'leader_len': len(leader) if leader else 0,
                'trailer': trailer,
                'trailer_len': len(trailer) if trailer else 0,
                'telo_start': telo_start,
                'telo_end': telo_end,
                'telo_purity': purity,
                'donor': row['donor'],
                'acceptor': row['acceptor'],
                'splice_type': row['splice_type'],
            }
            analyses.append(result)

    df_analysis = pd.DataFrame(analyses)
    print(f"  Extracted {len(df_analysis)} introns with terminal sequences")

    # ==============================================================================
    # ANALYSIS
    # ==============================================================================
    print("\n[4/5] Computing statistics...")

    leaders_data = []
    trailers_data = []

    for _, row in df_analysis.iterrows():
        if row['leader'] and len(row['leader']) > 0:
            leaders_data.append({
                'length': len(row['leader']),
                'nuc_freq': get_nuc_comp(row['leader']),
            })

        if row['trailer'] and len(row['trailer']) > 0:
            trailers_data.append({
                'length': len(row['trailer']),
                'nuc_freq': get_nuc_comp(row['trailer']),
            })

    df_leaders = pd.DataFrame(leaders_data)
    df_trailers = pd.DataFrame(trailers_data)

    print(f"\n5' LEADER STATISTICS:")
    if len(df_leaders) > 0:
        print(f"  Count: {len(df_leaders)}")
        print(f"  Length - Mean: {df_leaders['length'].mean():.1f}bp ± {df_leaders['length'].std():.1f}bp")
        print(f"  Length - Median: {df_leaders['length'].median():.0f}bp")
        print(f"  Length - Range: {df_leaders['length'].min()}-{df_leaders['length'].max()}bp")

    print(f"\n3' TRAILER STATISTICS:")
    if len(df_trailers) > 0:
        print(f"  Count: {len(df_trailers)}")
        print(f"  Length - Mean: {df_trailers['length'].mean():.1f}bp ± {df_trailers['length'].std():.1f}bp")
        print(f"  Length - Median: {df_trailers['length'].median():.0f}bp")
        print(f"  Length - Range: {df_trailers['length'].min()}-{df_trailers['length'].max()}bp")

    splice_map = map_splice_sites_in_ttaggg()
    print(f"\nTTAGGG SPLICE SITE MAPPING:")
    print(f"  Sequence: T(0)-T(1)-A(2)-G(3)-G(4)-G(5)")
    print(f"  GT (donor) at position: {splice_map['GT']}")
    print(f"  AG (acceptor) at position: {splice_map['AG']}")

    # ==============================================================================
    # FIGURE
    # ==============================================================================
    print("\n[5/5] Creating figure...")

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # Panel A: 5' Leader length
    ax_a = fig.add_subplot(gs[0, 0])
    if len(df_leaders) > 0:
        ax_a.hist(df_leaders['length'], bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax_a.axvline(df_leaders['length'].mean(), color='red', linestyle='--', linewidth=2,
                     label=f"Mean: {df_leaders['length'].mean():.0f}bp")
        ax_a.set_xlabel("Leader length (bp)", fontsize=11, fontweight='bold')
        ax_a.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_a.set_title("A. 5' Leader Length Distribution", fontsize=12, fontweight='bold')
        ax_a.legend(fontsize=9)
        ax_a.grid(axis='y', alpha=0.3)

    # Panel B: 3' Trailer length
    ax_b = fig.add_subplot(gs[0, 1])
    if len(df_trailers) > 0:
        ax_b.hist(df_trailers['length'], bins=20, color='coral', edgecolor='black', alpha=0.7)
        ax_b.axvline(df_trailers['length'].mean(), color='red', linestyle='--', linewidth=2,
                     label=f"Mean: {df_trailers['length'].mean():.0f}bp")
        ax_b.set_xlabel("Trailer length (bp)", fontsize=11, fontweight='bold')
        ax_b.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_b.set_title("B. 3' Trailer Length Distribution", fontsize=12, fontweight='bold')
        ax_b.legend(fontsize=9)
        ax_b.grid(axis='y', alpha=0.3)

    # Panel C: Nucleotide composition
    ax_c = fig.add_subplot(gs[1, 0])
    nts = ['A', 'T', 'G', 'C']
    x_pos = np.arange(len(nts))
    width = 0.25

    if len(df_leaders) > 0:
        leader_comp = [np.mean([f[nt] for f in df_leaders['nuc_freq']]) for nt in nts]
        ax_c.bar(x_pos - width, leader_comp, width, label="Leaders", color='steelblue', edgecolor='black')

    if len(df_trailers) > 0:
        trailer_comp = [np.mean([f[nt] for f in df_trailers['nuc_freq']]) for nt in nts]
        ax_c.bar(x_pos + width, trailer_comp, width, label="Trailers", color='coral', edgecolor='black')

    ax_c.set_xlabel("Nucleotide", fontsize=11, fontweight='bold')
    ax_c.set_ylabel("Frequency", fontsize=11, fontweight='bold')
    ax_c.set_title("C. Nucleotide Composition", fontsize=12, fontweight='bold')
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(nts)
    ax_c.legend(fontsize=10)
    ax_c.grid(axis='y', alpha=0.3)

    # Panel D: TTAGGG splice sites
    ax_d = fig.add_subplot(gs[1, 1])
    seq = "TTAGGG"
    for i, nt in enumerate(seq):
        ax_d.add_patch(Rectangle((i, 0.3), 0.8, 0.4, facecolor='lightgray', edgecolor='black', linewidth=2))
        ax_d.text(i + 0.4, 0.5, nt, ha='center', va='center', fontsize=14, fontweight='bold')

    for pos in splice_map['GT']:
        ax_d.plot([pos + 0.4, pos + 1.2], [0.15, 0.15], 'g-', linewidth=3)
        ax_d.text(pos + 0.8, 0.05, 'GT', ha='center', fontsize=9, color='green', fontweight='bold')

    for pos in splice_map['AG']:
        ax_d.plot([pos + 0.4, pos + 1.2], [0.85, 0.85], 'r-', linewidth=3)
        ax_d.text(pos + 0.8, 0.95, 'AG', ha='center', fontsize=9, color='red', fontweight='bold')

    ax_d.set_xlim(-0.5, 6.5)
    ax_d.set_ylim(0, 1)
    ax_d.set_aspect('equal')
    ax_d.axis('off')
    ax_d.set_title("D. Splice Sites in TTAGGG (positions 0-5)", fontsize=12, fontweight='bold')
    ax_d.legend(handles=[Line2D([0], [0], color='g', lw=3, label='GT'),
                         Line2D([0], [0], color='r', lw=3, label='AG')],
                loc='upper left', fontsize=9)

    # Panel E: Placeholder
    ax_e = fig.add_subplot(gs[2, 0])
    ax_e.text(0.5, 0.5, f"E. Boundary Sequences\n(n={len(df_leaders)} 5' boundaries)",
             ha='center', va='center', fontsize=11, transform=ax_e.transAxes,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_e.axis('off')

    # Panel F: Scatter plot
    ax_f = fig.add_subplot(gs[2, 1])
    if len(df_analysis) > 1:
        valid_df = df_analysis[(df_analysis['leader_len'] > 0) & (df_analysis['trailer_len'] > 0) &
                               (df_analysis['telo_purity'] > 0)]
        if len(valid_df) > 1:
            ax_f.scatter(valid_df['leader_len'], valid_df['telo_purity'],
                        c=valid_df['trailer_len'], cmap='viridis', s=100, alpha=0.6, edgecolor='black')
            ax_f.set_xlabel("Leader Length (bp)", fontsize=11, fontweight='bold')
            ax_f.set_ylabel("Telomeric Purity", fontsize=11, fontweight='bold')
            ax_f.set_title("F. Terminal Length vs Purity", fontsize=12, fontweight='bold')
            ax_f.grid(alpha=0.3)

    plt.suptitle("Extended Data Figure 14: Telotron Terminal Sequences",
                 fontsize=14, fontweight='bold', y=0.995)

    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {OUTPUT_FIG}")

    # ==============================================================================
    # SUMMARY
    # ==============================================================================
    print("\n" + "=" * 80)
    print("SUMMARY OF KEY FINDINGS")
    print("=" * 80)

    print(f"\nDataset Summary:")
    print(f"  Telomeric introns analyzed: {len(df_analysis)}")
    print(f"  From {df_analysis['mag_id'].nunique()} unique MAGs")
    print(f"  Introns with identifiable 5' leaders: {len(df_leaders)}")
    print(f"  Introns with identifiable 3' trailers: {len(df_trailers)}")

    if len(df_leaders) > 0:
        print(f"\n5' Leader Terminal Sequences:")
        print(f"  Length distribution:")
        print(f"    Mean ± SD: {df_leaders['length'].mean():.1f} ± {df_leaders['length'].std():.1f} bp")
        print(f"    Median: {df_leaders['length'].median():.0f} bp")
        print(f"    Min-Max: {df_leaders['length'].min()}-{df_leaders['length'].max()} bp")
        leader_nuc = {nt: [] for nt in 'ATGC'}
        for f in df_leaders['nuc_freq']:
            for nt in 'ATGC':
                leader_nuc[nt].append(f[nt])
        gc = np.mean(leader_nuc['G']) + np.mean(leader_nuc['C'])
        print(f"  Nucleotide composition:")
        print(f"    GC content: {gc:.1%}")
        print(f"    A: {np.mean(leader_nuc['A']):.1%}, T: {np.mean(leader_nuc['T']):.1%}")

    if len(df_trailers) > 0:
        print(f"\n3' Trailer Terminal Sequences:")
        print(f"  Length distribution:")
        print(f"    Mean ± SD: {df_trailers['length'].mean():.1f} ± {df_trailers['length'].std():.1f} bp")
        print(f"    Median: {df_trailers['length'].median():.0f} bp")
        print(f"    Min-Max: {df_trailers['length'].min()}-{df_trailers['length'].max()} bp")
        trailer_nuc = {nt: [] for nt in 'ATGC'}
        for f in df_trailers['nuc_freq']:
            for nt in 'ATGC':
                trailer_nuc[nt].append(f[nt])
        gc = np.mean(trailer_nuc['G']) + np.mean(trailer_nuc['C'])
        print(f"  Nucleotide composition:")
        print(f"    GC content: {gc:.1%}")
        print(f"    A: {np.mean(trailer_nuc['A']):.1%}, T: {np.mean(trailer_nuc['T']):.1%}")

    print(f"\nTTAGGG Repeat and Splice Sites:")
    print(f"  Canonical sequence: 5'-TTAGGG-3'")
    print(f"  Positions:          0-1-2-3-4-5")
    print(f"  Donor site (GT): positions {splice_map['GT']}")
    print(f"  Acceptor site (AG): positions {splice_map['AG']}")
    print(f"  ")
    print(f"  KEY FINDING: AG dinucleotide appears at positions 2-3 within TTAGGG")
    print(f"  This means acceptor splice sites (AG) can occur within telomeric repeats,")
    print(f"  enabling cryptic splice sites within the telomeric array itself.")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)

if __name__ == '__main__':
    main()
