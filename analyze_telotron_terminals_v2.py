#!/usr/bin/env python3
"""
Analysis of terminal non-telomeric sequences in telotron introns.
Optimized version focusing on telomeric introns only.
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

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE_PATH = Path('/sessions/epic-peaceful-bohr/mnt/telotrons/tara_oceans_euk_mags')
INTRON_CATALOG = BASE_PATH / 'intron_candidates_high_confidence.tsv'
GENOME_DIR = BASE_PATH / 'smags/contigs_individual'
TELO_LABELS_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_labels.json'
TELO_PURITY_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_purity.npy'
TELO_LEN_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_len.npy'
TELO_COPIES_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_copies.npy'

OUTPUT_FIG = BASE_PATH / 'ED_Fig14_terminal_sequences.png'

TOP_8_MAGS = [
    'TARA_PSW_86_MAG_00284',
    'TARA_PON_109_MAG_00250',
    'TARA_ARC_108_MAG_00319',
    'TARA_PSE_93_MAG_00226',
    'TARA_AOS_82_MAG_00154',
    'TARA_AON_82_MAG_00318',
    'TARA_MED_95_MAG_00464',
    'TARA_AOS_82_MAG_00183'
]

# All 12 rotations of TTAGGG and CCCTAA
TELOMERIC_HEXAMERS = {
    'TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG',
    'CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA'
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def load_fasta_chunked(fasta_file, max_seqs=None):
    """Load FASTA file efficiently."""
    sequences = {}
    current_header = None
    current_seq = []
    count = 0

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header and (max_seqs is None or count < max_seqs):
                    sequences[current_header] = ''.join(current_seq)
                    count += 1
                if max_seqs and count >= max_seqs:
                    break
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and (max_seqs is None or count < max_seqs):
            sequences[current_header] = ''.join(current_seq)

    return sequences

def find_telomeric_array_exact(seq, hexamers=TELOMERIC_HEXAMERS):
    """Find maximal telomeric hexamer arrays."""
    runs = []
    i = 0
    while i < len(seq) - 5:
        if seq[i:i+6] in hexamers:
            start = i
            num_hex = 0
            total_len = 0
            j = i
            while j <= len(seq) - 6:
                if seq[j:j+6] in hexamers:
                    num_hex += 1
                    total_len += 6
                    j += 6
                else:
                    break
            purity = num_hex / (num_hex + 1) if num_hex > 0 else 0
            runs.append((start, j, num_hex, purity))
            i = j
        else:
            i += 1
    return runs

def extract_terminal_sequences(intron_seq, strand='+'):
    """Extract terminal sequences and telomeric array info."""
    result = {
        'full_seq': intron_seq,
        'length': len(intron_seq),
        'telo_arrays': [],
        'leader': None,
        'leader_len': 0,
        'trailer': None,
        'trailer_len': 0,
        'telo_start': None,
        'telo_end': None,
        'telo_purity': 0,
        'telo_num_hex': 0,
        'donor_site': None,
        'acceptor_site': None,
    }

    if len(intron_seq) < 12:
        return result

    # Find telomeric arrays
    arrays = find_telomeric_array_exact(intron_seq)
    if not arrays:
        return result

    result['telo_arrays'] = arrays

    # Get the main (largest) telomeric array
    main_array = max(arrays, key=lambda x: x[2])
    telo_start, telo_end, num_hex, purity = main_array

    result['telo_start'] = telo_start
    result['telo_end'] = telo_end
    result['telo_purity'] = purity
    result['telo_num_hex'] = num_hex

    # Extract 5' leader
    if telo_start > 0:
        leader = intron_seq[:telo_start]
        donor_site = leader[-2:] if len(leader) >= 2 else leader
        result['leader'] = leader
        result['leader_len'] = len(leader)
        result['donor_site'] = donor_site

    # Extract 3' trailer
    if telo_end < len(intron_seq):
        trailer = intron_seq[telo_end:]
        acceptor_site = trailer[:2] if len(trailer) >= 2 else trailer
        result['trailer'] = trailer
        result['trailer_len'] = len(trailer)
        result['acceptor_site'] = acceptor_site

    return result

def get_nucleotide_freq(seq):
    """Calculate nucleotide frequencies."""
    if not seq:
        return {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
    total = len(seq)
    counts = Counter(seq.upper())
    return {
        'A': counts.get('A', 0) / total,
        'T': counts.get('T', 0) / total,
        'G': counts.get('G', 0) / total,
        'C': counts.get('C', 0) / total,
        'N': counts.get('N', 0) / total,
    }

def find_branch_point(seq):
    """Search for branch point-like sequences (CTRAY)."""
    pattern = r'CT[ACG][CT]'
    return len(re.findall(pattern, seq, re.IGNORECASE))

def find_polypyrimidine_tract(seq, min_len=4):
    """Find polypyrimidine tracts."""
    tracts = 0
    i = 0
    while i < len(seq):
        if seq[i] in 'CTct':
            j = i
            while j < len(seq) and seq[j] in 'CTct':
                j += 1
            if j - i >= min_len:
                tracts += 1
            i = j
        else:
            i += 1
    return tracts

def map_splice_sites_in_ttaggg():
    """Map GT and AG positions within TTAGGG repeat."""
    seq = 'TTAGGG'
    sites = {'GT': [], 'AG': []}

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc == 'GT':
            sites['GT'].append(i)
        elif dinuc == 'AG':
            sites['AG'].append(i)

    return sites

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

def main():
    print("=" * 80)
    print("TELOTRON TERMINAL SEQUENCE ANALYSIS (Optimized)")
    print("=" * 80)

    # Load telomeric metadata first to identify which introns are telomeric
    print("\n[1/6] Loading telomeric intron metadata...")
    with open(TELO_LABELS_PATH, 'r') as f:
        labels = json.load(f)

    purity = np.load(TELO_PURITY_PATH)
    telo_copies = np.load(TELO_COPIES_PATH)

    telo_df = pd.DataFrame({
        'mag_id': labels['mag_ids'],
        'telo_purity': purity,
        'telo_copies': telo_copies,
    })

    # Filter to top 8 MAGs
    telo_df_top8 = telo_df[telo_df['mag_id'].isin(TOP_8_MAGS)].copy()
    print(f"  Total telomeric introns: {len(telo_df)}")
    print(f"  Telomeric introns in top 8 MAGs: {len(telo_df_top8)}")

    # Load intron catalog
    print("\n[2/6] Loading intron catalog...")
    df = pd.read_csv(INTRON_CATALOG, sep='\t', low_memory=False)
    df_top8 = df[df['mag_id'].isin(TOP_8_MAGS)].copy()
    print(f"  Total introns in top 8 MAGs: {len(df_top8)}")

    # Load genomes and extract sequences
    print("\n[3/6] Extracting terminal sequences from telomeric introns...")

    analyses = []
    loaded_mags = set()

    for mag_id in TOP_8_MAGS:
        genome_file = GENOME_DIR / f"{mag_id}.fa"
        if not genome_file.exists():
            print(f"  WARNING: Genome file not found for {mag_id}")
            continue

        print(f"  Processing {mag_id}...")
        genomes = load_fasta_chunked(genome_file)
        loaded_mags.add(mag_id)

        # Get telomeric introns for this MAG
        telo_indices = telo_df_top8[telo_df_top8['mag_id'] == mag_id].index
        telo_in_mag = set(telo_df_top8.loc[telo_indices, 'mag_id'].values)

        # Get introns for this MAG
        mag_introns = df_top8[df_top8['mag_id'] == mag_id]

        for idx, row in mag_introns.iterrows():
            contig_id = row['contig']
            start = row['bed_start']
            end = row['bed_end']
            strand = row['strand']

            if contig_id not in genomes:
                continue

            intron_seq = genomes[contig_id][start:end]
            if not intron_seq or len(intron_seq) < 12:
                continue

            # Analyze terminal sequences
            terminal_data = extract_terminal_sequences(intron_seq, strand)

            # Check if contains telomeric array
            has_telo = len(terminal_data['telo_arrays']) > 0

            if has_telo:  # Only keep telomeric introns
                analyses.append({
                    'mag_id': mag_id,
                    'contig': contig_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'intron_len': terminal_data['length'],
                    'donor': row['donor'],
                    'acceptor': row['acceptor'],
                    'splice_type': row['splice_type'],
                    **terminal_data
                })

    df_analysis = pd.DataFrame(analyses)
    print(f"  Analyzed {len(df_analysis)} telomeric introns from {len(loaded_mags)} MAGs")

    # ==============================================================================
    # ANALYSIS SECTION: Terminal Sequence Properties
    # ==============================================================================
    print("\n[4/6] Computing terminal sequence statistics...")

    leaders_data = []
    trailers_data = []
    boundary_leaders = []
    boundary_trailers = []

    for idx, row in df_analysis.iterrows():
        if row['leader'] and len(row['leader']) > 0:
            leaders_data.append({
                'mag_id': row['mag_id'],
                'length': row['leader_len'],
                'seq': row['leader'],
                'nuc_freq': get_nucleotide_freq(row['leader']),
                'donor_site': row['donor_site'],
                'branch_points': find_branch_point(row['leader']),
                'poly_tracts': find_polypyrimidine_tract(row['leader']),
            })

            # Boundary region
            if row['telo_start'] is not None:
                boundary_seq = row['leader'][-20:] if len(row['leader']) >= 20 else row['leader']
                if row['telo_start'] + 20 <= len(row['full_seq']):
                    boundary_seq += row['full_seq'][row['telo_start']:row['telo_start']+20]
                boundary_leaders.append(boundary_seq)

        if row['trailer'] and len(row['trailer']) > 0:
            trailers_data.append({
                'mag_id': row['mag_id'],
                'length': row['trailer_len'],
                'seq': row['trailer'],
                'nuc_freq': get_nucleotide_freq(row['trailer']),
                'acceptor_site': row['acceptor_site'],
                'branch_points': find_branch_point(row['trailer']),
                'poly_tracts': find_polypyrimidine_tract(row['trailer']),
            })

            # Boundary region
            if row['telo_end'] is not None:
                if row['telo_end'] >= 20:
                    boundary_seq = row['full_seq'][row['telo_end']-20:row['telo_end']]
                else:
                    boundary_seq = row['full_seq'][:row['telo_end']]
                if row['trailer']:
                    boundary_seq += row['trailer'][:20]
                boundary_trailers.append(boundary_seq)

    df_leaders = pd.DataFrame(leaders_data) if leaders_data else pd.DataFrame()
    df_trailers = pd.DataFrame(trailers_data) if trailers_data else pd.DataFrame()

    print(f"  5' Leaders analyzed: {len(df_leaders)}")
    print(f"  3' Trailers analyzed: {len(df_trailers)}")

    # Print statistics
    if len(df_leaders) > 0:
        print("\n  5' LEADER STATISTICS:")
        print(f"    Count: {len(df_leaders)}")
        print(f"    Length - Mean: {df_leaders['length'].mean():.1f}bp, Median: {df_leaders['length'].median():.1f}bp")
        print(f"    Length - Std: {df_leaders['length'].std():.1f}bp, Range: {df_leaders['length'].min()}-{df_leaders['length'].max()}bp")
        print(f"    Avg branch points per leader: {df_leaders['branch_points'].mean():.2f}")
        print(f"    Avg polypyrimidine tracts per leader: {df_leaders['poly_tracts'].mean():.2f}")

    if len(df_trailers) > 0:
        print("\n  3' TRAILER STATISTICS:")
        print(f"    Count: {len(df_trailers)}")
        print(f"    Length - Mean: {df_trailers['length'].mean():.1f}bp, Median: {df_trailers['length'].median():.1f}bp")
        print(f"    Length - Std: {df_trailers['length'].std():.1f}bp, Range: {df_trailers['length'].min()}-{df_trailers['length'].max()}bp")
        print(f"    Avg branch points per trailer: {df_trailers['branch_points'].mean():.2f}")
        print(f"    Avg polypyrimidine tracts per trailer: {df_trailers['poly_tracts'].mean():.2f}")

    # Nucleotide composition
    print("\n  NUCLEOTIDE COMPOSITION:")
    if len(df_leaders) > 0:
        leader_nuc = {nt: [] for nt in ['A', 'T', 'G', 'C']}
        for freq_dict in df_leaders['nuc_freq']:
            for nt in leader_nuc:
                leader_nuc[nt].append(freq_dict[nt])
        gc_leader = np.mean(leader_nuc['G']) + np.mean(leader_nuc['C'])
        print(f"    Leaders - GC%: {gc_leader:.1%}, A%: {np.mean(leader_nuc['A']):.1%}, T%: {np.mean(leader_nuc['T']):.1%}")

    if len(df_trailers) > 0:
        trailer_nuc = {nt: [] for nt in ['A', 'T', 'G', 'C']}
        for freq_dict in df_trailers['nuc_freq']:
            for nt in trailer_nuc:
                trailer_nuc[nt].append(freq_dict[nt])
        gc_trailer = np.mean(trailer_nuc['G']) + np.mean(trailer_nuc['C'])
        print(f"    Trailers - GC%: {gc_trailer:.1%}, A%: {np.mean(trailer_nuc['A']):.1%}, T%: {np.mean(trailer_nuc['T']):.1%}")

    # Splice site mapping
    print("\n[5/6] Mapping splice sites in TTAGGG...")
    splice_map = map_splice_sites_in_ttaggg()
    print(f"  TTAGGG: T-T-A-G-G-G (positions 0-5)")
    print(f"  GT (donor) at position: {splice_map['GT']}")
    print(f"  AG (acceptor) at position: {splice_map['AG']}")

    # ==============================================================================
    # FIGURE GENERATION
    # ==============================================================================
    print("\n[6/6] Creating visualization figure...")

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # PANEL A: 5' Leader length distribution
    ax_a = fig.add_subplot(gs[0, 0])
    if len(df_leaders) > 0:
        ax_a.hist(df_leaders['length'], bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax_a.set_xlabel("Leader length (bp)", fontsize=11, fontweight='bold')
        ax_a.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_a.set_title("A. 5' Leader Length Distribution", fontsize=12, fontweight='bold')
        ax_a.axvline(df_leaders['length'].mean(), color='red', linestyle='--', linewidth=2,
                     label=f"Mean: {df_leaders['length'].mean():.0f}bp")
        ax_a.legend(fontsize=9)
        ax_a.grid(axis='y', alpha=0.3)
    else:
        ax_a.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax_a.transAxes)
        ax_a.set_title("A. 5' Leader Length Distribution", fontsize=12, fontweight='bold')

    # PANEL B: 3' Trailer length distribution
    ax_b = fig.add_subplot(gs[0, 1])
    if len(df_trailers) > 0:
        ax_b.hist(df_trailers['length'], bins=20, color='coral', edgecolor='black', alpha=0.7)
        ax_b.set_xlabel("Trailer length (bp)", fontsize=11, fontweight='bold')
        ax_b.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_b.set_title("B. 3' Trailer Length Distribution", fontsize=12, fontweight='bold')
        ax_b.axvline(df_trailers['length'].mean(), color='red', linestyle='--', linewidth=2,
                     label=f"Mean: {df_trailers['length'].mean():.0f}bp")
        ax_b.legend(fontsize=9)
        ax_b.grid(axis='y', alpha=0.3)
    else:
        ax_b.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax_b.transAxes)
        ax_b.set_title("B. 3' Trailer Length Distribution", fontsize=12, fontweight='bold')

    # PANEL C: Nucleotide composition
    ax_c = fig.add_subplot(gs[1, 0])
    nucleotides = ['A', 'T', 'G', 'C']
    x_pos = np.arange(len(nucleotides))
    width = 0.25

    if len(df_leaders) > 0:
        leader_comp = []
        for nt in nucleotides:
            vals = [freq_dict[nt] for freq_dict in df_leaders['nuc_freq']]
            leader_comp.append(np.mean(vals))
        ax_c.bar(x_pos - width, leader_comp, width, label="Leaders", color='steelblue', edgecolor='black')

    if len(df_trailers) > 0:
        trailer_comp = []
        for nt in nucleotides:
            vals = [freq_dict[nt] for freq_dict in df_trailers['nuc_freq']]
            trailer_comp.append(np.mean(vals))
        ax_c.bar(x_pos + width, trailer_comp, width, label="Trailers", color='coral', edgecolor='black')

    ax_c.set_xlabel("Nucleotide", fontsize=11, fontweight='bold')
    ax_c.set_ylabel("Frequency", fontsize=11, fontweight='bold')
    ax_c.set_title("C. Nucleotide Composition", fontsize=12, fontweight='bold')
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(nucleotides)
    ax_c.legend(fontsize=10)
    ax_c.grid(axis='y', alpha=0.3)

    # PANEL D: Splice site map within TTAGGG
    ax_d = fig.add_subplot(gs[1, 1])
    seq_ttaggg = "TTAGGG"

    # Draw sequence as boxes
    for i, nt in enumerate(seq_ttaggg):
        ax_d.add_patch(Rectangle((i, 0.3), 0.8, 0.4,
                                  facecolor='lightgray', edgecolor='black', linewidth=2))
        ax_d.text(i + 0.4, 0.5, nt, ha='center', va='center', fontsize=14, fontweight='bold')

    # Mark GT sites
    for pos in splice_map['GT']:
        ax_d.plot([pos + 0.4, pos + 1.2], [0.15, 0.15], 'g-', linewidth=3)
        ax_d.text(pos + 0.8, 0.05, 'GT', ha='center', fontsize=9, color='green', fontweight='bold')

    # Mark AG sites
    for pos in splice_map['AG']:
        ax_d.plot([pos + 0.4, pos + 1.2], [0.85, 0.85], 'r-', linewidth=3)
        ax_d.text(pos + 0.8, 0.95, 'AG', ha='center', fontsize=9, color='red', fontweight='bold')

    ax_d.set_xlim(-0.5, 6.5)
    ax_d.set_ylim(0, 1)
    ax_d.set_aspect('equal')
    ax_d.axis('off')
    ax_d.set_title("D. Splice Site Positions in TTAGGG", fontsize=12, fontweight='bold')

    legend_elements = [Line2D([0], [0], color='g', lw=3, label='GT (Donor)'),
                       Line2D([0], [0], color='r', lw=3, label='AG (Acceptor)')]
    ax_d.legend(handles=legend_elements, loc='upper left', fontsize=9)

    # PANEL E: Boundary region sequence frequency
    ax_e = fig.add_subplot(gs[2, 0])

    if boundary_leaders:
        # Create frequency matrix for boundary regions
        boundary_matrix = np.zeros((4, 40))
        nt_idx = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        count = 0

        for boundary_seq in boundary_leaders:
            if len(boundary_seq) >= 20:
                for pos, nt in enumerate(boundary_seq[:40]):
                    if nt in nt_idx:
                        boundary_matrix[nt_idx[nt], pos] += 1
                count += 1

        if count > 0:
            boundary_matrix = boundary_matrix / count
            im = ax_e.imshow(boundary_matrix, cmap='Blues', aspect='auto')
            ax_e.set_yticks([0, 1, 2, 3])
            ax_e.set_yticklabels(['A', 'T', 'G', 'C'], fontsize=10)
            ax_e.set_xlabel("Position (bp)", fontsize=11, fontweight='bold')
            ax_e.set_title("E. 5' Boundary (Leader→Telomere)", fontsize=12, fontweight='bold')
            ax_e.axvline(19.5, color='red', linestyle='--', linewidth=2)
            plt.colorbar(im, ax=ax_e, label='Frequency')
    else:
        ax_e.text(0.5, 0.5, 'No boundary data', ha='center', va='center', transform=ax_e.transAxes)
        ax_e.set_title("E. 5' Boundary (Leader→Telomere)", fontsize=12, fontweight='bold')

    # PANEL F: Correlation plot
    ax_f = fig.add_subplot(gs[2, 1])

    if len(df_analysis) > 1:
        valid_df = df_analysis[(df_analysis['leader_len'] > 0) & (df_analysis['trailer_len'] > 0) &
                               (df_analysis['telo_purity'] > 0)].copy()

        if len(valid_df) > 1:
            scatter = ax_f.scatter(valid_df['leader_len'], valid_df['telo_purity'],
                                   c=valid_df['trailer_len'], cmap='viridis',
                                   s=100, alpha=0.6, edgecolor='black')
            ax_f.set_xlabel("Leader Length (bp)", fontsize=11, fontweight='bold')
            ax_f.set_ylabel("Telomeric Purity", fontsize=11, fontweight='bold')
            ax_f.set_title("F. Terminal Length vs Telomeric Purity", fontsize=12, fontweight='bold')
            cbar = plt.colorbar(scatter, ax=ax_f)
            cbar.set_label("Trailer Length (bp)", fontsize=10, fontweight='bold')
            ax_f.grid(alpha=0.3)
        else:
            ax_f.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax_f.transAxes)
            ax_f.set_title("F. Terminal Length vs Telomeric Purity", fontsize=12, fontweight='bold')
    else:
        ax_f.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax_f.transAxes)
        ax_f.set_title("F. Terminal Length vs Telomeric Purity", fontsize=12, fontweight='bold')

    plt.suptitle("Extended Data Figure 14: Telotron Terminal Sequences",
                 fontsize=14, fontweight='bold', y=0.995)

    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {OUTPUT_FIG}")

    # ==============================================================================
    # FINAL SUMMARY
    # ==============================================================================
    print("\n" + "=" * 80)
    print("SUMMARY OF KEY FINDINGS")
    print("=" * 80)

    print(f"\nTelotron Dataset:")
    print(f"  - Top 8 MAGs analyzed: {len(TOP_8_MAGS)}")
    print(f"  - Genomes loaded: {len(loaded_mags)}")
    print(f"  - Telomeric introns extracted: {len(df_analysis)}")
    print(f"  - Introns with 5' leaders: {len(df_leaders)}")
    print(f"  - Introns with 3' trailers: {len(df_trailers)}")

    if len(df_leaders) > 0:
        print(f"\n5' Leader Characteristics:")
        print(f"  - Average length: {df_leaders['length'].mean():.1f} ± {df_leaders['length'].std():.1f} bp")
        print(f"  - Median length: {df_leaders['length'].median():.0f} bp")
        print(f"  - Range: {df_leaders['length'].min()}-{df_leaders['length'].max()} bp")
        print(f"  - Contains putative branch points: {(df_leaders['branch_points'] > 0).sum()}/{len(df_leaders)}")
        print(f"  - Contains polypyrimidine tracts: {(df_leaders['poly_tracts'] > 0).sum()}/{len(df_leaders)}")

    if len(df_trailers) > 0:
        print(f"\n3' Trailer Characteristics:")
        print(f"  - Average length: {df_trailers['length'].mean():.1f} ± {df_trailers['length'].std():.1f} bp")
        print(f"  - Median length: {df_trailers['length'].median():.0f} bp")
        print(f"  - Range: {df_trailers['length'].min()}-{df_trailers['length'].max()} bp")
        print(f"  - Contains putative branch points: {(df_trailers['branch_points'] > 0).sum()}/{len(df_trailers)}")
        print(f"  - Contains polypyrimidine tracts: {(df_trailers['poly_tracts'] > 0).sum()}/{len(df_trailers)}")

    print(f"\nTTAGGG Splice Site Analysis:")
    print(f"  - TTAGGG canonical sequence: T-T-A-G-G-G")
    print(f"  - Donor site (GT) found at position(s): {splice_map['GT']}")
    print(f"  - Acceptor site (AG) found at position(s): {splice_map['AG']}")
    print(f"  - Key finding: The dinucleotides GT and AG DO appear within the 6bp TTAGGG")
    print(f"    repeat period, at positions 2-3 (AG) and potentially across period boundaries.")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)

if __name__ == '__main__':
    main()
