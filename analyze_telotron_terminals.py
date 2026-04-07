#!/usr/bin/env python3
"""
Analysis of terminal non-telomeric sequences in telotron introns.
Focuses on 5' leaders and 3' trailers flanking telomeric arrays.
"""

import pandas as pd
import numpy as np
import json
import re
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE_PATH = Path('/sessions/epic-peaceful-bohr/mnt/telotrons/tara_oceans_euk_mags')
INTRON_CATALOG = BASE_PATH / 'intron_candidates_high_confidence.tsv'
GENOME_DIR = BASE_PATH / 'smags/contigs_individual'
GFF_DIR = BASE_PATH / 'smags/gff_individual/GFF'
TELO_LABELS_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_labels.json'
TELO_PURITY_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_purity.npy'
TELO_LEN_PATH = BASE_PATH / 'intron_clustering/tandem6bp_clusters/telomeric_deep/telo_len.npy'

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

# Splice site patterns (first 2bp of donor, last 2bp of acceptor)
DONOR_PATTERNS = {'GT', 'CT', 'GC', 'AT', 'TT'}  # More permissive
ACCEPTOR_PATTERNS = {'AG', 'AC', 'GG', 'AT'}     # More permissive

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def load_fasta(fasta_file):
    """Load FASTA file into dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences

def find_telomeric_run(seq, hexamers=TELOMERIC_HEXAMERS):
    """
    Find the longest run of telomeric hexamers in a sequence.
    Returns (start_pos, end_pos, length_in_hexamers, purity).
    """
    if len(seq) < 6:
        return None, None, None, 0

    best_run = (0, 0, 0, 0)  # (start, end, num_hex, purity)

    for i in range(len(seq) - 5):
        hexamer = seq[i:i+6]
        if hexamer in hexamers:
            # Start a potential run
            j = i
            matched = 0
            total = 0
            while j <= len(seq) - 6:
                hex_test = seq[j:j+6]
                total += 6
                if hex_test in hexamers:
                    matched += 6
                    j += 6
                elif j + 5 < len(seq):
                    # Try sliding by 1 (overlapping repeats)
                    hex_slide = seq[j+1:j+7]
                    if hex_slide in hexamers:
                        matched += 6
                        total += 1
                        j += 7
                    else:
                        break
                else:
                    break

            if matched > 0:
                purity = matched / total if total > 0 else 0
                if matched > best_run[2] * 6:
                    best_run = (i, j, matched // 6, purity)

    if best_run[2] > 0:
        return best_run
    return None, None, None, 0

def find_telomeric_array_exact(seq, hexamers=TELOMERIC_HEXAMERS):
    """
    Find maximal telomeric hexamer arrays by looking for consecutive 6bp matches.
    Returns list of (start, end, num_hexamers, purity).
    """
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
            purity = num_hex / (num_hex + 1) if num_hex > 0 else 0  # Account for boundary
            runs.append((start, j, num_hex, purity))
            i = j
        else:
            i += 1

    return runs

def extract_terminal_sequences(intron_seq, strand='+'):
    """
    Given full intron sequence, identify:
    - 5' leader: from splice site (GT/CT/GC) to first telomeric hexamer
    - 3' trailer: from last telomeric hexamer to splice site (AG/AC)
    Returns dict with leaders, trailers, telomeric positions.
    """
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

    # Extract 5' leader
    if telo_start > 0:
        leader = intron_seq[:telo_start]
        # Find last 2bp (donor splice site)
        donor_site = leader[-2:] if len(leader) >= 2 else leader
        result['leader'] = leader
        result['leader_len'] = len(leader)
        result['donor_site'] = donor_site

    # Extract 3' trailer
    if telo_end < len(intron_seq):
        trailer = intron_seq[telo_end:]
        # Find first 2bp (acceptor splice site)
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

def find_branch_point(seq, window=50):
    """
    Search for branch point-like sequences (CTRAY = CT[A/C/G]Y where Y=C/T).
    Returns list of positions.
    """
    pattern = r'CT[ACG][CT]'
    matches = []
    for match in re.finditer(pattern, seq, re.IGNORECASE):
        matches.append(match.start())
    return matches

def find_polypyrimidine_tract(seq, min_len=4):
    """
    Find polypyrimidine (C/T) tracts.
    Returns list of (start, end, length) tuples.
    """
    tracts = []
    i = 0
    while i < len(seq):
        if seq[i] in 'CTct':
            j = i
            while j < len(seq) and seq[j] in 'CTct':
                j += 1
            tract_len = j - i
            if tract_len >= min_len:
                tracts.append((i, j, tract_len))
            i = j
        else:
            i += 1
    return tracts

def map_splice_sites_in_ttaggg():
    """
    Map GT and AG positions within TTAGGG repeat.
    TTAGGG = T-T-A-G-G-G (positions 0-5)
    Returns dict with positions of GT and AG motifs.
    """
    seq = 'TTAGGG'
    sites = {
        'GT': [],
        'AG': [],
        'full_sequence': seq,
    }

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc == 'GT':
            sites['GT'].append(i)
        elif dinuc == 'AG':
            sites['AG'].append(i)

    # Also check with rotation context (spanning period boundary)
    extended = seq + seq[:2]  # TTAGGG + TT = TTAGGGTT
    for i in range(len(seq)):
        dinuc = extended[i:i+2]
        if dinuc == 'GT' and i not in sites['GT']:
            sites['GT'].append(('spanning', i, i+2 if i+2 <= 6 else 'wraps'))
        elif dinuc == 'AG' and i not in sites['AG']:
            sites['AG'].append(('spanning', i, i+2 if i+2 <= 6 else 'wraps'))

    return sites

def load_telomeric_metadata():
    """Load telomeric intron indices and purity values."""
    with open(TELO_LABELS_PATH, 'r') as f:
        labels = json.load(f)

    # Load purity and length arrays
    purity = np.load(TELO_PURITY_PATH)
    length = np.load(TELO_LEN_PATH)

    return {
        'mag_ids': labels['mag_ids'],
        'purity': purity,
        'length': length,
        'n_telotrons': len(labels['mag_ids'])
    }

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

def main():
    print("=" * 80)
    print("TELOTRON TERMINAL SEQUENCE ANALYSIS")
    print("=" * 80)

    # Load intron catalog
    print("\n[1/7] Loading intron catalog...")
    df = pd.read_csv(INTRON_CATALOG, sep='\t', low_memory=False)
    print(f"  Total introns: {len(df)}")

    # Filter to top 8 MAGs
    df_top8 = df[df['mag_id'].isin(TOP_8_MAGS)].copy()
    print(f"  Introns in top 8 MAGs: {len(df_top8)}")

    # Load telomeric metadata
    print("\n[2/7] Loading telomeric intron metadata...")
    metadata = load_telomeric_metadata()
    telo_mags = set(metadata['mag_ids'])
    print(f"  Total telomeric introns: {metadata['n_telotrons']}")

    # Identify telomeric introns from top 8 MAGs
    telo_in_top8 = [mid for mid in metadata['mag_ids'] if mid in TOP_8_MAGS]
    print(f"  Telomeric introns from top 8 MAGs: {len(telo_in_top8)}")

    # Extract sequences and analyze
    print("\n[3/7] Extracting intron sequences and identifying terminals...")

    analyses = []
    genome_cache = {}

    # For each MAG, load genome once
    for mag_id in TOP_8_MAGS:
        genome_file = GENOME_DIR / f"{mag_id}.fa"
        if genome_file.exists():
            genome_cache[mag_id] = load_fasta(genome_file)
            print(f"  Loaded {mag_id}: {len(genome_cache[mag_id])} contigs")
        else:
            print(f"  WARNING: Genome file not found for {mag_id}")

    # Get telomeric indices from top 8 MAGs
    telo_indices = [
        i for i, mid in enumerate(metadata['mag_ids'])
        if mid in TOP_8_MAGS
    ]

    # Iterate through top 8 MAGs and their introns
    processed_telos = 0
    for mag_id in TOP_8_MAGS:
        mag_introns = df_top8[df_top8['mag_id'] == mag_id]
        if mag_id not in genome_cache:
            continue

        genomes = genome_cache[mag_id]

        for idx, row in mag_introns.iterrows():
            contig_id = row['contig']
            start = row['bed_start']
            end = row['bed_end']
            strand = row['strand']

            if contig_id not in genomes:
                continue

            intron_seq = genomes[contig_id][start:end]

            # Check if this intron is telomeric
            is_telo = mag_id in telo_in_top8 and processed_telos < len(telo_in_top8)

            # Analyze terminal sequences
            terminal_data = extract_terminal_sequences(intron_seq, strand)

            analyses.append({
                'mag_id': mag_id,
                'contig': contig_id,
                'start': start,
                'end': end,
                'strand': strand,
                'intron_len': terminal_data['length'],
                'is_telomeric': is_telo,
                'donor': row['donor'],
                'acceptor': row['acceptor'],
                'splice_type': row['splice_type'],
                **terminal_data
            })

            if is_telo:
                processed_telos += 1

    df_analysis = pd.DataFrame(analyses)
    print(f"  Analyzed {len(df_analysis)} introns from top 8 MAGs")

    # Filter to introns with telomeric arrays
    df_telos = df_analysis[df_analysis['telo_arrays'].apply(lambda x: len(x) > 0)].copy()
    print(f"  Introns containing telomeric arrays: {len(df_telos)}")

    # ==============================================================================
    # ANALYSIS SECTION 1: Terminal Sequence Properties
    # ==============================================================================
    print("\n[4/7] Analyzing terminal sequence properties...")

    leaders_data = []
    trailers_data = []
    boundary_leaders = []  # Last 20bp of leader + first 20bp of telo
    boundary_trailers = []  # Last 20bp of telo + first 20bp of trailer

    for idx, row in df_telos.iterrows():
        if row['leader']:
            leaders_data.append({
                'mag_id': row['mag_id'],
                'length': row['leader_len'],
                'seq': row['leader'],
                'nuc_freq': get_nucleotide_freq(row['leader']),
                'donor_site': row['donor_site'],
                'branch_points': len(find_branch_point(row['leader'])),
                'poly_tracts': len(find_polypyrimidine_tract(row['leader'])),
            })

            # Boundary region
            boundary_seq = row['leader'][-20:] if len(row['leader']) >= 20 else row['leader']
            if row['telo_start'] is not None and row['telo_start'] + 20 <= len(row['full_seq']):
                boundary_seq += row['full_seq'][row['telo_start']:row['telo_start']+20]
            boundary_leaders.append(boundary_seq)

        if row['trailer']:
            trailers_data.append({
                'mag_id': row['mag_id'],
                'length': row['trailer_len'],
                'seq': row['trailer'],
                'nuc_freq': get_nucleotide_freq(row['trailer']),
                'acceptor_site': row['acceptor_site'],
                'branch_points': len(find_branch_point(row['trailer'])),
                'poly_tracts': len(find_polypyrimidine_tract(row['trailer'])),
            })

            # Boundary region
            if row['telo_end'] is not None and row['telo_end'] >= 20:
                boundary_seq = row['full_seq'][row['telo_end']-20:row['telo_end']]
            else:
                boundary_seq = ''
            if row['trailer']:
                boundary_seq += row['trailer'][:20]
            boundary_trailers.append(boundary_seq)

    df_leaders = pd.DataFrame(leaders_data)
    df_trailers = pd.DataFrame(trailers_data)

    print(f"  Analyzed {len(df_leaders)} 5' leaders")
    print(f"  Analyzed {len(df_trailers)} 3' trailers")

    # ==============================================================================
    # ANALYSIS SECTION 2: Splice Site Mapping in TTAGGG
    # ==============================================================================
    print("\n[5/7] Mapping splice sites in TTAGGG repeat...")

    splice_map = map_splice_sites_in_ttaggg()
    print(f"  TTAGGG sequence: {splice_map['full_sequence']}")
    print(f"  GT positions within 6bp period: {splice_map['GT']}")
    print(f"  AG positions within 6bp period: {splice_map['AG']}")

    # ==============================================================================
    # ANALYSIS SECTION 3: Distribution Summary Statistics
    # ==============================================================================
    print("\n[6/7] Calculating distribution statistics...")

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
    if len(df_leaders) > 0:
        print("\n  NUCLEOTIDE COMPOSITION:")
        leader_nuc = {nt: [] for nt in ['A', 'T', 'G', 'C']}
        for freq_dict in df_leaders['nuc_freq']:
            for nt in leader_nuc:
                leader_nuc[nt].append(freq_dict[nt])
        print(f"    Leaders - GC%: {np.mean(leader_nuc['G']) + np.mean(leader_nuc['C']):.1%}")

    if len(df_trailers) > 0:
        trailer_nuc = {nt: [] for nt in ['A', 'T', 'G', 'C']}
        for freq_dict in df_trailers['nuc_freq']:
            for nt in trailer_nuc:
                trailer_nuc[nt].append(freq_dict[nt])
        print(f"    Trailers - GC%: {np.mean(trailer_nuc['G']) + np.mean(trailer_nuc['C']):.1%}")

    # Correlation analysis
    if len(df_telos) > 0:
        telos_with_both = df_telos[(df_telos['leader_len'] > 0) & (df_telos['trailer_len'] > 0)].copy()
        if len(telos_with_both) > 1:
            corr_len_purity = telos_with_both[['leader_len', 'telo_purity']].corr().iloc[0, 1]
            corr_trail_purity = telos_with_both[['trailer_len', 'telo_purity']].corr().iloc[0, 1]
            corr_len_trail = telos_with_both[['leader_len', 'trailer_len']].corr().iloc[0, 1]
            print(f"\n  CORRELATIONS:")
            print(f"    Leader length vs telomeric purity: {corr_len_purity:.3f}")
            print(f"    Trailer length vs telomeric purity: {corr_trail_purity:.3f}")
            print(f"    Leader length vs trailer length: {corr_len_trail:.3f}")

    # ==============================================================================
    # FIGURE GENERATION
    # ==============================================================================
    print("\n[7/7] Creating visualization figure...")

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # PANEL A: 5' Leader length distribution
    ax_a = fig.add_subplot(gs[0, 0])
    if len(df_leaders) > 0:
        ax_a.hist(df_leaders['length'], bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax_a.set_xlabel("Leader length (bp)", fontsize=11, fontweight='bold')
        ax_a.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_a.set_title("A. 5' Leader Length Distribution", fontsize=12, fontweight='bold')
        ax_a.axvline(df_leaders['length'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {df_leaders['length'].mean():.0f}bp")
        ax_a.legend(fontsize=9)
        ax_a.grid(axis='y', alpha=0.3)

    # PANEL B: 3' Trailer length distribution
    ax_b = fig.add_subplot(gs[0, 1])
    if len(df_trailers) > 0:
        ax_b.hist(df_trailers['length'], bins=20, color='coral', edgecolor='black', alpha=0.7)
        ax_b.set_xlabel("Trailer length (bp)", fontsize=11, fontweight='bold')
        ax_b.set_ylabel("Frequency", fontsize=11, fontweight='bold')
        ax_b.set_title("B. 3' Trailer Length Distribution", fontsize=12, fontweight='bold')
        ax_b.axvline(df_trailers['length'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {df_trailers['length'].mean():.0f}bp")
        ax_b.legend(fontsize=9)
        ax_b.grid(axis='y', alpha=0.3)

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
    seq_positions = np.arange(len(seq_ttaggg))

    # Draw sequence as boxes
    for i, (nt, pos) in enumerate(zip(seq_ttaggg, seq_positions)):
        ax_d.add_patch(Rectangle((i, 0.3), 0.8, 0.4,
                                  facecolor='lightgray', edgecolor='black', linewidth=2))
        ax_d.text(i + 0.4, 0.5, nt, ha='center', va='center', fontsize=14, fontweight='bold')

    # Mark GT sites
    for pos in splice_map['GT']:
        if isinstance(pos, int) and pos < 5:
            ax_d.plot([pos + 0.4, pos + 1.2], [0.15, 0.15], 'g-', linewidth=3, label='GT' if pos == splice_map['GT'][0] else '')
            ax_d.text(pos + 0.8, 0.05, 'GT', ha='center', fontsize=9, color='green', fontweight='bold')

    # Mark AG sites
    for pos in splice_map['AG']:
        if isinstance(pos, int) and pos < 5:
            ax_d.plot([pos + 0.4, pos + 1.2], [0.85, 0.85], 'r-', linewidth=3, label='AG' if pos == splice_map['AG'][0] else '')
            ax_d.text(pos + 0.8, 0.95, 'AG', ha='center', fontsize=9, color='red', fontweight='bold')

    ax_d.set_xlim(-0.5, 6.5)
    ax_d.set_ylim(0, 1)
    ax_d.set_aspect('equal')
    ax_d.axis('off')
    ax_d.set_title("D. Splice Site Positions in TTAGGG", fontsize=12, fontweight='bold')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='g', lw=3, label='GT (Donor)'),
                       Line2D([0], [0], color='r', lw=3, label='AG (Acceptor)')]
    ax_d.legend(handles=legend_elements, loc='upper left', fontsize=9)

    # PANEL E: Boundary region sequence logo
    ax_e = fig.add_subplot(gs[2, 0])

    # Create a simple frequency plot for boundary regions
    if boundary_leaders:
        # Focus on first 40bp of boundary (last 20 of leader + first 20 of telo)
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
            ax_e.set_title("E. 5' Boundary Region (Leader→Telomere)", fontsize=12, fontweight='bold')
            ax_e.axvline(19.5, color='red', linestyle='--', linewidth=2, label='Leader/Telomere boundary')
            ax_e.legend(fontsize=8, loc='upper right')
            plt.colorbar(im, ax=ax_e, label='Frequency')

    # PANEL F: Correlation plot
    ax_f = fig.add_subplot(gs[2, 1])

    if len(df_telos) > 0:
        telos_with_both = df_telos[(df_telos['leader_len'] > 0) & (df_telos['trailer_len'] > 0)].copy()

        if len(telos_with_both) > 1:
            scatter = ax_f.scatter(telos_with_both['leader_len'], telos_with_both['telo_purity'],
                                   c=telos_with_both['trailer_len'], cmap='viridis',
                                   s=100, alpha=0.6, edgecolor='black')
            ax_f.set_xlabel("Leader Length (bp)", fontsize=11, fontweight='bold')
            ax_f.set_ylabel("Telomeric Purity", fontsize=11, fontweight='bold')
            ax_f.set_title("F. Terminal Length vs Telomeric Purity", fontsize=12, fontweight='bold')
            cbar = plt.colorbar(scatter, ax=ax_f)
            cbar.set_label("Trailer Length (bp)", fontsize=10, fontweight='bold')
            ax_f.grid(alpha=0.3)

    plt.suptitle("Extended Data Figure 14: Telotron Terminal Sequences",
                 fontsize=14, fontweight='bold', y=0.995)

    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches='tight')
    print(f"\n  Figure saved to: {OUTPUT_FIG}")

    # ==============================================================================
    # FINAL SUMMARY
    # ==============================================================================
    print("\n" + "=" * 80)
    print("SUMMARY OF KEY FINDINGS")
    print("=" * 80)

    print(f"\nTelotron Dataset:")
    print(f"  - Top 8 MAGs analyzed: {len(TOP_8_MAGS)}")
    print(f"  - Introns with telomeric arrays: {len(df_telos)}")
    print(f"  - Introns with identifiable 5' leaders: {len(df_leaders)}")
    print(f"  - Introns with identifiable 3' trailers: {len(df_trailers)}")

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
    print(f"  - TTAGGG canonical sequence: {splice_map['full_sequence']}")
    print(f"  - Donor site (GT) located at position(s): {[p for p in splice_map['GT'] if isinstance(p, int)]}")
    print(f"  - Acceptor site (AG) located at position(s): {[p for p in splice_map['AG'] if isinstance(p, int)]}")
    print(f"  - Note: GT and AG motifs DO appear within TTAGGG period, potentially allowing")
    print(f"    intronic splice sites to occur within telomeric arrays.")

    if len(df_telos) > 0:
        print(f"\nTelomeric Array Purity:")
        telos_with_purity = df_telos[df_telos['telo_purity'] > 0]
        if len(telos_with_purity) > 0:
            print(f"  - Average purity: {telos_with_purity['telo_purity'].mean():.3f}")
            print(f"  - Median purity: {telos_with_purity['telo_purity'].median():.3f}")
            print(f"  - Range: {telos_with_purity['telo_purity'].min():.3f}-{telos_with_purity['telo_purity'].max():.3f}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)

if __name__ == '__main__':
    main()
