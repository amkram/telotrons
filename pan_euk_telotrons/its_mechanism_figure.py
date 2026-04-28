#!/usr/bin/env python3
"""
ITS-mechanism summary figure: 6-panel dashboard.

Panel A: Class distribution by dataset (stacked bar)
Panel B: Linker length distribution (histogram, by class)
Panel C: TSD presence per class (bar)
Panel D: NR2C/F motif rate per class (bar) — TTI signal
Panel E: Survivorship: telomeric density per Mb in CDS / intron / intergenic
         across 5 haptophyte MAGs (grouped bar)
Panel F: Linker composition: fraction with telo content (bar)
"""

import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import defaultdict


def main():
    real_dir = Path(__file__).parent / "real_telotrons"

    # Load summaries
    with open(real_dir / "its_mechanism_v2.json") as f:
        v2 = json.load(f)
    with open(real_dir / "linker_analysis_summary.json") as f:
        la = json.load(f)
    with open(real_dir / "true_nonrepeat_linkers_summary.json") as f:
        nrl = json.load(f)
    with open(real_dir / "survivorship_test.json") as f:
        surv = json.load(f)

    # Load per-record TSV for histogram
    linker_lens_per_class = defaultdict(list)
    with open(real_dir / "its_mechanism_v2.tsv") as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            ll = int(row.get('linker_len', 0) or 0)
            cls = row.get('classification', '')
            if ll >= 5:
                linker_lens_per_class[cls].append(ll)

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # === Panel A: Class distribution ===
    ax = axes[0, 0]
    classes = ['SINGLE_G', 'SINGLE_C', 'CONVERGING_NOLINKER', 'CONVERGING_WITHLINKER',
               'DIVERGING_WITHLINKER', 'DIVERGING_NOLINKER', 'NONE']
    pan_counts = [v2['class_distribution_by_dataset']['pan_euk'].get(c, 0) for c in classes]
    tara_counts = [v2['class_distribution_by_dataset']['tara'].get(c, 0) for c in classes]
    x = np.arange(len(classes))
    ax.bar(x, tara_counts, color='steelblue', label='Tara Oceans (haptophyte MAGs)')
    ax.bar(x, pan_counts, bottom=tara_counts, color='orange', label='Pan-eukaryotic (validated)')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', '\n') for c in classes], rotation=0, fontsize=8)
    ax.set_ylabel('count')
    ax.set_yscale('symlog')
    ax.set_title('A. Orientation class distribution', loc='left', fontweight='bold')
    ax.legend(fontsize=8)
    for i, (t, p) in enumerate(zip(tara_counts, pan_counts)):
        if t + p > 0:
            ax.text(i, t + p, f"{t+p:,}", ha='center', va='bottom', fontsize=7)

    # === Panel B: Linker length histogram ===
    ax = axes[0, 1]
    cwl = linker_lens_per_class.get('CONVERGING_WITHLINKER', [])
    if cwl:
        bins = np.linspace(5, 200, 40)
        ax.hist(cwl, bins=bins, color='steelblue', alpha=0.7, edgecolor='k', linewidth=0.3)
        ax.axvspan(5, 30, alpha=0.15, color='red', label='Short (5-30bp)\nNHEJ scar')
        ax.axvspan(30, 120, alpha=0.15, color='gold', label='Medium (30-120bp)\nTERC-frag (Nergadze)')
        ax.axvspan(120, 200, alpha=0.15, color='green', label='Long (120+bp)\nCaptured frag')
        ax.set_xlim(5, 200)
        ax.set_xlabel('Linker length (bp)')
        ax.set_ylabel('count')
        ax.set_title(f'B. CONVERGING_WITHLINKER linker lengths (n={len(cwl):,})', loc='left', fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')

    # === Panel C: TSD presence ===
    ax = axes[0, 2]
    tsd_data = v2['tsd_presence_per_class']
    keep = ['SINGLE_G', 'SINGLE_C', 'CONVERGING_NOLINKER', 'CONVERGING_WITHLINKER']
    pcts = [tsd_data.get(c, {}).get('pct_with_tsd_3plus', 0) for c in keep]
    means = [tsd_data.get(c, {}).get('mean_tsd_when_present', 0) for c in keep]
    x = np.arange(len(keep))
    bars = ax.bar(x, pcts, color='coral', edgecolor='k')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', '\n') for c in keep], fontsize=8)
    ax.set_ylabel('% with TSD ≥3bp')
    ax.set_title('C. TSDs (target site duplications)', loc='left', fontweight='bold')
    ax.set_ylim(0, 12)
    for i, (p, m) in enumerate(zip(pcts, means)):
        ax.text(i, p + 0.3, f"{p:.1f}%\n(mean {m:.1f}bp)", ha='center', va='bottom', fontsize=7)

    # === Panel D: NR2C/F motif rate ===
    ax = axes[1, 0]
    nr_data = v2['nr2cf_halfsites_per_class']
    pcts = [nr_data.get(c, {}).get('pct_with_any', 0) for c in keep]
    means = [nr_data.get(c, {}).get('mean', 0) for c in keep]
    x = np.arange(len(keep))
    bars = ax.bar(x, pcts, color='mediumpurple', edgecolor='k')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', '\n') for c in keep], fontsize=8)
    ax.set_ylabel('% with ≥1 GGGTCA half-site\nin 200bp flanks')
    ax.set_title('D. NR2C/F motif scan (TTI mechanism test)', loc='left', fontweight='bold')
    ax.axhline(y=50, color='red', linestyle='--', linewidth=1, label='TTI prediction\n(should be high)')
    ax.set_ylim(0, 60)
    for i, (p, m) in enumerate(zip(pcts, means)):
        ax.text(i, p + 1, f"{p:.1f}%", ha='center', va='bottom', fontsize=7)
    ax.legend(fontsize=8)

    # === Panel E: Survivorship test ===
    ax = axes[1, 1]
    mags = list(surv.keys())
    short_names = [m.split('_MAG_')[0].replace('TARA_', '') + '...' for m in mags]
    cds_per_mb = [surv[m]['cds']['telo_kmers'] / (surv[m]['cds']['bp'] / 1e6) if surv[m]['cds']['bp'] > 0 else 0 for m in mags]
    intron_per_mb = [surv[m]['intron']['telo_kmers'] / (surv[m]['intron']['bp'] / 1e6) if surv[m]['intron']['bp'] > 0 else 0 for m in mags]
    interg_per_mb = [surv[m]['intergenic']['telo_kmers'] / (surv[m]['intergenic']['bp'] / 1e6) if surv[m]['intergenic']['bp'] > 0 else 0 for m in mags]
    x = np.arange(len(mags))
    w = 0.27
    ax.bar(x - w, cds_per_mb, w, label='CDS', color='lightcoral', edgecolor='k')
    ax.bar(x, intron_per_mb, w, label='Intron', color='steelblue', edgecolor='k')
    ax.bar(x + w, interg_per_mb, w, label='Intergenic', color='gray', edgecolor='k')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=30, ha='right', fontsize=7)
    ax.set_ylabel('Telomeric kmers per Mb')
    ax.set_title('E. Survivorship: telomeric density by compartment', loc='left', fontweight='bold')
    ax.legend(fontsize=8)
    # Add ratio annotations
    for i, (cds, intr) in enumerate(zip(cds_per_mb, intron_per_mb)):
        if cds > 0:
            ratio = intr / cds
            ax.text(i, intr * 1.5, f"{ratio:.0f}×", ha='center', va='bottom', fontsize=7,
                    fontweight='bold', color='darkblue')

    # === Panel F: linker composition ===
    ax = axes[1, 2]
    # From linker_analysis: orientation breakdown by length bucket
    buckets_order = ['short_5_30', 'medium_30_120', 'long_120plus']
    # Pure non-repeat % from nrl
    pure_pct = nrl['pure_nonrepeat'] / nrl['total_linkers'] * 100

    # Stacked bars per length bucket
    components = ['no_telo_pct', 'fwd_only_pct', 'rev_only_pct', 'both_pct']
    colors = ['darkred', 'steelblue', 'coral', 'gold']
    labels = ['No telo (true non-repeat)', 'Only TTAGGG kmers',
              'Only CCCTAA kmers', 'Both (mosaic)']
    data = {comp: [la['orientation_breakdown'][b][comp] for b in buckets_order]
            for comp in components}
    x = np.arange(len(buckets_order))
    bottom = np.zeros(len(buckets_order))
    for comp, color, label in zip(components, colors, labels):
        ax.bar(x, data[comp], bottom=bottom, color=color, edgecolor='k', linewidth=0.5, label=label)
        bottom += np.array(data[comp])
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace('_', '\n') for b in buckets_order], fontsize=8)
    ax.set_ylabel('% of linkers in bucket')
    ax.set_title(f'F. Linker composition (only {pure_pct:.1f}% are truly non-repeat)',
                 loc='left', fontweight='bold')
    ax.legend(fontsize=7, loc='center right')
    ax.set_ylim(0, 100)

    plt.suptitle('Telotron architecture: testing ITS-mechanism hypotheses\n'
                 '55,822 telomeric introns from 8 haptophyte MAGs + 186 pan-eukaryotic validated',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()

    out = real_dir / 'its_mechanism_summary.png'
    plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
