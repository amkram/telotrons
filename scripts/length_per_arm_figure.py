#!/usr/bin/env python
"""Plot per-arm burst length comparison: 1-arm vs 2-arm/2 (telomeric_bases)."""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Repo-relative default; TELOTRON_ROOT env override keeps this runnable off /scratch1.
ROOT = os.environ.get('TELOTRON_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'work/results/length_distribution_2026-06-07')
os.makedirs(OUT, exist_ok=True)
df = pd.read_csv(os.path.join(ROOT, 'work/results/final_telotron_set_architecture.tsv'), sep='\t', low_memory=False)
def lin(row):
    org=row['organism']
    if org=='Tara Oceans MAG': return 'PSW_MAG'
    if isinstance(org,str) and org.startswith('Eimeria'): return 'Eimeria'
    return 'Other'
df['lineage']=df.apply(lin,axis=1)

fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5))
for ax, lineage, motif_len in [(axes[0],'PSW_MAG',6), (axes[1],'Eimeria',7)]:
    one = df[(df['lineage']==lineage)&(df['architecture'].isin(['GT-F-AG','GT-R-AG']))]['telomeric_bases'].astype(float).values
    two = df[(df['lineage']==lineage)&(df['architecture']=='GT-F-R-AG')]['telomeric_bases'].astype(float).values
    two_per_arm = two/2.0
    linker_arch = df[(df['lineage']==lineage)&(df['architecture'].isin(['GT-R-linker-F-AG','GT-F-linker-R-AG']))]['telomeric_bases'].astype(float).values
    linker_per_arm = linker_arch/2.0
    bins = np.arange(0, 250, 7 if motif_len==7 else 6)
    if len(one):
        ax.hist(one, bins=bins, alpha=0.5, color='#1f78b4', density=True, label=f'1-arm (n={len(one)}) med={int(np.median(one))}bp')
    if len(two):
        ax.hist(two_per_arm, bins=bins, alpha=0.5, color='#e31a1c', density=True, label=f'2-arm /2 (n={len(two)}) med={int(np.median(two_per_arm))}bp')
    if len(linker_per_arm):
        ax.hist(linker_per_arm, bins=bins, alpha=0.4, color='#33a02c', density=True, label=f'linker /2 (n={len(linker_arch)}) med={int(np.median(linker_per_arm))}bp')
    ax.set_xlabel(f'per-arm telomeric_bases (bp; ticks every {motif_len}bp = 1 repeat)')
    ax.set_ylabel('density')
    ax.set_title(f'{lineage} (motif_len={motif_len})')
    ax.legend(fontsize=7, loc='upper right')
plt.tight_layout()
out = os.path.join(OUT, 'per_arm_burst_length.pdf')
plt.savefig(out)
print('wrote', out)
