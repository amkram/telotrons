#!/usr/bin/env python
"""
Length-distribution analysis: telotron intron_len by (lineage x architecture).

- Splits PSW Tara MAG vs Eimeria genus.
- Per (lineage, architecture) bucket: count, median, IQR, range, Hartigan dip stat + p.
- Linker-bearing vs no-linker comparison (Mann-Whitney) within each lineage.
- Convergent two-arm hypothesis: does GT-F-R-AG have a mode at ~2x GT-F-AG (or
  ~2x GT-R-AG)? Compute median ratios.
- Hypothesis on linker telotrons: are they longer by ~linker_len than the
  matching no-linker class?

Outputs go under work/results/length_distribution_2026-06-07/.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import diptest

def bh_fdr(pvals):
    """Benjamini-Hochberg q-values (NaN-safe, order-preserving)."""
    p = np.asarray(pvals, dtype=float); ok = ~np.isnan(p); idx = np.where(ok)[0]
    q = np.full(p.shape, np.nan)
    if idx.size == 0:
        return q
    ps = p[idx]; order = np.argsort(ps); m = idx.size
    ranked = ps[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    qi = np.empty(m); qi[order] = ranked; q[idx] = qi
    return q

# Repo-relative default; CLI/env override keeps this runnable off /scratch1.
ROOT = os.environ.get('TELOTRON_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(ROOT, 'work/results/length_distribution')
os.makedirs(OUTDIR, exist_ok=True)

TSV = os.path.join(ROOT, 'work/results/final_telotron_set_architecture.tsv')

df = pd.read_csv(TSV, sep='\t', low_memory=False)

# Lineage assignment
def lineage(row):
    org = row['organism']
    if org == 'Tara Oceans MAG':
        return 'PSW_MAG'
    if isinstance(org, str) and org.startswith('Eimeria'):
        return 'Eimeria'
    if org == 'Monocercomonoides exilis':
        return 'Monocercomonoides'
    return 'Other'

df['lineage'] = df.apply(lineage, axis=1)

# Focus architectures
ARCHS = ['GT-F-AG', 'GT-R-AG', 'GT-F-R-AG', 'GT-R-F-AG', 'GT-R-linker-F-AG', 'GT-F-linker-R-AG', 'Multi-junction', 'Other']

# ===== Per-bucket summary =====
rows = []
for lin in ['PSW_MAG', 'Eimeria']:
    for arch in ARCHS:
        sub = df[(df['lineage']==lin) & (df['architecture']==arch)]
        n = len(sub)
        if n == 0:
            rows.append({'lineage':lin, 'architecture':arch, 'n':0, 'median':np.nan, 'q25':np.nan, 'q75':np.nan, 'min':np.nan, 'max':np.nan, 'mean':np.nan, 'std':np.nan, 'dip_stat':np.nan, 'dip_p':np.nan})
            continue
        x = sub['intron_len'].astype(float).values
        med = float(np.median(x)); q25, q75 = float(np.percentile(x,25)), float(np.percentile(x,75))
        mn, mx = float(np.min(x)), float(np.max(x))
        # Operator precedence: `A, B if n>1 else 0.0` parses as `A, (B if n>1 else 0.0)`
        # which is a tuple, then a tuple/float mix on unpack; guard the whole line.
        mean = float(np.mean(x))
        std = float(np.std(x, ddof=1)) if n > 1 else 0.0
        # Hartigan dip on log10(len) to be scale-insensitive when distributions span orders
        # but report on raw too; choose raw because biological lengths are bounded short.
        if n >= 4:
            try:
                dip_stat, dip_p = diptest.diptest(x)
            except Exception as e:
                dip_stat, dip_p = (np.nan, np.nan)
        else:
            dip_stat, dip_p = (np.nan, np.nan)
        rows.append({'lineage':lin, 'architecture':arch, 'n':n, 'median':med, 'q25':q25, 'q75':q75,
                     'min':mn, 'max':mx, 'mean':mean, 'std':std, 'dip_stat':dip_stat, 'dip_p':dip_p})

summary = pd.DataFrame(rows)
summary['dip_q_bh'] = bh_fdr(summary['dip_p'].values)  # BH across all (lineage x arch) dip tests
summary.to_csv(os.path.join(OUTDIR,'summary_by_lineage_arch.tsv'), sep='\t', index=False)
print('Summary (dip_q_bh = BH-FDR over buckets):')
print(summary.to_string(index=False))

# ===== Linker vs no-linker MWU within lineage =====
def has_linker(a):
    return 'linker' in a if isinstance(a, str) else False

rows = []
for lin in ['PSW_MAG', 'Eimeria']:
    sub = df[df['lineage']==lin]
    linker = sub[sub['architecture'].apply(has_linker)]['intron_len'].astype(float).values
    nolink = sub[sub['architecture'].isin(['GT-F-R-AG','GT-F-AG','GT-R-AG'])]['intron_len'].astype(float).values
    if len(linker)>=5 and len(nolink)>=5:
        u, p = mannwhitneyu(linker, nolink, alternative='greater')
        rows.append({'lineage':lin, 'n_linker':len(linker), 'n_nolink':len(nolink),
                     'med_linker':float(np.median(linker)), 'med_nolink':float(np.median(nolink)),
                     'delta_med':float(np.median(linker)-np.median(nolink)),
                     'mwu_U':float(u), 'mwu_p_greater':float(p)})
    else:
        rows.append({'lineage':lin, 'n_linker':len(linker), 'n_nolink':len(nolink),
                     'med_linker':np.nan, 'med_nolink':np.nan, 'delta_med':np.nan, 'mwu_U':np.nan, 'mwu_p_greater':np.nan})

linker_cmp = pd.DataFrame(rows)
# NB: PSW_MAG counts are intervals from ONE genome (TARA_PSW_86_MAG_00284) and Eimeria pools 5
# assemblies, so these p's reflect within-assembly variance, not biological replication. Treat as
# descriptive; BH q corrects only for the number of tests run here.
linker_cmp['mwu_q_bh'] = bh_fdr(linker_cmp['mwu_p_greater'].values)
linker_cmp.to_csv(os.path.join(OUTDIR,'linker_vs_nolinker.tsv'), sep='\t', index=False)
print('\nLinker vs no-linker (intron length; mwu_q_bh = BH-FDR; single-genome PSW caveat):')
print(linker_cmp.to_string(index=False))

# ===== Convergent two-arm hypothesis =====
# Median ratio: GT-F-R-AG median / median of (GT-F-AG, GT-R-AG); should be ~2 under independent extensions.
rows = []
for lin in ['PSW_MAG','Eimeria']:
    sub = df[df['lineage']==lin]
    def m(arch):
        v = sub[sub['architecture']==arch]['intron_len'].values
        return float(np.median(v)) if len(v) else np.nan
    fr = m('GT-F-R-AG'); f = m('GT-F-AG'); r = m('GT-R-AG')
    one_arm_med = np.nanmedian([f, r])
    rows.append({'lineage':lin, 'med_GT-F-AG':f, 'med_GT-R-AG':r, 'med_one_arm':one_arm_med,
                 'med_GT-F-R-AG':fr, 'ratio_FR_over_onearm': fr/one_arm_med if one_arm_med else np.nan})

twoarm = pd.DataFrame(rows)
twoarm.to_csv(os.path.join(OUTDIR,'two_arm_hypothesis_ratio.tsv'), sep='\t', index=False)
print('\nTwo-arm hypothesis (median GT-F-R-AG / one-arm median):')
print(twoarm.to_string(index=False))

# ===== Linker offset hypothesis =====
# (linker_len + median_one_arm_arch)  vs   median_linker_arch
rows = []
for lin in ['PSW_MAG','Eimeria']:
    sub = df[df['lineage']==lin]
    # GT-R-linker-F-AG vs GT-R-AG + linker  AND vs GT-F-R-AG + linker
    for arch_link in ['GT-R-linker-F-AG','GT-F-linker-R-AG']:
        s = sub[sub['architecture']==arch_link]
        if len(s)==0: continue
        med_link_arch = float(np.median(s['intron_len']))
        med_link = float(np.median(pd.to_numeric(s['linker_len'], errors='coerce').dropna()))
        # corresponding no-linker
        if arch_link=='GT-R-linker-F-AG':
            base = sub[sub['architecture']=='GT-F-R-AG']['intron_len']
        else:
            base = sub[sub['architecture']=='GT-F-R-AG']['intron_len']
        med_base = float(np.median(base)) if len(base) else np.nan
        rows.append({'lineage':lin, 'arch_link':arch_link, 'n':len(s),
                     'med_link_arch':med_link_arch, 'med_linker_len':med_link,
                     'med_GT-F-R-AG':med_base,
                     'observed_delta': med_link_arch - med_base,
                     'predicted_delta_eq_linker_len': med_link})
linker_offset = pd.DataFrame(rows)
linker_offset.to_csv(os.path.join(OUTDIR,'linker_offset_hypothesis.tsv'), sep='\t', index=False)
print('\nLinker offset hypothesis:')
print(linker_offset.to_string(index=False))

# ===== Figures =====
# 2x1 panel: histogram per (lineage, arch). Log-x.
def plot_histograms():
    # Multi-junction + GT-R-F-AG included: both are bidirectional classes that can
    # single-handedly admit a confident bearer, and omitting them made these
    # figures read as a complete architecture census when they were not.
    archs_plot = ['GT-F-AG','GT-R-AG','GT-F-R-AG','GT-R-F-AG','GT-R-linker-F-AG','GT-F-linker-R-AG','Multi-junction','Other']
    lineages = ['PSW_MAG','Eimeria']
    fig, axes = plt.subplots(2, len(archs_plot), figsize=(3.3*len(archs_plot), 6), sharey=False)
    for i, lin in enumerate(lineages):
        for j, arch in enumerate(archs_plot):
            ax = axes[i,j]
            v = df[(df['lineage']==lin) & (df['architecture']==arch)]['intron_len'].astype(float).values
            if len(v)==0:
                ax.set_title(f'{lin}/{arch}\nn=0', fontsize=8)
                ax.axis('off'); continue
            v_clip = v[(v>0) & (v<2000)]
            ax.hist(v_clip, bins=40, color='#4a7ab8' if i==0 else '#b8584a', alpha=0.85)
            ax.axvline(np.median(v_clip), color='k', ls='--', lw=0.8)
            try:
                ds, dp = diptest.diptest(v) if len(v)>=4 else (np.nan,np.nan)
            except Exception:
                ds, dp = (np.nan, np.nan)
            ax.set_title(f'{lin}\n{arch}\nn={len(v)} med={int(np.median(v))} dip_p={dp:.3g}', fontsize=7)
            ax.set_xlabel('intron_len (bp)', fontsize=7)
            ax.tick_params(labelsize=6)
    plt.tight_layout()
    p = os.path.join(OUTDIR,'length_hist_by_lineage_arch.pdf')
    plt.savefig(p)
    plt.close()
    return p

# Overlaid KDE/hist comparison: PSW vs Eimeria for each architecture
def plot_overlay():
    archs_plot = ['GT-F-AG','GT-R-AG','GT-F-R-AG','GT-R-F-AG','GT-R-linker-F-AG','GT-F-linker-R-AG','Multi-junction']
    fig, axes = plt.subplots(1, len(archs_plot), figsize=(3.3*len(archs_plot), 3.2), sharey=False)
    for j, arch in enumerate(archs_plot):
        ax = axes[j]
        for color, lin in [('#4a7ab8','PSW_MAG'), ('#b8584a','Eimeria')]:
            v = df[(df['lineage']==lin) & (df['architecture']==arch)]['intron_len'].astype(float).values
            if len(v)<3: continue
            v_clip = v[(v>0) & (v<1500)]
            ax.hist(v_clip, bins=40, color=color, alpha=0.45, label=f'{lin} (n={len(v)}, med={int(np.median(v))})', density=True)
        ax.set_title(arch, fontsize=8)
        ax.legend(fontsize=6, loc='upper right')
        ax.set_xlabel('intron_len (bp)', fontsize=7)
        ax.tick_params(labelsize=6)
    plt.tight_layout()
    p = os.path.join(OUTDIR,'length_overlay_psw_vs_eimeria.pdf')
    plt.savefig(p); plt.close()
    return p

p1 = plot_histograms()
p2 = plot_overlay()
print(f'\nFigures: {p1}\n         {p2}')

# ===== PSW vs Eimeria per architecture: MWU =====
rows = []
for arch in ARCHS:
    psw = df[(df['lineage']=='PSW_MAG') & (df['architecture']==arch)]['intron_len'].astype(float).values
    eim = df[(df['lineage']=='Eimeria') & (df['architecture']==arch)]['intron_len'].astype(float).values
    if len(psw)>=5 and len(eim)>=5:
        u, p = mannwhitneyu(psw, eim, alternative='two-sided')
        rows.append({'architecture':arch, 'n_PSW':len(psw), 'n_Eim':len(eim),
                     'med_PSW':float(np.median(psw)), 'med_Eim':float(np.median(eim)),
                     'mwu_U':float(u), 'mwu_p':float(p)})
    else:
        rows.append({'architecture':arch, 'n_PSW':len(psw), 'n_Eim':len(eim),
                     'med_PSW':float(np.median(psw)) if len(psw) else np.nan,
                     'med_Eim':float(np.median(eim)) if len(eim) else np.nan,
                     'mwu_U':np.nan, 'mwu_p':np.nan})

cmp_lin = pd.DataFrame(rows)
cmp_lin['mwu_q_bh'] = bh_fdr(cmp_lin['mwu_p'].values)  # BH across the 7 per-architecture tests
cmp_lin.to_csv(os.path.join(OUTDIR,'psw_vs_eimeria_per_arch.tsv'), sep='\t', index=False)
print('\nPSW vs Eimeria per architecture (MWU two-sided; mwu_q_bh = BH-FDR over architectures):')
print(cmp_lin.to_string(index=False))

# ===== Per-arm burst-length figure =====
# 1-arm vs 2-arm/2 (telomeric_bases); per-arm burst-length comparison, ticked
# at the motif-length step so each tick = one repeat unit.
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5))
for ax, lin, motif_len in [(axes[0], 'PSW_MAG', 6), (axes[1], 'Eimeria', 7)]:
    sub = df[df['lineage'] == lin]
    one = sub[sub['architecture'].isin(['GT-F-AG', 'GT-R-AG'])]['telomeric_bases'].astype(float).values
    two = sub[sub['architecture'] == 'GT-F-R-AG']['telomeric_bases'].astype(float).values / 2.0
    link = sub[sub['architecture'].isin(['GT-R-linker-F-AG', 'GT-F-linker-R-AG'])]['telomeric_bases'].astype(float).values / 2.0
    bins = np.arange(0, 250, motif_len)
    if len(one):
        ax.hist(one, bins=bins, alpha=0.5, color='#1f78b4', density=True,
                label=f'1-arm (n={len(one)}) med={int(np.median(one))}bp')
    if len(two):
        ax.hist(two, bins=bins, alpha=0.5, color='#e31a1c', density=True,
                label=f'2-arm /2 (n={2*len(two)//1}) med={int(np.median(two))}bp')
    if len(link):
        ax.hist(link, bins=bins, alpha=0.4, color='#33a02c', density=True,
                label=f'linker /2 (n={2*len(link)//1}) med={int(np.median(link))}bp')
    ax.set_xlabel(f'per-arm telomeric_bases (bp; ticks every {motif_len}bp = 1 repeat)')
    ax.set_ylabel('density')
    ax.set_title(f'{lin} (motif_len={motif_len})')
    ax.legend(fontsize=7, loc='upper right')
plt.tight_layout()
p = os.path.join(OUTDIR, 'per_arm_burst_length.pdf')
plt.savefig(p); plt.close()
print(f'\nwrote {p}')
