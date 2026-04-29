#!/usr/bin/env python3
"""
Summary stats figure for §9f: aggregate counts + per-linker classification +
hit pident-vs-length scatter colored by lineage.
"""

import csv, glob, json, bisect
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

OUT = Path("/scratch1/alex/telotrons/pan_euk_telotrons/real_telotrons/eimeria_linker_blast/figures/summary.png")

# Build full Eimeria telotron interval map
by_acc_contig = defaultdict(lambda: defaultdict(list))
for f in sorted(glob.glob('/scratch1/alex/telotrons/pan_euk_telotrons/ultra_results/GCF_000499*_telotrons.tsv')):
    if 'Phaseolus' in f: continue
    with open(f) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for r in rdr:
            by_acc_contig[r['genome_acc']][r['contig']].append((int(r['start']), int(r['end'])))
for acc in by_acc_contig:
    for c in by_acc_contig[acc]:
        by_acc_contig[acc][c].sort()

EIM_ACCS = set(by_acc_contig.keys())
SARC_ACCS = {'GCF_000006565.2','GCF_002563875.1','GCF_000208865.1','GCA_000727475.1'}
CYCLO = 'GCF_002999335.1'
DEEP = {'GCF_000006425.1','GCF_900005855.1'}

def overlaps(acc, contig, lo, hi):
    intervals = by_acc_contig.get(acc, {}).get(contig, [])
    if not intervals: return False
    starts = [s for s,e in intervals]
    i = bisect.bisect_right(starts, hi)
    for j in range(max(0, i-2), i):
        s, e = intervals[j]
        if not (hi < s or lo > e): return True
    return False

# Per-linker classification + per-hit scatter data
per_lin = defaultdict(lambda: {'eim_in':0,'eim_out':0,'sarc':0,'cyclo':0,'deep':0})
scatter_pts = []  # (length, pident, lineage_class)
with open('/scratch1/alex/telotrons/pan_euk_telotrons/real_telotrons/eimeria_linker_blast/eimeria_linkers_vs_apicomplexa.blastn.tsv') as f:
    for line in f:
        p = line.rstrip().split('\t')
        if len(p) < 12: continue
        qseqid, sseqid = p[0], p[1]
        pident, length = float(p[2]), int(p[3])
        if pident < 80 or length < 30: continue
        qacc = qseqid.split('__',1)[0]
        sacc = sseqid.split('::',1)[0] if '::' in sseqid else sseqid
        if sacc == qacc: continue
        sstart, send = int(p[8]), int(p[9])
        s_lo, s_hi = min(sstart,send), max(sstart,send)
        contig = sseqid.split('::',1)[1] if '::' in sseqid else sseqid
        if sacc in EIM_ACCS:
            if overlaps(sacc, contig, s_lo, s_hi):
                per_lin[qseqid]['eim_in'] += 1
                lineage = 'eim_in_telotron'
            else:
                per_lin[qseqid]['eim_out'] += 1
                lineage = 'eim_out_telotron'
        elif sacc in SARC_ACCS:
            per_lin[qseqid]['sarc'] += 1
            lineage = 'sarcocystidae'
        elif sacc == CYCLO:
            per_lin[qseqid]['cyclo'] += 1
            lineage = 'cyclospora'
        elif sacc in DEEP:
            per_lin[qseqid]['deep'] += 1
            lineage = 'deep'
        else:
            continue
        scatter_pts.append((length, pident, lineage))

# Aggregate
n_with_any = sum(1 for v in per_lin.values() if any(v.values()))

# Lineage tally (linker counts)
n_eim = sum(1 for v in per_lin.values() if v['eim_in'] or v['eim_out'])
n_sarc = sum(1 for v in per_lin.values() if v['sarc'])
n_cyclo = sum(1 for v in per_lin.values() if v['cyclo'])
n_deep = sum(1 for v in per_lin.values() if v['deep'])

# Within Eimeria-hit linkers, classify by in-telotron vs out-telotron
n_only_in = sum(1 for v in per_lin.values() if v['eim_in']>0 and v['eim_out']==0)
n_only_out = sum(1 for v in per_lin.values() if v['eim_in']==0 and v['eim_out']>0)
n_mixed = sum(1 for v in per_lin.values() if v['eim_in']>0 and v['eim_out']>0)

# 3-panel layout
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], hspace=0.35, wspace=0.25)

# === Panel A: aggregate funnel ===
axA = fig.add_subplot(gs[0, 0])
axA.set_title("A. Linker→Apicomplexa BLAST aggregate", fontsize=11, fontweight="bold", loc="left")
labels = [
    f"759\npure non-repeat\nEimeria linkers",
    f"378\nwith ≥1 strong hit\n(≥80% id, ≥30 bp)",
    f"371\nin other Eimeria",
    f"30\nin Sarcocystidae",
    f"17\nin deep\nApicomplexa",
    f"9\nin Cyclospora",
]
counts = [759, 378, 371, 30, 17, 9]
colors = ["#cccccc", "#80a6c8", "#5c9b5c", "#c87f7f", "#a07ec0", "#c8a060"]
y = 0
for i, (lbl, n, c) in enumerate(zip(labels, counts, colors)):
    width = n / 759 * 6
    axA.barh(y, width, color=c, edgecolor="#444", height=0.7)
    axA.text(width + 0.1, y, lbl, va="center", ha="left", fontsize=9)
    y -= 1
axA.set_xlim(0, 9)
axA.set_ylim(-len(labels), 1)
axA.axis("off")

# === Panel B: per-linker classification within Eimeria hits ===
axB = fig.add_subplot(gs[0, 1])
axB.set_title("B. Where do cross-Eimeria hits land? (n=371 linkers with Eimeria hits)",
              fontsize=11, fontweight="bold", loc="left")
classes = [
    f"Hits ONLY at non-telotron\nlocations in sister Eimeria\n(ancestral DNA candidate)\nn={n_only_out}",
    f"Mixed: some hits in telotrons,\nsome plain DNA\nn={n_mixed}",
    f"Hits ONLY inside other\ntelotrons in sister Eimeria\n(shared telotron)\nn={n_only_in}",
]
sizes = [n_only_out, n_mixed, n_only_in]
colors_p = ["#5c9b5c", "#c8a060", "#7090c0"]
wedges, _ = axB.pie(sizes, labels=None, colors=colors_p,
                     startangle=90, wedgeprops={"edgecolor": "white", "linewidth": 2})
for w, lbl, sz in zip(wedges, classes, sizes):
    ang = (w.theta2 + w.theta1) / 2
    x = 1.25 * np.cos(np.radians(ang))
    y_ = 1.25 * np.sin(np.radians(ang))
    ha = "left" if x >= 0 else "right"
    axB.text(x, y_, lbl, ha=ha, va="center", fontsize=8.5)

# === Panel C: hit length vs pident scatter ===
axC = fig.add_subplot(gs[1, :])
axC.set_title("C. BLAST hit length vs identity, colored by lineage class",
              fontsize=11, fontweight="bold", loc="left")
class_colors = {
    'eim_in_telotron':  '#7090c0',
    'eim_out_telotron': '#5c9b5c',
    'cyclospora':       '#c8a060',
    'sarcocystidae':    '#c87f7f',
    'deep':             '#a07ec0',
}
class_labels = {
    'eim_in_telotron':  f'sister Eimeria, INSIDE telotron',
    'eim_out_telotron': f'sister Eimeria, plain DNA (ancestral)',
    'cyclospora':       'Cyclospora (Eimeriidae)',
    'sarcocystidae':    'Sarcocystidae (~150 My)',
    'deep':             'Plasmodium / Cryptosporidium (>400 My)',
}

# Group by class
from collections import defaultdict
grouped = defaultdict(list)
for length, pid, cls in scatter_pts:
    grouped[cls].append((length, pid))

# Subsample dense classes
import random
random.seed(0)
plot_order = ['deep', 'sarcocystidae', 'cyclospora', 'eim_out_telotron', 'eim_in_telotron']
for cls in plot_order:
    pts = grouped.get(cls, [])
    if not pts: continue
    if len(pts) > 3000:
        pts = random.sample(pts, 3000)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    n_pts = len(grouped[cls])
    axC.scatter(xs, ys, s=10, alpha=0.45, c=class_colors[cls],
                label=f"{class_labels[cls]} (n={n_pts} hits)",
                edgecolors='none')

axC.set_xscale("log")
axC.set_xlabel("BLAST hit alignment length (bp, log scale)", fontsize=10)
axC.set_ylabel("Percent identity", fontsize=10)
axC.set_xlim(28, 1200)
axC.set_ylim(78, 102)
axC.axhline(80, color="#888", linestyle=":", linewidth=0.5)
axC.axvline(30, color="#888", linestyle=":", linewidth=0.5)
axC.legend(loc="lower right", fontsize=8.5, framealpha=0.9)
axC.grid(True, alpha=0.3)

fig.suptitle(
    "§9f Eimeria pure non-repeat linker → 12-Apicomplexa-genome BLAST",
    fontsize=13, fontweight="bold", y=0.99,
)

plt.savefig(OUT, dpi=180, bbox_inches="tight")
plt.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
plt.close()
print(f"wrote {OUT}")
