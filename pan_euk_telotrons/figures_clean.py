#!/usr/bin/env python3
"""
Simplified sequence-level figures.

Outputs:
  fig_classes.png         5 orientation classes, 4 examples each
  fig_compare.png         hapt CONVERGING vs Eimeria DIVERGING
  fig_shared.png          5 shared-linker tandem pairs
  fig_linkers.png         linker examples (mosaic vs pure)
  fig_tandem.png          tandem cluster context
  fig_origins.png         linker→origin alignments with 200bp flanks
"""
import csv, json, re, random
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

from sequence_viewer import classify_bases, render_row, COLOR_BG

random.seed(42)


def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

FAMILIES = {
    'TTAGGG':  {'fwd': _rotations('TTAGGG'),  'rev': _rotations('CCCTAA')},
    'TTTAGGG': {'fwd': _rotations('TTTAGGG'), 'rev': _rotations('CCCTAAA')},
}


def find_runs(seq, kmers, min_run=10):
    n = len(seq); cov = bytearray(n)
    for k in kmers:
        klen = len(k); i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1: break
            for j in range(idx, min(idx+klen, n)): cov[j] = 1
            i = idx + 1
    runs = []; in_r = False; rs = 0
    for i in range(n):
        if cov[i] and not in_r: in_r=True; rs=i
        elif not cov[i] and in_r:
            in_r=False
            if i-rs >= min_run: runs.append((rs,i))
    if in_r and n-rs >= min_run: runs.append((rs,n))
    return runs


def classify_orientation(seq, family):
    sets = FAMILIES.get(family)
    if not sets: return 'NONE', None
    fwd = find_runs(seq, sets['fwd']); rev = find_runs(seq, sets['rev'])
    n = len(seq)
    fbp = sum(e-s for s,e in fwd); rbp = sum(e-s for s,e in rev)
    has_f = fbp >= 10 and fbp/n >= 0.05
    has_r = rbp >= 10 and rbp/n >= 0.05
    if has_f and not has_r: return 'SINGLE_G', None
    if has_r and not has_f: return 'SINGLE_C', None
    if not has_f and not has_r: return 'NONE', None
    fmax = max(fwd, key=lambda r: r[1]-r[0])
    rmax = max(rev, key=lambda r: r[1]-r[0])
    if fmax[0] < rmax[0]:
        gap = rmax[0] - fmax[1]
        return ('CONVERGING_NOLINKER' if gap <= 4 else 'CONVERGING_WITHLINKER',
                (fmax[1], rmax[0]) if gap > 4 else None)
    else:
        gap = fmax[0] - rmax[1]
        return ('DIVERGING_NOLINKER' if gap <= 4 else 'DIVERGING_WITHLINKER',
                (rmax[1], fmax[0]) if gap > 4 else None)


def has_telo_kmer(seq, family):
    kmers = FAMILIES[family]['fwd'] | FAMILIES[family]['rev']
    return any(k in seq for k in kmers)


def find_tsd(left, right, min_len=3, max_len=20):
    if not left or not right: return 0
    L = min(max_len, len(left), len(right))
    for tl in range(L, min_len-1, -1):
        if left[-tl:] == right[:tl]: return tl
    return 0


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


# ----------------- loaders -----------------

def load_tara():
    out = []
    p = Path('real_telotrons/_tara_oceans_telotrons.tsv')
    if not p.exists(): return out
    with open(p) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            if r.get('intron_seq'):
                r['_family'] = 'TTAGGG'
                out.append(r)
    return out


def load_eimeria_with_flanks():
    contigs = {}
    eim = []
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                if row.get('intron_seq'):
                    row['acc'] = acc
                    row['_family'] = 'TTTAGGG'
                    eim.append(row)
    for acc in sorted(set(r['acc'] for r in eim)):
        contigs[acc] = {}
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen: continue
        cur, parts = None, []
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contigs[acc][cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]
                    parts = []
                else: parts.append(line.strip().upper())
            if cur: contigs[acc][cur] = ''.join(parts)
    for r in eim:
        cs = contigs.get(r['acc'], {}).get(r['contig'])
        if not cs: continue
        try:
            s = int(r['start']); e = int(r['end'])
        except: continue
        r['left_flank_100bp'] = cs[max(0, s-200):s]
        r['right_flank_100bp'] = cs[e:e+200]
    return eim, contigs


# ----------------- core renderer -----------------

def draw_intron(ax, y, intron, left, right, family, char_w=1.0, flank_bp=40,
                fontsize=5.5, tsd_len=0, highlights=None, x0=0):
    """Render: [left] [splice5'] [intron] [splice3'] [right]"""
    lf = left[-flank_bp:] if left else ''
    rf = right[:flank_bp] if right else ''
    x = x0
    if lf:
        cls = ['flank']*len(lf)
        ov = [(len(lf)-tsd_len, len(lf), 'tsd')] if tsd_len else []
        render_row(ax, y, lf, cls, x_start=x, char_width=char_w,
                   highlight_ranges=ov, fontsize=fontsize)
        x += len(lf)*char_w
    ax.add_patch(Rectangle((x, y), char_w*0.6, 1.0, facecolor='black'))
    x += char_w*0.6
    cls = classify_bases(intron, family)
    render_row(ax, y, intron, cls, x_start=x, char_width=char_w,
               highlight_ranges=highlights or [], fontsize=fontsize)
    x += len(intron)*char_w
    ax.add_patch(Rectangle((x, y), char_w*0.6, 1.0, facecolor='black'))
    x += char_w*0.6
    if rf:
        cls = ['flank']*len(rf)
        ov = [(0, tsd_len, 'tsd')] if tsd_len else []
        render_row(ax, y, rf, cls, x_start=x, char_width=char_w,
                   highlight_ranges=ov, fontsize=fontsize)
        x += len(rf)*char_w
    return x


def add_legend_simple(ax, x, y, fontsize=8):
    items = [('fwd_canonical', 'TTAGGG/TTTAGGG'),
             ('rev_canonical', 'CCCTAA/CCCTAAA'),
             ('fwd_variant',   'fwd 1-mm'),
             ('rev_variant',   'rev 1-mm'),
             ('flank',         'flank/intergenic'),
             ('tsd',           'TSD'),
             ('shared',        'BLAST hit / shared')]
    for i, (cls, lbl) in enumerate(items):
        ax.add_patch(Rectangle((x + i*22, y), 2.5, 1, facecolor=COLOR_BG[cls],
                                edgecolor='black', linewidth=0.3))
        ax.text(x + i*22 + 3.5, y + 0.5, lbl, fontsize=fontsize, va='center')


# ----------------- FIGURES -----------------

def fig_classes():
    tara = load_tara()
    eim, _ = load_eimeria_with_flanks()
    by_class = defaultdict(list)
    for r in tara + eim:
        seq = r.get('intron_seq', '').upper()
        if not seq or not (40 <= len(seq) <= 250): continue
        if not r.get('left_flank_100bp') or not r.get('right_flank_100bp'): continue
        cls, link = classify_orientation(seq, r['_family'])
        if cls == 'NONE': continue
        by_class[cls].append({'rec': r, 'seq': seq, 'family': r['_family'], 'link': link})

    classes = ['SINGLE_G', 'SINGLE_C', 'CONVERGING_NOLINKER',
               'CONVERGING_WITHLINKER', 'DIVERGING_WITHLINKER']
    selected = {}
    for c in classes:
        cands = sorted(by_class.get(c, []), key=lambda x: len(x['seq']))
        if 'WITHLINKER' in c:
            cands = [x for x in cands if x['link'] and 15 <= x['link'][1]-x['link'][0] <= 60]
        if len(cands) > 4:
            step = len(cands) // 4
            cands = [cands[i*step] for i in range(4)]
        selected[c] = cands

    char_w = 1.0; flank_bp = 40
    max_intron = max((len(x['seq']) for c in selected.values() for x in c), default=200)
    row_w = flank_bp*2 + max_intron + 5
    n_rows = sum(len(v)+1 for v in selected.values())

    fig, ax = plt.subplots(figsize=(max(20, row_w*0.06), max(11, n_rows*0.4)))
    ax.set_xlim(-25, row_w+5); ax.set_ylim(n_rows+2, -1)
    ax.axis('off')

    cur_y = 0
    for c in classes:
        ax.text(-25, cur_y+0.5, c, ha='left', va='center', fontsize=11, fontweight='bold')
        cur_y += 1
        for ex in selected[c]:
            r = ex['rec']
            tsd = find_tsd(r['left_flank_100bp'][-flank_bp:], r['right_flank_100bp'][:flank_bp])
            draw_intron(ax, cur_y, ex['seq'], r['left_flank_100bp'], r['right_flank_100bp'],
                        ex['family'], char_w=char_w, flank_bp=flank_bp, fontsize=5.5, tsd_len=tsd)
            sp = r.get('species', r.get('mag', ''))
            sp = sp.replace('TARA_', '').split('_MAG_')[0]
            sp = sp.replace('Eimeria_', 'E.')
            ax.text(-25, cur_y+0.5, f"{sp[:14]} {len(ex['seq'])}bp",
                    ha='left', va='center', fontsize=7)
            cur_y += 1
        cur_y += 0.5

    add_legend_simple(ax, 0, n_rows+0.5, fontsize=8)
    plt.savefig('real_telotrons/fig_classes.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("fig_classes.png")


def fig_compare(n_each=8):
    tara = load_tara()
    eim, _ = load_eimeria_with_flanks()
    hapt = []
    for r in tara:
        seq = r.get('intron_seq', '').upper()
        if not seq or not (80 <= len(seq) <= 200): continue
        if not r.get('left_flank_100bp') or not r.get('right_flank_100bp'): continue
        cls, link = classify_orientation(seq, 'TTAGGG')
        if cls == 'CONVERGING_WITHLINKER' and link and 15 <= link[1]-link[0] <= 60:
            hapt.append({'rec': r, 'seq': seq})
    random.shuffle(hapt); hapt = hapt[:n_each]

    eim_div = []
    for r in eim:
        seq = r.get('intron_seq', '').upper()
        if not seq or not (80 <= len(seq) <= 200): continue
        if not r.get('left_flank_100bp') or not r.get('right_flank_100bp'): continue
        cls, link = classify_orientation(seq, 'TTTAGGG')
        if cls == 'DIVERGING_WITHLINKER' and link and 15 <= link[1]-link[0] <= 60:
            eim_div.append({'rec': r, 'seq': seq})
    random.shuffle(eim_div); eim_div = eim_div[:n_each]

    char_w = 1.0; flank_bp = 35
    col_max = max(len(x['seq']) for x in hapt + eim_div)
    col_w = flank_bp*2 + col_max
    gap = 30
    n_rows = max(len(hapt), len(eim_div))

    fig, ax = plt.subplots(figsize=(28, max(8, n_rows*0.6)+1))
    ax.set_xlim(-22, 2*col_w+gap+5); ax.set_ylim(n_rows+1.5, -2)
    ax.axis('off')

    ax.text(col_w/2, -1, "Haptophyte CONVERGING", ha='center', fontsize=12, fontweight='bold')
    eim_x = col_w + gap
    ax.text(eim_x + col_w/2, -1, "Eimeria DIVERGING", ha='center', fontsize=12, fontweight='bold')

    for i, ex in enumerate(hapt):
        r = ex['rec']
        tsd = find_tsd(r['left_flank_100bp'][-flank_bp:], r['right_flank_100bp'][:flank_bp])
        draw_intron(ax, i, ex['seq'], r['left_flank_100bp'], r['right_flank_100bp'],
                    'TTAGGG', char_w=char_w, flank_bp=flank_bp, fontsize=5.5, tsd_len=tsd)
        mag = r.get('mag', '').replace('TARA_', '').split('_MAG_')[0]
        ax.text(-22, i+0.5, f"{mag} {len(ex['seq'])}bp", ha='left', va='center', fontsize=7)

    for i, ex in enumerate(eim_div):
        r = ex['rec']
        tsd = find_tsd(r['left_flank_100bp'][-flank_bp:], r['right_flank_100bp'][:flank_bp])
        draw_intron(ax, i, ex['seq'], r['left_flank_100bp'], r['right_flank_100bp'],
                    'TTTAGGG', char_w=char_w, flank_bp=flank_bp, fontsize=5.5, tsd_len=tsd, x0=eim_x)
        sp = r.get('species', '').replace('Eimeria_', 'E.')
        ax.text(eim_x-22, i+0.5, f"{sp[:14]} {len(ex['seq'])}bp", ha='left', va='center', fontsize=7)

    add_legend_simple(ax, 0, n_rows+0.5, fontsize=8)
    plt.savefig('real_telotrons/fig_compare.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("fig_compare.png")


def fig_shared():
    pairs_path = Path('real_telotrons/shared_linker_pairs.json')
    if not pairs_path.exists(): return
    pairs = json.load(open(pairs_path))

    contigs = defaultdict(dict)
    for acc in set(p['acc'] for p in pairs):
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen: continue
        cur, parts = None, []
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contigs[acc][cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]
                    parts = []
                else: parts.append(line.strip().upper())
            if cur: contigs[acc][cur] = ''.join(parts)

    fig, axes = plt.subplots(len(pairs), 1, figsize=(22, len(pairs)*3.5))
    if len(pairs) == 1: axes = [axes]

    for idx, p in enumerate(pairs):
        ax = axes[idx]; ax.axis('off')
        t1 = p['src_telotron']; t2 = p['paralog_telotron']
        link = t1['intron_seq'][t1['linker_offset']:t1['linker_end_offset']]
        seq2 = t2['intron_seq']; L = len(link)
        best_pos, best_match = -1, -1
        for i in range(len(seq2)-L+1):
            m = sum(a==b for a, b in zip(seq2[i:i+L], link))
            if m > best_match: best_match = m; best_pos = i

        cs = contigs.get(p['acc'], {})
        cs1 = cs.get(t1['contig'], '')
        left1 = cs1[max(0, t1['start']-100):t1['start']]
        right1 = cs1[t1['end']:t1['end']+100]
        cs2 = cs.get(t2['contig'], '')
        left2 = cs2[max(0, t2['start']-100):t2['start']]
        right2 = cs2[t2['end']:t2['end']+100]

        same = t1['contig'] == t2['contig']
        gap = abs(t2['start']-t1['end']) if same and t2['start']>t1['end'] else (abs(t1['start']-t2['end']) if same else None)
        gap_str = f"{gap}bp apart" if gap else "different contig"

        char_w = 1.0; flank_bp = 50
        max_w = max(len(t1['intron_seq']), len(t2['intron_seq'])) + 100
        ax.set_xlim(-25, max_w+5); ax.set_ylim(7, -2)

        ax.text(0, -1, f"{p['acc']}  {t1['contig']}:{t1['start']}-{t1['end']} ↔ "
                       f"{t2['contig']}:{t2['start']}-{t2['end']}  ({gap_str}, {best_match}/{L} matches)",
                fontsize=10, fontweight='bold')

        h1 = {pos: 'shared' for pos in range(t1['linker_offset'], t1['linker_end_offset'])}
        ranges_1 = [(t1['linker_offset'], t1['linker_end_offset'], 'shared')]
        draw_intron(ax, 1, t1['intron_seq'], left1, right1, 'TTTAGGG',
                    char_w=char_w, flank_bp=flank_bp, fontsize=5.5, highlights=ranges_1)

        ranges_2 = [(best_pos, best_pos+L, 'shared')] if best_pos >= 0 else []
        draw_intron(ax, 4, t2['intron_seq'], left2, right2, 'TTTAGGG',
                    char_w=char_w, flank_bp=flank_bp, fontsize=5.5, highlights=ranges_2)

    plt.tight_layout()
    plt.savefig('real_telotrons/fig_shared.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("fig_shared.png")


def fig_linkers(n_per=8):
    tara = load_tara()
    eim, _ = load_eimeria_with_flanks()

    pure_eim, mosaic_eim, pure_tara, mosaic_tara = [], [], [], []
    for r in eim:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, link = classify_orientation(seq, 'TTTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if not (20 <= len(linker) <= 80): continue
        (mosaic_eim if has_telo_kmer(linker, 'TTTAGGG') else pure_eim).append({
            'rec': r, 'seq': seq, 'link': link, 'family': 'TTTAGGG'})

    for r in tara:
        seq = r.get('intron_seq', '').upper()
        if not seq: continue
        cls, link = classify_orientation(seq, 'TTAGGG')
        if not link: continue
        ls, le = link
        linker = seq[ls:le]
        if not (20 <= len(linker) <= 80): continue
        (mosaic_tara if has_telo_kmer(linker, 'TTAGGG') else pure_tara).append({
            'rec': r, 'seq': seq, 'link': link, 'family': 'TTAGGG'})

    for L in [pure_eim, mosaic_eim, pure_tara, mosaic_tara]: random.shuffle(L)
    pure_eim, mosaic_eim = pure_eim[:n_per], mosaic_eim[:n_per]
    pure_tara, mosaic_tara = pure_tara[:n_per], mosaic_tara[:n_per]

    sections = [
        ('Tara mosaic linker (91%)', mosaic_tara),
        ('Tara pure linker (8.6%)', pure_tara),
        ('Eimeria mosaic linker', mosaic_eim),
        ('Eimeria pure linker (26%)', pure_eim),
    ]
    total = sum(len(s[1])+1 for s in sections)
    max_intron = max((len(x['seq']) for s in sections for x in s[1]), default=200)

    fig, ax = plt.subplots(figsize=(20, max(10, total*0.45)))
    ax.set_xlim(-30, max_intron+50); ax.set_ylim(total+1, -1)
    ax.axis('off')

    cur_y = 0
    for title, exs in sections:
        ax.text(-30, cur_y+0.5, title, ha='left', va='center', fontsize=10, fontweight='bold')
        cur_y += 1
        for ex in exs:
            seq = ex['seq']; ls, le = ex['link']
            cls = classify_bases(seq, ex['family'])
            render_row(ax, cur_y, seq, cls, x_start=0, char_width=1.0, fontsize=5.5)
            ax.plot([ls, ls], [cur_y, cur_y+1], color='green', linewidth=1.2, linestyle='--')
            ax.plot([le, le], [cur_y, cur_y+1], color='green', linewidth=1.2, linestyle='--')
            ax.text(-30, cur_y+0.5, f"{len(seq)}bp link={le-ls}",
                    ha='left', va='center', fontsize=7)
            ax.text(max_intron+5, cur_y+0.5, seq[ls:le][:50],
                    ha='left', va='center', fontsize=6.5, family='monospace', color='darkgreen')
            cur_y += 1
        cur_y += 0.5

    plt.savefig('real_telotrons/fig_linkers.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("fig_linkers.png")


def fig_tandem():
    eim, _ = load_eimeria_with_flanks()
    by_contig = defaultdict(list)
    for r in eim:
        seq = r.get('intron_seq', '')
        if not seq: continue
        try:
            by_contig[(r['acc'], r['contig'])].append({
                'start': int(r['start']), 'end': int(r['end'])
            })
        except: pass
    interesting = []
    for k, telos in by_contig.items():
        if not (3 <= len(telos) <= 5): continue
        telos.sort(key=lambda x: x['start'])
        gaps = [telos[i+1]['start']-telos[i]['end'] for i in range(len(telos)-1)]
        if all(0 <= g <= 1000 for g in gaps):
            interesting.append((k, telos, sum(gaps)))
    interesting.sort(key=lambda x: x[2])
    selected = interesting[:5]
    if not selected: return

    fig, axes = plt.subplots(len(selected), 1, figsize=(22, len(selected)*3))
    if len(selected) == 1: axes = [axes]
    for idx, ((acc, contig), telos, _) in enumerate(selected):
        ax = axes[idx]; ax.axis('off')
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        cs = ''
        if gen:
            cur, parts = None, []
            with open(gen) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur == contig: cs = ''.join(parts); break
                        cur = line[1:].strip().split()[0]; parts = []
                    elif cur == contig:
                        parts.append(line.strip().upper())
                if cur == contig and not cs: cs = ''.join(parts)
        if not cs: continue

        vis_s = max(0, telos[0]['start']-50)
        vis_e = min(len(cs), telos[-1]['end']+50)
        full_seq = cs[vis_s:vis_e]
        full_cls = ['flank']*len(full_seq)
        for t in telos:
            ts, te = t['start']-vis_s, t['end']-vis_s
            seq_intron = full_seq[ts:te]
            ic = classify_bases(seq_intron, 'TTTAGGG')
            for i, c in enumerate(ic):
                if ts+i < len(full_cls): full_cls[ts+i] = c

        ax.set_xlim(-15, len(full_seq)+5); ax.set_ylim(3, -1.5)
        ax.text(0, -1, f"{contig}  ({len(telos)} telotrons in {vis_e-vis_s}bp)",
                fontsize=9, fontweight='bold')
        render_row(ax, 0, full_seq, full_cls, x_start=0, char_width=1.0, fontsize=5.5)
        for t in telos:
            ts, te = t['start']-vis_s, t['end']-vis_s
            ax.add_patch(Rectangle((ts-0.3, 0), 0.6, 1.0, facecolor='black'))
            ax.add_patch(Rectangle((te-0.3, 0), 0.6, 1.0, facecolor='black'))
            ax.text((ts+te)/2, 1.5, f"{t['end']-t['start']}bp",
                    ha='center', va='center', fontsize=7)
    plt.tight_layout()
    plt.savefig('real_telotrons/fig_tandem.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("fig_tandem.png")


def fig_origins():
    """Linker→origin alignments with substantial flanking context (200bp each side)."""
    pairs = json.load(open('real_telotrons/linker_origin_pairs.json'))

    # Sort: shared-linker first, then same-contig, then cross
    def sort_key(p):
        if p.get('origin_in_telotron'): return (0, -p['pident'])
        if p['source_contig'] == p['origin_contig']: return (1, -p['pident'])
        return (2, -p['pident'])
    pairs.sort(key=sort_key)

    # Need: source telotron seq + flanks, origin sequence with 200bp flanks
    contigs = defaultdict(dict)
    for acc in set(p['acc'] for p in pairs):
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if not gen: continue
        cur, parts = None, []
        with open(gen) as f:
            for line in f:
                if line.startswith('>'):
                    if cur: contigs[acc][cur] = ''.join(parts)
                    cur = line[1:].strip().split()[0]; parts = []
                else: parts.append(line.strip().upper())
            if cur: contigs[acc][cur] = ''.join(parts)

    # For each pair build the origin context with 200bp flanks
    FLANK = 200
    for p in pairs:
        cs = contigs.get(p['acc'], {}).get(p['origin_contig'], '')
        os_, oe_ = p['origin_start'] - 1, p['origin_end']  # 1-based to 0-based
        flank_l = min(FLANK, os_)
        flank_r = min(FLANK, len(cs) - oe_)
        p['_origin_full'] = cs[os_-flank_l:oe_+flank_r]
        p['_origin_hit_start'] = flank_l
        p['_origin_hit_end'] = flank_l + (oe_ - os_)
        p['_origin_left_flank'] = flank_l
        p['_origin_right_flank'] = flank_r

        # source: telotron sequence + 200bp flanks
        cs_src = contigs.get(p['acc'], {}).get(p['source_contig'], '')
        src_telo_s = p.get('src_telotron_start')
        src_telo_e = p.get('src_telotron_end')
        if src_telo_s and src_telo_e and cs_src:
            sf_l = min(FLANK, src_telo_s)
            sf_r = min(FLANK, len(cs_src) - src_telo_e)
            p['_src_full'] = cs_src[src_telo_s-sf_l:src_telo_e+sf_r]
            p['_src_telo_start'] = sf_l
            p['_src_telo_end'] = sf_l + (src_telo_e - src_telo_s)
            p['_src_linker_offset'] = (p['source_linker_start'] - src_telo_s) + sf_l
            p['_src_linker_end'] = (p['source_linker_end'] - src_telo_s) + sf_l
        else:
            p['_src_full'] = ''
            p['_src_telo_start'] = 0
            p['_src_telo_end'] = 0
            p['_src_linker_offset'] = 0
            p['_src_linker_end'] = 0

    # Render in batches of 6
    batch_size = 6
    for batch_idx, b_start in enumerate(range(0, len(pairs), batch_size)):
        batch = pairs[b_start:b_start+batch_size]
        # Layout: each pair gets 5 row units (header + source + spacer + origin + spacer)
        PH = 5
        max_w = max(max(len(p['_src_full']), len(p['_origin_full'])) for p in batch)
        max_w = max(max_w, 600)

        fig, ax = plt.subplots(figsize=(max(18, max_w*0.04+5), len(batch)*PH*0.42))
        ax.set_xlim(-30, max_w+5); ax.set_ylim(len(batch)*PH+1, -2)
        ax.axis('off')

        for idx, p in enumerate(batch):
            y0 = idx * PH
            same = p['source_contig'] == p['origin_contig']
            if p.get('origin_in_telotron'):
                cls_lbl = 'shared-linker'
            elif same: cls_lbl = 'same-contig'
            else: cls_lbl = 'cross-contig'

            if same:
                src_mid = (p['source_linker_start'] + p['source_linker_end'])//2
                org_mid = (p['origin_start'] + p['origin_end'])//2
                dist_str = f"{abs(org_mid-src_mid)}bp apart"
            else:
                dist_str = f"different contig"

            # Header — concise
            global_idx = b_start + idx + 1
            ax.text(-30, y0+0.5,
                    f"#{global_idx}  {p['source_contig']}:{p['source_linker_start']}-{p['source_linker_end']} → "
                    f"{p['origin_contig']}:{p['origin_start']}-{p['origin_end']}  "
                    f"[{cls_lbl}, {dist_str}, {p['pident']:.0f}% over {p['aln_length']}bp]",
                    ha='left', va='center', fontsize=8.5, fontweight='bold')

            # Row 1: source telotron with linker highlighted
            if p['_src_full']:
                src_seq = p['_src_full']
                cls = classify_bases(src_seq, 'TTTAGGG')
                # Make non-telotron flanks display as 'flank'
                for i in range(min(p['_src_telo_start'], len(cls))):
                    cls[i] = 'flank'
                for i in range(p['_src_telo_end'], len(cls)):
                    cls[i] = 'flank'
                # Highlight linker
                hl = [(p['_src_linker_offset'], p['_src_linker_end'], 'shared')]
                # Splice site bars
                ax.add_patch(Rectangle((p['_src_telo_start']-0.3, y0+1.5),
                                        0.6, 1, facecolor='black'))
                ax.add_patch(Rectangle((p['_src_telo_end']-0.3, y0+1.5),
                                        0.6, 1, facecolor='black'))
                render_row(ax, y0+1.5, src_seq, cls, x_start=0, char_width=1.0,
                           highlight_ranges=hl, fontsize=5.5)
                ax.text(-30, y0+2, "source", ha='left', va='center', fontsize=7, color='gray')
            else:
                ax.text(0, y0+2, "[source telotron sequence unavailable]", fontsize=7, color='red')

            # Row 2: origin with 200bp flanks, hit highlighted
            origin_seq = p['_origin_full']
            cls2 = classify_bases(origin_seq, 'TTTAGGG') if p.get('origin_in_telotron') \
                   else ['flank']*len(origin_seq)
            hl2 = [(p['_origin_hit_start'], p['_origin_hit_end'], 'shared')]
            render_row(ax, y0+3.5, origin_seq, cls2, x_start=0, char_width=1.0,
                       highlight_ranges=hl2, fontsize=5.5)
            ax.text(-30, y0+4, "origin", ha='left', va='center', fontsize=7, color='gray')

            # Mark contig edges if very close
            if p['_origin_left_flank'] < 50:
                ax.plot([0-0.5, 0-0.5], [y0+3.5, y0+4.5],
                        color='red', linewidth=2)
                ax.text(-3, y0+4, "5'", ha='right', va='center', fontsize=7, color='red')
            if p['_origin_right_flank'] < 50:
                ax.plot([len(origin_seq)+0.5, len(origin_seq)+0.5], [y0+3.5, y0+4.5],
                        color='red', linewidth=2)
                ax.text(len(origin_seq)+3, y0+4, "3'", fontsize=7, color='red')

        out = f'real_telotrons/fig_origins_batch{batch_idx+1}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(out)


if __name__ == '__main__':
    fig_classes()
    fig_compare()
    fig_shared()
    fig_linkers()
    fig_tandem()
    fig_origins()
