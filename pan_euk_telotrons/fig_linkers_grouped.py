#!/usr/bin/env python3
"""
Grouped linker figure: each panel shows one group of telotrons that share
the same linker, plus the originating genomic locus with substantial
flanking context.

Layout per group:
  Row 1+: each telotron in the group (full intron + 100bp flanks).
          Linker region highlighted in green.
  Row N:  the originating locus (intergenic/subtelomeric, or "in another
          telotron"). Shown with 200bp flanks. Aligned position highlighted.
"""
import csv, json, re
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from sequence_viewer import classify_bases, render_row, COLOR_BG


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def find_best_match(query, target):
    """Return (start, n_matches) for best ungapped placement of query in target.
    Tries both fwd and revcomp. Returns positive start for fwd, negative for rev."""
    L = len(query); n = len(target)
    if L > n: return None
    best = (-1, -1)  # start, matches
    for i in range(n - L + 1):
        m = sum(a == b for a, b in zip(query, target[i:i+L]))
        if m > best[1]: best = (i, m)
    rc = revcomp(query)
    rc_best = (-1, -1)
    for i in range(n - L + 1):
        m = sum(a == b for a, b in zip(rc, target[i:i+L]))
        if m > rc_best[1]: rc_best = (i, m)
    if rc_best[1] > best[1]:
        return rc_best[0], rc_best[1], 'rev'
    return best[0], best[1], 'fwd'


def main():
    groups_path = Path('real_telotrons/linker_groups.json')
    if not groups_path.exists():
        print("Run the grouping script first")
        return
    groups = json.load(open(groups_path))
    print(f"Loaded {len(groups)} groups")

    # Load contigs (per acc)
    accs = set()
    for g in groups:
        for m in g['members']:
            accs.add(m[0])
        if g['type'] == 'locus':
            accs.add(g['origin_acc'])

    contigs = defaultdict(dict)
    for acc in accs:
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
        print(f"  Loaded {acc}: {len(contigs[acc])} contigs")

    # Telotron lookup
    telo_lookup = {}  # (acc, contig, start, end) -> intron_seq
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                try:
                    key = (acc, row['contig'], int(row['start']), int(row['end']))
                    telo_lookup[key] = row.get('intron_seq', '').upper()
                except: pass

    # For each group, compute the canonical linker sequence (from first member)
    # and find its position in each member's intron and in the origin locus
    rendered_groups = []
    for g in groups:
        members = [tuple(m) for m in g['members']]
        # Get each telotron's intron seq + flanks
        member_recs = []
        for m in members:
            acc, contig, ts, te = m
            seq = telo_lookup.get((acc, contig, ts, te), '')
            if not seq: continue
            cs = contigs.get(acc, {}).get(contig, '')
            left = cs[max(0, ts-100):ts] if cs else ''
            right = cs[te:te+100] if cs else ''
            member_recs.append({
                'acc': acc, 'contig': contig, 'start': ts, 'end': te,
                'intron': seq, 'left': left, 'right': right,
            })
        if not member_recs: continue

        # Determine the canonical linker — for SHARED-LINKER groups, find
        # the linker region in the first telotron (intron has fwd/rev arrays
        # with non-repeat between them). For LOCUS groups, use the linker
        # region of the source telotron.

        # Try to find linker in first member: scan for non-repeat region
        first = member_recs[0]
        seq0 = first['intron']
        # Find longest non-telomeric stretch
        TT = {'TTTAGGG', 'TTAGGG', 'CCCTAAA', 'CCCTAA',
              'TAGGGTT', 'AGGGTTT', 'GGGTTTA', 'GGTTTAG', 'GTTTAGG', 'TTTAGGG',
              'AAACCCT', 'AACCCTA', 'ACCCTAA', 'CCCTAAA', 'CCTAAAC', 'CTAAACC'}
        # Build coverage
        n = len(seq0); cov = bytearray(n)
        for k in TT:
            i = 0
            while True:
                idx = seq0.find(k, i)
                if idx == -1: break
                for j in range(idx, idx+len(k)):
                    if j < n: cov[j] = 1
                i = idx + 1
        # Find longest run of 0
        best = (0, 0); s = None
        for i in range(n):
            if cov[i] == 0:
                if s is None: s = i
            else:
                if s is not None:
                    if i - s > best[1] - best[0]: best = (s, i)
                    s = None
        if s is not None and n - s > best[1] - best[0]: best = (s, n)
        link_start, link_end = best
        if link_end - link_start < 15: continue  # too short
        canonical_linker = seq0[link_start:link_end]

        # Find position of canonical linker in each member
        for mr in member_recs:
            pos, matches, strand = find_best_match(canonical_linker, mr['intron'])
            if pos < 0:
                mr['link_pos'] = (link_start, link_end)
                mr['link_matches'] = 0
                mr['link_strand'] = 'fwd'
            else:
                mr['link_pos'] = (pos, pos + len(canonical_linker))
                mr['link_matches'] = matches
                mr['link_strand'] = strand

        # Origin info
        origin_info = None
        if g['type'] == 'locus':
            org_contig = g['origin_contig']
            org_acc = g['origin_acc']
            cs = contigs.get(org_acc, {}).get(org_contig, '')
            # Use bin to bound region; refine to actual hit
            org_full_start = g['origin_bin_start']
            org_full_end = g['origin_bin_end']
            # Find the linker in the bin region
            region = cs[org_full_start:org_full_end] if cs else ''
            pos, matches, strand = find_best_match(canonical_linker, region) if region else (-1, 0, 'fwd')
            if pos >= 0:
                hit_genomic_start = org_full_start + pos
                hit_genomic_end = hit_genomic_start + len(canonical_linker)
                # Extend with 200bp flanks
                ext_l = min(200, hit_genomic_start)
                ext_r = min(200, len(cs) - hit_genomic_end) if cs else 0
                origin_info = {
                    'acc': org_acc, 'contig': org_contig,
                    'hit_start': hit_genomic_start, 'hit_end': hit_genomic_end,
                    'show_start': hit_genomic_start - ext_l,
                    'show_end': hit_genomic_end + ext_r,
                    'show_seq': cs[hit_genomic_start - ext_l : hit_genomic_end + ext_r] if cs else '',
                    'hit_offset': ext_l,
                    'hit_strand': strand,
                    'matches': matches,
                    'is_telotron': False,
                }

        rendered_groups.append({
            'type': g['type'],
            'members': member_recs,
            'canonical_linker': canonical_linker,
            'origin': origin_info,
        })

    print(f"\nRendering {len(rendered_groups)} groups")

    # Render in batches
    batch_size = 6
    for batch_idx, b_start in enumerate(range(0, len(rendered_groups), batch_size)):
        batch = rendered_groups[b_start:b_start+batch_size]
        # Compute layout
        # Each group needs: 1 header + len(members) + 1 origin + 1 spacer
        rows_per_group = lambda g: 2 + len(g['members']) + (1 if g['origin'] else 0)
        total_rows = sum(rows_per_group(g) for g in batch)
        max_w = 0
        for g in batch:
            for mr in g['members']:
                w = len(mr['left'][-100:]) + 1 + len(mr['intron']) + 1 + len(mr['right'][:100])
                max_w = max(max_w, w)
            if g['origin']:
                max_w = max(max_w, len(g['origin']['show_seq']))

        fig, ax = plt.subplots(figsize=(max(20, max_w*0.04+5), max(8, total_rows*0.45)))
        ax.set_xlim(-30, max_w+5); ax.set_ylim(total_rows+1, -1)
        ax.axis('off')

        cur_y = 0
        for g_idx, g in enumerate(batch):
            global_idx = b_start + g_idx + 1
            n_members = len(g['members'])

            # Header
            if g['type'] == 'shared':
                header = (f"#{global_idx}  Shared-linker group ({n_members} telotrons share "
                         f"{len(g['canonical_linker'])}bp linker)")
            else:
                origin = g['origin']
                if origin:
                    header = (f"#{global_idx}  {n_members} telotron(s) → genomic origin "
                             f"{origin['contig']}:{origin['hit_start']}-{origin['hit_end']}")
                else:
                    header = f"#{global_idx}  {n_members} telotron(s) (origin not localized)"
            ax.text(-30, cur_y+0.5, header, ha='left', va='center',
                    fontsize=10, fontweight='bold')
            cur_y += 1

            # Each telotron
            for mr in g['members']:
                lf = mr['left'][-100:]
                rf = mr['right'][:100]
                full_x = 0
                # Left flank
                if lf:
                    cls = ['flank']*len(lf)
                    render_row(ax, cur_y, lf, cls, x_start=full_x, char_width=1.0, fontsize=5.5)
                    full_x += len(lf)
                # 5' splice
                ax.add_patch(Rectangle((full_x, cur_y), 0.6, 1.0, facecolor='black'))
                full_x += 0.6
                # Intron with linker highlighted
                cls = classify_bases(mr['intron'], 'TTTAGGG')
                hl = [(mr['link_pos'][0], mr['link_pos'][1], 'shared')]
                render_row(ax, cur_y, mr['intron'], cls, x_start=full_x,
                           char_width=1.0, highlight_ranges=hl, fontsize=5.5)
                full_x += len(mr['intron'])
                # 3' splice
                ax.add_patch(Rectangle((full_x, cur_y), 0.6, 1.0, facecolor='black'))
                full_x += 0.6
                # Right flank
                if rf:
                    cls = ['flank']*len(rf)
                    render_row(ax, cur_y, rf, cls, x_start=full_x, char_width=1.0, fontsize=5.5)

                strand_mark = '←rc' if mr['link_strand'] == 'rev' else ''
                ax.text(-30, cur_y+0.5,
                        f"{mr['contig']}:{mr['start']}-{mr['end']} ({len(mr['intron'])}bp) "
                        f"link@{mr['link_pos'][0]}-{mr['link_pos'][1]} {strand_mark}",
                        ha='left', va='center', fontsize=7)
                cur_y += 1

            # Origin (if non-telotron locus)
            if g['origin']:
                origin = g['origin']
                cls = ['flank']*len(origin['show_seq'])
                hl = [(origin['hit_offset'], origin['hit_offset']+len(g['canonical_linker']), 'shared')]
                render_row(ax, cur_y, origin['show_seq'], cls, x_start=0,
                           char_width=1.0, highlight_ranges=hl, fontsize=5.5)
                strand_mark = '←rc' if origin['hit_strand'] == 'rev' else ''
                ax.text(-30, cur_y+0.5,
                        f"ORIGIN {origin['contig']}:{origin['show_start']}-{origin['show_end']} "
                        f"hit@{origin['hit_offset']}-{origin['hit_offset']+len(g['canonical_linker'])} "
                        f"{strand_mark}",
                        ha='left', va='center', fontsize=7, fontweight='bold')
                cur_y += 1
            cur_y += 1  # spacer

        out = f'real_telotrons/fig_linkers_grouped_batch{batch_idx+1}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  {out}")


if __name__ == '__main__':
    main()
