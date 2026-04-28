#!/usr/bin/env python3
"""Per-MAG linker→origin MSA figures for Tara haptophyte telotrons.
Mirrors fig_msa_per_species.py logic but for Tara MAGs."""
import csv, json, re
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from sequence_viewer import classify_bases, render_row


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def find_best_match(query, target):
    L = len(query); n = len(target)
    if L > n: return -1, -1, 'fwd'
    best = (-1, -1, 'fwd')
    for i in range(n - L + 1):
        m = sum(a == b for a, b in zip(query, target[i:i+L]))
        if m > best[1]: best = (i, m, 'fwd')
    rc = revcomp(query)
    for i in range(n - L + 1):
        m = sum(a == b for a, b in zip(rc, target[i:i+L]))
        if m > best[1]: best = (i, m, 'rev')
    return best


def parse_gff_genes(gff_path):
    genes = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 9: continue
            if p[2] != 'gene' and p[2] != 'mRNA': continue
            attrs = dict(re.findall(r'(\w+)=([^;]+)', p[8]))
            name = attrs.get('Name', attrs.get('ID', '?'))
            try:
                genes[p[0]].append({
                    'start': int(p[3]), 'end': int(p[4]),
                    'name': name,
                })
            except: pass
    return genes


def find_gene(gene_index, contig, pos):
    for g in gene_index.get(contig, []):
        if g['start'] <= pos <= g['end']:
            return g
    return None


def main():
    groups = json.load(open('real_telotrons/tara_linker_groups.json'))

    # Group by MAG, sort by group size (members)
    by_mag = defaultdict(list)
    for g in groups:
        by_mag[g['mag']].append(g)
    for m in by_mag:
        by_mag[m].sort(key=lambda x: -len(x['members']))

    # Cap to top 15 per MAG
    LIMIT = 15
    for m in by_mag:
        by_mag[m] = by_mag[m][:LIMIT]

    # Tara genome paths
    GENOME_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/contigs_individual')
    GFF_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/gff_individual/GFF')

    contigs = defaultdict(dict)
    gene_index = {}
    for mag in by_mag:
        fa = GENOME_DIR / f"{mag}.fa"
        if fa.exists():
            cur, parts = None, []
            with open(fa) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur: contigs[mag][cur] = ''.join(parts)
                        cur = line[1:].strip().split()[0]; parts = []
                    else: parts.append(line.strip().upper())
                if cur: contigs[mag][cur] = ''.join(parts)
        gff = GFF_DIR / f"{mag}.gmove.gff"
        if gff.exists():
            gene_index[mag] = parse_gff_genes(gff)

    # Telotron lookup (use _tara_oceans_telotrons.tsv)
    telo_lookup = {}
    contig_telotron_intervals = defaultdict(list)
    with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            try:
                key = (row['mag'], row['contig'], int(row['start']), int(row['end']))
                telo_lookup[key] = row.get('intron_seq', '').upper()
                contig_telotron_intervals[(row['mag'], row['contig'])].append(
                    (int(row['start']), int(row['end'])))
            except: pass

    def in_telotron(mag, contig, pos):
        for s, e in contig_telotron_intervals.get((mag, contig), []):
            if s <= pos <= e: return (s, e)
        return None

    FLANK = 200

    for mag, gs in by_mag.items():
        if not gs: continue
        rendered = []
        for g in gs:
            members = [tuple(m) for m in g['members']]
            first = members[0]
            seq0 = telo_lookup.get(first, '')
            if not seq0: continue

            # Find longest non-telomeric run = canonical linker
            n = len(seq0); cov = bytearray(n)
            for k in {'TTAGGG','TAGGGT','AGGGTT','GGGTTA','GGTTAG','GTTAGG',
                      'CCCTAA','CCTAAC','CTAACC','TAACCC','AACCCT','ACCCTA'}:
                i = 0
                while True:
                    idx = seq0.find(k, i)
                    if idx == -1: break
                    for j in range(idx, idx+6):
                        if j < n: cov[j] = 1
                    i = idx + 1
            best = (0, 0); s = None
            for i in range(n):
                if cov[i] == 0:
                    if s is None: s = i
                else:
                    if s is not None:
                        if i - s > best[1]-best[0]: best = (s, i)
                        s = None
            if s is not None and n - s > best[1]-best[0]: best = (s, n)
            if best[1] - best[0] < 15: continue
            canon = seq0[best[0]:best[1]]

            # Build rows
            rows = []

            origin_row = None
            if g['type'] == 'locus':
                org_contig = g['origin_contig']
                cs = contigs.get(mag, {}).get(org_contig, '')
                bin_s, bin_e = g['origin_bin_start'], g['origin_bin_end']
                region = cs[bin_s:bin_e] if cs else ''
                if region:
                    pos, m, strand = find_best_match(canon, region)
                    if pos >= 0:
                        gpos_s = bin_s + pos
                        gpos_e = gpos_s + len(canon)
                        if strand == 'rev':
                            ext_l = min(FLANK, len(cs) - gpos_e)
                            ext_r = min(FLANK, gpos_s)
                            ext_seq = revcomp(cs[gpos_s - ext_r : gpos_e + ext_l])
                            link_off = ext_l
                        else:
                            ext_l = min(FLANK, gpos_s)
                            ext_r = min(FLANK, len(cs) - gpos_e)
                            ext_seq = cs[gpos_s - ext_l : gpos_e + ext_r]
                            link_off = ext_l
                        gene = find_gene(gene_index.get(mag, {}), org_contig, gpos_s)
                        in_telo = in_telotron(mag, org_contig, gpos_s)
                        contig_len = len(cs)
                        d_to_end = min(gpos_s, contig_len - gpos_e)
                        annot = []
                        if gene: annot.append(f"gene={gene['name'][:14]}")
                        if in_telo: annot.append(f"INSIDE telotron[{in_telo[0]}-{in_telo[1]}]")
                        else: annot.append("intergenic")
                        if d_to_end < 5000:
                            annot.append(f"subtelo({d_to_end:,}bp)")
                        else:
                            annot.append(f"internal(end={d_to_end:,}bp)")
                        origin_row = {
                            'kind': 'ORIGIN',
                            'label': f"{org_contig[-25:]}:{gpos_s}-{gpos_e}{('(rc)' if strand=='rev' else '')}",
                            'annot': '; '.join(annot),
                            'seq': ext_seq, 'link_offset': link_off,
                            'link_length': len(canon), 'is_telotron_seq': bool(in_telo),
                        }

            if origin_row is None and g['type'] == 'locus':
                continue

            member_rows = []
            for m in members:
                mac, mc, ts, te = m
                seq = telo_lookup.get((mac, mc, ts, te), '')
                cs = contigs.get(mag, {}).get(mc, '')
                if not seq or not cs: continue
                pos, n_match, strand = find_best_match(canon, seq)
                flank = 100
                left = cs[max(0, ts-flank):ts]
                right = cs[te:te+flank]
                if strand == 'rev':
                    full = revcomp(right) + revcomp(seq) + revcomp(left)
                    link_pos = len(revcomp(right)) + (len(seq) - (pos + len(canon)))
                    intron_off = len(revcomp(right))
                else:
                    full = left + seq + right
                    link_pos = len(left) + pos
                    intron_off = len(left)
                annot = []
                gene = find_gene(gene_index.get(mag, {}), mc, ts)
                if gene: annot.append(f"gene={gene['name'][:14]}")
                annot.append(f"in_telotron[{ts}-{te}]")
                d_end = min(ts, len(cs)-te)
                if d_end < 5000: annot.append(f"subtelo({d_end:,}bp)")
                member_rows.append({
                    'kind': 'TELOTRON',
                    'label': f"{mc[-25:]}:{ts}-{te}{'(rc)' if strand=='rev' else ''}",
                    'annot': '; '.join(annot),
                    'seq': full,
                    'link_offset': link_pos,
                    'link_length': len(canon),
                    'intron_offset': intron_off,
                    'intron_length': len(seq),
                })

            if not member_rows: continue
            all_rows = ([origin_row] if origin_row else []) + member_rows
            max_left = max(r['link_offset'] for r in all_rows)
            max_right = max(len(r['seq']) - r['link_offset'] - r['link_length'] for r in all_rows)
            for r in all_rows:
                lp = max_left - r['link_offset']
                rp = (max_left + r['link_length'] + max_right) - (lp + len(r['seq']))
                r['_padded'] = '-'*lp + r['seq'] + '-'*rp
                r['_alink_s'] = max_left
                r['_alink_e'] = max_left + r['link_length']
                if r['kind'] == 'TELOTRON':
                    r['_aint_s'] = lp + r['intron_offset']
                    r['_aint_e'] = r['_aint_s'] + r['intron_length']

            rendered.append({
                'type': g['type'], 'rows': all_rows, 'canon': canon,
                'aligned_width': max_left + len(canon) + max_right,
            })

        if not rendered: continue
        total_rows = sum(2 + len(g['rows']) for g in rendered)
        max_w = min(max(g['aligned_width'] for g in rendered), 2500)
        fig_w = max(20, max_w * 0.04)
        fig_h = max(8, total_rows * 0.45)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(-50, max_w + 5)
        ax.set_ylim(total_rows + 1, -2)
        ax.axis('off')

        plt.suptitle(f"{mag}: linker→origin MSAs (top {len(rendered)} groups)",
                     fontsize=12, fontweight='bold')

        cur_y = 0
        for g_idx, g in enumerate(rendered):
            n_telo = sum(1 for r in g['rows'] if r['kind'] == 'TELOTRON')
            has_origin = any(r['kind'] == 'ORIGIN' for r in g['rows'])
            if has_origin:
                hdr = f"#{g_idx+1}  {n_telo} telotron(s) → origin (linker {len(g['canon'])}bp)"
            else:
                hdr = f"#{g_idx+1}  Shared linker between {n_telo} telotrons ({len(g['canon'])}bp)"
            ax.text(-50, cur_y+0.5, hdr, ha='left', va='center',
                    fontsize=9.5, fontweight='bold')
            cur_y += 1
            for r in g['rows']:
                seq = r['_padded']
                cls = ['flank']*len(seq)
                if r['kind'] == 'TELOTRON':
                    intron_seq = seq[r['_aint_s']:r['_aint_e']]
                    if '-' not in intron_seq:
                        ic = classify_bases(intron_seq, 'TTAGGG')
                        for j, c in enumerate(ic):
                            if r['_aint_s'] + j < len(cls):
                                cls[r['_aint_s'] + j] = c
                hl = [(r['_alink_s'], r['_alink_e'], 'shared')]
                bg_overrides = []
                for j, c in enumerate(seq):
                    if c == '-':
                        bg_overrides.append((j, j+1, (1,1,1)))
                kind_label = 'ORIGIN' if r['kind'] == 'ORIGIN' else 'telo'
                kind_color = 'darkred' if r['kind'] == 'ORIGIN' else 'black'
                fw = 'bold' if r['kind'] == 'ORIGIN' else 'normal'
                ax.text(-50, cur_y+0.5, f"{kind_label} {r['label']}",
                        ha='left', va='center', fontsize=7, color=kind_color, fontweight=fw)
                render_row(ax, cur_y, seq, cls, x_start=0, char_width=1.0,
                           highlight_ranges=hl, bg_override=bg_overrides, fontsize=5.5)
                ax.text(g['aligned_width']+3, cur_y+0.5, r['annot'],
                        ha='left', va='center', fontsize=6.5, color='gray')
                if r['kind'] == 'TELOTRON':
                    ax.add_patch(Rectangle((r['_aint_s']-0.4, cur_y), 0.4, 1.0, facecolor='black'))
                    ax.add_patch(Rectangle((r['_aint_e'], cur_y), 0.4, 1.0, facecolor='black'))
                cur_y += 1
            cur_y += 1

        out = f'real_telotrons/fig_msa_{mag}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Saved: {out}  ({len(rendered)} groups)")


if __name__ == '__main__':
    main()
