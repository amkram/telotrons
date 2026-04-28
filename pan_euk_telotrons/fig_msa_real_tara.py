#!/usr/bin/env python3
"""Real MSA per Tara MAG. Same logic as fig_msa_real.py but for Tara."""
import csv, json, re, subprocess, tempfile, os
from pathlib import Path
from collections import defaultdict, Counter
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


def parse_gff_full(gff_path):
    by_contig = defaultdict(lambda: {'genes': [], 'exons': defaultdict(list),
                                      'introns': []})
    mrna_to_gene = {}
    gene_info = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 9: continue
            try: s, e = int(p[3]), int(p[4])
            except: continue
            attrs = dict(re.findall(r'(\w+)=([^;]+)', p[8]))
            ftype = p[2]; contig = p[0]
            if ftype == 'gene':
                gid = attrs.get('ID', '')
                name = attrs.get('Name', attrs.get('gene', gid))
                gene_info[gid] = {'name': name, 'start': s, 'end': e, 'contig': contig}
                by_contig[contig]['genes'].append({'gid': gid, 'name': name, 'start': s, 'end': e})
            elif ftype in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                par = attrs.get('Parent', '')
                mrna_to_gene[mid] = par
            elif ftype in ('exon', 'CDS'):
                par = attrs.get('Parent', '')
                by_contig[contig]['exons'][par].append((s, e, ftype))
    for c, info in by_contig.items():
        for mid, exs in info['exons'].items():
            exs_sorted = sorted(set((s, e) for s, e, t in exs if t == 'exon'))
            if not exs_sorted:
                exs_sorted = sorted(set((s, e) for s, e, t in exs if t == 'CDS'))
            if len(exs_sorted) < 2: continue
            for i in range(len(exs_sorted) - 1):
                istart = exs_sorted[i][1] + 1
                iend = exs_sorted[i+1][0] - 1
                if iend > istart:
                    gid = mrna_to_gene.get(mid, '')
                    gname = gene_info.get(gid, {}).get('name', '?')
                    info['introns'].append((istart, iend, gname))
    return by_contig


def annotate_position(by_contig, contig, pos, in_telotron_func=None):
    info = by_contig.get(contig, {})
    in_intron = None
    for s, e, gn in info.get('introns', []):
        if s <= pos <= e:
            in_intron = (s, e, gn); break
    in_gene = None
    for g in info.get('genes', []):
        if g['start'] <= pos <= g['end']:
            in_gene = g; break
    annot = []
    if in_telotron_func:
        t = in_telotron_func(contig, pos)
        if t: annot.append(f"in_telotron[{t[0]}-{t[1]}]")
    if in_intron:
        annot.append(f"in_intron({in_intron[2]})")
    elif in_gene:
        annot.append(f"in_exon({in_gene['name']})")
    else:
        annot.append("intergenic")
    return annot


def run_mafft(seqs):
    if len(seqs) == 1: return seqs
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for name, s in seqs:
            f.write(f">{name}\n{s}\n")
        in_path = f.name
    try:
        result = subprocess.run(['mafft', '--auto', '--quiet', in_path],
                                 capture_output=True, text=True, timeout=120)
        if result.returncode != 0: return seqs
        out = []
        cur = None; parts = []
        for line in result.stdout.split('\n'):
            if line.startswith('>'):
                if cur: out.append((cur, ''.join(parts).upper()))
                cur = line[1:].strip(); parts = []
            else: parts.append(line.strip())
        if cur: out.append((cur, ''.join(parts).upper()))
        return out
    finally:
        os.unlink(in_path)


def main():
    groups = json.load(open('real_telotrons/tara_linker_groups.json'))

    # Group by mag, sort by group size, cap to top 15
    by_mag = defaultdict(list)
    for g in groups:
        by_mag[g['mag']].append(g)
    LIMIT = 15
    for m in by_mag:
        by_mag[m].sort(key=lambda x: -len(x['members']))
        by_mag[m] = by_mag[m][:LIMIT]

    GENOME_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/contigs_individual')
    GFF_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/gff_individual/GFF')

    contigs = defaultdict(dict)
    gff_data = {}
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
            gff_data[mag] = parse_gff_full(gff)

    telo_lookup = {}
    contig_telo = defaultdict(list)
    with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            try:
                key = (row['mag'], row['contig'], int(row['start']), int(row['end']))
                telo_lookup[key] = row.get('intron_seq', '').upper()
                contig_telo[(row['mag'], row['contig'])].append((int(row['start']), int(row['end'])))
            except: pass

    LINKER_FLANK = 80

    for mag, gs in by_mag.items():
        if not gs: continue
        def in_telo(contig, pos):
            for s, e in contig_telo.get((mag, contig), []):
                if s <= pos <= e: return (s, e)
            return None

        rendered = []
        for g in gs:
            members = [tuple(m) for m in g['members']]
            first = members[0]
            seq0 = telo_lookup.get(first, '')
            if not seq0: continue

            # Canonical linker = longest non-telomeric run (TTAGGG family)
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
                        if i - s > best[1] - best[0]: best = (s, i)
                        s = None
            if s is not None and n - s > best[1] - best[0]: best = (s, n)
            if best[1] - best[0] < 15: continue
            canon = seq0[best[0]:best[1]]

            seqs_to_align = []
            row_meta = []

            # ORIGIN row
            if g['type'] == 'locus':
                org_contig = g['origin_contig']
                cs = contigs.get(mag, {}).get(org_contig, '')
                bin_s, bin_e = g['origin_bin_start'], g['origin_bin_end']
                region = cs[bin_s:bin_e] if cs else ''
                if region:
                    pos, _, strand = find_best_match(canon, region)
                    if pos >= 0:
                        gpos_s = bin_s + pos
                        gpos_e = gpos_s + len(canon)
                        if strand == 'rev':
                            ext_l = min(LINKER_FLANK, len(cs) - gpos_e)
                            ext_r = min(LINKER_FLANK, gpos_s)
                            ext_seq = revcomp(cs[gpos_s - ext_r : gpos_e + ext_l])
                            link_off = ext_l
                        else:
                            ext_l = min(LINKER_FLANK, gpos_s)
                            ext_r = min(LINKER_FLANK, len(cs) - gpos_e)
                            ext_seq = cs[gpos_s - ext_l : gpos_e + ext_r]
                            link_off = ext_l
                        annot = annotate_position(gff_data.get(mag, {}), org_contig, gpos_s, in_telo)
                        contig_len = len(cs)
                        d_to_end = min(gpos_s, contig_len - gpos_e)
                        if d_to_end < 5000:
                            annot.append(f"subtelo({d_to_end:,}bp)")
                        else:
                            annot.append(f"end={d_to_end:,}bp")
                        seqs_to_align.append((f"ORIGIN_{gpos_s}", ext_seq))
                        row_meta.append({
                            'label': f"ORIGIN  {org_contig[-25:]}:{gpos_s}-{gpos_e}{(' (rc)' if strand=='rev' else '')}",
                            'kind': 'ORIGIN',
                            'annot': '; '.join(annot),
                            'link_in_unaligned': link_off,
                        })

            if g['type'] == 'locus' and not seqs_to_align: continue

            for m in members:
                _, mc, ts, te = m
                seq = telo_lookup.get((mag, mc, ts, te), '')
                if not seq: continue
                pos, _, strand = find_best_match(canon, seq)
                if pos < 0: continue
                if strand == 'rev':
                    s_in_telo = len(seq) - (pos + len(canon))
                    seq_ori = revcomp(seq)
                    pos_ori = s_in_telo
                else:
                    seq_ori = seq; pos_ori = pos
                ext_l = min(LINKER_FLANK, pos_ori)
                ext_r = min(LINKER_FLANK, len(seq_ori) - (pos_ori + len(canon)))
                ext_seq = seq_ori[pos_ori - ext_l : pos_ori + len(canon) + ext_r]
                link_off = ext_l
                annot = annotate_position(gff_data.get(mag, {}), mc, ts, in_telo)
                seqs_to_align.append((f"telo_{mc[-12:]}_{ts}{'rc' if strand=='rev' else ''}", ext_seq))
                row_meta.append({
                    'label': f"telo  {mc[-25:]}:{ts}-{te}{(' (rc)' if strand=='rev' else '')}",
                    'kind': 'TELOTRON',
                    'annot': '; '.join(annot),
                    'link_in_unaligned': link_off,
                })

            if len(seqs_to_align) < 2: continue
            aligned = run_mafft(seqs_to_align)
            if len(aligned) != len(seqs_to_align): continue

            # Find aligned linker positions in the reference (first row)
            ref_aligned = aligned[0][1]
            ref_unaligned_to_aligned = []
            for i, c in enumerate(ref_aligned):
                if c != '-':
                    ref_unaligned_to_aligned.append(i)
            ref_link_off = row_meta[0]['link_in_unaligned']
            if ref_link_off < len(ref_unaligned_to_aligned):
                aligned_link_start = ref_unaligned_to_aligned[ref_link_off]
                ref_link_end_pos = ref_link_off + len(canon)
                if ref_link_end_pos <= len(ref_unaligned_to_aligned):
                    aligned_link_end = ref_unaligned_to_aligned[ref_link_end_pos - 1] + 1
                else:
                    aligned_link_end = len(ref_aligned)
            else:
                aligned_link_start = 0; aligned_link_end = len(canon)

            rendered.append({
                'type': g['type'], 'aligned': aligned, 'meta': row_meta,
                'canon': canon,
                'alink_s': aligned_link_start, 'alink_e': aligned_link_end,
            })

        if not rendered: continue
        total_rows = sum(2 + len(g['aligned']) for g in rendered)
        max_w = max(len(g['aligned'][0][1]) for g in rendered)
        max_w = min(max_w, 2500)
        fig, ax = plt.subplots(figsize=(max(20, max_w*0.04), max(8, total_rows*0.45)))
        ax.set_xlim(-50, max_w+5); ax.set_ylim(total_rows+1, -2)
        ax.axis('off')
        plt.suptitle(f"{mag}: linker MSA per group (top {len(rendered)} groups)",
                     fontsize=12, fontweight='bold')

        cur_y = 0
        for g_idx, g in enumerate(rendered):
            n_telo = sum(1 for m in g['meta'] if m['kind'] == 'TELOTRON')
            has_origin = any(m['kind'] == 'ORIGIN' for m in g['meta'])
            if has_origin:
                hdr = f"#{g_idx+1}  {n_telo} telotron(s) + origin (linker {len(g['canon'])}bp)"
            else:
                hdr = f"#{g_idx+1}  Shared linker among {n_telo} telotrons ({len(g['canon'])}bp)"
            ax.text(-50, cur_y+0.5, hdr, ha='left', va='center', fontsize=10, fontweight='bold')
            cur_y += 1

            ref_seq = g['aligned'][0][1]
            consensus = []
            for col in range(len(ref_seq)):
                col_chars = [seq[col] for _, seq in g['aligned']
                             if col < len(seq) and seq[col] != '-']
                if col_chars:
                    consensus.append(Counter(col_chars).most_common(1)[0][0])
                else: consensus.append('-')

            for (label, seq), meta in zip(g['aligned'], g['meta']):
                cls = ['flank']*len(seq)
                bg_overrides = []
                for j, c in enumerate(seq):
                    if c == '-':
                        bg_overrides.append((j, j+1, (1, 1, 1)))
                    elif c != consensus[j] and consensus[j] != '-':
                        bg_overrides.append((j, j+1, (1.0, 0.7, 0.7)))
                hl = [(g['alink_s'], g['alink_e'], 'shared')]
                kind_color = 'darkred' if meta['kind'] == 'ORIGIN' else 'black'
                fw = 'bold' if meta['kind'] == 'ORIGIN' else 'normal'
                ax.text(-50, cur_y+0.5, meta['label'][:50],
                        ha='left', va='center', fontsize=7, color=kind_color, fontweight=fw)
                render_row(ax, cur_y, seq, cls, x_start=0, char_width=1.0,
                           highlight_ranges=hl, bg_override=bg_overrides, fontsize=5.5)
                ax.text(max_w+3, cur_y+0.5, meta['annot'],
                        ha='left', va='center', fontsize=6.5, color='gray')
                cur_y += 1
            cur_y += 1

        out = f'real_telotrons/fig_msa_real_{mag}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Saved: {out}  ({len(rendered)} groups)")


if __name__ == '__main__':
    main()
