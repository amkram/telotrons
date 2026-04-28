#!/usr/bin/env python3
"""
Real MSA figure per species. For each linker group, run mafft to compute
a full multiple alignment of:
  1. Origin sequence (with flanks) — placed at TOP
  2. All linker sequences from group telotrons (with flanks)
For singleton groups (1 telotron + origin), MSA is just origin + that linker.
For shared-linker groups (no separate origin), MSA is the linkers themselves.
Show gaps/mismatches in the MSA. Annotate: gene context, exon/intron, telotron.
"""
import csv, json, re, subprocess, tempfile, os
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

from sequence_viewer import classify_bases, render_row, COLOR_BG


ACC_NAME = {
    'GCF_000499385.1': 'Eimeria_necatrix',
    'GCF_000499425.1': 'Eimeria_acervulina',
    'GCF_000499545.2': 'Eimeria_tenella',
    'GCF_000499605.1': 'Eimeria_maxima',
}


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
    """Parse GFF for genes, mRNAs, exons → introns."""
    by_contig = defaultdict(lambda: {'genes': [], 'exons': defaultdict(list),
                                      'mrna_contig': {}, 'mrna_strand': {}})
    mrna_to_gene = {}
    gene_info = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 9: continue
            ftype = p[2]
            contig = p[0]
            try:
                s, e = int(p[3]), int(p[4])
            except: continue
            attrs = dict(re.findall(r'(\w+)=([^;]+)', p[8]))
            if ftype == 'gene':
                gid = attrs.get('ID', '')
                name = attrs.get('Name', attrs.get('gene', gid))
                biotype = attrs.get('gene_biotype', '?')
                gene_info[gid] = {'name': name, 'biotype': biotype,
                                   'start': s, 'end': e, 'strand': p[6],
                                   'contig': contig}
                by_contig[contig]['genes'].append({'gid': gid, 'name': name,
                                                    'start': s, 'end': e})
            elif ftype in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                par = attrs.get('Parent', '')
                mrna_to_gene[mid] = par
                by_contig[contig]['mrna_contig'][mid] = contig
                by_contig[contig]['mrna_strand'][mid] = p[6]
            elif ftype in ('exon', 'CDS'):
                par = attrs.get('Parent', '')
                by_contig[contig]['exons'][par].append((s, e, ftype))
    # Compute introns per mRNA
    for c, info in by_contig.items():
        info['introns'] = []
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
    return by_contig, gene_info


def annotate_position(by_contig, contig, pos, in_telotron_func=None):
    """Return list of annotation strings for a position."""
    info = by_contig.get(contig, {})
    # In annotated intron?
    in_intron = None
    for s, e, gname in info.get('introns', []):
        if s <= pos <= e:
            in_intron = (s, e, gname); break
    # In gene?
    in_gene = None
    for g in info.get('genes', []):
        if g['start'] <= pos <= g['end']:
            in_gene = g; break
    annot = []
    if in_telotron_func:
        t = in_telotron_func(contig, pos)
        if t:
            annot.append(f"in_telotron[{t[0]}-{t[1]}]")
    if in_intron:
        annot.append(f"in_intron({in_intron[2]})")
    elif in_gene:
        # In gene but not in intron → exon (or untranslated region)
        annot.append(f"in_exon({in_gene['name']})")
    else:
        annot.append("intergenic")
    return annot


def run_mafft(seqs):
    """Run mafft on a list of (name, seq) pairs. Returns aligned (name, seq) pairs."""
    if len(seqs) == 1:
        return seqs
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for name, s in seqs:
            f.write(f">{name}\n{s}\n")
        in_path = f.name
    try:
        result = subprocess.run(['mafft', '--auto', '--quiet', in_path],
                                 capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return seqs
        # Parse FASTA output
        out = []
        cur_name = None; cur_parts = []
        for line in result.stdout.split('\n'):
            if line.startswith('>'):
                if cur_name:
                    out.append((cur_name, ''.join(cur_parts).upper()))
                cur_name = line[1:].strip()
                cur_parts = []
            else:
                cur_parts.append(line.strip())
        if cur_name:
            out.append((cur_name, ''.join(cur_parts).upper()))
        return out
    finally:
        os.unlink(in_path)


def main():
    groups = json.load(open('real_telotrons/linker_groups.json'))

    # Group by acc (host species)
    by_acc = defaultdict(list)
    for g in groups:
        acc = g['members'][0][0]
        by_acc[acc].append(g)

    # Load Eimeria genomes + GFFs
    contigs = defaultdict(dict)
    gff_data = {}
    gene_info_global = {}
    for acc in by_acc:
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if gen:
            cur, parts = None, []
            with open(gen) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur: contigs[acc][cur] = ''.join(parts)
                        cur = line[1:].strip().split()[0]; parts = []
                    else: parts.append(line.strip().upper())
                if cur: contigs[acc][cur] = ''.join(parts)
        gff = next(Path(f'genomes/{acc}').rglob('*.gff'), None)
        if gff:
            by_contig, gene_info = parse_gff_full(gff)
            gff_data[acc] = by_contig
            gene_info_global[acc] = gene_info
            print(f"  {acc}: {sum(len(v['genes']) for v in by_contig.values())} genes, "
                  f"{sum(len(v['introns']) for v in by_contig.values())} introns")

    # Telotron lookup
    telo_lookup = {}
    contig_telo = defaultdict(list)
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
                    contig_telo[(acc, row['contig'])].append((int(row['start']), int(row['end'])))
                except: pass

    def make_in_telotron(acc):
        def f(contig, pos):
            for s, e in contig_telo.get((acc, contig), []):
                if s <= pos <= e: return (s, e)
            return None
        return f

    FLANK = 200
    LINKER_FLANK = 80

    for acc, gs in by_acc.items():
        species = ACC_NAME.get(acc, acc)
        in_telo_fn = make_in_telotron(acc)
        rendered_groups = []

        for g in gs:
            members = [tuple(m) for m in g['members']]
            first = members[0]
            seq0 = telo_lookup.get(first, '')
            if not seq0: continue

            # Find canonical linker = longest non-telomeric run
            n = len(seq0); cov = bytearray(n)
            for k in {'TTTAGGG','TTAGGG','CCCTAAA','CCCTAA','TAGGGTT','AGGGTTT',
                      'GGGTTTA','GGTTTAG','GTTTAGG','AAACCCT','AACCCTA','ACCCTAA',
                      'CCCTAAA','CCTAAAC','CTAAACC'}:
                i = 0
                while True:
                    idx = seq0.find(k, i)
                    if idx == -1: break
                    for j in range(idx, idx+len(k)):
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

            # Build sequences to MSA: origin (if exists) + each member's linker±flank
            seqs_to_align = []  # (label, sequence)
            row_metadata = []  # (label, kind, annot, link_offset_in_seq)

            # ORIGIN
            if g['type'] == 'locus':
                org_acc = g['origin_acc']; org_contig = g['origin_contig']
                cs = contigs.get(org_acc, {}).get(org_contig, '')
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
                        annot = annotate_position(gff_data.get(org_acc, {}),
                                                    org_contig, gpos_s, in_telo_fn)
                        contig_len = len(cs)
                        d_to_end = min(gpos_s, contig_len - gpos_e)
                        if d_to_end < 5000:
                            annot.append(f"subtelo({d_to_end:,}bp)")
                        else:
                            annot.append(f"end={d_to_end:,}bp")
                        seqs_to_align.append((f"ORIGIN_{org_contig}_{gpos_s}{('rc' if strand=='rev' else '')}", ext_seq))
                        row_metadata.append({
                            'label': f"ORIGIN  {org_contig}:{gpos_s}-{gpos_e}{(' (rc)' if strand=='rev' else '')}",
                            'kind': 'ORIGIN',
                            'annot': '; '.join(annot),
                            'link_in_unaligned': link_off,
                        })

            # Skip locus groups where we couldn't localize origin
            if g['type'] == 'locus' and not seqs_to_align:
                continue

            # Members
            for m in members:
                mac, mc, ts, te = m
                seq = telo_lookup.get((mac, mc, ts, te), '')
                if not seq: continue
                pos, _, strand = find_best_match(canon, seq)
                if pos < 0: continue

                # Extract linker±flank
                if strand == 'rev':
                    s_in_telo = len(seq) - (pos + len(canon))
                    seq_ori = revcomp(seq)
                    pos_ori = s_in_telo
                else:
                    seq_ori = seq
                    pos_ori = pos

                ext_l = min(LINKER_FLANK, pos_ori)
                ext_r = min(LINKER_FLANK, len(seq_ori) - (pos_ori + len(canon)))
                ext_seq = seq_ori[pos_ori - ext_l : pos_ori + len(canon) + ext_r]
                link_off = ext_l

                annot = []
                if mac in gff_data:
                    a = annotate_position(gff_data[mac], mc, ts, in_telo_fn)
                    annot.extend(a)
                else:
                    annot.append(f"in_telotron[{ts}-{te}]")

                seqs_to_align.append((f"telo_{mc}_{ts}_{te}{'rc' if strand=='rev' else ''}", ext_seq))
                row_metadata.append({
                    'label': f"telo  {mc}:{ts}-{te}{(' (rc)' if strand=='rev' else '')}",
                    'kind': 'TELOTRON',
                    'annot': '; '.join(annot),
                    'link_in_unaligned': link_off,
                })

            if len(seqs_to_align) < 2: continue

            # Run MAFFT
            aligned = run_mafft(seqs_to_align)
            if len(aligned) != len(seqs_to_align):
                # fallback: pad
                continue

            # Find aligned linker positions
            # The reference (first input) had linker at position [link_off, link_off+len(canon)]
            # in aligned coords, find where positions map
            ref_aligned = aligned[0][1]
            ref_unaligned_to_aligned = []
            for i, c in enumerate(ref_aligned):
                if c != '-':
                    ref_unaligned_to_aligned.append(i)
            ref_link_off = row_metadata[0]['link_in_unaligned']
            if ref_link_off < len(ref_unaligned_to_aligned):
                aligned_link_start = ref_unaligned_to_aligned[ref_link_off]
                aligned_link_end_pos = ref_link_off + len(canon)
                if aligned_link_end_pos <= len(ref_unaligned_to_aligned):
                    aligned_link_end = ref_unaligned_to_aligned[aligned_link_end_pos - 1] + 1
                else:
                    aligned_link_end = len(ref_aligned)
            else:
                aligned_link_start = 0; aligned_link_end = len(canon)

            rendered_groups.append({
                'type': g['type'],
                'aligned': aligned,
                'metadata': row_metadata,
                'canon': canon,
                'aligned_link_start': aligned_link_start,
                'aligned_link_end': aligned_link_end,
            })

        if not rendered_groups: continue

        # Render: stack all groups for this species
        total_rows = sum(2 + len(g['aligned']) for g in rendered_groups)
        max_w = max(len(g['aligned'][0][1]) for g in rendered_groups)
        fig_w = max(20, max_w * 0.04)
        fig_h = max(8, total_rows * 0.45)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(-50, max_w + 5)
        ax.set_ylim(total_rows + 1, -2)
        ax.axis('off')

        plt.suptitle(f"{species}: linker MSA per group ({len(rendered_groups)} groups)",
                     fontsize=12, fontweight='bold')

        cur_y = 0
        for g_idx, g in enumerate(rendered_groups):
            n = len(g['aligned'])
            n_telo = sum(1 for m in g['metadata'] if m['kind'] == 'TELOTRON')
            has_origin = any(m['kind'] == 'ORIGIN' for m in g['metadata'])
            if has_origin:
                hdr = f"#{g_idx+1}  {n_telo} telotron(s) + origin (linker {len(g['canon'])}bp)"
            else:
                hdr = f"#{g_idx+1}  Shared linker among {n_telo} telotrons ({len(g['canon'])}bp)"
            ax.text(-50, cur_y+0.5, hdr, ha='left', va='center',
                    fontsize=10, fontweight='bold')
            cur_y += 1

            # Compute consensus-based mismatch detection
            # For each column, find majority base; mismatches = differ from majority
            ref_seq = g['aligned'][0][1]
            n_seqs = len(g['aligned'])
            consensus = []
            for col in range(len(ref_seq)):
                col_chars = [seq[col] for _, seq in g['aligned'] if col < len(seq) and seq[col] != '-']
                if col_chars:
                    from collections import Counter
                    c = Counter(col_chars)
                    consensus.append(c.most_common(1)[0][0])
                else:
                    consensus.append('-')

            for (label, seq), meta in zip(g['aligned'], g['metadata']):
                # Build per-base classes
                cls = []
                bg_overrides = []
                for j, c in enumerate(seq):
                    if c == '-':
                        cls.append('flank')  # we'll override bg to white below
                        bg_overrides.append((j, j+1, (1, 1, 1)))
                    elif c != consensus[j] and consensus[j] != '-':
                        cls.append('flank')
                        bg_overrides.append((j, j+1, (1.0, 0.7, 0.7)))  # mismatch = pink
                    else:
                        cls.append('flank')

                # Highlight linker region with green overlay
                hl = [(g['aligned_link_start'], g['aligned_link_end'], 'shared')]

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

        out = f'real_telotrons/fig_msa_real_{species}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Saved: {out}  ({len(rendered_groups)} groups)")


if __name__ == '__main__':
    main()
