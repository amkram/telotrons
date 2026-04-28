#!/usr/bin/env python3
"""
One figure per Eimeria species. For each species, every linker group is shown:
  Top row of group: origin (preferred non-telotron locus) with 200bp flanks.
  Below: each derived telotron's linker region with 200bp flanks.
All rows in a group are aligned by the linker (anchored MSA).
Each row is annotated with: gene context (if any), in-telotron status, position.
"""
import csv, json, re
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
    'GCF_000499745.2': 'Eimeria_mitis',
}


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def find_best_match(query, target):
    """Return (start, n_matches, strand) for best ungapped placement of query in target."""
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
    """Return per-contig list of (start, end, name, strand)."""
    genes = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 9: continue
            if p[2] != 'gene': continue
            attrs = dict(re.findall(r'(\w+)=([^;]+)', p[8]))
            name = attrs.get('Name', attrs.get('gene', '?'))
            biotype = attrs.get('gene_biotype', '?')
            try:
                genes[p[0]].append({
                    'start': int(p[3]), 'end': int(p[4]),
                    'strand': p[6], 'name': name, 'biotype': biotype,
                })
            except: pass
    return genes


def find_gene(gene_index, contig, pos):
    for g in gene_index.get(contig, []):
        if g['start'] <= pos <= g['end']:
            return g
    return None


def main():
    groups = json.load(open('real_telotrons/linker_groups.json'))

    # Group by accession (the "host" species of the source telotrons)
    by_acc = defaultdict(list)
    for g in groups:
        # group's host acc — use first member
        acc = g['members'][0][0]
        by_acc[acc].append(g)

    # Load contigs and gene index per acc
    contigs = defaultdict(dict)
    gene_index = {}
    for acc in by_acc:
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if gen:
            cur, parts = None, []
            with open(gen) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur: contigs[acc][cur] = ''.join(parts)
                        cur = line[1:].strip().split()[0]
                        parts = []
                    else: parts.append(line.strip().upper())
                if cur: contigs[acc][cur] = ''.join(parts)
        gff = next(Path(f'genomes/{acc}').rglob('*.gff'), None)
        if gff:
            gene_index[acc] = parse_gff_genes(gff)
            print(f"  {acc}: {len(contigs[acc])} contigs, "
                  f"{sum(len(v) for v in gene_index[acc].values())} genes")

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

    # Telotron interval index for fast in-telotron checks
    contig_telotron_intervals = defaultdict(list)
    for (acc, c, s, e) in telo_lookup:
        contig_telotron_intervals[(acc, c)].append((s, e))

    def in_telotron(acc, contig, pos):
        for s, e in contig_telotron_intervals.get((acc, contig), []):
            if s <= pos <= e:
                return (s, e)
        return None

    FLANK = 200
    PIN_LINKER_START = 250  # column in MSA where linker start is anchored

    # Build group records
    species_data = {}
    for acc, gs in by_acc.items():
        rendered = []
        for g in gs:
            members = g['members']
            # Compute canonical linker from first member
            first_acc, first_contig, ts, te = members[0]
            seq0 = telo_lookup.get((first_acc, first_contig, ts, te), '')
            if not seq0: continue

            # Linker = longest non-telomeric run in seq0
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
            best = (0, 0); s_run = None
            for i in range(n):
                if cov[i] == 0:
                    if s_run is None: s_run = i
                else:
                    if s_run is not None:
                        if i - s_run > best[1] - best[0]: best = (s_run, i)
                        s_run = None
            if s_run is not None and n - s_run > best[1] - best[0]: best = (s_run, n)
            if best[1] - best[0] < 15: continue
            canon = seq0[best[0]:best[1]]

            # Build rows
            rows = []

            # Origin row (preferred non-telotron locus)
            origin_row = None
            if g['type'] == 'locus':
                org_acc = g['origin_acc']
                org_contig = g['origin_contig']
                cs = contigs.get(org_acc, {}).get(org_contig, '')
                bin_s = g['origin_bin_start']; bin_e = g['origin_bin_end']
                region = cs[bin_s:bin_e] if cs else ''
                if region:
                    pos, m, strand = find_best_match(canon, region)
                    if pos >= 0:
                        gpos_s = bin_s + pos
                        gpos_e = gpos_s + len(canon)
                        if strand == 'rev':
                            # We want to show the strand where the linker is fwd
                            # so revcomp the surrounding sequence
                            ext_l = min(FLANK, len(cs) - gpos_e)
                            ext_r = min(FLANK, gpos_s)
                            ext_seq_raw = cs[gpos_s - ext_r:gpos_e + ext_l]
                            ext_seq = revcomp(ext_seq_raw)
                            link_off = ext_l
                        else:
                            ext_l = min(FLANK, gpos_s)
                            ext_r = min(FLANK, len(cs) - gpos_e)
                            ext_seq = cs[gpos_s - ext_l:gpos_e + ext_r]
                            link_off = ext_l

                        # Annotate position
                        gene = find_gene(gene_index.get(org_acc, {}), org_contig, gpos_s)
                        in_telo = in_telotron(org_acc, org_contig, gpos_s)
                        # Distance to contig end
                        contig_len = len(cs)
                        d_to_end = min(gpos_s, contig_len - gpos_e)
                        annot = []
                        if gene:
                            annot.append(f"gene={gene['name'][:14]}")
                        if in_telo:
                            annot.append(f"INSIDE telotron[{in_telo[0]}-{in_telo[1]}]")
                        else:
                            annot.append("intergenic")
                        if d_to_end < 5000:
                            annot.append(f"subtelo({d_to_end:,}bp from end)")
                        else:
                            annot.append(f"internal(end={d_to_end:,}bp)")

                        origin_row = {
                            'kind': 'ORIGIN',
                            'label': f"{org_contig}:{gpos_s}-{gpos_e}{('(rc)' if strand=='rev' else '')}",
                            'annot': '; '.join(annot),
                            'seq': ext_seq,
                            'link_offset': link_off,
                            'link_length': len(canon),
                            'family': 'TTTAGGG',
                            'is_telotron_seq': bool(in_telo),
                        }

            # If no origin row but type is shared, no separate origin section.
            # If no origin row and locus failed, skip this group.

            if origin_row is None and g['type'] == 'locus':
                continue

            # Member rows
            member_rows = []
            for m in members:
                mac, mc, ms, me = m
                seq = telo_lookup.get((mac, mc, ms, me), '')
                cs = contigs.get(mac, {}).get(mc, '')
                if not seq or not cs: continue

                # Find linker position in this telotron
                pos, n_match, strand = find_best_match(canon, seq)

                # Show full intron + 100bp flanks on each side
                flank = 100
                left = cs[max(0, ms-flank):ms]
                right = cs[me:me+flank]

                # If the canonical match is on rc strand, revcomp the whole thing
                if strand == 'rev':
                    full = revcomp(right) + revcomp(seq) + revcomp(left)
                    link_pos = len(revcomp(right)) + (len(seq) - (pos + len(canon)))
                else:
                    full = left + seq + right
                    link_pos = len(left) + pos

                # Annotation
                annot = []
                gene = find_gene(gene_index.get(mac, {}), mc, ms)
                if gene:
                    annot.append(f"gene={gene['name'][:14]}")
                annot.append(f"in_telotron[{ms}-{me}]")
                # is the telotron itself near a contig end?
                d_end = min(ms, len(cs) - me)
                if d_end < 5000:
                    annot.append(f"subtelo({d_end:,}bp)")

                member_rows.append({
                    'kind': 'TELOTRON',
                    'label': f"{mc}:{ms}-{me}{'(rc)' if strand=='rev' else ''}",
                    'annot': '; '.join(annot),
                    'seq': full,
                    'link_offset': link_pos,
                    'link_length': len(canon),
                    'family': 'TTTAGGG',
                    'intron_offset': len(left) if strand == 'fwd' else len(revcomp(right)),
                    'intron_length': len(seq),
                    'is_telotron_seq': True,
                })

            if not member_rows: continue

            # Determine MSA anchor
            # All rows should have their linker start at the same column
            # Determine pad needed for each
            all_rows = []
            if origin_row: all_rows.append(origin_row)
            all_rows.extend(member_rows)

            max_left_pad = max(r['link_offset'] for r in all_rows)
            max_right_pad = max(len(r['seq']) - r['link_offset'] - r['link_length']
                                  for r in all_rows)

            # Pre-pad each row with gaps
            for r in all_rows:
                left_pad = max_left_pad - r['link_offset']
                right_pad = (max_left_pad + r['link_length'] + max_right_pad) - \
                            (left_pad + len(r['seq']))
                r['_padded_seq'] = '-' * left_pad + r['seq'] + '-' * right_pad
                r['_aligned_link_start'] = max_left_pad
                r['_aligned_link_end'] = max_left_pad + r['link_length']
                # Also adjust intron offset if telotron
                if r['kind'] == 'TELOTRON':
                    r['_padded_intron_start'] = left_pad + r['intron_offset']
                    r['_padded_intron_end'] = r['_padded_intron_start'] + r['intron_length']

            rendered.append({
                'type': g['type'],
                'rows': all_rows,
                'canon': canon,
                'aligned_width': max_left_pad + len(canon) + max_right_pad,
            })

        species_data[acc] = rendered
        print(f"  {acc}: {len(rendered)} groups rendered (after filtering)")

    # Render one figure per species
    for acc, rendered in species_data.items():
        if not rendered:
            print(f"  Skip {acc}: no groups")
            continue
        species = ACC_NAME.get(acc, acc)
        # Layout: each group has 1 header + N rows + 1 spacer
        total_rows = sum(2 + len(g['rows']) for g in rendered)
        max_w = max(g['aligned_width'] for g in rendered)
        max_w = min(max_w, 2000)  # cap

        fig_w = max(20, max_w * 0.04)
        fig_h = max(8, total_rows * 0.45)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(-50, max_w + 5)
        ax.set_ylim(total_rows + 1, -2)
        ax.axis('off')

        plt.suptitle(f"{species}: linker→origin MSAs ({len(rendered)} groups)",
                     fontsize=13, fontweight='bold')

        cur_y = 0
        for g_idx, g in enumerate(rendered):
            n_telotrons = sum(1 for r in g['rows'] if r['kind'] == 'TELOTRON')
            has_origin = any(r['kind'] == 'ORIGIN' for r in g['rows'])
            if has_origin:
                hdr = f"#{g_idx+1}  {n_telotrons} telotron(s) → genomic origin (linker {len(g['canon'])}bp)"
            else:
                hdr = f"#{g_idx+1}  Shared linker between {n_telotrons} telotrons ({len(g['canon'])}bp)"
            ax.text(-50, cur_y + 0.5, hdr, ha='left', va='center',
                    fontsize=9.5, fontweight='bold')
            cur_y += 1

            # Render each row
            for r in g['rows']:
                seq = r['_padded_seq']
                # Build classes:
                # - gaps as 'flank' (light)
                # - within telotron-intron region: classify_bases
                # - outside telotron (= flanking exon for telotron rows OR genomic context for origin): 'flank'
                # - within linker: highlight as 'shared'
                cls = ['flank'] * len(seq)
                if r['kind'] == 'TELOTRON':
                    intron_s = r['_padded_intron_start']
                    intron_e = r['_padded_intron_end']
                    # Classify only the intron portion (no gaps inside it)
                    intron_seq = seq[intron_s:intron_e]
                    intron_clean = intron_seq.replace('-', '')
                    # Edge case: if there are gaps in the intron, our classification
                    # won't match exactly. Skip if so.
                    if len(intron_clean) == len(intron_seq):
                        ic = classify_bases(intron_seq, 'TTTAGGG')
                        for j, c in enumerate(ic):
                            if intron_s + j < len(cls):
                                cls[intron_s + j] = c
                else:
                    # Origin row — if it's in another telotron, classify the linker neighborhood
                    if r['is_telotron_seq']:
                        # Find a window around the linker that's a telomeric array
                        win_s = max(0, r['_aligned_link_start'] - 100)
                        win_e = min(len(seq), r['_aligned_link_end'] + 100)
                        win_seq = seq[win_s:win_e].replace('-', '')
                        if win_seq:
                            wc = classify_bases(win_seq, 'TTTAGGG')
                            # We'd need to map back through gaps which is messy.
                            # Simplification: just leave as flank.
                            pass

                # Mark gaps explicitly (no color)
                # Mark linker
                hl = [(r['_aligned_link_start'], r['_aligned_link_end'], 'shared')]

                # Render
                # Use space char for gaps
                seq_display = seq
                # We can't easily tell render_row to skip gaps; let's draw gaps as white
                # by overriding bg for gap positions
                bg_overrides = []
                for j, c in enumerate(seq_display):
                    if c == '-':
                        bg_overrides.append((j, j+1, (1.0, 1.0, 1.0)))

                # Label
                kind_label = 'ORIGIN' if r['kind'] == 'ORIGIN' else 'telo'
                kind_color = 'darkred' if r['kind'] == 'ORIGIN' else 'black'
                ax.text(-50, cur_y + 0.5,
                        f"{kind_label} {r['label']}",
                        ha='left', va='center', fontsize=7.5, color=kind_color, fontweight='bold' if r['kind']=='ORIGIN' else 'normal')

                render_row(ax, cur_y, seq_display, cls, x_start=0, char_width=1.0,
                           highlight_ranges=hl, bg_override=bg_overrides, fontsize=5.5)

                # Annotation text on the right
                ax.text(g['aligned_width'] + 3, cur_y + 0.5,
                        r['annot'], ha='left', va='center', fontsize=6.5, color='gray')

                # If telotron, draw splice site bars at intron boundaries
                if r['kind'] == 'TELOTRON':
                    ax.add_patch(Rectangle((r['_padded_intron_start'] - 0.4, cur_y),
                                            0.4, 1.0, facecolor='black'))
                    ax.add_patch(Rectangle((r['_padded_intron_end'], cur_y),
                                            0.4, 1.0, facecolor='black'))

                cur_y += 1
            cur_y += 1  # spacer

        out = f'real_telotrons/fig_msa_{species}.png'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  Saved: {out}")


if __name__ == '__main__':
    main()
