#!/usr/bin/env python3
"""
Strict-filtered per-species/per-MAG MSA figures with 250bp flanks.

Drops:
  - N-rich sequences (>5% N)
  - linkers that contain >1 exact telomeric kmer or ≥10% 1-mismatch coverage
  - canonical linker length < 25bp

Layout:
  - One figure per species/MAG
  - For each linker group: real MAFFT MSA of origin + linkers (250bp flanks)
  - Origin row at top (red bold)
  - Mismatches: pink. Gaps: white. Linker region: green overlay.
  - Annotations: in_intron(gene) | in_exon(gene) | intergenic | in_telotron[s-e]
                 + subtelo(Xbp) or end=Xbp.
"""
import csv, json, re, subprocess, tempfile, os
from pathlib import Path
from collections import defaultdict, Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from sequence_viewer import classify_bases, render_row


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


def find_canonical_linker(seq, family_kmers):
    """Longest run of bases not covered by any exact telomeric kmer in family_kmers."""
    n = len(seq); cov = bytearray(n)
    for k in family_kmers:
        klen = len(k); i = 0
        while True:
            idx = seq.find(k, i)
            if idx == -1: break
            for j in range(idx, idx+klen):
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
    if best[1] - best[0] < 25: return None
    return best


FAMILY_KMERS = {
    'TTAGGG':  set(),
    'TTTAGGG': set(),
}
for b in ['TTAGGG', 'CCCTAA']:
    for i in range(len(b)):
        FAMILY_KMERS['TTAGGG'].add(b[i:] + b[:i])
for b in ['TTTAGGG', 'CCCTAAA']:
    for i in range(len(b)):
        FAMILY_KMERS['TTTAGGG'].add(b[i:] + b[:i])


# Render figure for one species
def render_species_fig(species_label, groups, contigs, gff_data, telo_lookup,
                        contig_telo, family, flank_bp=250, out_path=None):
    rendered = []
    family_kmers = FAMILY_KMERS[family]

    def in_telo(contig, pos):
        for s, e in contig_telo.get(contig, []):
            if s <= pos <= e: return (s, e)
        return None

    for g in groups:
        members = [tuple(m) for m in g['members']]
        first = members[0]
        # first member format depends on type — last 3 are contig, start, end
        first_acc = first[0]
        first_contig = first[1]
        first_ts = first[2]; first_te = first[3]
        seq0 = telo_lookup.get((first_acc, first_contig, first_ts, first_te), '')
        if not seq0: continue
        # Canonical linker
        link_range = find_canonical_linker(seq0, family_kmers)
        if not link_range: continue
        canon = seq0[link_range[0]:link_range[1]]

        seqs_to_align = []
        row_meta = []

        # ORIGIN row
        if g['type'] == 'locus':
            org_acc = g.get('origin_acc') or g.get('mag')
            org_contig = g['origin_contig']
            cs = contigs.get(org_acc, {}).get(org_contig, '')
            bin_s, bin_e = g['origin_bin_start'], g['origin_bin_end']
            region = cs[bin_s:bin_e] if cs else ''
            if region:
                pos, _, strand = find_best_match(canon, region)
                if pos >= 0:
                    gpos_s = bin_s + pos
                    gpos_e = gpos_s + len(canon)
                    if strand == 'rev':
                        ext_l = min(flank_bp, len(cs) - gpos_e)
                        ext_r = min(flank_bp, gpos_s)
                        ext_seq = revcomp(cs[gpos_s - ext_r : gpos_e + ext_l])
                        link_off = ext_l
                    else:
                        ext_l = min(flank_bp, gpos_s)
                        ext_r = min(flank_bp, len(cs) - gpos_e)
                        ext_seq = cs[gpos_s - ext_l : gpos_e + ext_r]
                        link_off = ext_l
                    # Drop if Ns dominate
                    n_frac = ext_seq.count('N') / len(ext_seq) if ext_seq else 1
                    if n_frac > 0.05: continue
                    annot = annotate_position(gff_data.get(org_acc, {}), org_contig, gpos_s, in_telo)
                    contig_len = len(cs)
                    d_to_end = min(gpos_s, contig_len - gpos_e)
                    annot.append(f"subtelo({d_to_end:,}bp)" if d_to_end < 5000
                                 else f"end={d_to_end:,}bp")
                    seqs_to_align.append((f"ORIGIN_{gpos_s}", ext_seq))
                    row_meta.append({
                        'label': f"ORIGIN  {org_contig}:{gpos_s}-{gpos_e}{(' (rc)' if strand=='rev' else '')}",
                        'kind': 'ORIGIN',
                        'annot': '; '.join(annot),
                        'link_in_unaligned': link_off,
                    })

        if g['type'] == 'locus' and not seqs_to_align: continue

        for m in members:
            mac, mc, ts, te = m
            seq = telo_lookup.get((mac, mc, ts, te), '')
            if not seq: continue
            pos, _, strand = find_best_match(canon, seq)
            if pos < 0: continue
            cs = contigs.get(mac, {}).get(mc, '')
            if not cs: continue
            if strand == 'rev':
                left_raw = cs[max(0, ts-flank_bp):ts] if cs else ''
                right_raw = cs[te:te+flank_bp] if cs else ''
                whole = left_raw + seq + right_raw
                whole_rc = revcomp(whole)
                rev_link_s = len(whole) - (len(left_raw) + pos + len(canon))
                ext_l = min(flank_bp, rev_link_s)
                ext_r = min(flank_bp, len(whole) - rev_link_s - len(canon))
                ext_seq = whole_rc[rev_link_s - ext_l : rev_link_s + len(canon) + ext_r]
                link_off = ext_l
                # intron position in displayed (rc): ext_l - (len(seq) - pos - len(canon))
                intron_start_disp = ext_l - (len(seq) - pos - len(canon))
                intron_end_disp = intron_start_disp + len(seq)
            else:
                gpos_s = ts + pos
                gpos_e = gpos_s + len(canon)
                ext_l = min(flank_bp, gpos_s)
                ext_r = min(flank_bp, len(cs) - gpos_e)
                ext_seq = cs[gpos_s - ext_l : gpos_e + ext_r]
                link_off = ext_l
                intron_start_disp = ext_l - pos
                intron_end_disp = intron_start_disp + len(seq)

            n_frac = ext_seq.count('N') / len(ext_seq) if ext_seq else 1
            if n_frac > 0.05: continue

            annot = []
            if mac in gff_data:
                annot.extend(annotate_position(gff_data[mac], mc, ts, in_telo))
            else:
                t = in_telo(mc, ts)
                if t: annot.append(f"in_telotron[{t[0]}-{t[1]}]")

            seqs_to_align.append((f"telo_{mc}_{ts}{'rc' if strand=='rev' else ''}", ext_seq))
            row_meta.append({
                'label': f"telo  {mc}:{ts}-{te}{(' (rc)' if strand=='rev' else '')}",
                'kind': 'TELOTRON',
                'annot': '; '.join(annot),
                'link_in_unaligned': link_off,
                'intron_start_unaligned': intron_start_disp,
                'intron_end_unaligned': intron_end_disp,
            })

        if len(seqs_to_align) < 2: continue
        aligned = run_mafft(seqs_to_align)
        if len(aligned) != len(seqs_to_align): continue

        # Map unaligned position -> aligned column for each row
        def unaln_to_aln(aln_str, unaln_pos):
            """Given aligned sequence (with '-'), return aligned column index
            corresponding to unaligned base index unaln_pos. Clamps to bounds."""
            if unaln_pos < 0: return 0
            cnt = 0
            for i, c in enumerate(aln_str):
                if c != '-':
                    if cnt == unaln_pos: return i
                    cnt += 1
            return len(aln_str)

        # Find aligned linker columns (use first row as reference)
        ref_aligned = aligned[0][1]
        ref_link_off = row_meta[0]['link_in_unaligned']
        aligned_link_start = unaln_to_aln(ref_aligned, ref_link_off)
        aligned_link_end = unaln_to_aln(ref_aligned, ref_link_off + len(canon))

        # Compute aligned intron boundaries for each TELOTRON row
        for (name, aln_str), meta in zip(aligned, row_meta):
            if meta['kind'] == 'TELOTRON':
                meta['intron_start_aligned'] = unaln_to_aln(aln_str, meta['intron_start_unaligned'])
                meta['intron_end_aligned'] = unaln_to_aln(aln_str, meta['intron_end_unaligned'])

        rendered.append({
            'type': g['type'], 'aligned': aligned, 'meta': row_meta,
            'canon': canon,
            'alink_s': aligned_link_start, 'alink_e': aligned_link_end,
        })

    if not rendered: return False

    # Render
    total_rows = sum(2 + len(g['aligned']) for g in rendered)
    max_w = max(len(g['aligned'][0][1]) for g in rendered)
    fig_w = max(20, max_w * 0.04)
    fig_h = max(8, total_rows * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-50, max_w + 5); ax.set_ylim(total_rows + 1, -2)
    ax.axis('off')

    plt.suptitle(f"{species_label}: linker MSA per group ({len(rendered)} groups, ±{flank_bp}bp flanks)",
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
            col_chars = [s[col] for _, s in g['aligned']
                         if col < len(s) and s[col] != '-']
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
                       highlight_ranges=hl, bg_override=bg_overrides, fontsize=4.5)
            # Telotron intron boundary bars
            if meta['kind'] == 'TELOTRON':
                i_s = meta.get('intron_start_aligned', 0)
                i_e = meta.get('intron_end_aligned', len(seq))
                # 5' splice site bar (intron start)
                if 0 <= i_s <= len(seq):
                    ax.add_patch(Rectangle((i_s - 0.4, cur_y), 0.7, 1.0,
                                             facecolor='black', edgecolor='none',
                                             zorder=5))
                # 3' splice site bar (intron end)
                if 0 <= i_e <= len(seq):
                    ax.add_patch(Rectangle((i_e - 0.4, cur_y), 0.7, 1.0,
                                             facecolor='black', edgecolor='none',
                                             zorder=5))
            ax.text(max_w+3, cur_y+0.5, meta['annot'],
                    ha='left', va='center', fontsize=6.5, color='gray')
            cur_y += 1
        cur_y += 1

    plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def main():
    # === Eimeria ===
    eim_groups = json.load(open('real_telotrons/strict_linker_groups_eimeria.json'))
    by_acc_eim = defaultdict(list)
    for g in eim_groups:
        first_acc = g['members'][0][0]
        by_acc_eim[first_acc].append(g)

    # Load Eimeria data
    eim_contigs = defaultdict(dict)
    eim_gff = {}
    for acc in by_acc_eim:
        gen = next(Path(f'genomes/{acc}').rglob('*.fna'), None)
        if gen:
            cur, parts = None, []
            with open(gen) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur: eim_contigs[acc][cur] = ''.join(parts)
                        cur = line[1:].strip().split()[0]; parts = []
                    else: parts.append(line.strip().upper())
                if cur: eim_contigs[acc][cur] = ''.join(parts)
        gff = next(Path(f'genomes/{acc}').rglob('*.gff'), None)
        if gff:
            eim_gff[acc] = parse_gff_full(gff)

    eim_telo_lookup = {}
    eim_contig_telo = defaultdict(list)
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                seq = row.get('intron_seq', '').upper()
                if not seq: continue
                # Drop N-rich introns
                if seq.count('N') / len(seq) > 0.05: continue
                try:
                    key = (acc, row['contig'], int(row['start']), int(row['end']))
                    eim_telo_lookup[key] = seq
                    eim_contig_telo[(acc, row['contig'])].append((int(row['start']), int(row['end'])))
                except: pass

    for acc, gs in by_acc_eim.items():
        species = ACC_NAME.get(acc, acc)
        # Build per-acc contig_telo subset
        per_contig = defaultdict(list)
        for (a, c), v in eim_contig_telo.items():
            if a == acc: per_contig[c] = v
        out = f'real_telotrons/fig_msa_strict_{species}.png'
        ok = render_species_fig(species, gs, eim_contigs, eim_gff,
                                  eim_telo_lookup, per_contig, 'TTTAGGG',
                                  flank_bp=250, out_path=out)
        if ok: print(f"Saved: {out}")
        else: print(f"  No valid groups for {species}")

    # === Tara ===
    tara_groups = json.load(open('real_telotrons/strict_linker_groups_tara.json'))
    by_mag = defaultdict(list)
    for g in tara_groups:
        by_mag[g['members'][0][0]].append(g)

    GENOME_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/contigs_individual')
    GFF_DIR = Path('/Users/alex/telotrons/tara_oceans_euk_mags/smags/gff_individual/GFF')

    tara_contigs = defaultdict(dict)
    tara_gff = {}
    for mag in by_mag:
        fa = GENOME_DIR / f"{mag}.fa"
        if fa.exists():
            cur, parts = None, []
            with open(fa) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur: tara_contigs[mag][cur] = ''.join(parts)
                        cur = line[1:].strip().split()[0]; parts = []
                    else: parts.append(line.strip().upper())
                if cur: tara_contigs[mag][cur] = ''.join(parts)
        gff = GFF_DIR / f"{mag}.gmove.gff"
        if gff.exists():
            tara_gff[mag] = parse_gff_full(gff)

    tara_telo_lookup = {}
    tara_contig_telo = defaultdict(list)
    with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            seq = row.get('intron_seq', '').upper()
            if not seq: continue
            if seq.count('N') / len(seq) > 0.05: continue
            try:
                key = (row['mag'], row['contig'], int(row['start']), int(row['end']))
                tara_telo_lookup[key] = seq
                tara_contig_telo[(row['mag'], row['contig'])].append(
                    (int(row['start']), int(row['end'])))
            except: pass

    for mag, gs in by_mag.items():
        per_contig = defaultdict(list)
        for (m, c), v in tara_contig_telo.items():
            if m == mag: per_contig[c] = v
        out = f'real_telotrons/fig_msa_strict_{mag}.png'
        ok = render_species_fig(mag, gs, tara_contigs, tara_gff,
                                  tara_telo_lookup, per_contig, 'TTAGGG',
                                  flank_bp=250, out_path=out)
        if ok: print(f"Saved: {out}")
        else: print(f"  No valid groups for {mag}")


if __name__ == '__main__':
    main()
