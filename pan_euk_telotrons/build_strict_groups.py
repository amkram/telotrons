#!/usr/bin/env python3
"""Build linker groups from strict-filtered BLAST results."""
import csv, json, re
from pathlib import Path
from collections import defaultdict


def build_eimeria_groups():
    # Load Eimeria telotron records
    contig_telos = defaultdict(list)
    for f in Path('ultra_results').glob('*Eimeria*.tsv'):
        m = re.match(r'(GCF_\d+\.\d+)_', f.name)
        if not m: continue
        acc = m.group(1)
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                try:
                    contig_telos[(acc, row['contig'])].append(
                        (int(row['start']), int(row['end'])))
                except: pass

    def find_telotron(acc, contig, start, end):
        for s, e in contig_telos.get((acc, contig), []):
            if start >= s - 50 and end <= e + 50:
                return (s, e)
        return None

    all_hits = []
    for hf in Path('blast_out/eimeria_strict').glob('*_hits.tsv'):
        acc = hf.name.replace('_hits.tsv', '')
        with open(hf) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 14: continue
                qseq, sseq = p[0], p[1]
                pident = float(p[2]); aln_len = int(p[3])
                qstart, qend = int(p[6]), int(p[7])
                sstart, send = int(p[8]), int(p[9])
                qlen = int(p[12])
                qcov = (qend - qstart + 1) / qlen
                qid = qseq.split('|')
                if len(qid) < 5: continue
                src_contig = qid[2]
                try: src_s, src_e = (int(x) for x in qid[3].split('-'))
                except: continue
                mtch = re.search(r'(NW_\d+\.\d+)', sseq)
                sseq_contig = mtch.group(1) if mtch else sseq
                self_hit = (sseq_contig == src_contig and
                           min(sstart, send) <= src_e + 200 and
                           max(sstart, send) >= src_s - 200)
                if self_hit: continue
                if qcov < 0.80 or pident < 85: continue
                all_hits.append({
                    'acc': acc, 'src_contig': src_contig,
                    'src_start': src_s, 'src_end': src_e,
                    'origin_contig': sseq_contig,
                    'origin_start': min(sstart, send),
                    'origin_end': max(sstart, send),
                })
    print(f"Eimeria non-self HQ hits: {len(all_hits)}")

    groups = defaultdict(set)
    for h in all_hits:
        src_telo = None
        for s, e in contig_telos.get((h['acc'], h['src_contig']), []):
            if s <= h['src_start'] and e >= h['src_end']:
                src_telo = (s, e); break
        if not src_telo: continue
        dest = find_telotron(h['acc'], h['origin_contig'],
                              h['origin_start'], h['origin_end'])
        a = (h['acc'], h['src_contig'], src_telo[0], src_telo[1])
        if dest:
            b = (h['acc'], h['origin_contig'], dest[0], dest[1])
            canon = tuple(sorted([a, b]))
            groups[('SHARED', canon)].add(a); groups[('SHARED', canon)].add(b)
        else:
            bin_id = (h['acc'], h['origin_contig'],
                       h['origin_start'] // 1000, h['origin_end'] // 1000)
            groups[('LOCUS',) + bin_id].add(a)

    out = []
    for gid, members in groups.items():
        if gid[0] == 'SHARED':
            out.append({'type': 'shared', 'members': [list(m) for m in members]})
        else:
            out.append({
                'type': 'locus',
                'origin_acc': gid[1], 'origin_contig': gid[2],
                'origin_bin_start': gid[3]*1000,
                'origin_bin_end': (gid[4]+1)*1000,
                'members': [list(m) for m in members],
            })
    return out


def build_tara_groups():
    tara = defaultdict(list)
    with open('real_telotrons/_tara_oceans_telotrons.tsv') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            try:
                tara[(row['mag'], row['contig'])].append(
                    (int(row['start']), int(row['end'])))
            except: pass

    def find_telotron(mag, contig, start, end):
        for s, e in tara.get((mag, contig), []):
            if start >= s - 50 and end <= e + 50: return (s, e)
        return None

    all_hits = []
    for hf in Path('blast_out/tara_strict').glob('*_hits.tsv'):
        mag = hf.name.replace('_hits.tsv', '')
        with open(hf) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 14: continue
                qseq, sseq = p[0], p[1]
                pident = float(p[2]); aln_len = int(p[3])
                qstart, qend = int(p[6]), int(p[7])
                sstart, send = int(p[8]), int(p[9])
                qlen = int(p[12])
                qcov = (qend - qstart + 1) / qlen
                qid = qseq.split('|')
                if len(qid) < 5: continue
                src_contig = qid[2]
                try: src_s, src_e = (int(x) for x in qid[3].split('-'))
                except: continue
                sseq_contig = sseq
                self_hit = (sseq_contig == src_contig and
                           min(sstart, send) <= src_e + 200 and
                           max(sstart, send) >= src_s - 200)
                if self_hit: continue
                if qcov < 0.80 or pident < 85: continue
                all_hits.append({
                    'mag': mag, 'src_contig': src_contig,
                    'src_start': src_s, 'src_end': src_e,
                    'origin_contig': sseq_contig,
                    'origin_start': min(sstart, send),
                    'origin_end': max(sstart, send),
                })
    print(f"Tara non-self HQ hits: {len(all_hits)}")

    groups = defaultdict(set)
    for h in all_hits:
        src_telo = None
        for s, e in tara.get((h['mag'], h['src_contig']), []):
            if s <= h['src_start'] and e >= h['src_end']:
                src_telo = (s, e); break
        if not src_telo: continue
        dest = find_telotron(h['mag'], h['origin_contig'],
                              h['origin_start'], h['origin_end'])
        a = (h['mag'], h['src_contig'], src_telo[0], src_telo[1])
        if dest:
            b = (h['mag'], h['origin_contig'], dest[0], dest[1])
            canon = tuple(sorted([a, b]))
            groups[('SHARED', canon)].add(a); groups[('SHARED', canon)].add(b)
        else:
            bin_id = (h['mag'], h['origin_contig'],
                       h['origin_start'] // 1000, h['origin_end'] // 1000)
            groups[('LOCUS',) + bin_id].add(a)

    out = []
    for gid, members in groups.items():
        if gid[0] == 'SHARED':
            out.append({'type': 'shared',
                       'mag': list(members)[0][0],
                       'members': [list(m) for m in members]})
        else:
            out.append({
                'type': 'locus',
                'mag': gid[1],
                'origin_contig': gid[2],
                'origin_bin_start': gid[3]*1000,
                'origin_bin_end': (gid[4]+1)*1000,
                'members': [list(m) for m in members],
            })
    return out


if __name__ == '__main__':
    eim_groups = build_eimeria_groups()
    print(f"Eimeria groups: {len(eim_groups)} "
          f"(shared={sum(1 for g in eim_groups if g['type']=='shared')}, "
          f"locus={sum(1 for g in eim_groups if g['type']=='locus')})")
    with open('real_telotrons/strict_linker_groups_eimeria.json', 'w') as f:
        json.dump(eim_groups, f, indent=1)

    tara_groups = build_tara_groups()
    print(f"Tara groups: {len(tara_groups)} "
          f"(shared={sum(1 for g in tara_groups if g['type']=='shared')}, "
          f"locus={sum(1 for g in tara_groups if g['type']=='locus')})")
    with open('real_telotrons/strict_linker_groups_tara.json', 'w') as f:
        json.dump(tara_groups, f, indent=1)
