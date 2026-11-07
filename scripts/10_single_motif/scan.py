"""Score every annotated intron of one genome against every primitive motif of period 1-8."""

import argparse
import collections
import csv
import gzip
import json
import time
from pathlib import Path

from annotation import annotated_introns, fasta
from repeats import Scanner, canonical

TELOMERIC = [
    'TTAGG',
    'TCAGG',
    'TTGGG',
    'TTAGGG',
    'TTAGGC',
    'TTGGGG',
    'TTTAGGG',
    'TTTTAGGG',
    'TTTTGGGG',
    'TATTAGGG',
    'TTATTGGG',
]
# intron length bin edges (bp)
LENGTH_BOUNDS = [50, 100, 200, 500, 1000, 2000, 10000]
HIT_FIELDS = [
    'locus_id',
    'motif',
    'period',
    'screen_bp',
    'tract_bp',
    'screen_pass',
    'tract_pass',
    'intron_bp',
    'is_telomeric',
    'is_lineage_motif',
]
CANDIDATE_FIELDS = [
    'index',
    'locus_id',
    'seqid',
    'start',
    'end',
    'strand',
    'intron_bp',
    'gc',
    'length_bin',
    'gc_bin',
    'ambiguous_bp',
]


def write_rows(path, rows, fields):
    with (gzip.open if str(path).endswith('.gz') else open)(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--library', required=True, help='compiled repeats.c')
    args = parser.parse_args()
    with open(args.manifest) as handle:
        meta = next(
            row for row in csv.DictReader(handle, delimiter='\t') if row['genome_id'] == args.genome
        )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    scanner = Scanner(args.library)
    telomeric_ids = {scanner.ids[canonical(motif)] for motif in TELOMERIC}
    lineage_id = scanner.ids[canonical(meta['lineage_motif'])]

    def screened(motif_id):
        # the screen measure applies to motifs of period >= 6 and to telomere motifs
        return len(scanner.motifs[motif_id]) >= 6 or motif_id in telomeric_ids

    # intron intervals and their strands per contig, from the GFF
    introns = collections.defaultdict(dict)
    for contig, start, end, strand in annotated_introns(meta['gff']):
        introns[contig].setdefault((start, end), set()).add(strand)

    totals = collections.Counter(annotated_intervals=sum(map(len, introns.values())))
    denominators = collections.Counter()
    motif_counts = collections.Counter()
    seen = set()
    candidates = []
    next_progress = 100000

    # stream the assembly, writing hits and flanks as we go
    with gzip.open(out / 'hits.tsv.gz', 'wt', compresslevel=5) as hits_handle, gzip.open(
        out / 'flanks.fa.gz', 'wt', compresslevel=5
    ) as flanks_handle:
        hit_writer = csv.DictWriter(
            hits_handle, fieldnames=HIT_FIELDS, delimiter='\t', lineterminator='\n'
        )
        hit_writer.writeheader()
        for contig, sequence in fasta(meta['fasta']):
            if contig not in introns:
                continue
            seen.add(contig)
            for (start, end), strands in sorted(introns[contig].items()):
                # assessable introns: inside the contig (not wrapping a circular organelle),
                # >= 30 bp, with 12 bp flanks
                if not 0 <= start < end <= len(sequence):
                    totals['outside_sequence'] += 1
                    continue
                length = end - start
                if length < 30:
                    totals['shorter_than_30'] += 1
                    continue
                if start < 12 or end + 12 > len(sequence):
                    totals['short_flank'] += 1
                    continue

                # denominator by length bin and GC bin (fifths)
                intron = sequence[start:end]
                length_bin = sum(length >= bound for bound in LENGTH_BOUNDS)
                acgt = sum(intron.count(base) for base in 'ACGT')
                gc = (intron.count('G') + intron.count('C')) / acgt if acgt else 0
                gc_bin = min(4, int(gc * 5))
                denominators[length_bin, gc_bin] += 1
                totals['assessable_introns'] += 1
                totals['ambiguous_introns'] += acgt < length

                # motifs over 50% of the intron
                selected = [
                    (motif_id, screen_bp, tract_bp)
                    for motif_id, screen_bp, tract_bp in scanner.scan(intron)
                    if (screened(motif_id) and screen_bp * 2 > length) or tract_bp * 2 > length
                ]
                if not selected:
                    continue
                left = sequence[max(0, start - 200) : start].encode()
                right = sequence[end : end + 200].encode()

                passing = []
                for motif_id, screen_bp, tract_bp in selected:
                    motif = scanner.motifs[motif_id]
                    measures = [('screen', screen_bp), ('tract', tract_bp)]
                    for method, covered in measures:
                        if method == 'screen' and not screened(motif_id):
                            continue
                        for threshold in [50, 85]:
                            if covered * 100 > threshold * length:
                                motif_counts[motif, method, threshold, 'before_flanks'] += 1

                    # junction (< 6 of 12 bp at each end) and 200 bp flank (<= 20%) repeat filters
                    screen_pass = False
                    tract_pass = False
                    if screened(motif_id) and screen_bp * 2 > length:
                        left_junction = scanner.coverage(left[-12:], motif_id)[1]
                        right_junction = scanner.coverage(right[:12], motif_id)[1]
                        flank_bp = (
                            scanner.coverage(left, motif_id)[0]
                            + scanner.coverage(right, motif_id)[0]
                        )
                        screen_pass = max(
                            left_junction, right_junction
                        ) < 6 and flank_bp * 5 <= len(left) + len(right)
                    if tract_bp * 2 > length:
                        left_tract = scanner.tract(left, motif_id)
                        right_tract = scanner.tract(right, motif_id)
                        tract_pass = max(left_tract[2], right_tract[1]) < 6 and (
                            (left_tract[0] + right_tract[0]) * 5 <= len(left) + len(right)
                        )
                    if not (screen_pass or tract_pass):
                        continue

                    passing.append((motif_id, screen_bp, tract_bp, screen_pass, tract_pass))
                    for method, covered in measures:
                        if not (screen_pass if method == 'screen' else tract_pass):
                            continue
                        if method == 'screen' and not screened(motif_id):
                            continue
                        for threshold in [50, 85]:
                            if covered * 100 > threshold * length:
                                motif_counts[motif, method, threshold, 'after_flanks'] += 1
                if not passing:
                    continue

                # record the candidate, its 250 bp flanks and its hits
                locus_id = f'{args.genome}|{contig}:{start + 1}-{end}'
                index = len(candidates)
                candidates.append(
                    dict(
                        index=index,
                        locus_id=locus_id,
                        seqid=contig,
                        start=start + 1,
                        end=end,
                        strand=next(iter(strands)) if len(strands) == 1 else '.',
                        intron_bp=length,
                        gc=gc,
                        length_bin=length_bin,
                        gc_bin=gc_bin,
                        ambiguous_bp=length - acgt,
                    )
                )
                upstream = sequence[max(0, start - 250) : start]
                downstream = sequence[end : end + 250]
                flanks_handle.write(f'>L{index}|U\n{upstream}\n>L{index}|D\n{downstream}\n')
                for motif_id, screen_bp, tract_bp, screen_pass, tract_pass in passing:
                    hit_writer.writerow(
                        dict(
                            locus_id=locus_id,
                            motif=scanner.motifs[motif_id],
                            period=len(scanner.motifs[motif_id]),
                            screen_bp=screen_bp,
                            tract_bp=tract_bp,
                            screen_pass=int(screen_pass),
                            tract_pass=int(tract_pass),
                            intron_bp=length,
                            is_telomeric=int(motif_id in telomeric_ids),
                            is_lineage_motif=int(motif_id == lineage_id),
                        )
                    )
            if totals['assessable_introns'] >= next_progress:
                print(args.genome, totals['assessable_introns'], 'introns', flush=True)
                next_progress = totals['assessable_introns'] + 100000

    # every contig with introns must be in the assembly
    assert set(introns) == seen, sorted(set(introns) - seen)[:10]

    # candidates, denominators, motif counts and the run summary
    write_rows(out / 'candidates.tsv.gz', candidates, CANDIDATE_FIELDS)
    write_rows(
        out / 'denominators.tsv',
        [
            dict(length_bin=length_bin, gc_bin=gc_bin, introns=count)
            for (length_bin, gc_bin), count in sorted(denominators.items())
        ],
        ['length_bin', 'gc_bin', 'introns'],
    )
    write_rows(
        out / 'motifs.tsv',
        [
            dict(motif=motif, method=method, threshold=threshold, stage=stage, count=count)
            for (motif, method, threshold, stage), count in sorted(motif_counts.items())
        ],
        ['motif', 'method', 'threshold', 'stage', 'count'],
    )
    totals['repeat_candidate_intervals'] = len(candidates)
    summary = dict(genome=meta, counts=dict(totals), seconds=time.time() - started, complete=True)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    scanner.close()
    print(args.genome, dict(totals), round(time.time() - started, 1), 'seconds', flush=True)


if __name__ == '__main__':
    main()
