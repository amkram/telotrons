"""Group the >85% single-motif loci of one genome into duplicates by BLAST of both 250 bp flanks."""

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path

MINIMUM_IDENTITY = 95
MINIMUM_ALIGNMENT_BP = 100
FLANK_BP = 250


def read(path):
    with (gzip.open if str(path).endswith('.gz') else open)(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def primary(hit):
    """>85% of the intron under the screen rule or the tract rule, after the flank filters."""
    length = int(hit['intron_bp'])
    screened = int(hit['period']) >= 6 or int(hit['is_telomeric'])
    screen = int(hit['screen_pass']) and screened and int(hit['screen_bp']) * 100 > 85 * length
    tract = int(hit['tract_pass']) and int(hit['tract_bp']) * 100 > 85 * length
    return screen or tract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', required=True, help='scan.py output of one genome')
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    directory = Path(args.directory)

    # job scratch for the BLAST files, under the caller's TMPDIR
    scratch = Path(os.environ['TMPDIR']).resolve(strict=True) / f'deduplicate_{directory.name}'
    scratch.mkdir(exist_ok=True)

    loci = {hit['locus_id'] for hit in read(directory / 'hits.tsv.gz') if primary(hit)}
    rows = [row for row in read(directory / 'candidates.tsv.gz') if row['locus_id'] in loci]
    strands = {int(row['index']): row['strand'] for row in rows}

    # union-find over candidate indices, one root per duplicate group
    parent = {index: index for index in strands}

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first, second):
        first = root(first)
        second = root(second)
        if first != second:
            parent[max(first, second)] = min(first, second)

    # flanks of the primary loci, minus-strand ones reverse-complemented with U and D swapped
    flanks = scratch / 'flanks.fa'
    alignments = scratch / 'alignments.tsv'
    edges = set()
    with gzip.open(directory / 'flanks.fa.gz', 'rt') as source, flanks.open('w') as out:
        keep = False
        for line in source:
            if line.startswith('>'):
                index, side = line[2:].strip().split('|')
                index = int(index)
                keep = index in strands
                reverse = keep and strands[index] == '-'
                if reverse:
                    line = f'>L{index}|' + ('D' if side == 'U' else 'U') + '\n'
            elif reverse:
                line = line.strip().translate(str.maketrans('ACGT', 'TGCA'))[::-1] + '\n'
            if keep:
                out.write(line)

    # all-against-all BLAST of the flanks
    if len(strands) > 1:
        database = str(scratch / 'flanks')
        subprocess.run(
            ['makeblastdb', '-in', str(flanks), '-dbtype', 'nucl', '-out', database],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            [
                'blastn',
                '-query',
                str(flanks),
                '-db',
                database,
                '-out',
                str(alignments),
                '-outfmt',
                '6 qseqid sseqid pident length',
                '-evalue',
                '1e-10',
                '-word_size',
                '11',
                '-dust',
                'no',
                '-max_target_seqs',
                str(2 * len(strands)),
                '-num_threads',
                str(args.threads),
            ],
            check=True,
        )

        # two loci are duplicates when both U and D flanks match (>= 95% identity over >= 100 bp)
        upstream = set()
        downstream = set()
        with alignments.open() as handle:
            for line in handle:
                query, subject, identity, length = line.split()
                query_index, query_side = query[1:].split('|')
                subject_index, subject_side = subject[1:].split('|')
                pair = tuple(sorted([int(query_index), int(subject_index)]))
                if query_side != subject_side or pair[0] == pair[1]:
                    continue
                if float(identity) < MINIMUM_IDENTITY or int(length) < MINIMUM_ALIGNMENT_BP:
                    continue
                (upstream if query_side == 'U' else downstream).add(pair)
        edges = upstream & downstream
        for first, second in edges:
            union(first, second)

    # group of each locus and the run summary
    with (directory / 'groups.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['locus_id', 'group'])
        for row in rows:
            writer.writerow([row['locus_id'], root(int(row['index']))])
    summary = dict(
        primary_intervals=len(rows),
        groups=len({root(index) for index in strands}),
        paired_flank_edges=len(edges),
        minimum_identity=MINIMUM_IDENTITY,
        minimum_alignment_bp=MINIMUM_ALIGNMENT_BP,
        flank_bp=FLANK_BP,
        max_target_seqs=2 * len(strands),
    )
    (directory / 'dedup.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(directory.name, summary, flush=True)
    shutil.rmtree(scratch)


if __name__ == '__main__':
    main()
