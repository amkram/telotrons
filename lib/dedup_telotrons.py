"""Collapse loci whose upstream and downstream flanks both match by blastn.

A match spans >= MIN_LEN bp at >= MIN_PIDENT % identity.
"""

from collections import defaultdict
import subprocess

from _common import rc, slug

MIN_PIDENT = 95.0
MIN_LEN = 100
FLANK = 250


def locus_flanks(row, chrom):
    """Up to FLANK bp upstream and downstream of a 1-based locus, in transcript orientation."""
    start, end = int(row['start']), int(row['end'])
    left = chrom[max(0, start - 1 - FLANK) : start - 1]
    right = chrom[end : min(len(chrom), end + FLANK)]

    # minus strand: swap sides and reverse-complement
    if row['strand'] == '-':
        return rc(right), rc(left)
    return left, right


def write_flanks(path, rows, contigs):
    """Write flanks as L<row>|U (upstream) and L<row>|D (downstream), True when any was written."""
    written = False
    with open(path, 'w') as handle:
        for i, row in enumerate(rows):
            chrom = contigs.get(row['seqid'], '')
            if not chrom:
                continue
            upstream, downstream = locus_flanks(row, chrom)
            if upstream:
                handle.write(f'>L{i}|U\n{upstream}\n')
                written = True
            if downstream:
                handle.write(f'>L{i}|D\n{downstream}\n')
                written = True
    return written


def blast_flanks(fasta_path, prefix, threads):
    """Path of the tabular all-against-all blastn hits between flanks."""
    database = f'{prefix}_db'
    hits = f'{prefix}_blast.tsv'
    subprocess.run(
        ['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl', '-out', database],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # no -max_target_seqs: the order-of-encounter cap drops true partners in ~3000-locus genomes
    subprocess.run(
        [
            'blastn',
            '-query',
            fasta_path,
            '-db',
            database,
            '-out',
            hits,
            '-outfmt',
            '6 qseqid sseqid pident length',
            '-evalue',
            '1e-10',
            '-word_size',
            '11',
            '-dust',
            'no',
            '-num_threads',
            str(threads),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return hits


def flank_pairs(hits_path):
    """{'U': pairs, 'D': pairs}: locus index pairs (lower first) whose same-side flanks match."""
    pairs = {'U': set(), 'D': set()}
    with open(hits_path) as handle:
        for line in handle:
            query, subject, pident, length = line.rstrip('\n').split('\t')
            query_index, query_side = query.split('|')
            subject_index, subject_side = subject.split('|')
            query_index = int(query_index[1:])
            subject_index = int(subject_index[1:])

            # only same-side hits between different loci above both cutoffs count
            if query_side != subject_side or query_index == subject_index:
                continue
            if float(pident) < MIN_PIDENT or int(length) < MIN_LEN:
                continue
            pairs[query_side].add(tuple(sorted((query_index, subject_index))))
    return pairs


def union_find_groups(node_count, edges):
    """Member lists of the connected components, in order of their first member."""
    parent = list(range(node_count))

    # root lookup with path halving
    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for first, second in edges:
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[root_second] = root_first

    groups = {}
    for node in range(node_count):
        groups.setdefault(find(node), []).append(node)
    return list(groups.values())


def dedup_groups(genome_id, rows, contigs, scratch, threads):
    """Duplicate-group id per locus, shared by loci whose two flanks both match."""
    prefix = f'{scratch}/{slug(genome_id)}'
    fasta_path = f'{prefix}_flanks.fa'
    edges = []
    if len(rows) > 1 and write_flanks(fasta_path, rows, contigs):
        pairs = flank_pairs(blast_flanks(fasta_path, prefix, threads))

        # both flanks required: multigene families often share >=95% identity
        # on one side only (Ling 2007)
        edges = sorted(pairs['U'] & pairs['D'])
        one_side_only = len((pairs['U'] | pairs['D']) - (pairs['U'] & pairs['D']))
        print(
            f'    [dedup] both-flank edges: {len(edges)}; '
            f'one-side-only (potential under-collapse): {one_side_only}',
            flush=True,
        )

    group_of = [0] * len(rows)
    for group_id, members in enumerate(union_find_groups(len(rows), edges)):
        for member in members:
            group_of[member] = group_id
    return group_of


def choose_representatives(rows, group_of):
    """Row index of each row's group representative: standard tier, longest intron, coordinates."""
    members = defaultdict(list)
    for i, group in enumerate(group_of):
        members[group].append(i)
    representative = {
        group: min(
            indices,
            key=lambda i: (
                rows[i]['tier'] != 'standard',
                -int(rows[i]['intron_len']),
                rows[i]['seqid'],
                int(rows[i]['start']),
                rows[i]['strand'],
            ),
        )
        for group, indices in members.items()
    }
    return [representative[group] for group in group_of]
