"""FASTA, GFF, gene-mask and interval helpers shared by the ITS scans."""

import bisect, collections, csv, gzip
from urllib.parse import unquote

# feature types in the whole-model gene mask, besides every type whose name contains 'RNA'
MODEL_TYPES = ('gene', 'pseudogene', 'mRNA', 'transcript', 'exon', 'CDS')

Feature = collections.namedtuple(
    'Feature', 'contig kind start end strand phase feature_id parents attributes'
)


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_tsv(path):
    """Rows of a plain or gzipped tab-separated file, as dicts."""
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def rows_by_genome(path):
    """Rows of a tab-separated file keyed by their genome_id column."""
    return {row['genome_id']: row for row in read_tsv(path)}


def write_tsv(path, rows):
    """Write dict rows to a tab-separated file, gzipped when the name ends in .gz."""
    rows = list(rows)
    assert rows, path
    if str(path).endswith('.gz'):
        handle = gzip.open(path, 'wt', compresslevel=3)
    else:
        handle = open(path, 'w')
    with handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'
        )
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    """(name, uppercase sequence) per FASTA record, named by the header's first word."""
    name = None
    parts = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                # a new header closes the previous record
                if name is not None:
                    yield name, ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                if name is None and line.strip():
                    raise ValueError(f'Sequence before the first FASTA header in {path}')
                parts.append(line.strip())

    if name is not None:
        yield name, ''.join(parts).upper()


def count_acgt(sequence):
    """Number of A, C, G and T bases."""
    return sum(sequence.count(base) for base in 'ACGT')


def gff_features(path):
    """Feature rows of a GFF3 file, with 0-based half-open coordinates."""
    with open_text(path) as handle:
        for line in handle:
            # stop at embedded sequence, skip comments and blank lines
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Malformed annotation row in {path}: {line!r}')
            # GFF 1-based inclusive to 0-based half-open
            start = int(fields[3]) - 1
            end = int(fields[4])
            if not 0 <= start < end:
                raise ValueError(f'Invalid feature coordinates in {path}: {line!r}')

            # column 9 key=value pairs, Parent split on ',' before unquoting
            attributes = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            parents = tuple(
                unquote(name) for name in attributes.get('Parent', '').split(',') if name
            )
            yield Feature(
                contig=fields[0],
                kind=fields[2],
                start=start,
                end=end,
                strand=fields[6],
                phase=fields[7],
                feature_id=unquote(attributes.get('ID', '')),
                parents=parents,
                attributes=attributes,
            )


def gene_mask(path):
    """Merged whole-model genic intervals per contig, each parent's full span (introns) included."""
    intervals = collections.defaultdict(list)
    # full span of each exon/CDS parent, keyed by (contig, strand, parent)
    spans = {}

    for row in gff_features(path):
        if row.kind in MODEL_TYPES or 'RNA' in row.kind:
            intervals[row.contig].append((row.start, row.end))

        # widen each parent's span to cover all of its exon/CDS rows
        if row.kind in ('exon', 'CDS'):
            for parent in row.parents:
                key = row.contig, row.strand, parent
                low, high = spans.get(key, (row.start, row.end))
                spans[key] = min(low, row.start), max(high, row.end)

    for (contig, strand, parent), span in spans.items():
        intervals[contig].append(span)
    return {contig: merge(pieces) for contig, pieces in intervals.items()}


def merge(intervals):
    """Sorted union of half-open intervals, touching ones joined."""
    result = []
    for start, end in sorted(set(intervals)):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def overlaps(start, end, mask):
    """Whether [start, end) overlaps a merged, sorted mask."""
    # the last mask interval starting before end must end after start
    i = bisect.bisect_left(mask, (end,))
    return i > 0 and mask[i - 1][1] > start


def subtract(intervals, mask):
    """Parts of each interval outside a merged, sorted mask."""
    result = []
    for start, end in intervals:
        # only the last mask interval starting at or before start can contain start
        i = max(0, bisect.bisect_right(mask, (start, 10**30)) - 1)
        cursor = start

        # walk the overlapping mask intervals, keeping the uncovered pieces
        while i < len(mask) and mask[i][0] < end:
            low, high = mask[i]
            if high > cursor:
                if low > cursor:
                    result.append((cursor, min(low, end)))
                cursor = max(cursor, high)
            i += 1

        # uncovered tail after the last mask interval
        if cursor < end:
            result.append((cursor, end))
    return result
