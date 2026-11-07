"""Annotated introns from a GFF, and FASTA records, for the single-motif scan."""

import collections
import gzip
from urllib.parse import unquote

# features that are pieces of a transcript rather than its parents
SEGMENTS = ('exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR')


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def attributes(text):
    """ID and Parent IDs of a GFF attribute column."""
    raw = dict(item.split('=', 1) for item in text.split(';') if '=' in item)
    # split Parent before unquoting: IDs may hold escaped commas
    parents = tuple(unquote(parent) for parent in raw.get('Parent', '').split(',') if parent)
    return unquote(raw.get('ID', '')), parents


def merge(intervals):
    """Union of half-open intervals, sorted."""
    merged = []
    for start, end in sorted(set(intervals)):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def gaps(intervals):
    """Gaps between the merged blocks of one model: its introns."""
    blocks = merge(intervals)
    return [(left[1], right[0]) for left, right in zip(blocks, blocks[1:]) if left[1] < right[0]]


def annotated_introns(path):
    """(contig, start0, end, strand) of every intron of a stranded annotated model."""
    exons = collections.defaultdict(list)
    cds = collections.defaultdict(list)
    feature_type = {}
    child_types = []
    skipped = 0

    # one pass over the GFF: feature types, parent links and exon/CDS segments
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Malformed annotation: {path}')
            contig, kind, strand = fields[0], fields[2], fields[6]
            start, end = int(fields[3]) - 1, int(fields[4])
            # skip records with start >= end (organellar trans-splicing)
            if not 0 <= start < end:
                skipped += 1
                continue
            feature_id, parents = attributes(fields[8])
            if feature_id:
                key = (contig, strand, feature_id)
                if key in feature_type and feature_type[key] != kind:
                    raise ValueError(f'ID used for different feature types: {key}')
                feature_type[key] = kind
            if feature_id and parents and kind not in SEGMENTS:
                child_types.extend(((contig, strand, parent), kind) for parent in parents)
            if kind in ('exon', 'CDS'):
                segments = exons if kind == 'exon' else cds
                for parent in parents:
                    segments[contig, strand, parent].append((start, end))
    if skipped:
        print(f'annotation: skipped {skipped} records with start>=end in {path}', flush=True)

    # genes with transcripts: their own segments duplicate the transcripts' introns
    genes_with_transcripts = {
        parent
        for parent, kind in child_types
        if feature_type.get(parent) in ('gene', 'pseudogene')
        and (kind == 'transcript' or 'RNA' in kind)
    }

    # introns from each model's exons, or its CDS when it has no exons
    introns = set()
    for key in set(exons) | set(cds):
        contig, strand, _ = key
        if key in genes_with_transcripts or strand not in ('+', '-'):
            continue
        for start, end in gaps(exons.get(key) or cds[key]):
            introns.add((contig, start, end, strand))
    return introns


def fasta(path):
    """(name, upper-case sequence) of each FASTA record."""
    name = None
    parts = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                if name is None and line.strip():
                    raise ValueError('Sequence before FASTA header')
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts).upper()
