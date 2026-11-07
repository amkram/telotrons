"""Unique annotated introns: gaps between merged exon segments of each model, else CDS."""

from collections import defaultdict
import gzip
from urllib.parse import unquote


def merged(blocks):
    """Union of overlapping or abutting half-open blocks."""
    union = []
    for start, end in sorted(set(blocks)):
        if union and start <= union[-1][1]:
            union[-1] = (union[-1][0], max(end, union[-1][1]))
        else:
            union.append((start, end))
    return union


def parse_gff(path):
    """Feature types and parents keyed by (contig, strand, ID), and exon/CDS blocks per parent."""
    feature_type = {}
    parents = defaultdict(set)
    segments = defaultdict(lambda: defaultdict(list))
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            # stop at an embedded FASTA section, skip comments and blank lines
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('Malformed GFF row')
            contig, kind, strand = fields[0], fields[2], fields[6]
            # GFF 1-based closed to 0-based half-open
            start, end = int(fields[3]) - 1, int(fields[4])
            if not 0 <= start < end:
                raise ValueError('Invalid feature coordinates')

            # attributes: ID and the comma-separated Parent list
            attributes = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            feature_id = unquote(attributes.get('ID', ''))
            parent_ids = [
                unquote(value) for value in attributes.get('Parent', '').split(',') if value
            ]
            if feature_id:
                key = (contig, strand, feature_id)
                if key in feature_type and feature_type[key] != kind:
                    raise ValueError('Inconsistent feature type for ID')
                feature_type[key] = kind
                parents[key].update((contig, strand, parent) for parent in parent_ids)
            if kind in ('exon', 'CDS'):
                for parent in parent_ids:
                    segments[(contig, strand, parent)][kind].append((start, end))
    return feature_type, parents, segments


def derive(path):
    """Set of (contig, start0, end0, strand) introns of a GFF, 0-based half-open."""
    feature_type, parents, segments = parse_gff(path)

    # segments listed directly under a gene that also has transcripts are skipped
    genes_with_transcripts = {
        parent
        for key, key_parents in parents.items()
        for parent in key_parents
        if feature_type.get(parent) in ('gene', 'pseudogene')
        and (feature_type.get(key) == 'transcript' or 'RNA' in feature_type.get(key, ''))
    }

    # one intron per gap between consecutive merged blocks of each model
    introns = set()
    for key, parts in segments.items():
        contig, strand, _model = key
        if strand not in ('+', '-') or key in genes_with_transcripts:
            continue
        # exons define the model when present, else CDS
        blocks = merged(parts['exon'] if parts.get('exon') else parts['CDS'])
        for upstream, downstream in zip(blocks, blocks[1:]):
            introns.add((contig, upstream[1], downstream[0], strand))
    return introns
