"""CDS models from a GFF, their exact intron gaps, and selected contigs from a FASTA."""

import gzip
from collections import defaultdict
from urllib.parse import unquote


def open_text(path):
    """Open a plain or gzipped text file."""
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def gff_attributes(column):
    """GFF column 9 as a dict, URL-decoded."""
    pairs = (item.split('=', 1) for item in column.split(';') if '=' in item)
    return {key: unquote(value) for key, value in pairs}


def load_cds(gff):
    """Transcript ID -> CDS model (strand, protein ID, contig, segments, phases, genes)."""
    cds = {}
    parents = {}
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            attributes = gff_attributes(fields[8])

            # every feature's parents, for the gene lookup below
            feature_id = attributes.get('ID')
            if feature_id:
                parents[feature_id] = attributes.get('Parent', '').split(',')
            if fields[2] != 'CDS':
                continue

            for tid in attributes.get('Parent', '').split(','):
                if not tid:
                    continue
                # the first CDS row creates the model and later rows must agree on strand and contig
                model = cds.setdefault(
                    tid,
                    dict(
                        strand=fields[6],
                        prot=attributes.get('protein_id', tid),
                        seqid=fields[0],
                        segs=[],
                        phase={},
                    ),
                )
                assert model['strand'] == fields[6] and model['seqid'] == fields[0], (
                    gff,
                    tid,
                    'inconsistent CDS parent',
                )
                segment = (int(fields[3]), int(fields[4]))
                if segment not in model['segs']:
                    model['segs'].append(segment)
                model['phase'][segment] = int(fields[7]) if fields[7] in ('0', '1', '2') else 0

    # a transcript's genes are its parents
    for tid, model in cds.items():
        model['genes'] = [gene for gene in parents.get(tid, []) if gene]
    return cds


def exact_gaps(cds):
    """1-based closed CDS gap (seqid, start, end, strand) -> [(transcript, residue, phase)]."""
    gaps = defaultdict(list)
    for tid, model in cds.items():
        minus = model['strand'] == '-'
        segments = sorted(model['segs'], reverse=minus)

        # coding bases before each gap, less the start phase
        coding_bp = -model['phase'][segments[0]]
        for upstream, downstream in zip(segments, segments[1:]):
            coding_bp += upstream[1] - upstream[0] + 1
            if minus:
                gap = (downstream[1] + 1, upstream[0] - 1)
            else:
                gap = (upstream[1] + 1, downstream[0] - 1)
            gaps[(model['seqid'], *gap, model['strand'])].append((tid, *divmod(coding_bp, 3)))
    return gaps


def read_subset(path, wanted):
    """Upper-case sequences of the wanted contigs, all of which must be present."""
    sequences = {}
    name = None
    parts = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name in wanted:
                    sequences[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            elif name in wanted:
                parts.append(line.strip())
        if name in wanted:
            sequences[name] = ''.join(parts).upper()
    assert wanted <= set(sequences), (path, wanted - set(sequences))
    return sequences
