"""Extract the miniprot-predicted model of each comparator pair as a non-independent tier."""

import json
from collections import defaultdict

from exact_position import OUT, GenomeSequence, check_scratch, genome_files, read_json_gz, sha
from exact_position import write_json_gz


def read_predictions(path, protein_ids):
    """mRNA rows of the wanted proteins, and their CDS rows, from a miniprot GFF."""
    mrna = {}
    cds = defaultdict(list)
    for line in open(path):
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            continue
        attributes = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
        if fields[2] == 'mRNA' and attributes['Target'].split()[0] in protein_ids:
            mrna[attributes['ID']] = (fields, attributes)
        elif fields[2] == 'CDS' and attributes.get('Parent') in mrna:
            cds[attributes['Parent']].append((fields, attributes))
    return mrna, cds


def matching_prediction(locus, pair, mrna, cds):
    """The single Rank=1 mRNA that reproduces the stored hit: (mRNA ID, CDS rows)."""
    hit = pair['miniprot_hit']
    found = []
    for mrna_id, (fields, attributes) in mrna.items():
        target = attributes['Target'].split()
        same_hit = (
            target[0] == locus['protein_id']
            and int(attributes['Rank']) == 1
            and fields[0] == hit['contig']
            and fields[6] == hit['strand']
            and int(target[1]) == hit['qmin']
            and int(target[2]) == hit['qmax']
        )
        if not same_hit:
            continue
        segments = sorted((int(row[3]), int(row[4])) for row, _ in cds[mrna_id])
        if segments == sorted(tuple(seg) for seg in hit['segs']):
            found.append((mrna_id, cds[mrna_id]))
    assert len(found) == 1, (locus['locus_id'], pair['ortholog'], len(found))
    return found[0]


def predicted_model(mrna_id, items, protein_id, contig, strand, genome):
    """Predicted CDS in transcript order with their DNA, and the introns between them."""
    exons = []
    for fields, attributes in sorted(
        items, key=lambda item: int(item[0][3]), reverse=strand == '-'
    ):
        start, end = int(fields[3]), int(fields[4])
        exons.append(
            dict(
                start=start,
                end=end,
                phase=fields[7],
                sequence=genome.fetch_stranded(contig, start - 1, end, strand),
                query=attributes['Target'],
                attributes=attributes,
            )
        )

    introns = []
    for left, right in zip(exons, exons[1:]):
        if strand == '-':
            start, end = right['end'] + 1, left['start'] - 1
        else:
            start, end = left['end'] + 1, right['start'] - 1
        introns.append(dict(start=start, end=end, bp=end - start + 1))
    return dict(
        parent=mrna_id,
        protein_id=protein_id,
        product='miniprot prediction (not independent annotation)',
        cds=exons,
        introns=introns,
    )


def main():
    check_scratch()
    evidence = read_json_gz('evidence.json.gz')
    annotations = {row['key']: row for row in read_json_gz('annotations.json.gz')['rows']}

    # evidence pairs grouped by miniprot GFF, whose hash must match the evidence
    pairs_by_gff = defaultdict(list)
    for locus in evidence['loci']:
        for pair in locus['pairs']:
            pairs_by_gff[pair['miniprot_source']].append((locus, pair))

    by_genome = defaultdict(list)
    sources = {}
    for path, pairs in pairs_by_gff.items():
        sources[path] = sha(path)
        assert sources[path] == evidence['source_hashes'][path]
        mrna, cds = read_predictions(path, {locus['protein_id'] for locus, _ in pairs})
        for locus, pair in pairs:
            by_genome[pair['ortholog']].append(
                (locus, pair, matching_prediction(locus, pair, mrna, cds))
            )

    # per comparator genome: the annotation request, with the prediction as its only model
    rows = []
    for genome_id, items in by_genome.items():
        fasta, _ = genome_files(genome_id)
        assert fasta
        genome = GenomeSequence(fasta, {pair['miniprot_hit']['contig'] for _, pair, _ in items})
        for locus, pair, (mrna_id, cds_rows) in items:
            hit = pair['miniprot_hit']
            request = annotations[locus['locus_id'] + '|' + genome_id]
            row = {k: v for k, v in request.items() if k != 'models'}
            row['models'] = [
                predicted_model(
                    mrna_id, cds_rows, locus['protein_id'], hit['contig'], hit['strand'], genome
                )
            ]
            row['annotation_available'] = True
            row['evidence_type'] = 'miniprot_prediction'
            row['source'] = pair['miniprot_source']
            rows.append(row)
        genome.close()
        print(genome_id, len(items), 'predicted models', flush=True)

    write_json_gz('predictions.json.gz', dict(rows=rows, sources=sources))
    complete = dict(
        complete=True,
        requests=len(rows),
        output_sha256=sha(OUT / 'predictions.json.gz'),
        script_sha256=sha(__file__),
        evidence_sha256=sha(OUT / 'evidence.json.gz'),
    )
    (OUT / 'prediction_extraction_complete.json').write_text(json.dumps(complete, indent=2) + '\n')


if __name__ == '__main__':
    main()
