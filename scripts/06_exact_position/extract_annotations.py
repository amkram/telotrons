"""Extract the annotated coding models, with DNA, spanning each focal and comparator site."""

import json
from collections import defaultdict
from pathlib import Path

from exact_position import (
    ORTHOLOGY,
    OUT,
    RESOLVED_STATES,
    GenomeSequence,
    check_scratch,
    genome_files,
    read_json_gz,
    sha,
    write_json_gz,
)
from _common import rc, read_aln
from hostlib import open_text

FLANK_BP = 30  # genomic bases kept either side of each requested site
INTRON_DNA_WINDOW_BP = 100  # model introns this close to the site keep their DNA


def focal_proteins(loci):
    """Focal protein of every locus, read back from its protein MSA."""
    proteins = {}
    for locus in loci:
        genome_id = locus['genome_id']
        directory = ORTHOLOGY / 'recall/primary/orthology' / genome_id / 'protein_msa'
        msa = read_aln(directory / f"{genome_id}__{locus['protein_id']}.aln.fasta")
        protein = ''.join(seq for name, seq in msa if name.split('|')[0] == genome_id)
        protein = protein.replace('-', '')
        assert protein, locus['locus_id']
        proteins[locus['locus_id']] = protein
    return proteins


def site_requests(loci):
    """Per genome, one request per focal intron and per resolved comparator site."""
    requests = defaultdict(list)
    for locus in loci:
        tid = locus['locus_id']
        catalogue = locus['catalogue']
        requests[locus['genome_id']].append(
            dict(
                key=tid + '|FOCAL',
                locus_id=tid,
                role='focal',
                seqid=catalogue['seqid'],
                strand=catalogue['strand'],
                start=int(catalogue['start']),
                end=int(catalogue['end']),
                expected_intron=locus['focal_intron'],
            )
        )

        # comparators with a resolved state and a miniprot hit
        for pair in locus['pairs']:
            if not pair['state'] or pair['state']['state'] not in RESOLVED_STATES:
                continue
            hit = pair.get('miniprot_hit')
            if not hit:
                continue

            # comparator site: the recovered intron, else the predicted one, else the codon
            # of the focal residue
            predicted = pair.get('selected_predicted_intron')
            if pair.get('correction') and pair['correction'].get('recovered_intron'):
                predicted = pair['correction']['recovered_intron']
            if predicted:
                start, end = predicted['start'], predicted['end']
            else:
                residue = locus['intron_residue']
                exons = [cds for cds in hit['cds_q'] if cds[0] <= residue <= cds[1]]
                if not exons:
                    continue
                query_start, _, exon_start, exon_end, _ = exons[0]

                # approximate centre, used only to select models (position comes from the alignment)
                if hit['strand'] == '+':
                    start = exon_start + (residue - query_start) * 3
                else:
                    start = exon_end - (residue - query_start) * 3
                end = start
            requests[pair['ortholog']].append(
                dict(
                    key=tid + '|' + pair['ortholog'],
                    locus_id=tid,
                    role=pair['state']['state'],
                    seqid=hit['contig'],
                    strand=hit['strand'],
                    start=start,
                    end=end,
                    expected_intron=pair['intron'],
                )
            )
    return requests


def cds_models(gff, surfaces):
    """(contig, strand) -> [(parent, CDS rows, span start, span end)] on the requested surfaces."""
    models = defaultdict(list)
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'CDS' or (fields[0], fields[6]) not in surfaces:
                continue
            attributes = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            for parent in attributes.get('Parent', attributes.get('ID', '')).split(','):
                models[(fields[0], fields[6], parent)].append((fields, attributes))

    by_surface = defaultdict(list)
    for (contig, strand, parent), items in models.items():
        low = min(int(fields[3]) for fields, _ in items)
        high = max(int(fields[4]) for fields, _ in items)
        by_surface[(contig, strand)].append((parent, items, low, high))
    return by_surface


def spanning_models(request, by_surface, genome):
    """Models spanning the site: CDS in transcript order and the introns between them."""
    contig, strand = request['seqid'], request['strand']
    start, end = request['start'], request['end']
    models = []
    for parent, items, low, high in by_surface[(contig, strand)]:
        if not (low < start and high > end):
            continue
        items = sorted(items, key=lambda item: int(item[0][3]), reverse=strand == '-')
        cds = []
        for fields, _ in items:
            exon_start, exon_end = int(fields[3]), int(fields[4])
            seq = genome.fetch_stranded(contig, exon_start - 1, exon_end, strand)
            cds.append(dict(start=exon_start, end=exon_end, phase=fields[7], sequence=seq))

        introns = []
        for first, second in zip(cds, cds[1:]):
            if strand == '-':
                intron_start, intron_end = second['end'] + 1, first['start'] - 1
            else:
                intron_start, intron_end = first['end'] + 1, second['start'] - 1
            intron = dict(start=intron_start, end=intron_end, bp=intron_end - intron_start + 1)

            # DNA only for introns within INTRON_DNA_WINDOW_BP of the site, to bound output size
            near = intron_start <= end + INTRON_DNA_WINDOW_BP
            near = near and intron_end >= start - INTRON_DNA_WINDOW_BP
            if near:
                intron['sequence'] = genome.fetch_stranded(
                    contig, intron_start - 1, intron_end, strand
                )
            introns.append(intron)

        models.append(
            dict(
                parent=parent,
                protein_id=items[0][1].get('protein_id', ''),
                product=items[0][1].get('product', ''),
                cds=cds,
                introns=introns,
            )
        )
    return models


def main():
    check_scratch()
    evidence = read_json_gz('evidence.json.gz')
    proteins = focal_proteins(evidence['loci'])
    requests = site_requests(evidence['loci'])

    rows = []
    sources = {}
    for genome_id, genome_requests in requests.items():
        fasta, gff = genome_files(genome_id)
        assert fasta, genome_id
        genome = GenomeSequence(fasta, {request['seqid'] for request in genome_requests})
        by_surface = defaultdict(list)
        if gff:
            sources[gff] = sha(Path(gff))
            surfaces = {(request['seqid'], request['strand']) for request in genome_requests}
            by_surface = cds_models(gff, surfaces)

        # stored intron and FLANK_BP flanks on the request strand, and the models spanning it
        for request in genome_requests:
            contig, strand = request['seqid'], request['strand']
            start, end = request['start'], request['end']
            genomic = genome.fetch_stranded(contig, start - 1, end, strand)
            left = genome.fetch(contig, start - 1 - FLANK_BP, start - 1)
            right = genome.fetch(contig, end, end + FLANK_BP)
            if strand == '-':
                left, right = rc(right), rc(left)
            expected = request['expected_intron']
            rows.append(
                dict(
                    **request,
                    genome_id=genome_id,
                    focal_protein=proteins[request['locus_id']],
                    annotation_available=bool(gff),
                    genomic_sequence=genomic,
                    genomic_matches_stored_intron=genomic == expected if expected else None,
                    flank5=left,
                    flank3=right,
                    models=spanning_models(request, by_surface, genome),
                    reference=genome.provenance,
                )
            )

        n_models = sum(len(row['models']) for row in rows if row['genome_id'] == genome_id)
        print(genome_id, len(genome_requests), 'requests', n_models, 'models', flush=True)
        genome.close()

    write_json_gz('annotations.json.gz', dict(rows=rows, sources=sources))
    complete = dict(
        complete=True,
        requests=len(rows),
        genomes=len(requests),
        evidence_sha256=sha(OUT / 'evidence.json.gz'),
        annotations_sha256=sha(OUT / 'annotations.json.gz'),
        script_sha256=sha(__file__),
    )
    (OUT / 'annotation_extraction_complete.json').write_text(json.dumps(complete, indent=2) + '\n')
    print('COMPLETE', len(rows), 'annotation requests', flush=True)


if __name__ == '__main__':
    main()
