"""Map one genome's catalogued introns onto exact CDS gaps and translate the host proteins."""

import hashlib
import json
from collections import Counter

from common import WORK, dump, header, parse_stage, read, sha, task_genome, use_scratch, write
from hostlib import exact_gaps, load_cds, read_subset
from telotron_ortholog_align import protein_of
from _common import rc

MIN_PROTEIN_AA = 20
MAX_X_FRACTION = 0.1  # unknown residues allowed in a usable translation
AUDIT_FIELDS = [
    'locus_id',
    'genome_id',
    'tier',
    'status',
    'matching_CDS_models',
    'gene_ids',
    'chosen_tx',
    'protein_id',
    'protein_length',
    'intron_residue',
    'intron_phase',
    'alternative_models_json',
]


def candidate_models(row, models, gaps, contig):
    """Every CDS model whose exact gap is this intron, translated and checked for use."""
    start, end, strand = int(row['start']), int(row['end']), row['strand']
    if strand not in ('+', '-'):
        return []
    candidates = []
    for tx_id, residue, phase in gaps.get((row['seqid'], start, end, strand), []):
        model = models[tx_id]
        protein = protein_of(model, contig)

        # usable: long enough, no internal stop, few X, junction inside the protein
        usable = (
            len(protein) >= MIN_PROTEIN_AA
            and '*' not in protein
            and protein.count('X') <= MAX_X_FRACTION * len(protein)
            and 0 <= residue <= len(protein)
        )
        candidates.append(
            dict(
                tx_id=tx_id,
                protein_id=model['prot'],
                protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
                protein_length=len(protein),
                intron_residue=residue,
                intron_phase=phase,
                genes=model['genes'],
                usable=usable,
                sequence=protein,
            )
        )
    return candidates


def mapping_status(strand, candidates, genes):
    """Status named by the first failing check, or 'mapped' when all pass."""
    if strand not in ('+', '-'):
        return 'unresolved_strand'
    if not candidates:
        return 'no_exact_CDS_gap'
    if not any(c['usable'] for c in candidates):
        return 'invalid_CDS_translation'
    if len(genes) > 1:
        return 'ambiguous_overlapping_genes'
    return 'mapped'


def main():
    use_scratch()
    stage = parse_stage(__doc__)
    meta = task_genome()
    genome_id = meta['genome_id']
    out_dir = WORK / stage / 'host' / genome_id
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort = WORK / 'inputs' / f'{stage}.tsv'
    rows = [row for row in read(cohort) if row['genome_id'] == genome_id]

    models = load_cds(meta['gff'])
    gaps = exact_gaps(models)
    contigs = read_subset(meta['fasta'], {row['seqid'] for row in rows})

    audit = []
    mapped = []
    proteins = {}
    for row in rows:
        # intron DNA on the locus strand must match the catalogue
        start, end, strand = int(row['start']), int(row['end']), row['strand']
        dna = contigs[row['seqid']][start - 1 : end]
        if strand == '-':
            dna = rc(dna)
        assert len(dna) == int(row['intron_len'])
        if strand in ('+', '-'):
            assert hashlib.sha256(dna.encode()).hexdigest() == row['intron_sha256'], row['locus_id']

        candidates = candidate_models(row, models, gaps, contigs[row['seqid']])
        genes = {gene for c in candidates for gene in c['genes']}
        status = mapping_status(strand, candidates, genes)
        record = dict(
            locus_id=row['locus_id'],
            genome_id=genome_id,
            tier=stage,
            status=status,
            matching_CDS_models=len(candidates),
            gene_ids=';'.join(sorted(genes)),
            chosen_tx='',
            protein_id='',
            protein_length='',
            intron_residue='',
            intron_phase='',
            alternative_models_json=json.dumps(
                [{k: v for k, v in c.items() if k != 'sequence'} for c in candidates]
            ),
        )

        # mapped: the longest usable protein, ties broken by transcript ID
        if status == 'mapped':
            usable = [c for c in candidates if c['usable']]
            chosen = min(usable, key=lambda c: (-c['protein_length'], c['tx_id']))
            protein_id = chosen['protein_id']
            protein = chosen['sequence']
            collision = protein_id in proteins and proteins[protein_id] != protein
            assert not collision, (genome_id, protein_id, 'protein identifier collision')
            proteins[protein_id] = protein
            record.update(
                chosen_tx=chosen['tx_id'],
                protein_id=protein_id,
                protein_length=len(protein),
                intron_residue=chosen['intron_residue'],
                intron_phase=chosen['intron_phase'],
            )
            mapped.append(dict(row, tx_id=chosen['tx_id'], gene_id=';'.join(sorted(genes))))
        audit.append(record)

    # the mapped catalogue keeps the cohort's columns, plus the host transcript and genes
    write(out_dir / 'host_audit.tsv', audit, AUDIT_FIELDS)
    fields = header(cohort)
    fields += [name for name in ('tx_id', 'gene_id') if name not in fields]
    write(out_dir / 'mapped_catalogue.tsv', mapped, fields)

    statuses = Counter(record['status'] for record in audit)
    dump(
        out_dir / 'complete.json',
        dict(
            complete=True,
            stage=stage,
            genome_id=genome_id,
            loci=len(rows),
            mapped=len(mapped),
            proteins=len(proteins),
            status=dict(statuses),
            cohort_sha256=sha(cohort),
            scratch=use_scratch(),
            outputs={
                p.name: sha(p)
                for p in out_dir.iterdir()
                if p.is_file() and p.name != 'complete.json'
            },
        ),
    )
    print('HOST COMPLETE', stage, genome_id, len(rows), len(mapped), statuses, flush=True)


if __name__ == '__main__':
    main()
