"""Join RNA, orthology and exact-position evidence onto the catalogue (missing stays blank)."""

import collections
import json

from common import PACKAGE, REPO, WORK, read_tsv, write_tsv


def by_locus(path):
    """Rows keyed by locus_id (each locus must appear once)."""
    rows = read_tsv(path)
    keyed = {row['locus_id']: row for row in rows}
    assert len(keyed) == len(rows), path
    return keyed


def tier_summary(loci):
    """Splice-site, architecture, RNA and closest-comparator counts for one tier."""
    assayed = [locus for locus in loci if locus['rna_status'] == 'assayed']
    return dict(
        n=len(loci),
        GT_AG=sum(locus['splice_class'] == 'GT-AG' for locus in loci),
        GC_AG=sum(locus['splice_class'] == 'GC-AG' for locus in loci),
        architectures=dict(collections.Counter(locus['architecture'] for locus in loci)),
        classified=sum(locus['architecture'] != 'mixed_unclear' for locus in loci),
        RNA=dict(
            assayed=len(assayed),
            supported=sum(locus['rna_exact'] >= 1 for locus in assayed),
            support_ge3=sum(locus['rna_exact'] >= 3 for locus in assayed),
            unassayed_with_dataset=sum(locus['rna_status'] == 'not_assayed' for locus in loci),
        ),
        orthology_states=dict(
            collections.Counter(locus['closest_state'] or 'NOT_EVALUATED' for locus in loci)
        ),
    )


def orthology_summary(loci):
    """Closest-comparator outcomes of primary loci."""
    absent = sum(locus['closest_state'] == 'NO_INTRON' for locus in loci)
    ordinary = sum(locus['ordinary_exact'] for locus in loci)
    prediction_only = sum(locus['ordinary_prediction_only'] for locus in loci)

    # the reported denominator: intron absent plus annotation-supported ordinary introns
    return dict(
        absent=absent,
        ordinary_exact=ordinary,
        denominator=absent + ordinary,
        prediction_only=prediction_only,
        sensitivity_denominator=absent + ordinary + prediction_only,
    )


def main():
    catalogue = read_tsv(WORK / 'catalogue/catalogue.tsv.gz')
    rna = by_locus(WORK / '04_rna/RNA_loci.tsv')
    closest = by_locus(WORK / '05_orthology/results/closest_comparators.tsv')
    exact_position = by_locus(WORK / '06_exact_position/combined_locus_results.tsv')
    genomes_with_rna = {run['genome_id'] for run in read_tsv(REPO / 'config/rna_runs.tsv')}

    for locus in catalogue:
        locus_id = locus['locus_id']

        # RNA junction reads (no_dataset: the genome has no RNA-seq run)
        assayed = rna.get(locus_id, {}).get('assayed') == '1'
        locus['rna_exact'] = int(rna[locus_id]['exact_junction_reads']) if assayed else ''
        if assayed:
            locus['rna_status'] = 'assayed'
        elif locus['genome_id'] in genomes_with_rna:
            locus['rna_status'] = 'not_assayed'
        else:
            locus['rna_status'] = 'no_dataset'

        # closest comparator and the exact position/phase audit of ordinary introns
        comparison = closest.get(locus_id, {})
        locus['closest_state'] = comparison.get('closest_state', '')
        audit = exact_position.get(locus_id, {}).get('combined_status', 'not_audited')
        locus['exact_phase_audit'] = audit
        locus['ordinary_exact'] = int(audit == 'annotation_corroborated_exact')
        locus['ordinary_prediction_only'] = int(audit == 'prediction_only_exact')
        if locus['ordinary_exact']:
            assert locus['closest_state'] == 'ORDINARY', locus_id
        locus['closest_identity'] = comparison.get('closest_identity', '')
        locus['closest_comparators'] = comparison.get('closest_comparators', '')

    # summaries of the primary (standard, >85%) and additional (>50% to 85%) tiers
    primary = [locus for locus in catalogue if locus['tier'] == 'standard']
    additional = [locus for locus in catalogue if locus['tier'] == 'additional']
    lineages = sorted({locus['clade'] for locus in primary})
    summary = dict(
        primary=tier_summary(primary),
        additional=tier_summary(additional),
        primary_orthology=orthology_summary(primary),
        primary_orthology_by_lineage={
            lineage: orthology_summary([locus for locus in primary if locus['clade'] == lineage])
            for lineage in lineages
        },
    )

    # TSD totals over all primary loci (stage 3)
    totals = next(
        row for row in read_tsv(WORK / '03_tsd/TSD_summary.tsv') if row['group_type'] == 'all'
    )
    assert int(totals['loci']) == len(primary)
    summary['TSD'] = dict(
        loci=int(totals['loci']), TSD=int(totals['TSD']), percent=float(totals['percent'])
    )

    write_tsv(PACKAGE / 'results/catalogue_with_evidence.tsv.gz', catalogue)
    (PACKAGE / 'results/summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
