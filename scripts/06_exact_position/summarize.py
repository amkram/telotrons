"""Per-pair, per-locus, per-lineage and per-genome exact position results and a summary."""

import json
from collections import Counter

from exact_position import ORTHOLOGY, OUT, primary_closest, read_json_gz, sha, write_tsv

EXACT = 'exact_position_phase'
EXACT_NEAR_GAP = 'exact_direct_site_near_alignment_gap'


def pair_row(locus, pair, focal, check):
    """One comparator pair: focal and comparator statuses and the stored-boundary sites."""
    sites = [
        dict(parent=m['parent'], phase_source=m.get('phase_source'), qualified=m['qualified'], **z)
        for m in check['models']
        for z in m.get('introns', [])
        if z['stored_boundary']
    ]
    return dict(
        locus_id=locus['locus_id'],
        genome_id=locus['genome_id'],
        clade=locus['closest']['clade'],
        comparator=pair['ortholog'],
        identity=pair['state']['protein_ident'],
        focal_status=focal['status'],
        comparator_status=check['status'],
        original_annotated_phase=check['original_annotated_phase'],
        inferred_phase_only=check['inferred_phase_only'],
        raw_query_offset=(pair.get('selected_predicted_intron') or {}).get('delta', ''),
        annotation_sites=json.dumps(sites, separators=(',', ':')),
    )


def locus_status(focal, checks):
    """(status, corroborated, direct-site sensitivity, exact checks), first matching rule wins."""
    exact = [check for check in checks if check['status'] == EXACT]
    direct = [check for check in checks if check['status'] in (EXACT, EXACT_NEAR_GAP)]
    focal_ok = focal['status'] == EXACT
    supported = focal_ok and bool(exact)

    # sensitivity: a direct site next to an alignment gap also counts
    direct_supported = focal['status'] in (EXACT, EXACT_NEAR_GAP) and bool(direct)
    if supported:
        status = 'corroborated_exact_position_phase'
    elif direct_supported:
        status = 'exact_direct_site_but_nearby_gap'
    elif not focal_ok:
        status = 'focal_boundary_or_phase_unresolved'
    elif any(check['status'] == 'different_position_or_phase' for check in checks):
        status = 'mapped_position_or_phase_mismatch'
    elif any(check['status'] == 'exact_site_different_genomic_boundary' for check in checks):
        status = 'genomic_boundary_discrepancy'
    else:
        status = 'comparator_annotation_or_alignment_unresolved'
    return status, supported, direct_supported, exact


def percent_absent(absent, ordinary):
    """Closest-absent loci as a percentage of absent plus ordinary, None when both are 0."""
    return 100 * absent / (absent + ordinary) if absent + ordinary else None


def main():
    # evidence and model reviews, indexed by request key
    evidence = read_json_gz('evidence.json.gz')
    checks = {row['key']: row for row in read_json_gz('model_review.json.gz')['rows']}

    # per locus: one pair row per comparator, then a locus status
    loci = []
    pairs = []
    for locus in evidence['loci']:
        tid = locus['locus_id']
        focal = checks[tid + '|FOCAL']
        locus_checks = [checks[tid + '|' + pair['ortholog']] for pair in locus['pairs']]
        for pair, check in zip(locus['pairs'], locus_checks):
            pairs.append(pair_row(locus, pair, focal, check))

        status, supported, direct_supported, exact = locus_status(focal, locus_checks)
        loci.append(
            dict(
                locus_id=tid,
                genome_id=locus['genome_id'],
                organism=locus['organism'],
                clade=locus['closest']['clade'],
                focal_residue=locus['intron_residue'],
                focal_phase=locus['intron_phase'],
                focal_status=focal['status'],
                status=status,
                corroborated=supported,
                direct_site_sensitivity=direct_supported,
                tied_comparators=len(locus_checks),
                exact_comparators=len(exact),
                all_tied_comparators_exact=bool(exact) and len(exact) == len(locus_checks),
                both_phases_annotated=(
                    supported
                    and focal['original_annotated_phase']
                    and any(check['original_annotated_phase'] for check in exact)
                ),
                comparator_statuses=';'.join(check['status'] for check in locus_checks),
                rna_exact_junction_reads=(locus.get('rna') or {}).get(
                    'exact_junction_reads', 'not_assayed'
                ),
            )
        )

    write_tsv('locus_results.tsv', loci)
    write_tsv('pair_results.tsv', pairs)

    # reviewed loci must be exactly the closest-ORDINARY primary loci
    baseline = primary_closest()
    ordinary = {row['locus_id'] for row in baseline if row['closest_state'] == 'ORDINARY'}
    assert ordinary == {row['locus_id'] for row in loci}

    # per-lineage and per-genome tables before and after the audit
    lineages = []
    for clade in sorted({row['clade'] for row in baseline}):
        before = [row for row in baseline if row['clade'] == clade]
        reviewed = [row for row in loci if row['clade'] == clade]
        absent = sum(row['closest_state'] == 'NO_INTRON' for row in before)
        corroborated = sum(row['corroborated'] for row in reviewed)
        lineages.append(
            dict(
                clade=clade,
                primary_loci=len(before),
                closest_absent_unchanged=absent,
                ordinary_before=len(reviewed),
                corroborated_exact=corroborated,
                unconfirmed_ordinary=len(reviewed) - corroborated,
                original_other=len(before) - absent - len(reviewed),
                other_after=len(before) - absent - corroborated,
                conditional_absent_before=percent_absent(absent, len(reviewed)),
                conditional_absent_after=percent_absent(absent, corroborated),
                near_gap_sensitivity_additions=sum(
                    row['direct_site_sensitivity'] and not row['corroborated'] for row in reviewed
                ),
                status_counts=json.dumps(
                    dict(Counter(r['status'] for r in reviewed)), sort_keys=True
                ),
            )
        )
    write_tsv('lineage_changes.tsv', lineages)

    genomes = []
    for genome_id in sorted({row['genome_id'] for row in baseline}):
        before = [row for row in baseline if row['genome_id'] == genome_id]
        reviewed = [row for row in loci if row['genome_id'] == genome_id]
        absent = sum(row['closest_state'] == 'NO_INTRON' for row in before)
        corroborated = sum(row['corroborated'] for row in reviewed)
        genomes.append(
            dict(
                genome_id=genome_id,
                primary_loci=len(before),
                closest_absent_unchanged=absent,
                ordinary_before=len(reviewed),
                corroborated_exact=corroborated,
                unconfirmed_ordinary=len(reviewed) - corroborated,
                other_after=len(before) - absent - corroborated,
                status_counts=json.dumps(
                    dict(Counter(r['status'] for r in reviewed)), sort_keys=True
                ),
            )
        )
    write_tsv('genome_changes.tsv', genomes)

    # overall summary: closest-state totals before the audit and two sensitivities after it
    absent = sum(row['closest_state'] == 'NO_INTRON' for row in baseline)
    n_ordinary = len(loci)
    other = len(baseline) - absent - n_ordinary
    confirmed = sum(row['corroborated'] for row in loci)
    direct = sum(row['direct_site_sensitivity'] for row in loci)

    def totals(ordinary_count):
        return dict(
            absent=absent,
            ordinary=ordinary_count,
            other=other + n_ordinary - ordinary_count,
            conditional_absent_percent=percent_absent(absent, ordinary_count),
        )

    inputs = (
        OUT / 'evidence.json.gz',
        OUT / 'annotations.json.gz',
        OUT / 'model_review.json.gz',
        ORTHOLOGY / 'results/closest_comparators.tsv',
    )
    summary = dict(
        complete=True,
        loci=len(loci),
        pairs=len(pairs),
        locus_statuses=dict(Counter(row['status'] for row in loci)),
        pair_statuses=dict(Counter(row['comparator_status'] for row in pairs)),
        focal_statuses=dict(Counter(row['focal_status'] for row in loci)),
        confirmed=confirmed,
        both_phases_annotated=sum(row['both_phases_annotated'] for row in loci),
        confirmed_with_all_ties_exact=sum(
            row['corroborated'] and row['all_tied_comparators_exact'] for row in loci
        ),
        direct_site_sensitivity=direct,
        ordinary_unconfirmed=n_ordinary - confirmed,
        original=totals(n_ordinary),
        conservative_sensitivity=totals(confirmed),
        direct_site_sensitivity_counts=totals(direct),
        input_sha256={str(path): sha(path) for path in inputs},
        script_sha256=sha(__file__),
        interpretation=(
            'No new absent calls. Unconfirmed ordinary loci become unresolved in this audit '
            'sensitivity, not de novo gains. Unchanged absent calls and other composition states '
            'were not revalidated. Production release unchanged.'
        ),
    )
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
