"""One status per locus from the annotation, prediction and alternative-intron reviews."""

import json
from collections import Counter, defaultdict

from exact_position import OUT, read_json_gz, read_tsv, sha, write_tsv

UNQUALIFIED_ANNOTATION = (
    'annotation_alignment_or_phase_unqualified',
    'no_overlapping_coding_model',
)


def combined_status(locus, predictions, annotations, telomeric):
    """Annotation first, then non-conflicting predictions (first matching rule wins)."""
    exact = [p for p in predictions if p['status'] == 'exact_position_phase']

    # exact predictions where the annotation was unqualified or absent
    nonconflicting = [p for p in exact if p['original_annotation_status'] in UNQUALIFIED_ANNOTATION]
    if locus['corroborated'] == 'True':
        return 'annotation_corroborated_exact', exact
    if nonconflicting:
        return 'prediction_only_exact', exact
    if exact:
        return 'prediction_annotation_conflict', exact
    if locus['locus_id'] in telomeric:
        return 'exact_alternative_is_telomeric', exact
    boundary = any(p['status'] == 'exact_site_different_genomic_boundary' for p in predictions)
    if locus['status'] == 'genomic_boundary_discrepancy' or boundary:
        return 'alternative_boundary_disagreement', exact
    if locus['status'] == 'exact_direct_site_but_nearby_gap':
        return 'exact_direct_site_near_gap', exact
    if any(r['status'] == 'different_position_or_phase' for r in predictions + annotations):
        return 'position_or_phase_mismatch', exact
    return 'annotation_or_alignment_unresolved', exact


def percent_absent(absent, ordinary):
    """Closest-absent loci as a percentage of absent plus ordinary, None when both are 0."""
    return 100 * absent / (absent + ordinary) if absent + ordinary else None


def main():
    # annotation and prediction reviews, alternative-intron states, locus results
    annotation_review = {row['key']: row for row in read_json_gz('model_review.json.gz')['rows']}
    predictions = defaultdict(list)
    for row in read_json_gz('prediction_review.json.gz')['rows']:
        predictions[row['locus_id']].append(row)
    alternatives = json.loads((OUT / 'alternative_composition.json').read_text())
    telomeric = {r['locus_id'] for r in alternatives if r['state_by_focal_motif'] == 'TELOTRON'}
    rows = read_tsv(OUT / 'locus_results.tsv')

    for row in rows:
        locus_predictions = predictions[row['locus_id']]
        annotations = [annotation_review[p['key']] for p in locus_predictions]
        status, exact = combined_status(row, locus_predictions, annotations, telomeric)
        row['combined_status'] = status
        row['prediction_statuses'] = ';'.join(p['status'] for p in locus_predictions)
        row['exact_prediction_pairs'] = len(exact)
        row['exact_alternative_telomeric'] = row['locus_id'] in telomeric
    write_tsv('combined_locus_results.tsv', rows)

    # per-lineage counts and conditional absence percentages
    lineages = []
    before = {row['clade']: row for row in read_tsv(OUT / 'lineage_changes.tsv')}
    for clade in sorted(before):
        counts = Counter(row['combined_status'] for row in rows if row['clade'] == clade)
        n_loci = sum(counts.values())
        absent = int(before[clade]['closest_absent_unchanged'])
        corroborated = counts['annotation_corroborated_exact']
        with_predictions = corroborated + counts['prediction_only_exact']
        lineages.append(
            dict(
                clade=clade,
                ordinary_before=n_loci,
                annotation_corroborated=corroborated,
                prediction_only=counts['prediction_only_exact'],
                position_or_phase_mismatch=counts['position_or_phase_mismatch'],
                exact_alternative_telomeric=counts['exact_alternative_is_telomeric'],
                remaining_uncertain=(
                    n_loci
                    - corroborated
                    - counts['prediction_only_exact']
                    - counts['position_or_phase_mismatch']
                    - counts['exact_alternative_is_telomeric']
                ),
                closest_absent_unchanged=absent,
                conditional_absent_original=float(before[clade]['conditional_absent_before']),
                conditional_absent_annotation_only=percent_absent(absent, corroborated),
                conditional_absent_including_prediction=percent_absent(absent, with_predictions),
            )
        )
    write_tsv('combined_lineage_results.tsv', lineages)

    # overall counts against the closest-state totals of lineage_changes.tsv
    counts = Counter(row['combined_status'] for row in rows)
    corroborated = counts['annotation_corroborated_exact']
    with_predictions = corroborated + counts['prediction_only_exact']
    primary_loci = sum(int(row['primary_loci']) for row in before.values())
    absent = sum(int(row['closest_absent_unchanged']) for row in before.values())
    both_annotated = sum(row['both_phases_annotated'] == 'True' for row in rows)

    def totals(ordinary):
        return dict(
            absent=absent,
            ordinary=ordinary,
            other=primary_loci - absent - ordinary,
            conditional_absent_percent=percent_absent(absent, ordinary),
        )

    inputs = [
        OUT / 'locus_results.tsv',
        OUT / 'model_review.json.gz',
        OUT / 'prediction_review.json.gz',
        OUT / 'alternative_composition.json',
    ]
    summary = dict(
        complete=True,
        loci=len(rows),
        comparison_pairs=sum(int(row['tied_comparators']) for row in rows),
        categories=dict(counts),
        annotated_phases_both=both_annotated,
        at_least_one_inferred_phase=corroborated - both_annotated,
        original=totals(len(rows)),
        annotation_corroborated_only=totals(corroborated),
        include_nonconflicting_prediction_only=totals(with_predictions),
        alternative_exact_telomeric_loci=sorted(telomeric),
        no_automatic_absence_reclassifications=True,
        production_release_unchanged=True,
        source_sha256={str(path): sha(path) for path in inputs},
        script_sha256=sha(__file__),
        limits=(
            'Frozen original closest-eligible choices; no reranking to more distant comparisons. '
            'Absence calls, other present-state cohorts and evolutionary polarization were not '
            'revalidated. Other includes unresolved and shared-telomeric loci. The conditional '
            'absence percentages are not gain rates.'
        ),
    )
    (OUT / 'final_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
