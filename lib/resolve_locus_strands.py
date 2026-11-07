"""One record per physical locus, without inventing transcript-strand evidence."""

from collections import defaultdict


def locus_key(row):
    """(genome, contig, start, end) of a physical locus."""
    return row['genome_id'], row['seqid'], str(row['start']), str(row['end'])


def resolve_strands(rows, curated_strands):
    """One record per locus: its only strand record, else its curated strand, else unresolved."""
    records = defaultdict(list)
    for row in rows:
        records[locus_key(row)].append(row)

    resolved = []
    for key in sorted(records):
        candidates = records[key]

        # one record per strand, all with the same repeat coverage
        assert len({row['strand'] for row in candidates}) == len(candidates)
        assert len({int(row['covered_bp']) for row in candidates}) == 1
        curated = [
            row for row in candidates if curated_strands.get(row['locus_id']) == row['strand']
        ]

        # the strand_resolution labels are those of the published catalogue
        unresolved = False
        if len(candidates) == 1:
            chosen = candidates[0]
            reason = 'single_strand_record'
        elif curated:
            chosen = curated[0]
            reason = 'preserve_previous_catalogue_strand'
        else:
            # a display strand only, strand='.' below blocks a biological call
            chosen = min(candidates, key=lambda row: row['strand'])
            reason = 'unresolved_conflicting_strands'
            unresolved = True

        locus = dict(
            chosen,
            strand_resolution=reason,
            input_strand_options=';'.join(sorted(row['strand'] for row in candidates)),
            opposite_strand_records_removed=len(candidates) - 1,
            display_strand=chosen['strand'],
            display_splice_class=chosen['splice_class'],
        )

        # unresolved loci lose every strand-dependent call
        if unresolved:
            locus['strand'] = '.'
            locus['splice_class'] = 'unresolved_strand'
            locus['architecture'] = 'mixed_unclear'
            locus['architecture_reason'] = ';'.join(
                reason_part for reason_part in [locus['architecture_reason'], reason] if reason_part
            )
            locus['clean_depends_on_tsd_adjustment'] = 0
        resolved.append(locus)
    return resolved
