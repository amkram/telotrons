"""Closest comparators before (recall) and after the absence correction, for both tiers and their
union (expanded), under results/absence_correction/: pair_decisions.tsv, closest_state_changes.tsv,
overall_before_after.tsv and lineage_before_after.tsv.
"""

from collections import Counter

from closest import TOP_N, closest_selection, states_by_locus
from common import STAGES, WORK, genomes, orthology_dir, read, write

DECISION_FIELDS = ['locus', 'ortholog', 'protein', 'old_state', 'new_state', 'reason']


def closest_summary(states, host_status):
    """Closest state, top identity and top-TOP_N counts of one locus."""
    status, eligible, best, top = closest_selection(states, host_status)
    return dict(
        state=status,
        identity=float(best[0]['protein_ident']) if best else '',
        eligible=len(eligible),
        top10_absent=sum(state['state'] == 'NO_INTRON' for state in top),
        top10_ordinary=sum(state['state'] == 'ORDINARY' for state in top),
        at_least10=len(eligible) >= TOP_N,
    )


def before_after(rows, prefix):
    """Absent, ordinary and resolved counts and percentages of one side of the comparison."""
    counts = Counter(row[prefix + 'state'] for row in rows)
    absent = counts['NO_INTRON']
    ordinary = counts['ORDINARY']
    resolved = absent + ordinary
    top10_absent_only = sum(
        row[prefix + 'top10_absent'] > 0 and row[prefix + 'top10_ordinary'] == 0 for row in rows
    )
    return {
        prefix + 'absent': absent,
        prefix + 'ordinary': ordinary,
        prefix + 'resolved': resolved,
        prefix + 'unresolved': len(rows) - resolved,
        prefix + 'absent_percent': 100 * absent / resolved if resolved else '',
        prefix + 'resolved_percent': 100 * resolved / len(rows) if rows else '',
        prefix + 'top10_absent_no_ordinary': top10_absent_only,
        prefix + 'at_least10': sum(row[prefix + 'at_least10'] for row in rows),
    }


def summarize(name, key, value, rows):
    """One before/after summary row of a set of loci."""
    summary = dict(set=name, **{key: value}, loci=len(rows))
    summary.update(before_after(rows, 'before_'))
    summary.update(before_after(rows, 'after_'))
    summary['changed_closest_loci'] = sum(row['before_state'] != row['after_state'] for row in rows)
    return summary


def main():
    loci = {}
    decisions = []
    for stage in STAGES:
        catalogue = read(WORK / 'inputs' / f'{stage}.tsv')
        loci[stage] = []

        # genome by genome (genomes.json order), loci in catalogue order
        for meta in genomes():
            genome_id = meta['genome_id']
            host_audit = read(WORK / stage / 'host' / genome_id / 'host_audit.tsv')
            host = {row['locus_id']: row['status'] for row in host_audit}
            before = states_by_locus(orthology_dir(stage, genome_id, 'recall') / 'states.tsv')
            corrected_dir = orthology_dir(stage, genome_id, 'corrected')
            after = states_by_locus(corrected_dir / 'states.tsv')
            for decision in read(corrected_dir / 'pair_decisions.tsv'):
                decisions.append(
                    dict(stage=stage, genome_id=genome_id, clade=meta['clade'], **decision)
                )
            for locus in catalogue:
                if locus['genome_id'] != genome_id:
                    continue
                tid = locus['locus_id']
                row = dict(
                    stage=stage,
                    genome_id=genome_id,
                    organism=meta['organism'],
                    clade=meta['clade'],
                    locus_id=tid,
                )
                row.update(
                    {'before_' + k: v for k, v in closest_summary(before[tid], host[tid]).items()}
                )
                row.update(
                    {'after_' + k: v for k, v in closest_summary(after[tid], host[tid]).items()}
                )
                loci[stage].append(row)
    loci['expanded'] = loci['primary'] + loci['additional']

    # overall and per-lineage summaries, and the loci whose closest state changed
    overall = []
    lineages = []
    for name, rows in loci.items():
        overall.append(summarize(name, 'scope', 'all', rows))
        for clade in sorted({row['clade'] for row in rows}):
            lineages.append(
                summarize(name, 'clade', clade, [r for r in rows if r['clade'] == clade])
            )
    changes = [row for row in loci['expanded'] if row['before_state'] != row['after_state']]

    out = WORK / 'results/absence_correction'
    write(out / 'overall_before_after.tsv', overall)
    write(out / 'lineage_before_after.tsv', lineages)
    write(out / 'closest_state_changes.tsv', changes, list(loci['expanded'][0]))
    write(out / 'pair_decisions.tsv', decisions, ['stage', 'genome_id', 'clade'] + DECISION_FIELDS)
    for row in overall:
        print(
            row['set'],
            'absent',
            row['before_absent'],
            '->',
            row['after_absent'],
            'ordinary',
            row['before_ordinary'],
            '->',
            row['after_ordinary'],
            'changed',
            row['changed_closest_loci'],
        )


if __name__ == '__main__':
    main()
