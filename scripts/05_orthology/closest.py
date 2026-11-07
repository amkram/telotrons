"""Closest eligible comparator of every catalogued locus in both tiers (closest_comparators.tsv)."""

import json
from collections import Counter, defaultdict

from common import STAGES, WORK, genomes, orthology_dir, read, write

ELIGIBLE = ('NO_INTRON', 'ORDINARY')
TOP_N = 10  # comparators summarized by the top10_* columns, plus ties at the boundary


def states_by_locus(path):
    """State rows of one genome, grouped by telotron locus."""
    states = defaultdict(list)
    for state in read(path):
        states[state['telotron']].append(state)
    return states


def closest_selection(states, host_status):
    """(closest state, eligible states by identity, tied best, top TOP_N with boundary ties)."""
    eligible = sorted(
        (state for state in states if state['state'] in ELIGIBLE),
        key=lambda state: (-float(state['protein_ident']), state['ortholog']),
    )
    if not eligible:
        status = 'NO_ELIGIBLE_COMPARATOR' if host_status == 'mapped' else 'HOST_UNRESOLVED'
        return status, eligible, [], []

    # every comparator tied at the top identity counts, and conflicting tied states give TIED_STATES
    best_identity = float(eligible[0]['protein_ident'])
    best = [state for state in eligible if float(state['protein_ident']) == best_identity]
    best_states = {state['state'] for state in best}
    status = next(iter(best_states)) if len(best_states) == 1 else 'TIED_STATES'
    cutoff = float(eligible[min(TOP_N, len(eligible)) - 1]['protein_ident'])
    top = [state for state in eligible if float(state['protein_ident']) >= cutoff]
    return status, eligible, best, top


def closest_row(locus, clade, host, states):
    """Closest-comparator summary of one locus."""
    status, eligible, best, top = closest_selection(states, host['status'])
    return dict(
        locus_id=locus['locus_id'],
        genome_id=locus['genome_id'],
        clade=clade,
        tier=locus['tier'],
        architecture=locus['architecture'],
        host_status=host['status'],
        closest_state=status,
        closest_identity=best[0]['protein_ident'] if best else '',
        closest_tied_comparators=len(best),
        eligible_comparators=len(eligible),
        top10_including_boundary_ties=len(top),
        top10_NO_INTRON=sum(state['state'] == 'NO_INTRON' for state in top),
        top10_ORDINARY=sum(state['state'] == 'ORDINARY' for state in top),
        at_least10_eligible=int(len(eligible) >= TOP_N),
        any_ORDINARY=int(any(state['state'] == 'ORDINARY' for state in states)),
        any_NO_INTRON=int(any(state['state'] == 'NO_INTRON' for state in states)),
        any_TELOTRON=int(any(state['state'] == 'TELOTRON' for state in states)),
        protein_id=host['protein_id'],
        intron_residue=host['intron_residue'],
        intron_phase=host['intron_phase'],
        closest_comparators=';'.join(state['ortholog'] for state in best),
        all_states=len(states),
    )


def main():
    clades = {meta['genome_id']: meta['clade'] for meta in genomes()}
    closest = []

    # one row per catalogued locus in catalogue order, reading each genome once
    for stage in STAGES:
        loaded = {}
        for locus in read(WORK / 'inputs' / f'{stage}.tsv'):
            genome_id = locus['genome_id']
            if genome_id not in loaded:
                host_audit = read(WORK / stage / 'host' / genome_id / 'host_audit.tsv')
                reviewed = orthology_dir(stage, genome_id, 'reviewed') / 'states.tsv'
                host = {row['locus_id']: row for row in host_audit}
                loaded[genome_id] = host, states_by_locus(reviewed)
            host, states = loaded[genome_id]
            locus_id = locus['locus_id']
            closest.append(closest_row(locus, clades[genome_id], host[locus_id], states[locus_id]))
    write(WORK / 'results/closest_comparators.tsv', closest)

    summary = {
        tier: dict(Counter(row['closest_state'] for row in closest if row['tier'] == tier))
        for tier in sorted({row['tier'] for row in closest})
    }
    print(json.dumps(dict(loci=len(closest), by_tier=summary), indent=1))


if __name__ == '__main__':
    main()
