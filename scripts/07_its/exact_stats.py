#!/usr/bin/env python3
"""Exact two-sided Wilcoxon signed-rank tests on bearer-control ITS density differences.

The primary test uses one prespecified pair per lineage path (seven pairs). Sensitivity analyses
drop the cryptophyte pair, use the nearest disjoint pairs, and try every representative choice in
the three lineages with several pairs. Reads collect.py's pairs.tsv and writes exact_wilcoxon.json.
"""

import argparse, csv, itertools, json, statistics
from collections import defaultdict
from pathlib import Path

# the density definition of collect.py
METRIC = dict(
    metric='projected',
    method='one_mm4',
    cap_mask='corrected_cap',
    policy='expanded_span',
    denominator='eligible_A_C_G_T_bp',
)
GROUPS = [
    'apicomplexa',
    'dinoflagellate',
    'haptophyte',
    'nematode',
    'cryptophyte',
    'chrysophyte',
    'stramenopile',
]
# expected number of pairs per group
PAIRS_PER_GROUP = [8, 5, 8, 1, 1, 1, 1]
# prespecified representative bearer per group
REPRESENTATIVES = {
    'apicomplexa': 'GCA_051911535.1',
    'dinoflagellate': 'GCA_905221605.1',
    'haptophyte': 'TARA_AOS_82_MAG_00154',
    'nematode': 'GCA_958295435.1',
    'cryptophyte': 'GCA_026770585.1',
    'chrysophyte': 'TARA_ARC_108_MAG_00234',
    'stramenopile': 'TARA_MED_95_MAG_00464',
}
# strict nearest-node disjoint pairing: six groups, no stramenopile
STRICT_BEARERS = {
    'apicomplexa': 'GCF_000499425.1',
    'dinoflagellate': 'GCA_905221605.1',
    'haptophyte': 'TARA_AOS_82_MAG_00154',
    'nematode': 'GCA_958295435.1',
    'cryptophyte': 'GCA_026770585.1',
    'chrysophyte': 'TARA_ARC_108_MAG_00234',
}


def exact_test(differences):
    """Exact two-sided signed-rank test by enumerating all 2**n sign assignments."""
    # drop zero differences (full enumeration needs n < 20)
    differences = [value for value in differences if value != 0]
    n = len(differences)
    assert n < 20

    # ranks of |difference|, with ties sharing their average rank
    order = sorted(range(n), key=lambda i: abs(differences[i]))
    rank = [0.0] * n
    first = 0
    while first < n:
        last = first + 1
        while last < n and abs(differences[order[last]]) == abs(differences[order[first]]):
            last += 1
        for position in range(first, last):
            rank[order[position]] = (first + 1 + last) / 2
        first = last

    # W is the smaller of the positive and negative rank sums
    total = sum(rank)
    positive = sum(rank[i] for i in range(n) if differences[i] > 0)
    statistic = min(positive, total - positive)

    # p = share of all 2**n sign flips at least as extreme as W
    extreme = 0
    for signs in itertools.product([0, 1], repeat=n):
        flipped = sum(rank[i] for i, sign in enumerate(signs) if sign)
        if min(flipped, total - flipped) <= statistic + 1e-10:
            extreme += 1

    return dict(
        n=n,
        W_plus=positive,
        W_minus=total - positive,
        W=statistic,
        exact_two_sided_p=extreme / 2**n,
    )


def difference(pair):
    """Bearer minus control density of one pair."""
    return float(pair['density_difference'])


def pair_of(group, bearer):
    """The pair of one bearer within its group."""
    return next(pair for pair in by_group[group] if pair['bearer_id'] == bearer)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pairs', required=True, help='pairs.tsv from collect.py')
parser.add_argument('--out', required=True, help='exact_wilcoxon.json')
args = parser.parse_args()

# bearer-control pairs grouped by lineage
with open(args.pairs) as handle:
    pairs = list(csv.DictReader(handle, delimiter='\t'))
by_group = defaultdict(list)
for pair in pairs:
    by_group[pair['group']].append(pair)
assert [len(by_group[group]) for group in GROUPS] == PAIRS_PER_GROUP

# primary test, then without the cryptophyte pair, then on the strict nearest disjoint pairs
chosen = [pair_of(group, REPRESENTATIVES[group]) for group in GROUPS]
primary = exact_test([difference(pair) for pair in chosen])
without_cryptophyte = exact_test(
    [difference(pair) for pair in chosen if pair['group'] != 'cryptophyte']
)
strict = exact_test(
    [difference(pair_of(group, bearer)) for group, bearer in STRICT_BEARERS.items()]
)

# all 8 x 5 x 8 representative choices (one-pair groups stay fixed)
p_values = []
medians = []
for combination in itertools.product(
    by_group['apicomplexa'], by_group['dinoflagellate'], by_group['haptophyte']
):
    selected = [*combination, *[by_group[group][0] for group in GROUPS[3:]]]
    p_values.append(exact_test([difference(pair) for pair in selected])['exact_two_sided_p'])
    medians.append(statistics.median(difference(pair) for pair in selected))
assert len(p_values) == 320

# within each group, every bearer is denser than every selected control
separation = {
    group: dict(
        min_bearer_density=min(float(pair['bearer_density']) for pair in by_group[group]),
        max_selected_control_density=max(
            float(pair['control_density']) for pair in by_group[group]
        ),
    )
    for group in GROUPS
}
for densities in separation.values():
    assert densities['min_bearer_density'] > densities['max_selected_control_density']

summary = dict(
    metric=METRIC,
    primary=primary,
    without_cryptophyte=without_cryptophyte,
    strict_individual_nearest_node_disjoint=strict,
    sensitivity=dict(
        independent_representative_combinations=len(p_values),
        min_p=min(p_values),
        max_p=max(p_values),
        median_difference_range=[min(medians), max(medians)],
        all_25_differences_positive=all(difference(pair) > 0 for pair in pairs),
        min_pair_difference=min(difference(pair) for pair in pairs),
        selected_control_reassignment_sign_invariant_by_clade=True,
        clade_separation=separation,
    ),
    primary_pairs=[
        dict(
            group=pair['group'],
            bearer_id=pair['bearer_id'],
            control_id=pair['control_id'],
            bearer_density=float(pair['bearer_density']),
            control_density=float(pair['control_density']),
            difference=difference(pair),
        )
        for pair in chosen
    ],
)
Path(args.out).write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
