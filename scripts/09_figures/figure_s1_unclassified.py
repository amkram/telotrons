"""Supplementary Figure S1: classified loci and the first rule unclassified loci fail, by tier."""

import collections
import gzip

from style import *

catalogue_path = Path(__file__).resolve().parents[2] / 'work/catalogue/catalogue.tsv.gz'
with gzip.open(catalogue_path, 'rt') as handle:
    catalogue = list(csv.DictReader(handle, delimiter='\t'))

# display categories, top to bottom, starting with classified loci
CATEGORIES = [
    'Classified (R, F, F–R, F–L–R, R–L–F)',
    'Repeat sequence outside the arrays',
    'More than two arrays',
    'Two same-orientation arrays, separated',
    'Array ends >8 bp from an intron boundary',
    'Opposing arrays with a gap that is not a linker',
    'Single or convergent arrays at ≤85% coverage',
    'No array of three units, or ambiguous bases',
]


def first_reason(locus):
    return next((reason for reason in locus['architecture_reason'].split(';') if reason), '')


def category(locus):
    # classified loci, else the category of the first rule failed
    if locus['architecture'] in ARCH:
        return CATEGORIES[0]
    reason = first_reason(locus)
    if reason == 'repeat_islands_outside_clean_arrays':
        return CATEGORIES[1]
    if reason == 'more_than_two_arrays':
        return CATEGORIES[2]
    if reason == 'separated_arrays_of_same_orientation':
        return CATEGORIES[3]
    if reason.endswith('not_within_8bp_of_intron_or_tsd_boundary'):
        return CATEGORIES[4]
    if reason in (
        'interarray_gap_without_qualifying_linker',
        'degenerate_repeat_extension_leaves_lt18bp_core',
    ):
        return CATEGORIES[5]
    if reason.endswith('requires_strictly_gt85_percent_coverage'):
        return CATEGORIES[6]
    return CATEGORIES[7]


# loci per category and tier, with the catch-all category split by reason
counts = {
    tier: collections.Counter(category(locus) for locus in catalogue if locus['tier'] == tier)
    for tier in ['standard', 'additional']
}
catch_all = {
    tier: collections.Counter(
        first_reason(locus) or 'none'
        for locus in catalogue
        if locus['tier'] == tier and category(locus) == CATEGORIES[7]
    )
    for tier in counts
}

# source table and counts record
table = []
for tier, name in [('standard', 'primary'), ('additional', 'additional')]:
    total = sum(counts[tier].values())
    for label in CATEGORIES:
        table.append(
            dict(
                tier=name,
                category=label,
                loci=counts[tier][label],
                percent=100 * counts[tier][label] / total,
                tier_total=total,
            )
        )
with open(PACKAGE / 'tables/FigureS1_unclassified_reasons.tsv', 'w') as handle:
    writer = csv.DictWriter(handle, list(table[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(table)
record = dict(
    counts={tier: dict(tally) for tier, tally in counts.items()},
    other_detail={tier: dict(tally) for tier, tally in catch_all.items()},
)
(PACKAGE / 'results/figure_s1_reasons.json').write_text(json.dumps(record, indent=1) + '\n')

# one bar panel per tier, side by side
fig = plt.figure(figsize=(WIDTH, 78 * MM))
for j, (tier, title) in enumerate(
    [('standard', 'Primary (>85% coverage)'), ('additional', 'Additional (>50% to ≤85%)')]
):
    ax = fig.add_axes([0.36 if j == 0 else 0.70, 0.14, 0.24, 0.70])
    total = sum(counts[tier].values())
    values = [counts[tier][label] for label in CATEGORIES]
    # classified in blue, every failure category in grey
    colors = ['#2A78D6'] + ['#b8b8b8'] * (len(CATEGORIES) - 1)
    ax.barh(range(len(CATEGORIES)), values, color=colors, height=0.62, edgecolor='white', lw=0.3)

    # count and percent labels, with x room left for them
    x_max = max(values) * 1.42
    ax.set_xlim(0, x_max)
    for i, value in enumerate(values):
        ax.text(
            value + x_max * 0.015,
            i,
            f'{value:,} ({100 * value / total:.1f}%)' if value else '0',
            va='center',
            fontsize=6.4,
            color=INK,
        )

    # categories run top to bottom, named in panel A only
    ax.set_ylim(len(CATEGORIES) - 0.4, -0.6)
    ax.set_yticks(range(len(CATEGORIES)))
    ax.set_yticklabels(CATEGORIES if j == 0 else [], fontsize=6.8)
    ax.set_title(f'{title}\n{total:,} loci', fontsize=7.5, pad=6)
    axes_style(ax)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Loci (n)', fontsize=7)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    panel(fig, 'AB'[j], 0.04 if j == 0 else 0.655, 0.97)
save(fig, 'FigureS1_unclassified_reasons')
print(json.dumps(record['counts'], indent=1))
