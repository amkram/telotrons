"""Figure 1 (lineage tree) and Supplementary Figure S2 (genome tree) with their tables.

Deep nodes are TimeTree ages (inputs/backbone_dated.nwk). Each lineage is a polytomy at its
published crown age (inputs/crown_ages.tsv). Rows keep the order of the *_topology.nwk files.
"""

from style import *
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from statistics import mean
from decimal import Decimal, ROUND_HALF_UP


def one_decimal(value):
    # halves rounded up rather than to even
    return str(Decimal(str(value)).quantize(Decimal('.1'), rounding=ROUND_HALF_UP))


CROWN = {row['lineage']: float(row['crown_Ma']) for row in read_inputs('crown_ages.tsv')}


def lineage(clade):
    # the chrysophyte sits inside the stramenopile polytomy
    return 'stramenopile' if clade == 'chrysophyte' else clade


def dated_tree(order_file, clade_of):
    """Backbone with each lineage tip replaced by its rows at the crown age, and node ages in Ma."""
    order = [tip.name for tip in Phylo.read(PACKAGE / 'inputs' / order_file, 'newick').get_terminals()]
    tree = Phylo.read(PACKAGE / 'inputs/backbone_dated.nwk', 'newick')
    for tip in tree.get_terminals():
        crown = CROWN[tip.name]
        members = [name for name in order if lineage(clade_of[name]) == tip.name]
        tip.branch_length -= crown
        tip.clades = [Clade(branch_length=crown, name=name) for name in members]
        tip.name = None
    depth = tree.depths()
    root_age = max(depth.values())
    return tree, {node: root_age - d for node, d in depth.items()}


def draw_tree(ax, tree, heights, ys, left, width):
    """Draw the tree with x by node age and inner nodes at their children's mean y."""

    def x_of(node):
        return left + width * (1 - heights[node] / heights[tree.root])

    def y_of(node):
        if node not in ys:
            ys[node] = mean(y_of(child) for child in node.clades)
        return ys[node]

    y_of(tree.root)

    def draw(node):
        if node.clades:
            ax.plot(
                [x_of(node)] * 2,
                [min(ys[child] for child in node.clades), max(ys[child] for child in node.clades)],
                c='#999999',
                lw=0.55,
            )
            for child in node.clades:
                ax.plot([x_of(node), x_of(child)], [ys[child]] * 2, c='#999999', lw=0.55)
                draw(child)

    draw(tree.root)
    return x_of


def time_axis(ax, tree, ages, x_of, y, size, label_every=500):
    """Age axis under the tree, with ticks every 500 Ma and labels every label_every Ma."""
    left, right = x_of(tree.root), x_of(tree.get_terminals()[0])
    ax.plot([left, right], [y] * 2, c=RULE, lw=0.5)
    for ma in range(0, int(ages[tree.root]) + 1, 500):
        x = right - (right - left) * ma / ages[tree.root]
        ax.plot([x] * 2, [y, y - 0.9], c=RULE, lw=0.5)
        if ma % label_every == 0:
            ax.text(x, y - 1.6, f'{ma:,}', ha='center', va='top', fontsize=size, color=MUTED)
    ax.text((left + right) / 2, y - 5.2, 'Age (Ma)', ha='center', va='top', fontsize=size, color=MUTED)


def legend(fig, y):
    # bearer/control key on the left, lineage colours on the right
    fig.legend(
        [
            Line2D([], [], marker='o', ls='', c=INK, markersize=3.5),
            Line2D([], [], marker='s', ls='', c='#888', markerfacecolor='white', markersize=3.5),
        ],
        [
            'Bearer: ≥5 introns with >85% repeats',
            'Control: no intron with >50% repeats under the lineage motif',
        ],
        loc='lower left',
        bbox_to_anchor=(0.025, y),
        frameon=False,
        fontsize=6.6,
        handletextpad=0.4,
        labelspacing=0.3,
    )
    fig.legend(
        [Line2D([], [], marker='o', ls='', c=COLORS[i], markersize=3.5) for i in range(7)],
        [display_name(name) for name in NAMES],
        loc='lower right',
        bbox_to_anchor=(0.985, y),
        ncol=3,
        frameon=False,
        fontsize=6.6,
        columnspacing=1.0,
        handletextpad=0.3,
        labelspacing=0.3,
    )


def tip_marker(ax, x, y, row):
    # filled circles for bearers, open squares for controls
    bearer = row['role'] == 'bearer'
    color = COLORS[CLADES.index(row['clade'])]
    ax.scatter(
        x,
        y,
        s=12,
        marker='o' if bearer else 's',
        facecolors=color if bearer else 'white',
        edgecolors=color,
        lw=0.6,
        zorder=3,
    )


# Figure 1: lineage tree, one row per bearer group or control set
CONTROL_LABELS = {
    'nematode_control': 'Pristionchus',
    'cryptophyte_control': 'Cryptophyta MAG',
    'haptophyte_control': 'Other haptophytes (n = 8)',
    'GCF_000186865.1': 'Aureococcus',
    'GCA_001750085.1': 'Fragilariopsis',
    'dinoflagellate_control': 'Other dinoflagellates (n = 5)',
    'apicomplexa_control': 'Other coccidians (n = 8)',
}
rows = {row['row_id']: row for row in read('Figure1_source_data.tsv')}
tree, ages = dated_tree('main_topology.nwk', {k: row['clade'] for k, row in rows.items()})
tips = tree.get_terminals()
HEIGHT = 106
fig = plt.figure(figsize=(180 * MM, HEIGHT * MM))
# full-figure axes in millimetres
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, 180), ylim=(0, HEIGHT))
ax.axis('off')
TOP, STEP = 84, 4.7
ys = {tip: TOP - i * STEP for i, tip in enumerate(tips)}
x_of = draw_tree(ax, tree, ages, ys, 5, 58)
time_axis(ax, tree, ages, x_of, ys[tips[-1]] - 6, 6.2)
# labels go through text() so the Figure S2 overlap check sees them
labels = []


def text(x, y, string, **kw):
    labels.append(ax.text(x, y, string, va='center', **kw))


LABEL, MOTIF, GT85, GT50, ABSENT, TSD = 65.5, 109, 127, 142, 155, 178
for tip in tips:
    row = rows[tip.name]
    y = ys[tip]
    bearer = row['role'] == 'bearer'
    ink = INK if bearer else MUTED
    tip_marker(ax, x_of(tip), y, row)
    name = CONTROL_LABELS.get(row['row_id'], row['label'].replace(' (n = 1)', ''))
    text(LABEL, y, display_name(name), fontsize=7, color=ink)
    text(MOTIF, y, row['dominant_motif'], ha='center', fontsize=7, color=ink)
    text(GT85, y, f"{int(row['gt85']):,}", ha='right', fontsize=7, color=ink)
    text(GT50, y, f"{int(row['gt50']):,}", ha='right', fontsize=7, color=ink)
    denominator = int(row['ortholog_denominator'])
    absent = int(row['closest_intron_absent'])
    absent_percent = one_decimal(100 * absent / denominator) if denominator else '—'
    text(ABSENT, y, absent_percent, ha='right', fontsize=7, color=ink)
    if denominator:
        text(ABSENT + 1, y, f'({absent:,}/{denominator:,})', fontsize=6.1, color=MUTED)
    tsd_percent = one_decimal(row['TSD_percent_gt85']) if bearer else '—'
    text(TSD, y, tsd_percent, ha='right', fontsize=7, color=ink)

# column headers and the two group rules
HEADER = TOP + 7.5
for x, header in [
    (MOTIF, 'Repeat\nmotif'),
    (GT85 - 4, '>85%'),
    (GT50 - 4, '>50%'),
    (ABSENT + 5, 'Intron-absent\n% (n/N)'),
    (TSD - 3, 'TSDs\n(%)'),
]:
    text(x, HEADER, header, ha='center', fontsize=7)
text((GT85 + GT50) / 2 - 4, HEADER + 7.5, 'Telotrons (n)', ha='center', fontsize=7.3)
text((ABSENT - 8 + TSD) / 2, HEADER + 7.5, 'Primary (>85%)', ha='center', fontsize=7.3)
ax.plot([GT85 - 12, GT50 + 1], [HEADER + 4.5] * 2, c=RULE, lw=0.5)
ax.plot([ABSENT - 8, TSD + 1], [HEADER + 4.5] * 2, c=RULE, lw=0.5)
legend(fig, 0.012)
save(fig, 'Figure1_telotron_tree')

# Figure S2: genome tree, one row per genome
rows = {row['genome_id']: row for row in read('FigureS2_genomes.tsv')}
tree, ages = dated_tree('genomes_topology.nwk', {k: row['clade'] for k, row in rows.items()})
tips = tree.get_terminals()
fig = plt.figure(figsize=(307 * MM, 307 * MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, 307), ylim=(0, 307))
ax.axis('off')
labels = []
ys = {tip: 276 - i * 4.8 for i, tip in enumerate(tips)}
x_of = draw_tree(ax, tree, ages, ys, 7, 24)
time_axis(ax, tree, ages, x_of, ys[tips[-1]] - 5, 6.6, label_every=1000)

# column headers
text(7, 294, 'Genome', fontsize=8)
text(146, 298, 'Telotrons (n)', ha='center', fontsize=8)
text(237, 298, 'Primary (>85%)', ha='center', fontsize=8)
for x, header in [
    (119, 'Repeat\nmotif'),
    (136, '>85%'),
    (155, '>50%'),
    (181, 'Intron length (bp)\nmedian [range]'),
    (210, 'Included in\ncomparison (%)'),
    (241, 'Intron-absent\n% (n/N)'),
    (268, 'TSDs\n(%)'),
    (291, 'TSD length (bp)\nmedian [range]'),
]:
    text(x, 287, header, ha='center', fontsize=7.2)
ax.plot([7, 301], [281] * 2, c=RULE, lw=0.5)


def median_range(row, prefix):
    return (
        f"{float(row[prefix + '_median_bp']):g} "
        f"[{int(row[prefix + '_min_bp']):,}–{int(row[prefix + '_max_bp']):,}]"
    )


previous_clade = None
for tip in tips:
    row = rows[tip.name]
    y = ys[tip]
    bearer = row['role'] == 'bearer'
    ink = INK if bearer else MUTED
    # light rule between lineages
    if previous_clade and previous_clade != row['clade']:
        ax.plot([33, 301], [y + 2.4] * 2, c='#dddddd', lw=0.4)
    previous_clade = row['clade']
    tip_marker(ax, x_of(tip), y, row)

    # genome label: add the id tail to Tara MAG names, italicise binomials
    label = row['label'].replace('Sister_Phaeocystis', 'Phaeocystis-related').split(' (')[0]
    if row['genome_id'].startswith('TARA_') and row['genome_id'][-5:] not in label:
        label += ' ' + row['genome_id'][-5:]
    if 'MAG' not in label and len(label.split()) >= 2 and not label.startswith('Stramenopile'):
        words = label.split()
        label = r'$\it{' + r'\ '.join(words[:2]) + r'}$' + ' ' + ' '.join(words[2:])
    else:
        label = display_name(label)
    text(33, y, label, fontsize=7.2, color=ink)
    text(119, y, row['motif'], ha='center', fontsize=7.3, color=ink)
    for x, key in [(141, 'gt85'), (161, 'gt50')]:
        text(x, y, f"{int(row[key]):,}", ha='right', fontsize=7.3, color=ink)

    # bearer-only columns show a dash for controls
    intron_lengths = median_range(row, 'intron') if bearer else '—'
    text(181, y, intron_lengths, ha='center', fontsize=6.8, color=ink)
    included = one_decimal(row['ortholog_included_percent']) if bearer else '—'
    text(216, y, included, ha='right', fontsize=7.3, color=ink)
    denominator = int(row['ortholog_denominator'])
    absent = int(row['closest_intron_absent'])
    absent_percent = one_decimal(100 * absent / denominator) if denominator else '—'
    text(235, y, absent_percent, ha='right', fontsize=7.3, color=ink)
    if denominator:
        text(237, y, f'({absent:,}/{denominator:,})', fontsize=6.6, color=MUTED)
    tsd_percent = one_decimal(row['TSD_percent']) if bearer else '—'
    text(272, y, tsd_percent, ha='right', fontsize=7.3, color=ink)
    tsd_lengths = median_range(row, 'TSD') if int(row['TSD_n']) else '—'
    text(291, y, tsd_lengths, ha='center', fontsize=7, color=ink)
legend(fig, 0.015)

# record any pair of overlapping labels for QA
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
bounds = [(label.get_text(), label.get_window_extent(renderer)) for label in labels]
overlaps = [
    (first, second)
    for i, (first, box) in enumerate(bounds)
    for second, other in bounds[i + 1 :]
    if box.overlaps(other)
]
write_qa('tree_text_overlaps.json', overlaps)
save(fig, 'FigureS2_genome_tree')
