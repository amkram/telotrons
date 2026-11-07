"""Supplementary Figure S3: introns >85% one repeat motif, telomeric or any other, in bearers
and controls (stage 10 scans of config/single_motif_genomes.tsv)."""

import collections
import csv
import gzip
import random
import statistics

from style import *

REPO = Path(__file__).resolve().parents[2]
SCANS = REPO / 'work/10_single_motif'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def tsv(path):
    with gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def within_one_edit(a, b):
    """Equal, or one substitution, insertion or deletion apart."""
    if a == b:
        return True
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) <= 1
    if abs(len(a) - len(b)) != 1:
        return False
    short, long = sorted((a, b), key=len)
    return any(long[:i] + long[i + 1 :] == short for i in range(len(long)))


def rotations(motif):
    return {motif[i:] + motif[:i] for i in range(len(motif))}


lineage_motif = {row['genome_id']: row['lineage_motif'] for row in tsv(REPO / 'config/single_motif_genomes.tsv')}
rows = []
for genome in read('FigureS2_genomes.tsv'):
    genome_id = genome['genome_id']
    scan = SCANS / genome_id
    assessable = json.loads((scan / 'summary.json').read_text())['counts']['assessable_introns']
    group = {row['locus_id']: row['group'] for row in tsv(scan / 'groups.tsv')}
    motif = lineage_motif[genome_id]
    lineage_rotations = rotations(motif) | rotations(motif.translate(COMPLEMENT)[::-1])

    def telomere_like(candidate):
        # any rotation within one edit of a lineage motif rotation
        return any(within_one_edit(x, y) for x in rotations(candidate) for y in lineage_rotations)

    # best motif per duplicate group: highest coverage, then lineage, telomeric, shorter motif
    best = {}
    for hit in tsv(scan / 'hits.tsv.gz'):
        length = int(hit['intron_bp'])
        # the rule of deduplicate.py: the screen measure for periods 6-8 and telomere motifs,
        # otherwise the tract measure
        screened = int(hit['period']) >= 6 or int(hit['is_telomeric'])
        if screened and int(hit['screen_pass']) and int(hit['screen_bp']) * 100 > 85 * length:
            covered, passed = int(hit['screen_bp']), 1
        else:
            covered, passed = int(hit['tract_bp']), int(hit['tract_pass'])
        if not passed or covered * 100 <= 85 * length:
            continue
        unit = group.get(hit['locus_id'], hit['locus_id'])
        key = (
            covered / length,
            int(hit['is_lineage_motif']),
            int(hit['is_telomeric']),
            -len(hit['motif']),
            hit['motif'],
        )
        if unit not in best or key > best[unit][0]:
            best[unit] = (key, hit)

    # other motifs split by period: short (1-4 bp) and long (5-8 bp)
    classes = collections.Counter()
    others = {'short': collections.Counter(), 'long': collections.Counter()}
    for _, hit in best.values():
        if int(hit['is_lineage_motif']):
            classes['lineage'] += 1
        elif int(hit['is_telomeric']) or telomere_like(hit['motif']):
            classes['telomere_like'] += 1
        else:
            size = 'short' if len(hit['motif']) < 5 else 'long'
            classes[size] += 1
            others[size][hit['motif']] += 1
    top = {size: (counts.most_common(1) or [('', 0)])[0] for size, counts in others.items()}
    rows.append(
        dict(
            genome_id=genome_id,
            label=genome['label'],
            clade=genome['clade'],
            role=genome['role'],
            lineage_motif=motif,
            assessable_introns=assessable,
            single_motif_gt85=sum(classes.values()),
            lineage_motif_introns=classes['lineage'],
            telomere_like_introns=classes['telomere_like'],
            other_short_introns=classes['short'],
            top_other_short_motif=top['short'][0],
            other_long_introns=classes['long'],
            top_other_long_motif=top['long'][0],
            catalogue_primary=int(genome['gt85']),
        )
    )

with open(PACKAGE / 'tables/FigureS3_single_motif_introns.tsv', 'w') as handle:
    writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)

# six columns: telomeric, short other and long other motifs, in bearers then controls; each
# value carries the motif shown on its label (the lineage motif, or the most frequent other motif)
columns = []
for role, title in (('bearer', 'Bearers'), ('control', 'Controls')):
    members = [row for row in rows if row['role'] == role]
    telomeric = [row['lineage_motif_introns'] + row['telomere_like_introns'] for row in members]
    columns.append((f'{title}\ntelomeric motif', members, telomeric, [r['lineage_motif'] for r in members]))
    for size, span in (('short', '<5 bp'), ('long', '5–8 bp')):
        columns.append(
            (
                f'{title}\nany other motif\n{span}',
                members,
                [row[f'other_{size}_introns'] for row in members],
                [row[f'top_other_{size}_motif'] for row in members],
            )
        )
summary = {
    title.replace('\n', ' '): dict(genomes=len(values), median=statistics.median(values), total=sum(values))
    for title, _, values, _ in columns
}
(PACKAGE / 'results').mkdir(exist_ok=True)
(PACKAGE / 'results/figure_s3_summary.json').write_text(json.dumps(summary, indent=1) + '\n')

LABELED = 5  # genomes labeled per column, highest first
fig = plt.figure(figsize=(170 * MM, 90 * MM))
ax = fig.add_axes([0.09, 0.26, 0.89, 0.71])
jitter = random.Random(7)
points = []  # (column, x, value, motif) of every plotted genome
for x, (title, members, values, motifs) in enumerate(columns):
    ax.boxplot(
        [values],
        positions=[x],
        widths=0.46,
        showfliers=False,
        whis=(0, 100),
        patch_artist=True,
        boxprops=dict(facecolor='#f2f2f2', edgecolor='#999999', lw=0.6),
        medianprops=dict(color=INK, lw=1.0),
        whiskerprops=dict(color='#999999', lw=0.6),
        capprops=dict(color='#999999', lw=0.6),
        zorder=1,
    )
    # one point per genome: squares for MAGs, circles for isolates
    for row, value, motif in zip(members, values, motifs):
        px = x + (jitter.random() - 0.5) * 0.38
        points.append((x, px, value, motif))
        ax.scatter(
            px,
            value,
            s=16,
            marker='s' if row['genome_id'].startswith('TARA_') else 'o',
            facecolors=COLORS[CLADES.index(row['clade'])],
            edgecolors='white',
            lw=0.4,
            zorder=3,
        )
    median = statistics.median(values)
    ax.text(x - 0.27, median, f'{median:g}', ha='right', va='center', fontsize=6.8, color=INK)
ax.set_yscale('symlog', linthresh=1, linscale=0.5)
ax.set_ylim(-0.4, 6000)
ax.set_yticks([0, 1, 10, 100, 1000])
ax.set_yticklabels(['0', '1', '10', '100', '1,000'])
ax.set_xticks(range(len(columns)))
ax.set_xticklabels([title for title, _, _, _ in columns], fontsize=7)
ax.set_xlim(-0.6, 5.9)
ax.axvline(2.5, color='#cccccc', lw=0.6, ls=':')
ax.set_ylabel('Introns >85% one repeat motif per genome (n)', fontsize=7.5)
axes_style(ax)

# motif labels right of each column for its top genomes, spread at least GAP points apart;
# medians are labeled on the left
GAP = 7.5
to_display = ax.transData.transform
to_data = ax.transData.inverted().transform
for column in range(len(columns)):
    top = sorted((p for p in points if p[0] == column and p[2] > 0), key=lambda p: -p[2])[:LABELED]
    side = 1
    label_x = column + side * 0.30
    previous = None
    for _, px, value, motif in top:
        y = to_display((label_x, value))[1]
        if previous is not None and previous - y < GAP:
            y = previous - GAP
        previous = y
        label_y = to_data((0, y))[1]
        ax.plot([label_x - side * 0.02, px], [label_y, value], color='#bbbbbb', lw=0.4, zorder=2)
        ax.text(
            label_x,
            label_y,
            motif,
            ha='left' if side == 1 else 'right',
            va='center',
            fontsize=5.6,
            fontfamily='DejaVu Sans Mono',
            color=INK,
            bbox=dict(facecolor='white', edgecolor='none', pad=0.3),
            zorder=4,
        )
fig.legend(
    [Line2D([], [], marker='o', ls='', c=COLORS[i], markersize=4) for i in range(7)],
    [display_name(name) for name in NAMES],
    loc='lower center',
    bbox_to_anchor=(0.55, 0.0),
    ncol=4,
    frameon=False,
    fontsize=6,
    handletextpad=0.3,
    columnspacing=1.0,
    labelspacing=0.3,
)
save(fig, 'FigureS3_single_motif_introns')
