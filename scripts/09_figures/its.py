"""Figure 3: ITS density of bearers and controls, with the tested lineage pairs."""

from style import *
import math
import statistics

points = read('figure_points.tsv')
stats = json.loads((PACKAGE / 'tables/exact_wilcoxon.json').read_text())


def plot_y(density):
    # linear from 0 to 0.1 tracts/Mb (0 to 0.5 on the plot), log10 above
    return 0.5 * density / 0.1 if density <= 0.1 else 0.5 + math.log10(density / 0.1)


fig = plt.figure(figsize=(180 * MM, 110 * MM))
ax = fig.add_axes([0.10, 0.18, 0.49, 0.77])
side = fig.add_axes([0.65, 0.16, 0.34, 0.79])
side.axis('off')

# place points in two columns, nudging sideways until none overlap
coordinates = {}
for role, center in [('bearer', 0), ('control', 1)]:
    placed = []
    column = sorted(
        (point for point in points if point['role'] == role),
        key=lambda point: (float(point['density']), point['genome_id']),
    )
    for point in column:
        y = plot_y(float(point['density']))
        # try the centre first, then alternately left and right
        offsets = [0] + [sign * j * 0.055 for j in range(1, 15) for sign in [-1, 1]]
        offset = next(
            offset
            for offset in offsets
            if all(
                math.hypot((center + offset - x) * 190, (y - other_y) * 65) > 7
                for x, other_y in placed
            )
        )
        coordinates[point['genome_id']] = (center + offset, y)
        placed.append((center + offset, y))

# lines joining each tested bearer-control pair (cryptophyte dashed)
for pair in stats['primary_pairs']:
    bearer = coordinates[pair['bearer_id']]
    control = coordinates[pair['control_id']]
    ax.plot(
        [bearer[0], control[0]],
        [bearer[1], control[1]],
        c=COLORS[CLADES.index(pair['group'])],
        lw=0.8,
        alpha=0.55,
        ls='--' if pair['group'] == 'cryptophyte' else '-',
        zorder=1,
    )

# points (squares for MAGs, circles for isolates) and bars at the medians
for point in points:
    x, y = coordinates[point['genome_id']]
    ax.scatter(
        x,
        y,
        s=14,
        marker='s' if point['assembly_class'] == 'MAG' else 'o',
        c=COLORS[CLADES.index(point['clade'])],
        edgecolors='white',
        lw=0.3,
        zorder=3,
    )
for role, x in [('bearer', 0), ('control', 1)]:
    median = statistics.median(float(point['density']) for point in points if point['role'] == role)
    ax.hlines(plot_y(median), x - 0.24, x + 0.24, color=INK, lw=1.4, zorder=4)

ax.set(xlim=(-0.45, 1.45), ylim=(-0.08, plot_y(2000) + 0.1))
ticks = [0, 0.1, 1, 10, 100, 1000]
ax.set_yticks([plot_y(tick) for tick in ticks], [str(tick) for tick in ticks])
ax.set_xticks([0, 1], ['Bearers\n(n = 25)', 'Controls\n(n = 25)'])
ax.set_ylabel('ITS tracts per Mb')
axes_style(ax)
ax.tick_params(axis='x', length=0, pad=5)
ax.grid(axis='y', color='#e5e5e5', lw=0.5, zorder=0)

# side panel: lineage key, marker key, median key and test results
side.text(0, 0.96, 'Lineage', fontsize=8, fontweight='bold', transform=side.transAxes)
for i, (name, color) in enumerate(zip(NAMES, COLORS)):
    y = 0.88 - i * 0.075
    side.scatter(0.025, y, s=18, color=color, transform=side.transAxes)
    side.text(0.09, y, display_name(name), va='center', fontsize=7.5, transform=side.transAxes)
for y, marker, label in [(0.29, 'o', 'Isolate'), (0.22, 's', 'MAG')]:
    side.scatter(0.025, y, s=18, marker=marker, c=MUTED, transform=side.transAxes)
    side.text(0.09, y, label, va='center', fontsize=7.5, transform=side.transAxes)
side.plot([0, 0.055], [0.15, 0.15], c=INK, lw=1.4, transform=side.transAxes)
side.text(0.09, 0.15, 'Median', va='center', fontsize=7.5, transform=side.transAxes)
side.text(
    0,
    0.055,
    'Exact Wilcoxon\n7 pairs: p = 0.0156\n6 pairs: p = 0.0313',
    fontsize=7.3,
    linespacing=1.4,
    va='top',
    transform=side.transAxes,
)
side.text(
    0,
    -0.075,
    'Six-pair test omits the cryptophyte\ncomparison (dashed line, 40.8% BUSCO).',
    fontsize=6.4,
    linespacing=1.4,
    va='top',
    transform=side.transAxes,
)

save(fig, 'Figure3_ITS_density')
write_qa('ITS_plot_coordinates.json', coordinates)
