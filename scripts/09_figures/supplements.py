"""Unnumbered supplementary figure: architecture shares of classified loci by coverage tier."""

from style import *

clades = read('clade_summary.tsv')
tiers = read('tier_summary.tsv')
fig = plt.figure(figsize=(WIDTH, 105 * MM))
bar_audit = []

# one panel per tier: architecture shares of each lineage's classified loci
for j, (tier, title) in enumerate(
    [('gt85', 'Primary (>85%)'), ('gt50_to85', 'Additional (>50% to ≤85%)')]
):
    left_edge = 0.205 if j == 0 else 0.65
    ax = fig.add_axes([left_edge, 0.26, 0.25, 0.54])
    total = next(row for row in tiers if row['tier'] == tier)
    for i, clade in enumerate(CLADES):
        row = next(row for row in clades if row['clade'] == clade and row['tier'] == tier)
        classified = int(row['clean'])

        # stacked segments, one per architecture
        left = 0
        assert sum(int(row[architecture]) for architecture in ARCH) == classified
        for architecture in ARCH:
            count = int(row[architecture])
            share = 100 * count / classified if classified else 0
            ax.barh(
                i,
                share,
                left=left,
                height=0.55,
                color=ARCH_COLORS[architecture],
                edgecolor='white',
                lw=0.35,
            )
            left += share
            bar_audit.append(
                dict(
                    tier=tier,
                    clade=clade,
                    architecture=architecture,
                    count=count,
                    clean=classified,
                    percent=share,
                )
            )

        # classified/all label beyond the bar
        ax.text(
            103,
            i,
            f"{classified:,}/{int(row['loci']):,}",
            ha='left',
            va='center',
            fontsize=7,
            clip_on=False,
        )

    ax.set_xlim(0, 100)
    ax.set_ylim(6.65, -0.65)
    ax.set_xticks([0, 50, 100])
    ax.set_yticks(range(7), [display_name(name) for name in NAMES] if j == 0 else [])
    axes_style(ax)
    ax.tick_params(axis='y', length=0, pad=5)
    ax.spines['left'].set_visible(False)
    ax.set_title(title, pad=28, fontweight='normal', fontsize=7.5)
    # tier total under the title
    ax.text(
        0.5,
        1.08,
        f"{int(total['clean']):,}/{int(total['loci']):,} classified "
        f"({100 * int(total['clean']) / int(total['loci']):.1f}%)",
        transform=ax.transAxes,
        ha='center',
        fontsize=7,
        color=MUTED,
    )
    ax.text(103, -0.95, 'Classified/all', fontsize=7, ha='left')
    panel(fig, 'A' if j == 0 else 'B', left_edge - 0.035, 0.945)

# shared x label and architecture legend
fig.text(0.55, 0.155, 'Share of classified telotrons (%)', ha='center', fontsize=7.5)
fig.legend(
    [Patch(fc=ARCH_COLORS[architecture]) for architecture in ARCH],
    [architecture.replace('-', '–') for architecture in ARCH],
    loc='lower center',
    bbox_to_anchor=(0.54, 0.04),
    ncol=6,
    frameon=False,
    fontsize=7.5,
    columnspacing=1.6,
    handlelength=1.2,
)
save(fig, 'FigureS1_architecture_by_tier')
write_qa('tier_bars.json', bar_audit)
