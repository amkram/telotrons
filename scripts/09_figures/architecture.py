"""Figure 2: architecture shares by lineage (A) and example loci with their closest orthologs (B)."""

from style import *
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.font_manager import FontProperties
from Bio.Seq import Seq
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from repeat_policy import repeat_labels

# inputs: architecture counts, curated example loci, their display audit and exon phases
bars = read('architecture_counts.tsv')
counts = {(row['clade'], row['architecture']): int(row['count']) for row in bars}
classified = {row['clade']: int(row['clean_denominator']) for row in bars}
# the caption cites 4,331 classified loci
assert sum(classified.values()) == 4331
examples = json.loads((PACKAGE / 'inputs/examples_matched_verified.json').read_text())
display = json.loads((PACKAGE / 'inputs/sequence_display_audit_matched.json').read_text())
phases = json.loads((PACKAGE / 'inputs/exon_phases_matched.json').read_text())['loci']

# 270 mm wide, in points, so every displayed base stays legible without new elisions
CANVAS_WIDTH = 270 * MM * 72
CANVAS_HEIGHT = 439  # points: panels A and B plus a two-row legend
fig = plt.figure(figsize=(CANVAS_WIDTH / 72, CANVAS_HEIGHT / 72))
# full-figure axes whose data units are points
canvas = fig.add_axes([0, 0, 1, 1])
canvas.set_xlim(0, CANVAS_WIDTH)
canvas.set_ylim(0, CANVAS_HEIGHT)
canvas.axis('off')

# sequence grid in points, with right exons ending 8 pt from the edge
PITCH = 3.92
BASE_FONT = 6.5
RIGHT_EXON_X = CANVAS_WIDTH - 8 - 12 * PITCH
# panel A bar origin and 100% length, in points
BAR_X = 94
BAR_WIDTH = 63


def text(x, y, string, size=7.5, **kw):
    kw.setdefault('va', 'center')
    return canvas.text(x, y, string, fontsize=size, **kw)


def panel_label(x, y, letter):
    text(x, y, letter, 10, fontweight='bold')


def genome_label(right, y, name):
    # native italics: mathtext shifts the baseline
    styled = display_name(name)
    if styled.startswith('$\\it{'):
        taxon, suffix = styled[5:].split('}$', 1)
        parts = [(taxon.replace('\\ ', ' '), 'italic'), (suffix, 'normal')]
    else:
        parts = [(styled, 'normal')]

    # right-align at `right` by measuring each styled part
    renderer = fig.canvas.get_renderer()
    widths = [
        renderer.get_text_width_height_descent(
            part, FontProperties(family='DejaVu Sans', style=style, size=7), False
        )[0]
        * 72
        / fig.dpi
        for part, style in parts
    ]
    x = right - sum(widths)
    for (part, style), width in zip(parts, widths):
        if part:
            text(x, y, part, 7, fontstyle=style, ha='left', va='baseline')
        x += width


def glyph(x, y, architecture):
    # one arrow per repeat arm (F rightward, R leftward), joined by L for linkers
    arms = [arm for arm in architecture.split('-') if arm != 'L']
    linker = 'L' in architecture
    if len(arms) == 1:
        spans = [(x, x + 41)]
    else:
        spans = [(x, x + 16 if linker else x + 20), (x + 25 if linker else x + 21, x + 41)]
    if linker:
        canvas.plot([x + 16, x + 25], [y, y], color=INK, lw=0.6)
        text(x + 20.5, y + 8, 'L', 6.5, ha='center')
    for (low, high), arm in zip(spans, arms):
        color = F_COLOR if arm == 'F' else R_COLOR
        start, dx = (low, high - low) if arm == 'F' else (high, low - high)
        canvas.add_patch(
            FancyArrow(
                start,
                y,
                dx,
                0,
                width=7,
                head_width=11,
                head_length=4,
                length_includes_head=True,
                facecolor=color + '25',
                edgecolor=color,
                lw=0.7,
            )
        )
        text((low + high) / 2, y, arm, 6.5, ha='center')


def tokens(architecture, name, intron):
    """Shown bases ('b', index) and omitted runs ('g', (start, end)) of one displayed intron."""
    if not intron:
        return []
    match = [
        row
        for row in display
        if row['architecture'] == architecture
        and row['name'].replace('MAG_', 'MAG ') == display_name(name)
    ]
    assert len(match) == 1, (architecture, name)
    shown = set(match[0]['shown_indices0'])
    omitted = {start: end for start, end in match[0]['elisions0']}
    out = []
    i = 0
    while i < len(intron):
        if i in shown:
            out.append(('b', i))
            i += 1
        else:
            assert i in omitted
            out.append(('g', (i, omitted[i])))
            i = omitted[i]
    # every base is either shown or inside an omitted run
    assert sum(1 if kind == 'b' else value[1] - value[0] for kind, value in out) == len(intron)
    return out


def cell_width(sequence_tokens):
    # display width in cells: one per base, four per omission
    return sum(1 if kind == 'b' else 4.0 for kind, _ in sequence_tokens)


panel_label(10, CANVAS_HEIGHT - 12, 'A')
text(125.5, CANVAS_HEIGHT - 12, 'Within lineage (%)', 7, ha='center')
panel_label(180, CANVAS_HEIGHT - 12, 'B')
layout_audit = []
codon_audit = []

# panels A and B: one horizontal band per architecture
band_top = CANVAS_HEIGHT - 29
band_heights = {'R': 76, 'F': 62, 'R-F': 38, 'F-R': 62, 'R-L-F': 62, 'F-L-R': 62}
for architecture in ARCH:
    band_height = band_heights[architecture]
    center = band_top - band_height / 2
    text(10, center, architecture.replace('-', '–'), 7.5, fontweight='bold')
    glyph(38, center, architecture)

    # panel A: each lineage's share of its classified loci, labelled with counts
    bars_top = center + 23.4
    for j, (clade, color) in enumerate(zip(CLADES, COLORS)):
        count = counts[clade, architecture]
        percent = 100 * count / classified[clade]
        length = BAR_WIDTH * percent / 100
        y = bars_top - j * 7.8
        if count:
            canvas.add_patch(Rectangle((BAR_X, y - 2.35), length, 4.7, fc=color, ec='none'))
            text(BAR_X + length + 2, y, f'{count:,}', 7)
    # bar axis only where loci exist (R-F has none and its band is too short for the rule)
    if any(counts[clade, architecture] for clade in CLADES):
        canvas.plot([BAR_X, BAR_X], [bars_top - 49.5, bars_top + 3], color=RULE, lw=0.5)
    # percentage scale, above the first band only
    if architecture == 'R':
        canvas.plot([BAR_X, BAR_X + BAR_WIDTH], [bars_top + 8] * 2, color=RULE, lw=0.5)
        for percent in [0, 50, 100]:
            tick_x = BAR_X + BAR_WIDTH * percent / 100
            canvas.plot([tick_x, tick_x], [bars_top + 6, bars_top + 8], color=RULE, lw=0.5)
            text(tick_x, bars_top + 14, str(percent), 6.5, ha='center')

    # panel B: the example locus and its comparisons, or a note
    example = next((item for item in examples if item['architecture'] == architecture), None)
    if example is None:
        text((180 + CANVAS_WIDTH) / 2, center, 'Not observed', 7.5, color=MUTED, ha='center')
        band_top -= band_height
        canvas.plot([10, CANVAS_WIDTH - 8], [band_top] * 2, color='#e5e5e5', lw=0.4)
        continue
    focal_locus = example['catalogue']
    context = example['context']
    sequence_rows = [dict(name=example['name'], **context, kind='focal')] + example['comparisons']
    group_span = max(
        cell_width(tokens(architecture, row['name'], row['intron'])) for row in sequence_rows
    )
    assert group_span == 89, (architecture, group_span)
    row_ys = [
        center + (len(sequence_rows) - 1) * 10 - 3 - 20 * j for j in range(len(sequence_rows))
    ]
    for row, y in zip(sequence_rows, row_ys):
        focal = row['kind'] == 'focal'
        intron = row['intron']
        intron_tokens = tokens(architecture, row['name'], intron)

        # widen the omissions so every row spans group_span cells
        omissions = sum(kind == 'g' for kind, _ in intron_tokens)
        padding = group_span - cell_width(intron_tokens) if intron else 0
        assert not padding or omissions, (architecture, row['name'], 'cannot pad a full sequence')
        omission_width = 4.0 + padding / omissions if omissions else 4.0
        left_exon_x = RIGHT_EXON_X - (14 + group_span) * PITCH
        intron_x = left_exon_x + 13 * PITCH
        assert left_exon_x >= 292, (architecture, left_exon_x)
        genome_label(283, y - 2, row['name'])

        # putative TSD bases either side of each junction, shaded on the focal row
        tsd_shown = focal and focal_locus['tsd_adjustment_used'] == '1'
        tsd_back = int(float(focal_locus['tsd_core_back'])) if tsd_shown else 0
        tsd_forward = int(float(focal_locus['tsd_core_fwd'])) if tsd_shown else 0

        def base(x, letter, color=INK, shade=False):
            if shade:
                canvas.add_patch(Rectangle((x - 0.2, y - 4), PITCH, 8, fc=TSD_COLOR, ec='none'))
            text(
                x,
                y - 2,
                letter,
                BASE_FONT,
                fontfamily='DejaVu Sans Mono',
                color=color,
                ha='left',
                va='baseline',
            )

        # exon flanks, muting comparison bases that differ from the focal locus
        for k, letter in enumerate(row['left']):
            differs = not focal and letter != context['left'][k]
            base(
                left_exon_x + k * PITCH,
                letter,
                MUTED if differs else INK,
                tsd_back > 0 and k >= 12 - tsd_back,
            )
        for k, letter in enumerate(row['right']):
            differs = not focal and letter != context['right'][k]
            base(RIGHT_EXON_X + k * PITCH, letter, MUTED if differs else INK, k < tsd_forward)

        # focal row: amino acids of the complete codons in each exon
        if focal:
            phase = phases[example['pid']]['intron_phase']
            codons = []
            for side, exon, start, end, origin in [
                ('left', row['left'], (12 - phase) % 3, 12 - phase, left_exon_x),
                ('right', row['right'], (3 - phase) % 3, 12, RIGHT_EXON_X),
            ]:
                for k in range(start, end - 2, 3):
                    amino_acid = str(Seq(exon[k : k + 3]).translate())
                    assert amino_acid != '*'
                    text(
                        origin + (k + 1.5) * PITCH,
                        y + 6,
                        amino_acid,
                        6.5,
                        fontfamily='DejaVu Sans Mono',
                        color=MUTED,
                        ha='center',
                        va='baseline',
                    )
                    codons.append(
                        dict(exon=side, start0=k, codon=exon[k : k + 3], amino_acid=amino_acid)
                    )
            codon_audit.append(
                dict(
                    pid=example['pid'], architecture=architecture, intron_phase=phase, codons=codons
                )
            )

        # intron: a rule if absent, else bases (coloured F/R on the focal row) and omission counts
        if not intron:
            canvas.plot([intron_x, RIGHT_EXON_X - 0.3 * PITCH], [y] * 2, color=RULE, lw=0.65)
        else:
            marks = repeat_labels(intron, [focal_locus['motif']]) if focal else '.' * len(intron)
            x = intron_x
            for kind, value in intron_tokens:
                if kind == 'g':
                    start, end = value
                    text(
                        x + omission_width / 2 * PITCH,
                        y - 2,
                        str(end - start),
                        6.5,
                        ha='center',
                        color=MUTED,
                        va='baseline',
                    )
                    x += omission_width * PITCH
                else:
                    if marks[value] in 'Ff':
                        color = F_COLOR
                    elif marks[value] in 'Rr':
                        color = R_COLOR
                    else:
                        color = INK
                    in_tsd = value < tsd_forward or value >= len(intron) - tsd_back
                    base(x, intron[value].lower(), color, focal and in_tsd)
                    x += PITCH
            assert abs(x - (RIGHT_EXON_X - PITCH)) < 1e-6

        layout_audit.append(
            dict(
                label_right_x_pt=283,
                label_baseline_pt=y - 2,
                architecture=architecture,
                name=row['name'],
                left_exon_x_pt=left_exon_x,
                right_exon_x_pt=RIGHT_EXON_X,
                group_span_cells=group_span,
                omission_width_cells=omission_width,
                dna_baseline_pt=y - 2,
                base_pitch_pt=PITCH,
                base_font_pt=BASE_FONT,
                y_pt=y,
                sequence_tokens=intron_tokens,
                left=row['left'],
                intron=intron,
                right=row['right'],
            )
        )
    band_top -= band_height
    if architecture != ARCH[-1]:
        canvas.plot([10, CANVAS_WIDTH - 8], [band_top] * 2, color='#e5e5e5', lw=0.4)

# two-row legend under both panels
lineage_y = band_top - 18
sequence_y = band_top - 34
text(10, lineage_y, 'Lineage', 7.5, fontweight='bold')
for j, (clade, name, color) in enumerate(zip(CLADES, NAMES, COLORS)):
    x = 52 + j * 100
    canvas.add_patch(Rectangle((x, lineage_y - 2.3), 6, 4.6, fc=color, ec='none'))
    text(x + 10, lineage_y, display_name(name) + f' ({classified[clade]:,})', 7)
text(10, sequence_y, 'Sequences', 7.5, fontweight='bold')
x = 66
for color, key in [(F_COLOR, 'F: G-rich repeats'), (R_COLOR, 'R: C-rich repeats'), (INK, 'L: linker')]:
    canvas.plot([x, x + 13], [sequence_y] * 2, color=color, lw=1.2)
    text(x + 19, sequence_y, key, 7)
    x += 110
canvas.add_patch(Rectangle((x, sequence_y - 3.5), 13, 7, fc=TSD_COLOR, ec='none'))
text(x + 19, sequence_y, 'Putative TSD', 7)
x += 100
text(x + 6.5, sequence_y, '35', 6.5, color=MUTED, ha='center')
text(x + 19, sequence_y, 'Bases omitted', 7)
x += 100
text(x + 6.5, sequence_y, 'Q E', 6.5, fontfamily='DejaVu Sans Mono', color=MUTED, ha='center')
text(x + 19, sequence_y, 'Amino acids', 7)

# figure, plus layout and translation audits
save(fig, 'Figure2_architecture_examples')
write_qa('sequence_layout.json', layout_audit)
write_qa('amino_acids.json', codon_audit)
