"""Figure 1, Figure S2 and Table S3 source tables, plus the curated and ITS inputs of stage 9."""

import collections
import json
import shutil
import statistics

from common import PACKAGE, REPO, WORK, read_tsv, write_tsv

TABLES = PACKAGE / 'tables'
# Figure 1 names of the bearer groups
BEARER_LABELS = {
    'nematode': 'Mesorhabditis',
    'cryptophyte': 'Baffinella',
    'haptophyte': 'Haptophytes',
    'chrysophyte': 'Chrysophyte',
    'stramenopile': 'Other stramenopile',
    'dinoflagellate': 'Symbiodinium',
    'apicomplexa': 'Eimeria',
}
# the two stramenopile controls have one Figure 1 row each
STRAMENOPILE_CONTROL_LABELS = {
    'GCA_001750085.1': 'Diatom controls',
    'GCF_000186865.1': 'Pelagophyte controls',
}


def median_range(values):
    """Median, minimum and maximum, or blanks for no values."""
    return (statistics.median(values), min(values), max(values)) if values else ('', '', '')


def figure1_row(row_id, label, genome_rows):
    """One Figure 1 row pooling the Figure S2 rows of its genomes."""
    primary_loci = sum(row['gt85'] for row in genome_rows)
    absent = sum(row['closest_intron_absent'] for row in genome_rows)
    denominator = sum(row['ortholog_denominator'] for row in genome_rows)
    tsds = sum(row['TSD_n'] for row in genome_rows)

    # the motif behind the most primary loci of the group
    loci_by_motif = collections.Counter()
    for row in genome_rows:
        loci_by_motif[row['motif']] += row['gt85']
    dominant_motif = loci_by_motif.most_common(1)[0][0] if primary_loci else '—'

    return dict(
        row_id=row_id,
        label=f'{label} (n = {len(genome_rows)})',
        clade=genome_rows[0]['clade'],
        role=genome_rows[0]['role'],
        genomes=len(genome_rows),
        genome_ids=json.dumps([row['genome_id'] for row in genome_rows]),
        dominant_motif=dominant_motif,
        gt85=primary_loci,
        gt50=sum(row['gt50'] for row in genome_rows),
        closest_intron_absent=absent,
        ordinary_exact=denominator - absent,
        ortholog_denominator=denominator,
        intron_absent_percent=100 * absent / denominator if denominator else '',
        putative_TSDs_gt85=tsds,
        TSD_percent_gt85=100 * tsds / primary_loci if primary_loci else '',
    )


# curated figure inputs (topologies, exon phases, Figure 2B examples, MAG BUSCO) for stages 8-9
(PACKAGE / 'inputs').mkdir(parents=True, exist_ok=True)
for source in sorted((REPO / 'config/figures').iterdir()):
    shutil.copy2(source, PACKAGE / 'inputs' / source.name)

# Figure 3 inputs from stage 7 (the ITS points are also Table S5)
TABLES.mkdir(parents=True, exist_ok=True)
for name in ['figure_points.tsv', 'pairs.tsv', 'exact_wilcoxon.json']:
    shutil.copy2(WORK / '07_its' / name, TABLES / name)
shutil.copy2(WORK / '07_its/figure_points.tsv', TABLES / 'Table_S5_ITS.tsv')

evidence = read_tsv(PACKAGE / 'results/catalogue_with_evidence.tsv.gz')
tsd = {row['locus_id']: row for row in read_tsv(WORK / '03_tsd/TSD_loci.tsv')}
genomes = read_tsv(REPO / 'config/genomes.tsv')
# lineages in the order they first appear in the genome configuration
lineages = list(dict.fromkeys(genome['lineage'] for genome in genomes))
its_points = {row['genome_id']: row for row in read_tsv(WORK / '07_its/figure_points.tsv')}
assert set(its_points) == {genome['genome_id'] for genome in genomes}

# Table S3: the reported closest-comparator state of each primary locus
primary_by_genome = collections.defaultdict(list)
comparisons = collections.defaultdict(collections.Counter)
table_s3 = []
for locus in evidence:
    if locus['tier'] != 'standard':
        continue
    genome_id = locus['genome_id']
    state = locus['closest_state']
    audit = locus['exact_phase_audit']
    primary_by_genome[genome_id].append(locus)

    # only absent introns and annotation-supported ordinary introns enter the denominator
    if state == 'NO_INTRON':
        reported = 'intron_absent'
    elif audit == 'annotation_corroborated_exact':
        reported = 'ordinary_exact'
    else:
        reported = 'excluded'
    comparisons[genome_id][reported] += 1
    comparisons[genome_id]['prediction_only'] += audit == 'prediction_only_exact'

    table_s3.append(
        dict(
            locus_id=locus['locus_id'],
            genome_id=genome_id,
            lineage=locus['clade'],
            tier='primary',
            closest_state=state,
            position_phase_audit=audit,
            reported_state=reported,
            denominator_included=int(reported != 'excluded'),
            rna_exact_junction_reads=locus['rna_exact'],
            rna_status=locus['rna_status'],
        )
    )
write_tsv(TABLES / 'Table_S3_primary_orthology.tsv', table_s3)

# Figure S2: one row per bearer and control genome
loci_per_genome = collections.Counter(locus['genome_id'] for locus in evidence)
genome_rows = []
for genome in genomes:
    genome_id = genome['genome_id']
    bearer = genome['role'] == 'bearer'
    primary = primary_by_genome[genome_id]
    counts = comparisons[genome_id]
    denominator = counts['intron_absent'] + counts['ordinary_exact']
    intron_lengths = median_range([int(locus['intron_len']) for locus in primary])

    # lengths of the putative TSDs at this genome's primary loci
    tsd_lengths = [
        int(tsd[locus['locus_id']]['core_bp'])
        for locus in primary
        if tsd[locus['locus_id']]['putative_TSD'] == '1'
    ]
    tsd_length_range = median_range(tsd_lengths)

    its = its_points[genome_id]
    genome_rows.append(
        dict(
            genome_id=genome_id,
            label=genome['label'],
            clade=genome['lineage'],
            assembly_class=genome['assembly_class'],
            role=genome['role'],
            motif=genome['motif'] if bearer else '—',
            gt85=len(primary),
            gt50=loci_per_genome[genome_id],
            intron_median_bp=intron_lengths[0],
            intron_min_bp=intron_lengths[1],
            intron_max_bp=intron_lengths[2],
            closest_intron_absent=counts['intron_absent'],
            ordinary_exact=counts['ordinary_exact'],
            ordinary_prediction_only=counts['prediction_only'],
            ortholog_denominator=denominator,
            ortholog_included_percent=100 * denominator / len(primary) if bearer else '',
            intron_absent_percent=(
                100 * counts['intron_absent'] / denominator if denominator else ''
            ),
            TSD_n=len(tsd_lengths) if bearer else 0,
            TSD_percent=100 * len(tsd_lengths) / len(primary) if bearer else '',
            TSD_median_bp=tsd_length_range[0],
            TSD_min_bp=tsd_length_range[1],
            TSD_max_bp=tsd_length_range[2],
            ITS_components=int(its['ITS_components']),
            eligible_ACGT_bp=int(its['eligible_ACGT_bp']),
            ITS_density=float(its['density']),
            BUSCO_percent=genome['BUSCO_percent'],
            pair_bearer_id=genome['pair_bearer_id'],
            pair_control_id=genome['pair_control_id'],
        )
    )

# lineage order, bearers before controls, then by name
genome_rows.sort(
    key=lambda row: (
        lineages.index(row['clade']),
        row['role'] != 'bearer',
        row['label'],
        row['genome_id'],
    )
)
write_tsv(TABLES / 'FigureS2_genomes.tsv', genome_rows)

# Figure 1: one row per lineage and role, except one row per stramenopile control
figure1 = []
for lineage in lineages:
    for role in ['bearer', 'control']:
        group = [row for row in genome_rows if row['clade'] == lineage and row['role'] == role]
        if not group or (lineage == 'stramenopile' and role == 'control'):
            continue
        label = BEARER_LABELS[lineage] if role == 'bearer' else 'Controls'
        figure1.append(figure1_row(f'{lineage}_{role}', label, group))
for row in genome_rows:
    if row['clade'] == 'stramenopile' and row['role'] == 'control':
        label = STRAMENOPILE_CONTROL_LABELS[row['genome_id']]
        figure1.append(figure1_row(row['genome_id'], label, [row]))
write_tsv(TABLES / 'Figure1_source_data.tsv', figure1)

# intron-absent fraction of the closest comparisons, per bearer group
for row in figure1:
    if row['role'] == 'bearer':
        print(row['row_id'], row['closest_intron_absent'], row['ortholog_denominator'])
absent = sum(row['closest_intron_absent'] for row in figure1)
denominator = sum(row['ortholog_denominator'] for row in figure1)
print(f'intron absent {absent}/{denominator} ({100 * absent / denominator:.1f}%)')
