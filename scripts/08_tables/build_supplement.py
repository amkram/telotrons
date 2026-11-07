"""Supplementary Tables S1-S3: TSVs for the docx builder and one xlsx workbook."""

import collections
import csv
import re

from openpyxl import Workbook
from openpyxl.utils import get_column_letter

from common import PACKAGE, REPO, WORK, read_tsv

OUT = PACKAGE / 'supplementary'
TSV_DIR = OUT / 'tables'
TSV_DIR.mkdir(parents=True, exist_ok=True)
# row order of lineages in Table S1
LINEAGE_ORDER = [
    'Haptophytes',
    'Symbiodinium',
    'Baffinella',
    'Chrysophyte',
    'Other stramenopile',
    'Eimeria',
    'Mesorhabditis',
]
SOURCE = {
    'NCBI RefSeq': 'RefSeq',
    'NCBI GenBank': 'GenBank',
    'Tara Oceans SMAGs-v1': 'Tara Oceans SMAGs',
}
INTRONS_FROM = {
    'annotated exons': 'exons',
    'host annotation': 'exons',
    'CDS blocks': 'CDS',
    'annotated exons and CDS blocks': 'exons and CDS',
}
LIBRARY = {
    'RNA-Seq PAIRED': 'paired RNA-seq',
    'RNA-Seq SINGLE': 'single RNA-seq',
    'EST PAIRED': 'paired EST',
}
# display labels of comparator states, position/phase audits and unclassified reasons
STATE = {
    'NO_INTRON': 'intron absent',
    'ORDINARY': 'ordinary intron',
    'NO_ELIGIBLE_COMPARATOR': 'no eligible comparator',
    'HOST_UNRESOLVED': 'host unresolved',
    'TIED_STATES': 'tied states',
    '': 'not evaluated',
}
AUDIT = {
    'not_audited': 'not audited',
    'annotation_corroborated_exact': 'same position and phase, annotation-supported',
    'prediction_only_exact': 'same position and phase, prediction only',
    'position_or_phase_mismatch': 'position or phase mismatch',
    'annotation_or_alignment_unresolved': 'unresolved',
    'exact_alternative_is_telomeric': 'exact match is telomeric',
    'exact_direct_site_near_gap': 'exact match near an assembly gap',
    'prediction_annotation_conflict': 'prediction conflicts with annotation',
    'alternative_boundary_disagreement': 'alternative boundaries disagree',
}
REASON = {
    'repeat_islands_outside_clean_arrays': 'repeat sequence outside the arrays',
    'more_than_two_arrays': 'more than two arrays',
    'separated_arrays_of_same_orientation': 'two same-orientation arrays, separated',
    'interarray_gap_without_qualifying_linker': 'opposing arrays with a gap that is not a linker',
    'degenerate_repeat_extension_leaves_lt18bp_core': (
        'opposing arrays with a gap that is not a linker'
    ),
    'ambiguous_sequence': 'ambiguous bases',
    'no_qualifying_array': 'no array of three units',
}
CLASSIFIED = ['R', 'F', 'R-F', 'F-R', 'R-L-F', 'F-L-R']


def write_rows(name, header, rows):
    with open(TSV_DIR / name, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def read_rows(name):
    with open(TSV_DIR / name) as handle:
        return list(csv.reader(handle, delimiter='\t'))


def decimal(value, places=1):
    """Fixed decimal places, or blank for a missing value."""
    try:
        return f'{float(value):.{places}f}'
    except (TypeError, ValueError):
        return ''


def number(value):
    """An int or float from numeric text, anything else unchanged."""
    try:
        return int(value) if str(value).lstrip('-').isdigit() else float(value)
    except (ValueError, TypeError):
        return value


def display_name(organism):
    """Organism name without screening notes."""
    return re.sub(r'\s*\((no clean sister|non-bearer|genus unassigned)\)', '', organism).strip()


# genome display names and roles from the supplement configuration
genome_config = read_tsv(REPO / 'config/supplement/genomes_v3.tsv')
NAME = {row['Assembly or MAG']: display_name(row['Organism']) for row in genome_config}
assert len(NAME) == 50
assert not any(
    'TARA_' in name or 'non-bearer' in name or 'clean' in name for name in NAME.values()
), NAME
ROLE = {
    row['Assembly or MAG']: 'Bearer' if row['Role'].startswith('Telotron') else 'Control'
    for row in genome_config
}


def table_s1_order(row):
    """Lineage order, bearers before controls, then name."""
    genome_id = row['Assembly or MAG']
    lineage = LINEAGE_ORDER.index(row['Lineage']) if row['Lineage'] in LINEAGE_ORDER else 99
    return (lineage, ROLE[genome_id] != 'Bearer', NAME[genome_id])


# Table S1 inputs: package tables, MAG BUSCO (Delmont et al. 2022) and annotated intron counts
table_s1 = {
    row['genome_id']: row for row in read_tsv(PACKAGE / 'tables/Table_S1_genomes_motifs.tsv')
}
figure_s2 = {row['genome_id']: row for row in read_tsv(PACKAGE / 'tables/FigureS2_genomes.tsv')}
busco = {
    row['genome_id']: row['BUSCO_completion'] for row in read_tsv(PACKAGE / 'inputs/smag_busco.tsv')
}
annotated_introns = {
    row['genome_id']: int(row['annotated_introns'])
    for row in read_tsv(REPO / 'config/supplement/annotated_introns.tsv')
}

# one row per genome in part A (identity, screen) and part B (results)
rows_a = []
rows_b = []
for config in sorted(genome_config, key=table_s1_order):
    genome_id = config['Assembly or MAG']
    package_row = table_s1[genome_id]
    genome = figure_s2[genome_id]
    mag_busco = busco[genome_id] if config['Class'] == 'MAG' else ''
    introns = annotated_introns[genome_id]
    bearer = ROLE[genome_id] == 'Bearer'

    # the supplement configuration must agree with the package tables
    if package_row['BUSCO_percent']:
        assert abs(float(package_row['BUSCO_percent']) - float(mag_busco)) < 0.051, genome_id
    assert package_row['gt85'] == config['>85% loci'] == genome['gt85'], genome_id
    assert package_row['gt50'] == config['>50% loci'] == genome['gt50'], genome_id

    rows_a.append(
        [
            config['Lineage'],
            NAME[genome_id],
            genome_id,
            SOURCE[config['Source']],
            INTRONS_FROM[config['Intron source']],
            ROLE[genome_id],
            config['Telomere motif'],
            f'{introns:,}' if introns else '',
            config['>85% loci'],
            config['>50% loci'],
            decimal(mag_busco) if mag_busco else '',
        ]
    )

    # controls show a dash for bearer-only results
    denominator = int(genome['ortholog_denominator'])
    absent = int(genome['closest_intron_absent'])
    rows_b.append(
        [
            config['Lineage'],
            NAME[genome_id],
            genome_id,
            f'{absent}/{denominator}' if bearer else '—',
            genome['TSD_n'] if bearer else '—',
            decimal(int(genome['eligible_ACGT_bp']) / 1e6, 2),
            decimal(genome['ITS_density'], 3),
            genome['pair_control_id'] if bearer else genome['pair_bearer_id'],
        ]
    )

write_rows(
    'S1_genomes_A.tsv',
    [
        'Lineage',
        'Genome',
        'Accession or MAG',
        'Source',
        'Introns from',
        'Role',
        'Telomere motif',
        'Annotated introns',
        '>85% loci',
        '>50% loci',
        'BUSCO (%)',
    ],
    rows_a,
)
write_rows(
    'S1_genomes_B.tsv',
    [
        'Lineage',
        'Genome',
        'Accession or MAG',
        'Closest ortholog lacks intron (n/compared)',
        'Putative TSDs',
        'Eligible sequence (Mb)',
        'ITS tracts per Mb',
        'Paired genome',
    ],
    rows_b,
)

# Table S3 (package Table S4): RNA-seq junction support per genome, then a total row
table_s4 = {
    row['genome_id']: row for row in read_tsv(PACKAGE / 'tables/Table_S4_RNA_accessions.tsv')
}
baseline = {row['genome_id']: row for row in read_tsv(WORK / '04_rna/ordinary_baseline.tsv')}
COUNT_KEYS = ['primary', 'primary_supported', 'additional', 'additional_supported']
rows_s3 = []
totals = collections.Counter()
for config in read_tsv(REPO / 'config/supplement/rna_v3.tsv'):
    if config['Organism'] == 'Total':
        continue
    genome_id = config['Assembly']
    rna = table_s4[genome_id]
    ordinary = baseline[genome_id]

    # additional loci are the expanded set minus the primary loci
    primary = int(rna['primary_loci'])
    primary_supported = int(rna['primary_ge1'])
    additional = int(rna['expanded_loci']) - primary
    additional_supported = int(rna['expanded_ge1']) - primary_supported
    counts = [primary, primary_supported, additional, additional_supported]
    for key, value in zip(COUNT_KEYS, counts):
        totals[key] += value
    totals['ordinary'] += int(ordinary['ordinary_sampled'])
    totals['ordinary_supported'] += int(ordinary['ordinary_ge1'])

    rows_s3.append(
        [
            NAME[genome_id],
            genome_id,
            config['RNA-seq runs (SRA/ENA)'],
            '; '.join(LIBRARY[library.strip()] for library in config['Library'].split(';')),
            *[f'{value:,}' for value in counts],
            decimal(100 * float(ordinary['ordinary_rate'])),
        ]
    )
rows_s3.append(
    [
        'Total',
        '',
        '',
        '',
        *[f'{totals[key]:,}' for key in COUNT_KEYS],
        decimal(100 * totals['ordinary_supported'] / totals['ordinary']),
    ]
)
write_rows(
    'S3_rna.tsv',
    [
        'Genome',
        'Accession',
        'RNA-seq runs',
        'Library',
        'Primary assayed',
        'Primary supported',
        'Additional assayed',
        'Additional supported',
        'Ordinary introns supported (%)',
    ],
    rows_s3,
)

# the xlsx workbook: one sheet per table
workbook = Workbook()
workbook.remove(workbook.active)


def cell_value(value):
    """Numeric text as a number, blanks and dashes as empty cells."""
    if value in ('', '—', None):
        return None
    text = str(value).replace(',', '')
    if re.fullmatch(r'-?\d+', text):
        return int(text)
    if re.fullmatch(r'-?\d+\.\d+', text):
        return float(text)
    return value


def add_sheet(name, header, rows):
    sheet = workbook.create_sheet(name)
    sheet.append(header)
    for row in rows:
        sheet.append([cell_value(value) for value in row])

    # thousands separators on integers, column widths of 10 to 36 from header length
    for column in sheet.iter_cols(min_row=2):
        for cell in column:
            if isinstance(cell.value, int):
                cell.number_format = '#,##0'
    for index, title in enumerate(header, 1):
        sheet.column_dimensions[get_column_letter(index)].width = min(max(len(title) + 2, 10), 36)
    sheet.freeze_panes = 'A2'


# S1 sheet: part A, then part B without its three identity columns and with the
# n/compared column split into two numeric columns
part_a = read_rows('S1_genomes_A.tsv')
part_b = read_rows('S1_genomes_B.tsv')
header_b = part_b[0][3:]
split_at = header_b.index('Closest ortholog lacks intron (n/compared)')
header_b = (
    header_b[:split_at] + ['Intron absent (loci)', 'Compared loci'] + header_b[split_at + 1 :]
)


def split_fraction(row):
    values = row[3:]
    absent, compared = values[split_at].split('/') if '/' in values[split_at] else ('', '')
    return values[:split_at] + [absent, compared] + values[split_at + 1 :]


add_sheet(
    'S1_genomes',
    part_a[0] + header_b,
    [row_a + split_fraction(row_b) for row_a, row_b in zip(part_a[1:], part_b[1:])],
)


def unclassified_reason(locus):
    """The first architecture rule an unclassified locus fails (blank when classified)."""
    if locus['architecture'] in CLASSIFIED:
        return ''
    first = next((reason for reason in locus['architecture_reason'].split(';') if reason), '')
    if first.endswith('not_within_8bp_of_intron_or_tsd_boundary'):
        return 'array ends >8 bp from an intron boundary'
    if first.endswith('requires_strictly_gt85_percent_coverage'):
        return 'single or convergent arrays at ≤85% coverage'
    return REASON.get(first, first.replace('_', ' '))


# S2 sheet: Table S2 with display labels and the stage 5 count of eligible comparators
eligible = {
    row['locus_id']: row['eligible_comparators']
    for row in read_tsv(WORK / '05_orthology/results/closest_comparators.tsv')
}
table_s2 = read_tsv(PACKAGE / 'tables/Table_S2_telotrons.tsv.gz')


def s2_row(locus):
    # TSD fields only for putative TSDs, RNA reads only where assayed
    tsd = locus['putative_TSD'] == '1'
    architecture = locus['architecture']
    return [
        locus['locus_id'],
        NAME[locus['genome_id']],
        locus['seqid'],
        number(locus['start']),
        number(locus['end']),
        locus['strand'],
        number(locus['intron_len']),
        locus['splice_class'],
        round(float(locus['repeat_fraction']), 4),
        'unclassified' if architecture == 'mixed_unclear' else architecture.replace('-', '–'),
        unclassified_reason(locus),
        number(locus['linker_start']),
        number(locus['linker_end']),
        number(locus['tsd_core_bp']) if tsd else '',
        locus['tsd_core_sequence'] if tsd else '',
        number(locus['rna_exact']) if locus['rna_status'] == 'assayed' else '',
        STATE[locus['closest_state']],
        number(eligible.get(locus['locus_id'], '')),
        AUDIT.get(locus['exact_phase_audit'], 'not audited'),
    ]


add_sheet(
    'S2_catalogue',
    [
        'Locus',
        'Genome',
        'Contig',
        'Start',
        'End',
        'Strand',
        'Intron length (bp)',
        'Splice boundaries',
        'Repeat coverage',
        'Architecture',
        'Why unclassified',
        'Linker start',
        'Linker end',
        'TSD length (bp)',
        'TSD sequence',
        'RNA junction reads',
        'Closest comparator',
        'Eligible comparators',
        'Position and phase audit',
    ],
    [s2_row(locus) for locus in table_s2],
)

# S3 sheet straight from its TSV, then save the workbook
rows = read_rows('S3_rna.tsv')
add_sheet('S3_rna', rows[0], rows[1:])
workbook.save(OUT / 'SUPPLEMENTARY_TABLES.xlsx')
print('sheets', workbook.sheetnames, 'S2 rows', len(table_s2))
