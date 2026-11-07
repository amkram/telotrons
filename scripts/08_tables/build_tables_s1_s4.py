"""Table S1 (genomes, lineage motif, locus counts) and Table S4 (RNA-seq runs, junction support)."""

from common import PACKAGE, REPO, WORK, read_tsv, write_tsv

# Table S1: the genome configuration with the Figure S2 locus counts
figure_s2 = {row['genome_id']: row for row in read_tsv(PACKAGE / 'tables/FigureS2_genomes.tsv')}
table_s1 = []
for genome in read_tsv(REPO / 'config/genomes.tsv'):
    counts = figure_s2[genome['genome_id']]
    bearer = genome['role'] == 'bearer'
    table_s1.append(
        dict(
            genome_id=genome['genome_id'],
            label=genome['label'],
            role=genome['role'],
            reporting_lineage=genome['lineage'],
            assembly_class=genome['assembly_class'],
            tested_motifs=genome['motif'],
            control_scope='not_applicable' if bearer else 'single_lineage_motif_screen',
            BUSCO_percent=genome['BUSCO_percent'],
            assembly_source=genome['assembly_source'],
            annotation_provider=genome['annotation_provider'],
            gt85=counts['gt85'],
            gt50=counts['gt50'],
        )
    )
write_tsv(PACKAGE / 'tables/Table_S1_genomes_motifs.tsv', table_s1)

# Table S4: assayed loci and junction support per genome with RNA-seq
assayed = [locus for locus in read_tsv(WORK / '04_rna/RNA_loci.tsv') if locus['assayed'] == '1']
table_s4 = []
for run in read_tsv(REPO / 'config/rna_runs.tsv'):
    loci = [locus for locus in assayed if locus['genome_id'] == run['genome_id']]
    primary = [locus for locus in loci if locus['tier'] == 'standard']
    table_s4.append(
        dict(
            genome_id=run['genome_id'],
            organism=run['organism'],
            runs=run['runs'],
            preset=run['preset'],
            reference='genomic assembly',
            primary_loci=len(primary),
            primary_ge1=sum(locus['ge1'] == '1' for locus in primary),
            primary_ge3=sum(locus['ge3'] == '1' for locus in primary),
            expanded_loci=len(loci),
            expanded_ge1=sum(locus['ge1'] == '1' for locus in loci),
            expanded_ge3=sum(locus['ge3'] == '1' for locus in loci),
        )
    )
write_tsv(PACKAGE / 'tables/Table_S4_RNA_accessions.tsv', table_s4)
print('Table S1:', len(table_s1), 'genomes; Table S4:', len(table_s4), 'genomes with RNA-seq')
