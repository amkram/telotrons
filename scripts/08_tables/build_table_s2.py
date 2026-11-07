"""Table S2: every catalogue locus with its evidence and putative TSD call."""

from common import PACKAGE, WORK, read_tsv, write_tsv

# a putative TSD: a flanking duplication core scoring at least 6
MINIMUM_TSD_SCORE = 6

loci = read_tsv(PACKAGE / 'results/catalogue_with_evidence.tsv.gz')
tsd = {row['locus_id']: row for row in read_tsv(WORK / '03_tsd/TSD_loci.tsv')}

table = []
for locus in loci:
    flag = str(int(int(locus['tsd_core_score']) >= MINIMUM_TSD_SCORE))
    stage3 = tsd.get(locus['locus_id'])
    # stage 3 calls primary loci by the same rule
    assert stage3 is None or stage3['putative_TSD'] == flag, locus['locus_id']

    # putative_TSD goes after the catalogue columns, before the joined evidence
    row = {}
    for key, value in locus.items():
        if key == 'rna_exact':
            row['putative_TSD'] = flag
        row[key] = value
    table.append(row)

write_tsv(PACKAGE / 'tables/Table_S2_telotrons.tsv.gz', table)
print('Table S2:', len(table), 'loci,', sum(row['putative_TSD'] == '1' for row in table), 'TSDs')
