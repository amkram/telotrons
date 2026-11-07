"""Write TELO_WORKDIR/inputs: primary.tsv, additional.tsv, genomes.json and scan_genomes.json."""

import argparse
import csv
import gzip

from common import REPO, WORK, dump, write
from _common import find_genome_files

TIER_FILES = {'standard': 'primary.tsv', 'additional': 'additional.tsv'}


def read_config(name):
    """Rows of a tab-separated table in config/."""
    with (REPO / 'config' / name).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def genome_entries(catalogue, refseq_dir, tara_dir):
    """genomes.json entries: files, clade, organism, assembly class and comparator panel."""
    assembly_class = {row['genome_id']: row['assembly_class'] for row in read_config('genomes.tsv')}
    labels = {}
    for row in catalogue:
        labels.setdefault(row['genome_id'], set()).add((row['clade'], row['organism']))

    entries = []
    for panel in read_config('comparator_panels.tsv'):
        genome_id = panel['focal_genome_id']
        assert len(labels.get(genome_id, ())) == 1, (genome_id, labels.get(genome_id))
        clade, organism = next(iter(labels[genome_id]))
        assert clade == panel['clade'], (genome_id, clade, panel['clade'])
        comparators = panel['comparators'].split(';')
        assert len(comparators) == int(panel['n_comparators']), genome_id

        # the same lookup the aligner uses for every focal genome
        fasta, gff = find_genome_files(genome_id, refseq_dir, tara_dir)
        assert fasta and gff, (genome_id, 'genome FASTA or GFF not found under data/raw')
        entries.append(
            dict(
                genome_id=genome_id,
                fasta=fasta,
                gff=gff,
                clade=clade,
                organism=organism,
                assembly_class=assembly_class[genome_id],
                comparators=comparators,
            )
        )
    return entries


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
    args = ap.parse_args()

    with gzip.open(args.catalogue, 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        catalogue = list(reader)
        fields = reader.fieldnames

    # catalogue rows by tier, all columns kept
    assert {row['tier'] for row in catalogue} <= set(TIER_FILES)
    for tier, name in TIER_FILES.items():
        rows = [row for row in catalogue if row['tier'] == tier]
        write(WORK / 'inputs' / name, rows, fields)
        print(name, len(rows), 'loci', flush=True)

    raw = REPO / 'data/raw'
    entries = genome_entries(catalogue, str(raw / 'refseq'), str(raw / 'tara'))
    dump(WORK / 'inputs/genomes.json', entries)

    # telomere motifs and group per genome, read by the state caller for its focal genome
    scan = [
        dict(genome_id=row['genome_id'], group=row['group'], motifs=row['motifs'].split(';'))
        for row in read_config('comparator_motifs.tsv')
    ]
    dump(WORK / 'inputs/scan_genomes.json', scan)
    print('genomes.json', len(entries), 'bearers; scan_genomes.json', len(scan), 'genomes')


if __name__ == '__main__':
    main()
