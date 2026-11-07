#!/usr/bin/env python3
"""Map one genome's RNA runs (config/rna_runs.tsv) and count exact junction reads at each of its
catalogue loci.

Writes work/04_rna/<genome>/<genome>.splice_sr.bam (unless --bam names an existing BAM) and
work/04_rna/<genome>/junctions.tsv. scripts/04c_collect_rna.py joins all genomes.
"""

import argparse, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))
from _common import find_genome_files
from audit_its_opportunity import read_tsv, rows_by_genome
from rnaseq_map import map_runs
from telotron_splice_junctions import count_junctions

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--genome-id', required=True)
parser.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
parser.add_argument('--config', default=str(REPO / 'config/rna_runs.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/04_rna'))
parser.add_argument('--threads', type=int, default=12)
parser.add_argument('--max-size', default='120g', help='prefetch download limit per run')
parser.add_argument('--bam', help='existing BAM (with .bai) to count from instead of mapping')
args = parser.parse_args()

# this genome's runs and catalogue loci
runs = rows_by_genome(args.config)[args.genome_id]
loci = [row for row in read_tsv(args.catalogue) if row['genome_id'] == args.genome_id]
assert loci, (args.genome_id, 'no catalogue loci')
folder = Path(args.out) / args.genome_id
bam = args.bam or folder / f'{args.genome_id}.splice_sr.bam'

# map the runs unless an existing BAM was given
if not args.bam:
    fasta_path, gff_path = find_genome_files(args.genome_id, args.refseq_dir, args.tara_dir)
    assert fasta_path, (args.genome_id, 'missing assembly under data/raw')
    map_runs(
        fasta_path,
        runs['runs'].split(';'),
        runs['preset'],
        bam,
        args.genome_id,
        threads=args.threads,
        max_size=args.max_size,
    )

# exact junction reads per locus
rows = count_junctions(bam, loci, folder / 'junctions.tsv')
expressed = sum(row[3] > 0 for row in rows)
supported = sum(row[4] > 0 for row in rows)
print(
    args.genome_id, len(rows), 'loci;', expressed, 'with reads;', supported, 'with junction reads'
)
