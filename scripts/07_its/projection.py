#!/usr/bin/env python3
"""Projected gene masks and ITS recount for one genome of config/genomes.tsv (bearer or control).

Donor candidates and the marker threshold come from config/its_donors.tsv. Donor genomes go under
data/raw like any other genome. Each candidate's proteome and 80 protein markers
are derived from its annotation; the three candidates whose markers align best are projected onto
the genome with miniprot, then every other eligible candidate. The recount removes raw-scan
components (scan.py) that overlap a projected gene span. Bearers' known telotrons, read from
--catalogue, check that no counted component overlaps one.
Writes work/07_its/proj/<genome>/summary.tsv (components and eligible A/C/G/T bp).
"""

import argparse, json, shutil, sys, tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'lib'))
from _common import find_genome_files, scratch_root
from audit_its_opportunity import read_tsv, rows_by_genome
from prepare_proteomes import write_proteome
from recount import recount
import run_projection

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--genome-id', required=True)
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--donor-table', default=str(REPO / 'config/its_donors.tsv'))
parser.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/07_its'))
parser.add_argument('--threads', type=int, default=4)
args = parser.parse_args()


def genome_files(genome_id):
    """Assembly and annotation of one genome under data/raw."""
    fasta_path, gff_path = find_genome_files(genome_id, args.refseq_dir, args.tara_dir)
    assert fasta_path and gff_path, (genome_id, 'missing assembly or annotation under data/raw')
    return fasta_path, gff_path


# this genome's row, its donor candidates and the marker threshold
genome = rows_by_genome(args.config)[args.genome_id]
donor_row = rows_by_genome(args.donor_table)[args.genome_id]
candidates = donor_row['donor_candidates'].split(';')
min_markers = int(donor_row['min_marker_hits'])
folder = Path(args.out) / 'proj' / args.genome_id
proteomes = folder / 'proteomes'
folder.mkdir(parents=True, exist_ok=True)

# every candidate's proteome and markers
for donor in candidates:
    write_proteome(donor, *genome_files(donor), proteomes)

# index the genome in job scratch (TMPDIR, never /tmp) with the pinned miniprot
run_projection.check_version()
fasta_path, gff_path = genome_files(args.genome_id)
scratch = Path(tempfile.mkdtemp(prefix=f'its_{args.genome_id}_', dir=scratch_root()))
index = scratch / 'genome.mpi'
run_projection.index_genome(fasta_path, index, args.threads, folder / 'index.log')

# choose donors by marker alignment
selection = run_projection.select_donors(
    candidates, proteomes, index, folder, scratch, args.threads, min_markers
)
(folder / 'donor_selection.json').write_text(json.dumps(selection, indent=2) + '\n')
print(args.genome_id, 'donors', selection['selected'], selection['additional'], flush=True)

# project the selected donors, then every other eligible donor
projections = [('proteomes', selection['selected'])]
if selection['additional']:
    projections.append(('additional_proteomes', selection['additional']))
paths = [
    run_projection.project(donors, proteomes, index, folder, name, scratch, args.threads)
    for name, donors in projections
]
shutil.rmtree(scratch)

# known telotrons of this genome, 0-based half-open (controls have none)
known = [
    (row['seqid'], int(row['start']) - 1, int(row['end']))
    for row in read_tsv(args.catalogue)
    if row['genome_id'] == args.genome_id
]

# recount the raw-scan components outside the projected gene spans
raw_folder = Path(args.out) / 'its_raw/per_genome'
recount(genome, fasta_path, gff_path, raw_folder, folder, paths, known)
