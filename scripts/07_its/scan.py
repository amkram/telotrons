#!/usr/bin/env python3
"""ITS tract scan of one genome of config/genomes.tsv (bearer or control).

Components: merged tracts of >= 4 motif units with <= 1 substitution per adjacent unit pair,
outside the whole-model gene mask and the terminal telomeric caps. Writes
work/07_its/its_raw/per_genome/<genome>.components.tsv.gz, .caps.tsv.gz and .json (count and
eligible A/C/G/T bp).
"""

import argparse, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'lib'))
from _common import find_genome_files
from audit_its_opportunity import rows_by_genome
import its_rescan, raw_scanner

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--genome-id', required=True)
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/07_its'))
args = parser.parse_args()

# this genome's lineage motif and its assembly and annotation
meta = rows_by_genome(args.config)[args.genome_id]
fasta_path, gff_path = find_genome_files(args.genome_id, args.refseq_dir, args.tara_dir)
assert fasta_path and gff_path, (args.genome_id, 'missing assembly or annotation')

# scan and report the count and density
folder = Path(args.out) / 'its_raw'
core = its_rescan.build(folder)
result = raw_scanner.run_genome(
    args.genome_id, meta['motif'], fasta_path, gff_path, folder / 'per_genome', core
)
print(
    'ITS',
    args.genome_id,
    meta['motif'],
    result['n_components'],
    result['eligible_acgt_bp'],
    result['components_per_acgt_mb'],
    flush=True,
)
