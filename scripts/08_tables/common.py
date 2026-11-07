"""Paths and TSV helpers shared by the stage 8 table scripts."""

import csv
import gzip
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WORK = REPO / 'work'
PACKAGE = Path(os.environ['TELO_PACKAGE'])


def read_tsv(path):
    """Rows of a tab-separated file, gzipped or not, as dictionaries."""
    with gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows):
    """Write dictionaries as a tab-separated file, gzipped when the name ends in .gz."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt') if str(path).endswith('.gz') else open(path, 'w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
