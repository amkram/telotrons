"""Exact position audit paths and helpers (reads TELO_WORKDIR, writes work/06_exact_position)."""

import csv
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / 'lib'), str(REPO / 'scripts/05_orthology')]

from _common import find_genome_files, rc  # noqa: E402
from indexed_reference import certified_cache, read_fai, read_range  # noqa: E402

OUT = REPO / 'work/06_exact_position'
ORTHOLOGY = Path(os.environ['TELO_WORKDIR']).resolve()  # stage 5 work directory
RNA_LOCI = REPO / 'work/04_rna/RNA_loci.tsv'
PRIMARY_TIER = 'standard'  # catalogue tier of the primary loci
RESOLVED_STATES = ('NO_INTRON', 'ORDINARY', 'TELOTRON', 'MIXED', 'AMBIGUOUS')
ELIGIBLE_STATES = ('NO_INTRON', 'ORDINARY')


def sha(path):
    """SHA-256 of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_scratch():
    """TMPDIR, TMP and TEMP must name one job scratch directory outside /tmp."""
    scratch = Path(os.environ['TMPDIR']).resolve(strict=True)
    assert not str(scratch).startswith(('/tmp', '/private/tmp')), scratch
    assert os.environ['TMP'] == os.environ['TEMP'] == str(scratch)
    return scratch


def read_tsv(path):
    """Rows of a (possibly gzipped) tab-separated table."""
    with gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(name, rows):
    """Write rows to OUT/name, with the first row's keys as columns."""
    with (OUT / name).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_json_gz(name):
    """A gzipped JSON file in OUT."""
    with gzip.open(OUT / name, 'rt') as handle:
        return json.load(handle)


def write_json_gz(name, value):
    """Write a gzipped JSON file in OUT."""
    with gzip.open(OUT / name, 'wt') as handle:
        json.dump(value, handle)


def primary_closest():
    """Closest-comparator rows of the primary loci (stage 5 results)."""
    rows = read_tsv(ORTHOLOGY / 'results/closest_comparators.tsv')
    return [row for row in rows if row['tier'] == PRIMARY_TIER]


def genome_files(genome_id):
    """(FASTA, GFF or None) of a genome under data/raw."""
    raw = REPO / 'data/raw'
    return find_genome_files(genome_id, str(raw / 'refseq'), str(raw / 'tara'))


class GenomeSequence:
    """Contigs of one genome, from its stage 5 reference cache when one exists, else the FASTA."""

    def __init__(self, source, contigs):
        self.handle = None
        self.seqs = {}
        cached = certified_cache(source)
        if cached:
            fasta, certificate = cached
            self.entries = read_fai(fasta)
            self.handle = fasta.open('rb')
            self.provenance = dict(
                cache=str(fasta), certificate=str(certificate), certificate_sha256=sha(certificate)
            )
            return

        source = Path(source).resolve()
        name = None
        with gzip.open(source, 'rt') if source.name.endswith('.gz') else source.open() as handle:
            for line in handle:
                if line.startswith('>'):
                    name = line[1:].split()[0]
                elif name in contigs:
                    self.seqs[name] = self.seqs.get(name, '') + line.strip().upper()
        self.provenance = dict(source=str(source), sha256=sha(source))

    def fetch(self, contig, start, end):
        """Bases [start, end) of a contig, clipped to the contig."""
        if not self.handle:
            return self.seqs[contig][max(0, start) : end]
        entry = self.entries[contig]
        start = max(0, start)
        end = min(entry[0], end)
        if end <= start:
            return ''
        return read_range(self.handle, entry, start, end)

    def fetch_stranded(self, contig, start, end, strand):
        """Bases [start, end) read on the given strand."""
        seq = self.fetch(contig, start, end)
        return rc(seq) if strand == '-' else seq

    def close(self):
        if self.handle:
            self.handle.close()
