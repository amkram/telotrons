"""Indexed genome copies in reference_cache/, read by byte offset to keep panels out of RAM."""

import fcntl
import gzip
import hashlib
import json
import os
import subprocess
from collections import OrderedDict
from pathlib import Path

from common import WORK, dump, use_scratch

CACHE = WORK / 'reference_cache'
MAX_OPEN_FILES = 16
open_files = OrderedDict()  # least recently used first


def source_signature(source):
    """Path, size and mtime of a source FASTA, and the cache directory they key."""
    path = Path(source).resolve(strict=True)
    stat = path.stat()
    signature = dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:24]
    return signature, CACHE / key


def certified_cache(source):
    """(cached FASTA, certificate) of a source FASTA, or None when no cache was built."""
    signature, directory = source_signature(source)
    certificate = directory / 'complete.json'
    if not certificate.exists():
        return None
    assert json.loads(certificate.read_text())['source'] == signature, certificate
    return directory / 'reference.fa', certificate


def read_fai(fasta):
    """samtools .fai index: contig -> (length, offset, bases per line, bytes per line)."""
    lines = Path(str(fasta) + '.fai').read_text().splitlines()
    return {fields[0]: tuple(map(int, fields[1:5])) for fields in (x.split('\t') for x in lines)}


def read_range(handle, entry, start, end):
    """Bases [start, end) of one contig of an open fixed-width FASTA, given its .fai entry."""
    _, offset, bases_per_line, bytes_per_line = entry
    first = offset + start // bases_per_line * bytes_per_line + start % bases_per_line
    last = offset + (end - 1) // bases_per_line * bytes_per_line + (end - 1) % bases_per_line + 1
    handle.seek(first)
    seq = handle.read(last - first).replace(b'\n', b'').replace(b'\r', b'').decode().upper()
    assert len(seq) == end - start, (start, end, len(seq))
    return seq


def prepare(source):
    """Cached, indexed copy of a FASTA, built once per source while other jobs wait on the lock."""
    use_scratch()
    signature, directory = source_signature(source)
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / 'reference.fa'
    done = directory / 'complete.json'

    with (directory / 'build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if done.exists():
            assert json.loads(done.read_text())['source'] == signature
            assert fasta.exists() and Path(str(fasta) + '.fai').exists()
            return fasta

        # copy (decompressing .gz) to a partial file, then rename it into place
        path = Path(signature['path'])
        partial = directory / f'reference.{os.getpid()}.partial'
        opener = gzip.open if path.name.endswith('.gz') else open
        with opener(path, 'rb') as source_file, partial.open('wb') as out:
            for block in iter(lambda: source_file.read(1024 * 1024), b''):
                out.write(block)

        # the source must not change while it is copied
        stat = path.stat()
        assert stat.st_size == signature['size'] and stat.st_mtime_ns == signature['mtime_ns']
        os.replace(partial, fasta)
        subprocess.run(['samtools', 'faidx', str(fasta)], check=True)
        dump(done, dict(source=signature, reference=str(fasta), complete=True))
    return fasta


class Reference:
    """Contigs of one cached genome. build=False reads an existing cache and never builds."""

    def __init__(self, source, build=True):
        if build:
            self.path = prepare(source)
        else:
            cached = certified_cache(source)
            assert cached, ('no reference cache', source)
            self.path = cached[0]
        self.entries = read_fai(self.path)
        self.lengths = {name: entry[0] for name, entry in self.entries.items()}

    def get(self, name, default=None):
        return Contig(self, name, self.lengths[name]) if name in self.lengths else default

    def __getitem__(self, name):
        return Contig(self, name, self.lengths[name])

    def fetch(self, name, start, end):
        """Bases [start, end) of one contig, keeping at most MAX_OPEN_FILES files open."""
        path = str(self.path)
        if path not in open_files:
            while len(open_files) >= MAX_OPEN_FILES:
                open_files.popitem(last=False)[1].close()
            open_files[path] = open(path, 'rb')
        open_files.move_to_end(path)
        if start >= end:
            return ''
        return read_range(open_files[path], self.entries[name], start, end)


class Contig:
    """One contig of a Reference, sliced like a string."""

    def __init__(self, reference, name, length):
        self.reference = reference
        self.name = name
        self.length = length

    def __len__(self):
        return self.length

    def __getitem__(self, index):
        if isinstance(index, slice):
            start, end, step = index.indices(self.length)
            assert step == 1
            return self.reference.fetch(self.name, start, end)
        if index < 0:
            index += self.length
        return self.reference.fetch(self.name, index, index + 1)


def miniprot_reference(source, threads):
    """miniprot index of a cached genome, built once under a lock."""
    fasta = prepare(source)
    index = fasta.with_suffix('.mpi')
    with (fasta.parent / 'miniprot.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not index.exists():
            partial = Path(str(index) + f'.{os.getpid()}.partial')
            with (fasta.parent / 'miniprot_index.log').open('w') as log:
                command = ['miniprot', '-t', str(threads), '-d', str(partial), str(fasta)]
                subprocess.run(command, stdout=log, stderr=log, check=True)
            os.replace(partial, index)
    return str(index)
