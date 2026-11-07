"""Paths, task selection and table helpers shared by the orthology scripts (stage 5)."""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

WORK = Path(os.environ['TELO_WORKDIR']).resolve()  # orthology work directory
REPO = Path(__file__).resolve().parents[2]
STAGES = ('primary', 'additional')  # inputs/<stage>.tsv: catalogue tiers standard, additional
RECALL_POLICY = 'actual_matches12_of20_v1'
MIN_FLANK_MATCHES = 12  # matching residues required in each 20-residue flank (recall policy)

# the aligner, the state caller and the shared sequence helpers live in lib/
sys.path.insert(0, str(REPO / 'lib'))


def read(path):
    """Rows of a tab-separated table."""
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def header(path):
    """Column names of a tab-separated table."""
    with Path(path).open() as handle:
        return next(csv.reader(handle, delimiter='\t'))


def write(path, rows, fields=None):
    """Write rows as a tab-separated table (fields default to the first row's keys)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def dump(path, value):
    """Write a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def sha(path):
    """SHA-256 of a file."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def use_scratch():
    """Refuse /tmp and route temporary files and tool caches into the job's TMPDIR."""
    scratch = Path(os.environ['TMPDIR']).resolve(strict=True)
    assert scratch != Path('/tmp'), scratch
    os.environ['TMP'] = os.environ['TEMP'] = str(scratch)
    for variable, subdir in [
        ('MAFFT_TMPDIR', 'mafft'),
        ('PIP_CACHE_DIR', 'pip'),
        ('MPLCONFIGDIR', 'mpl'),
        ('XDG_CACHE_HOME', 'cache'),
        ('JOBLIB_TEMP_FOLDER', 'joblib'),
        ('NUMBA_CACHE_DIR', 'numba'),
    ]:
        target = scratch / subdir
        target.mkdir(exist_ok=True)
        os.environ[variable] = str(target)
    return str(scratch)


def run(command, **kwargs):
    """Echo and run a command, failing on a non-zero exit."""
    print('RUN', *map(str, command), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def parse_stage(description):
    """The stage named on the command line (primary or additional)."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('stage', choices=STAGES)
    return parser.parse_args().stage


def genomes():
    """inputs/genomes.json: one entry per bearer genome."""
    return json.loads((WORK / 'inputs/genomes.json').read_text())


def task_genome():
    """genomes.json entry of this array task (SLURM_ARRAY_TASK_ID, else TELO_TASK_INDEX)."""
    index = int(os.environ.get('SLURM_ARRAY_TASK_ID', os.environ.get('TELO_TASK_INDEX')))
    return genomes()[index]


def orthology_dir(stage, genome_id, step=''):
    """One genome's aligner output, or a later step's: 'recall', 'corrected' or 'reviewed'."""
    return (WORK / step if step else WORK) / stage / 'orthology' / genome_id
