"""Skip a finished alignment when its inputs and output are unchanged."""

import hashlib
import json
import subprocess
from pathlib import Path

from common import dump, sha


def cached_msa_job(align):
    """Wrap an alignment job (path, sequences, is_protein) with an .inputs.json sidecar check."""
    version = subprocess.run(
        ['mafft', '--version'], capture_output=True, text=True, check=True
    ).stderr.strip()

    def cached(job):
        path, seqs, is_protein = job
        output = Path(path)
        sidecar = Path(str(path) + '.inputs.json')

        # inputs digest: sequences, alphabet flag and MAFFT version
        payload = json.dumps([seqs, is_protein, version], separators=(',', ':'))
        digest = hashlib.sha256(payload.encode()).hexdigest()
        if output.exists() and sidecar.exists():
            prior = json.loads(sidecar.read_text())
            if prior['input_sha256'] == digest and prior['output_sha256'] == sha(output):
                return
        align(job)
        dump(sidecar, dict(input_sha256=digest, output_sha256=sha(output), mafft_version=version))

    return cached
