"""Align one genome's mapped loci against its comparator panel and classify the introns."""

import json
import os
import sys
from pathlib import Path

from classify_states import CALLER, classify
from common import REPO, WORK, dump, orthology_dir, parse_stage, read, run, sha, task_genome
from common import use_scratch, write

ALIGNER = REPO / 'lib/telotron_ortholog_align.py'


def main():
    use_scratch()
    stage = parse_stage(__doc__)
    meta = task_genome()
    genome_id = meta['genome_id']
    out_dir = orthology_dir(stage, genome_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    # the host mapping of this genome must be complete
    host_dir = WORK / stage / 'host' / genome_id
    assert json.loads((host_dir / 'complete.json').read_text())['complete']
    catalogue = host_dir / 'mapped_catalogue.tsv'
    rows = read(catalogue)

    # no mapped loci: an empty states table
    if not rows:
        write(
            out_dir / 'states.tsv', [], ['telotron', 'focal', 'ortholog', 'state', 'protein_ident']
        )
        dump(
            out_dir / 'complete.json',
            dict(
                complete=True,
                genome_id=genome_id,
                stage=stage,
                mapped_loci=0,
                comparators=len(meta['comparators']),
                states=0,
            ),
        )
        return

    # align the mapped loci against the comparator panel
    command = [sys.executable, ALIGNER, '--final', catalogue, '--focal-ids', genome_id]
    command += ['--ortholog-ids', ','.join(meta['comparators'])]
    command += ['--refseq-dir', REPO / 'data/raw/refseq', '--tara-dir', REPO / 'data/raw/tara']
    command += ['--outdir', out_dir, '--threads', os.environ.get('SLURM_CPUS_PER_TASK', '8')]
    search_path = f"{Path(__file__).resolve().parent}:{REPO / 'lib'}"
    run(command, cwd=REPO, env=dict(os.environ, PYTHONPATH=search_path))

    # classify states, then check the outputs cover exactly the mapped loci
    composition = classify(stage, meta)
    loci = {row['locus_id'] for row in rows}
    states = read(out_dir / 'states.tsv')
    assert {state['telotron'] for state in states} <= loci
    aligned = {row['telotron'] for row in read(out_dir / 'ortholog_map.tsv')}
    assert aligned == loci, (genome_id, len(aligned), len(loci))

    dump(
        out_dir / 'complete.json',
        dict(
            complete=True,
            genome_id=genome_id,
            stage=stage,
            mapped_loci=len(rows),
            comparators=len(meta['comparators']),
            states=len(states),
            states_sha256=sha(out_dir / 'states.tsv'),
            map_sha256=sha(out_dir / 'ortholog_map.tsv'),
            aligner_sha256=sha(ALIGNER),
            caller_sha256=sha(CALLER),
            scratch=use_scratch(),
            composition=composition,
        ),
    )
    print('ORTHOLOGY COMPLETE', stage, genome_id, len(states), flush=True)


if __name__ == '__main__':
    main()
