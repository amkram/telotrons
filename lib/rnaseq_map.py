"""Download SRA/ENA runs and map them with minimap2 into one sorted, indexed BAM."""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from _common import scratch_root
from rna_evidence import atomic_path, bam_header, checked_output


def log(*args):
    """Progress message on stderr."""
    print('[rnaseq_map]', *args, file=sys.stderr, flush=True)


def run(command):
    """Log and run a command. A nonzero exit raises."""
    log('$', *command)
    subprocess.run(command, check=True)


def check_preset_version(version, preset):
    """Raise unless the minimap2 version supports the preset (splice:sr needs >= 2.29)."""
    # major.minor at the start of the version string
    match = re.match(r'(\d+)\.(\d+)', version)
    if not match:
        raise ValueError(f'Unrecognized minimap2 version: {version}')
    if preset == 'splice:sr' and tuple(map(int, match.groups())) < (2, 29):
        raise ValueError('splice:sr requires minimap2 >=2.29; do not substitute a long-read preset')


def read_batches(scratch, accession):
    """FASTQ batches of one run: the mate pair, the unpaired reads, or both."""
    root = Path(scratch)
    # fasterq-dump --split-3 names: _1/_2 for mates, no suffix for unpaired reads
    first, second, single = [root / f'{accession}{suffix}.fastq' for suffix in ('_1', '_2', '')]
    present = {path: path.exists() and path.stat().st_size > 0 for path in (first, second, single)}
    if present[first] != present[second]:
        raise ValueError(f'Incomplete read pair for {accession}')

    batches = []
    if present[first]:
        batches.append([str(first), str(second)])
    if present[single]:
        batches.append([str(single)])
    if not batches:
        raise ValueError(f'No nonempty FASTQs for {accession}')
    return batches


def map_and_sort(command, output, logfile):
    """Pipe the mapper into samtools sort, writing output atomically."""
    log('$', *command)
    with atomic_path(output) as temporary:
        with open(logfile, 'w') as mapper_log:
            mapper = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=mapper_log)
            # stop the mapper if the sorter cannot start
            try:
                sorter = subprocess.Popen(
                    ['samtools', 'sort', '-@', '4', '-m', '2G', '-o', temporary, '-'],
                    stdin=mapper.stdout,
                    stderr=subprocess.PIPE,
                )
            except BaseException:
                mapper.stdout.close()
                mapper.terminate()
                mapper.wait()
                raise
            # close our pipe copy so the mapper sees the sorter exit
            mapper.stdout.close()

            # wait for both, and on interruption stop whichever is still running
            try:
                _, sort_error = sorter.communicate()
                mapper_status = mapper.wait()
            except BaseException:
                for process in (mapper, sorter):
                    if process.poll() is None:
                        process.terminate()
                mapper.wait()
                sorter.wait()
                raise

            if mapper_status or sorter.returncode:
                raise RuntimeError(
                    f'Mapping pipeline failed: mapper={mapper_status}, sorter={sorter.returncode}; '
                    f'see {logfile}\n{sort_error.decode(errors="replace")}'
                )

        checked_output(['samtools', 'quickcheck', '-v', temporary])


def map_runs(genome_fasta, accessions, preset, output, sample, threads=12, max_size='120g'):
    """Map every run to the genome into output (BAM), output.bai and output.logs/."""
    # distinct, well-formed run accessions
    if (
        not accessions
        or len(set(accessions)) != len(accessions)
        or any(not re.fullmatch(r'[SED]RR\d+', accession) for accession in accessions)
    ):
        raise ValueError(f'Supply distinct SRR/ERR/DRR accessions: {accessions}')
    assert preset in ('splice:sr', 'splice', 'splice:hq'), preset
    check_preset_version(checked_output(['minimap2', '--version']).strip(), preset)
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    # fresh scratch per call, kept on failure for inspection
    scratch = Path(tempfile.mkdtemp(prefix='rnaseq_', dir=scratch_root()))
    log(f'Scratch retained on failure: {scratch}')
    batch_bams = []
    log_paths = []

    # download each run, then map and sort each read batch on its own
    for accession in accessions:
        run(['prefetch', '-O', str(scratch), '--max-size', max_size, accession])
        run(
            [
                'fasterq-dump',
                '--split-3',
                '-e',
                str(threads),
                '--temp',
                str(scratch / 'tmp'),
                '-O',
                str(scratch),
                str(scratch / accession / (accession + '.sra')),
            ]
        )
        for batch_index, reads in enumerate(read_batches(scratch, accession)):
            bam = scratch / f'{accession}.{batch_index}.bam'
            logfile = scratch / f'{accession}.{batch_index}.minimap2.log'
            read_group = f'@RG\\tID:{accession}\\tSM:{sample}'
            command = ['minimap2', '-ax', preset, '-ub', '--secondary=no', '-R', read_group]
            command += ['-t', str(threads), str(genome_fasta), *reads]
            map_and_sort(command, str(bam), str(logfile))
            batch_bams.append(str(bam))
            log_paths.append(logfile)

    # a single batch is copied, several are merged, and the BAM is indexed before publishing
    with atomic_path(output) as temporary:
        index = temporary + '.bai'
        try:
            if len(batch_bams) == 1:
                shutil.copyfile(batch_bams[0], temporary)
            else:
                run(['samtools', 'merge', '-f', '-@', str(threads), temporary, *batch_bams])
            run(['samtools', 'index', temporary])
            bam_header(temporary)
            Path(index).replace(str(output) + '.bai')
        finally:
            if Path(index).exists():
                Path(index).unlink()

    # keep the mapper logs beside the BAM, removing scratch only after success
    logdir = Path(str(output) + '.logs')
    logdir.mkdir(exist_ok=True)
    for logfile in log_paths:
        shutil.copyfile(logfile, logdir / logfile.name)
    shutil.rmtree(scratch)
    log(f'Wrote checked BAM and index: {output}')
