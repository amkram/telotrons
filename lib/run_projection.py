"""miniprot projection of donor proteomes onto one target genome, with marker-based donor choice."""

import gzip, shutil, subprocess
from pathlib import Path
from audit_its_opportunity import open_text

MINIPROT_VERSION = '0.18-r281'
# marker search: near-complete best hits only
MARKER_OPTIONS = ['--outs=0.99', '--outc=0.5']
# full-proteome projection: up to 30 secondary alignments per protein
PROJECTION_OPTIONS = ['--outs=0.7', '--outc=0.2', '-N', '30', '-p', '0.7']
# a marker hit counts with >= 50% query coverage and >= 35% identity
MARKER_COVERAGE = 0.5
MARKER_IDENTITY = 0.35
MARKERS_PER_DONOR = 80


def tags(text):
    """GFF column 9 as a dict."""
    return dict(item.split('=', 1) for item in text.split(';') if '=' in item)


def alignments(path):
    """miniprot --gff alignments: PAF fields, mRNA coordinates and CDS blocks."""
    # current: alignment filled from GFF rows, pending: PAF record awaiting its mRNA row
    current = None
    pending = None
    with open_text(path) as handle:
        for line in handle:
            # each ##PAF line starts a new alignment and closes the previous one
            if line.startswith('##PAF\t'):
                if current is not None:
                    yield current
                    current = None
                fields = line.rstrip('\n').split('\t')[1:]
                # unmapped query
                if fields[5] == '*':
                    pending = None
                    continue

                # TAG:TYPE:VALUE fields after the 12 PAF columns
                paf_tags = {}
                for tag in fields[12:]:
                    name, _, value = tag.split(':', 2)
                    paf_tags[name] = value
                pending = dict(
                    query=fields[0],
                    qlen=int(fields[1]),
                    qstart=int(fields[2]),
                    qend=int(fields[3]),
                    strand=fields[4],
                    seqid=fields[5],
                    tlen=int(fields[6]),
                    pstart=int(fields[7]),
                    pend=int(fields[8]),
                    score=int(paf_tags['AS']),
                )
            elif line.startswith('#'):
                continue
            else:
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9:
                    continue
                attributes = tags(fields[8])

                # the mRNA row completes the pending PAF record
                if fields[2] == 'mRNA':
                    assert current is None and pending is not None
                    current = pending
                    pending = None
                    current.update(
                        start0=int(fields[3]) - 1,
                        end0=int(fields[4]),
                        identity=float(attributes['Identity']),
                        id=attributes['ID'],
                        cds=[],
                    )
                    # the GFF row must agree with its PAF record
                    assert attributes['Target'].split()[0] == current['query']
                    assert fields[0] == current['seqid'] and fields[6] == current['strand']
                    assert 0 <= current['start0'] < current['end0'] <= current['tlen']
                    # GFF mRNA may include the terminal stop (<=3 bp) that PAF omits
                    assert 0 <= current['pstart'] - current['start0'] <= 3
                    assert 0 <= current['end0'] - current['pend'] <= 3

                # CDS rows belong to the current mRNA and lie inside it
                elif fields[2] == 'CDS':
                    assert current is not None and attributes['Parent'] == current['id']
                    start, end = int(fields[3]) - 1, int(fields[4])
                    assert current['start0'] <= start < end <= current['end0']
                    current['cds'].append((start, end))

    # flush the last alignment, every PAF record having had its mRNA
    if current is not None:
        yield current
    assert pending is None


def check_version():
    """Stop unless miniprot is the pinned version."""
    version = subprocess.check_output(['miniprot', '--version'], text=True).strip()
    assert version == MINIPROT_VERSION, version


def index_genome(fasta_path, index, threads, log):
    """Build the miniprot index of the target genome."""
    with open(log, 'w') as handle:
        command = ['miniprot', '-t', str(threads), '-d', str(index), str(fasta_path)]
        subprocess.run(command, stdout=handle, stderr=handle, check=True)


def align(index, query, options, threads, output, log):
    """Align query proteins to the indexed genome, streaming the GFF into a gzip file."""
    command = ['miniprot', '-t', str(threads), '-I', '-j', '1', '--gff', *options]
    command += [str(index), str(query)]
    with open(log, 'w') as errors, gzip.open(output, 'wb', compresslevel=3) as out:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        try:
            shutil.copyfileobj(process.stdout, out, 1 << 20)
        finally:
            process.stdout.close()
        assert process.wait() == 0, (command, log)


def pool(paths, output):
    """Concatenate protein FASTA files."""
    with open(output, 'wb') as out:
        for path in paths:
            with open(path, 'rb') as handle:
                shutil.copyfileobj(handle, out)


def select_donors(candidates, proteomes, index, folder, scratch, threads, min_markers):
    """Score candidate donors by marker identity x coverage and select the top three."""
    marker_query = Path(scratch) / 'markers.faa'
    pool([Path(proteomes) / f'{donor}.markers.faa' for donor in candidates], marker_query)
    markers = Path(folder) / 'markers.gff.gz'
    align(index, marker_query, MARKER_OPTIONS, threads, markers, Path(folder) / 'markers.log')

    # best identity x coverage per marker, among hits that pass both cutoffs
    best = {}
    for hit in alignments(markers):
        coverage = (hit['qend'] - hit['qstart']) / hit['qlen']
        if coverage >= MARKER_COVERAGE and hit['identity'] >= MARKER_IDENTITY:
            best[hit['query']] = max(best.get(hit['query'], 0), hit['identity'] * coverage)

    # donor score: mean over its 80 markers (unhit ones score 0), eligible with >= min_markers hits
    scores = []
    for donor in candidates:
        values = [score for query, score in best.items() if query.startswith(donor + '|')]
        scores.append(
            dict(
                donor=donor,
                n_markers_hit=len(values),
                score=sum(values) / MARKERS_PER_DONOR,
                eligible=len(values) >= min_markers,
            )
        )

    # top three eligible donors, ties broken by donor ID, the other eligible ones additional
    scores.sort(key=lambda row: (-row['score'], row['donor']))
    eligible = [row['donor'] for row in scores if row['eligible']]
    # no eligible donor is a failure, never a silently passed mask
    assert eligible, ('No eligible related donors', scores)
    return dict(candidates=scores, selected=eligible[:3], additional=eligible[3:])


def project(donors, proteomes, index, folder, name, scratch, threads):
    """Project the donors' full proteomes into folder/<name>.gff.gz and return its path."""
    query = Path(scratch) / 'donors.faa'
    pool([Path(proteomes) / f'{donor}.faa' for donor in donors], query)
    output = Path(folder) / f'{name}.gff.gz'
    align(index, query, PROJECTION_OPTIONS, threads, output, Path(folder) / f'{name}.log')

    # every projected alignment must carry CDS rows
    for alignment in alignments(output):
        assert alignment['cds'], alignment['id']
    return output
