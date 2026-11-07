"""Per-genome ITS scan: components outside the whole-model gene mask and the terminal caps."""

import contextlib, csv, gzip, json, time
from pathlib import Path
from audit_its_opportunity import count_acgt, fasta, gene_mask, overlaps, subtract
from its_rescan import cap_mask, components, scan, validate_tracts


@contextlib.contextmanager
def table(path, header):
    """A gzipped TSV writer with its header row."""
    with gzip.open(path, 'wt', compresslevel=5) as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        yield writer


def run_genome(genome_id, motif, fasta_path, gff_path, folder, core):
    """Scan one genome into folder: <genome>.components/.caps.tsv.gz and <genome>.json."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    last_report = started
    genic = gene_mask(gff_path)
    seen = set()
    scanned_bp = 0
    n_components = 0
    eligible_acgt = 0

    with (
        table(folder / f'{genome_id}.components.tsv.gz', ['seqid', 'start0', 'end0']) as kept,
        table(folder / f'{genome_id}.caps.tsv.gz', ['seqid', 'start0', 'end0']) as caps,
    ):
        for contig, sequence in fasta(fasta_path):
            assert contig not in seen, contig
            seen.add(contig)
            scanned_bp += len(sequence)
            genes = genic.get(contig, [])
            assert all(0 <= start < end <= len(sequence) for start, end in genes)

            # tracts, each rechecked, and the terminal caps
            tracts = scan(sequence, motif, core)
            validate_tracts(sequence, motif, tracts)
            cap = cap_mask(sequence, motif)
            caps.writerows([contig, start, end] for start, end in cap)

            # eligible A/C/G/T bp: outside the gene mask and the cap
            masked = genes + subtract(cap, genes)
            eligible_acgt += count_acgt(sequence) - sum(
                count_acgt(sequence[start:end]) for start, end in masked
            )

            # components that touch neither mask
            for start, end in components(tracts):
                if not overlaps(start, end, genes) and not overlaps(start, end, cap):
                    kept.writerow([contig, start, end])
                    n_components += 1

            # progress line at most every 30 s
            if time.monotonic() - last_report > 30:
                print(f'  {genome_id}: {scanned_bp / 1e6:.0f} Mb scanned', flush=True)
                last_report = time.monotonic()

    # no masked contig may be missing from the FASTA
    assert not set(genic) - seen, sorted(set(genic) - seen)[:10]
    assert eligible_acgt > 0

    result = dict(
        genome_id=genome_id,
        motif=motif,
        n_components=n_components,
        eligible_acgt_bp=eligible_acgt,
        components_per_acgt_mb=n_components * 1e6 / eligible_acgt,
    )
    (folder / f'{genome_id}.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f'DONE {genome_id}: {time.monotonic() - started:.1f}s', flush=True)
    return result
