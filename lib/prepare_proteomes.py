"""Nonredundant donor proteomes and protein markers derived from one genome's annotation."""

import collections, hashlib, math
from pathlib import Path
from audit_its_opportunity import fasta, gff_features


def entropy(protein):
    """Shannon entropy, in bits, of the residue composition."""
    total = len(protein)
    return -sum(
        (count / total) * math.log2(count / total)
        for count in collections.Counter(protein).values()
    )


def read_cds_models(gff_path):
    """CDS segments and translation tables per (contig, strand, transcript), and parent links."""
    models = {}
    # parents of each non-CDS feature, used to climb from transcript to gene
    parents = {}

    for row in gff_features(gff_path):
        if row.feature_id and row.parents and row.kind != 'CDS':
            parents[row.contig, row.strand, row.feature_id] = row.parents

        # keep stranded CDS rows that name a parent transcript
        if row.kind != 'CDS' or not row.parents or row.strand not in '+-':
            continue

        # (start, end, phase) per transcript, an unknown phase counted as 0
        phase = int(row.phase) if row.phase in '012' else 0
        for transcript in row.parents:
            model = models.setdefault(
                (row.contig, row.strand, transcript), dict(segments=set(), tables=set())
            )
            model['segments'].add((row.start, row.end, phase))
            model['tables'].add(int(row.attributes.get('transl_table', '1')))
    return models, parents


def translate(dna, strand, segments, table):
    """Protein of one CDS model: spliced in transcript order, trailing stops removed."""
    from Bio.Seq import Seq

    # minus-strand segments are reverse-complemented and joined in reverse order
    if strand == '-':
        segments = segments[::-1]
    parts = [
        dna[start:end] if strand == '+' else str(Seq(dna[start:end]).reverse_complement())
        for start, end, phase in segments
    ]
    # skip the first segment's phase offset and translate whole codons only
    coding = ''.join(parts)[segments[0][2] :]
    return str(Seq(coding[: len(coding) // 3 * 3]).translate(table=table)).rstrip('*')


def gene_of(key, parents):
    """Climb single-parent links from a transcript to its gene."""
    contig, strand, gene = key
    visited = set()
    while (contig, strand, gene) in parents and gene not in visited:
        visited.add(gene)
        above = parents[contig, strand, gene]
        if len(above) != 1:
            break
        gene = above[0]
    return gene


def derive(fasta_path, gff_path):
    """Clean (contig, strand, gene, transcript, protein) rows: longest per gene, no repeated sequence."""
    models, parents = read_cds_models(gff_path)
    by_contig = collections.defaultdict(list)
    for key in models:
        by_contig[key[0]].append(key)

    seen = set()
    proteins = []
    for contig, dna in fasta(fasta_path):
        assert contig not in seen, contig
        seen.add(contig)
        for key in by_contig.get(contig, []):
            model = models[key]
            segments = sorted(model['segments'])

            # reject overlapping segments and mixed translation tables
            if any(left[1] > right[0] for left, right in zip(segments, segments[1:])):
                continue
            assert all(0 <= start < end <= len(dna) for start, end, phase in segments)
            if len(model['tables']) != 1:
                continue
            protein = translate(dna, key[1], segments, next(iter(model['tables'])))

            # reject short, stop-containing, ambiguous and low-complexity proteins
            if len(protein) < 50 or '*' in protein:
                continue
            if protein.count('X') / len(protein) > 0.1 or entropy(protein) < 2.5:
                continue
            proteins.append((key[0], key[1], gene_of(key, parents), key[2], protein))

    # every annotated contig must be present in the FASTA
    assert not set(by_contig) - seen

    # longest clean isoform per gene, then one entry per exact sequence
    longest = {}
    for row in sorted(proteins, key=lambda row: (-len(row[4]), row[:4])):
        longest.setdefault(row[:3], row)
    unique = {}
    for row in sorted(longest.values(), key=lambda row: row[:4]):
        unique.setdefault(row[4], row)

    return sorted(unique.values(), key=lambda row: row[:4])


def write_fasta(path, named):
    """Write (query, ..., protein) rows as FASTA."""
    with open(path, 'w') as handle:
        for query, *rest in named:
            handle.write(f'>{query}\n{rest[-1]}\n')


def write_proteome(genome_id, fasta_path, gff_path, folder):
    """<genome>.faa (every clean protein) and <genome>.markers.faa (80 marker proteins)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rows = derive(fasta_path, gff_path)
    assert rows, genome_id
    # query IDs: genome ID plus a 1-based protein number
    named = [(f'{genome_id}|p{i:06d}', *row) for i, row in enumerate(rows, 1)]
    write_fasta(folder / f'{genome_id}.faa', named)

    # deterministic marker sample, blind to repeats, role and alignment outcome:
    # the first 80 mid-length, complex, unambiguous proteins in sequence-hash order
    candidates = [
        row
        for row in named
        if 200 <= len(row[-1]) <= 1000 and entropy(row[-1]) >= 3.5 and 'X' not in row[-1]
    ]
    markers = sorted(candidates, key=lambda row: hashlib.sha256(row[-1].encode()).hexdigest())[:80]
    assert len(markers) == 80, (genome_id, len(markers))
    write_fasta(folder / f'{genome_id}.markers.faa', markers)
    print(genome_id, len(named), 'proteins', flush=True)
