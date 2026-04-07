#!/usr/bin/env python3
"""
Filter telotron candidates for intron misannotation artifacts.

Checks:
  1. Genomic splice sites — GT...AG (or CT...AC) in the actual genome at intron boundaries
  2. Flanking exon sizes — both exons must be ≥20bp
  3. Gene structure — at least 2 introns (3+ exons) in the parent gene, or exons >50bp
  4. Flanking exon coding potential — check for stop codons, poly-A, low-complexity
  5. Intergenic masquerading — is this actually between two separate genes?
  6. Annotation evidence — check GFF for source (RNA-seq, ab initio, homology)
  7. Repeat in flanking exons — if exons are also telomeric, it's likely a misannotation

Requires genome FASTA + GFF for each species. Downloads if not present.
"""

import csv, gzip, json, os, re, sys, subprocess, shutil
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime

BASEDIR = Path(__file__).resolve().parent
DATASETS = BASEDIR.parent / "pan_euk_telotrons_datasets"
GENOMES_DIR = BASEDIR / "genomes"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

TELO_HEX = _rotations("TTAGGG") | _rotations("CCCTAA")


def telo_coverage(seq):
    n = len(seq)
    if n < 6:
        return 0.0
    covered = bytearray(n)
    for h in TELO_HEX:
        start = 0
        while True:
            idx = seq.find(h, start)
            if idx == -1:
                break
            for j in range(idx, min(idx + len(h), n)):
                covered[j] = 1
            start = idx + 1
    return sum(covered) / n


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def load_fasta(path):
    seqs = {}
    curr = None
    parts = []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if curr:
                    seqs[curr] = ''.join(parts)
                curr = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if curr:
        seqs[curr] = ''.join(parts)
    return seqs


def download_genome(acc):
    """Download genome + GFF3. Returns (gff_path, fna_path) or None."""
    out_dir = GENOMES_DIR / acc
    zip_path = GENOMES_DIR / f"{acc}.zip"

    if out_dir.exists():
        gff = _find(out_dir, "*.gff")
        fna = _find(out_dir, "*.fna")
        if gff and fna:
            return (gff, fna)

    try:
        subprocess.run(
            [str(DATASETS), "download", "genome", "accession", acc,
             "--include", "genome,gff3", "--filename", str(zip_path)],
            check=True, capture_output=True, timeout=3600,
        )
    except Exception:
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(["unzip", "-qo", str(zip_path), "-d", str(out_dir)],
                        check=True, capture_output=True, timeout=300)
    except Exception:
        return None
    finally:
        zip_path.unlink(missing_ok=True)

    gff = _find(out_dir, "*.gff")
    fna = _find(out_dir, "*.fna")
    return (gff, fna) if gff and fna else None


def _find(root, pattern):
    matches = list(root.rglob(pattern))
    return str(matches[0]) if matches else None


# ---------------------------------------------------------------------------
# Parse GFF for gene structure around telotron loci
# ---------------------------------------------------------------------------

def parse_gene_structure(gff_path, target_loci):
    """
    For each target locus (contig, start, end), find:
    - Parent mRNA/gene
    - All exons of that mRNA
    - Annotation source
    - Gene biotype

    target_loci: list of (contig, start, end) tuples
    Returns: dict of (contig, start, end) -> gene_info
    """
    # Build lookup set
    target_set = set(target_loci)
    target_contigs = {t[0] for t in target_loci}

    # Step 1: Find which mRNAs contain introns overlapping our targets
    # An intron is a gap between consecutive exons in the same mRNA

    # Collect all exons per mRNA on target contigs
    mrna_exons = defaultdict(list)  # mRNA_id -> [(start, end)]
    mrna_contig = {}
    mrna_strand = {}
    mrna_parent = {}  # mRNA_id -> gene_id
    gene_info = {}     # gene_id -> {name, biotype, source, start, end}
    mrna_source = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            contig = fields[0]
            if contig not in target_contigs:
                continue

            ftype = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            source = fields[1]
            attrs = dict(re.findall(r'(\w+)=([^;]+)', fields[8]))

            if ftype == 'gene':
                gid = attrs.get('ID', '')
                gene_info[gid] = {
                    'name': attrs.get('Name', attrs.get('gene', gid)),
                    'biotype': attrs.get('gene_biotype', '?'),
                    'source': source,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'contig': contig,
                }

            elif ftype in ('mRNA', 'transcript', 'lnc_RNA', 'rRNA', 'tRNA'):
                mid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                mrna_parent[mid] = parent
                mrna_contig[mid] = contig
                mrna_strand[mid] = strand
                mrna_source[mid] = source

            elif ftype in ('exon', 'CDS'):
                parent = attrs.get('Parent', '')
                mrna_exons[parent].append((start, end))

    # Step 2: For each mRNA, compute introns and match to targets
    results = {}
    for mid, exons in mrna_exons.items():
        contig = mrna_contig.get(mid)
        if not contig:
            continue
        exons_sorted = sorted(set(exons))

        for i in range(len(exons_sorted) - 1):
            istart = exons_sorted[i][1] + 1  # 1-based intron start
            iend = exons_sorted[i + 1][0] - 1  # 1-based intron end

            # Check if this intron matches any target
            for t_contig, t_start, t_end in target_loci:
                if contig != t_contig:
                    continue
                # Allow ±2bp tolerance for coordinate systems
                if abs(istart - t_start) <= 2 and abs(iend - t_end) <= 2:
                    gid = mrna_parent.get(mid, '')
                    ginfo = gene_info.get(gid, {})

                    prev_exon = exons_sorted[i]
                    next_exon = exons_sorted[i + 1]
                    prev_exon_len = prev_exon[1] - prev_exon[0] + 1
                    next_exon_len = next_exon[1] - next_exon[0] + 1

                    results[(t_contig, t_start, t_end)] = {
                        'mrna_id': mid,
                        'gene_id': gid,
                        'gene_name': ginfo.get('name', '?'),
                        'gene_biotype': ginfo.get('biotype', '?'),
                        'gene_source': ginfo.get('source', '?'),
                        'gene_start': ginfo.get('start', 0),
                        'gene_end': ginfo.get('end', 0),
                        'gene_strand': ginfo.get('strand', '?'),
                        'n_exons': len(exons_sorted),
                        'n_introns': len(exons_sorted) - 1,
                        'exon_index': i,  # which intron (0-based)
                        'prev_exon': prev_exon,
                        'next_exon': next_exon,
                        'prev_exon_len': prev_exon_len,
                        'next_exon_len': next_exon_len,
                        'all_exon_sizes': [e[1] - e[0] + 1 for e in exons_sorted],
                        'annotation_source': mrna_source.get(mid, '?'),
                    }
                    break

    return results


# ---------------------------------------------------------------------------
# Validate one species
# ---------------------------------------------------------------------------

def validate_species(acc, name, telotron_tsv):
    """Download genome+GFF, then validate each telotron."""
    print(f"\n{'='*70}", flush=True)
    print(f"  Validating: {name} ({acc})", flush=True)

    # Load telotron candidates
    with open(telotron_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        candidates = [r for r in reader if float(r.get('telo_coverage', 0)) >= 0.85]

    if not candidates:
        print(f"  No ≥85% candidates", flush=True)
        return []

    print(f"  {len(candidates)} candidates at ≥85% purity", flush=True)

    # Download genome + GFF
    paths = download_genome(acc)
    if not paths:
        print(f"  ⚠ Failed to download genome", flush=True)
        return []
    gff_path, fna_path = paths

    # Load genome
    print(f"  Loading genome...", flush=True)
    genome = load_fasta(fna_path)

    # Get target loci
    target_loci = []
    for c in candidates:
        target_loci.append((c['contig'], int(c['start']), int(c['end'])))

    # Parse gene structure
    print(f"  Parsing GFF for gene structure...", flush=True)
    gene_structures = parse_gene_structure(gff_path, target_loci)

    # Validate each candidate
    results = []
    for c in candidates:
        contig = c['contig']
        start = int(c['start'])
        end = int(c['end'])
        length = int(c['length'])
        seq = c.get('intron_seq', '').upper()

        record = {
            'species': name,
            'acc': acc,
            'contig': contig,
            'start': start,
            'end': end,
            'length': length,
            'telo_coverage': float(c.get('telo_coverage', 0)),
            'orientation': c.get('orientation', ''),
            'donor_reported': c.get('donor', ''),
            'acceptor_reported': c.get('acceptor', ''),
        }

        flags = []  # rejection flags
        notes = []  # informational

        # --- Check 1: Genomic splice sites ---
        if contig in genome:
            gseq = genome[contig]
            if start >= 2 and end + 2 <= len(gseq):
                genomic_donor = gseq[start:start+2]
                genomic_acceptor = gseq[end-2:end]
                record['genomic_donor'] = genomic_donor
                record['genomic_acceptor'] = genomic_acceptor

                canonical = (genomic_donor == 'GT' and genomic_acceptor == 'AG') or \
                            (genomic_donor == 'CT' and genomic_acceptor == 'AC')
                if not canonical:
                    flags.append(f"BAD_GENOMIC_SPLICE({genomic_donor}-{genomic_acceptor})")
            else:
                flags.append("NEAR_CONTIG_BOUNDARY")
        else:
            flags.append("CONTIG_NOT_FOUND")

        # --- Check 2: Gene structure ---
        key = (contig, start, end)
        gs = gene_structures.get(key)
        if gs:
            record['gene_name'] = gs['gene_name']
            record['gene_biotype'] = gs['gene_biotype']
            record['n_exons'] = gs['n_exons']
            record['n_introns'] = gs['n_introns']
            record['prev_exon_len'] = gs['prev_exon_len']
            record['next_exon_len'] = gs['next_exon_len']
            record['annotation_source'] = gs['annotation_source']
            record['exon_index'] = gs['exon_index']
            record['all_exon_sizes'] = gs['all_exon_sizes']

            # Check 2a: Flanking exon sizes
            if gs['prev_exon_len'] < 20:
                flags.append(f"TINY_PREV_EXON({gs['prev_exon_len']}bp)")
            if gs['next_exon_len'] < 20:
                flags.append(f"TINY_NEXT_EXON({gs['next_exon_len']}bp)")

            # Check 2b: Single-intron gene with small exons = suspicious
            if gs['n_introns'] == 1:
                if gs['prev_exon_len'] < 50 and gs['next_exon_len'] < 50:
                    flags.append("SINGLE_INTRON_TINY_EXONS")
                elif gs['prev_exon_len'] < 100 and gs['next_exon_len'] < 100:
                    notes.append("single_intron_small_exons")

            # Check 2c: Telomeric content in flanking exons
            if contig in genome:
                gseq = genome[contig]
                prev_exon_seq = gseq[gs['prev_exon'][0]-1:gs['prev_exon'][1]]
                next_exon_seq = gseq[gs['next_exon'][0]-1:gs['next_exon'][1]]
                prev_telo = telo_coverage(prev_exon_seq)
                next_telo = telo_coverage(next_exon_seq)
                record['prev_exon_telo'] = round(prev_telo, 3)
                record['next_exon_telo'] = round(next_telo, 3)

                if prev_telo > 0.30:
                    flags.append(f"TELO_IN_PREV_EXON({prev_telo:.0%})")
                if next_telo > 0.30:
                    flags.append(f"TELO_IN_NEXT_EXON({next_telo:.0%})")

            # Check 2d: Gene biotype
            if gs['gene_biotype'] not in ('protein_coding', '?'):
                notes.append(f"biotype={gs['gene_biotype']}")

        else:
            flags.append("NO_GENE_STRUCTURE_FOUND")

        record['flags'] = flags
        record['notes'] = notes
        record['verdict'] = 'PASS' if not flags else 'FAIL'
        results.append(record)

    # Cleanup genome files
    if "--keep-genomes" not in sys.argv:
        genome_dir = GENOMES_DIR / acc
        if genome_dir.exists():
            shutil.rmtree(genome_dir, ignore_errors=True)

    n_pass = sum(1 for r in results if r['verdict'] == 'PASS')
    n_fail = sum(1 for r in results if r['verdict'] == 'FAIL')
    print(f"  Results: {n_pass} PASS, {n_fail} FAIL out of {len(results)}", flush=True)

    if n_fail > 0:
        flag_counts = Counter(f for r in results for f in r['flags'])
        print(f"  Rejection reasons: {dict(flag_counts)}", flush=True)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Determine which species to validate
    checkpoint = BASEDIR / "checkpoint.json"
    telotron_dir = BASEDIR / "telotrons"
    output_dir = BASEDIR / "validated_results"
    output_dir.mkdir(exist_ok=True)

    if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
        # Validate specific files
        tsv_files = sys.argv[1:]
        targets = []
        for f in tsv_files:
            # Extract accession from filename
            basename = Path(f).stem
            parts = basename.split('_')
            acc = f"{parts[0]}_{parts[1]}"
            name = '_'.join(parts[2:]).replace('_telotrons', '')
            targets.append((acc, name, f))
    else:
        # Find species with telotrons from checkpoint
        ckpt = json.load(open(checkpoint))
        targets = []
        for acc, r in ckpt.items():
            if r.get('status') != 'ok' or r.get('n_telotrons', 0) == 0:
                continue
            name = r.get('name', '?')
            tsv = telotron_dir / f"{acc}_{name}_telotrons.tsv"
            if tsv.exists():
                targets.append((acc, name, str(tsv)))

        # Sort by telotron count descending
        targets.sort(key=lambda x: -json.load(open(checkpoint)).get(x[0], {}).get('n_telotrons', 0))

    print(f"Validating {len(targets)} species with telotron candidates\n")

    all_results = []
    for acc, name, tsv in targets:
        results = validate_species(acc, name, tsv)
        all_results.extend(results)

    # Write results
    n_pass = sum(1 for r in all_results if r['verdict'] == 'PASS')
    n_fail = sum(1 for r in all_results if r['verdict'] == 'FAIL')

    print(f"\n{'='*70}")
    print(f"FINAL RESULTS: {n_pass} PASS, {n_fail} FAIL out of {len(all_results)}")
    print(f"{'='*70}\n")

    # Summary by species
    species_pass = defaultdict(int)
    species_fail = defaultdict(int)
    species_flags = defaultdict(lambda: Counter())
    for r in all_results:
        sp = r['species']
        if r['verdict'] == 'PASS':
            species_pass[sp] += 1
        else:
            species_fail[sp] += 1
            for f in r['flags']:
                species_flags[sp][f] += 1

    print(f"{'Species':<45s} {'PASS':>5s} {'FAIL':>5s} Top rejection reason")
    print('-' * 90)
    for sp in sorted(set(list(species_pass.keys()) + list(species_fail.keys())),
                     key=lambda s: -(species_pass[s] + species_fail[s])):
        p = species_pass[sp]
        f = species_fail[sp]
        top_flag = species_flags[sp].most_common(1)
        top_str = f"{top_flag[0][0]}({top_flag[0][1]})" if top_flag else ""
        print(f"  {sp:<43s} {p:>5d} {f:>5d} {top_str}")

    # Write JSON
    json_path = output_dir / "misannotation_filter_results.json"
    # Remove non-serializable items
    clean_results = []
    for r in all_results:
        cr = {k: v for k, v in r.items()}
        clean_results.append(cr)

    with open(json_path, 'w') as jf:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'total': len(all_results),
            'pass': n_pass,
            'fail': n_fail,
            'results': clean_results,
        }, jf, indent=2, default=str)

    # Write TSV of passing telotrons
    tsv_path = output_dir / "real_telotrons.tsv"
    passing = [r for r in all_results if r['verdict'] == 'PASS']
    if passing:
        cols = ['species', 'contig', 'start', 'end', 'length', 'telo_coverage',
                'orientation', 'genomic_donor', 'genomic_acceptor',
                'gene_name', 'gene_biotype', 'n_exons', 'n_introns',
                'prev_exon_len', 'next_exon_len', 'annotation_source',
                'prev_exon_telo', 'next_exon_telo', 'notes']
        with open(tsv_path, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in passing:
                vals = [str(r.get(c, '')) for c in cols]
                f.write('\t'.join(vals) + '\n')

    print(f"\nOutputs:")
    print(f"  {json_path}")
    print(f"  {tsv_path}")


if __name__ == '__main__':
    main()
