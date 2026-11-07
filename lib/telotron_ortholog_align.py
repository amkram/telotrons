#!/usr/bin/env python3
"""Map telotron host proteins to comparator orthologs and call the homologous intron FILL/ABSENT."""

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

from Bio.Align import PairwiseAligner, substitution_matrices

from _common import ensure_tool_on_path, find_genome_files, rc, read_fasta, scratch_root
from hostlib import load_cds
from indexed_reference import Reference, miniprot_reference
from msa_cache import cached_msa_job

MIN_IDENT = 0.40  # miniprot identity floor (protein twilight zone)
MSA_MIN_IDENT = 0.0  # MSAs keep sub-floor hits (tagged |lowid), calls stay gated on MIN_IDENT
INTRON_TOL = 8  # residues of slack when matching an intron position
ABSENT_MARGIN = 15  # aligned exon residues needed on each side to call ABSENT
FLANK_AA = 20  # protein flank compared on each side (residues)
MAX_FLANK_MISMATCH = 5  # per flank, above this a FILL/ABSENT becomes POOR_FLANK
MIN_INTRONS_ELSEWHERE = 1  # retrocopy guard
ABSENT_WINDOW = 2000  # bp either side of the site scanned for N-gaps
MAX_N_FRAC = 0.05  # max N fraction around an ABSENT site (N-gaps mimic "no intron")
EDGE_BP = 1000  # minimum distance from either contig end (bp)
MAX_ALT_LOCI = 0  # near-best alternative loci tolerated (paralog guard)
ALT_LOCUS_COVERAGE = 0.9  # a hit covering this fraction of the best is an alternative locus
MOTIF = "TTTAGGG"  # fallback telomere motif for rows without one
MIN_MSA_INTRON_BP = 10  # shorter comparator introns stay out of the intron MSA

CODON = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',
    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',
    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'TAT': 'Y',
    'TAC': 'Y',
    'TAA': '*',
    'TAG': '*',
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'TGT': 'C',
    'TGC': 'C',
    'TGA': '*',
    'TGG': 'W',
    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',
}

# one global BLOSUM62 aligner shared by all flank comparisons
ALIGNER = PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -11
ALIGNER.extend_gap_score = -1
ALIGNER.mode = "global"

MAP_HEADER = (
    "telotron\tfocal\tgene\tprotein\tintron_residue\tortholog\tpct_identity\t"
    "call\tn_frameshift\tn_stop\tframeshift_at_junction\t"
    "ortho_intron_len\tortho_intron_telofrac\tlocal_ident\t"
    "n_alternative_loci\tflank5_mm\tflank3_mm\trbh_pass\tintron_phase\n"
)


def translate(seq):
    """Frame-0 translation, unknown codons as X."""
    seq = seq.upper()
    return "".join(CODON.get(seq[i : i + 3], 'X') for i in range(0, len(seq) - 2, 3))


def ordered_segments(model):
    """CDS segments in transcript order (descending coordinates on the minus strand)."""
    return sorted(model["segs"], key=lambda seg: -seg[0] if model["strand"] == "-" else seg[0])


def first_phase(model):
    """Phase of the model's first CDS segment in transcript order."""
    segments = ordered_segments(model)
    return model.get("phase", {}).get(segments[0], 0) if segments else 0


def protein_of(model, contig_seq):
    """Translation of a CDS model, starting on its first whole codon."""
    parts = []
    for start, end in ordered_segments(model):
        exon = contig_seq[start - 1 : end]
        parts.append(rc(exon) if model["strand"] == "-" else exon)

    # skip the first segment's phase bases so translation starts on a codon
    return translate("".join(parts)[first_phase(model) :]).rstrip("*")


def intron_residue(model, intron_start, intron_end):
    """(residue, phase) of the CDS gap equal to the 1-based inclusive intron, else None."""
    segments = ordered_segments(model)
    minus = model["strand"] == "-"

    # coding bases upstream of each gap, less the start phase
    coding_bp = -first_phase(model)
    for upstream, downstream in zip(segments, segments[1:]):
        if minus:
            gap = (downstream[1] + 1, upstream[0] - 1)
        else:
            gap = (upstream[1] + 1, downstream[0] - 1)
        coding_bp += upstream[1] - upstream[0] + 1
        if gap == (intron_start, intron_end):
            return divmod(coding_bp, 3)
    return None


def run_miniprot(genome_fa, query_faa, threads, out_gff):
    """miniprot GFF of query proteins against a genome, published only when complete."""
    # write via a partial file so a failed run is never reused as cache
    partial = f"{out_gff}.partial.{os.getpid()}"
    command = ["miniprot", "-t", str(threads), "--gff"]
    command += [miniprot_reference(genome_fa, threads), query_faa]
    try:
        with open(partial, "w") as gff, open(out_gff + ".log", "w") as log:
            subprocess.run(command, stdout=gff, stderr=log, check=True)
        os.replace(partial, out_gff)
    finally:
        if os.path.exists(partial):
            os.unlink(partial)


def read_miniprot_gff(gff):
    """mRNA records by ID, and their CDS segments (genome and query coordinates)."""
    mrna = {}
    cds = {}
    for line in open(gff):
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]
        target = re.search(r"Target=(\S+) (\d+) (\d+)", attributes)

        if fields[2] == "mRNA":
            mrna_id = re.search(r"ID=([^;]+)", attributes).group(1)
            mrna[mrna_id] = dict(
                prot=target.group(1),
                rank=int(re.search(r"Rank=(\d+)", attributes).group(1)),
                ident=float(re.search(r"Identity=([0-9.]+)", attributes).group(1)),
                contig=fields[0],
                strand=fields[6],
                qmin=int(target.group(2)),
                qmax=int(target.group(3)),
            )
            cds[mrna_id] = []
        elif fields[2] == "CDS":
            parent = re.search(r"Parent=([^;]+)", attributes).group(1)
            frameshifts = re.search(r"Frameshift=(\d+)", attributes)
            stops = re.search(r"StopCodon=(\d+)", attributes)
            identity = re.search(r"Identity=([0-9.]+)", attributes)

            # (genome start, end, query start, end, frameshifts, stops, identity)
            cds.setdefault(parent, []).append(
                (
                    int(fields[3]),
                    int(fields[4]),
                    int(target.group(2)),
                    int(target.group(3)),
                    int(frameshifts.group(1)) if frameshifts else 0,
                    int(stops.group(1)) if stops else 0,
                    float(identity.group(1)) if identity else 0.0,
                )
            )
    return mrna, cds


def count_alternative_loci(mrna):
    """Per protein, other hits covering >= ALT_LOCUS_COVERAGE of its best Rank=1 hit."""
    best_coverage = {}
    for hit in mrna.values():
        if hit["rank"] == 1:
            coverage = hit["qmax"] - hit["qmin"]
            best_coverage[hit["prot"]] = max(best_coverage.get(hit["prot"], 0), coverage)

    counts = Counter()
    for hit in mrna.values():
        best = best_coverage.get(hit["prot"], 0)
        if best <= 0:
            continue
        coverage = hit["qmax"] - hit["qmin"]
        if coverage >= ALT_LOCUS_COVERAGE * best and (hit["rank"] != 1 or coverage < best):
            counts[hit["prot"]] += 1
    return counts


def hit_record(hit, segments, n_alternative_loci):
    """One alignment with its segments, introns and frameshift/stop query spans."""
    minus = hit["strand"] == "-"
    in_query_order = sorted(segments, key=lambda seg: seg[2])

    # introns: (query residue, genomic start, end) between query-adjacent segments
    introns = []
    for upstream, downstream in zip(in_query_order, in_query_order[1:]):
        if minus:
            start, end = downstream[1] + 1, upstream[0] - 1
        else:
            start, end = upstream[1] + 1, downstream[0] - 1
        if end >= start:
            introns.append((upstream[3], start, end))

    return dict(
        contig=hit["contig"],
        strand=hit["strand"],
        qmin=hit["qmin"],
        qmax=hit["qmax"],
        ident=hit["ident"],
        segs=[(seg[0], seg[1]) for seg in sorted(segments, key=lambda seg: seg[0])],
        cds_q=[(seg[2], seg[3], seg[0], seg[1], seg[6]) for seg in in_query_order],
        introns=introns,
        _cov=hit["qmax"] - hit["qmin"],
        n_frameshift=sum(seg[4] for seg in in_query_order),
        n_stop=sum(seg[5] for seg in in_query_order),
        fs_spans=[(seg[2], seg[3]) for seg in in_query_order if seg[4]],
        stop_spans=[(seg[2], seg[3]) for seg in in_query_order if seg[5]],
        n_alternative_loci=n_alternative_loci,
    )


def parse_miniprot(gff):
    """protein -> its max-coverage Rank=1 miniprot alignment."""
    mrna, cds = read_miniprot_gff(gff)
    alternatives = count_alternative_loci(mrna)

    # the first Rank=1 hit of greatest coverage wins
    best = {}
    for mrna_id, hit in mrna.items():
        if hit["rank"] != 1:
            continue
        protein = hit["prot"]
        if protein in best and hit["qmax"] - hit["qmin"] <= best[protein]["_cov"]:
            continue
        best[protein] = hit_record(hit, cds[mrna_id], alternatives.get(protein, 0))
    return best


def local_ident(hit, residue, window):
    """Best identity of any aligned segment overlapping residue +/- window."""
    low, high = residue - window, residue + window
    return max(
        (
            identity
            for query_start, query_end, _, _, identity in hit["cds_q"]
            if min(query_end, high) >= max(query_start, low)
        ),
        default=0.0,
    )


def protein_alignment(focal_seq, other_seq):
    """Top global BLOSUM62 alignment of two proteins as (focal, other) gapped strings."""
    alignment = ALIGNER.align(focal_seq, other_seq)[0]
    return str(alignment[0]), str(alignment[1])


def flank_mismatches(focal_seq, ortho_seq, residue, flank_aa):
    """(mm5, mm3, n5, n3): mismatches and focal residues in the flank_aa windows either side."""
    if not focal_seq or not ortho_seq:
        return flank_aa, flank_aa, 0, 0

    aligned_focal, aligned_ortholog = protein_alignment(focal_seq, ortho_seq)
    low5, high5 = residue - flank_aa, residue
    low3, high3 = residue, residue + flank_aa
    focal_index = 0
    mismatch5 = mismatch3 = 0
    compared5 = compared3 = 0

    # walk the alignment by focal residue, ortholog gaps counted as mismatches
    for focal_char, ortholog_char in zip(aligned_focal, aligned_ortholog):
        if focal_char == "-":
            continue
        if low5 <= focal_index < high5:
            compared5 += 1
            mismatch5 += focal_char != ortholog_char
        elif low3 <= focal_index < high3:
            compared3 += 1
            mismatch3 += focal_char != ortholog_char
        focal_index += 1

    # uncompared focal window residues count as mismatches, so truncated hits cannot pass
    expected5 = max(0, min(high5, len(focal_seq)) - max(low5, 0))
    expected3 = max(0, min(high3, len(focal_seq)) - max(low3, 0))
    mismatch5 += max(0, expected5 - compared5)
    mismatch3 += max(0, expected3 - compared3)
    return mismatch5, mismatch3, expected5, expected3


def telomere_frac(seq, motif):
    """Best-phase positional identity to a motif tandem (random DNA scores ~0.3-0.4)."""
    seq = "".join(base for base in seq.upper() if base in "ACGT")
    if len(seq) < len(motif):
        return 0.0
    best = 0

    # best over both strands and every phase of the tandem
    for unit in (motif, rc(motif)):
        unit_len = len(unit)
        for phase in range(unit_len):
            tandem = (unit[phase:] + unit * (len(seq) // unit_len + 2))[: len(seq)]
            best = max(best, sum(1 for x, y in zip(seq, tandem) if x == y))
    return best / len(seq)


def nearby_introns(hit, residue):
    """(start, end) of the hit's introns within INTRON_TOL residues of the junction."""
    return [
        (start, end)
        for position, start, end in hit["introns"]
        if abs(position - residue) <= INTRON_TOL
    ]


def longest_intron(introns):
    """The longest of several (start, end) introns."""
    return max(introns, key=lambda intron: intron[1] - intron[0])


def intron_dna(contig, hit, start, end):
    """Comparator intron [start, end] on the hit's strand."""
    seq = contig[start - 1 : end]
    return rc(seq) if hit["strand"] == "-" else seq


def structurally_absent(hit, residue):
    """Residue inside one exon with ABSENT_MARGIN aligned frameshift-free residues each side."""
    positions = [position for position, _, _ in hit["introns"]]
    left = max([p for p in positions if p <= residue] + [hit["qmin"]])
    right = min([p for p in positions if p >= residue] + [hit["qmax"]])
    clean = not any(
        start <= residue + ABSENT_MARGIN and end >= residue - ABSENT_MARGIN
        for start, end in hit["fs_spans"] + hit["stop_spans"]
    )
    return residue - left >= ABSENT_MARGIN and right - residue >= ABSENT_MARGIN and clean


def introns_elsewhere(hit, residue):
    """Retrocopy guard: whether the hit has introns away from the junction."""
    elsewhere = sum(1 for position, _, _ in hit["introns"] if abs(position - residue) > INTRON_TOL)
    return elsewhere >= MIN_INTRONS_ELSEWHERE


def residue_position(hit, residue):
    """Genomic position of the residue's codon, else the hit's midpoint."""
    for query_start, query_end, start, end, _ in hit.get("cds_q", []):
        if query_start <= residue <= query_end:
            offset = (residue - query_start) * 3
            return end - offset if hit["strand"] == "-" else start + offset
    if hit.get("segs"):
        return (min(start for start, _ in hit["segs"]) + max(end for _, end in hit["segs"])) // 2
    return None


def site_clear(contig, position):
    """No N-gap excess within ABSENT_WINDOW and at least EDGE_BP from both contig ends."""
    window = contig[max(0, position - ABSENT_WINDOW) : min(len(contig), position + ABSENT_WINDOW)]
    window = window.upper()
    n_fraction = window.count("N") / max(1, len(window))
    far_from_ends = position >= EDGE_BP and len(contig) - position >= EDGE_BP
    return n_fraction <= MAX_N_FRAC and far_from_ends


def mafft(seqs, threads, is_protein):
    """MAFFT alignment of (name, sequence) pairs, fewer than two returned unaligned."""
    if len(seqs) < 2:
        return seqs

    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, dir=scratch_root()) as tf:
        for name, seq in seqs:
            tf.write(f">{name}\n{seq}\n")
        unaligned = tf.name

    alphabet = ["--amino"] if is_protein else ["--nuc"]
    command = ["mafft", "--auto", "--quiet", "--thread", str(threads)] + alphabet + [unaligned]
    try:
        out = subprocess.run(command, capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as error:
        os.unlink(unaligned)
        raise RuntimeError(f"MAFFT failed; refusing unaligned fallback: {error.stderr}") from error
    os.unlink(unaligned)

    # parse the aligned FASTA from stdout
    aligned = []
    name = None
    buf = []
    for line in out.splitlines():
        if line.startswith(">"):
            if name is not None:
                aligned.append((name, "".join(buf)))
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        aligned.append((name, "".join(buf)))
    return aligned


def write_fasta(path, seqs):
    """Write (name, sequence) pairs as FASTA."""
    with open(path, "w") as handle:
        for name, seq in seqs:
            handle.write(f">{name}\n{seq}\n")


def run_msa_job(job):
    """Align one (output path, sequences, is_protein) job single-threaded."""
    path, seqs, is_protein = job
    write_fasta(path, mafft(seqs, 1, is_protein=is_protein))


def load_catalogue(path, focal_ids):
    """Telotron rows of the focal genomes (mapped catalogue from host.py), keyed by genome."""
    by_genome = defaultdict(list)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["genome_id"] not in focal_ids:
                continue
            by_genome[row["genome_id"]].append(
                dict(
                    seqid=row["seqid"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    tx=row["tx_id"],
                    gene=row["gene_id"],
                    motif=(row.get("motif") or "").strip() or MOTIF,
                )
            )
    return by_genome


def map_hosts(rows_by_genome, focal_ids, refseq_dir, tara_dir):
    """Host proteins (once each) and telotrons placed on them as (residue, phase)."""
    proteins = {}
    telotrons = []
    for focal in focal_ids:
        fasta, gff = find_genome_files(focal, refseq_dir, tara_dir)
        if not fasta or not gff:
            print(f"[skip focal {focal}] missing fasta/gff", file=sys.stderr)
            continue
        contigs = read_fasta(fasta)
        models = load_cds(gff)

        # each telotron needs its host transcript's CDS gap to be exactly the intron
        for row in rows_by_genome.get(focal, []):
            model = models.get(row["tx"])
            if not model or not model["prot"] or row["seqid"] not in contigs:
                continue
            junction = intron_residue(model, row["start"], row["end"])
            if junction is None:
                continue
            residue, phase = junction
            contig = contigs[row["seqid"]]

            # record each host protein once, with its CDS span for the RBH check
            protein_id = model["prot"]
            if protein_id not in proteins:
                proteins[protein_id] = dict(
                    focal=focal,
                    seq=protein_of(model, contig),
                    contig=row["seqid"],
                    cds_lo=min(start for start, _ in model["segs"]),
                    cds_hi=max(end for _, end in model["segs"]),
                )

            # intron DNA oriented on the gene strand
            dna = contig[row["start"] - 1 : row["end"]]
            if model["strand"] == "-":
                dna = rc(dna)
            telotrons.append(
                dict(
                    tid=f"{focal}|{row['seqid']}:{row['start']}-{row['end']}",
                    focal=focal,
                    prot=protein_id,
                    gene=row["gene"],
                    res=residue,
                    phase=phase,
                    dna=dna,
                    motif=row["motif"],
                )
            )
    return proteins, telotrons


def forward_search(ortho_ids, query_faa, query_sha, work_dir, threads, refseq_dir, tara_dir):
    """Best miniprot hit per focal protein and an indexed reference, per comparator."""
    comparators = {}
    for oid in ortho_ids:
        fasta, _ = find_genome_files(oid, refseq_dir, tara_dir)
        if not fasta:
            print(f"[skip ortholog {oid}] missing fasta", file=sys.stderr)
            continue

        # the query SHA1 keys the cache, so a changed focal panel never reuses stale GFFs
        out_gff = os.path.join(work_dir, f"{oid}.{query_sha}.miniprot.gff")
        if not (os.path.exists(out_gff) and os.path.getsize(out_gff) > 0):
            run_miniprot(fasta, query_faa, threads, out_gff)
        comparators[oid] = dict(hits=parse_miniprot(out_gff), reference=Reference(fasta))
        n_aligned = len(comparators[oid]["hits"])
        print(f"  miniprot vs {oid}: {n_aligned} proteins aligned", file=sys.stderr)
    return comparators


def ortholog_proteins(comparators, protein_ids):
    """(comparator, protein) -> (frame-0 translation of the hit CDS, hit), frameshifts kept."""
    translated = {}
    for oid, comparator in comparators.items():
        for protein_id in protein_ids:
            hit = comparator["hits"].get(protein_id)
            if not hit:
                continue
            contig = comparator["reference"].get(hit["contig"], "")
            cds = "".join(contig[start - 1 : end] for start, end in hit["segs"])
            if hit["strand"] == "-":
                cds = rc(cds)
            translated[(oid, protein_id)] = (translate(cds).rstrip("*"), hit)
    return translated


def reciprocal_best_hits(
    ortho_ids, proteins, translated, focal_fastas, query_sha, work_dir, threads
):
    """(comparator, protein) -> whether the reverse hit overlaps the host protein's CDS span."""
    # group each comparator's translated hits by focal genome, skipping self-comparison
    groups = defaultdict(list)
    for oid in ortho_ids:
        for protein_id, protein in proteins.items():
            if oid == protein["focal"]:
                continue
            ortholog_seq, hit = translated.get((oid, protein_id), (None, None))
            if hit and hit["ident"] >= MIN_IDENT and ortholog_seq:
                groups[(oid, protein["focal"])].append((protein_id, ortholog_seq))

    # reverse search each group against its focal genome (cached)
    reciprocal = {}
    for (oid, focal), items in groups.items():
        if focal not in focal_fastas:
            continue
        stem = os.path.join(work_dir, f"rev.{oid}__{focal}.{query_sha}")
        write_fasta(f"{stem}.faa", items)
        if not (os.path.exists(f"{stem}.miniprot.gff") and os.path.getsize(f"{stem}.miniprot.gff")):
            run_miniprot(focal_fastas[focal], f"{stem}.faa", threads, f"{stem}.miniprot.gff")
        reverse_hits = parse_miniprot(f"{stem}.miniprot.gff")

        # reciprocal when the reverse hit overlaps the focal protein's CDS span
        for protein_id, _ in items:
            protein = proteins[protein_id]
            reverse_hit = reverse_hits.get(protein_id)
            overlaps = False
            if reverse_hit and reverse_hit["contig"] == protein["contig"]:
                low = min(start for start, _ in reverse_hit["segs"])
                high = max(end for _, end in reverse_hit["segs"])
                overlaps = low <= protein["cds_hi"] and high >= protein["cds_lo"]
            reciprocal[(oid, protein_id)] = overlaps
    return reciprocal


def call_pair(telotron, focal_seq, ortholog_seq, hit, contig, rbh_pass):
    """One telotron x comparator call: (map row fields, comparator intron DNA or None)."""
    residue = telotron["res"]

    # defaults for a pair without a usable alignment
    fields = dict(
        ident=0.0,
        n_frameshift=0,
        n_stop=0,
        fs_at_intron=False,
        ortho_intron_len=0.0,
        ortho_intron_telofrac=0.0,
        local_ident=0.0,
        n_alternative_loci=0,
        flank5_mm=FLANK_AA,
        flank3_mm=FLANK_AA,
        rbh_pass=0,
    )
    if not hit:
        fields["call"] = "NOALIGN"
        return fields, None
    if hit["ident"] < MIN_IDENT:
        fields["call"] = "NOCOV"
        return fields, None

    fields.update(
        ident=hit["ident"],
        n_frameshift=hit["n_frameshift"],
        n_stop=hit["n_stop"],
        n_alternative_loci=hit.get("n_alternative_loci", 0),
        local_ident=local_ident(hit, residue, ABSENT_MARGIN),
    )
    covered = hit["qmin"] - INTRON_TOL <= residue <= hit["qmax"] + INTRON_TOL
    fields["fs_at_intron"] = any(
        start - INTRON_TOL <= residue <= end + INTRON_TOL
        for start, end in hit["fs_spans"] + hit["stop_spans"]
    )

    # gates: conserved flanks, reciprocal best hit, no near-best paralog
    mismatch5, mismatch3, _, _ = flank_mismatches(focal_seq, ortholog_seq, residue, FLANK_AA)
    fields.update(flank5_mm=mismatch5, flank3_mm=mismatch3, rbh_pass=int(rbh_pass))
    flank_ok = mismatch5 <= MAX_FLANK_MISMATCH and mismatch3 <= MAX_FLANK_MISMATCH
    gate_ok = flank_ok and rbh_pass and fields["n_alternative_loci"] <= MAX_ALT_LOCI

    # structure is judged even when the gates fail, so such pairs surface as POOR_FLANK
    candidates = nearby_introns(hit, residue)
    present = absent = struct_absent = False
    intron_seq = None
    if candidates and gate_ok:
        # FILL: the longest candidate intron and its comparator DNA
        present = True
        intron_seq = intron_dna(contig, hit, *longest_intron(candidates))
        fields["ortho_intron_len"] = len(intron_seq)
        fields["ortho_intron_telofrac"] = telomere_frac(intron_seq, telotron["motif"])
    elif not candidates:
        # ABSENT: no intron near the junction, and the site is not an N-gap or contig edge
        struct_absent = structurally_absent(hit, residue)
        position = residue_position(hit, residue)
        site_ok = not contig or position is None or site_clear(contig, position)
        absent = struct_absent and gate_ok and introns_elsewhere(hit, residue) and site_ok

    # first match wins: FILL, ABSENT, POOR_FLANK, UNCERTAIN, NOCOV
    if present:
        fields["call"] = "FILL"
    elif absent:
        fields["call"] = "ABSENT"
    elif not gate_ok and (candidates or struct_absent):
        fields["call"] = "POOR_FLANK"
    elif covered:
        fields["call"] = "UNCERTAIN"
    else:
        fields["call"] = "NOCOV"
    return fields, intron_seq


def write_map(path, map_rows):
    """ortholog_map.tsv: one row per telotron x comparator."""
    with open(path, "w") as handle:
        handle.write(MAP_HEADER)
        for row in map_rows:
            handle.write(
                f"{row['tid']}\t{row['focal']}\t{row['gene']}\t{row['prot']}\t{row['res']}\t"
                f"{row['ortholog']}\t{row['ident']:.3f}\t{row['call']}\t{row['n_frameshift']}\t"
                f"{row['n_stop']}\t{int(row['fs_at_intron'])}\t"
                f"{row['ortho_intron_len']}\t{row['ortho_intron_telofrac']:.3f}\t"
                f"{row['local_ident']:.3f}\t{row['n_alternative_loci']}\t"
                f"{row['flank5_mm']}\t{row['flank3_mm']}\t{row['rbh_pass']}\t{row['phase']}\n"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--final", required=True, help="mapped catalogue (host.py)")
    ap.add_argument("--focal-ids", default="", help="comma telotron-positive genome_ids")
    ap.add_argument(
        "--ortholog-ids", required=True, help="comma genome_ids to search for orthologs"
    )
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    ensure_tool_on_path("miniprot", "mafft")

    # focal genomes to map and comparator genomes to search
    focal_ids = sorted(set(x for x in args.focal_ids.split(",") if x))
    if not focal_ids:
        raise SystemExit("no focal_ids (use --focal-ids)")
    print(f"telotron_ortholog_align: {len(focal_ids)} focal species", flush=True)
    ortho_ids = [x for x in args.ortholog_ids.split(",") if x]

    protein_msa_dir = os.path.join(args.outdir, "protein_msa")
    intron_dna_dir = os.path.join(args.outdir, "intron_dna")
    work_dir = os.path.join(args.outdir, "_work")
    for directory in (protein_msa_dir, intron_dna_dir, work_dir):
        os.makedirs(directory, exist_ok=True)

    # host proteins and telotron junctions
    rows_by_genome = load_catalogue(args.final, focal_ids)
    proteins, telotrons = map_hosts(rows_by_genome, focal_ids, args.refseq_dir, args.tara_dir)
    print(
        f"focal host proteins: {len(proteins)}; telotrons mapped: {len(telotrons)}", file=sys.stderr
    )
    query_faa = os.path.join(work_dir, "focal_host_proteins.faa")
    write_fasta(query_faa, [(protein_id, p["seq"]) for protein_id, p in proteins.items()])
    with open(query_faa, "rb") as handle:
        query_sha = hashlib.sha1(handle.read()).hexdigest()[:12]

    # forward search, ortholog translations and the reciprocal-best-hit check
    comparators = forward_search(
        ortho_ids, query_faa, query_sha, work_dir, args.threads, args.refseq_dir, args.tara_dir
    )
    translated = ortholog_proteins(comparators, proteins)
    focal_fastas = {}
    for focal in focal_ids:
        fasta, _ = find_genome_files(focal, args.refseq_dir, args.tara_dir)
        if fasta:
            focal_fastas[focal] = fasta
    reciprocal = reciprocal_best_hits(
        ortho_ids, proteins, translated, focal_fastas, query_sha, work_dir, args.threads
    )
    print(
        f"  RBH: {sum(reciprocal.values())}/{len(reciprocal)} ortholog ORF proteins reciprocal",
        file=sys.stderr,
    )

    # protein MSAs keep sub-floor hits (MSA_MIN_IDENT): every map row needs an MSA row
    msa_jobs = []
    for protein_id, protein in proteins.items():
        seqs = [(f"{protein['focal']}|{protein_id}|FOCAL|fs=0|stop=0", protein["seq"])]
        for oid in ortho_ids:
            ortholog_seq, hit = translated.get((oid, protein_id), (None, None))
            if oid == protein["focal"] or not hit or hit["ident"] < MSA_MIN_IDENT:
                continue
            if ortholog_seq:
                low = "|lowid" if hit["ident"] < MIN_IDENT else ""
                name = f"{oid}|id={hit['ident']:.2f}|fs={hit['n_frameshift']}"
                seqs.append((f"{name}|stop={hit['n_stop']}{low}", ortholog_seq))
        if len(seqs) >= 2:
            path = os.path.join(protein_msa_dir, f"{protein['focal']}__{protein_id}.aln.fasta")
            msa_jobs.append((path, seqs, True))

    # call every telotron x comparator pair and collect comparator intron DNA for the MSAs
    map_rows = []
    for telotron in telotrons:
        intron_seqs = [
            (f"{telotron['focal']}|TELOTRON|len={len(telotron['dna'])}", telotron["dna"])
        ]
        for oid in ortho_ids:
            if oid == telotron["focal"]:
                continue
            ortholog_seq, hit = translated.get((oid, telotron["prot"]), ("", None))
            contig = comparators[oid]["reference"].get(hit["contig"], "") if hit else ""
            rbh_pass = reciprocal.get((oid, telotron["prot"]), False)
            fields, intron_seq = call_pair(
                telotron, proteins[telotron["prot"]]["seq"], ortholog_seq, hit, contig, rbh_pass
            )
            map_rows.append(dict(telotron, ortholog=oid, **fields))
            if intron_seq and len(intron_seq) >= MIN_MSA_INTRON_BP:
                flag = ",FRAMESHIFT@junction" if fields["fs_at_intron"] else ""
                length = fields["ortho_intron_len"]
                telofrac = fields["ortho_intron_telofrac"]
                name = f"{oid}|intron|len={length}|telo={telofrac:.2f}{flag}"
                intron_seqs.append((name, intron_seq))

        # always write the focal telotron, even alone, under the name every reader opens
        stem = telotron["tid"].replace("|", "__").replace(":", "_")
        path = os.path.join(intron_dna_dir, f"{stem}.aln.fasta")
        if len(intron_seqs) == 1:
            write_fasta(path, intron_seqs)
        else:
            msa_jobs.append((path, intron_seqs, False))

    # reuse finished alignments whose inputs are unchanged (msa_cache)
    print(f"aligning {len(msa_jobs)} MSAs across {args.threads} workers...", file=sys.stderr)
    align = cached_msa_job(run_msa_job)
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        list(executor.map(align, msa_jobs))

    write_map(os.path.join(args.outdir, "ortholog_map.tsv"), map_rows)
    calls = Counter(row["call"] for row in map_rows)
    print(f"\nDONE. telotrons={len(telotrons)}  ortholog calls={len(map_rows)}", file=sys.stderr)
    print("  " + "  ".join(f"{call} {n}" for call, n in sorted(calls.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
