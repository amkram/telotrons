"""Shared sequence, table, FASTA, genome-location, scratch and orthology helpers."""

import csv
import glob
import gzip
import os
import re
import shutil

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq):
    """Reverse complement."""
    return seq.translate(_RC_TABLE)[::-1]


def slug(text):
    """Text with every character outside [A-Za-z0-9_] replaced by '_'."""
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


def open_maybe_gz(path):
    """Open a text file for reading, gzip-compressed when its name ends in .gz."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_tsv(path):
    """Rows of a TSV(.gz) file as dicts."""
    with open_maybe_gz(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    """Write dict rows to a TSV(.gz) file, columns in first-seen order. No rows, no file."""
    rows = list(rows)
    if not rows:
        return
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    """{record id: upper-case sequence}, keyed by the first word of each header."""
    seqs, seqid, chunks = {}, None, []
    with open_maybe_gz(path) as handle:
        for line in handle:
            line = line.rstrip()
            # a header closes the previous record
            if line.startswith(">"):
                if seqid is not None:
                    seqs[seqid] = "".join(chunks).upper()
                seqid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        if seqid is not None:
            seqs[seqid] = "".join(chunks).upper()
    return seqs


def _comparator_dir(refseq_dir):
    """Annotation-free comparator assemblies, beside the refseq directory."""
    return os.path.join(os.path.dirname(os.path.abspath(refseq_dir)), "comparators")


def _genbank_dir(refseq_dir):
    """GenBank downloads, beside the refseq directory."""
    return os.path.join(os.path.dirname(os.path.abspath(refseq_dir)), "genbank")


def find_genome_files(gid, refseq_dir, tara_dir):
    """(FASTA, GFF) paths of a genome, None for a missing file."""
    if gid.startswith("TARA"):
        fasta = glob.glob(f"{tara_dir}/Contigs/{gid}.fa*")
        # Tara MAG GFFs are named <gid>.gmove.gff
        gff = (
            glob.glob(f"{tara_dir}/Genes/GFF/{gid}*.gff*")
            or glob.glob(f"{tara_dir}/Annot/{gid}*.gff*")
            or glob.glob(f"{tara_dir}/**/{gid}*.gff*", recursive=True)
        )
    else:
        fasta = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.fna*", recursive=True)
        gff = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.gff*", recursive=True)

        # fill each missing side independently (partial refseq downloads)
        genbank = _genbank_dir(refseq_dir)
        if not fasta:
            fasta = glob.glob(f"{genbank}/{gid}/**/*_genomic.fna*", recursive=True)
        if not gff:
            gff = glob.glob(f"{genbank}/{gid}/**/*_genomic.gff*", recursive=True)

        # last resort: annotation-free comparators (no GFF), so an annotated copy always wins
        if not fasta:
            fasta = glob.glob(f"{_comparator_dir(refseq_dir)}/{gid}.fna*")

    # drop index sidecars, then prefer uncompressed FASTA
    fasta = [path for path in fasta if not path.endswith((".fai", ".gzi", ".bai", ".ndx"))]
    fasta = [path for path in fasta if path.endswith((".fna", ".fa", ".fasta"))] or fasta
    return (fasta[0] if fasta else None), (gff[0] if gff else None)


def ensure_tool_on_path(*tools):
    """Put each missing tool on PATH from the first snakemake environment that has it."""
    envs = os.path.expanduser("~/.snakemake-envs/*/bin")
    for tool in tools:
        if shutil.which(tool):
            continue
        for candidate in glob.glob(f"{envs}/{tool}"):
            os.environ["PATH"] = (
                os.path.dirname(candidate) + os.pathsep + os.environ.get("PATH", "")
            )
            break


# Never /tmp: SIGKILLed slurm jobs skip TemporaryDirectory cleanup and fill node roots.
def scratch_root():
    """$TMPDIR as a writable job scratch directory. Raises when unset or under /tmp."""
    scratch = (os.environ.get("TMPDIR") or "").strip()
    if not scratch:
        raise RuntimeError(
            "TMPDIR is unset; set it to a job-specific scratch directory (never /tmp)"
        )

    # test /tmp and its resolved target, since /tmp may be a symlink
    resolved = os.path.realpath(scratch)
    for tmp in {"/tmp", os.path.realpath("/tmp")}:
        if resolved == tmp or resolved.startswith(tmp + os.sep):
            raise RuntimeError(
                f"TMPDIR={scratch} resolves under /tmp; use a job-specific scratch directory"
            )

    os.makedirs(scratch, exist_ok=True)
    if not os.access(scratch, os.W_OK):
        raise OSError(f"{scratch} not writable")
    return scratch


def intron_observed(row, intron_seq):
    """Whether ortho_intron_len > 0 (it may be the string "0.0") and the intron is not all gaps."""
    try:
        length = float(row.get("ortho_intron_len") or 0)
    except (TypeError, ValueError):
        length = 0.0
    return length > 0 and bool((intron_seq or "").replace("-", ""))


def its_spans(seq, motif, min_units=4):
    """Non-overlapping ITS tracts (start, end, units, identity), <=1 substitution per unit pair."""
    from paired_repeat_tracts import tracts

    candidates = []
    for start, end, mismatches, _pattern in tracts(seq or '', motif, min_units):
        identity = 1 - mismatches / (end - start)
        candidates.append((start, end, (end - start) // len(motif), round(identity, 3)))

    # greedy left to right. At equal start the longest, then most identical, tract wins
    spans = []
    for span in sorted(candidates, key=lambda tract: (tract[0], -(tract[1] - tract[0]), -tract[3])):
        if not spans or span[0] >= spans[-1][1]:
            spans.append(span)
    return spans


def read_aln(path):
    """[(header, sequence)] of an aligned FASTA, empty when the file is missing."""
    records, name, chunks = [], None, []
    if not os.path.exists(path):
        return records

    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())

    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def locus_stub(tid):
    """'GCF_x|NW_1:100-200' -> 'GCF_x__NW_1_100-200' (the per-locus file stem)."""
    focal, rest = tid.split("|", 1)
    return f"{focal}__{rest.replace(':', '_')}"


def flank_gate(msa, focal, res, gate_aa, min_ident, slop):
    """{taxon: (ok, ident, n_indel)} over gate_aa focal residues each side of junction residue res

    ok needs identity >= min_ident, a full window and no indel beyond slop residues from it.
    """
    focal_row = next((seq for header, seq in msa if header.split("|")[0] == focal), None)
    if focal_row is None:
        return {}

    # junction: alignment column of focal residue res (0-based)
    junction, seen = None, 0
    for i, residue in enumerate(focal_row):
        if residue != "-":
            if seen == res:
                junction = i
                break
            seen += 1
    if junction is None:
        junction = len(focal_row)

    # widen to gate_aa focal residues each side, n5/n3 counting those found
    lo, count = junction, 0
    while lo > 0 and count < gate_aa:
        lo -= 1
        if focal_row[lo] != "-":
            count += 1
    n5 = count

    hi, count = junction, 0
    while hi < len(focal_row) and count < gate_aa:
        if focal_row[hi] != "-":
            count += 1
        hi += 1
    n3 = count

    # residues_to_junction[i]: focal residues from column i to the junction (slop is in residues)
    residues_to_junction, count = {}, 0
    for i in range(junction - 1, lo - 1, -1):
        if focal_row[i] != "-":
            count += 1
        residues_to_junction[i] = count
    count = 0
    for i in range(junction, hi):
        if focal_row[i] != "-":
            count += 1
        residues_to_junction[i] = count

    gates = {}
    for header, seq in msa:
        taxon = header.split("|")[0]
        if taxon == focal:
            continue

        # identity over shared residues, gaps outside slop counted as indels
        match = compared = indel = 0
        for i in range(lo, hi):
            # columns past the end of a shorter row count as gaps
            focal_residue, other = focal_row[i], (seq[i] if i < len(seq) else "-")
            if focal_residue == "-" and other == "-":
                continue
            if focal_residue == "-" or other == "-":
                if residues_to_junction.get(i, 99) > slop:
                    indel += 1
                continue
            compared += 1
            match += focal_residue == other
        ident = match / compared if compared else 0.0

        # a window truncated at either end of the focal protein fails
        full = n5 >= gate_aa and n3 >= gate_aa
        gates[taxon] = (bool(full and indel == 0 and ident >= min_ident), ident, indel)
    return gates


def repeat_spans(seq, motifs):
    """Per-base labels: F/R exact unit, f/r two-unit window with <=1 substitution, '.' otherwise."""
    from repeat_policy import repeat_labels

    return repeat_labels(seq, motifs)
