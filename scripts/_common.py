"""Shared helpers for the telotron pipeline scripts.

These were previously copy-pasted (in near-identical form) across ~16 scripts:
sequence/motif utilities, FASTA readers, genome-file location, gz decompression,
a tqdm fallback, and the figure palette. Each script is run as
`python scripts/<x>.py`, so `scripts/` is on sys.path[0] and `from _common
import ...` resolves with no package install.

The genome-FASTA finder is the one place callers historically diverged (some
raised FileNotFoundError, some returned None; some filtered .fai/.gzi sidecars,
some did not). The unified version filters sidecars always (strictly safer) and
takes a `required` flag to select raise-vs-None.
"""
import glob
import gzip
import os
import re
import shutil
import subprocess

# --------------------------------------------------------------------------- #
# Sequence / motif utilities
# --------------------------------------------------------------------------- #
_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    """Reverse complement (case-preserving; N maps to N)."""
    return s.translate(_RC_TABLE)[::-1]


# scan_telotrons.py historically force-uppercased; its inputs are already upper,
# so the case-preserving rc is equivalent. Kept as an alias for clarity.
revcomp = rc
reverse_complement = rc


def rotations(seq):
    """All cyclic rotations of an (upper-cased) motif."""
    m = seq.upper()
    return [m[i:] + m[:i] for i in range(len(m))]


def expand_motifs(motifs):
    """Every rotation and reverse complement of each motif, deduped + sorted."""
    out = set()
    for m in motifs:
        for r in rotations(m):
            out.add(r)
            out.add(rc(r))
    return sorted(out)


def all_variants(motif):
    """All rotations of a single motif and of its reverse complement, as a set."""
    return set(expand_motifs([motif]))


def slug(s):
    """Filesystem-safe slug: non-alphanumerics -> underscore."""
    return re.sub(r"[^A-Za-z0-9_]", "_", str(s))


_slug = slug  # legacy name used across plot scripts


# --------------------------------------------------------------------------- #
# FASTA I/O
# --------------------------------------------------------------------------- #
def open_maybe_gz(path):
    """Text-mode handle, transparently decompressing .gz."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """Read a (optionally gzipped) FASTA into {id: upper_seq}. IDs are split on
    whitespace (first token)."""
    seqs, sid, buf = {}, None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf).upper()
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            seqs[sid] = "".join(buf).upper()
    return seqs


load_fasta = read_fasta  # legacy name


def read_fasta_desc(path):
    """Like read_fasta but also returns descriptions: ({id: seq}, {id: description}).
    The description is the header text after the first whitespace token."""
    seqs, descs, sid, buf = {}, {}, None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf).upper()
                head = line[1:]
                sid = head.split()[0]
                descs[sid] = head[len(sid):].strip()
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            seqs[sid] = "".join(buf).upper()
    return seqs, descs


# --------------------------------------------------------------------------- #
# Genome file location (RefSeq tree + Tara Contigs/Annot layout)
# --------------------------------------------------------------------------- #
# GenBank-only (GCA_) assemblies land in a sibling dir. Any resolver that
# takes `refseq_dir` also checks the parallel `genbank_dir` (default:
# refseq_dir/../genbank). Callers can override via the env var
# TELOTRON_GENBANK_DIR when running scripts by hand.
def _genbank_dir(refseq_dir):
    override = os.environ.get("TELOTRON_GENBANK_DIR")
    if override:
        return override
    parent = os.path.dirname(os.path.abspath(refseq_dir))
    return os.path.join(parent, "genbank")


def find_genome_fasta(gid, refseq_dir, tara_dir, required=True):
    """Locate a genome's nucleotide FASTA.

    RefSeq (GCF_):   {refseq_dir}/{gid}/**/*_genomic.fna[.gz]
    GenBank (GCA_):  {genbank_dir}/{gid}/**/*_genomic.fna[.gz]  (sibling of refseq_dir)
    Tara:            {tara_dir}/Contigs/{gid}.fa[.gz]
    Index/aux sidecars are excluded and a plain FASTA is preferred over a .gz
    (miniprot etc. want the uncompressed file). Returns the first such match;
    if none, raises FileNotFoundError when required else returns None.
    """
    if gid.startswith("TARA"):
        hits = glob.glob(f"{tara_dir}/Contigs/{gid}.fa*")
    else:
        hits = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.fna*", recursive=True)
        if not hits:
            hits = glob.glob(f"{_genbank_dir(refseq_dir)}/{gid}/**/*_genomic.fna*", recursive=True)
    hits = [p for p in hits if not p.endswith((".fai", ".gzi", ".bai", ".ndx"))]
    plain = [p for p in hits if p.endswith((".fna", ".fa", ".fasta"))]
    gz = [p for p in hits if p.endswith(".gz")]
    pick = plain or gz or hits
    if pick:
        return pick[0]
    if required:
        raise FileNotFoundError(f"no genome FASTA for {gid}")
    return None


def find_genome_files(gid, refseq_dir, tara_dir):
    """Locate (fasta, gff) for a genome; either element may be None.

    Checks {refseq_dir}/{gid}/ then {genbank_dir}/{gid}/ (parallel dir) then
    Tara. Necessary for the GenBank front-end: GCA_-only assemblies live
    under data/raw/genbank/, not data/raw/refseq/.
    """
    if gid.startswith("TARA"):
        fa = glob.glob(f"{tara_dir}/Contigs/{gid}.fa*")
        # Tara MAG GFFs are named `<gid>.gmove.gff` (Genoscope gmove), so the
        # glob needs a `*` between gid and `.gff` — the prior `{gid}.gff*` never
        # matched them (left find_genome_files returning gff=None for every MAG).
        gff = (glob.glob(f"{tara_dir}/Genes/GFF/{gid}*.gff*")
               or glob.glob(f"{tara_dir}/Annot/{gid}*.gff*")
               or glob.glob(f"{tara_dir}/**/{gid}*.gff*", recursive=True))
    else:
        fa = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.fna*", recursive=True)
        gff = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.gff*", recursive=True)
        # Fill EITHER missing side from the parallel genbank dir independently.
        # Previously used `and` — a partial refseq download (only fna, no gff)
        # blocked the genbank fallback for the missing gff, silently returning
        # gff=None even though the file existed on disk.
        gb = _genbank_dir(refseq_dir)
        if not fa:
            fa = glob.glob(f"{gb}/{gid}/**/*_genomic.fna*", recursive=True)
        if not gff:
            gff = glob.glob(f"{gb}/{gid}/**/*_genomic.gff*", recursive=True)
    fa = [p for p in fa if not p.endswith((".fai", ".gzi", ".bai", ".ndx"))]
    fa = [p for p in fa if p.endswith((".fna", ".fa", ".fasta"))] or fa
    return (fa[0] if fa else None), (gff[0] if gff else None)


# --------------------------------------------------------------------------- #
# Decompression / indexing
# --------------------------------------------------------------------------- #
def maybe_decompress(src, tmpdir):
    """Return a plain path; if src is .gz, stream-decompress it into tmpdir."""
    if src.endswith(".gz"):
        out = os.path.join(tmpdir, os.path.basename(src)[:-3])
        with gzip.open(src, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return out
    return src


def ensure_tool_on_path(*tools, envs_glob="~/.snakemake-envs/*/bin"):
    """Make CLI tools discoverable when a script is run as a bare `python
    scripts/X.py` outside `snakemake --use-conda`. For each tool not already on
    PATH, search the local snakemake conda envs and prepend the bin dir that
    actually contains it. No-ops when the tool already resolves. Replaces the
    brittle hardcoded env-hash PATH hacks that broke on every env rebuild.
    """
    base = os.path.expanduser(envs_glob)
    for t in tools:
        if shutil.which(t):
            continue
        for cand in glob.glob(f"{base}/{t}"):
            d = os.path.dirname(cand)
            os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")
            break


def ensure_uncompressed_faidx(fa_path, work_dir):
    """Return path to an uncompressed, samtools-faidx'd FASTA (decompressing into
    work_dir if needed)."""
    if fa_path.endswith(".gz"):
        local = os.path.join(work_dir, os.path.basename(fa_path)[:-3])
        if not os.path.exists(local):
            with gzip.open(fa_path, "rb") as ih, open(local, "wb") as oh:
                shutil.copyfileobj(ih, oh)
        fa_path = local
    if not os.path.exists(fa_path + ".fai"):
        subprocess.run(["samtools", "faidx", fa_path], check=True)
    return fa_path


# --------------------------------------------------------------------------- #
# Progress bar (no hard tqdm dependency)
# --------------------------------------------------------------------------- #
try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback when tqdm absent
    def tqdm(it=None, **kw):
        return it if it is not None else iter(())
    tqdm.write = print  # type: ignore[attr-defined]


# --------------------------------------------------------------------------- #
# Figure palette + architecture labels (single source of truth)
# --------------------------------------------------------------------------- #
TEAL = "#20808D"
RUST = "#A84B2F"
GOLD = "#FFC553"
GREY = "#AAAAAA"

ORIENT_COLORS = {
    "Fwd/G-rich": TEAL,
    "Rev/C-rich": RUST,
    # Keys must match scan_telotrons.classify_orientation outputs
    # ("Fwd/G-rich", "Rev/C-rich", "Mixed", "None"). Previously this map
    # used "Bidirectional"/"Unknown", which never matched and silently
    # collapsed bidirectional + linker loci to GREY (per
    # verify_figures_visualization Finding F, 2026-06-04).
    "Mixed": GOLD,
    "None": GREY,
}

# Splice-architecture categories emitted by classify_telotron_architecture.py.
# SINGLE SOURCE OF TRUTH — this list must stay in sync with that script's
# classify() return values. Consumers importing ARCH_ORDER/BIDIR_ARCHS get the
# labels for free; a hand-copied list silently drops whole classes (this list
# omitted "Multi-junction" until 2026-07-16, so every per-architecture plot
# dropped that class without warning).
ARCH_ORDER = [
    "GT-F-AG", "GT-R-AG", "GT-F-R-AG", "GT-R-F-AG",
    "GT-F-linker-R-AG", "GT-R-linker-F-AG", "Multi-junction", "Other",
]
ARCH_COLORS = {
    "GT-F-AG":          TEAL,
    "GT-R-AG":          RUST,
    "GT-F-R-AG":        "#6A8D5A",
    "GT-R-F-AG":        "#4E7A6B",
    "GT-F-linker-R-AG": GOLD,
    "GT-R-linker-F-AG": "#8B5A99",
    "Multi-junction":   "#7A7A4A",
    "Other":            "#BBBBBB",
}

# Architectures in which BOTH arms of the telomeric repeat are present — the
# distinctive telomerase-mediated signature that gates confident_species'
# second admission path. GT-F-R-AG (convergent) and GT-R-F-AG (divergent) are
# kept as SEPARATE labels because dyad orientation is itself a result (~99% of
# Mixed loci are G-rich-5'), but both count as bidirectional.
BIDIR_ARCHS = frozenset({
    "GT-F-R-AG", "GT-R-F-AG",
    "GT-F-linker-R-AG", "GT-R-linker-F-AG",
    "Multi-junction",
})


# --------------------------------------------------------------------------- #
# Focal-panel accession <-> species map (single source of truth)
# --------------------------------------------------------------------------- #
# The 5 RefSeq Eimeria genomes + the 2 outgroup coccidia used as the focal
# panel for ortholog/architecture work. Ground truth verified from each
# assembly's FASTA description lines and the strain code in its filename:
#   GCF_000499385.1 / ENH001 -> Eimeria necatrix Houghton
#   GCF_000499425.1 / EAH001 -> Eimeria acervulina Houghton
#   GCF_000499545.2 / ETH001 -> Eimeria tenella Houghton
#   GCF_000499605.1 / EMW001 -> Eimeria maxima Weybridge
#   GCF_000499745.2 / EMH001 -> Eimeria mitis Houghton
# NB: Eimeria brunetti and Eimeria praecox are NOT in this panel; if those
# names appear in a downstream TSV the upstream script had a swapped map.
SPECIES_BY_ACC = {
    "GCF_000499385.1": "Eimeria_necatrix",
    "GCF_000499425.1": "Eimeria_acervulina",
    "GCF_000499545.2": "Eimeria_tenella",
    "GCF_000499605.1": "Eimeria_maxima",
    "GCF_000499745.2": "Eimeria_mitis",
    "GCF_000006565.2": "Toxoplasma_gondii",
    "GCF_000208865.1": "Neospora_caninum",
}

# Short 4-letter abbreviations used in locus_text rows and downstream TSVs.
SPECIES_ABBREV = {
    "Eimeria_necatrix":   "Enec",
    "Eimeria_acervulina": "Eace",
    "Eimeria_tenella":    "Eten",
    "Eimeria_maxima":     "Emax",
    "Eimeria_mitis":      "Emit",
    "Toxoplasma_gondii":  "Tgon",
    "Neospora_caninum":   "Ncan",
}
