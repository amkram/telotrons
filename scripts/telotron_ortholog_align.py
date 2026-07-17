#!/usr/bin/env python3
"""Telotron ortholog alignment.

For every telotron-containing gene in the focal (telotron-positive) genomes:
  1. find the orthologous locus in each panel genome (miniprot, intron-aware);
  2. align the gene in PROTEIN space (introns removed) across focal + orthologs (MAFFT);
  3. for orthologs that carry an intron at the homologous position, extract that
     intronic DNA and align the telotron's DNA (the telomeric array + linker) to it
     (MAFFT) -- testing whether the orthologous intron holds telomere-like sequence.

Inputs are only genome FASTA + GFF (resolved via _common.find_genome_files), so the
rule is self-contained on the pipeline's downloads. Outputs to --outdir:
  ortholog_map.tsv          per (telotron, ortholog) locus call: ortholog gene, %id,
                            intron-at-locus (fill), ortholog-intron length + telomere frac
  protein_msa/<focal>__<protein>.aln.fasta    protein MSA per host gene
  intron_dna/<telotron>.aln.fasta             telotron-vs-ortholog-intron DNA MSA
  summary.tsv               per telotron rollup; gene_summary.tsv per host gene

This is the reproducible form of analysis/telotron_origin/homologs/ (fill-vs-create).
"""
import argparse
import concurrent.futures as cf
import hashlib
import os
import re
import subprocess as sp
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import find_genome_files, read_fasta, rc, ensure_tool_on_path  # noqa: E402

CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(seq):
    seq = seq.upper()
    return "".join(CODON.get(seq[i:i + 3], 'X') for i in range(0, len(seq) - 2, 3))


def open_maybe_gz(path):
    import gzip
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def extract_flanks(contig, gs, ge, strand, flank):
    """Exonic DNA flanking an intron at 1-based inclusive genomic [gs, ge], returned
    transcript-oriented as (5'-flank, 3'-flank): the exon DNA abutting the donor and
    acceptor respectively, `flank` bp each."""
    up = contig[max(0, gs - 1 - flank):gs - 1]   # genomic upstream of intron
    dn = contig[ge:ge + flank]                    # genomic downstream of intron
    if strand == "-":
        return rc(dn), rc(up)
    return up, dn


# ── GFF: per-transcript CDS structure ───────────────────────────────────────
def load_cds(gff):
    """tx_id -> dict(strand, prot, segs=[(start,end)] genomic-ascending)."""
    tx = {}
    with open_maybe_gz(gff) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"Parent=([^;]+)", f[8])
            if not m:
                continue
            t = m.group(1)
            # default protein key = transcript id (covers ab-initio gmove GFFs with no
            # protein_id, e.g. Tara MAGs); overridden by protein_id below when present
            d = tx.setdefault(t, {"strand": f[6], "prot": t, "segs": []})
            d["segs"].append((int(f[3]), int(f[4])))
            p = re.search(r"protein_id=([^;]+)", f[8])
            if p:
                d["prot"] = p.group(1)
    return tx


def ordered_segs(d):
    """CDS segments in transcription order (5'->3')."""
    return sorted(d["segs"], key=lambda x: -x[0] if d["strand"] == "-" else x[0])


def protein_of(d, contig_seq):
    """Translate a transcript's CDS to protein from the contig sequence."""
    parts = []
    for s, e in ordered_segs(d):
        sub = contig_seq[s - 1:e]
        parts.append(rc(sub) if d["strand"] == "-" else sub)
    return translate("".join(parts)).rstrip("*")


def intron_residue(d, ts, te):
    """Protein residue at which the intron [ts,te] sits; None if no CDS gap matches.

    Returns (residue, phase) where phase ∈ {0,1,2} is the codon-position of the
    intron in the upstream exon (0 = intron between codons, 1 = after 1st base,
    2 = after 2nd base). Phase ≠ 0 means the intron splits a codon; the caller
    can use this to translate flanking exon DNA in the correct frame. Reports
    a stderr diagnostic when any preceding CDS exon length is not a multiple
    of 3, since the cumulative count then drifts (focal-genome frameshift or
    partial CDS at contig edge).
    """
    segs = ordered_segs(d)
    rev = d["strand"] == "-"
    cum = 0
    drift = False
    for i in range(len(segs) - 1):
        a, b = segs[i], segs[i + 1]
        gap = (b[1] + 1, a[0] - 1) if rev else (a[1] + 1, b[0] - 1)
        seg_len = a[1] - a[0] + 1
        if seg_len % 3 != 0:
            drift = True
        cum += seg_len
        if abs(gap[0] - ts) <= 2 and abs(gap[1] - te) <= 2:
            res, phase = divmod(cum, 3)
            if drift:
                print(f"[intron_residue] CDS exon length not %3 upstream of "
                      f"intron {ts}-{te} (strand {d['strand']}); residue/phase may drift",
                      file=sys.stderr)
            return res, phase
    return None


# ── miniprot ─────────────────────────────────────────────────────────────────
def run_miniprot(miniprot, genome_fa, query_faa, threads, out_gff):
    with open(out_gff, "w") as fh:
        sp.run([miniprot, "-t", str(threads), "--gff", genome_fa, query_faa],
               stdout=fh, stderr=sp.DEVNULL, check=True)


def parse_miniprot(gff):
    """protein -> best alignment dict(contig, strand, qmin, qmax, ident,
       cds=[(gs,ge)], introns=[(qpos, gstart, gend)])."""
    mrna = {}
    cds = {}
    for line in open(gff):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "mRNA":
            mid = re.search(r"ID=([^;]+)", f[8]).group(1)
            rank = int(re.search(r"Rank=(\d+)", f[8]).group(1))
            ident = float(re.search(r"Identity=([0-9.]+)", f[8]).group(1))
            tg = re.search(r"Target=(\S+) (\d+) (\d+)", f[8])
            mrna[mid] = dict(prot=tg.group(1), rank=rank, ident=ident,
                             contig=f[0], strand=f[6],
                             qmin=int(tg.group(2)), qmax=int(tg.group(3)))
            cds[mid] = []
        elif f[2] == "CDS":
            par = re.search(r"Parent=([^;]+)", f[8]).group(1)
            tg = re.search(r"Target=(\S+) (\d+) (\d+)", f[8])
            fsm = re.search(r"Frameshift=(\d+)", f[8])
            stm = re.search(r"StopCodon=(\d+)", f[8])
            idm = re.search(r"Identity=([0-9.]+)", f[8])   # per-exon local identity
            cds.setdefault(par, []).append(
                (int(f[3]), int(f[4]), int(tg.group(2)), int(tg.group(3)),
                 int(fsm.group(1)) if fsm else 0, int(stm.group(1)) if stm else 0,
                 float(idm.group(1)) if idm else 0.0))
    best = {}
    # n_alternative_loci: count of Rank>=2 hits per protein whose query coverage is
    # within 10% of Rank=1 (paralog candidates that the Rank-1-only rule silently
    # drops; Wolf & Koonin 2012 doi:10.1093/molbev/mss165 caution against
    # single-best-hit ortholog assignment). Also track Rank=1 alternatives that
    # lose the max-coverage tie-break in this same protein.
    alt_counts = {}
    best_cov = {}
    for mid, m in mrna.items():
        if m["rank"] != 1:
            continue
        cov = m["qmax"] - m["qmin"]
        best_cov[m["prot"]] = max(best_cov.get(m["prot"], 0), cov)
    for mid, m in mrna.items():
        bc = best_cov.get(m["prot"], 0)
        if bc <= 0:
            continue
        cov = m["qmax"] - m["qmin"]
        # Rank>=2 within-10%-coverage and Rank=1 alternates that lost the tie-break
        # both count as "alternative loci" worth flagging.
        if cov >= 0.9 * bc and (m["rank"] != 1 or cov < bc):
            alt_counts[m["prot"]] = alt_counts.get(m["prot"], 0) + 1
    for mid, m in mrna.items():
        if m["rank"] != 1:
            continue
        cov = m["qmax"] - m["qmin"]
        if m["prot"] in best and cov <= best[m["prot"]]["_cov"]:
            continue
        segs = sorted(cds[mid], key=lambda x: x[0])  # genomic ascending
        rev = m["strand"] == "-"
        tx_order = sorted(cds[mid], key=lambda x: x[2])  # query ascending = transcript order
        introns = []
        for i in range(len(tx_order) - 1):
            a, b = tx_order[i], tx_order[i + 1]
            qpos = a[3]                       # query residue at intron
            gs, ge = (b[1] + 1, a[0] - 1) if rev else (a[1] + 1, b[0] - 1)
            if ge >= gs:
                introns.append((qpos, gs, ge))
        # frameshifts/stops with their query ranges (NOT corrected — reported as-is)
        fs_spans = [(qs, qe) for _, _, qs, qe, fs, _, _ in tx_order if fs]
        stop_spans = [(qs, qe) for _, _, qs, qe, _, st, _ in tx_order if st]
        best[m["prot"]] = dict(contig=m["contig"], strand=m["strand"],
                               qmin=m["qmin"], qmax=m["qmax"], ident=m["ident"],
                               segs=[(s, e) for s, e, _, _, _, _, _ in segs],
                               cds_q=[(qs, qe, s, e, iden) for s, e, qs, qe, _, _, iden in tx_order],
                               introns=introns, _cov=cov,
                               n_frameshift=sum(x[4] for x in tx_order),
                               n_stop=sum(x[5] for x in tx_order),
                               fs_spans=fs_spans, stop_spans=stop_spans,
                               n_alternative_loci=alt_counts.get(m["prot"], 0))
    return best


def local_ident(a, r, window):
    """Best per-exon alignment identity (miniprot's own, so it is phase- and frameshift-aware)
    among the exons overlapping the +-window region around the homologous position r. Taking
    the MAX means a frameshift garbling one flank doesn't sink a real orthologue (its other,
    frame-stable exon still scores), while a spurious hit is low on EVERY exon -> filtered."""
    lo, hi = r - window, r + window
    return max((iden for qs, qe, gs, ge, iden in a["cds_q"] if min(qe, hi) >= max(qs, lo)),
               default=0.0)


# ── strict per-flank amino-acid mismatch gate ────────────────────────────────
# A global protein alignment of the focal host protein vs the ortholog ORF protein
# gives an unambiguous, frameshift-robust column map; we then walk `flank_aa` focal
# non-gap residues on each side of the homologous junction (focal residue `res`) and
# count columns where focal != ortholog (a gap counts as a mismatch). This REPLACES
# the lenient local_ident MAX gate as the primary FILL/ABSENT admission test: a real
# orthologous junction has near-identical flanks, a spurious / paralog / out-of-register
# hit does not. local_ident is still computed for back-compat reporting.
try:
    import warnings as _warnings
    _warnings.filterwarnings("ignore")
    from Bio.Align import PairwiseAligner as _PairwiseAligner, substitution_matrices as _submat
    _BLOSUM62 = _submat.load("BLOSUM62")
    _ALIGNER = _PairwiseAligner()
    _ALIGNER.substitution_matrix = _BLOSUM62
    _ALIGNER.open_gap_score = -11
    _ALIGNER.extend_gap_score = -1
    _ALIGNER.mode = "global"
    HAVE_BIOPYTHON = True
except Exception:  # pragma: no cover - biopython optional
    _ALIGNER = None
    HAVE_BIOPYTHON = False


def _nw_align(a, b):
    """Compact global Needleman-Wunsch fallback (BLOSUM62-free, identity scoring) used
    only when biopython is unavailable. Returns (aln_a, aln_b) gapped strings."""
    n, m = len(a), len(b)
    GAP = -4
    MATCH, MISMATCH = 2, -1
    # score + traceback matrices
    sc = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        sc[i][0] = i * GAP
    for j in range(1, m + 1):
        sc[0][j] = j * GAP
    tb = [[0] * (m + 1) for _ in range(n + 1)]  # 0=diag 1=up(gap in b) 2=left(gap in a)
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        tb[i][0] = 1
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = sc[i], sc[i - 1]
        trow = tb[i]
        for j in range(1, m + 1):
            d = prev[j - 1] + (MATCH if ai == b[j - 1] else MISMATCH)
            u = prev[j] + GAP
            l = row[j - 1] + GAP
            best = d
            t = 0
            if u > best:
                best, t = u, 1
            if l > best:
                best, t = l, 2
            row[j] = best
            trow[j] = t
    i, j = n, m
    oa, ob = [], []
    while i > 0 or j > 0:
        t = tb[i][j]
        if i > 0 and j > 0 and t == 0:
            oa.append(a[i - 1]); ob.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and t == 1:
            oa.append(a[i - 1]); ob.append("-"); i -= 1
        else:
            oa.append("-"); ob.append(b[j - 1]); j -= 1
    return "".join(reversed(oa)), "".join(reversed(ob))


def _mafft_pair_align(a, b):
    """Last-resort pairwise alignment via MAFFT (when biopython absent and the inline NW
    is undesirable for very long proteins). Returns (aln_a, aln_b)."""
    aln = mafft([("a", a), ("b", b)], 1, is_protein=True)
    d = dict(aln)
    return d.get("a", a), d.get("b", b)


def _global_protein_align(focal_seq, ortho_seq):
    """Global protein alignment -> (aln_focal, aln_ortho) equal-length gapped strings.
    Prefers biopython PairwiseAligner(global, BLOSUM62); falls back to inline NW; for very
    long sequences without biopython, MAFFT."""
    if HAVE_BIOPYTHON:
        try:
            aln = _ALIGNER.align(focal_seq, ortho_seq)[0]
            # Bio.Align.Alignment indexing yields the two gapped rows as strings.
            return str(aln[0]), str(aln[1])
        except Exception:
            pass
    if max(len(focal_seq), len(ortho_seq)) <= 4000:
        return _nw_align(focal_seq, ortho_seq)
    return _mafft_pair_align(focal_seq, ortho_seq)


def flank_mismatches(focal_seq, ortho_seq, res, flank_aa):
    """Count amino-acid mismatches in the 5' flank and the 3' flank SEPARATELY, over a
    window of `flank_aa` focal residues immediately flanking the homologous junction:
      5'-flank = focal residues [res-flank_aa, res)
      3'-flank = focal residues [res, res+flank_aa)
    A column is a mismatch when focal != ortholog (a gap in either counts as a mismatch).
    Returns (flank5_mm, flank3_mm, flank5_len, flank3_len) where *_len is the number of
    focal residues actually compared (< flank_aa near a terminus). Returns large mismatch
    counts if the alignment could not be produced so the locus fails the gate."""
    if not focal_seq or not ortho_seq:
        return flank_aa, flank_aa, 0, 0
    af, ao = _global_protein_align(focal_seq, ortho_seq)
    # walk alignment columns tracking the focal-residue index (0-based count of focal
    # non-gap residues consumed). The junction sits between focal residue index res-1
    # and res (intron_residue returns a 1-based-ish residue count = cumulative codons),
    # so 5'-flank = focal indices [res-flank_aa, res), 3'-flank = [res, res+flank_aa).
    lo5, hi5 = res - flank_aa, res
    lo3, hi3 = res, res + flank_aa
    fidx = 0
    mm5 = mm3 = 0
    n5 = n3 = 0
    for cf, co in zip(af, ao):
        if cf != "-":
            # this column carries focal residue with index fidx (0-based)
            if lo5 <= fidx < hi5:
                n5 += 1
                if cf != co:
                    mm5 += 1
            elif lo3 <= fidx < hi3:
                n3 += 1
                if cf != co:
                    mm3 += 1
            fidx += 1
        else:
            # gap in focal (insertion in ortholog): does not consume a focal residue, but
            # if it falls inside a flank window it does not penalise the focal-residue
            # count; the per-focal-residue mismatch accounting above is what the spec asks.
            pass
    # focal residues in a flank window that the alignment never emitted (focal shorter than
    # res+flank_aa, or terminus) are counted as missing -> mismatches, so a truncated /
    # non-covering ortholog cannot pass by simply having fewer compared positions.
    exp5 = max(0, min(hi5, len(focal_seq)) - max(lo5, 0))
    exp3 = max(0, min(hi3, len(focal_seq)) - max(lo3, 0))
    mm5 += max(0, exp5 - n5)
    mm3 += max(0, exp3 - n3)
    return mm5, mm3, exp5, exp3


def absent_flanks(octg, a, r, flank):
    """For an ABSENT orthologue (no intron at residue r): extract the contiguous exon DNA
    flanking the homologous junction (between codon r and r+1), transcript-oriented, so the
    intron-length deletion is visible against the focal's flank|intron|flank. Reuses
    extract_flanks with a zero-width intron at the junction position."""
    exon = next(((qs, qe, gs, ge) for qs, qe, gs, ge, _ in a["cds_q"] if qs <= r <= qe), None)
    if not exon:
        return "", ""
    qs, qe, gs, ge = exon
    if a["strand"] == "-":
        g = ge - 3 * (r - qs) - 2          # low-genomic base of codon r
        return extract_flanks(octg, g, g - 1, "-", flank)
    p = gs + 3 * (r - qs) + 2              # high-genomic (3') base of codon r
    return extract_flanks(octg, p + 1, p, "+", flank)


# ── telomere content of a sequence ───────────────────────────────────────────
def telomere_frac(seq, motif):
    """Best-phase fractional match of seq to a tandem (motif) on either strand."""
    seq = "".join(c for c in seq.upper() if c in "ACGT")
    if len(seq) < len(motif):
        return 0.0
    best = 0
    for unit in (motif, rc(motif)):
        L = len(unit)
        for ph in range(L):
            patt = (unit[ph:] + unit * (len(seq) // L + 2))[:len(seq)]
            best = max(best, sum(1 for a, b in zip(seq, patt) if a == b))
    return best / len(seq)


def mafft(seqs, threads, is_protein):
    """seqs: list of (name, seq). Returns aligned list of (name, seq)."""
    if len(seqs) < 2:
        return seqs
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        for n, s in seqs:
            tf.write(f">{n}\n{s}\n")
        inp = tf.name
    try:
        amino = ["--amino"] if is_protein else ["--nuc"]
        out = sp.run(["mafft", "--auto", "--quiet", "--thread", str(threads)] + amino + [inp],
                     capture_output=True, text=True, check=True).stdout
    except sp.CalledProcessError:
        os.unlink(inp)
        return seqs
    os.unlink(inp)
    aln = []
    name = None
    buf = []
    for line in out.splitlines():
        if line.startswith(">"):
            if name is not None:
                aln.append((name, "".join(buf)))
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        aln.append((name, "".join(buf)))
    return aln


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for n, s in seqs:
            fh.write(f">{n}\n{s}\n")


def run_msa_job(job):
    """One MAFFT alignment (single-threaded) written to disk; runs in a worker pool
    so many small alignments use the core budget concurrently."""
    path, seqs, is_protein = job
    write_fasta(path, mafft(seqs, 1, is_protein=is_protein))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--final", required=True, help="final_telotron_set.tsv")
    ap.add_argument("--focal-ids", default="",
                    help="comma telotron-positive genome_ids (optional if --confident-species given)")
    ap.add_argument("--confident-species", default="",
                    help="TSV from confident_species rule; genome_id column becomes the focal set")
    ap.add_argument("--ortholog-ids", required=True, help="comma genome_ids to search for orthologs")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-array-bp", type=int, default=50)
    ap.add_argument("--min-ident", type=float, default=0.40,
                    help="miniprot ortholog identity floor (raised from prior 0.20: Rost 1999 "
                         "Protein Eng 12:85 twilight zone; memory locus_text_patterns_2026-06-03 "
                         "documents ABS-bucket leakage just under the local-identity gate)")
    ap.add_argument("--intron-tol", type=int, default=8)
    ap.add_argument("--absent-margin", type=int, default=15,
                    help="residues of continuous aligned exon required on BOTH sides of the "
                         "homologous position to call the intron ABSENT (intron-length deletion)")
    ap.add_argument("--min-local-ident", type=float, default=0.40,
                    help="minimum best-per-exon alignment identity among exons within "
                         "+-absent_margin of the homologous position to make a confident "
                         "FILL/ABSENT call; below this the hit is treated as non-homologous "
                         "-> UNCERTAIN. Default 0.40 is documented to leak twilight-zone "
                         "hits into ABS (median ABS lid 0.357 per locus_text_patterns_2026-06-03); "
                         "tighten to 0.50 (NTI median) or 0.70 (deep-dive stringent) and consult "
                         "the stringency_bins.tsv emitted alongside summary.tsv.")
    ap.add_argument("--flank-bp", type=int, default=60, help="exonic DNA flank length to align")
    # ── strict per-flank aa-mismatch gate (primary FILL/ABSENT admission test) ──
    ap.add_argument("--flank-aa", type=int, default=20,
                    help="aa window each side of the homologous junction over which to count "
                         "focal-vs-ortholog protein mismatches (5'-flank = focal residues "
                         "[res-flank_aa,res); 3'-flank = [res,res+flank_aa)).")
    ap.add_argument("--max-flank-mismatch", type=int, default=5,
                    help="max aa mismatches allowed in EACH flank window; a FILL/ABSENT call "
                         "whose 5' or 3' flank exceeds this is downgraded to POOR_FLANK and "
                         "dropped from the clean set / locus_text. Replaces the lenient "
                         "local_ident MAX gate as the primary homology test (local_ident is "
                         "still computed/reported for back-compat).")
    # ── reciprocal-best-hit (RBH) paralog / gene-copy exclusion ──
    ap.add_argument("--rbh", dest="rbh", action="store_true", default=True,
                    help="require the assigned ortholog ORF protein's best reverse miniprot hit "
                         "vs the focal genome to land on the SAME focal gene (default on).")
    ap.add_argument("--no-rbh", dest="rbh", action="store_false",
                    help="disable the reciprocal-best-hit paralog/gene-copy check.")
    ap.add_argument("--max-alt-loci", type=int, default=0,
                    help="max n_alternative_loci (paralog-candidate co-equal hits of the focal "
                         "protein in the ortholog) tolerated for a clean FILL/ABSENT; above this "
                         "the call is downgraded (paralog ambiguity). Default 0.")
    ap.add_argument("--motif", default="TTTAGGG")
    ap.add_argument("--miniprot", default="miniprot")
    args = ap.parse_args()
    ensure_tool_on_path("miniprot", "mafft")

    import csv as _csv0
    focal_ids = set(x for x in args.focal_ids.split(",") if x)
    if args.confident_species and os.path.exists(args.confident_species):
        for r in _csv0.DictReader(open(args.confident_species), delimiter="\t"):
            focal_ids.add(r["genome_id"])
    focal_ids = sorted(focal_ids)
    if not focal_ids:
        raise SystemExit("no focal_ids (use --focal-ids and/or --confident-species)")
    print(f"telotron_ortholog_align: {len(focal_ids)} focal species", flush=True)
    ortho_ids = [x for x in args.ortholog_ids.split(",") if x]
    os.makedirs(args.outdir, exist_ok=True)
    pmsa_dir = os.path.join(args.outdir, "protein_msa")
    idna_dir = os.path.join(args.outdir, "intron_dna")
    flank_dir = os.path.join(args.outdir, "flank_dna")
    os.makedirs(pmsa_dir, exist_ok=True)
    os.makedirs(idna_dir, exist_ok=True)
    os.makedirs(flank_dir, exist_ok=True)
    work = os.path.join(args.outdir, "_work")
    os.makedirs(work, exist_ok=True)

    # ── read telotrons for the focal genomes ────────────────────────────────
    # Read by column NAME (csv.DictReader) rather than by position — the header
    # gains an extra `tx_ids`/`gene_ids` column after collapse_unique_loci, and
    # positional indexing (previously f[8]/f[9]) misaligns silently on collapsed
    # rows. Also prefer the comma-joined `tx_ids`/`gene_ids` when the singular
    # column is empty (collapsed rows only carry the multi-id form).
    # NB: motif is honoured per-row (not the CLI default) so genomes with non-TTTAGGG
    # canonical telomeres get the right tandem reference for telomere_frac.
    tel_by_focal = {}
    import csv as _csv
    with open(args.final) as _fh:
        for r in _csv.DictReader(_fh, delimiter="\t"):
            if r.get("genome_id") not in focal_ids:
                continue
            try:
                s, e = int(r["start"]), int(r["end"])
            except (KeyError, ValueError):
                continue
            if e - s < args.min_array_bp:
                continue
            row_motif = (r.get("motif") or args.motif).strip() or args.motif
            # Prefer singular column; fall back to first of the comma-joined form.
            tx = r.get("tx_id") or (r.get("tx_ids") or "").split(",")[0]
            gene = r.get("gene_id") or (r.get("gene_ids") or "").split(",")[0]
            tel_by_focal.setdefault(r["genome_id"], []).append(
                dict(seqid=r["seqid"], start=s, end=e, tx=tx, gene=gene, motif=row_motif))

    # ── per focal genome: build host proteins + telotron residue positions ──
    # focal protein record: prot_id -> dict(focal, gene, tx, seq, telotrons=[(tid,res,dna)])
    prot_rec = {}
    telotrons = []   # flat list: dict(tid, focal, prot, gene, res, dna, motif)
    for fid in focal_ids:
        fa, gff = find_genome_files(fid, args.refseq_dir, args.tara_dir)
        if not fa or not gff:
            print(f"[skip focal {fid}] missing fasta/gff", file=sys.stderr)
            continue
        contigs = dict(read_fasta(fa))
        cds = load_cds(gff)
        for t in tel_by_focal.get(fid, []):
            d = cds.get(t["tx"])
            if not d or not d["prot"] or t["seqid"] not in contigs:
                continue
            rp = intron_residue(d, t["start"], t["end"])
            if rp is None:
                continue
            res, phase = rp
            pid = d["prot"]
            if pid not in prot_rec:
                # focal CDS genomic span (for reciprocal-best-hit overlap test):
                # the contig + [min,max] of this transcript's CDS segments.
                _segs = d["segs"]
                prot_rec[pid] = dict(focal=fid, gene=t["gene"], tx=t["tx"],
                                     seq=protein_of(d, contigs[t["seqid"]]),
                                     contig=t["seqid"],
                                     cds_lo=min(s for s, _ in _segs),
                                     cds_hi=max(e for _, e in _segs))
            dna = contigs[t["seqid"]][t["start"] - 1:t["end"]]
            if d["strand"] == "-":
                dna = rc(dna)
            f5, f3 = extract_flanks(contigs[t["seqid"]], t["start"], t["end"], d["strand"], args.flank_bp)
            tid = f"{fid}|{t['seqid']}:{t['start']}-{t['end']}"
            telotrons.append(dict(tid=tid, focal=fid, prot=pid, gene=t["gene"],
                                  res=res, phase=phase, dna=dna, flank5=f5, flank3=f3,
                                  motif=t.get("motif", args.motif)))
    print(f"focal host proteins: {len(prot_rec)}; telotrons mapped: {len(telotrons)}", file=sys.stderr)

    # combined query of focal host proteins
    query_faa = os.path.join(work, "focal_host_proteins.faa")
    write_fasta(query_faa, [(p, r["seq"]) for p, r in prot_rec.items()])
    # SHA1 of the query is the miniprot cache key. Any change to the focal panel
    # (genome added/removed, telotron set re-derived, protein-id collision fixed)
    # changes the digest, so stale GFFs from a prior query are NOT silently reused.
    with open(query_faa, "rb") as _qf:
        query_sha = hashlib.sha1(_qf.read()).hexdigest()[:12]

    # ── miniprot focal host proteins vs each ortholog genome ────────────────
    # ortho[oid][prot] = parsed alignment
    ortho = {}
    for oid in ortho_ids:
        fa, gff = find_genome_files(oid, args.refseq_dir, args.tara_dir)
        if not fa:
            print(f"[skip ortholog {oid}] missing fasta", file=sys.stderr)
            continue
        out_gff = os.path.join(work, f"{oid}.{query_sha}.miniprot.gff")
        if not (os.path.exists(out_gff) and os.path.getsize(out_gff) > 0):
            run_miniprot(args.miniprot, fa, query_faa, args.threads, out_gff)  # else reuse cache
        ortho[oid] = dict(aln=parse_miniprot(out_gff), contigs=dict(read_fasta(fa)))
        print(f"  miniprot vs {oid}: {len(ortho[oid]['aln'])} proteins aligned", file=sys.stderr)

    # ── per host protein: protein MSA (focal + orthologs, introns removed) ──
    def ortho_protein(oid, pid):
        """Naive translation of the ortholog ORF: CDS exons concatenated and read in
        frame 0. Frameshifts are NOT corrected -- they surface as internal '*' stops
        and as the fs=/stop= counts in the header, so a mis-aligned intron is visible."""
        a = ortho[oid]["aln"].get(pid)
        if not a or a["ident"] < args.min_ident:
            return None, None
        seq = "".join(ortho[oid]["contigs"].get(a["contig"], "")[s - 1:e] for s, e in a["segs"])
        if a["strand"] == "-":
            seq = rc(seq)
        return translate(seq).rstrip("*"), a   # internal stops kept (frameshift/pseudogene damage shown)

    # ── reciprocal-best-hit (RBH) paralog/gene-copy exclusion ────────────────
    # For each ortholog locus assigned to a focal protein, take that ortholog ORF protein
    # and miniprot it BACK against the FOCAL genome; require its best (Rank=1, max-coverage)
    # hit to land on the SAME focal gene (same contig, overlapping the focal CDS span). If
    # the reverse best hit falls on a different focal locus, the assigned ortholog is a
    # paralog/gene copy and the call FAILS RBH (Wolf & Koonin 2012; Tatusov 1997 COG RBH).
    # Batched: one reverse miniprot per (ortholog genome, focal genome) pair, cached on disk.
    # rbh_ok[(oid, pid)] = bool. Missing key (e.g. --no-rbh) is treated as pass by callers.
    rbh_ok = {}
    if args.rbh:
        # focal genome FASTAs (the focal id is also in the ortholog panel as a cross-Eimeria
        # target, but we need the path explicitly for the reverse search).
        focal_fa = {}
        for fid in focal_ids:
            ffa, _ = find_genome_files(fid, args.refseq_dir, args.tara_dir)
            if ffa:
                focal_fa[fid] = ffa
        # group assigned ortholog ORF proteins by (oid, focal genome of pid)
        rev_groups = {}   # (oid, fid) -> list of (pid, ortho_orf_protein)
        for oid in ortho_ids:
            for pid, rr in prot_rec.items():
                fid = rr["focal"]
                if oid == fid:
                    continue                       # a genome vs itself is trivially reciprocal
                op, a = ortho_protein(oid, pid)
                if not op:
                    continue
                rev_groups.setdefault((oid, fid), []).append((pid, op))
        for (oid, fid), items in rev_groups.items():
            ffa = focal_fa.get(fid)
            if not ffa:
                continue
            rev_faa = os.path.join(work, f"rev.{oid}__{fid}.{query_sha}.faa")
            # query id = pid so the reverse parse keys back to the focal protein it came from
            write_fasta(rev_faa, items)
            rev_gff = os.path.join(work, f"rev.{oid}__{fid}.{query_sha}.miniprot.gff")
            if not (os.path.exists(rev_gff) and os.path.getsize(rev_gff) > 0):
                run_miniprot(args.miniprot, ffa, rev_faa, args.threads, rev_gff)
            rev_best = parse_miniprot(rev_gff)
            for pid, _op in items:
                rr = prot_rec[pid]
                rb = rev_best.get(pid)
                # RBH passes iff the reverse best hit lands on the focal protein's own
                # contig and overlaps its CDS span (same gene); else paralog -> fail.
                ok = False
                if rb and rb["contig"] == rr["contig"]:
                    rlo = min(s for s, _ in rb["segs"])
                    rhi = max(e for _, e in rb["segs"])
                    ok = (rlo <= rr["cds_hi"] and rhi >= rr["cds_lo"])  # span overlap
                rbh_ok[(oid, pid)] = ok
        print(f"  RBH: {sum(rbh_ok.values())}/{len(rbh_ok)} ortholog ORF proteins reciprocal",
              file=sys.stderr)

    # Cache of flank-mismatch results so telotrons sharing a host protein don't
    # realign the same focal/ortholog protein pair repeatedly.
    #
    # The key MUST include `res`. It used to be just (pid, oid), on the premise
    # that res is "a property of the host gene" — it is not: res is the residue
    # position of WHICHEVER INTRON is telotronic, and one transcript can carry
    # several telotrons. A gene with telotrons in intron 2 (res=40) and intron 7
    # (res=520) had the second one silently reuse the first's gate result,
    # measured 480 residues away. If residue 520 sits in a poorly-conserved
    # region whose true counts are (14, 11), it still inherited (1, 0) ->
    # flank_ok=True -> a FILL or ABSENT call the strict gate exists to reject
    # as POOR_FLANK.
    _flank_cache = {}

    msa_jobs = []   # (out_path, seqs, is_protein) -- aligned concurrently in a pool below
    for pid, r in prot_rec.items():
        seqs = [(f"{r['focal']}|{pid}|FOCAL|fs=0|stop=0", r["seq"])]
        for oid in ortho_ids:
            if oid == r["focal"]:
                continue
            op, a = ortho_protein(oid, pid)
            if op:
                seqs.append((f"{oid}|id={a['ident']:.2f}|fs={a['n_frameshift']}|stop={a['n_stop']}", op))
        if len(seqs) >= 2:
            msa_jobs.append((os.path.join(pmsa_dir, f"{r['focal']}__{pid}.aln.fasta"), seqs, True))

    # ── per telotron: does each ortholog have an intron at res? extract it ──
    map_rows = []
    for t in telotrons:
        intron_dna_seqs = [(f"{t['focal']}|TELOTRON|len={len(t['dna'])}", t["dna"])]
        flank5_seqs = [(f"{t['focal']}|FOCAL|5p", t["flank5"])]
        flank3_seqs = [(f"{t['focal']}|FOCAL|3p", t["flank3"])]
        for oid in ortho_ids:
            if oid == t["focal"]:
                continue
            a = ortho.get(oid, {}).get("aln", {}).get(t["prot"])
            present = covered = fs_at_intron = absent = False
            ident = oint_len = oint_telo = lid = 0.0
            n_fs = n_stop = n_alt = 0
            flank5_mm = flank3_mm = args.flank_aa     # default = full mismatch (fails gate)
            rbh_pass = (not args.rbh)                 # vacuously pass when RBH disabled
            flank_ok = gate_ok = False
            struct_cand = []                          # structural intron-at-junction (pre-gate)
            struct_absent = False                     # structural intron-deletion (pre-gate)
            oint_seq = None
            if a and a["ident"] >= args.min_ident:
                ident = a["ident"]
                n_fs, n_stop = a["n_frameshift"], a["n_stop"]
                n_alt = a.get("n_alternative_loci", 0)
                covered = a["qmin"] - args.intron_tol <= t["res"] <= a["qmax"] + args.intron_tol
                # frameshift / premature stop sitting AT the homologous junction =>
                # the intron-boundary alignment is bad; flag it, do NOT correct it.
                fs_at_intron = any(qs - args.intron_tol <= t["res"] <= qe + args.intron_tol
                                   for qs, qe in (a["fs_spans"] + a["stop_spans"]))
                # LOCAL flank-alignment identity around the homologous position (kept for
                # back-compat reporting only; no longer the primary admission test).
                lid = local_ident(a, t["res"], args.absent_margin)
                # STRICT PER-FLANK aa-MISMATCH GATE (primary FILL/ABSENT admission test).
                # Global-align focal vs ortholog ORF protein and count mismatches in the
                # 5' and 3' flank windows SEPARATELY; require BOTH <= max_flank_mismatch.
                _fk = (t["prot"], oid, t["res"])   # res is part of the value; see cache note
                if _fk not in _flank_cache:
                    op, _a = ortho_protein(oid, t["prot"])
                    if op:
                        _flank_cache[_fk] = flank_mismatches(
                            prot_rec[t["prot"]]["seq"], op, t["res"], args.flank_aa)
                    else:
                        _flank_cache[_fk] = (args.flank_aa, args.flank_aa, 0, 0)
                flank5_mm, flank3_mm, _f5n, _f3n = _flank_cache[_fk]
                flank_ok = (flank5_mm <= args.max_flank_mismatch and
                            flank3_mm <= args.max_flank_mismatch)
                # RBH + paralog-ambiguity (alt-loci) exclusion.
                rbh_pass = (not args.rbh) or rbh_ok.get((oid, t["prot"]), False)
                alt_ok = n_alt <= args.max_alt_loci
                # combined gate REPLACES the lenient lid MAX gate as the FILL/ABSENT test.
                gate_ok = flank_ok and rbh_pass and alt_ok
                # structural intron-at-junction candidates, computed independently of the
                # gate so a structurally-clear FILL/ABSENT that the gate rejects can be
                # surfaced as POOR_FLANK (rather than silently demoted to plain UNCERTAIN).
                struct_cand = [(gs, ge) for qpos, gs, ge in a["introns"]
                               if abs(qpos - t["res"]) <= args.intron_tol]
                cand = struct_cand if gate_ok else []
                if cand:
                    present = True
                    gs, ge = max(cand, key=lambda x: x[1] - x[0])
                    octg = ortho[oid]["contigs"].get(a["contig"], "")
                    oint_seq = octg[gs - 1:ge]
                    if a["strand"] == "-":
                        oint_seq = rc(oint_seq)
                    oint_len = len(oint_seq)
                    oint_telo = telomere_frac(oint_seq, t["motif"])
                    # flanking-exon DNA of the orthologous intron (transcript-oriented)
                    of5, of3 = extract_flanks(octg, gs, ge, a["strand"], args.flank_bp)
                    flag5 = ",FS@junc" if fs_at_intron else ""
                    if of5:
                        flank5_seqs.append((f"{oid}|id={ident:.2f}{flag5}", of5))
                    if of3:
                        flank3_seqs.append((f"{oid}|id={ident:.2f}{flag5}", of3))
                elif not struct_cand:
                    # ABSENT = intron-length DELETION (relative to the loci that have it):
                    # residue r sits inside a single continuous exon, with >= absent_margin
                    # residues aligned on BOTH sides (nearest intron / alignment edge is
                    # >= absent_margin away) i.e. the orthologue encodes straight THROUGH the
                    # homologous position with no intron-sized gap. Also requires the contiguity
                    # window to be frameshift/stop-free (a frameshift there means miniprot may be
                    # stitching through a disrupted intron, so contiguity can't be trusted). The
                    # strict per-flank gate (gate_ok) is the homology test (replaces lid).
                    M = args.absent_margin
                    qp = sorted(q for q, _, _ in a["introns"])
                    left = max([q for q in qp if q <= t["res"]] + [a["qmin"]])
                    right = min([q for q in qp if q >= t["res"]] + [a["qmax"]])
                    clean = not any(qs <= t["res"] + M and qe >= t["res"] - M
                                    for qs, qe in (a["fs_spans"] + a["stop_spans"]))
                    struct_absent = (t["res"] - left >= M and right - t["res"] >= M and clean)
                    absent = struct_absent and gate_ok
                    if absent:
                        # show the orthologue's contiguous flank DNA (the deletion) in the panels
                        octg = ortho[oid]["contigs"].get(a["contig"], "")
                        af5, af3 = absent_flanks(octg, a, t["res"], args.flank_bp)
                        if af5:
                            flank5_seqs.append((f"{oid}|id={ident:.2f}|ABSENT", af5))
                        if af3:
                            flank3_seqs.append((f"{oid}|id={ident:.2f}|ABSENT", af3))
            # POOR_FLANK = a structurally-clear FILL or ABSENT whose strict flank / RBH /
            # alt-loci gate FAILED -> NOT a clean call, excluded from the clean set & locus_text.
            poor_flank = (a is not None and a["ident"] >= args.min_ident and not gate_ok
                          and (struct_cand or struct_absent))
            call = ("FILL" if present else "ABSENT" if absent
                    else "POOR_FLANK" if poor_flank
                    else "UNCERTAIN" if covered else "NOCOV" if a else "NOALIGN")
            map_rows.append(dict(tid=t["tid"], focal=t["focal"], gene=t["gene"],
                                 prot=t["prot"], res=t["res"], ortholog=oid, ident=ident,
                                 call=call, n_frameshift=n_fs, n_stop=n_stop,
                                 fs_at_intron=fs_at_intron, ortho_intron_len=oint_len,
                                 ortho_intron_telofrac=oint_telo, local_ident=lid,
                                 n_alternative_loci=n_alt,
                                 flank5_mm=flank5_mm, flank3_mm=flank3_mm,
                                 rbh_pass=int(rbh_pass)))
            if oint_seq and len(oint_seq) >= 10:
                flag = ",FRAMESHIFT@junction" if fs_at_intron else ""
                intron_dna_seqs.append(
                    (f"{oid}|intron|len={oint_len}|telo={oint_telo:.2f}{flag}", oint_seq))
        safe = t["tid"].replace("|", "__").replace(":", "_")
        if len(intron_dna_seqs) >= 2:
            msa_jobs.append((os.path.join(idna_dir, f"{safe}.aln.fasta"), intron_dna_seqs, False))
        if len(flank5_seqs) >= 2:
            msa_jobs.append((os.path.join(flank_dir, f"{safe}.5p.aln.fasta"), flank5_seqs, False))
        if len(flank3_seqs) >= 2:
            msa_jobs.append((os.path.join(flank_dir, f"{safe}.3p.aln.fasta"), flank3_seqs, False))

    # ── run all MAFFT alignments concurrently (up to --threads cores) ────────
    print(f"aligning {len(msa_jobs)} MSAs across {args.threads} workers...", file=sys.stderr)
    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(run_msa_job, msa_jobs))

    # ── write tables ─────────────────────────────────────────────────────────
    with open(os.path.join(args.outdir, "ortholog_map.tsv"), "w") as fh:
        # NB: new columns (flank5_mm, flank3_mm, rbh_pass) are APPENDED at the end so
        # existing index-based consumers (telotron_ortholog_textdump.load_omap reads cols
        # 0-12) keep working unchanged.
        fh.write("telotron\tfocal\tgene\tprotein\tintron_residue\tortholog\tpct_identity\t"
                 "call\tn_frameshift\tn_stop\tframeshift_at_junction\t"
                 "ortho_intron_len\tortho_intron_telofrac\tlocal_ident\t"
                 "n_alternative_loci\tflank5_mm\tflank3_mm\trbh_pass\n")
        for r in map_rows:
            fh.write(f"{r['tid']}\t{r['focal']}\t{r['gene']}\t{r['prot']}\t{r['res']}\t"
                     f"{r['ortholog']}\t{r['ident']:.3f}\t{r['call']}\t{r['n_frameshift']}\t"
                     f"{r['n_stop']}\t{int(r['fs_at_intron'])}\t"
                     f"{r['ortho_intron_len']}\t{r['ortho_intron_telofrac']:.3f}\t{r['local_ident']:.3f}\t"
                     f"{r.get('n_alternative_loci', 0)}\t"
                     f"{r.get('flank5_mm', '')}\t{r.get('flank3_mm', '')}\t{r.get('rbh_pass', '')}\n")

    # per-telotron rollup
    from collections import defaultdict
    by_t = defaultdict(list)
    for r in map_rows:
        by_t[r["tid"]].append(r)
    with open(os.path.join(args.outdir, "summary.tsv"), "w") as fh:
        # n_poor_flank appended at the end (new strict-flank/RBH downgrade bucket); existing
        # columns/order unchanged for back-compat.
        fh.write("telotron\tfocal\tgene\tn_ortho_aligned\tn_fill\tn_fill_clean\t"
                 "n_fill_frameshift\tn_absent\tn_uncertain\tn_fill_telomeric\t"
                 "median_ortho_intron_telofrac\tn_poor_flank\n")
        n_fill = n_fill_clean = n_fill_fs = n_absent = n_unc = n_poor = 0
        for tid, rs in by_t.items():
            aligned = [r for r in rs if r["call"] in
                       ("FILL", "ABSENT", "UNCERTAIN", "NOCOV", "POOR_FLANK")]
            fill = [r for r in rs if r["call"] == "FILL"]
            fill_clean = [r for r in fill if not r["fs_at_intron"]]
            fill_fs = [r for r in fill if r["fs_at_intron"]]
            absent = [r for r in rs if r["call"] == "ABSENT"]
            # POOR_FLANK rolls into the uncertain bucket (not a clean call)
            unc = [r for r in rs if r["call"] in ("UNCERTAIN", "NOCOV", "POOR_FLANK")]
            poor = [r for r in rs if r["call"] == "POOR_FLANK"]
            telo = [r for r in fill if r["ortho_intron_telofrac"] >= 0.6]
            tfracs = sorted(r["ortho_intron_telofrac"] for r in fill)
            med = tfracs[len(tfracs) // 2] if tfracs else 0.0
            n_fill += len(fill); n_fill_clean += len(fill_clean)
            n_fill_fs += len(fill_fs); n_absent += len(absent); n_unc += len(unc)
            n_poor += len(poor)
            fh.write(f"{tid}\t{rs[0]['focal']}\t{rs[0]['gene']}\t{len(aligned)}\t"
                     f"{len(fill)}\t{len(fill_clean)}\t{len(fill_fs)}\t{len(absent)}\t{len(unc)}\t"
                     f"{len(telo)}\t{med:.3f}\t{len(poor)}\n")
    dec = n_fill + n_absent
    print(f"\nDONE. telotrons={len(by_t)}  ortholog calls={len(map_rows)}", file=sys.stderr)
    if dec:
        print(f"intron PRESENT (FILL) {n_fill} vs ABSENT/deletion {n_absent}  "
              f"(of confidently-callable {dec}); UNCERTAIN {n_unc} (incl POOR_FLANK {n_poor}); "
              f"FILL clean {n_fill_clean} / frameshift@junction {n_fill_fs}", file=sys.stderr)
    print(f"protein MSAs: {len(os.listdir(pmsa_dir))}; intron-DNA alignments: {len(os.listdir(idna_dir))}",
          file=sys.stderr)

    # ── stringency-binned table ─────────────────────────────────────────────
    # FILL/ABSENT calls are sensitive to the --min-local-ident gate
    # (locus_text_patterns_2026-06-03: "110/110 collapses to 17 at lid≥0.70").
    # Re-bucket FILL/ABSENT/UNCERTAIN at multiple lid thresholds and at the
    # two telo_frac thresholds in use across this script (0.6 here) and
    # telotron_ortholog_textdump.py (0.45) so downstream consumers can read
    # the conservative count without re-running.
    lid_bins = [0.40, 0.50, 0.70, 0.85]
    tf_bins = [0.45, 0.60]
    with open(os.path.join(args.outdir, "stringency_bins.tsv"), "w") as fh:
        fh.write("lid_threshold\ttelofrac_threshold\tn_fill\tn_absent\tn_uncertain\t"
                 "n_fill_telomeric\tn_locus_fill\tn_locus_absent\n")
        for lid_t in lid_bins:
            # at this stringency, calls with local_ident < lid_t become UNCERTAIN
            def _recall(r):
                if r["local_ident"] < lid_t and r["call"] in ("FILL", "ABSENT"):
                    return "UNCERTAIN"
                return r["call"]
            for tf_t in tf_bins:
                nfill = nabs = nunc = ntelo = 0
                loc_fill = set()
                loc_abs = set()
                for r in map_rows:
                    c = _recall(r)
                    if c == "FILL":
                        nfill += 1
                        loc_fill.add(r["tid"])
                        if r["ortho_intron_telofrac"] >= tf_t:
                            ntelo += 1
                    elif c == "ABSENT":
                        nabs += 1
                        loc_abs.add(r["tid"])
                    elif c in ("UNCERTAIN", "NOCOV", "POOR_FLANK"):
                        nunc += 1
                fh.write(f"{lid_t:.2f}\t{tf_t:.2f}\t{nfill}\t{nabs}\t{nunc}\t"
                         f"{ntelo}\t{len(loc_fill)}\t{len(loc_abs)}\n")
    print("stringency_bins.tsv written (lid x telofrac sensitivity grid)", file=sys.stderr)


if __name__ == "__main__":
    main()
