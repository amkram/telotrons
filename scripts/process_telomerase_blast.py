#!/usr/bin/env python3
"""Post-process telomerase-vs-genomes BLAST output into a per-query report.

For one genome and one query set (e.g. tert_prot.tsv):

  1. Parse the outfmt-6 HSP rows (the ASCII alignment blocks under each row
     are ignored — we'll do our own alignment).
  2. Cluster HSPs that hit the same subject locus into one entry per
     (qseqid, locus). Two HSPs are "same locus" if they're on the same
     sseqid and their subject spans overlap or sit within --merge-gap bp.
  3. Group the deduplicated entries by qseqid → one section per query.
  4. For each (qseqid, locus): extract a reference window of size
     ±query_length around the best HSP's span; align the FULL query to the
     window (Bio.Align semi-global: free end-gaps on the reference so the
     window can extend past the query); trim the alignment to the
     well-supported region (longest contiguous block whose match density in
     a sliding window stays above --trim-min-frac).

  5. Write a human-readable report and a deduplicated summary TSV.

Output (per (genome, label)):
  <indir>/<label>.report.txt   pretty-printed sections
  <indir>/<label>.dedup.tsv    one row per (qseqid, locus) with metadata
"""
import argparse
import gzip
import os
import sys
import textwrap
from collections import defaultdict

from Bio.Align import PairwiseAligner, substitution_matrices

from _common import rc, find_genome_fasta, read_fasta_desc as load_fasta

LABELS_PROT = {"tert_prot", "other_prot"}


def bare_seqid(sseqid):
    """Strip NCBI pipe-format prefix: 'ref|NW_001.1|' → 'NW_001.1', pass-through bare IDs."""
    if "|" in sseqid:
        parts = sseqid.split("|")
        return parts[1] if len(parts) >= 2 else sseqid
    return sseqid


def parse_blast_tsv(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line[0] in " \t":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                rows.append({
                    "qseqid": f[0], "sseqid": f[1],
                    "pident": float(f[2]), "length": int(f[3]),
                    "mismatch": int(f[4]), "gapopen": int(f[5]),
                    "qstart": int(f[6]), "qend": int(f[7]),
                    "sstart": int(f[8]), "send": int(f[9]),
                    "evalue": float(f[10]), "bitscore": float(f[11]),
                })
            except ValueError:
                continue
    return rows


def cluster_hsps(hsps, merge_gap):
    """Return list of clusters: {qseqid, sseqid, lo, hi, strand, hsps:[...]}.

    HSPs on the same (qseqid, sseqid) whose subject spans overlap or are
    within `merge_gap` bp are merged. Strand is the strand of the best HSP.
    """
    by_q_s = defaultdict(list)
    for h in hsps:
        s_lo = min(h["sstart"], h["send"])
        s_hi = max(h["sstart"], h["send"])
        strand = "+" if h["sstart"] <= h["send"] else "-"
        by_q_s[(h["qseqid"], h["sseqid"], strand)].append(
            {**h, "s_lo": s_lo, "s_hi": s_hi, "strand": strand}
        )
    clusters = []
    for (q, s, strand), hs in by_q_s.items():
        hs_sorted = sorted(hs, key=lambda x: x["s_lo"])
        cur = None
        for h in hs_sorted:
            if cur is None:
                cur = {"qseqid": q, "sseqid": s, "strand": strand,
                       "lo": h["s_lo"], "hi": h["s_hi"], "hsps": [h]}
            elif h["s_lo"] <= cur["hi"] + merge_gap:
                cur["hi"] = max(cur["hi"], h["s_hi"])
                cur["hsps"].append(h)
            else:
                clusters.append(cur)
                cur = {"qseqid": q, "sseqid": s, "strand": strand,
                       "lo": h["s_lo"], "hi": h["s_hi"], "hsps": [h]}
        if cur is not None:
            clusters.append(cur)
    return clusters


def best_hsp(cluster):
    return min(cluster["hsps"], key=lambda h: (h["evalue"], -h["bitscore"]))


def make_aligner(is_protein):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    if is_protein:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
    else:
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -6
        aligner.extend_gap_score = -1
    # Free end-gaps on the TARGET (reference window) — query has to align fully.
    aligner.target_end_open_gap_score = 0
    aligner.target_end_extend_gap_score = 0
    return aligner


def align_pair(aligner, query, ref):
    """Return aligned (q, r) strings. Best alignment only."""
    aln = aligner.align(query, ref)[0]
    return str(aln[0]), str(aln[1])


def match_line(q_aln, r_aln, is_protein):
    out = []
    blosum = substitution_matrices.load("BLOSUM62") if is_protein else None
    for q, r in zip(q_aln, r_aln):
        if q == "-" or r == "-":
            out.append(" ")
        elif q.upper() == r.upper():
            out.append("|")
        elif is_protein and blosum is not None:
            try:
                score = blosum[q.upper(), r.upper()]
            except KeyError:
                score = 0
            out.append("+" if score > 0 else " ")
        else:
            out.append(" ")
    return "".join(out)


def trim_alignment(q_aln, r_aln, match, win=15, min_frac=0.5):
    """Trim leading/trailing positions until a sliding window has ≥ min_frac
    matches (or '+' for protein). Returns trimmed (q, m, r) and (l_trim, r_trim).
    """
    L = len(match)
    if L < win:
        return q_aln, match, r_aln, 0, 0
    # Boolean per-position match (count '|' and '+' as matches).
    good = [1 if c in "|+" else 0 for c in match]
    thresh = int(min_frac * win)
    # left trim
    left = 0
    for i in range(L - win + 1):
        if sum(good[i:i + win]) >= thresh:
            left = i
            break
    else:
        return q_aln, match, r_aln, 0, 0
    # right trim
    right = L
    for j in range(L, win - 1, -1):
        if sum(good[j - win:j]) >= thresh:
            right = j
            break
    return q_aln[left:right], match[left:right], r_aln[left:right], left, L - right


def wrap_alignment(q_aln, m_aln, r_aln, q_label, r_label, width=80):
    """Format alignment block in BLAST-pairwise style."""
    lines = []
    L = len(q_aln)
    qpos = rpos = 0
    label_w = max(len(q_label), len(r_label), 5)
    for i in range(0, L, width):
        qchunk = q_aln[i:i + width]
        mchunk = m_aln[i:i + width]
        rchunk = r_aln[i:i + width]
        qadv = sum(1 for c in qchunk if c != "-")
        radv = sum(1 for c in rchunk if c != "-")
        lines.append(f"{q_label:<{label_w}} {qpos + 1:>6} {qchunk}  {qpos + qadv:>6}")
        lines.append(f"{'':<{label_w}}        {mchunk}")
        lines.append(f"{r_label:<{label_w}} {rpos + 1:>6} {rchunk}  {rpos + radv:>6}")
        lines.append("")
        qpos += qadv
        rpos += radv
    return "\n".join(lines)


def section_for_query(qseqid, qseq, qdesc, clusters, genome_seqs,
                      is_protein, aligner, trim_win, trim_min_frac):
    out = []
    qlen = len(qseq)
    out.append(f"\n{'=' * 80}")
    out.append(f"Query: {qseqid}    len={qlen} {'aa' if is_protein else 'nt'}"
               f"    {qdesc[:120]}")
    out.append(f"{'=' * 80}")
    for i, c in enumerate(("\n".join(textwrap.wrap(qseq, 80))).splitlines()):
        out.append("  " + c)
    out.append(f"\nHits: {len(clusters)} deduplicated locus(es) "
               f"({sum(len(c['hsps']) for c in clusters)} HSPs)\n")
    # Sort clusters by best evalue.
    clusters_sorted = sorted(clusters, key=lambda c: best_hsp(c)["evalue"])
    for idx, c in enumerate(clusters_sorted, 1):
        bh = best_hsp(c)
        sseqid = c["sseqid"]
        strand = c["strand"]
        lo, hi = c["lo"], c["hi"]
        contig = genome_seqs.get(bare_seqid(sseqid))
        if contig is None:
            out.append(f"  Hit {idx}: {sseqid}:{lo}-{hi} ({strand})  "
                       f"[contig not in FASTA — skipping alignment]")
            continue
        # Extended window: ±qlen (nt-space; for protein query treat as ±qlen*3)
        expand = qlen * (3 if is_protein else 1)
        ws = max(0, lo - 1 - expand)
        we = min(len(contig), hi + expand)
        win_seq = contig[ws:we]
        if strand == "-":
            win_seq = rc(win_seq)
        # For protein query the ref window is nucleotide → we can't directly
        # align AA query to nt window with PairwiseAligner. Compromise: 6-frame
        # translate the window and pick the frame that aligns best.
        if is_protein:
            translations = six_frame(win_seq)
            best_aln = None
            best_score = float("-inf")
            best_frame = None
            for frame, prot in translations:
                if len(prot) < 5:
                    continue
                try:
                    aln = aligner.align(qseq, prot)[0]
                except ValueError:
                    continue
                if aln.score > best_score:
                    best_score = aln.score
                    best_aln = (str(aln[0]), str(aln[1]), prot, frame)
            if best_aln is None:
                out.append(f"  Hit {idx}: {sseqid}:{lo}-{hi} ({strand})  "
                           f"[no translatable frame] best HSP "
                           f"pident={bh['pident']:.1f}% e={bh['evalue']:.1e}")
                continue
            q_aln, r_aln, _, frame = best_aln
            r_label = f"{sseqid}_f{frame}"
        else:
            try:
                aln = aligner.align(qseq, win_seq)[0]
            except ValueError:
                out.append(f"  Hit {idx}: skipped (aligner failure)")
                continue
            q_aln, r_aln = str(aln[0]), str(aln[1])
            r_label = sseqid

        m_aln = match_line(q_aln, r_aln, is_protein)
        tq, tm, tr_, lt, rt = trim_alignment(q_aln, r_aln, m_aln,
                                             win=trim_win, min_frac=trim_min_frac)
        out.append(
            f"  --- Hit {idx}/{len(clusters_sorted)} ---  "
            f"Subject: {sseqid}:{lo}-{hi} ({strand})  "
            f"window {ws+1}-{we}  ({we - ws} {'nt'})\n"
            f"  Best HSP: pident={bh['pident']:.1f}%  length={bh['length']}  "
            f"evalue={bh['evalue']:.1e}  bitscore={bh['bitscore']:.0f}  "
            f"qcov={(bh['qend'] - bh['qstart'] + 1)}/{qlen}\n"
            f"  Cluster: {len(c['hsps'])} HSP(s)  trimmed {lt}+{rt} of {len(q_aln)} positions"
        )
        out.append("")
        out.append(wrap_alignment(tq, tm, tr_, "Query", r_label, width=80))
    return "\n".join(out)


def six_frame(seq):
    """Yield (frame, protein) for all 6 frames. frame ∈ {+1,+2,+3,-1,-2,-3}."""
    out = []
    for sign, s in ((1, seq), (-1, rc(seq))):
        for off in range(3):
            sub = s[off:]
            sub = sub[: len(sub) - (len(sub) % 3)]
            prot = translate_dna(sub)
            out.append((sign * (off + 1), prot))
    return out


_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(seq):
    return "".join(_CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))


def write_summary_tsv(out_path, label, clusters_by_q):
    with open(out_path, "w") as fh:
        fh.write("qseqid\tsseqid\tstrand\tlo\thi\tn_hsps\tbest_pident\t"
                 "best_length\tbest_evalue\tbest_bitscore\tqcov\n")
        for q, clusters in clusters_by_q.items():
            for c in sorted(clusters, key=lambda c: best_hsp(c)["evalue"]):
                bh = best_hsp(c)
                qcov = bh["qend"] - bh["qstart"] + 1
                fh.write(f"{q}\t{c['sseqid']}\t{c['strand']}\t{c['lo']}\t{c['hi']}\t"
                         f"{len(c['hsps'])}\t{bh['pident']:.2f}\t{bh['length']}\t"
                         f"{bh['evalue']:.2e}\t{bh['bitscore']:.1f}\t{qcov}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-id", required=True)
    ap.add_argument("--blast-dir", required=True,
                    help="work/results/blast_telomerase_vs_genomes/{gid}/")
    ap.add_argument("--db-dir", required=True,
                    help="work/results/telomerase_db/")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--merge-gap", type=int, default=200,
                    help="HSPs within this many bp on the same subject = same locus")
    ap.add_argument("--trim-window", type=int, default=15)
    ap.add_argument("--trim-min-frac", type=float, default=0.4)
    args = ap.parse_args()

    gid = args.genome_id
    fa_src = find_genome_fasta(gid, args.refseq_dir, args.tara_dir, required=False)
    if not fa_src:
        print(f"no FASTA for {gid}; skipping", flush=True)
        return 0
    if fa_src.endswith(".gz"):
        # decompress in-place (one-shot) so load_fasta can read
        plain = fa_src[:-3]
        if not os.path.exists(plain):
            with gzip.open(fa_src, "rb") as ih, open(plain, "wb") as oh:
                import shutil
                shutil.copyfileobj(ih, oh)
        fa_src = plain
    print(f"[{gid}] loading FASTA {fa_src}", flush=True)
    genome_seqs, _ = load_fasta(fa_src)

    labels = [("tr_nucl", "tr.nucl.fa"),
              ("tert_nucl", "tert.nucl.fa"),
              ("tert_prot", "tert.prot.fa"),
              ("other_nucl", "other.nucl.fa"),
              ("other_prot", "other.prot.fa")]
    for label, qfa_name in labels:
        tsv = os.path.join(args.blast_dir, f"{label}.tsv")
        if not os.path.exists(tsv):
            continue
        qfa = os.path.join(args.db_dir, qfa_name)
        q_seqs, q_descs = load_fasta(qfa)
        is_protein = label in LABELS_PROT
        aligner = make_aligner(is_protein)

        hsps = parse_blast_tsv(tsv)
        if not hsps:
            print(f"  {label}: no hits", flush=True)
            continue
        clusters = cluster_hsps(hsps, args.merge_gap)
        # Group by qseqid for output sections.
        by_q = defaultdict(list)
        for c in clusters:
            by_q[c["qseqid"]].append(c)
        # Sort queries by their best evalue across hits.
        q_order = sorted(by_q.keys(),
                         key=lambda q: min(best_hsp(c)["evalue"] for c in by_q[q]))

        report_path = os.path.join(args.blast_dir, f"{label}.report.txt")
        with open(report_path, "w") as out:
            out.write(f"Telomerase BLAST report — genome {gid} — {label}\n")
            out.write(f"  {len(hsps)} HSPs → {len(clusters)} deduplicated loci "
                      f"across {len(by_q)} queries with at least one hit\n")
            for q in q_order:
                if q not in q_seqs:
                    continue
                out.write(section_for_query(
                    q, q_seqs[q], q_descs.get(q, ""),
                    by_q[q], genome_seqs, is_protein, aligner,
                    args.trim_window, args.trim_min_frac))
        print(f"  {label}: {len(by_q)} queries, {len(clusters)} loci → "
              f"{report_path}", flush=True)
        write_summary_tsv(os.path.join(args.blast_dir, f"{label}.dedup.tsv"),
                          label, by_q)


if __name__ == "__main__":
    sys.exit(main())
