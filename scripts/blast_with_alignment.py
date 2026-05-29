#!/usr/bin/env python3
"""Run blastn with pairwise output, filter self-matches, emit TSV + alignment.

Output for each non-self HSP:
    <tab-separated fields>\n
    <alignment block (Score/Identities/Strand + Query/midline/Sbjct)>\n
    \n

Columns: qseqid sseqid pident length mismatch gapopen qstart qend
         sstart send evalue bitscore qcovs qlen
"""
import argparse
import re
import subprocess
import sys

HEADER = "\t".join([
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qlen",
])


# ── self-match filter (same logic as filter_blast_selfmatches.py) ─────────────

def parse_query_locus(qid):
    parts = qid.split("|")
    if len(parts) < 3:
        return None
    try:
        start, end = parts[2].split("-")
        return parts[1], int(start), int(end)
    except ValueError:
        return None


def is_self(qid, sid, sstart, send):
    parsed = parse_query_locus(qid)
    if parsed is None:
        return False
    q_seqid, locus_start, locus_end = parsed
    # sid may be bare seqid or genome_id|seqid
    s_seqid = sid.split("|")[-1] if "|" in sid else sid
    if s_seqid != q_seqid:
        return False
    s0, s1 = min(sstart, send), max(sstart, send)
    return not (s1 < locus_start or s0 > locus_end)


# ── pairwise BLAST output parser ──────────────────────────────────────────────

def count_gap_opens(seq):
    """Number of gap-open events (transitions from non-gap to gap)."""
    in_gap = False
    opens = 0
    for ch in seq:
        if ch == "-":
            if not in_gap:
                opens += 1
                in_gap = True
        else:
            in_gap = False
    return opens


def parse_blast_pairwise(stream):
    """Yield dicts with keys: qid, qlen, sid, evalue, bitscore,
    identities, aln_len, gaps, qstart, qend, sstart, send,
    q_segs, mid_segs, s_segs (alignment text segments).
    """
    qid = qlen = sid = None
    sid_full = ""          # full >line text for display
    # per-HSP state
    hsp = None

    def flush_hsp():
        if hsp and hsp.get("q_segs"):
            yield_hsp(hsp)

    pending = []

    def yield_hsp(h):
        pending.append(dict(h))

    lines = iter(stream)
    in_query_id = False
    qid_buf = []

    for line in lines:
        line = line.rstrip("\n")

        # ── query ID (may wrap) ───────────────────────────────────────────────
        if line.startswith("Query="):
            if hsp:
                yield_hsp(hsp)
                hsp = None
            qid_buf = [line[6:].strip()]
            in_query_id = True
            qlen = None
            sid = None
            continue

        if in_query_id:
            if re.match(r"^Length=\d+", line):
                qlen = int(line.split("=")[1])
                in_query_id = False
                qid = " ".join(qid_buf).replace("\n", "")
                # collapse wrapped query id (BLAST wraps at ~60 chars)
                qid = re.sub(r"\s+", "", qid)
                qid_buf = []
            elif line.strip():
                qid_buf.append(line.strip())
            continue

        # ── new subject ───────────────────────────────────────────────────────
        if line.startswith(">"):
            if hsp:
                yield_hsp(hsp)
                hsp = None
            # bare seqid is first word after ">"
            sid = line[1:].split()[0]
            sid_full = line[1:].strip()
            continue

        # ── new HSP ───────────────────────────────────────────────────────────
        m = re.match(r"^ Score = ([\d.e+\-]+) bits \((\d+)\),\s+Expect = (\S+)", line)
        if m:
            if hsp:
                yield_hsp(hsp)
            hsp = {
                "qid": qid, "qlen": qlen, "sid": sid,
                "bitscore": float(m.group(1)), "raw_score": int(m.group(2)),
                "evalue": m.group(3),
                "identities": None, "aln_len": None, "gaps": 0,
                "qstart": None, "qend": None,
                "sstart": None, "send": None,
                "q_segs": [], "mid_segs": [], "s_segs": [],
                "align_lines": [line],
            }
            continue

        if hsp is None:
            continue

        hsp["align_lines"].append(line)

        # ── Identities / Gaps ─────────────────────────────────────────────────
        m = re.match(r"^ Identities = (\d+)/(\d+)", line)
        if m:
            hsp["identities"] = int(m.group(1))
            hsp["aln_len"] = int(m.group(2))
            continue

        m = re.match(r"^ Gaps = (\d+)/(\d+)", line)
        if m:
            hsp["gaps"] = int(m.group(1))
            continue

        # ── alignment rows ────────────────────────────────────────────────────
        m = re.match(r"^Query\s+(\d+)\s+(\S+)\s+(\d+)", line)
        if m:
            s, seq, e = int(m.group(1)), m.group(2), int(m.group(3))
            hsp["q_segs"].append((s, seq, e))
            if hsp["qstart"] is None:
                hsp["qstart"] = s
            hsp["qend"] = e
            continue

        m = re.match(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line)
        if m:
            s, seq, e = int(m.group(1)), m.group(2), int(m.group(3))
            hsp["s_segs"].append((s, seq, e))
            if hsp["sstart"] is None:
                hsp["sstart"] = s
            hsp["send"] = e
            continue

        # capture midlines (line between Query and Sbjct)
        if hsp["q_segs"] and len(hsp["mid_segs"]) < len(hsp["q_segs"]):
            hsp["mid_segs"].append(line)

    if hsp:
        yield_hsp(hsp)

    yield from pending


def format_alignment(hsp):
    """Reconstruct the alignment block from parsed segments."""
    lines = []
    # header lines stored during parse
    for l in hsp.get("align_lines", []):
        if l.strip().startswith("Score"):
            lines.append(l)
        elif l.strip().startswith("Identities"):
            lines.append(l)
        elif l.strip().startswith("Gaps"):
            lines.append(l)
        elif l.strip().startswith("Strand"):
            lines.append(l)
    lines.append("")
    for i, (qs, qseq, qe) in enumerate(hsp["q_segs"]):
        lines.append(f"Query  {qs:<6} {qseq}  {qe}")
        if i < len(hsp["mid_segs"]):
            lines.append(hsp["mid_segs"][i])
        ss, sseq, se = hsp["s_segs"][i]
        lines.append(f"Sbjct  {ss:<6} {sseq}  {se}")
        lines.append("")
    return "\n".join(lines)


def hsp_to_tsv(hsp):
    qid = hsp["qid"]
    sid = hsp["sid"]
    aln_len = hsp["aln_len"] or 0
    identities = hsp["identities"] or 0
    gaps = hsp["gaps"]
    pident = round(100.0 * identities / aln_len, 3) if aln_len else 0
    mismatch = aln_len - identities - gaps
    # count gap opens from alignment sequences
    q_full = "".join(s for _, s, _ in hsp["q_segs"])
    s_full = "".join(s for _, s, _ in hsp["s_segs"])
    gapopen = count_gap_opens(q_full) + count_gap_opens(s_full)
    qstart, qend = hsp["qstart"], hsp["qend"]
    sstart, send = hsp["sstart"], hsp["send"]
    qlen = hsp["qlen"] or 0
    qcovs = round(abs(qend - qstart + 1) * 100 / qlen, 1) if qlen else 0
    return "\t".join(str(x) for x in [
        qid, sid, pident, aln_len, mismatch, gapopen,
        qstart, qend, sstart, send, hsp["evalue"], hsp["bitscore"],
        qcovs, qlen,
    ])


# ── main ──────────────────────────────────────────────────────────────────────

def run_blast(args, blast_args):
    cmd = [
        "blastn",
        "-db", blast_args.db,
        "-query", blast_args.query,
        "-outfmt", "0",
        "-num_threads", str(blast_args.threads),
        "-word_size", str(blast_args.word_size),
        "-evalue", str(blast_args.evalue),
        "-qcov_hsp_perc", str(blast_args.qcov_hsp_perc),
        "-perc_identity", str(blast_args.perc_identity),
    ]
    if blast_args.dust == "no":
        cmd += ["-dust", "no"]
    if not blast_args.soft_masking:
        cmd += ["-soft_masking", "false"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return proc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-db", required=True)
    ap.add_argument("-query", required=True)
    ap.add_argument("-out", default="-")
    ap.add_argument("-num_threads", dest="threads", type=int, default=8)
    ap.add_argument("-word_size", type=int, default=7)
    ap.add_argument("-evalue", type=float, default=1000)
    ap.add_argument("-qcov_hsp_perc", type=float, default=95)
    ap.add_argument("-perc_identity", type=float, default=98)
    ap.add_argument("-dust", default="no")
    ap.add_argument("-soft_masking", action="store_false", default=False)
    ap.add_argument("--no-self-filter", action="store_true")
    args = ap.parse_args()

    out = open(args.out, "w") if args.out != "-" else sys.stdout
    out.write("# " + HEADER + "\n")

    proc = run_blast(None, args)
    kept = dropped = 0

    for hsp in parse_blast_pairwise(proc.stdout):
        if hsp["qstart"] is None:
            continue
        qid, sid = hsp["qid"], hsp["sid"]
        sstart, send = hsp["sstart"], hsp["send"]
        if not args.no_self_filter and is_self(qid, sid, sstart, send):
            dropped += 1
            continue
        kept += 1
        out.write(hsp_to_tsv(hsp) + "\n")
        out.write(format_alignment(hsp) + "\n")

    proc.wait()
    if args.out != "-":
        out.close()
    print(f"kept {kept}  dropped {dropped} self-matches", file=sys.stderr)


if __name__ == "__main__":
    main()
