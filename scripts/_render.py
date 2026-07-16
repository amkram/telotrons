"""Shared helpers for the HTML preview scripts.

`render_candidates_html.py` and `render_confident_report.py` both need to
colour telomere-motif hits (exact rotations + 1-mismatch, forward and
reverse-complement) inside flanked intron sequence with a boundary marker
between flank and intron. This module owns the shared code so the two
scripts differ only in their per-row layout and header chrome.
"""
from __future__ import annotations
import gzip
import html


RC_TAB = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(RC_TAB)[::-1]


def rotations(m: str) -> set:
    return {m[i:] + m[:i] for i in range(len(m))}


def load_fasta(path: str) -> dict:
    """Load a FASTA (plain or .gz) into {seqid: uppercase_seq}."""
    opener = gzip.open if path.endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf).upper()
    return seqs


def build_kmer_sets(motif: str):
    """Return (fwd_exact, fwd_1mm, rev_exact, rev_1mm) sets of kmers of length len(motif)."""
    L = len(motif)
    fwd = rotations(motif)
    rev = rotations(rc(motif))

    def mismatches(kms):
        out = set()
        for km in kms:
            for i in range(L):
                for b in "ACGT":
                    if b != km[i]:
                        out.add(km[:i] + b + km[i + 1:])
        return out

    fwd_1mm = mismatches(fwd) - fwd - rev
    rev_1mm = mismatches(rev) - rev - fwd
    return fwd, fwd_1mm, rev, rev_1mm


# Position-class rank — higher wins when two candidate rotations overlap
# on the same base (exact beats 1-mismatch; G-orientation beats C).
CLS_RANK = {"ge": 5, "ce": 4, "gm": 3, "cm": 2}


def classify_positions(seq: str, motif: str):
    """Per-position class: 'ge' (G exact), 'ce' (C exact), 'gm' (G 1mm), 'cm' (C 1mm), or None."""
    n = len(seq)
    L = len(motif)
    if L == 0 or n < L:
        return [None] * n
    fwd_e, fwd_m, rev_e, rev_m = build_kmer_sets(motif)
    cls = [None] * n
    for i in range(n - L + 1):
        km = seq[i:i + L]
        if km in fwd_e:
            new = "ge"
        elif km in rev_e:
            new = "ce"
        elif km in fwd_m:
            new = "gm"
        elif km in rev_m:
            new = "cm"
        else:
            continue
        for j in range(i, i + L):
            if cls[j] is None or CLS_RANK[new] > CLS_RANK[cls[j]]:
                cls[j] = new
    return cls


def render_seq(seq: str, motif: str, up_len: int, dn_len: int) -> str:
    """HTML colored spans with intron-boundary markers ('⎢')."""
    n = len(seq)
    intron_start = up_len
    intron_end = n - dn_len
    cls = classify_positions(seq, motif)
    eff = []
    for i in range(n):
        c = cls[i]
        if c is None:
            c = "fl" if (i < intron_start or i >= intron_end) else "in"
        eff.append(c)
    parts = []
    buf, cur = [], None
    for i in range(n):
        if i == intron_start:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        if i == intron_end and intron_end > intron_start:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        c = eff[i]
        if c != cur:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf = []
            cur = c
        buf.append(html.escape(seq[i]))
    if buf:
        parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
    return "".join(parts)
