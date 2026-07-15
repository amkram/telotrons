#!/usr/bin/env python3
"""Per-telotron-locus TEXT dump of the ortholog alignments.

For each telotron locus, writes <outdir>/locus_text/<locus>.txt with four sections:

  UNALIGNED  DNA : left-flank exon, intron, right-flank exon (raw, per species)
  ALIGNED    DNA : left+right flank concatenated (MAFFT) ; intron (MAFFT)
  UNALIGNED  aa  : left-flank / intron(forced) / right-flank translations
  ALIGNED    aa  : left+right flank  DNA->translate->align (MAFFT)

Reuses the alignment outputs of telotron_ortholog_align.py (flank_dna/*.5p|3p,
intron_dna/*) — degapped for the unaligned sections, re-aligned for the concat
sections. IDs are fixed-width; aligned blocks are Clustal-style (interleaved, 60/line).
"""
import argparse
import concurrent.futures as cf
import os
import re
import subprocess as sp
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import ensure_tool_on_path, find_genome_files, open_maybe_gz  # noqa: E402

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
IDW = 22       # fixed-width id column
WRAP = 60      # residues / line
OPEN, CLOSE = "<color>", "</color>"   # wrap sequence letters for the IUPAC-color extension


def translate(dna, frame=0):
    dna = dna.upper()
    return "".join(CODON.get(dna[i:i + 3], "X") for i in range(frame, len(dna) - 2, 3))


def best_frame(dna):
    """frame (0/1/2) minimising internal stops (excluding a single trailing stop)."""
    best, bf = None, 0
    for f in range(3):
        aa = translate(dna, f).rstrip("*")
        s = aa.count("*")
        if best is None or s < best:
            best, bf = s, f
    return bf


def read_aln(path):
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def degap(s):
    return s.replace("-", "")


def mafft(rows, threads=1, amino=False):
    """rows: list of (id, seq). Returns aligned list of (id, seq) preserving input order."""
    if len(rows) < 2:
        return rows
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        for i, (_, s) in enumerate(rows):
            tf.write(f">{i}\n{s}\n")
        inp = tf.name
    try:
        out = sp.run(["mafft", "--auto", "--quiet", "--thread", str(threads),
                      "--amino" if amino else "--nuc", inp],
                     capture_output=True, text=True, check=True).stdout
    except sp.CalledProcessError:
        os.unlink(inp)
        return rows
    os.unlink(inp)
    aln, idx, buf = {}, None, []
    for line in out.splitlines():
        if line.startswith(">"):
            if idx is not None:
                aln[idx] = "".join(buf)
            idx, buf = int(line[1:]), []
        else:
            buf.append(line.strip())
    if idx is not None:
        aln[idx] = "".join(buf)
    return [(rows[i][0], aln.get(i, "")) for i in range(len(rows))]


CATNAME = {"NTI": "intron present (non-telomeric)",
           "TEL": "telotron present (telomeric intron)",
           "ABS": "intron absent (continuous CDS)",
           "UNC": "uncertain (no/low coverage or frameshift at junction)"}
LOCUS_CAT = {"NTI": "intron_present_nontelo", "TEL": "telotron_present",
             "ABS": "intron_absent", "UNC": "uncertain"}


def categorize(r, thr):
    """Per-ortholog category, taken straight from the pipeline's structural call:
      FILL    -> TEL (telomeric intron) / NTI (non-telomeric intron); UNC if the junction
                 alignment carries a frameshift (intron-position unreliable)
      ABSENT  -> ABS  (intron-length deletion: orthologue encodes straight through)
      else    -> UNC  (UNCERTAIN / NOCOV / NOALIGN — not confidently callable)"""
    if r["call"] == "FILL":                       # intron present (frameshift is a boundary
        return "TEL" if r["telofrac"] >= thr else "NTI"   # note, not absence)
    if r["call"] == "ABSENT":
        return "ABS"
    return "UNC"


def load_omap(path):
    """tid -> list of dict(oid, ident, call, fs, n_frameshift, n_stop, intron_len, telofrac)."""
    omap = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "telotron" or len(f) < 13:
            continue
        omap.setdefault(f[0], []).append(
            dict(oid=f[5], ident=float(f[6]), call=f[7], n_frameshift=int(f[8]),
                 n_stop=int(f[9]), fs=(f[10] == "1"), intron_len=int(float(f[11])),
                 telofrac=float(f[12])))
    return omap


def load_focal_cds(gff):
    """seqid -> list of (start, end, strand, phase) for CDS features."""
    by = {}
    if not gff:
        return by
    for line in open_maybe_gz(gff):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        ph = 0 if f[7] in (".", "") else int(f[7])
        by.setdefault(f[0], []).append((int(f[3]), int(f[4]), f[6], ph))
    return by


def intron_phase(by, seqid, s, e):
    """GFF phase of the downstream (3') exon = intron phase (0/1/2); None if not found."""
    for st, en, strand, ph in by.get(seqid, []):
        if strand == "+" and abs(st - (e + 1)) <= 2:
            return ph
        if strand == "-" and abs(en - (s - 1)) <= 2:
            return ph
    return None


def label(gid, focal_gid, code):
    return gid + ("*" if gid == focal_gid else "/" + code)


def emit_unaligned(out, rows):
    """rows: list of (id, seq). One unwrapped line per row: fixed id + full sequence.
    The sequence (not the id) is wrapped in <color>…</color> for the IUPAC-color extension."""
    for rid, seq in rows:
        if seq:
            out.write(f"{rid:<{IDW}.{IDW}} {OPEN} {seq} {CLOSE}\n")
        else:
            out.write(f"{rid:<{IDW}.{IDW}} (none)\n")
    out.write("\n")


def emit_aligned(out, rows):
    """rows: list of (id, seq) all same length. One unwrapped line per row, with a single
    full-width position ruler on top. Sequences are wrapped in <color>…</color>; the ruler is
    shifted by len(OPEN) so its columns still line up with the colored sequence."""
    if not rows:
        out.write("(no sequences)\n\n"); return
    W = max(len(s) for _, s in rows)
    out.write(" " * (IDW + 1 + len(OPEN) + 1) + _ruler(0, W) + "\n")
    for rid, seq in rows:
        out.write(f"{rid:<{IDW}.{IDW}} {OPEN} {seq} {CLOSE}\n")
    out.write("\n")


def mark_boundary(rows, ungapped_boundary, focal_idx=0):
    """Insert a '|' column at the exon-exon junction (intron boundary). The column is
    located in the focal row after `ungapped_boundary` non-gap chars, then a '|' is
    inserted at that column in every row so it stays vertically aligned."""
    if not rows:
        return rows
    foc = rows[focal_idx][1]
    cnt, col = 0, len(foc)
    for i, ch in enumerate(foc):
        if cnt == ungapped_boundary:
            col = i
            break
        if ch != "-":
            cnt += 1
    return [(rid, s[:col] + "|" + s[col:]) for rid, s in rows]


def _ruler(start, end):
    s = [" "] * (end - start)
    for p in range(start, end):
        if (p + 1) % 10 == 0:
            tag = str(p + 1)
            for j, c in enumerate(tag):
                pos = p - start - len(tag) + 1 + j
                if 0 <= pos < len(s):
                    s[pos] = c
    return "".join(s)


def order_species(rows5, focal_gid):
    """focal first, then orthologs by descending id= in the 5p header."""
    def idof(nm):
        m = re.search(r"id=([0-9.]+)", nm)
        return float(m.group(1)) if m else 1.0
    gids = []
    meta = {}
    for nm, _ in rows5:
        g = nm.split("|")[0]
        meta[g] = nm
        gids.append(g)
    others = sorted([g for g in gids if g != focal_gid], key=lambda g: -idof(meta[g]))
    return ([focal_gid] + others) if focal_gid in gids else gids


def process(tid, gene, residue, ortho_dir, out_dir, omap, thr, ph, panel):
    import io
    from collections import Counter
    safe = tid.replace("|", "__").replace(":", "_")
    p5 = os.path.join(ortho_dir, "flank_dna", f"{safe}.5p.aln.fasta")
    p3 = os.path.join(ortho_dir, "flank_dna", f"{safe}.3p.aln.fasta")
    pin = os.path.join(ortho_dir, "intron_dna", f"{safe}.aln.fasta")
    if not (os.path.exists(p5) and os.path.exists(p3)):
        return None
    focal_gid = tid.split("|")[0]
    a5 = {nm.split("|")[0]: seq for nm, seq in read_aln(p5)}
    a3 = {nm.split("|")[0]: seq for nm, seq in read_aln(p3)}
    aintron = {nm.split("|")[0]: seq for nm, seq in read_aln(pin)} if os.path.exists(pin) else {}

    # ---- per-ortholog category (from ortholog_map.tsv), restricted to the comparison PANEL ----
    rows_map = [r for r in omap.get(tid, []) if r["oid"] in panel]
    cat = {focal_gid: "foc"}
    info = {}
    for r in rows_map:
        c = categorize(r, thr)
        cat[r["oid"]] = c
        info[r["oid"]] = r | {"cat": c}
    counts = Counter(v["cat"] for v in info.values())
    # LOCUS category from the panel:
    #   intron_present_nontelo : a panel ortholog retains a non-telomeric intron (ancestral
    #                            intron, filled by the telomeric array)
    #   telotron_present       : no non-telo intron, but a telomeric intron is shared in the panel
    #   intron_absent          : a panel ortholog shows an intron-length deletion, none retains it
    #   uncertain              : the panel state can't be determined
    #
    # Each per-ortholog vote is WEIGHTED by its protein-level identity ("ident" from
    # ortholog_map.tsv). This replaces the old `any()` cascade, in which a single weak
    # NTI call from one panel ortholog could override multiple high-confidence TEL
    # calls (verify_ortholog_alignment.md new_concern 6 / FIXES_PRIORITIZED #25).
    # Memory `locus_text_patterns_2026-06-03`: "110/110 collapses to 17 at lid>=0.70"
    # is the direct evidence that locus categories were stringency-sensitive.
    w_nti = sum(r["ident"] for r in rows_map
                if r["call"] == "FILL" and r["telofrac"] < thr)
    w_tel = sum(r["ident"] for r in rows_map
                if r["call"] == "FILL" and r["telofrac"] >= thr)
    w_abs = sum(r["ident"] for r in rows_map if r["call"] == "ABSENT")
    if w_nti == 0 and w_tel == 0 and w_abs == 0:
        locus_cat = LOCUS_CAT["UNC"]
    elif w_nti >= w_tel and w_nti >= w_abs:
        locus_cat = LOCUS_CAT["NTI"]
    elif w_tel >= w_abs:
        locus_cat = LOCUS_CAT["TEL"]
    else:
        locus_cat = LOCUS_CAT["ABS"]
    # Preserve the legacy "any() saw at least one" booleans for the per-locus header
    # (so the user can see the underlying vote distribution at a glance).
    nti = w_nti > 0
    tel = w_tel > 0
    absent = w_abs > 0

    # species WITH sequence (focal + panel orthologs that have an intron), grouped by category
    rank = {"foc": 0, "NTI": 1, "TEL": 2, "UNC": 3}
    sp_order = sorted([g for g in a5.keys() if g == focal_gid or g in panel],
                      key=lambda g: (rank.get(cat.get(g, "UNC"), 4),
                                     -info.get(g, {}).get("ident", 1.0)))
    L = {g: degap(a5.get(g, "")) for g in sp_order}
    R = {g: degap(a3.get(g, "")) for g in sp_order}
    I = {g: degap(aintron.get(g, "")) for g in sp_order}
    lab = {g: label(g, focal_gid, cat.get(g, "UNC")) for g in sp_order}

    buf = io.StringIO()
    bar = "=" * 78
    buf.write(f"{bar}\n{gene}    {tid}    intron @ residue {residue}\n")
    buf.write(f"LOCUS CATEGORY: {locus_cat}\n{bar}\n\n")

    # ---- ortholog category table (ALL orthologs, incl. those without sequence) ----
    buf.write("## ortholog categories ##\n")
    buf.write(f"{'ortholog':<22}{'category':<50}{'%id':>5}{'intron_bp':>10}{'telofrac':>9}\n")
    buf.write(f"{focal_gid + '*':<22}{'focal (telotron)':<50}{'':>5}{'':>10}{'':>9}\n")
    for r in rows_map:
        c = info[r["oid"]]["cat"]
        ilen = r["intron_len"] if r["call"] == "FILL" else 0
        tf = f"{r['telofrac']:.2f}" if r["call"] == "FILL" else ""
        cn = CATNAME[c] + (" [fs@junction]" if (r["call"] == "FILL" and r["fs"]) else "")
        buf.write(f"{r['oid'] + '/' + c:<22}{cn:<50}{r['ident']:>5.2f}"
                  f"{(ilen or ''):>10}{tf:>9}\n")
    buf.write(f"\nid suffix: /NTI non-telo intron · /TEL telotron · /ABS intron absent · "
              f"/UNC uncertain  (telofrac>={thr} => TEL)\n")
    buf.write("(focal + FILL show flank|intron|flank; ABS show contiguous flanks = the deletion "
              "(no intron); UNCERTAIN / no-coverage carry no sequence)\n\n")

    # ---- 1. UNALIGNED DNA ----
    buf.write("### UNALIGNED — DNA ###\n\n-- left flank (5' exon) --\n")
    emit_unaligned(buf, [(lab[g], L[g]) for g in sp_order])
    buf.write("-- intron --\n")
    emit_unaligned(buf, [(lab[g], I[g]) for g in sp_order if I[g]])
    buf.write("-- right flank (3' exon) --\n")
    emit_unaligned(buf, [(lab[g], R[g]) for g in sp_order])

    # ---- 2. ALIGNED DNA  ('|' = intron boundary / splice junction) ----
    buf.write("### ALIGNED — DNA ###\n\n-- left+right flank concatenated (MAFFT;  | = intron boundary) --\n")
    aln_dna = mafft([(lab[g], L[g] + R[g]) for g in sp_order if (L[g] or R[g])], amino=False)
    emit_aligned(buf, mark_boundary(aln_dna, len(L[focal_gid])))
    buf.write("-- intron (MAFFT) --\n")
    # re-align the intron over THIS panel's species only (reusing the full-set alignment and
    # dropping columns would leave all-gap columns where excluded species had insertions)
    emit_aligned(buf, mafft([(lab[g], I[g]) for g in sp_order if I[g]], amino=False))

    # ---- 3. UNALIGNED aa  (flanks only; TRUE intron phase, frame not forced, stops kept) ----
    # Phase falls back to 0 when the focal GFF lacks an unambiguous phase tag for the
    # downstream exon; we now mark this explicitly with phase_assumed so consumers
    # (HTML reviewer, curator's column-conservation count) can stratify on it instead
    # of silently treating phase-0-assumed rows as known-phase
    # (verify_ortholog_alignment.md new_concern 4 / FIXES_PRIORITIZED #64).
    phase_assumed = ph is None
    phx = 0 if phase_assumed else ph
    split_up = (3 - phx) % 3                 # upstream bases in the split junction codon
    note = f"intron phase {ph}" if not phase_assumed else "intron phase unknown (assumed 0; phase_assumed=True)"
    buf.write(f"### UNALIGNED — aa ###\n(flanks translated in the gene reading frame; {note}; "
              f"phase 1/2 splits the junction codon; stops shown, not stripped)\n\n")
    def offL(g):
        return (len(L[g]) - split_up) % 3
    buf.write("-- left flank --\n")
    emit_unaligned(buf, [(lab[g], translate(L[g], offL(g))) for g in sp_order if L[g]])
    buf.write("-- right flank --\n")
    emit_unaligned(buf, [(lab[g], translate(R[g], phx)) for g in sp_order if R[g]])

    # ---- 4. ALIGNED aa  ('|' = intron boundary; codon split when phase != 0) ----
    # Header string is load-bearing: scripts/curate_locus_text.py:FOCAL_BLOCK regex and
    # scripts/wrap_locus_text.py:SECTION constant both key off this exact prefix.
    # Do NOT rename without updating those consumers + their _assert_schema guards.
    buf.write("### ALIGNED — aa flanks (focal frame) + DNA intron ###\n"
              "\n-- left+right flank: DNA->translate->align (MAFFT;  | = intron boundary) --\n")
    aln_aa = mafft([(lab[g], translate(L[g] + R[g], offL(g))) for g in sp_order if (L[g] or R[g])], amino=True)
    emit_aligned(buf, mark_boundary(aln_aa, (len(L[focal_gid]) - offL(focal_gid) + 2) // 3))

    # Second aa block: each species is translated in its own best (stop-minimising) frame.
    # When the orthologous intron sits at a different phase than the focal (inserted vs
    # ancestral), the focal-frame translation above puts the ortholog L-flank in the
    # wrong frame; this block recovers the real protein for visual inspection
    # (verify_ortholog_alignment.md new_concern 4). The consumer (curate_locus_text.py)
    # keys on "(focal frame)" only, so this block is informational and does not feed
    # the conservation count.
    buf.write("\n-- left+right flank: DNA->translate->align (MAFFT; best frame per species) --\n")
    aln_aa_bf = mafft([(lab[g], translate(L[g] + R[g], best_frame(L[g] + R[g])))
                       for g in sp_order if (L[g] or R[g])], amino=True)
    emit_aligned(buf, aln_aa_bf)

    sub = os.path.join(out_dir, locus_cat)
    os.makedirs(sub, exist_ok=True)
    with open(os.path.join(sub, f"{safe}.txt"), "w") as fh:
        fh.write(buf.getvalue())
    return dict(tid=tid, gene=gene, locus_cat=locus_cat,
                NTI=counts.get("NTI", 0), TEL=counts.get("TEL", 0),
                ABS=counts.get("ABS", 0), UNC=counts.get("UNC", 0),
                phase=("" if ph is None else ph), phase_assumed=int(phase_assumed))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ortho-dir", required=True)
    ap.add_argument("--out-dir", default=None, help="default <ortho-dir>/locus_text")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-loci", type=int, default=0)
    ap.add_argument("--telo-frac", type=float, default=0.45,
                    help="ortho-intron telomere-fraction >= this => 'telotron present' (TEL)")
    ap.add_argument("--refseq-dir", default="data/raw/refseq")
    ap.add_argument("--tara-dir", default="data/raw/tara")
    ap.add_argument("--panel-ids", default="",
                    help="comma genome_ids forming the comparison panel; restricts which "
                    "orthologs are shown AND decides the locus category (default: all in the map)")
    args = ap.parse_args()
    all_oids = {r["oid"] for rows in load_omap(os.path.join(args.ortho_dir, "ortholog_map.tsv")).values()
                for r in rows}
    panel = set(x for x in args.panel_ids.split(",") if x) or all_oids
    ensure_tool_on_path("mafft")
    out_dir = args.out_dir or os.path.join(args.ortho_dir, "locus_text")
    import shutil
    for d in set(LOCUS_CAT.values()):                 # clear stale category folders
        shutil.rmtree(os.path.join(out_dir, d), ignore_errors=True)
    os.makedirs(out_dir, exist_ok=True)

    omap = load_omap(os.path.join(args.ortho_dir, "ortholog_map.tsv"))
    loci = {}
    for line in open(os.path.join(args.ortho_dir, "ortholog_map.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] == "telotron" or len(f) < 5:
            continue
        loci.setdefault(f[0], dict(gene=f[2], residue=f[4]))
    items = sorted(loci.items(), key=lambda kv: (kv[1]["gene"], kv[0]))
    if args.max_loci:
        items = items[:args.max_loci]

    # true intron phase per locus from the focal GFF (so aa translation isn't frame-forced)
    cds_cache = {}
    for g in {tid.split("|")[0] for tid in loci}:
        _, gff = find_genome_files(g, args.refseq_dir, args.tara_dir)
        cds_cache[g] = load_focal_cds(gff)
    phase = {}
    for tid in loci:
        g, rest = tid.split("|", 1)
        seqid, coords = rest.split(":")
        s, e = (int(x) for x in coords.split("-"))
        phase[tid] = intron_phase(cds_cache.get(g, {}), seqid, s, e)

    results = []
    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(process, tid, d["gene"], d["residue"], args.ortho_dir, out_dir,
                          omap, args.telo_frac, phase.get(tid), panel)
                for tid, d in items]
        for fu in cf.as_completed(futs):
            r = fu.result()
            if r:
                results.append(r)
                if len(results) % 50 == 0:
                    print(f"  {len(results)} locus text files...", file=sys.stderr)

    # manifest + locus-category counts.
    # phase_assumed=1 marks loci where the focal GFF had no unambiguous phase tag
    # for the downstream exon and we silently defaulted to phase 0 (only correct ~1/3
    # of the time); downstream consumers can stratify on this column rather than
    # treating phase-0-assumed rows as known-phase (FIXES_PRIORITIZED #64).
    from collections import Counter
    with open(os.path.join(out_dir, "locus_categories.tsv"), "w") as fh:
        fh.write("telotron\tgene\tlocus_category\tn_NTI\tn_TEL\tn_ABS\tn_UNC\tphase\tphase_assumed\n")
        for r in sorted(results, key=lambda x: (x["locus_cat"], x["gene"])):
            fh.write(f"{r['tid']}\t{r['gene']}\t{r['locus_cat']}\t{r['NTI']}\t{r['TEL']}\t"
                     f"{r['ABS']}\t{r['UNC']}\t{r['phase']}\t{r['phase_assumed']}\n")
    lc = Counter(r["locus_cat"] for r in results)
    print(f"wrote {len(results)} locus text files to {out_dir}/<category>/", file=sys.stderr)
    print("locus categories:", dict(lc), file=sys.stderr)


if __name__ == "__main__":
    main()
