#!/usr/bin/env python3
"""locus_text-style view for NON-GENIC interstitial telomeric arrays, using
DNA-level orthology (250 bp flank BLAST) instead of protein orthology.

Focal arrays are admitted as ITS only if they meet the field gold standard
(Azzalin/Nergadze/Santagostino/Giulotto): an intrachromosomal array of >=4 telomeric
repeats with <1 mismatch per unit (== identity to a perfect tandem (TTTAGGG)n repeat
>= 6/7, over fwd & revcomp phases). This replaces an earlier ad-hoc coverage filter
that admitted degenerate microsatellite noise.

For each qualifying focal interstitial array (interstitial_arrays.tsv, genome_id == focal):
  1. extract  L = 250 bp upstream | A = array | R = 250 bp downstream
  2. blastn the L and R flanks against a panel of other Eimeria genomes
  3. an ortholog locus = a target contig carrying a colinear L-hit and R-hit on
     the same strand; the region between them is that genome's "array slot"
  4. classify: array_conserved (all orthologs telomeric) / array_mixed (some keep
     array, others empty/diverged — surfaced rather than collapsed) /
     array_absent (slot empty)
  5. MAFFT-align L, slot, R separately and emit  flank | array | flank  rows,
     wrapped in <color>…</color>, grouped into category folders.

Output: work/results/interstitial_orthologs/locus_text/<category>/<focal_locus>.txt
"""
from __future__ import annotations
import argparse, re, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import SPECIES_ABBREV, SPECIES_BY_ACC
from linker_segmentation import find_arrays            # reuse audited array finder

REPO = Path(__file__).resolve().parents[1]
ARRAYS_TSV = REPO / "work/results/interstitial_arrays.tsv"
OUT = REPO / "work/results/interstitial_orthologs"
FLANK = 250
MAXGAP = 1500          # max bp between paired L/R hits (ortholog slot size).
                       # Aligned with the focal-side `telo_block(maxext=1500)` per-side
                       # extension window so the focal and ortholog slot windows are
                       # symmetric (verifier review_interstitial_ITS M finding 2026-06-04).
MIN_HIT_LEN = 120      # blast hit must cover this much of the 250 bp flank
MIN_PIDENT = 75.0
EMPTY_SLOT_MAX = 15    # ortholog slot <= this (flanks ~abut) => clean array-absent (gain/loss);
                       # a larger non-telomeric slot = diverged orthologous position, not a loss
# ── ITS gold standard (Azzalin/Nergadze/Santagostino/Giulotto; e.g. Nergadze 2007,
#    Santagostino 2020): an interstitial telomeric sequence = an intrachromosomal array
#    of >=4 telomeric repeats with <1 mismatch per unit (their blastn query is (TTAGGG)4
#    with short-query auto-adjust off). For Eimeria the unit is TTTAGGG. "<1 mismatch/unit"
#    == identity to a perfect tandem repeat >= 6/7. ──
TELO_MOTIF = "TTTAGGG"
MIN_ITS_UNITS = 4              # >=4 telomeric repeats
MIN_ITS_IDENT = 6.0 / 7.0     # <1 mismatch per unit

FOCAL = "GCF_000499385.1"          # Eimeria necatrix (per _common.SPECIES_BY_ACC)
# Sister-Eimeria panel for orthology calls. Pulled from SPECIES_BY_ACC (the
# single source of truth) rather than re-hardcoded, so future edits to the
# panel map can't drift between scripts.
PANEL = {gid: SPECIES_BY_ACC[gid] for gid in (
    "GCF_000499425.1",            # Eimeria acervulina
    "GCF_000499605.1",            # Eimeria maxima
    "GCF_000499545.2",            # Eimeria tenella
    "GCF_000499745.2",            # Eimeria mitis
)}
SP = {FOCAL: SPECIES_BY_ACC[FOCAL], **PANEL}
ABBR = {SP[gid]: SPECIES_ABBREV[SP[gid]] for gid in (FOCAL, *PANEL)}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fna(gid):
    return next((REPO / "data/raw/refseq" / gid).glob("*_genomic.fna"))


def load_fasta(path):
    d, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                d[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip().upper())
    if name:
        d[name] = "".join(buf)
    return d


def rc(s):
    return s.translate(COMP)[::-1]


def tandem_identity(seq):
    """ITS gold-standard score: best identity of `seq` to a PERFECT tandem (TTTAGGG)n
    repeat over all 7 phases x {fwd, revcomp}. An ITS requires len >= 4 units AND
    identity >= 6/7 (= <1 mismatch per unit). This is the Giulotto-lab criterion;
    it cleanly rejects degenerate microsatellite noise that find_arrays' fuzzy
    <=1mm-per-unit matching admitted."""
    s = seq.upper()
    L = len(s)
    if L < 7:
        return 0.0
    best = 0
    for base in (TELO_MOTIF, rc(TELO_MOTIF)):
        for ph in range(len(base)):
            rot = base[ph:] + base[:ph]
            ref = (rot * (L // len(rot) + 2))[:L]
            m = sum(1 for a, b in zip(s, ref) if a == b)
            if m > best:
                best = m
    return best / L


TELO_HEPTAMERS = ({TELO_MOTIF[i:] + TELO_MOTIF[:i] for i in range(len(TELO_MOTIF))} |
                  {rc(TELO_MOTIF)[i:] + rc(TELO_MOTIF)[:i] for i in range(len(TELO_MOTIF))})


def _telo_cov(s):
    """boolean coverage of exact telomeric heptamers over uppercase string s."""
    n = len(s)
    cov = bytearray(n)
    for i in range(n - 6):
        if s[i:i + 7] in TELO_HEPTAMERS:
            cov[i:i + 7] = b"\x01" * 7
    return cov


def telo_block(seq, s, e, gap=24, maxext=1500):
    """Extend the called array [s,e] to the full CONTIGUOUS telomeric block — merge exact
    telomeric heptamers separated by <= `gap` bp (the Giulotto 'merge/cluster adjacent
    hits' step). This pulls adjacent telomeric repeats OUT of the flanks so the 250 bp
    flanks are clean unique DNA that MAFFT can align without low-complexity gaps. Capped
    at +/- maxext so subtelomeric runs don't blow up. Returns block (start, end)."""
    su = seq.upper()
    lo, hi = max(0, s - maxext), min(len(su), e + maxext)
    cov = _telo_cov(su[lo:hi])
    n = len(cov)
    runs, i = [], 0
    while i < n:
        if cov[i]:
            j = i
            while j < n and cov[j]:
                j += 1
            runs.append([i, j]); i = j
        else:
            i += 1
    merged = []
    for r in runs:
        if merged and r[0] - merged[-1][1] <= gap:
            merged[-1][1] = r[1]
        else:
            merged.append(list(r))
    cs, ce, bs, be = s - lo, e - lo, s - lo, e - lo
    for m in merged:
        if m[0] <= ce and m[1] >= cs:          # the merged block overlapping the called array
            bs, be = min(bs, m[0]), max(be, m[1])
    return lo + bs, lo + be


def flank_telo_frac(s):
    """fraction of a flank covered by exact telomeric heptamers (flags subtelomeric loci)."""
    return sum(_telo_cov(s.upper())) / max(1, len(s))


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw)


def mafft(seqs):
    """seqs: dict id->str (may include ''). Return dict id->aligned (all-gap for '')."""
    ne = {k: v for k, v in seqs.items() if v}
    if not ne:
        return {k: "" for k in seqs}
    if len(ne) == 1:
        k = next(iter(ne)); w = len(ne[k])
    else:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
            for k, v in ne.items():
                tf.write(f">{k}\n{v}\n")
            tmp = tf.name
        out = run(["mafft", "--quiet", "--auto", "--thread", "1", tmp]).stdout
        Path(tmp).unlink()
        aln, name, buf = {}, None, []
        for line in out.splitlines():
            if line.startswith(">"):
                if name:
                    aln[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name:
            aln[name] = "".join(buf)
        w = len(next(iter(aln.values())))
        ne = aln
    return {k: ne.get(k, "-" * w) for k in seqs}


def focal_arrays(limit=None):
    rows = []
    with open(ARRAYS_TSV) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[ix["genome_id"]] != FOCAL:
                continue
            rows.append({"seqid": c[ix["seqid"]], "start": int(c[ix["start"]]),
                         "end": int(c[ix["end"]]), "motif": c[ix["motif"]]})
            if limit and len(rows) >= limit:
                break
    return rows


def build_db(workdir, panel_seqs):
    cat = workdir / "panel.fa"
    with cat.open("w") as out:
        for gid, seqs in panel_seqs.items():
            for name, seq in seqs.items():
                out.write(f">{gid}::{name}\n{seq}\n")
    run(["makeblastdb", "-in", str(cat), "-dbtype", "nucl", "-out", str(workdir / "panel")])
    return workdir / "panel"


def blast_flanks(queries, db, workdir):
    qfa = workdir / "flanks.fa"
    with qfa.open("w") as out:
        for qid, seq in queries:
            out.write(f">{qid}\n{seq}\n")
    res = run(["blastn", "-query", str(qfa), "-db", str(db), "-task", "blastn",
               "-evalue", "1e-5", "-max_target_seqs", "50", "-num_threads", "8",
               "-outfmt", "6 qseqid sseqid pident length sstart send evalue bitscore"])
    hits = defaultdict(list)   # qid -> list of hit dicts
    for line in res.stdout.splitlines():
        q, s, pid, ln, ss, se, ev, bs = line.split("\t")
        pid, ln, ss, se, bs = float(pid), int(ln), int(ss), int(se), float(bs)
        if ln < MIN_HIT_LEN or pid < MIN_PIDENT:
            continue
        gid, contig = s.split("::", 1)
        lo, hi, strand = (ss, se, "+") if ss < se else (se, ss, "-")
        hits[q].append({"gid": gid, "contig": contig, "lo": lo - 1, "hi": hi,
                        "strand": strand, "pid": pid, "bs": bs})
    return hits


def best_per_genome(hit_list):
    best = {}
    for h in hit_list or []:
        key = (h["gid"], h["contig"], h["strand"])
        if key not in best or h["bs"] > best[key]["bs"]:
            best[key] = h
    return best


def orthologs(idx, Lhits, Rhits, panel_seqs):
    """Pair L/R hits per (genome,contig,strand) into ortholog slots."""
    out = []
    Lb, Rb = best_per_genome(Lhits), best_per_genome(Rhits)
    for key, L in Lb.items():
        R = Rb.get(key)
        if not R:
            continue
        gid, contig, strand = key
        seq = panel_seqs[gid][contig]
        if strand == "+":
            if not (L["hi"] <= R["lo"] and R["lo"] - L["hi"] <= MAXGAP):
                continue
            Lf, slot, Rf = seq[L["lo"]:L["hi"]], seq[L["hi"]:R["lo"]], seq[R["lo"]:R["hi"]]
        else:
            if not (R["hi"] <= L["lo"] and L["lo"] - R["hi"] <= MAXGAP):
                continue
            Lf, slot, Rf = rc(seq[L["lo"]:L["hi"]]), rc(seq[R["hi"]:L["lo"]]), rc(seq[R["lo"]:R["hi"]])
        # Slot admission must match focal admission: Giulotto gold standard
        # (>=4 telomeric units AND identity to perfect tandem >= 6/7).
        # `find_arrays` alone admits 3-copy/<=1mm noise and inflates
        # array_conserved counts (verifier review_interstitial_ITS Finding 2,
        # FIXES_PRIORITIZED.md #26).
        slot_up = slot.upper()
        has_array = (
            bool(find_arrays(slot_up))
            and len(slot_up) >= MIN_ITS_UNITS * len(TELO_MOTIF)
            and tandem_identity(slot_up) >= MIN_ITS_IDENT
        )
        out.append({"gid": gid, "contig": contig, "strand": strand,
                    "lpid": L["pid"], "rpid": R["pid"], "slot_len": len(slot),
                    "L": Lf, "A": slot, "R": Rf, "has_array": has_array})
    out.sort(key=lambda o: -(o["lpid"] + o["rpid"]))
    return out


def ortho_status(o):
    if o["has_array"]:
        return "array"
    return "empty" if o["slot_len"] <= EMPTY_SLOT_MAX else "diverged"


def categorize(orthos):
    """array_conserved (>=2 orthologs keep a telomeric array, or sole ortholog does) /
    array_mixed (some keep array, others empty/diverged -> mixed evidence; surfaced rather
    than silently collapsed to a single label per verifier review_interstitial_ITS M
    finding 2026-06-04) / array_absent (best ortholog has a near-empty slot -> clean gain
    or loss) / array_diverged (orthologous position holds substantial non-telomeric
    sequence -> NOT a clean loss) / no_ortholog."""
    if not orthos:
        return "no_ortholog"
    n_array = sum(1 for o in orthos if o["has_array"])
    if n_array == len(orthos):
        return "array_conserved"
    if n_array >= 1:
        # quorum: a sole supporting ortholog is still consistent with conservation;
        # 2+ orthologs disagreeing with 1+ array-bearing ortholog is mixed evidence
        if len(orthos) == 1:
            return "array_conserved"
        return "array_mixed"
    best = max(orthos, key=lambda o: o["lpid"] + o["rpid"])
    return "array_absent" if best["slot_len"] <= EMPTY_SLOT_MAX else "array_diverged"


OPEN, CLOSE, IDW = "<color>", "</color>", 20


def ruler(width, pad):
    s = [" "] * width
    for p in range(width):
        if (p + 1) % 10 == 0:
            for j, ch in enumerate(str(p + 1)):
                pos = p - len(str(p + 1)) + 1 + j
                if 0 <= pos < width:
                    s[pos] = ch
    return " " * pad + "".join(s).rstrip()


SUFFIX = {"array": "/ARR", "empty": "/EMPTY", "diverged": "/DIVERGED"}


def emit(focal_rec, orthos, path):
    rows = [("Emit*", focal_rec["L"], focal_rec["A"], focal_rec["R"], True)]
    for o in orthos:
        rows.append((ABBR[SP[o["gid"]]] + SUFFIX[ortho_status(o)],
                     o["L"], o["A"], o["R"], o["has_array"]))
    cat = categorize(orthos)
    Ls = mafft({r[0]: r[1] for r in rows})       # flanks ARE homologous -> align
    Rs = mafft({r[0]: r[3] for r in rows})
    # array column: only cross-align telomeric (homologous) slots; show non-array slots RAW
    telo = {r[0]: r[2] for r in rows if r[4]}
    Aa = mafft(telo) if len(telo) >= 2 else dict(telo)
    raw = {r[0]: r[2] for r in rows if not r[4]}
    wA = max([len(v) for v in list(Aa.values()) + list(raw.values())] + [1])
    As = {r[0]: (Aa.get(r[0], raw.get(r[0], "")) + "-" * (wA - len(Aa.get(r[0], raw.get(r[0], "")))))
          for r in rows}
    # Place the | at the INFERRED INSERTION boundary, not just the array edges. For a
    # gain (absent/diverged) the focal sequence the ortholog lacks can extend past the
    # repeat array into the "flanks"; the insertion = where the ortholog's flanks abut
    # (the edges of its alignment gap), read off the best (top) ortholog.
    GAP = set("-. ")
    Lw, Aw = len(Ls["Emit*"]), len(As["Emit*"])
    if cat in ("array_absent", "array_diverged") and len(rows) > 1:
        oL, oR = Ls[rows[1][0]], Rs[rows[1][0]]       # best ortholog (orthos sorted by flank id)
        left_cut = 1 + max((i for i, c in enumerate(oL) if c not in GAP), default=Lw - 1)
        right_cut = Lw + Aw + next((i for i, c in enumerate(oR) if c not in GAP), 0)
    else:
        left_cut, right_cut = Lw, Lw + Aw             # shared array -> array boundaries
    glued = {}
    for r in rows:
        full = Ls[r[0]] + As[r[0]] + Rs[r[0]]
        glued[r[0]] = full[:left_cut] + "|" + full[left_cut:right_cut] + "|" + full[right_cut:]
    W = max(len(v) for v in glued.values())
    lines = [
        "=" * 78,
        f"interstitial array  {FOCAL}|{focal_rec['seqid']}:{focal_rec['start']}-{focal_rec['end']}"
        f"  array_len={len(focal_rec['A'])}  motif={focal_rec['motif']}",
        f"CATEGORY: {cat}   ({len(orthos)} DNA ortholog slot(s))"
        + ("   [SUBTELOMERIC: flanks still telomeric-rich, alignment unreliable]"
           if focal_rec.get("subtelo") else ""),
        "=" * 78, "",
        "## orthologs ##",
        f"{'genome':6} {'contig':20} str  L%id R%id slot_bp array?",
        f"{'Emit*':6} {focal_rec['seqid']:20} +    100  100  {len(focal_rec['A']):>7} (focal)",
    ]
    for o in orthos:
        lines.append(f"{ABBR[SP[o['gid']]]:6} {o['contig']:20} {o['strand']:<4} "
                     f"{o['lpid']:>4.0f} {o['rpid']:>4.0f} {o['slot_len']:>7} "
                     f"{ortho_status(o)}")
    lines += ["", "### ALIGNED — flank | insertion | flank (MAFFT;  | = inferred insertion "
              "boundary: focal seq the ortholog lacks, may exceed the repeat array) ###", "",
              ruler(W, IDW + 1 + len(OPEN) + 1)]
    for tag, *_ in rows:
        lines.append(f"{tag:<{IDW}.{IDW}} {OPEN} {glued[tag]} {CLOSE}")
    (OUT / "locus_text" / cat).mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    # inferred insertion = focal sequence inside the | boundaries (array + co-inserted DNA)
    record = None
    if cat in ("array_absent", "array_diverged"):
        ins = glued["Emit*"].split("|")[1].replace("-", "")
        if ins:
            record = {"id": f"{focal_rec['seqid']}_{focal_rec['start']}-{focal_rec['end']}",
                      "cat": cat, "seq": ins, "ins_len": len(ins),
                      "array_len": len(focal_rec["A"]), "subtelo": bool(focal_rec.get("subtelo")),
                      "ortho": ABBR[SP[orthos[0]["gid"]]]}
    return cat, record


PANEL_SEQS = {}          # set in main, inherited by forked workers (COW)
WITH_NONE = False


def write_insertion_msa(insertions, fa_path, txt_path):
    """MAFFT-align the inferred insertions and present them locus_text-style
    (UPPER=telomeric array, lower=co-inserted non-array DNA). A-rich insertions
    are reverse-complemented first so every array reads G-rich/TTTAGGG and the
    cores align. FASTA is one line per sequence (no wrapping)."""
    def a_rich(s):                                   # C-strand (CCCTAAA) array?
        up = [c for c in s if c.isupper()]
        return sum(c == "C" for c in up) > sum(c == "G" for c in up)
    flip = {r["id"]: a_rich(r["seq"]) for r in insertions}
    oriented = {r["id"]: (rc(r["seq"]) if flip[r["id"]] else r["seq"]) for r in insertions}
    aln_raw = mafft(oriented)
    aln = {}                                 # MAFFT lowercases -> restore array/co-inserted case
    for k, a in aln_raw.items():
        src, i, out = oriented[k], 0, []
        for c in a:
            if c.isalpha():
                out.append(src[i]); i += 1
            else:
                out.append(c)
        aln[k] = "".join(out)
    with fa_path.open("w") as fh:
        for r in insertions:
            fh.write(f">{r['id']} cat={r['cat']} ins_len={r['ins_len']} array_len={r['array_len']}"
                     f"{' rc' if flip[r['id']] else ''}{' subtelo' if r['subtelo'] else ''}\n"
                     f"{aln[r['id']]}\n")
    W = max(len(v) for v in aln.values())
    rows = sorted(insertions, key=lambda r: (r["cat"] != "array_absent", -r["ins_len"]))

    def label(r):
        return ("ABS" if r["cat"] == "array_absent" else "DIV") + ("*" if r["subtelo"] else " ") \
               + ("r" if flip[r["id"]] else " ") + " " + r["id"]
    labels = [label(r) for r in rows]
    idw = min(40, max(len(x) for x in labels))
    pad = idw + 1 + len(OPEN) + 1
    tens = "".join(str((c // 10) % 10) if (c + 1) % 10 == 0 else " " for c in range(W))
    L = ["=" * 78,
         f"INSERTED-SEQUENCE MSA   (n={len(rows)} inferred ITS insertions, MAFFT)",
         "each row = focal sequence the ortholog lacks (the | insertion span):",
         "  UPPERCASE = telomeric array · lowercase = co-inserted non-array DNA",
         "  ABS = clean array_absent gain · DIV = array_diverged · * = subtelomeric (less certain)",
         "  r = A-rich insertion reverse-complemented so its array reads G-rich",
         "rows sorted by category then insertion length",
         "=" * 78, "", "### ALIGNED — inserted sequences ###", "", " " * pad + tens]
    for r, lab in zip(rows, labels):
        L.append(f"{lab:<{idw}.{idw}} {OPEN} {aln[r['id']]} {CLOSE}")
    txt_path.write_text("\n".join(L) + "\n")


def worker(args):
    rec, Lh, Rh = args
    orthos = orthologs(rec["i"], Lh, Rh, PANEL_SEQS)
    if not orthos and not WITH_NONE:
        return "no_ortholog(skipped)", None
    name = f"{FOCAL}__{rec['seqid']}_{rec['start']}-{rec['end']}.txt"
    cat = categorize(orthos)
    _, record = emit(rec, orthos, OUT / "locus_text" / cat / name)
    return cat, record


def main():
    global PANEL_SEQS, WITH_NONE, MIN_ITS_UNITS
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="process only N focal arrays (testing)")
    ap.add_argument("--with-no-ortholog", action="store_true",
                    help="also write locus files for arrays with zero DNA orthologs")
    ap.add_argument("--threads", type=int, default=16)
    # ITS gold standard (Nergadze 2007 / Santagostino 2020 / Giulotto). Default
    # is the literature value (4); exposing the flag makes the choice auditable
    # and lets the Snakefile pass config.yaml interstitial.min_its_units (single
    # source of truth — threshold map #6 at top of config.yaml). DO NOT relax
    # without updating the paper claim.
    ap.add_argument("--min-its-units", type=int, default=MIN_ITS_UNITS,
                    help="minimum telomeric-repeat units to count as an ITS (Giulotto-lab gold standard)")
    args = ap.parse_args()
    WITH_NONE = args.with_no_ortholog
    MIN_ITS_UNITS = args.min_its_units

    OUT.mkdir(parents=True, exist_ok=True)
    print("loading genomes…", file=sys.stderr)
    focal_seqs = load_fasta(fna(FOCAL))
    PANEL_SEQS = {gid: load_fasta(fna(gid)) for gid in PANEL}
    arrays = focal_arrays(args.limit)
    print(f"{len(arrays)} focal interstitial arrays", file=sys.stderr)

    recs, queries, n_noise, n_subtelo = [], [], 0, 0
    for i, a in enumerate(arrays):
        seq = focal_seqs.get(a["seqid"])
        if not seq:
            continue
        Acore = seq[a["start"]:a["end"]]        # gold-standard test on the CALLED core
        if len(Acore) < MIN_ITS_UNITS * len(TELO_MOTIF) or tandem_identity(Acore) < MIN_ITS_IDENT:
            n_noise += 1                        # fails ITS gold standard (>=4 units, <1mm/unit)
            continue
        bs, be = telo_block(seq, a["start"], a["end"])   # absorb adjacent telomeric repeats
        if bs < FLANK or be + FLANK > len(seq):
            continue
        L, A, R = seq[bs-FLANK:bs], seq[bs:be], seq[be:be+FLANK]
        subtelo = flank_telo_frac(L) > 0.25 or flank_telo_frac(R) > 0.25
        n_subtelo += subtelo
        recs.append({**a, "i": i, "start": bs, "end": be, "L": L, "A": A, "R": R, "subtelo": subtelo})
        queries.append((f"{i}_L", L)); queries.append((f"{i}_R", R))
    print(f"dropped {n_noise} arrays failing ITS gold standard (>=4 TTTAGGG units, <1mm/unit); "
          f"{len(recs)} ITS arrays (telomeric-block-extended, clean unique flanks; "
          f"{n_subtelo} flagged subtelomeric); blasting {len(queries)} flanks…", file=sys.stderr)

    from concurrent.futures import ProcessPoolExecutor
    counts = defaultdict(int)
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        hits = blast_flanks(queries, build_db(wd, PANEL_SEQS), wd)
        for c in ("array_conserved", "array_mixed", "array_absent", "array_diverged", "no_ortholog"):
            (OUT / "locus_text" / c).mkdir(parents=True, exist_ok=True)
        tasks = [(rec, hits.get(f"{rec['i']}_L"), hits.get(f"{rec['i']}_R")) for rec in recs]
        print(f"processing {len(tasks)} loci on {args.threads} workers…", file=sys.stderr)
        insertions = []
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            for cat, rec in ex.map(worker, tasks, chunksize=8):
                counts[cat] += 1
                if rec:
                    insertions.append(rec)
    # FASTA of the inferred insertions (array + co-inserted DNA the ortholog lacks)
    fa = OUT / "inserted_seqs.fasta"
    with fa.open("w") as fh:
        for r in sorted(insertions, key=lambda r: (r["cat"], -r["ins_len"])):
            fh.write(f">{FOCAL}|{r['id']} cat={r['cat']} ins_len={r['ins_len']} "
                     f"array_len={r['array_len']} extra={r['ins_len'] - r['array_len']} "
                     f"ortho={r['ortho']}{' subtelo' if r['subtelo'] else ''}\n{r['seq']}\n")
    print(f"\nwrote {len(insertions)} inferred insertions -> {fa}", file=sys.stderr)
    if insertions:
        write_insertion_msa(insertions, OUT / "inserted_seqs.msa.fasta", OUT / "inserted_seqs.msa.txt")
        print(f"wrote insertion MSA -> {OUT}/inserted_seqs.msa.txt", file=sys.stderr)
    print("\ncategory counts:")
    for k in sorted(counts):
        print(f"  {k:28s} {counts[k]}")


if __name__ == "__main__":
    main()
