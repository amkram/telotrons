#!/usr/bin/env python3
"""Per-genome scan: derive introns via GenomeTools, motif hits via seqkit, compose final TSV.

Heavy lifting is delegated:
  - intron coordinates from GFF       : `gt gff3 -retainids -addintrons`
  - telomeric motif occurrences (BED) : `seqkit locate --bed --pattern-file`
Project-specific composition (orientation, splice signals, terminal-motif inference) stays here.
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

from _common import revcomp, rotations, open_maybe_gz, read_fasta, tqdm

LOCI_HEADER = [
    "genome_id", "organism", "group", "source",
    "seqid", "start", "end", "strand",
    "tx_id", "gene_id", "intron_index", "intron_len",
    "contig_len", "distance_to_end",
    "motif", "telomeric_bases", "telomeric_frac", "flank_telomeric_frac",
    "donor", "acceptor", "splice_class",
    "canonical_strand", "splice_canonical",
    "orientation", "fwd_hits", "rev_hits",
    "first40", "last40",
    "terminal_motif", "terminal_motif_bases", "contig_both_ends_capped",
]
SUMMARY_HEADER = [
    "genome_id", "organism", "group", "source",
    "introns_scanned", "telotron_candidates", "contigs",
    "terminal_motif", "terminal_motif_bases", "status",
]


def write_pattern_file(path, motifs):
    """seqkit --pattern-file: each forward rotation gets its own line, all named after the base motif.
    seqkit searches both strands by default, so revcomps are handled implicitly."""
    with open(path, "w") as f:
        for m in motifs:
            m = m.upper()
            for i, rot in enumerate(rotations(m)):
                f.write(f">{m}\n{rot}\n")


def run_gt_introns(gff_path, out_path):
    """Decompress if gz; synthesize exon rows from CDS when no exons (e.g. Tara gmove);
    pipe to `gt gff3 -addintrons`. gt derives introns from exons only."""
    with open_maybe_gz(gff_path) as f:
        body = f.read()
    has_exon = any(
        line.split("\t", 3)[2] == "exon"
        for line in body.splitlines()
        if line and not line.startswith("#") and line.count("\t") >= 7
    )
    if not has_exon:
        extra = []
        for line in body.splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 9 and parts[2] == "CDS":
                parts[2] = "exon"
                extra.append("\t".join(parts))
        if extra:
            body = body + "\n" + "\n".join(extra) + "\n"

    # gt -tidy emits best-effort output then may bail with exit 1 on a late fatal
    # error (e.g. strand inconsistency between paralogs). Don't drop the genome
    # over it — keep whatever was written; load_introns will be empty if nothing did.
    with open(out_path, "w") as out:
        subprocess.run(["gt", "gff3", "-retainids", "-addintrons", "-tidy", "-"],
                       check=False, input=body, text=True,
                       stdout=out, stderr=subprocess.DEVNULL)


def parse_gff_attrs(field):
    out = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def load_introns(gt_gff_path):
    """Read intron rows from gt output; resolve Parent → transcript → gene via the same file."""
    transcripts, gene_of = {}, {}
    intron_rows = []
    with open(gt_gff_path) as f:
        for line in f:
            if not line.strip() or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start, end, _sc, strand, _ph, attr = parts
            attrs = parse_gff_attrs(attr)
            if ftype in ("mRNA", "transcript"):
                tid = attrs.get("ID")
                if tid:
                    transcripts[tid] = strand
                    gene_of[tid] = attrs.get("Parent", "")
            elif ftype == "intron":
                parent = attrs.get("Parent", "").split(",")[0]
                intron_rows.append([seqid, int(start), int(end), strand or transcripts.get(parent, "."),
                                    parent, gene_of.get(parent, "")])
    if not intron_rows:
        return pd.DataFrame(columns=["seqid", "start", "end", "strand", "tx_id", "gene_id",
                                     "intron_index", "intron_len"])
    df = pd.DataFrame(intron_rows, columns=["seqid", "start", "end", "strand", "tx_id", "gene_id"])
    df = df.sort_values(["tx_id", "start"]).reset_index(drop=True)
    df["intron_index"] = df.groupby("tx_id").cumcount() + 1
    df["intron_len"] = df.end - df.start + 1
    return df


def run_seqkit_locate(fasta_path, patterns_path, out_bed):
    """BED6: seqid, start, end, pattern_name (base motif), score, strand.
    -j 1 because parallelism comes from the ProcessPool around us; without this,
    seqkit's default 4 threads × N workers thrashes the box."""
    with open(out_bed, "w") as out:
        subprocess.run(["seqkit", "locate", "-j", "1", "--bed", "--ignore-case",
                        "--pattern-file", patterns_path, fasta_path],
                       check=True, stdout=out, stderr=subprocess.DEVNULL)


def _merged_coverage_per_group(df, group_cols):
    """Within each group (already sorted by start), sum merged interval coverage.
    Vectorized via cumulative-max-end: within a sorted group, new bp contributed by
    each row = max(0, cum_end[i] - max(ov_start[i], cum_end[i-1])). Returns one row
    per group with column `merged_bp`."""
    df = df.copy()
    g = df.groupby(group_cols, sort=False)
    df["_cum_end"] = g["ov_end"].cummax()
    # shift(1) inside a groupby stays within group; first row of each group → NaN.
    df["_prev_end"] = g["_cum_end"].shift(1)
    # For the first row in each group, no prior coverage → fill with ov_start so
    # new_bp = cum_end - ov_start (the full extent of that interval).
    df["_prev_end"] = df["_prev_end"].fillna(df["ov_start"])
    df["_new_bp"] = (df["_cum_end"] - df[["ov_start", "_prev_end"]].max(axis=1)).clip(lower=0)
    return (df.groupby(group_cols, sort=False)["_new_bp"].sum()
              .reset_index().rename(columns={"_new_bp": "merged_bp"}))


def _bedtools_intersect_intron_hits(bed_a_path, bed_b_path, out_path):
    with open(out_path, "w") as out:
        subprocess.run(["bedtools", "intersect", "-a", bed_a_path, "-b", bed_b_path,
                        "-wa", "-wb", "-sorted"],
                       check=True, stdout=out, stderr=subprocess.DEVNULL)


_INT_COLS = ["iseqid", "istart", "iend", "row", "ilen", "istrand",
             "hseqid", "hstart", "hend", "motif", "hscore", "hstrand"]


def _read_intersect(path):
    if not os.path.getsize(path):
        return pd.DataFrame(columns=_INT_COLS)
    return pd.read_csv(path, sep="\t", header=None, names=_INT_COLS)


def _sort_bed_inplace(path):
    subprocess.run(["sort", "-k1,1", "-k2,2n", "-o", path, path], check=True)


def compute_intron_motif_metrics(introns_df, hits_bed_path, contig_lens, tmp_dir):
    """Replace per-intron Python loops with one bedtools intersect + vectorized pandas.

    Returns (best_motif_df, flank_df, strand_df) keyed by `row`. ~30-100× faster than the
    previous per-intron pandas loop on multi-100k-intron genomes.
    """
    # Intron BED, 0-based half-open, with row id as the BED name field.
    introns_bed = os.path.join(tmp_dir, "introns.bed")
    pd.DataFrame({
        "seqid": introns_df.seqid,
        "start0": (introns_df.start - 1).astype(int),
        "end": introns_df.end.astype(int),
        "row": introns_df.row.astype(int),
        "ilen": introns_df.intron_len.astype(int),
        "strand": introns_df.strand,
    }).to_csv(introns_bed, sep="\t", header=False, index=False)
    _sort_bed_inplace(introns_bed)

    hits_sorted = os.path.join(tmp_dir, "hits.sorted.bed")
    subprocess.run(f"sort -k1,1 -k2,2n {hits_bed_path} > {hits_sorted}",
                   shell=True, check=True)

    # === Intron × hit overlaps ===
    inter_path = os.path.join(tmp_dir, "intersect.tsv")
    _bedtools_intersect_intron_hits(introns_bed, hits_sorted, inter_path)
    inter = _read_intersect(inter_path)

    if inter.empty:
        best_motif = pd.DataFrame(columns=["row", "motif", "telomeric_bases", "telomeric_frac"])
        strand = pd.DataFrame(columns=["row", "fwd_hits", "rev_hits"])
    else:
        inter["ov_start"] = inter[["istart", "hstart"]].max(axis=1)
        inter["ov_end"] = inter[["iend", "hend"]].min(axis=1)

        # Per-(row, motif) merged covered bp.
        inter_sorted = inter.sort_values(["row", "motif", "ov_start"])
        cov = _merged_coverage_per_group(inter_sorted, ["row", "motif"])
        cov = cov.merge(introns_df[["row", "intron_len"]], on="row")
        cov["frac"] = cov.merged_bp / cov.intron_len
        best_motif = cov.loc[cov.groupby("row").frac.idxmax(),
                             ["row", "motif", "merged_bp", "frac"]]
        best_motif.columns = ["row", "motif", "telomeric_bases", "telomeric_frac"]

        # Strand-adjusted fwd/rev counts (no merging — raw hit counts).
        fwd_mask = ((inter.istrand != "-") & (inter.hstrand == "+")) | \
                   ((inter.istrand == "-") & (inter.hstrand == "-"))
        strand = inter.assign(fwd=fwd_mask).groupby("row").agg(
            fwd_hits=("fwd", "sum"), total_hits=("fwd", "size")
        ).reset_index()
        strand["rev_hits"] = strand.total_hits - strand.fwd_hits
        strand = strand[["row", "fwd_hits", "rev_hits"]]

    # === Flank ±50 bp, excluding the intron itself ===
    clen = introns_df.seqid.map(contig_lens).fillna(0).astype(int)
    flanks = pd.concat([
        pd.DataFrame({  # left flank
            "seqid": introns_df.seqid,
            "start0": (introns_df.start - 51).clip(lower=0).astype(int),
            "end": (introns_df.start - 1).clip(lower=0).astype(int),
            "row": introns_df.row.astype(int),
        }),
        pd.DataFrame({  # right flank
            "seqid": introns_df.seqid,
            "start0": introns_df.end.astype(int),
            "end": (introns_df.end + 50).clip(upper=clen).astype(int),
            "row": introns_df.row.astype(int),
        }),
    ], ignore_index=True)
    flanks = flanks[flanks.end > flanks.start0].copy()
    flanks["flen"] = flanks.end - flanks.start0

    if flanks.empty:
        flank_df = pd.DataFrame(columns=["row", "flank_telomeric_frac"])
        return best_motif, flank_df, strand

    flanks_bed = os.path.join(tmp_dir, "flanks.bed")
    flanks[["seqid", "start0", "end", "row", "flen"]].assign(strand=".").to_csv(
        flanks_bed, sep="\t", header=False, index=False)
    _sort_bed_inplace(flanks_bed)

    finter_path = os.path.join(tmp_dir, "flank_intersect.tsv")
    _bedtools_intersect_intron_hits(flanks_bed, hits_sorted, finter_path)
    finter = _read_intersect(finter_path)

    if finter.empty:
        flank_df = pd.DataFrame({"row": introns_df.row, "flank_telomeric_frac": 0.0})
    else:
        finter["ov_start"] = finter[["istart", "hstart"]].max(axis=1)
        finter["ov_end"] = finter[["iend", "hend"]].min(axis=1)
        finter_sorted = finter.sort_values(["row", "motif", "ov_start"])
        fcov = _merged_coverage_per_group(finter_sorted, ["row", "motif"])
        flank_len = flanks.groupby("row").flen.sum().reset_index().rename(columns={"flen": "total_flank_len"})
        fcov = fcov.merge(flank_len, on="row")
        fcov["frac"] = fcov.merged_bp / fcov.total_flank_len
        flank_df = fcov.loc[fcov.groupby("row").frac.idxmax(), ["row", "frac"]]
        flank_df.columns = ["row", "flank_telomeric_frac"]

    return best_motif, flank_df, strand


def classify_orientation(fwd, rev):
    if fwd == 0 and rev == 0:
        return "None"
    if fwd >= 2 and fwd > 2 * rev:
        return "Fwd/G-rich"
    if rev >= 2 and rev > 2 * fwd:
        return "Rev/C-rich"
    return "Mixed"


def terminal_motif(hits, contig_lens, window=500, min_frac=0.5):
    """Dominant motif that actually forms a telomeric ARRAY at a contig end.
    A contig end (first or last `window` bp) only votes if it is ≥`min_frac` covered
    by rotations of a single base motif. Fully vectorized — the previous per-contig
    Python loop was 60%+ of total scan time on big genomes."""
    if hits.empty:
        return "", 0
    h = hits.copy()
    h["clen"] = h.seqid.map(contig_lens)
    h = h.dropna(subset=["clen"])
    if h.empty:
        return "", 0
    h["clen"] = h["clen"].astype(int)

    # Left end: hit overlaps [0, window). Right end: hit overlaps [clen-window, clen).
    left = h[h.start < window].copy()
    left["end_lo"] = 0
    left["end_hi"] = window
    left["end_id"] = left.seqid.astype(str) + ":L"

    right = h[h.end > (h.clen - window)].copy()
    right["end_lo"] = (right.clen - window).clip(lower=0).astype(int)
    right["end_hi"] = right.clen
    right["end_id"] = right.seqid.astype(str) + ":R"

    near = pd.concat([left, right], ignore_index=True)
    if near.empty:
        return "", 0

    near["ov_start"] = near[["start", "end_lo"]].max(axis=1)
    near["ov_end"] = near[["end", "end_hi"]].min(axis=1)
    near = near.rename(columns={"name": "motif"})
    near = near.sort_values(["end_id", "motif", "ov_start"])
    merged = _merged_coverage_per_group(near, ["end_id", "motif"])

    # Per-end window length (constant across rows of an end_id).
    win_len = (near.groupby("end_id", sort=False)
                   .agg(end_lo=("end_lo", "first"), end_hi=("end_hi", "first"))
                   .reset_index())
    win_len["window_len"] = (win_len.end_hi - win_len.end_lo).clip(lower=1)
    merged = merged.merge(win_len[["end_id", "window_len"]], on="end_id")
    merged["frac"] = merged.merged_bp / merged.window_len

    qualifying = merged[merged.frac >= min_frac]
    if qualifying.empty:
        return "", 0
    summed = qualifying.groupby("motif")["merged_bp"].sum().sort_values(ascending=False)
    return str(summed.index[0]), int(summed.iloc[0])


def genome_paths(row, refseq_dir, tara_dir):
    gid, _org, group, source = row
    if source == "genoscope" or group == "tara":
        fa_candidates = [p for p in glob.glob(f"{tara_dir}/Contigs/{gid}.fa*")
                         if not p.endswith((".fai", ".gzi"))]
        gff_candidates = [p for p in glob.glob(f"{tara_dir}/Genes/GFF/{gid}.gmove.gff*")
                          if not p.endswith((".tbi", ".idx"))]
        return fa_candidates[0], gff_candidates[0]
    fa_paths = (glob.glob(f"{refseq_dir}/{gid}/*/*/*_genomic.fna*")
                + glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.fna*", recursive=True))
    gff_paths = (glob.glob(f"{refseq_dir}/{gid}/**/genomic.gff*", recursive=True)
                 + glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.gff*", recursive=True))
    fa = next(p for p in fa_paths if not p.endswith((".fai", ".gzi")))
    gff = next(p for p in gff_paths if not p.endswith((".tbi", ".idx")))
    return fa, gff


def scan_one(args):
    row, opts = args
    gid, org, group, source = row
    try:
        fa, gff = genome_paths(row, opts["refseq_dir"], opts["tara_dir"])
        seqs = read_fasta(fa)
        contig_lens = {sid: len(s) for sid, s in seqs.items()}

        with tempfile.TemporaryDirectory(prefix=f"telo_{gid}_") as tmp:
            gff_with_introns = os.path.join(tmp, "with_introns.gff3")
            hits_bed = os.path.join(tmp, "hits.bed")
            patterns = os.path.join(tmp, "patterns.fa")

            run_gt_introns(gff, gff_with_introns)
            introns = load_introns(gff_with_introns)

            write_pattern_file(patterns, opts["motifs"])
            run_seqkit_locate(fa, patterns, hits_bed)
            if os.path.getsize(hits_bed):
                hits = pd.read_csv(hits_bed, sep="\t", header=None,
                                   names=["seqid", "start", "end", "name", "score", "strand"])
            else:
                hits = pd.DataFrame(columns=["seqid", "start", "end", "name", "score", "strand"])

            if introns.empty:
                return {"summary": [gid, org, group, source, 0, 0, len(seqs), "", 0, "OK"],
                        "introns": [], "loci": []}

            introns = introns.reset_index(drop=True)
            introns["row"] = introns.index

            # One bedtools intersect → coverage, flank coverage, strand counts. Replaces
            # three per-intron Python loops; ~30-100× faster on big genomes.
            best_motif, flank_df, strand_df = compute_intron_motif_metrics(
                introns, hits_bed, contig_lens, tmp,
            )

        introns = introns.merge(best_motif, on="row", how="left")
        introns["telomeric_bases"] = introns.telomeric_bases.fillna(0).astype(int)
        introns["telomeric_frac"] = introns.telomeric_frac.fillna(0.0)
        introns["motif"] = introns.motif.fillna("")

        introns = introns.merge(flank_df, on="row", how="left")
        introns["flank_telomeric_frac"] = pd.to_numeric(
            introns.flank_telomeric_frac, errors="coerce"
        ).fillna(0.0)

        introns = introns.merge(strand_df, on="row", how="left").fillna(
            {"fwd_hits": 0, "rev_hits": 0}
        )
        introns["fwd_hits"] = introns.fwd_hits.astype(int)
        introns["rev_hits"] = introns.rev_hits.astype(int)
        introns["orientation"] = [
            classify_orientation(f, r) for f, r in zip(introns.fwd_hits, introns.rev_hits)
        ]

        # Substring extracts: donor, acceptor, first40, last40 (spliced orientation).
        donors, acceptors, first40s, last40s, contig_lens_col, dists = [], [], [], [], [], []
        for r in introns.itertuples(index=False):
            chrom = seqs.get(r.seqid, "")
            genomic = chrom[r.start - 1:r.end]
            spliced = genomic if r.strand == "+" else revcomp(genomic)
            donors.append(spliced[:2])
            acceptors.append(spliced[-2:])
            first40s.append(spliced[:40])
            last40s.append(spliced[-40:])
            clen = len(chrom)
            contig_lens_col.append(clen)
            dists.append(min(r.start - 1, clen - r.end))
        introns["donor"] = donors
        introns["acceptor"] = acceptors
        introns["splice_class"] = [f"{d}-{a}" for d, a in zip(donors, acceptors)]
        # Which strand carries a canonical GT..AG (independent of the host-gene
        # strand). Loci canonical only on the opposite strand surface here as
        # CT-AC etc.; this column makes that explicit instead of silent.
        canon_strand, canon_ok = [], []
        for d, a in zip(donors, acceptors):
            fwd = (d == "GT" and a == "AG")
            rev = (d == "CT" and a == "AC")  # revcomp of GT..AG
            canon_strand.append("+" if fwd else "-" if rev else ".")
            canon_ok.append(fwd or rev)
        introns["canonical_strand"] = canon_strand
        introns["splice_canonical"] = canon_ok
        introns["first40"] = first40s
        introns["last40"] = last40s
        introns["contig_len"] = contig_lens_col
        introns["distance_to_end"] = dists

        # Genome-level terminal motif: prefer literature-curated assignment;
        # fall back to direct scan only when canonical mapping is unset.
        canonical = opts.get("canonical", {}).get(gid, "")
        if canonical:
            term_motif = canonical
            # Bases = how many bp of this motif actually appear within 500 bp of any contig end.
            if hits.empty:
                term_bases = 0
            else:
                h_motif = hits[hits.name == canonical]
                clen_series = h_motif.seqid.map(contig_lens)
                near = h_motif[(h_motif.start < 500) | (h_motif.end > clen_series - 500)]
                term_bases = int((near.end - near.start).sum())
        else:
            term_motif, term_bases = terminal_motif(hits, contig_lens)
        introns["terminal_motif"] = term_motif
        introns["terminal_motif_bases"] = term_bases

        # Per-contig: does the genome's terminal motif appear within 500 bp of
        # BOTH ends? The subtelomeric distance test must restrict to such
        # fully-capped contigs — fragmented contigs lack a real "chromosome end"
        # and inflate distance-to-end (manuscript Section 3 / review).
        capped = set()
        if term_motif and not hits.empty:
            hm = hits[hits.name == term_motif]
            for seqid, g in hm.groupby("seqid"):
                clen = contig_lens.get(seqid, 0)
                if clen and (g.start < 500).any() and (g.end > clen - 500).any():
                    capped.add(seqid)
        introns["contig_both_ends_capped"] = introns.seqid.isin(capped)

        # Manifest fields.
        introns["genome_id"] = gid
        introns["organism"] = org
        introns["group"] = group
        introns["source"] = source

        out = introns[LOCI_HEADER].values.tolist()
        loci = introns[
            (introns.telomeric_frac >= opts["min_repeat_frac"])
            & (introns.flank_telomeric_frac <= opts["max_flank_repeat_frac"])
            & (introns.motif != "")
        ][LOCI_HEADER].values.tolist()
        return {
            "summary": [gid, org, group, source, len(introns), len(loci),
                        len(seqs), term_motif, term_bases, "OK"],
            "introns": out,
            "loci": loci,
        }

    except subprocess.CalledProcessError as e:
        return {"summary": [gid, org, group, source, 0, 0, 0, "", 0, f"ERROR:tool-{e.returncode}"],
                "introns": [], "loci": []}
    except Exception as e:
        return {"summary": [gid, org, group, source, 0, 0, 0, "", 0, f"ERROR:{e}"],
                "introns": [], "loci": []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--canonical-motifs", default="",
                    help="TSV (genome_id, motif) of literature-curated telomere motifs. "
                         "When set for a genome, OVERRIDES the contig-end terminal_motif scan.")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--min-repeat-frac", type=float, default=0.85)
    ap.add_argument("--max-flank-repeat-frac", type=float, default=0.50)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    for tool in ("gt", "seqkit", "bedtools", "sort"):
        if not shutil.which(tool):
            sys.exit(f"required tool not on PATH: {tool}")

    motifs = [m for m in args.motifs.split(",") if m]
    rows = []
    with open(args.manifest) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append([r["genome_id"], r["organism"], r["group"], r["source"]])

    canonical = {}
    if args.canonical_motifs:
        with open(args.canonical_motifs) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["motif"]:
                    canonical[r["genome_id"]] = r["motif"]

    opts = {
        "refseq_dir": args.refseq_dir,
        "tara_dir": args.tara_dir,
        "motifs": motifs,
        "min_repeat_frac": args.min_repeat_frac,
        "max_flank_repeat_frac": args.max_flank_repeat_frac,
        "canonical": canonical,
    }

    os.makedirs(os.path.dirname(args.loci), exist_ok=True)
    with open(args.loci, "w", newline="") as lf, \
         open(args.introns, "w", newline="") as inf, \
         open(args.summary, "w", newline="") as sf:
        loci_w = csv.writer(lf, delimiter="\t")
        intron_w = csv.writer(inf, delimiter="\t")
        summary_w = csv.writer(sf, delimiter="\t")
        loci_w.writerow(LOCI_HEADER)
        intron_w.writerow(LOCI_HEADER)
        summary_w.writerow(SUMMARY_HEADER)

        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = [ex.submit(scan_one, (r, opts)) for r in rows]
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="scanning genomes", unit="genome"):
                res = fut.result()
                summary_w.writerow(res["summary"])
                intron_w.writerows(res["introns"])
                loci_w.writerows(res["loci"])
                sf.flush(); inf.flush(); lf.flush()


if __name__ == "__main__":
    main()
