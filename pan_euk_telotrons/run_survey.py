#!/usr/bin/env python3
"""
Pan-eukaryotic telotron survey — full NCBI RefSeq scan.

Downloads ALL annotated eukaryotic genomes from NCBI, extracts introns,
screens for telomeric content (≥85% purity), and runs TERC boundary analysis.

Design:
  - 15-core parallelism (one genome per worker)
  - Streaming FASTA: only contigs with genes are loaded, one at a time
  - Disk-conscious: raw genome files cleaned up after each genome
  - Resumable: checkpoint.json tracks per-accession state
  - Genome size cap: skip genomes >5 GB to avoid memory/disk blowout
"""

import json, os, re, sys, shutil, subprocess, threading, time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASEDIR    = Path(__file__).resolve().parent
DATASETS   = BASEDIR.parent / "pan_euk_telotrons_datasets"
ACCESSION_TSV = BASEDIR / "all_euk_accessions.tsv"
WORKERS    = 15
CHECKPOINT = BASEDIR / "checkpoint.json"
LOCK       = threading.Lock()

MAX_GENOME_BP   = 5_000_000_000   # skip genomes >5 Gbp
DOWNLOAD_TIMEOUT = 1800           # 30 min per genome download
UNZIP_TIMEOUT    = 300            # 5 min unzip

DIR_GENOMES   = BASEDIR / "genomes"
DIR_TELOTRONS = BASEDIR / "telotrons"
DIR_RESULTS   = BASEDIR / "results"

# ---------------------------------------------------------------------------
# Telomeric hexamer sets
# ---------------------------------------------------------------------------
def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

TELO_HEX_6 = _rotations("TTAGGG") | _rotations("CCCTAA")

VARIANT_SETS = {
    "TTTAGGG":  _rotations("TTTAGGG")  | _rotations("CCCTAAA"),
    "TTGGGG":   _rotations("TTGGGG")   | _rotations("CCCCAA"),
    "TTTGGG":   _rotations("TTTGGG")   | _rotations("CCCAAA"),
    "TTAGG":    _rotations("TTAGG")    | _rotations("CCTAA"),
    "TTTTAGGG": _rotations("TTTTAGGG") | _rotations("CCCTAAAA"),
}

ALL_SCAN_SETS = {"TTAGGG": TELO_HEX_6, **VARIANT_SETS}

# ---------------------------------------------------------------------------
# Load genome list from TSV
# ---------------------------------------------------------------------------
def load_genome_list():
    """Load accession list. Returns [(acc, organism_name, genome_bp), ...]"""
    genomes = []
    with open(ACCESSION_TSV) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            acc = parts[0]
            name = parts[1].replace(" ", "_").replace("/", "_")
            try:
                bp = int(parts[4])
            except (ValueError, IndexError):
                bp = 0
            if acc:
                genomes.append((acc, name, bp))
    return genomes

# ---------------------------------------------------------------------------
# Checkpoint (thread-safe)
# ---------------------------------------------------------------------------
def load_checkpoint():
    if CHECKPOINT.exists():
        try:
            return json.loads(CHECKPOINT.read_text())
        except json.JSONDecodeError:
            return {}
    return {}

def mark_done(acc, status_dict):
    with LOCK:
        ckpt = load_checkpoint()
        ckpt[acc] = status_dict
        tmp = CHECKPOINT.with_suffix(".tmp")
        tmp.write_text(json.dumps(ckpt))
        tmp.rename(CHECKPOINT)

# ---------------------------------------------------------------------------
# Streaming FASTA: only load contigs that have genes
# ---------------------------------------------------------------------------
def load_fasta_subset(fna_path, needed_contigs):
    """Load only the contigs in needed_contigs set. Memory-efficient."""
    seqs = {}
    curr = None
    parts = []
    keep = False
    with open(fna_path) as f:
        for line in f:
            if line.startswith(">"):
                if keep and curr:
                    seqs[curr] = "".join(parts)
                curr = line[1:].strip().split()[0]
                parts = []
                keep = curr in needed_contigs
            elif keep:
                parts.append(line.strip().upper())
    if keep and curr:
        seqs[curr] = "".join(parts)
    return seqs

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def download_genome(acc):
    out_dir = DIR_GENOMES / acc
    zip_path = DIR_GENOMES / f"{acc}.zip"

    if out_dir.exists():
        gff = _find_file(out_dir, "*.gff")
        fna = _find_file(out_dir, "*.fna")
        if gff and fna:
            return (gff, fna)

    try:
        subprocess.run(
            [str(DATASETS), "download", "genome", "accession", acc,
             "--include", "genome,gff3", "--filename", str(zip_path)],
            check=True, capture_output=True, timeout=DOWNLOAD_TIMEOUT,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        zip_path.unlink(missing_ok=True)
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["unzip", "-qo", str(zip_path), "-d", str(out_dir)],
            check=True, capture_output=True, timeout=UNZIP_TIMEOUT,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        shutil.rmtree(out_dir, ignore_errors=True)
        return None
    finally:
        zip_path.unlink(missing_ok=True)

    gff = _find_file(out_dir, "*.gff")
    fna = _find_file(out_dir, "*.fna")
    return (gff, fna) if gff and fna else None


def _find_file(root, pattern):
    matches = list(Path(root).rglob(pattern))
    return str(matches[0]) if matches else None

# ---------------------------------------------------------------------------
# Parse GFF → per-contig exon structure (no genome needed)
# ---------------------------------------------------------------------------
def parse_gff_exons(gff_path):
    """Returns {contig: {mrna: [(start,end)]}} and {mrna: strand}."""
    contig_mrna_exons = defaultdict(lambda: defaultdict(list))
    mrna_strand = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in ("exon", "CDS"):
                continue
            attrs = dict(re.findall(r"(\w+)=([^;]+)", parts[8]))
            parent = attrs.get("Parent", "")
            contig = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            contig_mrna_exons[contig][parent].append((start, end))
            mrna_strand[parent] = strand

    return contig_mrna_exons, mrna_strand

# ---------------------------------------------------------------------------
# Extract introns from parsed structure + loaded contigs
# ---------------------------------------------------------------------------
def extract_introns_from_contigs(contig_mrna_exons, mrna_strand, genome_seqs):
    introns = []
    for contig, mrna_exons in contig_mrna_exons.items():
        if contig not in genome_seqs:
            continue
        seq = genome_seqs[contig]
        for mrna, exons in mrna_exons.items():
            exons.sort()
            strand = mrna_strand.get(mrna, "+")
            for i in range(len(exons) - 1):
                istart = exons[i][1]
                iend = exons[i + 1][0]
                if iend <= istart:
                    continue
                length = iend - istart
                if length < 20 or length > 100000:
                    continue
                intron_seq = seq[istart:iend]
                if not intron_seq:
                    continue
                introns.append({
                    "contig": contig, "start": istart, "end": iend,
                    "length": length, "strand": strand, "mrna": mrna,
                    "donor": intron_seq[:2], "acceptor": intron_seq[-2:],
                    "seq": intron_seq,
                })
    return introns

# ---------------------------------------------------------------------------
# Screen + boundary analysis (unchanged logic)
# ---------------------------------------------------------------------------
def telo_coverage(seq, hexset):
    n = len(seq)
    if n == 0:
        return 0.0
    covered = bytearray(n)
    for h in hexset:
        start = 0
        while True:
            idx = seq.find(h, start)
            if idx == -1:
                break
            end = min(idx + len(h), n)
            for j in range(idx, end):
                covered[j] = 1
            start = idx + 1
    return sum(covered) / n


def strand_ratio(seq):
    fwd_pats = ["TTAGGG", "TAGGGT", "AGGGTT", "GGGTTA", "GGTTAG", "GTTAGG"]
    rev_pats = ["CCCTAA", "CCTAAC", "CTAACC", "TAACCC", "AACCCT", "ACCCTA"]
    fwd = sum(len(re.findall(p, seq)) for p in fwd_pats)
    rev = sum(len(re.findall(p, seq)) for p in rev_pats)
    total = fwd + rev
    if total == 0:
        return 0.5, "unknown"
    frac = fwd / total
    if frac > 0.7:
        return frac, "coding-strand"
    elif frac < 0.3:
        return frac, "template-strand"
    return frac, "converging"


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def screen_introns(introns):
    """Screen all introns against all repeat types. ≥85% purity."""
    telotrons = []
    for intron in introns:
        seq = intron["seq"]
        best_cov = 0.0
        best_repeat = None

        for rtype, hexset in ALL_SCAN_SETS.items():
            cov = telo_coverage(seq, hexset)
            if cov > best_cov:
                best_cov = cov
                best_repeat = rtype

        if best_cov < 0.85:
            continue

        frac, orient = strand_ratio(seq)
        norm_seq = seq
        if seq[:2] == "CT" and seq[-2:] == "AC":
            norm_seq = rc(seq)

        donor_8 = norm_seq[:8] if len(norm_seq) >= 8 else norm_seq
        acceptor_8 = norm_seq[-8:] if len(norm_seq) >= 8 else norm_seq

        telotrons.append({
            **intron,
            "telo_coverage": best_cov,
            "repeat_type": best_repeat,
            "fwd_frac": frac,
            "orientation": orient,
            "donor_8": donor_8,
            "acceptor_8": acceptor_8,
            "norm_seq": norm_seq,
        })
    return telotrons


def analyze_boundaries(telotrons, genome_id):
    stats = {
        "genome_id": genome_id, "n_total": 0, "n_gt_ag": 0,
        "donor_4": Counter(), "acceptor_4": Counter(),
        "by_orient": defaultdict(lambda: {"n": 0, "acceptor_4": Counter()}),
    }
    for t in telotrons:
        stats["n_total"] += 1
        seq = t["norm_seq"]
        if seq[:2] != "GT" or seq[-2:] != "AG":
            continue
        stats["n_gt_ag"] += 1
        if len(seq) >= 8:
            d4 = seq[2:6]
            a4 = seq[-6:-2]
            stats["donor_4"][d4] += 1
            stats["acceptor_4"][a4] += 1
            orient = t["orientation"]
            stats["by_orient"][orient]["n"] += 1
            stats["by_orient"][orient]["acceptor_4"][a4] += 1
    return stats

# ---------------------------------------------------------------------------
# Per-genome worker
# ---------------------------------------------------------------------------
def process_one_genome(acc, name, genome_bp):
    result = {"acc": acc, "name": name, "genome_bp": genome_bp,
              "status": "failed", "error": None}

    if genome_bp > MAX_GENOME_BP:
        result["error"] = f"skipped_too_large_{genome_bp/1e9:.1f}GB"
        return result

    try:
        paths = download_genome(acc)
        if paths is None:
            result["error"] = "download_failed"
            return result
        gff_path, fna_path = paths

        # Parse GFF first (small memory)
        contig_mrna_exons, mrna_strand = parse_gff_exons(gff_path)
        needed_contigs = set(contig_mrna_exons.keys())

        # Load only gene-bearing contigs
        genome_seqs = load_fasta_subset(fna_path, needed_contigs)

        # Extract introns
        introns = extract_introns_from_contigs(contig_mrna_exons, mrna_strand, genome_seqs)
        del genome_seqs  # free memory immediately
        result["n_introns"] = len(introns)

        # Screen
        telotrons = screen_introns(introns)
        del introns
        result["n_telotrons"] = len(telotrons)

        # Only write output files if telotrons found
        if telotrons:
            out_path = DIR_TELOTRONS / f"{acc}_{name}_telotrons.tsv"
            _write_telotrons_tsv(telotrons, out_path, acc, name)

            boundary_stats = analyze_boundaries(telotrons, f"{acc}_{name}")
            result["boundary"] = _boundary_summary(boundary_stats)

            detail_path = DIR_RESULTS / f"{acc}_{name}_boundary.json"
            detail_path.write_text(json.dumps({
                "genome_id": f"{acc}_{name}",
                "n_telotrons": len(telotrons),
                "n_gt_ag": boundary_stats["n_gt_ag"],
                "top_donor_4": boundary_stats["donor_4"].most_common(5),
                "top_acceptor_4": boundary_stats["acceptor_4"].most_common(5),
                "orientation_breakdown": {
                    o: {"n": d["n"], "top_acceptor": d["acceptor_4"].most_common(3)}
                    for o, d in boundary_stats["by_orient"].items()
                },
            }, indent=2))

        result["status"] = "ok"

    except Exception as e:
        result["error"] = str(e)

    finally:
        genome_dir = DIR_GENOMES / acc
        if genome_dir.exists():
            shutil.rmtree(genome_dir, ignore_errors=True)
        (DIR_GENOMES / f"{acc}.zip").unlink(missing_ok=True)

    return result


def _write_telotrons_tsv(telotrons, path, acc, name):
    cols = ["contig", "start", "end", "length", "strand", "mrna",
            "donor", "acceptor", "telo_coverage", "repeat_type",
            "fwd_frac", "orientation", "donor_8", "acceptor_8"]
    with open(path, "w") as f:
        f.write("genome_acc\tgenome_name\t" + "\t".join(cols) + "\tintron_seq\n")
        for t in telotrons:
            vals = [str(t.get(c, "")) for c in cols]
            f.write(f"{acc}\t{name}\t" + "\t".join(vals) + f"\t{t['seq']}\n")


def _boundary_summary(stats):
    total_a = sum(stats["acceptor_4"].values())
    aacc_n = stats["acceptor_4"].get("AACC", 0)
    aacc_pct = 100 * aacc_n / total_a if total_a > 0 else 0
    orient_aacc = {}
    for o in ["template-strand", "converging", "coding-strand"]:
        d = stats["by_orient"][o]
        ot = sum(d["acceptor_4"].values())
        orient_aacc[o] = 100 * d["acceptor_4"].get("AACC", 0) / ot if ot > 0 else 0
    top_a = stats["acceptor_4"].most_common(1)
    top_d = stats["donor_4"].most_common(1)
    return {
        "n_gt_ag": stats["n_gt_ag"],
        "top_donor": top_d[0] if top_d else None,
        "top_acceptor": top_a[0] if top_a else None,
        "AACC_AG_pct": round(aacc_pct, 1),
        "orient_AACC": {k: round(v, 1) for k, v in orient_aacc.items()},
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for d in [DIR_GENOMES, DIR_TELOTRONS, DIR_RESULTS]:
        d.mkdir(parents=True, exist_ok=True)

    all_genomes = load_genome_list()
    ckpt = load_checkpoint()
    todo = [(a, n, bp) for a, n, bp in all_genomes if a not in ckpt]

    # Sort: small genomes first for fast initial throughput
    todo.sort(key=lambda x: x[2])

    print(f"=== Pan-eukaryotic telotron survey ===")
    print(f"  {len(all_genomes)} total genomes, {len(ckpt)} already done, {len(todo)} remaining")
    print(f"  Skipping genomes > {MAX_GENOME_BP/1e9:.0f} GB")
    print(f"  Workers: {WORKERS}")
    print()

    if not todo:
        print("All genomes processed.")
        _write_final_report(all_genomes, ckpt)
        return

    done_count = len(ckpt)
    total_count = len(all_genomes)
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(process_one_genome, a, n, bp): (a, n)
                   for a, n, bp in todo}

        for fut in as_completed(futures):
            acc, name = futures[fut]
            try:
                result = fut.result()
            except Exception as e:
                result = {"acc": acc, "name": name, "status": "crashed", "error": str(e)}

            mark_done(acc, result)
            done_count += 1

            status = result.get("status", "?")
            n_int = result.get("n_introns", 0)
            n_tel = result.get("n_telotrons", 0)
            err = result.get("error", "")

            elapsed = time.time() - t0
            rate = (done_count - len(ckpt)) / elapsed * 3600 if elapsed > 0 else 0

            if n_tel > 0:
                tag = f"HIT({n_tel})"
            elif status == "ok":
                tag = "ok"
            else:
                tag = err[:18] if err else "fail"

            print(f"  [{done_count:>4}/{total_count}] [{tag:.<20s}] "
                  f"{name:.<45s} introns={n_int:>9,}  "
                  f"telotrons={n_tel:>3}  ({rate:.0f}/hr)")

    ckpt = load_checkpoint()
    _write_final_report(all_genomes, ckpt)


def _write_final_report(all_genomes, ckpt):
    summary_path = DIR_RESULTS / "pan_euk_summary.tsv"
    hits_path = DIR_RESULTS / "pan_euk_hits.tsv"
    report_path = DIR_RESULTS / "pan_euk_report.md"

    # Summary TSV (all genomes)
    with open(summary_path, "w") as f:
        f.write("accession\tspecies\tgenome_bp\tstatus\ttotal_introns\t"
                "telotrons_85pct\tn_gt_ag\ttop_acceptor\ttop_donor\n")
        for acc, name, bp in all_genomes:
            r = ckpt.get(acc, {})
            bnd = r.get("boundary", {})
            top_a = bnd.get("top_acceptor", "") if bnd else ""
            top_d = bnd.get("top_donor", "") if bnd else ""
            f.write(f"{acc}\t{name}\t{bp}\t{r.get('status','missing')}\t"
                    f"{r.get('n_introns',0)}\t{r.get('n_telotrons',0)}\t"
                    f"{bnd.get('n_gt_ag',0) if bnd else 0}\t{top_a}\t{top_d}\n")

    # Hits-only TSV
    hits = [(acc, name, bp, ckpt[acc]) for acc, name, bp in all_genomes
            if acc in ckpt and ckpt[acc].get("n_telotrons", 0) > 0]

    with open(hits_path, "w") as f:
        f.write("accession\tspecies\tgenome_bp\tn_introns\tn_telotrons\t"
                "n_gt_ag\ttop_acceptor\ttop_donor\tAACC_AG_pct\n")
        for acc, name, bp, r in sorted(hits, key=lambda x: -x[3].get("n_telotrons", 0)):
            bnd = r.get("boundary", {})
            top_a = bnd.get("top_acceptor", "")
            top_d = bnd.get("top_donor", "")
            f.write(f"{acc}\t{name}\t{bp}\t{r.get('n_introns',0)}\t"
                    f"{r.get('n_telotrons',0)}\t"
                    f"{bnd.get('n_gt_ag',0) if bnd else 0}\t"
                    f"{top_a}\t{top_d}\t{bnd.get('AACC_AG_pct',0) if bnd else 0}\n")

    # Markdown report
    ok = sum(1 for v in ckpt.values() if v.get("status") == "ok")
    fail = sum(1 for v in ckpt.values() if v.get("status") != "ok")
    total_introns = sum(v.get("n_introns", 0) for v in ckpt.values())
    total_telotrons = sum(v.get("n_telotrons", 0) for v in ckpt.values())

    lines = [
        "# Pan-Eukaryotic Telotron Survey — All NCBI RefSeq Eukaryotes\n\n",
        f"**Date**: {__import__('datetime').datetime.now().isoformat()[:10]}\n",
        f"**Genomes in catalog**: {len(all_genomes)}\n",
        f"**Successfully processed**: {ok}\n",
        f"**Failed/skipped**: {fail}\n",
        f"**Total introns scanned**: {total_introns:,}\n",
        f"**Telotrons found (≥85% purity)**: {total_telotrons}\n",
        f"**Species with telotrons**: {len(hits)}\n\n",
    ]

    if hits:
        lines.append("## Species with Telotrons (≥85% telomeric purity)\n\n")
        lines.append("| Species | Accession | Genome (Mbp) | Introns | Telotrons | GT-AG | Top 3' motif |\n")
        lines.append("|---|---|---:|---:|---:|---:|---|\n")
        for acc, name, bp, r in sorted(hits, key=lambda x: -x[3].get("n_telotrons", 0)):
            bnd = r.get("boundary", {})
            top_a = bnd.get("top_acceptor")
            top_a_str = f"{top_a[0]}({top_a[1]})" if top_a else "—"
            lines.append(
                f"| {name.replace('_',' ')} | {acc} | {bp/1e6:.0f} "
                f"| {r.get('n_introns',0):,} | {r.get('n_telotrons',0)} "
                f"| {bnd.get('n_gt_ag',0) if bnd else 0} | {top_a_str} |\n"
            )
    else:
        lines.append("## No telotrons found at ≥85% purity.\n")

    report_path.write_text("".join(lines))
    print(f"\n=== Report:     {report_path} ===")
    print(f"=== Summary:    {summary_path} ===")
    print(f"=== Hits only:  {hits_path} ===")
    print(f"=== {ok} processed, {total_telotrons} telotrons in {len(hits)} species ===")


if __name__ == "__main__":
    main()
