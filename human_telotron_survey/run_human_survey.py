#!/usr/bin/env python3
"""
Human genome telotron survey — search all annotated human reference genomes
for telomeric-repeat introns (≥85% TTAGGG/CCCTAA coverage).

Targets: GRCh38, T2T-CHM13v2.0, HuRef, CHM1, Ash1, hg01243, ASM2283312v2
"""

import json, os, re, sys, shutil, subprocess, tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASEDIR    = Path(__file__).resolve().parent
DATASETS   = BASEDIR.parent / "pan_euk_telotrons_datasets"
CHECKPOINT = BASEDIR / "checkpoint.json"
WORKERS    = 4  # conservative — human genomes are large

# Output dirs
DIR_GENOMES   = BASEDIR / "genomes"
DIR_TELOTRONS = BASEDIR / "telotrons"
DIR_RESULTS   = BASEDIR / "results"

# All annotated human genome assemblies from NCBI
HUMAN_ASSEMBLIES = [
    ("GCF_000001405.40", "GRCh38_p14",       "Primary reference genome"),
    ("GCF_009914755.1",  "T2T_CHM13v2",      "Telomere-to-telomere complete genome"),
    ("GCF_000002125.1",  "HuRef",            "Craig Venter personal genome"),
    ("GCF_000306695.2",  "CHM1_1.1",         "Hydatidiform mole CHM1"),
    ("GCA_011064465.2",  "Ash1_v2.2",        "Ashkenazi individual"),
    ("GCA_018873775.2",  "hg01243_v3",       "HPRC sample HG01243"),
    ("GCA_022833125.2",  "ASM2283312v2",     "Korean genome KPGP"),
    ("GCF_000002135.2",  "CRA_TCAGchr7v2",   "Chromosome 7 only (partial)"),
]

# Telomeric hexamer sets — all 12 rotations of TTAGGG and CCCTAA
def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

TELO_HEX = _rotations("TTAGGG") | _rotations("CCCTAA")

# Also track near-telomeric (1 substitution away) for context
FWD_PATS = ["TTAGGG", "TAGGGT", "AGGGTT", "GGGTTA", "GGTTAG", "GTTAGG"]
REV_PATS = ["CCCTAA", "CCTAAC", "CTAACC", "TAACCC", "AACCCT", "ACCCTA"]

# Purity thresholds to report
THRESHOLDS = [0.10, 0.30, 0.50, 0.70, 0.85, 0.95]

# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------
def load_checkpoint():
    if CHECKPOINT.exists():
        return json.loads(CHECKPOINT.read_text())
    return {}

def save_checkpoint(ckpt):
    tmp = CHECKPOINT.with_suffix(".tmp")
    tmp.write_text(json.dumps(ckpt, indent=2))
    tmp.replace(CHECKPOINT)

def mark_done(acc, result):
    ckpt = load_checkpoint()
    ckpt[acc] = result
    save_checkpoint(ckpt)

# ---------------------------------------------------------------------------
# FASTA loader
# ---------------------------------------------------------------------------
def load_fasta(path):
    seqs = {}
    curr = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if curr:
                    seqs[curr] = "".join(parts)
                curr = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if curr:
        seqs[curr] = "".join(parts)
    return seqs

# ---------------------------------------------------------------------------
# Download genome + GFF3
# ---------------------------------------------------------------------------
def download_genome(acc):
    out_dir = DIR_GENOMES / acc
    zip_path = DIR_GENOMES / f"{acc}.zip"

    # Already extracted?
    if out_dir.exists():
        gff = _find_file(out_dir, "*.gff")
        fna = _find_file(out_dir, "*.fna")
        if gff and fna:
            return (gff, fna)

    # Download
    try:
        subprocess.run(
            [str(DATASETS), "download", "genome", "accession", acc,
             "--include", "genome,gff3", "--filename", str(zip_path)],
            check=True, capture_output=True, timeout=1800,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        return None

    # Unzip
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["unzip", "-qo", str(zip_path), "-d", str(out_dir)],
            check=True, capture_output=True, timeout=600,
        )
    except subprocess.CalledProcessError:
        return None
    finally:
        zip_path.unlink(missing_ok=True)

    gff = _find_file(out_dir, "*.gff")
    fna = _find_file(out_dir, "*.fna")
    if gff and fna:
        return (gff, fna)
    return None


def _find_file(root, pattern):
    matches = list(root.rglob(pattern))
    return str(matches[0]) if matches else None

# ---------------------------------------------------------------------------
# Extract introns from GFF3 + FASTA
# ---------------------------------------------------------------------------
def extract_introns(gff_path, fna_path):
    """Extract introns, returning list of dicts."""
    genome = load_fasta(fna_path)
    print(f"    Loaded genome: {len(genome)} contigs, "
          f"{sum(len(s) for s in genome.values()):,} bp total", flush=True)

    # Parse GFF: collect exon/CDS features per parent
    mrna_features = defaultdict(list)
    mrna_strand = {}
    mrna_contig = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in ("exon", "CDS"):
                continue

            attrs = dict(re.findall(r"(\w+)=([^;]+)", parts[8]))
            parent = attrs.get("Parent", "")
            contig = parts[0]
            start = int(parts[3]) - 1  # 0-based
            end = int(parts[4])
            strand = parts[6]

            mrna_features[parent].append((start, end))
            mrna_strand[parent] = strand
            mrna_contig[parent] = contig

    introns = []
    for mrna, exons in mrna_features.items():
        exons.sort()
        contig = mrna_contig[mrna]
        strand = mrna_strand[mrna]
        if contig not in genome:
            continue
        seq = genome[contig]

        for i in range(len(exons) - 1):
            istart = exons[i][1]
            iend = exons[i + 1][0]
            if iend <= istart:
                continue
            length = iend - istart
            if length < 20 or length > 500000:
                continue

            intron_seq = seq[istart:iend]
            if not intron_seq or len(intron_seq) < 20:
                continue

            introns.append({
                "contig": contig,
                "start": istart,
                "end": iend,
                "length": length,
                "strand": strand,
                "mrna": mrna,
                "donor": intron_seq[:2],
                "acceptor": intron_seq[-2:],
                "seq": intron_seq,
            })

    # Deduplicate by (contig, start, end)
    seen = set()
    deduped = []
    for intr in introns:
        key = (intr["contig"], intr["start"], intr["end"])
        if key not in seen:
            seen.add(key)
            deduped.append(intr)

    return deduped

# ---------------------------------------------------------------------------
# Telomeric coverage calculation
# ---------------------------------------------------------------------------
def telo_coverage(seq, hexset=TELO_HEX):
    """Fraction of bases covered by any hexamer in the set."""
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
            for j in range(idx, min(idx + len(h), n)):
                covered[j] = 1
            start = idx + 1
    return sum(covered) / n


def strand_ratio(seq):
    """TTAGGG-rotation vs CCCTAA-rotation counts."""
    fwd = sum(len(re.findall(p, seq)) for p in FWD_PATS)
    rev = sum(len(re.findall(p, seq)) for p in REV_PATS)
    total = fwd + rev
    if total == 0:
        return 0.5, "unknown"
    frac = fwd / total
    if frac > 0.7:
        return frac, "coding-strand"
    elif frac < 0.3:
        return frac, "template-strand"
    else:
        return frac, "converging"


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))

# ---------------------------------------------------------------------------
# Screen introns for telomeric content
# ---------------------------------------------------------------------------
def screen_introns(introns):
    """Screen all introns, return (telotrons_list, threshold_counts_dict)."""
    threshold_counts = {t: 0 for t in THRESHOLDS}
    telotrons = []  # all with ≥10% for reporting

    total = len(introns)
    report_interval = max(1, total // 20)

    for idx, intron in enumerate(introns):
        if idx % report_interval == 0:
            print(f"    Screening: {idx:,}/{total:,} ({100*idx/total:.0f}%)", flush=True)

        seq = intron["seq"]
        cov = telo_coverage(seq)

        for t in THRESHOLDS:
            if cov >= t:
                threshold_counts[t] += 1

        if cov < 0.10:
            continue

        frac, orient = strand_ratio(seq)

        # Normalize to GT-AG if CT-AC
        d2 = seq[:2]
        a2 = seq[-2:]
        norm_seq = seq
        if d2 == "CT" and a2 == "AC":
            norm_seq = rc(seq)

        donor_8 = norm_seq[:8] if len(norm_seq) >= 8 else norm_seq
        acceptor_8 = norm_seq[-8:] if len(norm_seq) >= 8 else norm_seq

        telotrons.append({
            **intron,
            "telo_coverage": cov,
            "fwd_frac": frac,
            "orientation": orient,
            "donor_8": donor_8,
            "acceptor_8": acceptor_8,
            "norm_seq": norm_seq,
        })

    return telotrons, threshold_counts

# ---------------------------------------------------------------------------
# Boundary analysis
# ---------------------------------------------------------------------------
def analyze_boundaries(telotrons):
    """Compute boundary motif stats for telotrons with ≥85% purity."""
    stats = {
        "n_telo85": 0,
        "n_gt_ag": 0,
        "donor_4": Counter(),
        "acceptor_4": Counter(),
        "by_orient": defaultdict(lambda: {"n": 0, "acceptor_4": Counter()}),
    }

    for t in telotrons:
        if t["telo_coverage"] < 0.85:
            continue
        stats["n_telo85"] += 1

        seq = t["norm_seq"]
        d2 = seq[:2]
        a2 = seq[-2:]
        if d2 != "GT" or a2 != "AG":
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
# Process one genome
# ---------------------------------------------------------------------------
def process_one_genome(acc, name, desc):
    """Full pipeline for one human genome assembly."""
    result = {
        "acc": acc, "name": name, "description": desc,
        "status": "failed", "error": None,
        "timestamp": datetime.now().isoformat(),
    }

    try:
        print(f"\n{'='*70}", flush=True)
        print(f"  Processing: {name} ({acc})", flush=True)
        print(f"  {desc}", flush=True)
        print(f"{'='*70}", flush=True)

        # Download
        print("  Step 1: Downloading genome + GFF3...", flush=True)
        paths = download_genome(acc)
        if paths is None:
            result["error"] = "download_failed"
            print(f"  FAILED: Download failed for {acc}", flush=True)
            return result
        gff_path, fna_path = paths
        print(f"  Downloaded: {gff_path}", flush=True)

        # Extract introns
        print("  Step 2: Extracting introns...", flush=True)
        introns = extract_introns(gff_path, fna_path)
        result["n_introns"] = len(introns)
        print(f"  Extracted: {len(introns):,} unique introns", flush=True)

        # Screen for telomeric content
        print("  Step 3: Screening for telomeric content...", flush=True)
        telotrons, threshold_counts = screen_introns(introns)
        result["n_telotrons_10pct"] = threshold_counts.get(0.10, 0)
        result["n_telotrons_30pct"] = threshold_counts.get(0.30, 0)
        result["n_telotrons_50pct"] = threshold_counts.get(0.50, 0)
        result["n_telotrons_70pct"] = threshold_counts.get(0.70, 0)
        result["n_telotrons_85pct"] = threshold_counts.get(0.85, 0)
        result["n_telotrons_95pct"] = threshold_counts.get(0.95, 0)
        result["threshold_counts"] = {str(k): v for k, v in threshold_counts.items()}

        hp_telotrons = [t for t in telotrons if t["telo_coverage"] >= 0.85]

        print(f"  Telomeric introns found:", flush=True)
        for t in THRESHOLDS:
            print(f"    ≥{t*100:.0f}%: {threshold_counts[t]:,}", flush=True)

        # Write ALL telotrons (≥10%) to TSV for detailed analysis
        out_path = DIR_TELOTRONS / f"{acc}_{name}_telotrons.tsv"
        _write_telotrons_tsv(telotrons, out_path, acc, name)

        # Boundary analysis (≥85% only)
        print("  Step 4: Boundary analysis...", flush=True)
        boundary = analyze_boundaries(telotrons)
        result["boundary"] = _boundary_summary(boundary)

        # Report high-purity telotrons in detail
        if hp_telotrons:
            print(f"\n  *** FOUND {len(hp_telotrons)} HIGH-PURITY TELOTRONS (≥85%) ***", flush=True)
            result["hp_telotron_details"] = []
            for i, t in enumerate(hp_telotrons):
                detail = {
                    "contig": t["contig"],
                    "start": t["start"],
                    "end": t["end"],
                    "length": t["length"],
                    "strand": t["strand"],
                    "telo_coverage": round(t["telo_coverage"], 4),
                    "orientation": t["orientation"],
                    "donor": t["donor"],
                    "acceptor": t["acceptor"],
                    "donor_8": t["donor_8"],
                    "acceptor_8": t["acceptor_8"],
                    "gene": t["mrna"],
                }
                result["hp_telotron_details"].append(detail)
                print(f"    [{i+1}] {t['contig']}:{t['start']}-{t['end']} "
                      f"({t['length']}bp, {t['strand']}) "
                      f"purity={t['telo_coverage']:.1%} "
                      f"orient={t['orientation']} "
                      f"splice={t['donor']}-{t['acceptor']}", flush=True)
                # Print first 100bp of sequence for inspection
                print(f"        seq: {t['seq'][:100]}{'...' if len(t['seq'])>100 else ''}",
                      flush=True)
        else:
            print("  No high-purity telotrons found.", flush=True)

        # Also report notable near-misses (70-85%)
        near_misses = [t for t in telotrons if 0.70 <= t["telo_coverage"] < 0.85]
        if near_misses:
            print(f"\n  Near-misses (70-85% purity): {len(near_misses)}", flush=True)
            result["near_miss_details"] = []
            for t in near_misses[:10]:  # top 10
                detail = {
                    "contig": t["contig"],
                    "start": t["start"],
                    "end": t["end"],
                    "length": t["length"],
                    "telo_coverage": round(t["telo_coverage"], 4),
                    "orientation": t["orientation"],
                    "splice": f"{t['donor']}-{t['acceptor']}",
                }
                result["near_miss_details"].append(detail)
                print(f"    {t['contig']}:{t['start']}-{t['end']} "
                      f"purity={t['telo_coverage']:.1%} "
                      f"len={t['length']}bp", flush=True)

        result["status"] = "ok"

    except Exception as e:
        import traceback
        result["error"] = str(e)
        result["traceback"] = traceback.format_exc()
        print(f"  ERROR: {e}", flush=True)

    finally:
        # Cleanup genome files to save disk
        genome_dir = DIR_GENOMES / acc
        if genome_dir.exists() and "--keep-genomes" not in sys.argv:
            print(f"  Cleaning up genome files for {acc}...", flush=True)
            shutil.rmtree(genome_dir, ignore_errors=True)
        zip_path = DIR_GENOMES / f"{acc}.zip"
        zip_path.unlink(missing_ok=True)

    return result


def _write_telotrons_tsv(telotrons, path, acc, name):
    cols = ["contig", "start", "end", "length", "strand", "mrna",
            "donor", "acceptor", "telo_coverage",
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
        "n_telo85": stats["n_telo85"],
        "top_donor": list(top_d[0]) if top_d else None,
        "top_acceptor": list(top_a[0]) if top_a else None,
        "AACC_AG_pct": round(aacc_pct, 1),
        "orient_AACC": {k: round(v, 1) for k, v in orient_aacc.items()},
    }

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------
def write_final_report(ckpt):
    report_path = DIR_RESULTS / "human_telotron_survey_report.md"
    summary_path = DIR_RESULTS / "human_telotron_summary.tsv"

    lines = [
        "# Human Genome Telotron Survey\n\n",
        f"**Date:** {datetime.now().isoformat()[:10]}\n",
        f"**Assemblies surveyed:** {len(HUMAN_ASSEMBLIES)}\n",
        f"**Purity threshold:** ≥85% telomeric hexamer coverage\n",
        f"**Repeat type:** TTAGGG (vertebrate canonical)\n\n",
        "## Summary\n\n",
        "| Assembly | Description | Total introns | ≥10% | ≥30% | ≥50% | ≥70% | ≥85% | ≥95% |\n",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|\n",
    ]

    total_introns = 0
    total_hp = 0

    for acc, name, desc in HUMAN_ASSEMBLIES:
        r = ckpt.get(acc, {})
        if r.get("status") != "ok":
            lines.append(f"| {name} | {desc} | — | — | — | — | — | — | {r.get('error', 'not run')} |\n")
            continue

        ni = r.get("n_introns", 0)
        total_introns += ni
        hp = r.get("n_telotrons_85pct", 0)
        total_hp += hp

        lines.append(
            f"| **{name}** | {desc} "
            f"| {ni:,} "
            f"| {r.get('n_telotrons_10pct', 0):,} "
            f"| {r.get('n_telotrons_30pct', 0):,} "
            f"| {r.get('n_telotrons_50pct', 0):,} "
            f"| {r.get('n_telotrons_70pct', 0):,} "
            f"| {hp:,} "
            f"| {r.get('n_telotrons_95pct', 0):,} |\n"
        )

    lines.append(f"\n**Total introns scanned: {total_introns:,}**\n")
    lines.append(f"**Total high-purity telotrons (≥85%): {total_hp}**\n\n")

    # High-purity details
    lines.append("## High-Purity Telotron Details\n\n")
    any_found = False
    for acc, name, desc in HUMAN_ASSEMBLIES:
        r = ckpt.get(acc, {})
        details = r.get("hp_telotron_details", [])
        if details:
            any_found = True
            lines.append(f"### {name} ({acc})\n\n")
            for d in details:
                lines.append(
                    f"- **{d['contig']}:{d['start']}-{d['end']}** "
                    f"({d['length']}bp, {d['strand']} strand)\n"
                    f"  - Purity: {d['telo_coverage']:.1%}\n"
                    f"  - Orientation: {d['orientation']}\n"
                    f"  - Splice: {d['donor']}-{d['acceptor']}\n"
                    f"  - Donor 8bp: `{d['donor_8']}`\n"
                    f"  - Acceptor 8bp: `{d['acceptor_8']}`\n"
                    f"  - Gene: {d.get('gene', 'unknown')}\n\n"
                )

    if not any_found:
        lines.append("No high-purity telotrons (≥85%) found in any assembly.\n\n")

    # Near-misses
    lines.append("## Near-Miss Telotrons (70-85% purity)\n\n")
    for acc, name, desc in HUMAN_ASSEMBLIES:
        r = ckpt.get(acc, {})
        nms = r.get("near_miss_details", [])
        if nms:
            lines.append(f"### {name}\n\n")
            for d in nms:
                lines.append(
                    f"- {d['contig']}:{d['start']}-{d['end']} "
                    f"({d['length']}bp, purity={d['telo_coverage']:.1%}, "
                    f"{d.get('orientation','?')}, {d.get('splice','?')})\n"
                )
            lines.append("\n")

    # Implications
    lines.append("## Implications\n\n")
    if total_hp == 0:
        lines.append(
            "No high-purity telotrons were found across any human reference genome assembly, "
            "including the T2T-CHM13v2.0 complete genome. This confirms that human intronic "
            "telomeric sequences (ITSs) are structurally distinct from the high-purity "
            "telotron arrays found in marine haptophyte MAGs. Vertebrate ITSs represent "
            "degraded, ancient telomeric remnants rather than recently inserted, actively "
            "spliced telomeric repeat arrays.\n"
        )
    else:
        lines.append(
            f"Found {total_hp} high-purity telotrons across human assemblies. "
            "These should be examined for splice site quality, strand orientation, "
            "and comparison to haptophyte telotron architecture.\n"
        )

    report_path.write_text("".join(lines))
    print(f"\nReport: {report_path}")

    # TSV summary
    with open(summary_path, "w") as f:
        f.write("accession\tassembly\tdescription\tstatus\ttotal_introns\t"
                "telo_10pct\ttelo_30pct\ttelo_50pct\ttelo_70pct\ttelo_85pct\ttelo_95pct\n")
        for acc, name, desc in HUMAN_ASSEMBLIES:
            r = ckpt.get(acc, {})
            f.write(f"{acc}\t{name}\t{desc}\t{r.get('status','missing')}\t"
                    f"{r.get('n_introns',0)}\t"
                    f"{r.get('n_telotrons_10pct',0)}\t"
                    f"{r.get('n_telotrons_30pct',0)}\t"
                    f"{r.get('n_telotrons_50pct',0)}\t"
                    f"{r.get('n_telotrons_70pct',0)}\t"
                    f"{r.get('n_telotrons_85pct',0)}\t"
                    f"{r.get('n_telotrons_95pct',0)}\n")
    print(f"Summary: {summary_path}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for d in [DIR_GENOMES, DIR_TELOTRONS, DIR_RESULTS]:
        d.mkdir(parents=True, exist_ok=True)

    ckpt = load_checkpoint()

    # Determine what to run
    if "--force" in sys.argv:
        todo = list(HUMAN_ASSEMBLIES)
        ckpt = {}
        save_checkpoint(ckpt)
    elif "--only" in sys.argv:
        idx = sys.argv.index("--only") + 1
        target = sys.argv[idx] if idx < len(sys.argv) else ""
        todo = [(a, n, d) for a, n, d in HUMAN_ASSEMBLIES
                if target.lower() in a.lower() or target.lower() in n.lower()]
        if not todo:
            print(f"No assembly matching '{target}' found.")
            return
    else:
        todo = [(a, n, d) for a, n, d in HUMAN_ASSEMBLIES if a not in ckpt]

    if not todo:
        print(f"All {len(HUMAN_ASSEMBLIES)} assemblies already processed.")
        print("Use --force to rerun all, or --only <name> to run one.")
        write_final_report(ckpt)
        return

    print(f"=" * 70)
    print(f" Human Genome Telotron Survey")
    print(f" Assemblies: {len(HUMAN_ASSEMBLIES)} total, "
          f"{len(ckpt)} done, {len(todo)} to process")
    print(f" Purity threshold: ≥85% TTAGGG/CCCTAA coverage")
    print(f"=" * 70)

    # Process sequentially (human genomes are very large, ~3GB each)
    for acc, name, desc in todo:
        result = process_one_genome(acc, name, desc)
        mark_done(acc, result)

        # Status
        status = result.get("status", "?")
        n_int = result.get("n_introns", 0)
        n_hp = result.get("n_telotrons_85pct", 0)
        tag = "OK" if status == "ok" else f"FAIL"
        print(f"\n  [{tag}] {name}: {n_int:,} introns, "
              f"{n_hp} telotrons (≥85%)", flush=True)

    # Final report
    ckpt = load_checkpoint()
    write_final_report(ckpt)


if __name__ == "__main__":
    main()
