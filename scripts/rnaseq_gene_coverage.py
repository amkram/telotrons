#!/usr/bin/env python3
"""Build a gene-level RNA-seq coverage TSV for one species.

Generalises the per-species E. necatrix recipe (work/necatrix_rnaseq/run_pipeline.sh) so the
expression-vs-telotron arm is reproducible from a Snakemake rule instead of a hand-run shell script
and an ephemeral /tmp output. Steps, per accession:

  prefetch + fasterq-dump (--split-3)  ->  minimap2 -ax splice  ->  samtools sort/merge  ->  samtools bedcov

Uses SPLICE-AWARE alignment (`minimap2 -ax splice -ub --secondary=no`) so exon-junction reads are
placed as gapped alignments instead of soft-clipped by the short-read preset. This matters directly
for the telotron_expression analysis: intron-rich telotron-host genes have many exon-exon junctions,
and short-read alignment would systematically depress their coverage — producing a spurious
"host genes have lower expression" that is really the alignment artifact of the exact confound
(intron count) the downstream OLS tries to control for.

Output is BED6 + a 7th column = summed read depth over the gene interval (samtools bedcov), one row per
`gene` feature in the GFF. Consumers (telotron_expression.py) compute expression = depth_sum / gene_length.

Large intermediates (FASTQ, per-run BAM) are removed after the merged BAM + bedcov are produced unless
--keep-intermediates is given. Idempotent: skips the final TSV if it already exists and is non-empty
(use --force to rebuild).
"""
import argparse
import os
import re
import subprocess
import sys


def log(*a):
    print("[rnaseq_gene_coverage]", *a, file=sys.stderr, flush=True)


def run(cmd, **kw):
    log("$", cmd if isinstance(cmd, str) else " ".join(cmd))
    subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)


def genes_bed_from_gff(gff, out_bed):
    """Write a BED6 of every `gene` feature (0-based start), sorted by the caller."""
    n = 0
    with open(out_bed, "w") as o:
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            if not m:
                continue
            o.write(f"{f[0]}\t{int(f[3]) - 1}\t{f[4]}\t{m.group(1)}\t0\t{f[6]}\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, help="short species key, e.g. eten / necatrix")
    ap.add_argument("--genome", required=True, help="genome FASTA (.fna[.gz])")
    ap.add_argument("--gff", required=True, help="genome GFF (for gene intervals)")
    ap.add_argument("--srr", required=True, help="comma-separated SRA run accession(s), e.g. SRR7153211,SRR7153212")
    ap.add_argument("--out", required=True, help="output gene-coverage TSV")
    ap.add_argument("--workdir", default=None, help="scratch dir for FASTQ/BAM (default: dir of --out)")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--max-size", default="30g", help="prefetch --max-size")
    ap.add_argument("--keep-intermediates", action="store_true")
    ap.add_argument("--force", action="store_true", help="rebuild even if --out exists")
    a = ap.parse_args()

    if os.path.exists(a.out) and os.path.getsize(a.out) > 0 and not a.force:
        log(f"{a.out} exists and is non-empty; skipping (use --force to rebuild)")
        return

    wd = a.workdir or os.path.dirname(os.path.abspath(a.out)) or "."
    scratch = os.path.join(wd, f"_{a.species}_rnaseq")
    os.makedirs(scratch, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    srrs = [s.strip() for s in a.srr.split(",") if s.strip()]

    per_run_bams = []
    for srr in srrs:
        bam = os.path.join(scratch, f"{srr}.bam")
        if os.path.exists(bam) and os.path.getsize(bam) > 0 and not a.force:
            log(f"{bam} exists; reusing")
            per_run_bams.append(bam)
            continue
        log(f"prefetch {srr}")
        run(["prefetch", "-O", scratch, "--max-size", a.max_size, srr])
        log(f"fasterq-dump {srr}")
        run(["fasterq-dump", "--split-3", "-e", str(a.threads),
             "--temp", os.path.join(scratch, "tmp"), "-O", scratch,
             os.path.join(scratch, srr, f"{srr}.sra")])
        r1 = os.path.join(scratch, f"{srr}_1.fastq")
        r2 = os.path.join(scratch, f"{srr}_2.fastq")
        single = os.path.join(scratch, f"{srr}.fastq")
        if os.path.exists(r1) and os.path.exists(r2):
            reads = [r1, r2]
        elif os.path.exists(single):
            reads = [single]
        else:
            raise SystemExit(f"no FASTQ produced for {srr} in {scratch}")
        log(f"minimap2 -ax splice {srr} -> sort")
        # SPLICE-AWARE preset: reads spanning exon-exon junctions get gapped
        # alignments (essential for intron-rich telotron-host genes).
        # `-ub` = unstranded (both-strand splice-site discovery); correct for
        # the standard-Illumina PE dUTP-free libraries in config['rnaseq'].
        # Use `-uf` only for stranded libraries (dUTP/directional). `--secondary=no`
        # keeps bedcov from double-counting multimappers.
        run(f"minimap2 -ax splice -ub --secondary=no -t {a.threads} {a.genome} {' '.join(reads)} "
            f"2>{scratch}/{srr}.mm2.log | samtools sort -@ 4 -m 2G -o {bam} -")
        run(["samtools", "index", bam])
        per_run_bams.append(bam)
        if not a.keep_intermediates:
            for f in reads:
                try:
                    os.remove(f)
                except OSError:
                    pass
            run(f"rm -rf {scratch}/{srr} {scratch}/tmp")

    # Merge per-run BAMs (or use the single one) for bedcov.
    if len(per_run_bams) == 1:
        merged = per_run_bams[0]
    else:
        merged = os.path.join(scratch, f"{a.species}.merged.bam")
        run(["samtools", "merge", "-f", "-@", str(a.threads), merged] + per_run_bams)
        run(["samtools", "index", merged])

    bed = os.path.join(scratch, f"{a.species}.genes.bed")
    bed_sorted = os.path.join(scratch, f"{a.species}.genes.sorted.bed")
    n = genes_bed_from_gff(a.gff, bed)
    log(f"genes: {n}")
    run(f"sort -k1,1 -k2,2n {bed} > {bed_sorted}")
    run(f"samtools bedcov {bed_sorted} {merged} > {a.out}")
    with open(a.out) as fh:
        rows = sum(1 for _ in fh)
    log(f"bedcov done: {rows} genes -> {a.out}")

    if not a.keep_intermediates and len(per_run_bams) > 1:
        # We made a separate merged BAM; the per-run BAMs are now redundant. (When there is a single
        # run, merged IS that BAM, so we keep it — re-bedcov is cheap, re-aligning is not.)
        for b in per_run_bams:
            for f in (b, b + ".bai"):
                try:
                    os.remove(f)
                except OSError:
                    pass


if __name__ == "__main__":
    main()
