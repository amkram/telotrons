# Telotron survey pipeline: scan eukaryotic introns (NCBI RefSeq + Tara SMAGs)
# for telomeric-repeat arrays, filter, analyze, plot, and package outputs.


import os

configfile: "config.yaml"


TARA_BASE = config["tara_base"]
REFSEQ_URL = config["refseq_url"]
GENBANK_URL = config.get("genbank_url",
    "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt")
THREADS = int(config["threads"])

# Downloaded assemblies: union sentinel = both dirs present. Rules that scan
# genomes (scan_all, extract_telotron_fasta, telotron_ortholog_align, etc.)
# require both so the union is a single input list.
ASSEMBLIES_DONE = ["data/raw/refseq/.done", "data/raw/genbank/.done"]
ENV = config["env"]
REFSEQ_GROUPS = "|".join(config["refseq_groups"])
ACCESSIONS = config.get("accessions") or []
CANONICAL = config.get("canonical_telomere_motifs", {}) or {}
CANONICAL_BY_GENOME = CANONICAL.get("by_genome", {}) or {}
CANONICAL_BY_GROUP = CANONICAL.get("by_group", {}) or {}

# CLI accepts either a comma-separated string or a YAML list.
_motifs = config["telomere_motifs"]
TELOMERE_MOTIFS = ",".join(_motifs) if isinstance(_motifs, list) else _motifs

SCAN_MIN_FRAC = config["scan"]["min_repeat_frac"]
SCAN_MAX_FLANK_FRAC = config["scan"]["max_flank_repeat_frac"]
SCAN_MIN_INTRON_LEN = config["scan"].get("min_intron_len", 30)
FILTER_MIN_FRAC = config["filter"]["min_repeat_frac"]
FILTER_BIDIR_FRAC = config["filter"]["bidir_min_repeat_frac"]
FILTER_BIDIR_HITS = config["filter"]["bidir_min_hits"]
# Consolidated "what counts as a telotron" thresholds — see top of config.yaml
# for the full threshold map. These keys default to the script-internal values
# when the config block is missing (backward-compatible).
ARCH_MIN_ARRAY = int(config.get("architecture", {}).get("min_array_bp", 18))
INTERSTITIAL_MIN_ARRAY_LEN = int(config.get("interstitial", {}).get("min_array_len", 30))
INTERSTITIAL_MIN_ITS_UNITS = int(config.get("interstitial", {}).get("min_its_units", 4))
FILTER_FLAGS = " ".join(
    flag
    for flag, on in [
        ("--collapse-unique-loci", config["filter"]["collapse_unique_loci"]),
        ("--require-terminal-motif-match", config["filter"]["require_terminal_motif_match"]),
        ("--require-canonical-splice", config["filter"]["require_canonical_splice"]),
    ]
    if on
)

MEME_ENV = "envs/meme.yaml"

# Sentinels (path = ".../<name>/.done" or the canonical output file).
TERT_DEEP_HOMOLOGY_TSV          = "work/results/tert_deep_homology/confirmed_tert.tsv"
TELOTRON_ORTHO_SENTINEL         = "work/results/telotron_orthologs/.done"
TELOTRON_ORTHO_LOCI_PDF         = "work/results/figures/telotron_ortholog_loci.pdf"
TELOTRON_ORTHO_V2               = "work/results/telotron_orthologs_v2"
TELOTRON_ORTHO_TEXT             = TELOTRON_ORTHO_V2 + "/locus_text/.textdump.done"

# Per-locus 5'/3' end orientation. Hybrid (GC, CG) = G-rich one end + C-rich
# the other (strand switch within an intron). Interstitial arrays are
# single-direction (GG or CC only). Non-telotron control introns aren't stratified.
TELOTRON_CAT = {
    "GG": ("architecture", "GT-F-AG"),
    "CC": ("architecture", "GT-R-AG"),
    "GC": ("architecture", "GT-F-R-AG,GT-F-linker-R-AG"),
    "CG": ("architecture", "GT-R-linker-F-AG"),
}
INTERSTITIAL_CAT = {c: ("cat5_3", c) for c in ("GG", "GC", "CG", "CC")}


# Default target = core survey through `package`. Sub-arms are invoked by
# rule name (e.g. `snakemake analysis_arm`, `snakemake find_tert`,
# `snakemake linker_analysis`). `analysis_arm` is the only surviving aggregate.
rule all:
    input:
        "work/results/telotron_pipeline_outputs.zip",


# Build the unified genome manifest from NCBI GenBank + RefSeq assembly
# summaries + Tara SMAGs index. Every ANNOTATED eukaryote is included; RefSeq
# (GCF_) wins where a genome has both a GCF_ and a GCA_ record (RefSeq
# annotations are curated). If config["accessions"] is non-empty, subset to
# just those genome_ids.
rule manifests:
    output:
        refseq="work/manifests/refseq_euk.tsv",
        genbank="work/manifests/genbank_euk.tsv",
        tara="work/manifests/tara_mags.tsv",
        all="work/manifests/all_genomes.tsv",
    params:
        accessions=ACCESSIONS,
    shell:
        r"""
        mkdir -p manifests
        # RefSeq (curated) — group (col 25), ftp_path (col 20).
        curl -L --fail -s {REFSEQ_URL} \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1 || /^#/ {{next}} \
              $25 ~ /^({REFSEQ_GROUPS})$/ && $20!="na" \
              {{print $1,$8,$25,$20,"refseq"}}' \
          > {output.refseq}

        # GenBank (all annotated eukaryotes) — annotation_provider (col 33) != "na"
        # ensures a GFF exists. Exclude rows that have a paired GCF (col 18): those
        # are already represented in the RefSeq manifest.
        curl -L --fail -s {GENBANK_URL} \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1 || /^#/ {{next}} \
              $25 ~ /^({REFSEQ_GROUPS})$/ && $20!="na" && $33!="na" && ($18=="na" || $18=="") \
              {{print $1,$8,$25,$20,"genbank"}}' \
          > {output.genbank}

        curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
          | tar -tzf - \
          | sed -n 's#.*/\(TARA_.*MAG_[0-9][0-9]*\)\.gmove\.gff#\1#p' \
          | sort -u \
          | awk 'BEGIN{{OFS="\t"}} {{print $1,"Tara Oceans MAG","tara","genoscope","tara"}}' \
          > {output.tara}

        printf "genome_id\torganism\tgroup\tftp_path\tsource\n" > {output.all}
        cat {output.refseq} {output.genbank} {output.tara} >> {output.all}

        # Optional accession whitelist: keep header + matching genome_ids.
        if [ -n "{params.accessions}" ]; then
            printf '%s\n' {params.accessions} > work/manifests/.accessions.txt
            for f in {output.all} {output.refseq} {output.genbank} {output.tara}; do
                awk 'NR==FNR{{keep[$1];next}} FNR==1 || $1 in keep' \
                    work/manifests/.accessions.txt "$f" > "$f".subset
                mv "$f".subset "$f"
            done
        fi
        """


# Tara SMAGs ship as two monolithic ~50 GB tarballs (contigs + GFF). Stream-extract only
# the MAGs listed in work/manifests/tara_mags.tsv so the test set doesn't blow up disk.
rule tara_archives:
    input:
        "work/manifests/tara_mags.tsv",
    output:
        fna="data/raw/tara/.fna.done",
        gff="data/raw/tara/.gff.done",
    shell:
        r"""
        mkdir -p data/raw/tara
        mapfile -t MAGS < <(awk -F'\t' '{{print $1}}' {input})
        if [ "${{#MAGS[@]}}" -eq 0 ]; then
            touch {output.fna} {output.gff}; exit 0
        fi
        FNA_PATS=(); GFF_PATS=()
        for m in "${{MAGS[@]}}"; do
            FNA_PATS+=("Contigs/${{m}}.*")
            GFF_PATS+=("Genes/GFF/${{m}}.gmove.*")
        done
        # Stream both tarballs concurrently — independent, no shared output paths.
        ( curl -L --fail -s {TARA_BASE}/SMAGs_contigs_individual.fna.tar.gz \
            | tar --wildcards -xzf - -C data/raw/tara "${{FNA_PATS[@]}}" ) &
        ( curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
            | tar --wildcards -xzf - -C data/raw/tara "${{GFF_PATS[@]}}" ) &
        wait
        touch {output.fna} {output.gff}
        """


# Derive per-genome FNA and GFF URLs from every assembly's FTP path. Union
# writes: acc\tsource\tfna_url\tgff_url. Source-tagged so downloads land under
# Fan out per-genome curls. Derives URLs from refseq_euk + genbank_euk
# manifests inline (the standalone `assembly_urls` rule was inlined 2026-07-14).
# NCBI rate-limits aggressive parallel pulls, so cap the inner pool at 8
# regardless of THREADS. curl retries handle transient 503s; `exit 255` aborts
# xargs on persistent failure. GCF_ accessions land in data/raw/refseq/
# (unchanged); GCA_-only in data/raw/genbank/.
rule download_assemblies:
    input:
        refseq="work/manifests/refseq_euk.tsv",
        genbank="work/manifests/genbank_euk.tsv",
    output:
        refseq_done=touch("data/raw/refseq/.done"),
        genbank_done=touch("data/raw/genbank/.done"),
    threads: 8
    shell:
        r"""
        mkdir -p data/raw/refseq data/raw/genbank
        cat {input.refseq} {input.genbank} \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} {{
              sub(/\/$/, "", $4);
              n=split($4,a,"/"); base=a[n];
              root=($5=="refseq"?"data/raw/refseq":"data/raw/genbank");
              print root"\t"$1"\t"$4"/"base"_genomic.fna.gz";
              print root"\t"$1"\t"$4"/"base"_genomic.gff.gz"
            }}' \
          | xargs -P 8 -n3 sh -c \
              'mkdir -p $0/$1 && curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s "$2" -o $0/$1/$(basename "$2") || exit 255'
        touch {output.refseq_done} {output.genbank_done}
        """


# Canonical per-genome telomere motif TSV, built once from config.
# Consumed by scan_one_genome (and downstream length/interstitial arms).
rule canonical_motifs:
    input:
        manifest="work/manifests/all_genomes.tsv",
    output:
        "work/manifests/canonical_motifs.tsv",
    params:
        by_genome=CANONICAL_BY_GENOME,
        by_group=CANONICAL_BY_GROUP,
    run:
        import csv
        with open(input.manifest) as fh, open(output[0], "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["genome_id", "motif"])
            for r in csv.DictReader(fh, delimiter="\t"):
                m = params.by_genome.get(r["genome_id"]) or params.by_group.get(r["group"], "")
                w.writerow([r["genome_id"], m or ""])


# Checkpoint: emit the genome list from the manifest so downstream rules can
# fan out one sbatch job per genome. Snakemake re-plans the DAG once this
# completes.
checkpoint scan_list:
    input:
        manifest="work/manifests/all_genomes.tsv",
    output:
        "work/results/scan/scan_list.txt",
    run:
        import csv
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(input.manifest) as fh, open(output[0], "w") as out:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.write(r["genome_id"] + "\n")


# Per-genome scan — one sbatch job per genome under slurm; slurm scatters
# them across nodes. 4 CPUs / 8 GB mem per genome by default; big vertebrate
# genomes get bumped via set-resources in profiles/slurm/config.yaml.
# Intron motif coverage, splice signals, orientation, distance-to-end.
# Broad sieve (low --min-repeat-frac) so linker-architecture loci survive;
# strict cutoffs live in filter_final. --max-flank-repeat-frac 0.25 +
# --min-intron-len 30 reject misannotated subtelomere introns and 1-bp
# degenerate GFF intron annotations.
rule scan_one_genome:
    input:
        manifest="work/manifests/all_genomes.tsv",
        canonical="work/manifests/canonical_motifs.tsv",
        tara=["data/raw/tara/.fna.done", "data/raw/tara/.gff.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        loci=temp("work/results/scan/per_genome/{gid}.loci.tsv"),
        introns=temp("work/results/scan/per_genome/{gid}.introns.tsv"),
        summary=temp("work/results/scan/per_genome/{gid}.summary.tsv"),
    threads: 4
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/scan/per_genome
        python scripts/scan_telotrons.py \
            --manifest {input.manifest} \
            --single-genome {wildcards.gid} \
            --canonical-motifs {input.canonical} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --min-repeat-frac {SCAN_MIN_FRAC} --max-flank-repeat-frac {SCAN_MAX_FLANK_FRAC} \
            --min-intron-len {SCAN_MIN_INTRON_LEN} \
            --bidir-min-hits {FILTER_BIDIR_HITS} \
            --threads {threads} \
            --loci {output.loci} --introns {output.introns} --summary {output.summary}
        """


# RERUNS: --forcerun on the wildcarded WORKER (scan_one_genome), not on the
# aggregator (scan_all). Forcing scan_all alone only re-concats stale shards.
def _scan_per_genome_outputs(wildcards):
    """Post-checkpoint aggregation: list every per-genome TSV to expand over."""
    with open(checkpoints.scan_list.get().output[0]) as fh:
        gids = [line.strip() for line in fh if line.strip()]
    return {
        "loci":    expand("work/results/scan/per_genome/{gid}.loci.tsv",    gid=gids),
        "introns": expand("work/results/scan/per_genome/{gid}.introns.tsv", gid=gids),
        "summary": expand("work/results/scan/per_genome/{gid}.summary.tsv", gid=gids),
    }


# Aggregator: concat the per-genome TSVs into the pipeline-level TSVs.
# Runs on the driver node; cheap (< 1 min for ~11k genomes). Uses a Python
# run-block so snakemake's dry-run doesn't try to expand {input.loci[0]}
# before the scan_list checkpoint materialises.
rule scan_all:
    input:
        unpack(_scan_per_genome_outputs),
    output:
        loci="work/results/all_telotron_loci.tsv",
        introns="work/results/all_introns_scanned.tsv",
        summary="work/results/all_species_raw_summary.tsv",
    threads: 1
    run:
        os.makedirs("work/results", exist_ok=True)
        for key, out in [("loci", output.loci),
                         ("introns", output.introns),
                         ("summary", output.summary)]:
            parts = input[key]
            # Find the first non-empty shard for the header — an early empty
            # temp() file (from a genome with no introns) would otherwise
            # produce a headerless merged TSV and misalign every downstream
            # pandas read WITHOUT raising.
            header = None
            for p in parts:
                try:
                    with open(p) as fin:
                        header = fin.readline()
                        if header:
                            break
                except FileNotFoundError:
                    continue
            with open(out, "w") as fout:
                if header:
                    fout.write(header)
                for p in parts:
                    try:
                        with open(p) as fin:
                            first = fin.readline()
                            # Only strip if it looks like the header we wrote.
                            if first and first != header:
                                fout.write(first)
                            for line in fin:
                                fout.write(line)
                    except FileNotFoundError:
                        continue


# ── Broad any-repeat scan ─────────────────────────────────────────────────
# Sibling of scan_all: finds introns dominated by ANY short tandem repeat
# unit (k=4..10 by default), not just the configured telomere motifs.
# Flags matches to config telomere motifs for downstream curation. Reuses
# the scan_list checkpoint so slurm fans out one sbatch per genome.
REPEAT_SCAN = config.get("repeat_scan", {}) or {}
REPEAT_K_MIN = int(REPEAT_SCAN.get("k_min", 4))
REPEAT_K_MAX = int(REPEAT_SCAN.get("k_max", 10))
REPEAT_MIN_FRAC = float(REPEAT_SCAN.get("min_frac", 0.5))
REPEAT_MIN_INTRON_LEN = int(REPEAT_SCAN.get("min_intron_len", 30))


rule scan_repeat_introns_one:
    input:
        manifest="work/manifests/all_genomes.tsv",
        tara=["data/raw/tara/.fna.done", "data/raw/tara/.gff.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        tsv=temp("work/results/repeat_scan/per_genome/{gid}.tsv"),
    threads: 4
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/repeat_scan/per_genome
        python scripts/scan_repeat_introns.py \
            --manifest {input.manifest} \
            --single-genome {wildcards.gid} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --telomere-motifs {TELOMERE_MOTIFS} \
            --k-min {REPEAT_K_MIN} --k-max {REPEAT_K_MAX} \
            --min-frac {REPEAT_MIN_FRAC} --min-intron-len {REPEAT_MIN_INTRON_LEN} \
            --threads {threads} \
            --out {output.tsv}
        """


def _repeat_scan_per_genome_outputs(wildcards):
    with open(checkpoints.scan_list.get().output[0]) as fh:
        gids = [line.strip() for line in fh if line.strip()]
    return expand("work/results/repeat_scan/per_genome/{gid}.tsv", gid=gids)


# Aggregator: concat per-genome any-repeat TSVs. Uses the same header-safe
# concat as scan_all to survive empty shards from zero-intron genomes.
rule scan_repeat_introns:
    input:
        parts=_repeat_scan_per_genome_outputs,
    output:
        tsv="work/results/all_repeat_introns.tsv",
    threads: 1
    run:
        os.makedirs("work/results", exist_ok=True)
        header = None
        for p in input.parts:
            try:
                with open(p) as fin:
                    header = fin.readline()
                    if header:
                        break
            except FileNotFoundError:
                continue
        with open(output.tsv, "w") as fout:
            if header:
                fout.write(header)
            for p in input.parts:
                try:
                    with open(p) as fin:
                        first = fin.readline()
                        if first and first != header:
                            fout.write(first)
                        for line in fin:
                            fout.write(line)
                except FileNotFoundError:
                    continue


# Tighten candidates → final set. --require-terminal-motif-match enforces that
# the intronic motif matches the genome's actual telomere motif (kills cross-motif noise).
rule filter_final:
    input:
        loci="work/results/all_telotron_loci.tsv",
        introns="work/results/all_introns_scanned.tsv",
        summary="work/results/all_species_raw_summary.tsv",
    output:
        final="work/results/final_telotron_set.tsv",
        species="work/results/final_species_summary.tsv",
        neg="work/results/final_negative_controls.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/filter_final_set.py \
            --loci {input.loci} --introns {input.introns} --summary {input.summary} \
            --min-repeat-frac {FILTER_MIN_FRAC} \
            --bidir-min-repeat-frac {FILTER_BIDIR_FRAC} \
            --bidir-min-hits {FILTER_BIDIR_HITS} \
            {FILTER_FLAGS} \
            --final {output.final} --species {output.species} --negatives {output.neg}
        """


# Collapse doubly-assembled exons / close gene copies by aligning per-locus
# 250 bp upstream + 250 bp downstream flanks within each species (all-vs-all
# blastn). Two loci are duplicates if either flank shares >=100 bp at >=90%
# identity; the longest-intron member of each component is kept.
rule dedup_telotrons:
    input:
        final="work/results/final_telotron_set.tsv",
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        final="work/results/final_telotron_set_dedup.tsv",
        log="work/results/dedup_log.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/dedup_telotrons.py \
            --final {input.final} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --out-final {output.final} --out-log {output.log} \
            --threads {threads}
        """


# Classify each (deduplicated) telotron locus into one of:
# GT-F-AG, GT-R-AG, GT-F-R-AG, GT-F-linker-R-AG, GT-R-linker-F-AG, Other.
# Emits the long-form boundary-kmer-by-architecture TSV used by the per-species
# figures (donor / acceptor / linker_left / linker_right 6-mers).
rule classify_architecture:
    input:
        final="work/results/final_telotron_set_dedup.tsv",
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        loci="work/results/final_telotron_set_architecture.tsv",
        kmers="work/results/boundary_kmers_by_architecture.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/classify_telotron_architecture.py \
            --final {input.final} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --min-array {ARCH_MIN_ARRAY} \
            --out-loci {output.loci} --out-kmers {output.kmers}
        """


# Species-level confident-bearer set. A species is admitted when it either has
# >=3 telotrons passing filter_final OR at least one bidirectional architecture
# (GT-F-R-AG or a linker variant — a distinctive telomerase-mediated signature).
# All downstream analyses (gene-class, expression, nucleosome, ortholog panels)
# key off this file so new bearer species flow through automatically.
_CONFIDENT = config.get("confident_species", {}) or {}
CONFIDENT_MIN_N = int(_CONFIDENT.get("min_n", 3))
CONFIDENT_MIN_BIDIR = int(_CONFIDENT.get("min_bidir", 2))

rule confident_species:
    input:
        arch="work/results/final_telotron_set_architecture.tsv",
        manifest="work/manifests/all_genomes.tsv",
    output:
        "work/results/confident_species.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/confident_species.py \
            --arch {input.arch} --manifest {input.manifest} \
            --min-n {CONFIDENT_MIN_N} --min-bidir {CONFIDENT_MIN_BIDIR} \
            --out {output}
        """


# Boundary k-mer enrichment, length-matched distance-to-end test, architecture summary.
rule analyze:
    input:
        final="work/results/final_telotron_set_dedup.tsv",
        introns="work/results/all_introns_scanned.tsv",
        species="work/results/final_species_summary.tsv",
        neg="work/results/final_negative_controls.tsv",
    output:
        kmers="work/results/boundary_kmer_enrichment.tsv",
        dist="work/results/distance_to_end.tsv",
        arch="work/results/architecture_summary.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/analyze_telotrons.py \
            --final {input.final} --introns {input.introns} --species {input.species} \
            --boundary-kmers {output.kmers} --distance {output.dist} --architecture {output.arch} \
            --threads {threads}
        """


# Per-species × architecture FASTA + flanked-text extracts. Wildcarded on
# {set} (telotron | non_telotron); one rule, one script — the two extract
# arms (was extract_telotron_fasta + extract_non_telotron_fasta) are the
# same pipeline against different input tables.
# Flanked lines: [LEFT100] [INTRON] [RIGHT100], or for linker archs
# [LEFT100] [ARRAY1] [LINKER] [ARRAY2] [RIGHT100]. One subdir per species,
# one file per architecture.
rule extract_fasta:
    input:
        table=lambda w: ("work/results/final_telotron_set_architecture.tsv"
                          if w.set == "telotron"
                          else "work/results/non_telotron_controls.tsv"),
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        touch("work/results/{set}_fasta/.done"),
    wildcard_constraints:
        set=r"telotron|non_telotron",
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/{wildcards.set}_fasta work/results/{wildcards.set}_flanked
        python scripts/extract_telotron_fasta.py \
            --final {input.table} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --fasta-dir work/results/{wildcards.set}_fasta \
            --flanked-dir work/results/{wildcards.set}_flanked
        touch {output}
        """


# Interstitial telomeric arrays: non-terminal (>=5 kb from contig end),
# non-genic, non-intronic. Exclude (annotated genes ∪ introns) ∪ 6-frame
# ORF mask (work/results/masks/{gid}.bed, ORFs >= 450 nt). Emits 5'/3'
# boundary 6- and 12-mers.
rule find_interstitial_arrays:
    input:
        manifest="work/results/all_species_raw_summary.tsv",
        canonical="work/manifests/canonical_motifs.tsv",
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        "work/results/interstitial_arrays.tsv",
    params:
        mask_dir="work/results/masks",
        min_orf=450,
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        # 6-frame ORF mask: excludes unannotated coding regions from the
        # interstitial-array set. Was a separate `make_unannotated_masks` rule
        # until 2026-07-14; folded inline (one rule, one deliverable).
        mkdir -p {params.mask_dir}
        python scripts/make_unannotated_mask.py \
            --manifest {input.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir {params.mask_dir} \
            --min-orf-nt {params.min_orf}
        python scripts/find_interstitial_arrays.py \
            --manifest {input.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --canonical-motifs {input.canonical} \
            --mask-dir {params.mask_dir} \
            --min-array-len {INTERSTITIAL_MIN_ARRAY_LEN} \
            --out {output} --threads {threads}
        """


# STREME: --dna, w 4-12, 10 motifs, p<0.05, deterministic seed.
STREME_ARGS = "--dna --minw 4 --maxw 12 --nmotifs 10 --thresh 0.05 --seed 1"


# ── TERT by deep homology (miniprot + Pfam domain architecture) ─────────────
# BLAST/tblastn fail on apicomplexan TERT. Field-standard route: intron-aware
# miniprot with apicomplexan TERT seeds, confirmed by RBD (PF12009) UPSTREAM of
# RT (PF00078) -- rules out Eimeria retroelement RTs. Iterative re-seeding from
# confirmed proteins. Toxo gondii = positive control (sister; TERT unannotated).
# Outgroup / control targets (telotron-negative sisters). Union'd with
# work/results/confident_species.tsv at rule runtime so ALL confident bearers
# get searched for TERT, not just a hardcoded panel.
TERT_OUTGROUP_IDS = config.get("tert_outgroup_ids") or [
    "GCF_000006565.2",   # Toxoplasma gondii
    "GCF_000208865.1",   # Neospora caninum
    "GCF_002999335.1",   # Cyclospora cayetanensis
    "GCF_002563875.1",   # Besnoitia besnoiti
    "GCF_002600585.1",   # Cystoisospora suis
]


rule fetch_tert_seeds_hmms:
    output:
        seeds="work/results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="work/results/tert_deep_homology/refs/PF12009.hmm",
        rt="work/results/tert_deep_homology/refs/PF00078.hmm",
    shell:
        r"""
        python scripts/fetch_tert_seeds_hmms.py \
            --outdir work/results/tert_deep_homology/refs
        """


rule find_tert:
    # Kept monolithic (single sbatch, ~16 cpus) — iterate 3 expands the seed
    # set from cross-genome hits between rounds, so per-genome shards would
    # lose the iteration. Target set is small (~15-25 confident bearers +
    # outgroups) so single-node is fine; wall-clock ~ 60 min at 16 cpus.
    # If a shardable variant is ever needed use `find_tert_deep_homology.py
    # --single-genome {gid} --iterate 1` (flag already wired).
    input:
        seeds="work/results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="work/results/tert_deep_homology/refs/PF12009.hmm",
        rt="work/results/tert_deep_homology/refs/PF00078.hmm",
        confident="work/results/confident_species.tsv",
        assemblies=ASSEMBLIES_DONE,
        tara="data/raw/tara/.fna.done",
    output:
        tsv=TERT_DEEP_HOMOLOGY_TSV,
        faa="work/results/tert_deep_homology/all_confirmed.faa",
    params:
        outgroups=",".join(TERT_OUTGROUP_IDS),
    threads: 16
    conda:
        ENV
    shell:
        r"""
        # Target set = every confident bearer species + fixed outgroup panel.
        # Auto-scales with the confident-species set — no hardcoded per-species list.
        python scripts/find_tert_deep_homology.py \
            --confident-species {input.confident} \
            --outgroup-ids {params.outgroups} \
            --seeds {input.seeds} --trbd-hmm {input.trbd} --rt-hmm {input.rt} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir work/results/tert_deep_homology \
            --iterate 3 --threads {threads}
        """


# ── Telotron ortholog alignment ──────────────────────────────────────────────
# Per telotron-host gene: intron-aware miniprot maps the ortholog in each panel
# genome, aligns it in protein space, then DNA-aligns the telotron vs the
# orthologous intron. Resolves fill-vs-create. Frameshifts/stops at the homologous
# junction are flagged (not corrected).
_TELO_ORTHO = config.get("telotron_ortholog", {}) or {}
_TELO_ORTHO_FOCAL_IDS = _TELO_ORTHO.get("focal_ids", []) or []
_TELO_ORTHO_PANEL_IDS = _TELO_ORTHO.get("ortholog_ids", []) or []
_TELO_ORTHO_FOCAL = ",".join(_TELO_ORTHO_FOCAL_IDS)
_TELO_ORTHO_PANEL = ",".join(_TELO_ORTHO_PANEL_IDS)
# panel-split for the canonical locus_text: within-Eimeria = focal genomes,
# outgroup = the telotron-negative sisters (ortholog_ids minus the cross-Eimeria).
_TELO_WITHIN_PANEL = ",".join(_TELO_ORTHO_FOCAL_IDS)
_TELO_OUTGROUP_PANEL = ",".join([x for x in _TELO_ORTHO_PANEL_IDS if x not in set(_TELO_ORTHO_FOCAL_IDS)])


# Telotron-host-gene ortholog alignment: miniprot-map the ortholog in each
# panel genome, protein-align it, DNA-align telotron vs orthologous intron to
# resolve fill-vs-create. Emits the compiled per-locus PDF in the same rule
# (was 2 rules: align + plot; merged 2026-07-14).
# Focal-species list is derived from confident_species.tsv; wrap in a
# checkpoint so the DAG can fan out one sbatch per bearer.
checkpoint telotron_ortho_focal_list:
    input:
        confident="work/results/confident_species.tsv",
    output:
        "work/results/telotron_orthologs/focal_list.txt",
    run:
        import csv, os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(input.confident) as fh, open(output[0], "w") as out:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.write(r["genome_id"] + "\n")


# Per-focal-species ortholog run — one sbatch per bearer, iterates the
# config-driven ortholog panel and produces a per-focal directory of
# alignments/text-dumps under `work/results/telotron_orthologs/per_focal/`.
rule telotron_orthologs_one:
    input:
        final="work/results/final_telotron_set_architecture.tsv",
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        sentinel=touch("work/results/telotron_orthologs/per_focal/{fid}/.done"),
    params:
        panel=_TELO_ORTHO_PANEL,
        min_array=_TELO_ORTHO.get("min_array_bp", 50),
        # See main-rule comment: 0.40 kills ABS-bucket leakage.
        min_id=_TELO_ORTHO.get("min_identity", 0.40),
        tol=_TELO_ORTHO.get("intron_tol", 8),
    threads: 8
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/telotron_orthologs/per_focal/{wildcards.fid}
        python scripts/telotron_ortholog_align.py \
            --final {input.final} \
            --focal-ids {wildcards.fid} \
            --ortholog-ids {params.panel} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir work/results/telotron_orthologs/per_focal/{wildcards.fid} \
            --min-array-bp {params.min_array} --min-ident {params.min_id} \
            --intron-tol {params.tol} --threads {threads}
        """


def _telotron_ortho_focal_outputs(wildcards):
    with open(checkpoints.telotron_ortho_focal_list.get().output[0]) as fh:
        fids = [line.strip() for line in fh if line.strip()]
    return expand("work/results/telotron_orthologs/per_focal/{fid}/.done", fid=fids)


# Aggregator: merge every per-focal output into the flat directory layout
# downstream rules expect (`work/results/telotron_orthologs/<...>`), then
# render the combined PDF.
rule telotron_orthologs:
    input:
        per_focal=_telotron_ortho_focal_outputs,
    output:
        sentinel=touch(TELOTRON_ORTHO_SENTINEL),
        pdf=TELOTRON_ORTHO_LOCI_PDF,
    threads: 2
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/telotron_orthologs work/results/figures
        # Flatten per-focal subdirs into the top-level dir. Each per_focal
        # subdir is disjoint (keyed by focal genome_id inside the filenames),
        # so a plain cp -rn preserves outputs from prior runs.
        for d in work/results/telotron_orthologs/per_focal/*/; do
            [ -d "$d" ] || continue
            cp -rn "$d"/* work/results/telotron_orthologs/ 2>/dev/null || true
        done
        python scripts/plot_telotron_ortholog_loci.py \
            --ortho-dir work/results/telotron_orthologs \
            --out {output.pdf} --window 28 --max-rows 12
        """


# Definitely-non-telotron intron control set: introns with telomeric_frac < 10%
# from the same positive species, sampled per genome. Produced TSV mirrors the
# architecture-table schema so extract_telotron_fasta.py can ingest it.
rule build_non_telotron_controls:
    input:
        introns="work/results/all_introns_scanned.tsv",
        final="work/results/final_telotron_set.tsv",
    output:
        tsv="work/results/non_telotron_controls.tsv",
    conda: ENV
    shell:
        r"""
        python scripts/build_non_telotron_controls.py \
            --introns {input.introns} --final {input.final} \
            --out {output.tsv} --max-frac 0.01 \
            --all-genomes
        """


# Bundle all final TSVs + figures + manifest into a single zip for sharing.
PACKAGE_INPUTS = [
    "work/results/final_telotron_set.tsv",
    "work/results/final_telotron_set_dedup.tsv",
    "work/results/final_telotron_set_architecture.tsv",
    "work/results/final_species_summary.tsv",
    "work/results/final_negative_controls.tsv",
    "work/results/confident_species.tsv",
    "work/results/boundary_kmer_enrichment.tsv",
    "work/results/boundary_kmers_by_architecture.tsv",
    "work/results/distance_to_end.tsv",
    "work/results/architecture_summary.tsv",
    "work/results/dedup_log.tsv",
    "work/results/interstitial_arrays.tsv",
    "work/manifests/all_genomes.tsv",
]


rule package:
    input:
        PACKAGE_INPUTS,
    output:
        "work/results/telotron_pipeline_outputs.zip",
    shell:
        "zip -qj {output} {input}"


# ── Architecture / linker analyses ─────────────────────────────────────────
# Downstream of telotron_orthologs. Textdump rule (below) produces the
# per-locus text-view layout that linker_analysis + mask_telotron_arrays
# consume. Invoke individually — no aggregate target.

V2_LOCUS_TEXT = TELOTRON_ORTHO_TEXT   # textdump sentinel


# Per-locus text-view of the ortholog alignments, panel-split. The linker
# and mask rules read this v2 directory layout directly (was killed in the
# 4th deslop then restored 2026-07-14 when both downstream rules turned
# out to depend on it).
rule telotron_ortholog_textdump:
    input:
        TELOTRON_ORTHO_SENTINEL,
    output:
        touch(TELOTRON_ORTHO_TEXT),
    params:
        v2=TELOTRON_ORTHO_V2,
        within=_TELO_WITHIN_PANEL,
        outgroup=_TELO_OUTGROUP_PANEL,
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/telotron_ortholog_textdump.py \
            --ortho-dir work/results/telotron_orthologs --threads {threads} \
            --panel-ids {params.within} \
            --out-dir {params.v2}/locus_text/within_eimeria
        python scripts/telotron_ortholog_textdump.py \
            --ortho-dir work/results/telotron_orthologs --threads {threads} \
            --panel-ids {params.outgroup} \
            --out-dir {params.v2}/locus_text/outgroup
        """


# Linker analysis: segment each telotron into arrays + linkers, then cluster
# the linker sequences by 7-mer Jaccard to quantify cross-locus recurrence.
# One DAG node covering both steps (previously two rules).
rule linker_analysis:
    input:
        V2_LOCUS_TEXT,
    output:
        seg="work/results/mechanism_deepdive/linker_segmentation.tsv",
        seg_jsonl="work/results/mechanism_deepdive/linker_segmentation.jsonl",
        arch="work/results/mechanism_deepdive/architecture_per_locus.tsv",
        clusters="work/results/mechanism_deepdive/linker_clusters.tsv",
        recurrence="work/results/mechanism_deepdive/linker_recurrence_summary.txt",
    conda:
        ENV
    shell:
        r"""
        python scripts/linker_segmentation.py
        python scripts/cluster_linkers.py
        """


rule mask_telotron_arrays:
    input:
        V2_LOCUS_TEXT,
    output:
        "work/results/telotron_masked/telotron_masked.fasta",
        "work/results/telotron_masked/telotron_masked.msa.fasta",
        "work/results/telotron_masked/telotron_masked_msa.txt",
        "work/results/telotron_masked/telotron_masked_msa.html",
    conda:
        ENV
    shell:
        r"""
        python scripts/mask_telotron_arrays.py
        """


# ── Assembly-based telomere boundary detection ──────────────────────────────
# Field-standard motif sliding-window scan on assembled long-read contigs.
# Reports the most-internal motif position contiguous with each contig end.
# Config block `telomere_boundaries` sets per-species assembly + motif.
_TELOBND = config.get("telomere_boundaries", {}) or {}
_TELOBND_ASM = _TELOBND.get("assemblies", {}) or {}


# ── Long-read data + Telogator2 allele-specific telomere-length analysis ────
# data/raw/longread/manifest.tsv is the source of truth. Telogator2 consumes
# the nanopore subset, anchored to the custom Eimeria subtelomere reference at
# work/results/telogator2_ref/. Replaces the retired TARPON pipeline.
# Targets: `longread_data` (download all), `telogator2` (run on every ONT run).

import csv as _csv

LR_MANIFEST = "data/raw/longread/manifest.tsv"

def _read_lr_manifest():
    import os
    if not os.path.exists(LR_MANIFEST):
        return []
    with open(LR_MANIFEST) as fh:
        return list(_csv.DictReader(fh, delimiter="\t"))

LR_RUNS = _read_lr_manifest()
TG2_RUNS = [r for r in LR_RUNS if r["platform"] == "nanopore"]

# accession -> list of local fastq paths (handles multi-file runs e.g. SRR24971026)
_TG2_FILES = {}
for _r in TG2_RUNS:
    _TG2_FILES.setdefault(_r["accession"], []).append("data/raw/longread/" + _r["fastq_local"])

# accession -> species (one row per accession is fine, multi-file runs share species)
_TG2_SPECIES = {r["accession"]: r["species"] for r in TG2_RUNS}

# URL lookup keyed by local destination path
_LR_URL_BY_DEST = {("data/raw/longread/" + r["fastq_local"]): ("https://" + r["fastq_url"]) for r in LR_RUNS}




# ── Downstream analysis arm (was hand-run; now wired) ───────────────────────
# These rules cover the previously DAG-disconnected analysis scripts whose inputs
# ARE regenerable from the core pipeline (the architecture table, the controls
# table, the RefSeq GFFs, the nucleosome manifests). They write to the script's
# own output paths (mostly analysis/ + work/results/...). Scripts that depend on
# non-regenerable external data (ONT/Hi-C/RNA-seq reads, the work/old/ archive,
# hand-built telogator2 refs) are NOT wired — they are documented under the
# `analysis_external_inputs` note below, because a rule with an input no other
# rule can produce would only ever error.

ARCH_TSV = "work/results/final_telotron_set_architecture.tsv"

# RNA-seq gene-coverage arm (config-driven; SRA accessions in config["rnaseq"]).
_RNASEQ = config.get("rnaseq", {}) or {}
_RNASEQ_SP = _RNASEQ.get("species", {}) or {}
RNASEQ_THREADS = int(_RNASEQ.get("threads", THREADS))
RNASEQ_COV = {sp: f"work/results/rnaseq/{sp}_gene_cov.tsv" for sp in _RNASEQ_SP}


# Per-species RNA-seq -> gene depth-sum TSV (samtools bedcov). Replaces the hand-run
# work/necatrix_rnaseq/run_pipeline.sh and the ephemeral /tmp/eten_gene_cov.tsv.
rule rnaseq_gene_coverage:
    input:
        genome=lambda w: _RNASEQ_SP[w.species]["genome"],
        gff=lambda w: _RNASEQ_SP[w.species]["gff"],
        assemblies=ASSEMBLIES_DONE,
    output:
        "work/results/rnaseq/{species}_gene_cov.tsv",
    params:
        srr=lambda w: ",".join(_RNASEQ_SP[w.species]["srr"]),
        maxsize=lambda w: _RNASEQ_SP[w.species].get("max_size", "30g"),
    threads: RNASEQ_THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/rnaseq_gene_coverage.py \
            --species {wildcards.species} \
            --genome {input.genome} --gff {input.gff} \
            --srr {params.srr} --max-size {params.maxsize} \
            --threads {threads} --out {output}
        """


# Insertion-site composition / periodicity panel + within-gene control.
# One rule: (1) build per-locus flank FASTAs + control FASTAs; (2) compute
# composition / 10-bp WW periodicity / CpG panel with telomere-masked flanks
# and BH-FDR; (3) run the within-gene control (telotron introns vs same-gene
# sibling introns vs random non-host, separates gene-class confound from
# local targeting). The old NuPoP occupancy arm was retired 2026-07-14
# (linker interpretation was artifact per memory
# telotron_nucleosome_nupop_2026-06-08); the signals here are the ones that
# held up.
rule nucleosome_analysis:
    input:
        arch=ARCH_TSV,
        controls="work/results/non_telotron_controls.tsv",
        tara=["data/raw/tara/.fna.done"],
        assemblies=ASSEMBLIES_DONE,
    output:
        summary="work/results/nucleosome/nucleosome_feature_summary.png",
        withingene="work/results/nucleosome/withingene_control.txt",
        manifest="work/results/nucleosome/manifest.tsv",
        control_manifest="work/results/nucleosome/control_manifest.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/nucleosome_inputs.py \
            --table {input.arch} --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --out work/results/nucleosome
        python scripts/nucleosome_control_inputs.py \
            --controls {input.controls} --telo-manifest {output.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara --out work/results/nucleosome
        python scripts/nucleosome_features.py
        python scripts/nucleosome_withingene.py | tee {output.withingene}
        """




# Length distribution by (lineage x architecture) — BH-FDR corrected, single-MAG
# caveat in the output. Per-arm burst-length subfigure is emitted from the same
# TSVs by the same script.
rule length_distribution_by_arch:
    input:
        ARCH_TSV,
    output:
        summary="work/results/length_distribution/summary_by_lineage_arch.tsv",
        linker="work/results/length_distribution/linker_vs_nolinker.tsv",
        per_arch="work/results/length_distribution/psw_vs_eimeria_per_arch.tsv",
        hist="work/results/length_distribution/length_hist_by_lineage_arch.pdf",
        per_arm="work/results/length_distribution/per_arm_burst_length.pdf",
    conda:
        ENV
    shell:
        r"""
        TELOTRON_ROOT="$(pwd)" python scripts/length_distribution_by_arch_lineage.py
        """


# Expression-vs-telotron-presence figure. Gene-level coverage comes from the
# reproducible rnaseq_gene_coverage rule; the optional locus-level panel reads
# data/raw/rnaseq_splice_2026/per_locus_counts.tsv (manual build) and is skipped
# gracefully when absent. One consolidated script (was 2: presence + final).
rule telotron_expr_figures:
    input:
        eten="work/results/rnaseq/eten_gene_cov.tsv",
        necatrix="work/results/rnaseq/necatrix_gene_cov.tsv",
        arch=ARCH_TSV,
        assemblies=ASSEMBLIES_DONE,
    output:
        "work/results/figures/telotron_expression.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/telotron_expression.py
        """


# Host-gene-class dissection (host vs disjoint non-host intron-count + gene length;
# per-intron rate vs gene intron count + 5'->3' position). One consolidated
# analysis for the intron-rich long-gene bias (was 2 rules: telotron_gene_bias
# + telotron_per_intron, merged 2026-07-14).
rule telotron_gene_class:
    input:
        arch=ARCH_TSV,
        confident="work/results/confident_species.tsv",
        assemblies=ASSEMBLIES_DONE,
    output:
        "work/results/figures/telotron_gene_class.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/telotron_gene_class.py
        """


# Aggregate target for the wired analysis arm.
rule analysis_arm:
    input:
        rules.nucleosome_analysis.output,
        rules.telotron_gene_class.output,
        rules.length_distribution_by_arch.output,
        rules.telotron_expr_figures.output,
#
# detect_telomere_boundaries.py is fully CLI-parameterised and driven by the
# telomere_boundaries rule.
