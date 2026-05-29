# Telotron survey pipeline: scan eukaryotic introns (NCBI RefSeq + Tara SMAGs)
# for telomeric-repeat arrays, filter, analyze, plot, and package outputs.


import os

configfile: "config.yaml"


TARA_BASE = config["tara_base"]
REFSEQ_URL = config["refseq_url"]
THREADS = int(config["threads"])
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
FILTER_MIN_FRAC = config["filter"]["min_repeat_frac"]
FILTER_BIDIR_FRAC = config["filter"]["bidir_min_repeat_frac"]
FILTER_BIDIR_HITS = config["filter"]["bidir_min_hits"]
FILTER_FLAGS = " ".join(
    flag
    for flag, on in [
        ("--collapse-unique-loci", config["filter"]["collapse_unique_loci"]),
        ("--require-terminal-motif-match", config["filter"]["require_terminal_motif_match"]),
        ("--require-canonical-splice", config["filter"]["require_canonical_splice"]),
    ]
    if on
)

FIGURES = [
    "results/figures/telotron_counts.png",
]
TERMINAL_DENSITY_SENTINEL = "results/figures/terminal_motif/.done"
ARCH_KMER_SENTINEL = "results/figures/boundary_kmers_by_arch/.done"
INTERSTITIAL_KMER_SENTINEL = "results/figures/interstitial_boundary_kmers/.done"
INTERSTITIAL_SPLICE_KMER_SENTINEL = "results/figures/interstitial_boundary_kmers_c_rich_splice_candidates/.done"
INTERSTITIAL_LOGO_SENTINEL = "results/figures/interstitial_boundary_logos/.done"
ARRAY_LEN_DIST_PNG = "results/figures/array_length_distribution.png"
ARRAY_LEN_DIST_MIN40_PNG = "results/figures/array_length_distribution_min40.png"
PIPELINE_STAGES_SENTINEL = "results/figures/pipeline_stages/.done"
EXTRACT_FASTA_SENTINEL = "results/telotron_fasta/.done"
MSA_SENTINEL = "results/msa_regions/.done"
STREME_TELO_SENTINEL = "results/streme/telotrons/streme.xml"
STREME_NONTELO_SENTINEL = "results/streme/non_telo_introns/streme.xml"
STREME_LINKER_SENTINEL = "results/streme/linkers/streme.xml"
STREME_BRANCHPOINT_SENTINEL = "results/streme/branchpoint/.done"
FIMO_BRANCHPOINT_SENTINEL   = "results/fimo/branchpoint/.done"
MEME_ENV = "envs/meme.yaml"
TELOMERASE_DB_SENTINEL = "results/telomerase_db/.blastdb.done"
TELOMERASE_BLAST_SENTINEL = "results/blast_telomerase_vs_genomes/.done"
TELOMERASE_BLAST_REPORT_SENTINEL = "results/blast_telomerase_vs_genomes/.report.done"
TERT_DEEP_HOMOLOGY_TSV = "results/tert_deep_homology/confirmed_tert.tsv"

TELOTRON_SPLICE_LOGO_SENTINEL = "results/figures/telotron_splice_logos/.done"
CTRL_SPLICE_LOGO_SENTINEL    = "results/figures/non_telotron_splice_logos/.done"
CTRL_BOUNDARY_KMER_SENTINEL  = "results/figures/non_telotron_boundary_kmers/.done"
CTRL_FLANKED_SENTINEL        = "results/non_telotron_fasta/.done"
TELOTRON_BOUNDARY_KMER_SENTINEL = "results/figures/telotron_boundary_kmers/.done"
COMPOSITE_KMER_SENTINEL  = "results/figures/composite_boundary_kmers/.done"
COMPOSITE_LOGO_SENTINEL  = "results/figures/composite_boundary_logos/.done"

# Per-locus 5'/3' end orientation categories (GG / GC / CG / CC). Hybrid
# telotrons (GC, CG) arise when a single intron contains a G-rich array on one
# end and a C-rich array on the other — i.e. strand switching within the same
# intron. Interstitial arrays are single-direction (GG or CC only). Non-telotron
# control introns have no telomeric repeat and are not stratified.
TELOTRON_CAT = {
    "GG": ("architecture", "GT-F-AG"),
    "CC": ("architecture", "GT-R-AG"),
    "GC": ("architecture", "GT-F-R-AG,GT-F-linker-R-AG"),
    "CG": ("architecture", "GT-R-linker-F-AG"),
}
INTERSTITIAL_CAT = {
    "GG": ("cat5_3", "GG"),
    "GC": ("cat5_3", "GC"),
    "CG": ("cat5_3", "CG"),
    "CC": ("cat5_3", "CC"),
}
COMPOSITE_5P3P_KMER_SENTINEL = "results/figures/composite_boundary_kmers_5p3p/.done"
COMPOSITE_5P3P_LOGO_SENTINEL = "results/figures/composite_boundary_logos_5p3p/.done"


# Default target: the core survey through the packaged zip. `package` already
# pulls the final TSVs, analysis tables, interstitial arrays, linker-BLAST hits,
# the combined boundary-kmer PDF, and the headline figures, so the default
# `snakemake` reproduces the survey without the exploratory branches below.
# Run those on demand by target name:
#   snakemake sequences        # per-species FASTA extracts + MAFFT MSAs
#   snakemake figures_extra    # boundary-kmer / logo / composite / density figures
#   snakemake motif_discovery  # STREME + FIMO de-novo motif discovery
#   snakemake telomerase_search# telomerase BLAST + TERT deep-homology search
#   snakemake everything       # core + every branch above (the old default)
rule all:
    input:
        "results/telotron_pipeline_outputs.zip",


rule sequences:
    input:
        EXTRACT_FASTA_SENTINEL,
        MSA_SENTINEL,


rule figures_extra:
    input:
        TERMINAL_DENSITY_SENTINEL,
        ARCH_KMER_SENTINEL,
        INTERSTITIAL_KMER_SENTINEL,
        INTERSTITIAL_SPLICE_KMER_SENTINEL,
        INTERSTITIAL_LOGO_SENTINEL,
        ARRAY_LEN_DIST_PNG,
        ARRAY_LEN_DIST_MIN40_PNG,
        PIPELINE_STAGES_SENTINEL,
        TELOTRON_SPLICE_LOGO_SENTINEL,
        CTRL_SPLICE_LOGO_SENTINEL,
        CTRL_BOUNDARY_KMER_SENTINEL,
        TELOTRON_BOUNDARY_KMER_SENTINEL,
        COMPOSITE_KMER_SENTINEL,
        COMPOSITE_LOGO_SENTINEL,
        COMPOSITE_5P3P_KMER_SENTINEL,
        COMPOSITE_5P3P_LOGO_SENTINEL,


rule motif_discovery:
    input:
        STREME_TELO_SENTINEL,
        STREME_NONTELO_SENTINEL,
        STREME_LINKER_SENTINEL,
        STREME_BRANCHPOINT_SENTINEL,
        FIMO_BRANCHPOINT_SENTINEL,


rule telomerase_search:
    input:
        TELOMERASE_DB_SENTINEL,
        TELOMERASE_BLAST_SENTINEL,
        TELOMERASE_BLAST_REPORT_SENTINEL,
        TERT_DEEP_HOMOLOGY_TSV,


# Core survey + every exploratory branch (the pre-refactor default target).
rule everything:
    input:
        rules.all.input,
        rules.sequences.input,
        rules.figures_extra.input,
        rules.motif_discovery.input,
        rules.telomerase_search.input,


# Build the unified genome manifest from RefSeq assembly summary + Tara SMAGs index.
# Keeps only annotated euk lineages (group, col 25) with a non-"na" FTP path (col 20).
# If config["accessions"] is non-empty, subset to just those genome_ids.
rule manifests:
    output:
        refseq="manifests/refseq_euk.tsv",
        tara="manifests/tara_mags.tsv",
        all="manifests/all_genomes.tsv",
    params:
        accessions=ACCESSIONS,
    shell:
        r"""
        mkdir -p manifests
        curl -L --fail -s {REFSEQ_URL} \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1 || /^#/ {{next}} \
              $25 ~ /^({REFSEQ_GROUPS})$/ && $20!="na" \
              {{print $1,$8,$25,$20}}' \
          > {output.refseq}

        curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
          | tar -tzf - \
          | sed -n 's#.*/\(TARA_.*MAG_[0-9][0-9]*\)\.gmove\.gff#\1#p' \
          | sort -u \
          | awk 'BEGIN{{OFS="\t"}} {{print $1,"Tara Oceans MAG","tara","genoscope"}}' \
          > {output.tara}

        printf "genome_id\torganism\tgroup\tsource\n" > {output.all}
        cat {output.refseq} {output.tara} >> {output.all}

        # Optional accession whitelist: keep header + matching genome_ids.

        # (subset is applied below; canonical_motifs rule consumes the final {output.all})
        if [ -n "{params.accessions}" ]; then
            printf '%s\n' {params.accessions} > manifests/.accessions.txt
            awk 'NR==FNR{{keep[$1];next}} FNR==1 || $1 in keep' \
                manifests/.accessions.txt {output.all} > {output.all}.subset
            mv {output.all}.subset {output.all}
            awk 'NR==FNR{{keep[$1];next}} $1 in keep' \
                manifests/.accessions.txt {output.refseq} > {output.refseq}.subset
            mv {output.refseq}.subset {output.refseq}
            awk 'NR==FNR{{keep[$1];next}} $1 in keep' \
                manifests/.accessions.txt {output.tara} > {output.tara}.subset
            mv {output.tara}.subset {output.tara}
        fi
        """


# Tara SMAGs ship as two monolithic ~50 GB tarballs (contigs + GFF). Stream-extract only
# the MAGs listed in manifests/tara_mags.tsv so the test set doesn't blow up disk.
rule tara_archives:
    input:
        "manifests/tara_mags.tsv",
    output:
        fna="raw/tara/.fna.done",
        gff="raw/tara/.gff.done",
    shell:
        r"""
        mkdir -p raw/tara
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
            | tar --wildcards -xzf - -C raw/tara "${{FNA_PATS[@]}}" ) &
        ( curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
            | tar --wildcards -xzf - -C raw/tara "${{GFF_PATS[@]}}" ) &
        wait
        touch {output.fna} {output.gff}
        """


# Derive per-genome FNA and GFF URLs from the RefSeq FTP path.
rule refseq_urls:
    input:
        "manifests/refseq_euk.tsv",
    output:
        "manifests/refseq_urls.tsv",
    shell:
        r"""
        awk -F'\t' 'BEGIN{{OFS="\t"}} {{
          sub(/\/$/, "", $4);
          n=split($4,a,"/"); base=a[n];
          print $1, $4"/"base"_genomic.fna.gz", $4"/"base"_genomic.gff.gz"
        }}' {input} > {output}
        """


# Fan out per-genome curls. NCBI rate-limits aggressive parallel pulls, so cap the
# inner pool at 8 regardless of THREADS. curl retries handle transient 503s.
# `exit 255` aborts xargs on a persistent failure after retries.
rule download_refseq:
    input:
        "manifests/refseq_urls.tsv",
    output:
        touch("raw/refseq/.done"),
    threads: 8
    shell:
        r"""
        mkdir -p raw/refseq
        awk -F'\t' '{{print $1"\t"$2"\n"$1"\t"$3}}' {input} \
          | xargs -P 8 -n2 sh -c \
              'mkdir -p raw/refseq/$0 && curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s "$1" -o raw/refseq/$0/$(basename "$1") || exit 255'
        touch {output}
        """


# Build a `genome_id → motif` TSV from the curated literature mapping in config.yaml.
# Scan stage uses this in preference to contig-end scanning.
rule canonical_motifs:
    input:
        manifest="manifests/all_genomes.tsv",
    output:
        "manifests/canonical_motifs.tsv",
    params:
        by_genome=CANONICAL_BY_GENOME,
        by_group=CANONICAL_BY_GROUP,
    run:
        import csv
        with open(input.manifest) as f, open(output[0], "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["genome_id", "motif"])
            for r in csv.DictReader(f, delimiter="\t"):
                gid = r["genome_id"]
                grp = r["group"]
                motif = params.by_genome.get(gid)
                if motif is None:
                    motif = params.by_group.get(grp, "")
                w.writerow([gid, motif or ""])


# The core scan: per-intron motif coverage, splice signals, orientation, distance-to-end.
# Broad sieve: keep --min-repeat-frac low so DIVERGING+LINKER architectures (where a
# non-telomeric linker dilutes the per-intron telomeric fraction) survive into the
# loci TSV. Strict cutoffs live downstream in filter_final.
# --max-flank-repeat-frac 0.50 still rejects misannotated telomere arrays (predictors
# sometimes split a long telomere into fake exon/intron structure).
rule scan_all:
    input:
        manifest="manifests/all_genomes.tsv",
        canonical="manifests/canonical_motifs.tsv",
        tara=["raw/tara/.fna.done", "raw/tara/.gff.done"],
        refseq="raw/refseq/.done",
    output:
        loci="results/all_telotron_loci.tsv",
        introns="results/all_introns_scanned.tsv",
        summary="results/all_species_raw_summary.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p results
        python scripts/scan_telotrons.py \
            --manifest {input.manifest} \
            --canonical-motifs {input.canonical} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --min-repeat-frac {SCAN_MIN_FRAC} --max-flank-repeat-frac {SCAN_MAX_FLANK_FRAC} \
            --threads {threads} \
            --loci {output.loci} --introns {output.introns} --summary {output.summary}
        """


# Tighten candidates → final set. --require-terminal-motif-match enforces that
# the intronic motif matches the genome's actual telomere motif (kills cross-motif noise).
rule filter_final:
    input:
        loci="results/all_telotron_loci.tsv",
        introns="results/all_introns_scanned.tsv",
        summary="results/all_species_raw_summary.tsv",
    output:
        final="results/final_telotron_set.tsv",
        species="results/final_species_summary.tsv",
        neg="results/final_negative_controls.tsv",
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
        final="results/final_telotron_set.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        final="results/final_telotron_set_dedup.tsv",
        log="results/dedup_log.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/dedup_telotrons.py \
            --final {input.final} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --out-final {output.final} --out-log {output.log} \
            --threads {threads}
        """


# Classify each (deduplicated) telotron locus into one of:
# GT-F-AG, GT-R-AG, GT-F-R-AG, GT-F-linker-R-AG, GT-R-linker-F-AG, Other.
# Emits the long-form boundary-kmer-by-architecture TSV used by the per-species
# figures (donor / acceptor / linker_left / linker_right 6-mers).
rule classify_architecture:
    input:
        final="results/final_telotron_set_dedup.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        loci="results/final_telotron_set_architecture.tsv",
        kmers="results/boundary_kmers_by_architecture.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/classify_telotron_architecture.py \
            --final {input.final} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --out-loci {output.loci} --out-kmers {output.kmers}
        """


# Boundary k-mer enrichment, length-matched distance-to-end test, architecture summary.
rule analyze:
    input:
        final="results/final_telotron_set_dedup.tsv",
        introns="results/all_introns_scanned.tsv",
        species="results/final_species_summary.tsv",
        neg="results/final_negative_controls.tsv",
    output:
        kmers="results/boundary_kmer_enrichment.tsv",
        dist="results/distance_to_end.tsv",
        arch="results/architecture_summary.tsv",
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


# Four headline figures (counts, orientation, boundary k-mers, distance scatter).
rule figures:
    input:
        species="results/final_species_summary.tsv",
        final="results/final_telotron_set.tsv",
        kmers="results/boundary_kmer_enrichment.tsv",
        dist="results/distance_to_end.tsv",
        arch="results/architecture_summary.tsv",
    output:
        FIGURES,
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures
        python scripts/plot_telotrons.py \
            --species {input.species} --final {input.final} \
            --kmers {input.kmers} --distance {input.dist} --architecture {input.arch} \
            --outdir results/figures
        """


# Per-species "where do telomere repeats sit?" plots — one figure per species,
# one panel per top-N contig, x = genomic coordinate, y = rolling hit density of the
# species's identified terminal motif. Telotron candidate loci overlaid as ticks.
rule terminal_motif_figures:
    input:
        summary="results/all_species_raw_summary.tsv",
        loci="results/all_telotron_loci.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        touch(TERMINAL_DENSITY_SENTINEL),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures/terminal_motif
        python scripts/plot_terminal_motif_density.py \
            --summary {input.summary} --loci {input.loci} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --outdir results/figures/terminal_motif \
            --threads {threads}
        """


# Per-species, per-architecture FASTA + flanked-text extracts.
# Flanked lines: [LEFT100] [INTRON] [RIGHT100], or for linker archs
# [LEFT100] [ARRAY1] [LINKER] [ARRAY2] [RIGHT100]. One subdir per species,
# one file per architecture.
rule extract_telotron_fasta:
    input:
        arch="results/final_telotron_set_architecture.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        touch(EXTRACT_FASTA_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/telotron_fasta results/telotron_flanked
        python scripts/extract_telotron_fasta.py \
            --final {input.arch} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --fasta-dir results/telotron_fasta \
            --flanked-dir results/telotron_flanked
        touch {output}
        """


# BLAST linker sequences from linker-architecture telotrons (1) against their own
# species genome and (2) against a concatenated all-species DB. Outputs hit TSVs
# annotated with the query genome_id so cross-species recurrences are findable.
rule blast_linkers:
    input:
        arch="results/final_telotron_set_architecture.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        own="results/linker_blast_hits_own_genome.tsv",
        all="results/linker_blast_hits_all_genomes.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/blast_linkers
        python scripts/blast_linkers.py \
            --arch-tsv {input.arch} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --workdir results/blast_linkers \
            --out-own {output.own} --out-all {output.all} \
            --threads {threads}
        """


# Per-species × per-architecture MAFFT MSAs of the locus regions
# (upstream_50 | intron|arm1[+linker+arm2]|arm1[+arm2] | downstream_50).
# Both raw and homopolymer-compressed alignments; combined.aln.txt collapses
# everything into one fixed-width view with spaces at region boundaries.
rule msa_telotron_regions:
    input:
        arch="results/final_telotron_set_architecture.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        touch(MSA_SENTINEL),
    threads: 4
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/msa_regions
        python scripts/msa_telotron_regions.py \
            --arch-tsv {input.arch} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir results/msa_regions --threads {threads}
        touch {output}
        """


# UniProt SwissProt eukaryote-rich proteome (~90 MB) — used by miniprot to
# mask homology-detectable CDS in genomes whose published GFF misses things
# (including pseudogenes, which miniprot keeps as Frameshift=/StopCodon= hits).
# To override, place your own proteome at manifests/proteome.faa.gz before
# running and Snakemake will skip the download.
rule download_proteome:
    output:
        "manifests/proteome.faa.gz",
    shell:
        r"""
        mkdir -p manifests
        curl -sSL --retry 5 --retry-delay 3 -o {output} \
            https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        """


# Per-genome closest-relative proteomes (same-genus RefSeq siblings, up to
# 5 each). Layered on top of the SwissProt base for miniprot homology mapping
# so divergent lineage-specific genes still register. TARA MAGs and genomes
# with no same-genus annotated relative get empty placeholders (→ SwissProt
# alone). Per-sibling cache lives under manifests/close_proteomes/_cache/.
rule fetch_close_proteomes:
    input:
        manifest="results/all_species_raw_summary.tsv",
    output:
        touch("manifests/close_proteomes/.done"),
    params:
        outdir="manifests/close_proteomes",
    shell:
        r"""
        python scripts/fetch_close_proteomes.py \
            --manifest {input.manifest} \
            --outdir {params.outdir}
        """


# Per-genome ORF mask: six-frame ATG→stop scan, ORFs >= min_orf_nt.
# Output is a sorted+merged BED per genome under results/masks/. Used by
# find_interstitial_arrays to subtract on top of the annotated gene/intron mask.
# Cutoff rationale: on a shuffled E. necatrix genome ORFs >= 300 nt cover 10.6%
# of bp (noise), ORFs >= 450 nt cover 1.9% (noise floor), so 450 is the default.
rule make_unannotated_masks:
    input:
        manifest="results/all_species_raw_summary.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        touch("results/masks/.done"),
    params:
        outdir="results/masks",
        min_orf=450,
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p {params.outdir}
        python scripts/make_unannotated_mask.py \
            --manifest {input.manifest} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir {params.outdir} \
            --min-orf-nt {params.min_orf}
        """


# All interstitial telomeric repeat arrays — i.e. arrays that are
# non-terminal (>=5 kb from contig end), non-genic, and non-intronic.
# Exclusion set = (annotated genes ∪ gt-derived introns) ∪ ORF mask
# (results/masks/{gid}.bed: six-frame ORFs >= 450 nt).
# seqkit locate → bedtools merge → bedtools intersect -v vs that exclusion BED.
# Emits 5'/3' boundary 6- and 12-mers.
rule find_interstitial_arrays:
    input:
        manifest="results/all_species_raw_summary.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
        masks="results/masks/.done",
    output:
        "results/interstitial_arrays.tsv",
    params:
        mask_dir="results/masks",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/find_interstitial_arrays.py \
            --manifest {input.manifest} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --mask-dir {params.mask_dir} \
            --out {output} --threads {threads}
        """


# Per-species figure: top boundary 6-mers (GT-XXXX donor, XXXX-AG acceptor,
# and 4+2 linker-edge 6-mers) split by architecture.
rule plot_boundary_kmers_by_arch:
    input:
        kmers="results/boundary_kmers_by_architecture.tsv",
        species="results/final_species_summary.tsv",
    output:
        touch(ARCH_KMER_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures/boundary_kmers_by_arch
        python scripts/plot_boundary_kmers_by_arch.py \
            --kmers {input.kmers} --species {input.species} \
            --outdir results/figures/boundary_kmers_by_arch
        touch {output}
        """


# Putative-telotron splice candidates: C-rich-on-plus interstitial arrays
# bracketed by a GT (donor) immediately upstream or AG (acceptor) immediately
# downstream on the plus strand — consistent with an unannotated antisense
# (mRNA-CCCTAA) telotron inside a plus-strand gene.
rule filter_interstitial_splice_candidates:
    input:
        arrays="results/interstitial_arrays.tsv",
    output:
        tsv="results/interstitial_arrays_c_rich_splice_candidates.tsv",
    conda: ENV
    shell:
        r"""
        python scripts/filter_interstitial_splice_candidates.py \
            --arrays {input.arrays} --out {output.tsv}
        """


rule plot_interstitial_splice_candidate_kmers:
    input:
        arrays="results/interstitial_arrays_c_rich_splice_candidates.tsv",
    output:
        touch(INTERSTITIAL_SPLICE_KMER_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/interstitial_boundary_kmers_c_rich_splice_candidates
        python scripts/plot_interstitial_boundary_kmers.py \
            --arrays {input.arrays} \
            --outdir results/figures/interstitial_boundary_kmers_c_rich_splice_candidates
        """


# Per-species figure: top 5' and 3' boundary 6-mers and 12-mers of interstitial
# (non-genic, non-intronic, non-terminal) telomeric repeat arrays.
rule plot_interstitial_boundary_kmers:
    input:
        arrays="results/interstitial_arrays.tsv",
    output:
        touch(INTERSTITIAL_KMER_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures/interstitial_boundary_kmers
        python scripts/plot_interstitial_boundary_kmers.py \
            --arrays {input.arrays} \
            --outdir results/figures/interstitial_boundary_kmers
        touch {output}
        """


# Per-species sequence logos (information bits, logomaker) of the 5'-flank
# + first repeat unit and last repeat unit + 3'-flank for interstitial arrays.
rule plot_interstitial_boundary_logos:
    input:
        arrays="results/interstitial_arrays.tsv",
    output:
        touch(INTERSTITIAL_LOGO_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures/interstitial_boundary_logos
        python scripts/plot_interstitial_boundary_logos.py \
            --arrays {input.arrays} \
            --outdir results/figures/interstitial_boundary_logos \
            --flank-len 10
        touch {output}
        """


rule plot_array_length_distribution:
    """Per-species histograms: ITS array_len | telotron telomeric_bases |
    non-telotron intron_len. Eimeria + MAGs only."""
    input:
        interstitial="results/interstitial_arrays.tsv",
        telotrons="results/final_telotron_set.tsv",
        non_telotrons="results/non_telotron_controls.tsv",
    output:
        ARRAY_LEN_DIST_PNG,
    conda:
        ENV
    shell:
        r"""
        python scripts/plot_array_length_distribution.py \
            --interstitial {input.interstitial} \
            --telotrons {input.telotrons} \
            --non-telotrons {input.non_telotrons} \
            --out {output}
        """


rule plot_array_length_distribution_min40:
    """Same but drops arrays/telotrons/introns < 40 bp."""
    input:
        interstitial="results/interstitial_arrays.tsv",
        telotrons="results/final_telotron_set.tsv",
        non_telotrons="results/non_telotron_controls.tsv",
    output:
        ARRAY_LEN_DIST_MIN40_PNG,
    conda:
        ENV
    shell:
        r"""
        python scripts/plot_array_length_distribution.py \
            --interstitial {input.interstitial} \
            --telotrons {input.telotrons} \
            --non-telotrons {input.non_telotrons} \
            --min-len 40 \
            --out {output}
        """


# Composite figure (PDF + tall PNG): per-species rows, intron architecture
# boundary panel on the left, interstitial array boundary panel on the right.
rule plot_combined_boundary_kmers:
    input:
        species="results/final_species_summary.tsv",
        intron_dir=ARCH_KMER_SENTINEL,
        interst_dir=INTERSTITIAL_KMER_SENTINEL,
    output:
        pdf="results/figures/boundary_kmers_combined.pdf",
        png="results/figures/boundary_kmers_combined.png",
    conda:
        ENV
    shell:
        r"""
        python scripts/plot_combined_boundary_kmers.py \
            --species {input.species} \
            --intron-dir results/figures/boundary_kmers_by_arch \
            --interstitial-dir results/figures/interstitial_boundary_kmers \
            --out-pdf {output.pdf} --out-png {output.png}
        """


# Pipeline-stage summary figures (manifests, scan, filter funnel, distance,
# architecture).
rule plot_pipeline_stages:
    input:
        manifest="manifests/all_genomes.tsv",
        raw="results/all_species_raw_summary.tsv",
        sp="results/final_species_summary.tsv",
        dist="results/distance_to_end.tsv",
        arch="results/final_telotron_set_architecture.tsv",
    output:
        touch(PIPELINE_STAGES_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p results/figures/pipeline_stages
        python scripts/plot_pipeline_stages.py \
            --manifests {input.manifest} \
            --raw-summary {input.raw} \
            --species-final {input.sp} \
            --distance {input.dist} \
            --architecture-loci {input.arch} \
            --outdir results/figures/pipeline_stages
        touch {output}
        """


# ── STREME (MEME suite) de novo motif discovery ──────────────────────────────
# Run STREME separately on three positive-sequence sets:
#   1. telotron intron sequences (full intron, all architectures)
#   2. non-telomeric introns       (intron_telomeric_frac < 0.10, random sample)
#   3. linker sequences            (middle non-telomeric block from linker-arch telotrons)
# STREME generates its own shuffled control. Inputs are built once by build_streme_inputs.py.

rule build_streme_inputs:
    input:
        telotrons_done=EXTRACT_FASTA_SENTINEL,
        introns="results/all_introns_scanned.tsv",
        linkers="results/blast_linkers/linker_queries/_all_linkers.fa",
    output:
        telo="results/streme_inputs/telotrons.fa",
        non_telo="results/streme_inputs/non_telo_introns.fa",
        linkers="results/streme_inputs/linkers.fa",
    params:
        n_non_telo=5000,
        non_telo_max_frac=0.10,
    conda: ENV
    shell:
        r"""
        python scripts/build_streme_inputs.py \
            --telotron-fasta-dir results/telotron_fasta \
            --introns-tsv {input.introns} \
            --linkers-fa {input.linkers} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir results/streme_inputs \
            --n-non-telo {params.n_non_telo} \
            --non-telo-max-frac {params.non_telo_max_frac}
        """


# STREME parameters: DNA mode, motif width 4-12, 10 motifs per run, p-value 0.05.
# STREME shuffles the input to build its own background.
STREME_ARGS = "--dna --minw 4 --maxw 12 --nmotifs 10 --thresh 0.05"

rule streme_telotrons:
    input: "results/streme_inputs/telotrons.fa"
    output: STREME_TELO_SENTINEL
    conda: MEME_ENV
    shell:
        r"""
        streme {STREME_ARGS} --p {input} --oc results/streme/telotrons
        """

rule streme_non_telo_introns:
    input: "results/streme_inputs/non_telo_introns.fa"
    output: STREME_NONTELO_SENTINEL
    conda: MEME_ENV
    shell:
        r"""
        streme {STREME_ARGS} --p {input} --oc results/streme/non_telo_introns
        """

rule streme_linkers:
    input: "results/streme_inputs/linkers.fa"
    output: STREME_LINKER_SENTINEL
    conda: MEME_ENV
    shell:
        r"""
        streme {STREME_ARGS} --p {input} --oc results/streme/linkers
        """


# ── Branchpoint motif discovery ─────────────────────────────────────────────
# Restrict to the last 80 bp of each canonical GT-AG intron (where the branch-
# point sits, -18 to -40 nt from the AG). Discriminative STREME with a length-
# matched mid-intron window as the explicit negative removes the AT-richness/
# simple-repeat noise that drowned out branchpoint signal in the full-intron
# run. Two clade pools (yeast+apicomplexa "tight" expecting YACTAAC; everything
# else "loose" expecting YNYTRAY) keep distinct consensuses from blurring.
def _branchpoint_genome_ids(min_introns=50):
    """Genomes with >= min_introns qualifying GT-AG introns in all_introns_scanned.tsv."""
    path = "results/all_introns_scanned.tsv"
    if not os.path.exists(path):
        return []
    counts = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        try:
            gi = h.index("genome_id")
            fi = h.index("telomeric_frac")
            li = h.index("intron_len")
            si = h.index("splice_class")
        except ValueError:
            return []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                if (f[si] == "GT-AG"
                        and float(f[fi]) < 0.10
                        and int(f[li]) >= 400):
                    counts[f[gi]] = counts.get(f[gi], 0) + 1
            except (ValueError, IndexError):
                continue
    return [gid for gid, n in counts.items() if n >= min_introns]


rule build_branchpoint_inputs_one:
    """Extract 3'-end branchpoint windows (pos) and mid-intron controls (neg)
    for a single genome, then dustmask."""
    input:
        introns="results/all_introns_scanned.tsv",
        tara="raw/tara/.fna.done",
        refseq="raw/refseq/.done",
    output:
        pos="results/streme_inputs/branchpoint_{gid}_pos.fa",
        neg="results/streme_inputs/branchpoint_{gid}_neg.fa",
    conda: ENV
    shell:
        r"""
        python scripts/build_branchpoint_inputs.py \
            --introns {input.introns} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir results/streme_inputs \
            --genome-id {wildcards.gid} \
            --pos-window 200 --min-intron-len 400 \
            --n-per-species 3000
        """


# Branchpoint-targeted STREME per species: 6-9 bp window, discriminative.
STREME_BP_ARGS = "--dna --minw 6 --maxw 9 --nmotifs 5 --thresh 0.05 --objfun de"

rule streme_branchpoint_one:
    input:
        pos="results/streme_inputs/branchpoint_{gid}_pos.fa",
        neg="results/streme_inputs/branchpoint_{gid}_neg.fa",
    output:
        "results/streme/branchpoint/{gid}/streme.xml",
    conda: MEME_ENV
    shell:
        r"""
        streme {STREME_BP_ARGS} --p {input.pos} --n {input.neg} \
            --oc results/streme/branchpoint/{wildcards.gid}
        """


rule streme_branchpoint:
    input:
        lambda wc: expand(
            "results/streme/branchpoint/{gid}/streme.xml",
            gid=_branchpoint_genome_ids(),
        ),
    output:
        touch(STREME_BRANCHPOINT_SENTINEL),


rule fimo_branchpoint_one:
    input:
        pos="results/streme_inputs/branchpoint_{gid}_pos.fa",
        pwms="scripts/branchpoint_pwms.meme",
    output:
        "results/fimo/branchpoint/{gid}/fimo.tsv",
    conda: MEME_ENV
    shell:
        r"""
        fimo --thresh 1e-3 --oc results/fimo/branchpoint/{wildcards.gid} \
            {input.pwms} {input.pos}
        """


rule fimo_branchpoint:
    input:
        lambda wc: expand(
            "results/fimo/branchpoint/{gid}/fimo.tsv",
            gid=_branchpoint_genome_ids(),
        ),
    output:
        touch(FIMO_BRANCHPOINT_SENTINEL),


# ── Telomerase Database (telomerase.us) reference BLAST dbs ─────────────────
# Scrape category pages for GenBank accessions, fetch FASTAs via NCBI E-utils,
# and build per-category BLAST databases. Nucleotide-only for TR (RNA
# component); both nucleotide and protein for TERT and Other (accessory
# proteins). HTML pages are cached under raw/telomerase.us/ so re-runs don't
# re-hit the site. Set NCBI_API_KEY in your shell to bump the E-utils rate.
rule fetch_telomerase_db:
    output:
        tr_nucl="results/telomerase_db/tr.nucl.fa",
        tert_nucl="results/telomerase_db/tert.nucl.fa",
        tert_prot="results/telomerase_db/tert.prot.fa",
        other_nucl="results/telomerase_db/other.nucl.fa",
        other_prot="results/telomerase_db/other.prot.fa",
    shell:
        r"""
        mkdir -p raw/telomerase.us results/telomerase_db
        python scripts/fetch_telomerase_db.py \
            --cache-dir raw/telomerase.us \
            --outdir results/telomerase_db
        """


# BLAST the curated telomerase reference set (telomerase.us) AS QUERY against
# each of the surveyed genomes. Permissive params (low word size, -evalue 10,
# masking disabled) so short divergent hits survive. Each HSP is written as
# the standard outfmt 6 12-column tabular row followed by an ASCII pairwise
# alignment (qseq / match line / sseq).
#
# Output: results/blast_telomerase_vs_genomes/{genome_id}/{tr_nucl,tert_nucl,
# tert_prot,other_nucl,other_prot}.tsv
# Telomerase-ortholog blast targets (config: blast_genome_ids): the 5 Eimeria
# spp. + 3 telotron-positive Tara MAGs. An explicit list rather than a prefix
# filter on the survey summary (which wrongly caught GCF_000499845.2 /
# Phaseolus vulgaris and may omit the Tara MAGs).
_BLAST_GENOME_IDS = config.get("blast_genome_ids") or [
    "GCF_000499385.1",   # Eimeria necatrix
    "GCF_000499425.1",   # Eimeria acervulina
    "GCF_000499545.2",   # Eimeria tenella
    "GCF_000499605.1",   # Eimeria maxima
    "GCF_000499745.2",   # Eimeria mitis
    "TARA_AON_82_MAG_00313",
    "TARA_MED_95_MAG_00464",
    "TARA_PSW_86_MAG_00284",
]


# Per-genome BLAST: each genome's BLASTs run as an independent Snakemake job,
# so 22 genomes can fan out across the cluster. Each job builds a temporary
# nucleotide BLAST db for its genome and runs all 5 query sets against it.
rule blast_telomerase_vs_one_genome:
    input:
        tr_nucl="results/telomerase_db/tr.nucl.fa",
        tert_nucl="results/telomerase_db/tert.nucl.fa",
        tert_prot="results/telomerase_db/tert.prot.fa",
        other_nucl="results/telomerase_db/other.nucl.fa",
        other_prot="results/telomerase_db/other.prot.fa",
        refseq="raw/refseq/.done",
        tara="raw/tara/.fna.done",
    output:
        touch("results/blast_telomerase_vs_genomes/{gid}/.done"),
    threads: 4
    shell:
        r"""
        mkdir -p results/blast_telomerase_vs_genomes/{wildcards.gid}
        python scripts/blast_telomerase_vs_genomes.py \
            --genome-id {wildcards.gid} \
            --db-dir results/telomerase_db \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir results/blast_telomerase_vs_genomes/{wildcards.gid} \
            --threads {threads}
        """


# Aggregator: depends on a per-genome .done sentinel for every surveyed genome.
rule blast_telomerase_vs_genomes:
    input:
        lambda wc: expand(
            "results/blast_telomerase_vs_genomes/{gid}/.done",
            gid=_BLAST_GENOME_IDS,
        ),
    output:
        touch(TELOMERASE_BLAST_SENTINEL),
    shell:
        "touch {output}"


# Per-genome post-processor: dedup HSPs at the same locus, split into per-query
# sections, and re-align the full query to a ±qlen reference window using
# Biopython's semi-global aligner (free reference end-gaps). Output is one
# pretty-printed report per query set plus a deduplicated summary TSV.
rule process_telomerase_blast_one_genome:
    input:
        blast_done="results/blast_telomerase_vs_genomes/{gid}/.done",
        tr_nucl="results/telomerase_db/tr.nucl.fa",
        tert_nucl="results/telomerase_db/tert.nucl.fa",
        tert_prot="results/telomerase_db/tert.prot.fa",
        other_nucl="results/telomerase_db/other.nucl.fa",
        other_prot="results/telomerase_db/other.prot.fa",
        refseq="raw/refseq/.done",
        tara="raw/tara/.fna.done",
    output:
        touch("results/blast_telomerase_vs_genomes/{gid}/.report.done"),
    conda: ENV
    shell:
        r"""
        python scripts/process_telomerase_blast.py \
            --genome-id {wildcards.gid} \
            --blast-dir results/blast_telomerase_vs_genomes/{wildcards.gid} \
            --db-dir results/telomerase_db \
            --refseq-dir raw/refseq --tara-dir raw/tara
        """


rule process_telomerase_blast:
    input:
        lambda wc: expand(
            "results/blast_telomerase_vs_genomes/{gid}/.report.done",
            gid=_BLAST_GENOME_IDS,
        ),
    output:
        touch(TELOMERASE_BLAST_REPORT_SENTINEL),
    shell:
        "touch {output}"



# ── TERT by deep homology (miniprot + Pfam domain architecture) ─────────────
# BLAST/tblastn fail on apicomplexan TERT (~25% id, multi-exon, unannotated
# across Coccidia). This rule recovers it the field-standard way: intron-aware
# miniprot with phylogenetically-local apicomplexan TERT seeds, then confirm by
# domain architecture -- a real TERT carries the telomerase RBD (Pfam PF12009)
# upstream of the RT domain (PF00078); abundant Eimeria retroelement RTs carry
# the RT but never the RBD. Iterative rounds add confirmed proteins back as
# within-clade seeds to rescue fragmented/divergent genomes. Toxoplasma is
# included as a positive control (sister taxon; its TERT is also unannotated).
TERT_GENOME_IDS = config.get("tert_genome_ids") or [
    "GCF_000499385.1",   # Eimeria necatrix
    "GCF_000499425.1",   # Eimeria acervulina
    "GCF_000499545.2",   # Eimeria tenella
    "GCF_000499605.1",   # Eimeria maxima
    "GCF_000499745.2",   # Eimeria mitis
    "GCF_000006565.2",   # Toxoplasma gondii (sister-taxon positive control)
]


rule fetch_tert_seeds_hmms:
    output:
        seeds="results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="results/tert_deep_homology/refs/PF12009.hmm",
        rt="results/tert_deep_homology/refs/PF00078.hmm",
    shell:
        r"""
        python scripts/fetch_tert_seeds_hmms.py \
            --outdir results/tert_deep_homology/refs
        """


rule find_tert:
    input:
        seeds="results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="results/tert_deep_homology/refs/PF12009.hmm",
        rt="results/tert_deep_homology/refs/PF00078.hmm",
        refseq="raw/refseq/.done",
        tara="raw/tara/.fna.done",
    output:
        tsv=TERT_DEEP_HOMOLOGY_TSV,
        faa="results/tert_deep_homology/all_confirmed.faa",
    params:
        gids=",".join(TERT_GENOME_IDS),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/find_tert_deep_homology.py \
            --genome-ids {params.gids} \
            --seeds {input.seeds} --trbd-hmm {input.trbd} --rt-hmm {input.rt} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --outdir results/tert_deep_homology \
            --iterate 3 --threads {threads}
        """


rule telomerase_db_blastdbs:
    input:
        tr_nucl="results/telomerase_db/tr.nucl.fa",
        tert_nucl="results/telomerase_db/tert.nucl.fa",
        tert_prot="results/telomerase_db/tert.prot.fa",
        other_nucl="results/telomerase_db/other.nucl.fa",
        other_prot="results/telomerase_db/other.prot.fa",
    output:
        touch(TELOMERASE_DB_SENTINEL),
    shell:
        r"""
        cd results/telomerase_db
        for f in tr.nucl.fa tert.nucl.fa other.nucl.fa; do
            [ -s "$f" ] && makeblastdb -in "$f" -dbtype nucl -parse_seqids -out "${{f%.fa}}" || echo "skip empty $f"
        done
        for f in tert.prot.fa other.prot.fa; do
            [ -s "$f" ] && makeblastdb -in "$f" -dbtype prot -parse_seqids -out "${{f%.fa}}" || echo "skip empty $f"
        done
        """


# ── Splice-signal sequence logos and non-telotron control set ───────────────
# Telotron splice-signal logos: ±10 bp of the GT (donor) and AG (acceptor).
# Reads the flanked .txt files produced by extract_telotron_fasta.
rule plot_telotron_splice_logos:
    input:
        EXTRACT_FASTA_SENTINEL,
    output:
        touch(TELOTRON_SPLICE_LOGO_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/telotron_splice_logos
        python scripts/plot_splice_signal_logos.py \
            --flanked-dir results/telotron_flanked \
            --outdir results/figures/telotron_splice_logos
        touch {output}
        """


# Definitely-non-telotron intron control set: introns with telomeric_frac < 10%
# from the same positive species, sampled per genome. Produced TSV mirrors the
# architecture-table schema so extract_telotron_fasta.py can ingest it.
rule build_non_telotron_controls:
    input:
        introns="results/all_introns_scanned.tsv",
        final="results/final_telotron_set.tsv",
    output:
        tsv="results/non_telotron_controls.tsv",
    conda: ENV
    shell:
        r"""
        python scripts/build_non_telotron_controls.py \
            --introns {input.introns} --final {input.final} \
            --out {output.tsv} --max-frac 0.10 --n-per-species 5000
        """


# Per-species and per-architecture FASTAs + flanked .txt for the control set,
# under results/non_telotron_fasta and results/non_telotron_flanked (architecture
# is always "control"). Reuses the telotron extractor.
rule extract_non_telotron_fasta:
    input:
        controls="results/non_telotron_controls.tsv",
        tara=["raw/tara/.fna.done"],
        refseq="raw/refseq/.done",
    output:
        touch(CTRL_FLANKED_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/non_telotron_fasta results/non_telotron_flanked
        python scripts/extract_telotron_fasta.py \
            --final {input.controls} \
            --refseq-dir raw/refseq --tara-dir raw/tara \
            --fasta-dir results/non_telotron_fasta \
            --flanked-dir results/non_telotron_flanked
        """


# Splice-signal logos for the control set (same script as telotron version).
rule plot_non_telotron_splice_logos:
    input:
        CTRL_FLANKED_SENTINEL,
    output:
        touch(CTRL_SPLICE_LOGO_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/non_telotron_splice_logos
        python scripts/plot_splice_signal_logos.py \
            --flanked-dir results/non_telotron_flanked \
            --outdir results/figures/non_telotron_splice_logos
        touch {output}
        """


# Boundary 6-mer / 12-mer rank plots for the control set, same layout as the
# interstitial boundary k-mer figures. Uses first40/last40 already in the TSV.
rule plot_non_telotron_boundary_kmers:
    input:
        controls="results/non_telotron_controls.tsv",
    output:
        touch(CTRL_BOUNDARY_KMER_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/non_telotron_boundary_kmers
        python scripts/plot_intron_boundary_kmers.py \
            --introns {input.controls} \
            --outdir results/figures/non_telotron_boundary_kmers \
            --label "non-telotron control (telomeric_frac<0.10)"
        """


# Per-species 2x2 boundary k-mer figures for the FINAL telotron set, using the
# same plot_intron_boundary_kmers script that produces the control figures.
# Filenames slugged as "{gid}_{organism}.png" so they line up with the splice
# and interstitial figure dirs for compositing.
rule plot_telotron_boundary_kmers:
    input:
        final="results/final_telotron_set.tsv",
    output:
        touch(TELOTRON_BOUNDARY_KMER_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/telotron_boundary_kmers
        python scripts/plot_intron_boundary_kmers.py \
            --introns {input.final} \
            --outdir results/figures/telotron_boundary_kmers \
            --label "telotrons"
        """


# Composite per-species panels: interstitial | non-telotron intron | telotron,
# side by side, for both the boundary k-mer figures and the sequence logos.
rule plot_composite_boundary_kmers:
    input:
        interst=INTERSTITIAL_KMER_SENTINEL,
        non_telo=CTRL_BOUNDARY_KMER_SENTINEL,
        telo=TELOTRON_BOUNDARY_KMER_SENTINEL,
    output:
        touch(COMPOSITE_KMER_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/composite_boundary_kmers
        python scripts/plot_composite_per_species.py \
            --outdir results/figures/composite_boundary_kmers \
            --panel "results/figures/interstitial_boundary_kmers:interstitial arrays" \
            --panel "results/figures/non_telotron_boundary_kmers:non-telotron control introns" \
            --panel "results/figures/telotron_boundary_kmers:telotrons"
        """


rule plot_composite_boundary_logos:
    input:
        interst=INTERSTITIAL_LOGO_SENTINEL,
        non_telo=CTRL_SPLICE_LOGO_SENTINEL,
        telo=TELOTRON_SPLICE_LOGO_SENTINEL,
    output:
        touch(COMPOSITE_LOGO_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/composite_boundary_logos
        python scripts/plot_composite_per_species.py \
            --outdir results/figures/composite_boundary_logos \
            --panel "results/figures/interstitial_boundary_logos:interstitial array boundary" \
            --panel "results/figures/non_telotron_splice_logos:non-telotron intron splice signal" \
            --panel "results/figures/telotron_splice_logos:telotron splice signal"
        """


# ── 5'/3' end-orientation stratified figures ───────────────────────────────
# For each telotron locus we classify the 5' and 3' ends independently (from
# architecture); for interstitial arrays the array orientation defines both
# ends. The final compiled figure shows every category side by side per species.

rule plot_interstitial_boundary_kmers_by_cat:
    input:
        arrays="results/interstitial_arrays.tsv",
    output:
        touch("results/figures/interstitial_boundary_kmers_{cat}/.done"),
    conda: ENV
    wildcard_constraints:
        cat="GG|GC|CG|CC",
    params:
        col=lambda w: INTERSTITIAL_CAT[w.cat][0],
        val=lambda w: INTERSTITIAL_CAT[w.cat][1],
    shell:
        r"""
        mkdir -p results/figures/interstitial_boundary_kmers_{wildcards.cat}
        python scripts/plot_interstitial_boundary_kmers.py \
            --arrays {input.arrays} \
            --outdir results/figures/interstitial_boundary_kmers_{wildcards.cat} \
            --filter-col {params.col} --filter-value "{params.val}"
        """


rule plot_interstitial_boundary_logos_by_cat:
    input:
        arrays="results/interstitial_arrays.tsv",
    output:
        touch("results/figures/interstitial_boundary_logos_{cat}/.done"),
    conda: ENV
    wildcard_constraints:
        cat="GG|GC|CG|CC",
    params:
        col=lambda w: INTERSTITIAL_CAT[w.cat][0],
        val=lambda w: INTERSTITIAL_CAT[w.cat][1],
    shell:
        r"""
        mkdir -p results/figures/interstitial_boundary_logos_{wildcards.cat}
        python scripts/plot_interstitial_boundary_logos.py \
            --arrays {input.arrays} \
            --outdir results/figures/interstitial_boundary_logos_{wildcards.cat} \
            --flank-len 10 \
            --filter-col {params.col} --filter-value "{params.val}"
        """


rule plot_telotron_boundary_kmers_by_cat:
    input:
        final="results/final_telotron_set_architecture.tsv",
    output:
        touch("results/figures/telotron_boundary_kmers_{cat}/.done"),
    conda: ENV
    wildcard_constraints:
        cat="GG|GC|CG|CC",
    params:
        col=lambda w: TELOTRON_CAT[w.cat][0],
        val=lambda w: TELOTRON_CAT[w.cat][1],
    shell:
        r"""
        mkdir -p results/figures/telotron_boundary_kmers_{wildcards.cat}
        python scripts/plot_intron_boundary_kmers.py \
            --introns {input.final} \
            --outdir results/figures/telotron_boundary_kmers_{wildcards.cat} \
            --label "telotrons 5'/3' = {wildcards.cat}" \
            --filter-col {params.col} --filter-value "{params.val}"
        """


rule plot_telotron_splice_logos_by_cat:
    input:
        sentinel=EXTRACT_FASTA_SENTINEL,
        ann="results/final_telotron_set_architecture.tsv",
    output:
        touch("results/figures/telotron_splice_logos_{cat}/.done"),
    conda: ENV
    wildcard_constraints:
        cat="GG|GC|CG|CC",
    params:
        col=lambda w: TELOTRON_CAT[w.cat][0],
        val=lambda w: TELOTRON_CAT[w.cat][1],
    shell:
        r"""
        mkdir -p results/figures/telotron_splice_logos_{wildcards.cat}
        python scripts/plot_splice_signal_logos.py \
            --flanked-dir results/telotron_flanked \
            --outdir results/figures/telotron_splice_logos_{wildcards.cat} \
            --annotation-tsv {input.ann} \
            --filter-col {params.col} --filter-value "{params.val}"
        """


# One compiled side-by-side per-species figure containing every 5'/3' category.
# Layout (left to right):
#   interstitial GG │ interstitial CC │ non-telotron control │
#   telotron GG │ telotron GC │ telotron CG │ telotron CC
rule plot_composite_5p3p_boundary_kmers:
    input:
        i_gg="results/figures/interstitial_boundary_kmers_GG/.done",
        i_gc="results/figures/interstitial_boundary_kmers_GC/.done",
        i_cg="results/figures/interstitial_boundary_kmers_CG/.done",
        i_cc="results/figures/interstitial_boundary_kmers_CC/.done",
        ctrl=CTRL_BOUNDARY_KMER_SENTINEL,
        t_gg="results/figures/telotron_boundary_kmers_GG/.done",
        t_gc="results/figures/telotron_boundary_kmers_GC/.done",
        t_cg="results/figures/telotron_boundary_kmers_CG/.done",
        t_cc="results/figures/telotron_boundary_kmers_CC/.done",
    output:
        touch(COMPOSITE_5P3P_KMER_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/composite_boundary_kmers_5p3p
        python scripts/plot_composite_per_species.py \
            --outdir results/figures/composite_boundary_kmers_5p3p \
            --panel "results/figures/interstitial_boundary_kmers_GG:interstitial 5'=G 3'=G" \
            --panel "results/figures/interstitial_boundary_kmers_GC:interstitial 5'=G 3'=C (hybrid)" \
            --panel "results/figures/interstitial_boundary_kmers_CG:interstitial 5'=C 3'=G (hybrid)" \
            --panel "results/figures/interstitial_boundary_kmers_CC:interstitial 5'=C 3'=C" \
            --panel "results/figures/non_telotron_boundary_kmers:non-telotron control" \
            --panel "results/figures/telotron_boundary_kmers_GG:telotron 5'=G 3'=G" \
            --panel "results/figures/telotron_boundary_kmers_GC:telotron 5'=G 3'=C (hybrid)" \
            --panel "results/figures/telotron_boundary_kmers_CG:telotron 5'=C 3'=G (hybrid)" \
            --panel "results/figures/telotron_boundary_kmers_CC:telotron 5'=C 3'=C"
        """


rule plot_composite_5p3p_boundary_logos:
    input:
        i_gg="results/figures/interstitial_boundary_logos_GG/.done",
        i_gc="results/figures/interstitial_boundary_logos_GC/.done",
        i_cg="results/figures/interstitial_boundary_logos_CG/.done",
        i_cc="results/figures/interstitial_boundary_logos_CC/.done",
        ctrl=CTRL_SPLICE_LOGO_SENTINEL,
        t_gg="results/figures/telotron_splice_logos_GG/.done",
        t_gc="results/figures/telotron_splice_logos_GC/.done",
        t_cg="results/figures/telotron_splice_logos_CG/.done",
        t_cc="results/figures/telotron_splice_logos_CC/.done",
    output:
        touch(COMPOSITE_5P3P_LOGO_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p results/figures/composite_boundary_logos_5p3p
        python scripts/plot_composite_per_species.py \
            --outdir results/figures/composite_boundary_logos_5p3p \
            --panel "results/figures/interstitial_boundary_logos_GG:interstitial 5'=G 3'=G" \
            --panel "results/figures/interstitial_boundary_logos_GC:interstitial 5'=G 3'=C (hybrid)" \
            --panel "results/figures/interstitial_boundary_logos_CG:interstitial 5'=C 3'=G (hybrid)" \
            --panel "results/figures/interstitial_boundary_logos_CC:interstitial 5'=C 3'=C" \
            --panel "results/figures/non_telotron_splice_logos:non-telotron control splice signal" \
            --panel "results/figures/telotron_splice_logos_GG:telotron 5'=G 3'=G" \
            --panel "results/figures/telotron_splice_logos_GC:telotron 5'=G 3'=C (hybrid)" \
            --panel "results/figures/telotron_splice_logos_CG:telotron 5'=C 3'=G (hybrid)" \
            --panel "results/figures/telotron_splice_logos_CC:telotron 5'=C 3'=C"
        """


# Bundle all final TSVs + figures + manifest into a single zip for sharing.
rule package:
    input:
        "results/final_telotron_set.tsv",
        "results/final_telotron_set_dedup.tsv",
        "results/final_telotron_set_architecture.tsv",
        "results/final_species_summary.tsv",
        "results/final_negative_controls.tsv",
        "results/boundary_kmer_enrichment.tsv",
        "results/boundary_kmers_by_architecture.tsv",
        "results/distance_to_end.tsv",
        "results/architecture_summary.tsv",
        "results/dedup_log.tsv",
        "results/interstitial_arrays.tsv",
        "results/linker_blast_hits_own_genome.tsv",
        "results/linker_blast_hits_all_genomes.tsv",
        "results/figures/boundary_kmers_combined.pdf",
        *FIGURES,
        "manifests/all_genomes.tsv",
    output:
        "results/telotron_pipeline_outputs.zip",
    shell:
        "zip -qj {output} {input}"
